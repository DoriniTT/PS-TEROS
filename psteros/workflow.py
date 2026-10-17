"""Public WorkGraph builder for surface-energy calculations."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from psteros.backends import add_cp2k_task, add_qe_task, add_vasp_task
from psteros.backends.qe import as_aiida_structure
from psteros.config import Cp2kCalculationConfig, QeCalculationConfig, SurfaceWorkflowConfig


def _refuse_cp2k(*configs: SurfaceWorkflowConfig) -> None:
    """The QE/VASP builders branch on ``backend == "vasp"``; keep a CP2K recipe out of their QE branch."""

    if any(config.backend == "cp2k" for config in configs):
        raise ValueError(
            "CP2K recipes are built with psteros.build_cp2k_surface_workgraph or "
            "psteros.build_cp2k_relax_static_workgraph"
        )


def _check_labels(structures: Mapping[str, Any]) -> None:
    if not structures:
        raise ValueError("at least one labelled structure is required")
    invalid = [label for label in structures if not label or not label.replace("_", "").replace("-", "").isalnum()]
    if invalid:
        raise ValueError(f"invalid calculation labels: {invalid}")


def build_surface_workgraph(
    structures: Mapping[str, Any],
    config: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build or submit a serial AiiDA WorkGraph for labelled structures.

    The same typed surface recipe supports VASP (the central backend) and
    Quantum ESPRESSO.  Each label becomes one backend task.  The number of
    calculations running at once is ``ExecutionPolicy.max_concurrent_jobs``
    (one by default).

    Parameters
    ----------
    structures
        Mapping from stable calculation labels to pymatgen ``Structure``,
        AiiDA ``StructureData``, or an AiiDA node PK.
    config
        A typed backend and execution recipe.
    submit
        Submit the completed graph when true; otherwise return it unsubmitted.
    """

    if not structures:
        raise ValueError("at least one labelled structure is required")
    invalid = [label for label in structures if not label or not label.replace("_", "").replace("-", "").isalnum()]
    if invalid:
        raise ValueError(f"invalid calculation labels: {invalid}")
    _refuse_cp2k(config)

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=config.name)
    if config.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = config.execution.max_concurrent_jobs
    for label, structure in structures.items():
        override = config.role_overrides.get(label)
        if config.backend == "qe":
            task = add_qe_task(
                workgraph,
                label=label,
                structure=structure,
                config=config.calculation,
                execution=config.execution,
                override=override,
            )
            workgraph.outputs.__setattr__(f"{label}_parameters", task.outputs.output_parameters)
            workgraph.outputs.__setattr__(f"{label}_structure", task.outputs.output_structure)
            workgraph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
        else:
            task = add_vasp_task(
                workgraph,
                label=label,
                structure=structure,
                config=config.calculation,
                execution=config.execution,
                override=override,
            )
            workgraph.outputs.__setattr__(f"{label}_misc", task.outputs.misc)
            workgraph.outputs.__setattr__(f"{label}_structure", task.outputs.structure)
            workgraph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
    if submit:
        workgraph.submit()
    return workgraph


def _relaxation_check(incar: Mapping[str, Any], label: str, stage: str) -> None:
    """Refuse a VASP relaxation that does not move ions, or a static SCF that does."""

    values = {str(key).upper(): value for key, value in incar.items()}
    moves = int(values.get("IBRION", -1)) >= 0 and int(values.get("NSW", 0)) > 0
    if stage == "relax" and not moves:
        raise ValueError(f"{label}: the relaxation INCAR needs IBRION >= 0 and NSW > 0")
    if stage == "static" and moves:
        raise ValueError(f"{label}: the static INCAR must not relax (use NSW = 0 or IBRION = -1)")


def build_relax_static_workgraph(
    structures: Mapping[str, Any],
    relaxation: SurfaceWorkflowConfig,
    static: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build a serial relaxation-to-static WorkGraph for each structure, with VASP or QE.

    Each static calculation starts from the structure relaxed by the
    preceding relaxation; the static energy is the one to analyse (a
    relaxation with a changing cell has a basis-set error). Both recipes
    share one backend and one execution policy, whose
    ``max_concurrent_jobs`` limits the calculations running at once.

    Graph outputs per label: ``<label>_relaxed_structure`` and, for VASP,
    ``<label>_relax_misc`` and ``<label>_static_misc`` (energies), for QE
    ``<label>_relax_parameters`` and ``<label>_static_parameters``, plus the
    ``_retrieved`` folders of both stages. :func:`psteros.read_vasp_results`
    and :func:`psteros.read_qe_results` read them.

    With VASP, every relaxation INCAR (recipe plus override) must move the
    ions (``IBRION >= 0``, ``NSW > 0``) and every static INCAR must not.
    With QE the relaxations run as
    :class:`~psteros.backends.qe_workchains.PwRelaxStageWorkChain`, which
    accepts a ``vc-relax`` whose final SCF exceeded the thresholds (exit 501);
    psteros must then be installed in the environment of the AiiDA daemon.
    """

    if not structures:
        raise ValueError("at least one labelled structure is required")
    _refuse_cp2k(relaxation, static)
    if relaxation.backend != static.backend:
        raise ValueError("relaxation and static recipes must use the same backend")
    if relaxation.execution != static.execution:
        raise ValueError("relaxation and static recipes must share one execution policy")
    for label in structures:
        if not label or not label.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"invalid calculation label: {label!r}")
    if relaxation.backend == "vasp":
        from psteros.backends.vasp import vasp_incar

        for label in structures:
            _relaxation_check(vasp_incar(relaxation.calculation, relaxation.role_overrides.get(label)), label, "relax")
            _relaxation_check(vasp_incar(static.calculation, static.role_overrides.get(label)), label, "static")

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=f"{relaxation.name}_relax_static")
    if relaxation.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = relaxation.execution.max_concurrent_jobs
    for label, source_structure in structures.items():
        if relaxation.backend == "vasp":
            from psteros.backends.vasp import as_vasp_structure, has_custom_kinds

            initial_structure = as_vasp_structure(source_structure)
            relax_task = add_vasp_task(
                workgraph, label=f"{label}_relax", structure=initial_structure,
                config=relaxation.calculation, execution=relaxation.execution,
                override=relaxation.role_overrides.get(label),
            )
            static_structure = relax_task.outputs.structure
            if has_custom_kinds(initial_structure):
                # aiida-vasp returns the relaxed structure with element names only; a kind such as
                # the pseudo-hydrogen H0p75 must be restored or the static step has no POTCAR for it.
                from psteros.backends.vasp_kinds import restore_vasp_kinds

                static_structure = workgraph.add_task(
                    restore_vasp_kinds, name=f"{label}_restore_kinds",
                    relaxed=relax_task.outputs.structure, template=initial_structure,
                ).outputs.result
            static_task = add_vasp_task(
                workgraph, label=f"{label}_static", structure=static_structure,
                kinds_structure=initial_structure, config=static.calculation, execution=static.execution,
                override=static.role_overrides.get(label),
            )
            relaxed, energies = relax_task.outputs.structure, ("misc", "misc")
        else:
            initial_structure = as_aiida_structure(source_structure)
            relax_task = add_qe_task(
                workgraph, label=f"{label}_relax", structure=initial_structure,
                config=relaxation.calculation, execution=relaxation.execution,
                override=relaxation.role_overrides.get(label), relax_stage=True,
            )
            static_task = add_qe_task(
                workgraph, label=f"{label}_static", structure=relax_task.outputs.output_structure,
                pseudo_structure=initial_structure, config=static.calculation, execution=static.execution,
                override=static.role_overrides.get(label),
            )
            relaxed, energies = relax_task.outputs.output_structure, ("parameters", "output_parameters")
        name, port = energies
        workgraph.outputs.__setattr__(f"{label}_relaxed_structure", relaxed)
        workgraph.outputs.__setattr__(f"{label}_relax_{name}", getattr(relax_task.outputs, port))
        workgraph.outputs.__setattr__(f"{label}_relax_retrieved", relax_task.outputs.retrieved)
        workgraph.outputs.__setattr__(f"{label}_static_{name}", getattr(static_task.outputs, port))
        workgraph.outputs.__setattr__(f"{label}_static_retrieved", static_task.outputs.retrieved)
    if submit:
        workgraph.submit()
    return workgraph


def build_qe_relax_static_workgraph(
    structures: Mapping[str, Any],
    relaxation: SurfaceWorkflowConfig,
    static: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """:func:`build_relax_static_workgraph` restricted to Quantum ESPRESSO recipes."""

    if not structures:
        raise ValueError("at least one labelled structure is required")
    if relaxation.backend != "qe" or static.backend != "qe":
        raise TypeError("build_qe_relax_static_workgraph requires QE recipes")
    if not isinstance(relaxation.calculation, QeCalculationConfig):
        raise TypeError("relaxation requires QeCalculationConfig")
    if not isinstance(static.calculation, QeCalculationConfig):
        raise TypeError("static requires QeCalculationConfig")
    return build_relax_static_workgraph(structures, relaxation, static, submit=submit)


def _cp2k_recipe(config: SurfaceWorkflowConfig, builder: str) -> Cp2kCalculationConfig:
    if config.backend != "cp2k" or not isinstance(config.calculation, Cp2kCalculationConfig):
        raise TypeError(f"{builder} requires CP2K recipes (backend='cp2k', Cp2kCalculationConfig)")
    return config.calculation


def build_cp2k_surface_workgraph(
    structures: Mapping[str, Any],
    config: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build or submit a WorkGraph with one CP2K calculation per labelled structure.

    The run type is the recipe's ``GLOBAL/RUN_TYPE`` (``ENERGY_FORCE`` when
    left out), so this is the CP2K counterpart of
    :func:`build_surface_workgraph` for static calculations. At most
    ``ExecutionPolicy.max_concurrent_jobs`` calculations run at once.

    Graph outputs per label: ``<label>_parameters`` (aiida-cp2k output, energy
    in Hartree), ``<label>_structure`` (when CP2K writes one) and
    ``<label>_retrieved``, which holds the force and density-of-states files
    of an ``ENERGY_FORCE`` run. :func:`psteros.read_cp2k_results`,
    :func:`psteros.read_cp2k_forces` and :func:`psteros.read_cp2k_dos` read
    them in eV and eV/A.
    """

    _check_labels(structures)
    calculation = _cp2k_recipe(config, "build_cp2k_surface_workgraph")

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=config.name)
    if config.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = config.execution.max_concurrent_jobs
    for label, structure in structures.items():
        task = add_cp2k_task(
            workgraph, label=label, structure=structure, config=calculation,
            execution=config.execution, override=config.role_overrides.get(label),
        )
        workgraph.outputs.__setattr__(f"{label}_parameters", task.outputs.output_parameters)
        workgraph.outputs.__setattr__(f"{label}_structure", task.outputs.output_structure)
        workgraph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
    if submit:
        workgraph.submit()
    return workgraph


def build_cp2k_relax_static_workgraph(
    structures: Mapping[str, Any],
    relaxation: SurfaceWorkflowConfig,
    static: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build a relaxation-to-static WorkGraph for each structure with CP2K.

    The CP2K counterpart of :func:`build_relax_static_workgraph`: each static
    calculation starts from the structure of the preceding relaxation, and its
    energy is the one to analyse. Every relaxation (recipe plus override)
    must have ``GLOBAL/RUN_TYPE`` ``GEO_OPT`` or ``CELL_OPT``, every static
    one ``ENERGY_FORCE`` or ``ENERGY`` (or none). ``fixed_sites`` of the
    relaxation overrides become ``FIXED_ATOMS``. Both recipes share one
    execution policy.

    Graph outputs per label, in the pattern of the QE builder:
    ``<label>_relaxed_structure``, ``<label>_relax_parameters``,
    ``<label>_static_parameters`` and the ``_retrieved`` folders of both
    stages; :func:`psteros.read_cp2k_results` reads them.
    """

    from psteros.backends.cp2k import RELAX_RUN_TYPES, STATIC_RUN_TYPES, cp2k_run_type

    _check_labels(structures)
    relax_recipe = _cp2k_recipe(relaxation, "build_cp2k_relax_static_workgraph")
    static_recipe = _cp2k_recipe(static, "build_cp2k_relax_static_workgraph")
    if relaxation.execution != static.execution:
        raise ValueError("relaxation and static recipes must share one execution policy")
    for label in structures:
        relax_type = cp2k_run_type(relax_recipe, relaxation.role_overrides.get(label))
        if relax_type not in RELAX_RUN_TYPES:
            raise ValueError(f"{label}: the CP2K relaxation needs GLOBAL/RUN_TYPE GEO_OPT or CELL_OPT, got {relax_type}")
        static_type = cp2k_run_type(static_recipe, static.role_overrides.get(label))
        if static_type not in STATIC_RUN_TYPES:
            raise ValueError(f"{label}: the CP2K static calculation must not move atoms, got RUN_TYPE {static_type}")

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=f"{relaxation.name}_relax_static")
    if relaxation.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = relaxation.execution.max_concurrent_jobs
    for label, source_structure in structures.items():
        initial_structure = as_aiida_structure(source_structure)
        relax_task = add_cp2k_task(
            workgraph, label=f"{label}_relax", structure=initial_structure, config=relax_recipe,
            execution=relaxation.execution, override=relaxation.role_overrides.get(label),
        )
        static_task = add_cp2k_task(
            workgraph, label=f"{label}_static", structure=relax_task.outputs.output_structure,
            kinds_structure=initial_structure, config=static_recipe, execution=static.execution,
            override=static.role_overrides.get(label),
        )
        workgraph.outputs.__setattr__(f"{label}_relaxed_structure", relax_task.outputs.output_structure)
        workgraph.outputs.__setattr__(f"{label}_relax_parameters", relax_task.outputs.output_parameters)
        workgraph.outputs.__setattr__(f"{label}_relax_retrieved", relax_task.outputs.retrieved)
        workgraph.outputs.__setattr__(f"{label}_static_parameters", static_task.outputs.output_parameters)
        workgraph.outputs.__setattr__(f"{label}_static_retrieved", static_task.outputs.retrieved)
    if submit:
        workgraph.submit()
    return workgraph
