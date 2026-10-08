"""Public WorkGraph builder for surface-energy calculations."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from psteros.backends import add_qe_task, add_vasp_task
from psteros.backends.qe import as_aiida_structure
from psteros.config import QeCalculationConfig, SurfaceWorkflowConfig


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
            from psteros.backends.vasp import as_vasp_structure

            initial_structure = as_vasp_structure(source_structure)
            relax_task = add_vasp_task(
                workgraph, label=f"{label}_relax", structure=initial_structure,
                config=relaxation.calculation, execution=relaxation.execution,
                override=relaxation.role_overrides.get(label),
            )
            static_task = add_vasp_task(
                workgraph, label=f"{label}_static", structure=relax_task.outputs.structure,
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
