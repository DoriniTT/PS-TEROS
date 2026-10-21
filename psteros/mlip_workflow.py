"""WorkGraph of the MLIP route to a surface phase diagram.

:func:`build_mlip_surface_workgraph` runs the route (``docs/MLIP_SURFACES.md``):

1. one DFT static calculation (energy and forces) of every structure, plus
   optional extra training structures such as :func:`psteros.rattled_copies`;
2. a replay fine-tune of a MACE foundation model on those labels (one model
   per seed of :class:`~psteros.MaceFineTuneConfig`; or of a UMA model with
   :class:`~psteros.UmaFineTuneConfig`, or of a NequIP or Allegro model with
   :class:`~psteros.NequipFineTuneConfig`), on a GPU computer;
3. the model of the first seed relaxes every structure (cell fixed);
4. ``active_learning_rounds`` times: one DFT static at each relaxed geometry
   is added to the labels, the model is fine-tuned again from the foundation
   model on all labels and relaxes the structures again from the start;
5. with ``dft_check``, one DFT static at each final geometry;
6. with ``vibrations``, the harmonic modes of every structure from the
   deployed model at its relaxed geometry (one more GPU job).

The energies of one diagram come from one source (see
:func:`psteros.check_one_energy_source`): ``"mlip_dft"``, the DFT statics of
step 5 (the default), or ``"mlip"``, the model at its own geometries.
:func:`read_mlip_energies` returns one of
them for every label; :func:`read_mlip_results` the relaxations.

:class:`MlipSurfaceStudy` is the route for your own slabs of a binary
compound from structures to the phase diagram: it names the calculations,
fixes the slab centres, builds the graph and turns its energies into a
:class:`psteros.SurfacePhaseDiagram`. :class:`MlipRoute` runs the slabs that
psteros cuts with the right electron count (:class:`psteros.ChargeNeutralSurfaceStudy`,
:class:`psteros.PolarSurfaceStudy`), so it serves semiconductors and
insulators as well as metals.

Why the metal check: in a metal the electron count of a slab needs no care,
any cut of the bulk will do. A slab of a semiconductor or an insulator must
be charge-neutral (or passivated), which pymatgen's ``SlabGenerator`` alone
does not ensure. Hand-cut slabs therefore need a metallic bulk unless
``allow_non_metallic=True``; the slabs of :class:`MlipRoute` are cut with
electron counting and need no check.
"""

from __future__ import annotations

import json
import re
import tempfile
import warnings
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride, ExecutionPolicy, SurfaceWorkflowConfig
from psteros.mlip import MaceFineTuneConfig, MlipRelaxation, MlipRelaxConfig
from psteros.nequip import NequipFineTuneConfig
from psteros.uma import UmaFineTuneConfig
from psteros.surface_study import ChargeNeutralStudyResult
from psteros.vibrations_workflow import VibrationsConfig

MLIP_ENERGY_SOURCES = ("mlip", "mlip_dft")

#: Warning of a graph whose bulk is not checked to be a metal.
NON_METALLIC_WARNING = (
    "allow_non_metallic=True: the bulk is not checked to be a metal. The slabs of a semiconductor or an insulator "
    "need the right electron count (charge-neutral or passivated), which pymatgen's SlabGenerator alone does not "
    "ensure; psteros.MlipRoute with a ChargeNeutralSurfaceStudy or a PolarSurfaceStudy cuts such slabs"
)


@dataclass(frozen=True)
class MlipJobConfig:
    """Where the MACE (or UMA, NequIP, Allegro) fine-tunes and relaxations run.

    ``code_label`` is an AiiDA code whose executable is the Python of an
    environment with ase, torch and mace-torch on the GPU computer that holds
    ``MaceFineTuneConfig.foundation_model`` (for a :class:`psteros.UmaFineTuneConfig`,
    an environment with fairchem-core instead; for a :class:`psteros.NequipFineTuneConfig`,
    one with nequip and nequip-allegro); ``execution`` gives the scheduler
    settings of one job (one GPU; ``with_mpi=False``, the jobs run one Python
    process). Each fine-tune and each relaxation is one job; a relaxation
    links the remote folder of its fine-tune, where the models stay, so that
    folder is never cleaned. ``device`` is ``"cuda"`` or ``"cpu"``. The jobs
    write their stderr (MACE's log, torch warnings) into stdout
    (``redirect_stderr``), unless ``execution.extra_options`` says otherwise.
    """

    code_label: str
    execution: ExecutionPolicy
    device: str = "cuda"

    def __post_init__(self) -> None:
        if not self.code_label or not self.code_label.strip():
            raise ValueError("MlipJobConfig: code_label must name an AiiDA code")
        if self.execution.with_mpi:
            raise ValueError("MlipJobConfig: the MACE jobs run one Python process; use ExecutionPolicy(with_mpi=False)")
        if self.device not in ("cuda", "cpu"):
            raise ValueError(f"MlipJobConfig: device must be 'cuda' or 'cpu', got {self.device!r}")


def _pymatgen(label: str, structure: Any) -> Any:
    """A plain pymatgen ``Structure`` of a pymatgen structure, ``StructureData`` or PK known at build time."""

    from pymatgen.core import Structure

    if isinstance(structure, int):
        from aiida import orm

        structure = orm.load_node(structure)
    if hasattr(structure, "get_pymatgen_structure"):
        structure = structure.get_pymatgen_structure()
    if not isinstance(structure, Structure):
        raise TypeError(f"{label}: give a pymatgen Structure, a StructureData or its PK (known when the graph is built)")
    return Structure(structure.lattice, structure.species, structure.frac_coords)


_PSEUDO_HYDROGEN_KIND = re.compile(r"H\d+p\d+")  # H0p75, H1p25, ...: the kinds of psteros.PseudoHydrogen


def _hydrogen_kinds(structure: Any) -> set[str]:
    """The hydrogen kinds of a structure: its pseudo-hydrogen kind names (``H0p75``, ...) and ``"H"`` for hydrogen.

    ``structure`` is a pymatgen ``Structure`` (kinds in the ``kind_name`` site property), a ``StructureData`` or its
    PK. Hydrogen kinds with any other name (``H``, ``H1``, ...) are hydrogen.
    """

    if isinstance(structure, int):
        from aiida import orm

        structure = orm.load_node(structure)
    if hasattr(structure, "get_kind"):
        names = [site.kind_name for site in structure.sites if structure.get_kind(site.kind_name).symbol == "H"]
    else:
        kinds = structure.site_properties.get("kind_name") or [None] * len(structure)
        names = [kind for site, kind in zip(structure, kinds) if site.specie.symbol == "H"]
    return {name if name and _PSEUDO_HYDROGEN_KIND.fullmatch(name) else "H" for name in names}


def _check_hydrogen_kinds(structures: Mapping[str, Any]) -> None:
    """Refuse hydrogen together with pseudo-hydrogens: the models see elements, not kinds, so both are one species
    ``H`` and can sit in the same environments (an H and an H0p75 on an As atom), with different DFT labels.

    Several pseudo-hydrogen kinds are allowed: each passivates one element, so the model tells them apart by their
    neighbour. The error names the labels and their hydrogen kinds.
    """

    found = {label: _hydrogen_kinds(structure) for label, structure in structures.items()}
    kinds = set().union(*found.values())
    if "H" in kinds and kinds != {"H"}:
        detail = "; ".join(f"{label}: {', '.join(sorted(item))}" for label, item in found.items() if item)
        raise ValueError(
            "the MLIP route cannot hold hydrogen and pseudo-hydrogens together: its models see elements, not kinds, "
            f"so both would be one species H ({detail}). Run this study with the DFT route"
        )


def _replay_elements(structures: Mapping[str, Any], plain: Sequence[Any]) -> list[str]:
    """The elements that the replay structures may hold: those of the study, without hydrogen when every hydrogen
    of ``structures`` is a pseudo-hydrogen.

    The models see elements, not kinds, so a pseudo-hydrogen is hydrogen for them; replay structures with real
    hydrogen, labelled by the foundation model, would train that species on hydrogen chemistry that no structure of
    a passivated study has. ``plain`` are the structures as plain pymatgen structures.
    """

    elements = {site.specie.symbol for item in plain for site in item}
    hydrogen = set().union(*(_hydrogen_kinds(structure) for structure in structures.values()))
    if hydrogen and "H" not in hydrogen:
        elements.discard("H")
    return sorted(elements)


def _check_static(static: SurfaceWorkflowConfig, labels: Sequence[str]) -> None:
    """Every static calculation must give forces and must not move atoms."""

    if static.backend == "cp2k":
        from psteros.backends.cp2k import cp2k_run_type

        for label in labels:
            run_type = cp2k_run_type(static.calculation, static.role_overrides.get(label))
            if run_type != "ENERGY_FORCE":
                raise ValueError(f"{label}: the CP2K statics of the MLIP route need RUN_TYPE ENERGY_FORCE, got {run_type}")
    elif static.backend == "vasp":
        from psteros.backends.vasp import vasp_incar
        from psteros.workflow import _relaxation_check

        for label in labels:
            _relaxation_check(vasp_incar(static.calculation, static.role_overrides.get(label)), label, "static")
    else:
        for label in labels:
            override = static.role_overrides.get(label)
            control = dict(dict(static.calculation.parameters).get("CONTROL", {}))
            if override and override.parameters:
                control.update(dict(override.parameters).get("CONTROL", {}))
            if control.get("calculation", "scf") != "scf":
                raise ValueError(f"{label}: the QE statics of the MLIP route need calculation 'scf', got {control['calculation']!r}")


def _static_override(static: SurfaceWorkflowConfig, label: str) -> CalculationOverride:
    """The override of a static: no fixed sites, and forces printed by QE."""

    override = static.role_overrides.get(label)
    override = replace(override, fixed_sites=()) if override else CalculationOverride()
    if static.backend == "qe":
        parameters = {namelist: dict(values) for namelist, values in dict(override.parameters or {}).items()}
        parameters.setdefault("CONTROL", {})["tprnfor"] = True
        override = replace(override, parameters=parameters)
    return override


def _add_static(workgraph: Any, static: SurfaceWorkflowConfig, name: str, label: str, structure: Any, reference: Any):
    """Add one DFT static of ``label``; return its (parsed output, retrieved folder) sockets."""

    override = _static_override(static, label)
    if static.backend == "cp2k":
        from psteros.backends.cp2k import add_cp2k_task

        node = add_cp2k_task(workgraph, label=name, structure=structure, kinds_structure=reference,
                             config=static.calculation, execution=static.execution, override=override)
        return node.outputs.output_parameters, node.outputs.retrieved
    if static.backend == "vasp":
        from psteros.backends.vasp import add_vasp_task

        node = add_vasp_task(workgraph, label=name, structure=structure, kinds_structure=reference,
                             config=static.calculation, execution=static.execution, override=override)
        return node.outputs.misc, node.outputs.retrieved
    from psteros.backends.qe import add_qe_task

    node = add_qe_task(workgraph, label=name, structure=structure, pseudo_structure=reference,
                       config=static.calculation, execution=static.execution, override=override)
    return node.outputs.output_parameters, node.outputs.retrieved


def build_mlip_surface_workgraph(
    structures: Mapping[str, Any],
    static: SurfaceWorkflowConfig,
    fine_tune: MaceFineTuneConfig | UmaFineTuneConfig | NequipFineTuneConfig,
    relax: MlipRelaxConfig = MlipRelaxConfig(),
    *,
    job: MlipJobConfig,
    metallic_bulk: str,
    fixed_sites: Mapping[str, Sequence[int]] | None = None,
    training_structures: Mapping[str, Any] | None = None,
    active_learning_rounds: int = 1,
    dft_check: bool = True,
    vibrations: VibrationsConfig | None = None,
    allow_non_metallic: bool = False,
    submit: bool = False,
) -> Any:
    """Build or submit the WorkGraph of the MLIP route (see the module docstring).

    Parameters
    ----------
    structures
        Labels -> structures (pymatgen ``Structure``, ``StructureData`` or PK)
        of every energy of the diagram: the slabs and the references (bulk,
        metal, gas molecule in a box), the references already relaxed with
        the DFT recipe. Every structure is relaxed by the model (cell fixed).
    static
        DFT static recipe of every calculation (QE, VASP or CP2K): forces are
        needed (CP2K ``ENERGY_FORCE``; QE runs with ``tprnfor``), and
        ``role_overrides`` apply by label (the extra training structures use
        the override of their own label, if any).
    fine_tune, relax
        The MACE recipe and the relaxation settings. With replay, the replay
        geometries are chosen from ``fine_tune.replay_file`` (a local file)
        when the graph is built, for the elements of ``structures``.
    job
        Code and scheduler of the GPU jobs (:class:`MlipJobConfig`).
    metallic_bulk
        Label of the bulk that must be a metal (partly filled states in its
        static); the first fine-tune stops otherwise, unless
        ``allow_non_metallic``.
    fixed_sites
        Labels -> zero-based sites that the model keeps in place, e.g.
        :func:`psteros.central_sites` of each slab.
    training_structures
        Extra labels -> structures with a DFT static in the training set and no
        relaxation, e.g. :func:`psteros.rattled_copies` of the slabs.
    active_learning_rounds
        Number of rounds of step 4 (1 by default; 0 skips them).
    dft_check
        Run the DFT statics of step 5 (needed for the ``"mlip_dft"`` energies).
    vibrations
        A :class:`psteros.VibrationsConfig`: the deployed model of the last
        round gives the Gamma-point harmonic modes of every structure at its
        relaxed geometry, by the finite differences of
        :func:`psteros.build_vibrations_workgraph` (one GPU job). The
        ``fixed_sites`` are not displaced (their region is counted as bulk);
        a bulk is displaced in its ``supercells`` entry; ``molecules`` lose
        their zero modes. ``None`` (the default) adds nothing.
    allow_non_metallic
        Do not check that ``metallic_bulk`` is a metal. For a semiconductor or
        an insulator, the slabs must then have the right electron count
        (a warning says so); see :class:`MlipRoute`.
    submit
        Submit the graph when true; otherwise return it unsubmitted.

    Graph outputs (all new; no other builder changes): per label
    ``<label>_static_parameters`` and ``<label>_static_retrieved`` (step 1),
    ``<label>_mlip_relaxed_structure`` and ``<label>_mlip_energy`` (eV, the
    first seed's model of the last round), with ``dft_check``
    ``<label>_mlip_dft_energy`` (eV), ``<label>_mlip_dft_parameters`` and
    ``<label>_mlip_dft_retrieved``; per round ``r`` (0 is the first fine-tune)
    ``mlip_round<r>_finetune_summary``, ``mlip_round<r>_model_folder`` and
    ``mlip_round<r>_relax_results`` (every seed's relaxations); with
    ``vibrations`` ``<label>_mlip_vibrations`` (read by
    :func:`psteros.read_vibrations`) and ``mlip_vibrations_results``.
    """

    if not isinstance(allow_non_metallic, bool):
        raise TypeError("allow_non_metallic must be True or False")
    if allow_non_metallic:
        warnings.warn(NON_METALLIC_WARNING, UserWarning, stacklevel=2)
    return _mlip_workgraph(
        structures, static, fine_tune, relax, job=job, metallic_bulk=metallic_bulk, metal_check=not allow_non_metallic,
        fixed_sites=fixed_sites, training_structures=training_structures, active_learning_rounds=active_learning_rounds,
        dft_check=dft_check, vibrations=vibrations, submit=submit,
    )


def _mlip_workgraph(
    structures: Mapping[str, Any],
    static: SurfaceWorkflowConfig,
    fine_tune: MaceFineTuneConfig | UmaFineTuneConfig | NequipFineTuneConfig,
    relax: MlipRelaxConfig,
    *,
    job: MlipJobConfig,
    metallic_bulk: str,
    metal_check: bool,
    fixed_sites: Mapping[str, Sequence[int]] | None,
    training_structures: Mapping[str, Any] | None,
    active_learning_rounds: int,
    dft_check: bool,
    vibrations: VibrationsConfig | None,
    submit: bool,
) -> Any:
    """The graph of :func:`build_mlip_surface_workgraph`; ``metal_check=False`` skips the metal check silently."""

    from psteros.mlip import select_replay_structures
    from psteros.workflow import _check_labels

    _check_labels(structures)
    extra = dict(training_structures or {})
    if extra:
        _check_labels(extra)
    _check_hydrogen_kinds({**structures, **extra})
    if not isinstance(job, MlipJobConfig):
        raise TypeError("job must be an MlipJobConfig")
    families = {MaceFineTuneConfig: "mace", UmaFineTuneConfig: "uma", NequipFineTuneConfig: "nequip"}
    if type(fine_tune) not in families or not isinstance(relax, MlipRelaxConfig):
        raise TypeError("fine_tune must be a MaceFineTuneConfig, a UmaFineTuneConfig or a NequipFineTuneConfig and "
                        "relax an MlipRelaxConfig")
    family = families[type(fine_tune)]
    # the head the deployed models run: MACE's fine-tuned head (relax.head), or the UMA task the fine-tune kept;
    # a NequIP-framework model has none
    head = fine_tune.task if family == "uma" else relax.head
    suffix = ".nequip.zip" if family == "nequip" else ".model"  # the file the fine-tune job writes per seed
    if isinstance(active_learning_rounds, bool) or not isinstance(active_learning_rounds, int) or active_learning_rounds < 0:
        raise ValueError("active_learning_rounds must be a non-negative integer")
    if metallic_bulk not in structures:
        raise ValueError(f"metallic_bulk {metallic_bulk!r} is not a label of structures")
    overlap = sorted(set(extra).intersection(structures))
    if overlap:
        raise ValueError(f"training_structures repeat labels of structures: {overlap}")
    frame_labels = list(structures) + list(extra)
    taken = sorted(
        f"{label}_al{round_}" for label in structures for round_ in range(1, active_learning_rounds + 1)
        if f"{label}_al{round_}" in frame_labels
    )
    if taken:
        raise ValueError(f"labels of the active-learning frames are already used: {taken}")
    fixed_sites = {label: tuple(sorted({int(site) for site in sites})) for label, sites in dict(fixed_sites or {}).items()}
    unknown = sorted(set(fixed_sites).difference(structures))
    if unknown:
        raise ValueError(f"fixed_sites names labels that are not in structures: {unknown}")
    _check_static(static, frame_labels)
    plain = {label: _pymatgen(label, structure) for label, structure in structures.items()}
    plain_extra = {label: _pymatgen(label, structure) for label, structure in extra.items()}
    if fine_tune.foundation_name is not None:
        from psteros.mace_models import foundation_level_warnings
        from psteros.nequip import nequip_level_warnings
        from psteros.uma import uma_level_warnings

        present = {site.specie.symbol for item in (*plain.values(), *plain_extra.values()) for site in item}
        if family == "uma":
            messages = uma_level_warnings(fine_tune.foundation_name, fine_tune.task, static.calculation, present)
        elif family == "nequip":
            messages = nequip_level_warnings(fine_tune.foundation_name, static.calculation, present)
        else:
            messages = foundation_level_warnings(fine_tune.foundation_name, fine_tune.foundation_head,
                                                 static.calculation, present)
        for message in messages:
            warnings.warn(f"{type(fine_tune).__name__}: {message}", UserWarning, stacklevel=3)
    for label, sites in fixed_sites.items():
        if any(site < 0 or site >= len(plain[label]) for site in sites):
            raise ValueError(f"{label}: fixed_sites must be within [0, {len(plain[label])})")
    displaced: dict[str, tuple[int, ...]] = {}
    if vibrations is not None:
        if not isinstance(vibrations, VibrationsConfig):
            raise TypeError("vibrations must be a VibrationsConfig")
        unknown = sorted(
            (set(vibrations.displaced_sites) | set(vibrations.supercells) | set(vibrations.molecules)).difference(structures)
        )
        if unknown:
            raise ValueError(f"VibrationsConfig names labels that are not in structures: {unknown}")
        for label, structure in plain.items():
            override = CalculationOverride(fixed_sites=fixed_sites.get(label, ()))
            displaced[label] = vibrations.sites_to_displace(label, len(structure), override)

    from aiida import orm
    from aiida_workgraph import WorkGraph, shelljob

    from psteros.backends.mlip_tasks import mlip_dft_energy, mlip_relaxed, mlip_training_folder
    from psteros.backends.qe import as_aiida_structure
    from psteros.mlip import prepare_relax_job

    code = orm.load_code(job.code_label)
    replay_nodes: dict[str, Any] = {}
    replay_record: dict[str, Any] = {}
    if fine_tune.replay:
        from ase.io import write

        elements = _replay_elements({**structures, **extra}, [*plain.values(), *plain_extra.values()])
        train, valid, replay_record = select_replay_structures(
            fine_tune.replay_file, elements=elements, min_atoms=fine_tune.replay_min_atoms,
            max_atoms=fine_tune.replay_max_atoms, validation_every=fine_tune.replay_validation_every,
        )
        with tempfile.TemporaryDirectory() as scratch:
            for key, atoms in (("replay_train", train), ("replay_valid", valid)):
                path = Path(scratch) / f"{key}_geometries.extxyz"
                write(path, list(atoms), format="extxyz")
                replay_nodes[key] = orm.SinglefileData(file=str(path), filename=path.name)

    workgraph = WorkGraph(name=f"{static.name}_mlip")
    if static.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = static.execution.max_concurrent_jobs
    start = {label: as_aiida_structure(structure) for label, structure in structures.items()}
    start.update({label: as_aiida_structure(structure) for label, structure in extra.items()})

    # inputs of the training folder: structure, parsed output and retrieved folder of every frame
    frames: dict[str, Any] = {}
    for label, structure in start.items():
        output, retrieved = _add_static(workgraph, static, f"{label}_static", label, structure, structure)
        frames.update({f"structure_{label}": structure, f"output_{label}": output, f"retrieved_{label}": retrieved})
        workgraph.outputs.__setattr__(f"{label}_static_parameters", output)
        workgraph.outputs.__setattr__(f"{label}_static_retrieved", retrieved)
    labels_of_frames = list(frame_labels)

    # MACE logs to stderr, and aiida-shell fails a job whose stderr is not empty (exit 410): merge it into stdout
    execution = {"redirect_stderr": True, **job.execution.scheduler_options()}
    deployed = f"seed{fine_tune.seeds[0]}"
    recipe = json.loads(json.dumps(asdict(fine_tune)))
    if recipe["foundation_name"] is None:  # left out: the settings of the graph stay as they were
        del recipe["foundation_name"]
    final: dict[str, Any] = {}
    for round_ in range(active_learning_rounds + 1):
        name = f"{static.name}_round{round_}"
        options = {
            "backend": static.backend, "name": name, "device": job.device, "metallic_bulk": metallic_bulk,
            "frames": list(labels_of_frames), "fine_tune": recipe,
            "replay_selection": replay_record,
        }
        if not metal_check:
            options["metal_check"] = False
        if family != "mace":  # left out for MACE: the settings of MACE graphs stay as they were
            options["family"] = family
        settings = orm.Dict(options)
        folder = workgraph.add_task(
            mlip_training_folder, name=f"mlip_round{round_}_training", settings=settings,
            inputs={**frames, **replay_nodes},
        )
        finetune = workgraph.add_task(
            shelljob, name=f"mlip_round{round_}_finetune", command=code, arguments=["-u", "finetune.py"],
            nodes={"job": folder.outputs.result}, outputs=["summary.json"],
            metadata={"label": f"{name}_finetune", "options": dict(execution)},
        )
        models = {f"seed{seed}": (f"finetune/models/{name}_seed{seed}{suffix}", head) for seed in fine_tune.seeds}
        with tempfile.TemporaryDirectory() as scratch:
            relax_folder = orm.FolderData(tree=str(prepare_relax_job(
                Path(scratch) / "relax", plain, models, relax, fixed_sites=fixed_sites, device=job.device,
                family=family)))
        relaxation = workgraph.add_task(
            shelljob, name=f"mlip_round{round_}_relax", command=code, arguments=["-u", "relax.py"],
            nodes={"job": relax_folder, "finetune": finetune.outputs.remote_folder}, filenames={"finetune": "finetune"},
            outputs=["results.json", "relaxed"],
            metadata={"label": f"{name}_relax", "options": {**execution, "use_symlinks": True}},
        )
        workgraph.outputs.__setattr__(f"mlip_round{round_}_finetune_summary", finetune.outputs.summary_json)
        workgraph.outputs.__setattr__(f"mlip_round{round_}_model_folder", finetune.outputs.remote_folder)
        workgraph.outputs.__setattr__(f"mlip_round{round_}_relax_results", relaxation.outputs.results_json)
        last = round_ == active_learning_rounds
        for label in structures:
            relaxed = workgraph.add_task(
                mlip_relaxed, name=f"{label}_mlip_round{round_}_relaxed", structure=start[label],
                results=relaxation.outputs.results_json, relaxed=relaxation.outputs.relaxed,
                settings=orm.Dict({"label": label, "model": deployed}),
            )
            if not last:
                frame = f"{label}_al{round_ + 1}"
                output, retrieved = _add_static(workgraph, static, f"{frame}_static", label,
                                                relaxed.outputs.structure, start[label])
                frames.update({f"structure_{frame}": relaxed.outputs.structure, f"output_{frame}": output,
                               f"retrieved_{frame}": retrieved})
                labels_of_frames.append(frame)
                continue
            final[label] = relaxed.outputs.structure
            workgraph.outputs.__setattr__(f"{label}_mlip_relaxed_structure", relaxed.outputs.structure)
            workgraph.outputs.__setattr__(f"{label}_mlip_energy", relaxed.outputs.energy)
            if dft_check:
                output, retrieved = _add_static(workgraph, static, f"{label}_mlip_dft", label,
                                                relaxed.outputs.structure, start[label])
                energy = workgraph.add_task(
                    mlip_dft_energy, name=f"{label}_mlip_dft_energy", output=output,
                    settings=orm.Dict({"backend": static.backend, "label": label}),
                )
                workgraph.outputs.__setattr__(f"{label}_mlip_dft_energy", energy.outputs.result)
                workgraph.outputs.__setattr__(f"{label}_mlip_dft_parameters", output)
                workgraph.outputs.__setattr__(f"{label}_mlip_dft_retrieved", retrieved)
    if vibrations is not None:
        _add_mlip_vibrations(workgraph, vibrations, final, displaced, plain, code=code, execution=execution,
                             model=(models[deployed][0], head), deployed=deployed, relax=relax, job=job,
                             finetune=finetune, name=static.name, family=family)
    if submit:
        workgraph.submit()
    return workgraph


def _add_mlip_vibrations(workgraph, vibrations, structures, displaced, plain, *, code, execution, model, deployed,
                         relax, job, finetune, name, family="mace") -> None:
    """The GPU job of the harmonic modes and one modes task per label (``<label>_mlip_vibrations``)."""

    from aiida import orm
    from aiida_workgraph import shelljob

    from psteros.backends.mlip_tasks import mlip_harmonic_modes, mlip_vibrations_folder
    from psteros.backends.vibration_tasks import repeat_structure

    targets, sites, sizes = {}, {}, {}
    for label, structure in structures.items():
        size = vibrations.supercells.get(label, (1, 1, 1))
        sizes[label] = size[0] * size[1] * size[2]
        sites[label] = list(displaced[label])
        if sizes[label] > 1:
            structure = workgraph.add_task(
                repeat_structure, name=f"{label}_mlip_supercell", structure=structure, size=orm.List(list(size))
            ).outputs.result
            sites[label] = list(range(len(plain[label]) * sizes[label]))
        targets[label] = structure
    options = {"model": list(model), "displaced_sites": sites, "device": job.device,
               "displacement_angstrom": vibrations.displacement_angstrom, "default_dtype": relax.default_dtype}
    if family != "mace":
        options["family"] = family
    folder = workgraph.add_task(
        mlip_vibrations_folder, name="mlip_vibrations_job", settings=orm.Dict(options),
        structures={f"structure_{label}": structure for label, structure in targets.items()},
    )
    forces = workgraph.add_task(
        shelljob, name="mlip_vibrations", command=code, arguments=["-u", "vibrations.py"],
        nodes={"job": folder.outputs.result, "finetune": finetune.outputs.remote_folder}, filenames={"finetune": "finetune"},
        outputs=["results.json", "vibrations"],
        metadata={"label": f"{name}_vibrations", "options": {**execution, "use_symlinks": True}},
    )
    workgraph.outputs.__setattr__("mlip_vibrations_results", forces.outputs.results_json)
    for label, structure in targets.items():
        modes = workgraph.add_task(
            mlip_harmonic_modes, name=f"{label}_mlip_vibrations", structure=structure, forces=forces.outputs.vibrations,
            settings=orm.Dict({"label": label, "model": deployed, "displaced_sites": sites[label],
                               "displacement_angstrom": vibrations.displacement_angstrom,
                               "supercell_size": sizes[label], "molecule": label in vibrations.molecules}),
        )
        workgraph.outputs.__setattr__(f"{label}_mlip_vibrations", modes.outputs.result)


def _task_node(graph: Any, name: str) -> Any:
    """The newest process that ``graph`` ran for the task ``name`` (None if it has not run)."""

    from psteros.backends.cp2k import _children

    return _children(graph).get(name)


def _relax_rounds(graph: Any) -> list[int]:
    import re

    from psteros.backends.cp2k import _children

    pattern = re.compile(r"mlip_round(\d+)_relax")
    return sorted(int(match.group(1)) for name in _children(graph) if (match := pattern.fullmatch(name)))


def read_mlip_results(
    graph: Any,
    labels: Sequence[str],
    *,
    seed: int | None = None,
    round: int | None = None,
) -> dict[str, MlipRelaxation]:
    """Relaxations of a :func:`build_mlip_surface_workgraph` graph, as :class:`~psteros.MlipRelaxation`.

    ``graph`` is the WorkGraph process node (``orm.load_node(pk)``); it may
    still run or have failed. ``round`` picks the fine-tune round (default:
    the last one that has a finished relaxation job) and ``seed`` the model
    (default: the first seed, the one the graph deploys).
    """

    import io

    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor

    rounds = _relax_rounds(graph)
    if round is None:
        finished = [r for r in rounds if getattr(_task_node(graph, f"mlip_round{r}_relax"), "is_finished_ok", False)]
        if not finished:
            raise ValueError("no relaxation job of this graph has finished")
        round = finished[-1]
    node = _task_node(graph, f"mlip_round{round}_relax")
    if node is None or not node.is_finished_ok:
        raise ValueError(f"the relaxation job of round {round} has not finished")
    results = json.loads(node.outputs.results_json.get_content())
    if seed is None:
        model = next(iter(results["models"]))
    else:
        model = f"seed{seed}"
        if model not in results["models"]:
            raise ValueError(f"no model of seed {seed} in round {round}; models: {sorted(results['models'])}")
    rows = results["models"][model]["structures"]
    missing = [label for label in labels if label not in rows]
    if missing:
        raise ValueError(f"no relaxation of {missing} in round {round}")
    output = {}
    for label in labels:
        row = rows[label]
        content = node.outputs.relaxed.base.repository.get_object_content(f"{model}/{label}.extxyz")
        atoms = read(io.StringIO(content), format="extxyz")
        output[label] = MlipRelaxation(
            structure=AseAtomsAdaptor.get_structure(atoms), energy_ev=float(row["energy_ev"]),
            max_force_ev_per_angstrom=float(row["max_force_ev_per_angstrom"]), steps=int(row["steps"]),
            converged=bool(row["converged"]),
        )
    return output


def read_mlip_energies(graph: Any, labels: Sequence[str], *, energy_source: str = "mlip_dft") -> dict[str, float]:
    """Energies (eV) of ``labels`` from one source of a :func:`build_mlip_surface_workgraph` graph.

    ``"mlip_dft"``: the DFT statics at the final geometries (``dft_check``);
    ``"mlip"``: the deployed model at its own geometries in the last round.
    Every label comes from the same source, so slabs and references can go
    into one diagram; a label without that energy raises ``ValueError``.
    """

    if energy_source not in MLIP_ENERGY_SOURCES:
        raise ValueError(f"energy_source must be one of {MLIP_ENERGY_SOURCES}, got {energy_source!r}")
    if energy_source == "mlip":
        rounds = _relax_rounds(graph)
        last = _task_node(graph, f"mlip_round{rounds[-1]}_relax") if rounds else None
        if last is None or not last.is_finished_ok:
            raise ValueError("the relaxation job of the last round has not finished")
        relaxations = read_mlip_results(graph, labels, round=rounds[-1])
        return {label: relaxations[label].energy_ev for label in labels}
    energies, missing = {}, []
    for label in labels:
        node = _task_node(graph, f"{label}_mlip_dft_energy")
        if node is None or not node.is_finished_ok:
            missing.append(label)
            continue
        energies[label] = float(node.outputs.result.value)
    if missing:
        raise ValueError(f"no DFT energy at the model geometry for {missing} (dft_check off, or not finished)")
    return energies


def _relaxed_structures(graph: Any, labels: Sequence[str]) -> dict[str, Any]:
    """The structures of the last relaxation round of a graph that may still run or have failed (``StructureData``)."""

    rounds = _relax_rounds(graph)
    if not rounds:
        raise ValueError("no relaxation job of this graph has run")
    structures, missing = {}, []
    for label in labels:
        node = _task_node(graph, f"{label}_mlip_round{rounds[-1]}_relaxed")
        if node is None or not node.is_finished_ok:
            missing.append(label)
            continue
        structures[label] = node.outputs.structure
    if missing:
        raise ValueError(f"no relaxed structure of {missing} in round {rounds[-1]}")
    return structures


def _with_molecules(vibrations: VibrationsConfig | None, molecules: Sequence[str]) -> VibrationsConfig | None:
    """``vibrations`` with the gas references among its molecules (their rotations are zero modes)."""

    if vibrations is None:
        return None
    if not isinstance(vibrations, VibrationsConfig):
        raise TypeError("vibrations must be a VibrationsConfig")
    return replace(vibrations, molecules=tuple(dict.fromkeys([*vibrations.molecules, *molecules])))


def _check_half_width(value: Any, slabs: Mapping[str, Any]) -> None:
    """``fixed_half_width_angstrom``: positive, ``None``, or a mapping of slab labels to either."""

    values = dict(value) if isinstance(value, Mapping) else {None: value}
    unknown = sorted(label for label in values if label is not None and label not in slabs)
    if unknown:
        raise ValueError(f"fixed_half_width_angstrom names labels that are not slabs: {unknown}")
    if any(width is not None and not width > 0 for width in values.values()):
        raise ValueError("fixed_half_width_angstrom must be positive (or None)")


def _central_sites(slabs: Mapping[str, Any], half_width: Any) -> dict[str, tuple[int, ...]]:
    """The sites of each slab within its half-width of the mid-plane (a mapping gives one width per slab)."""

    from psteros.structures import central_sites

    widths = dict(half_width) if isinstance(half_width, Mapping) else {label: half_width for label in slabs}
    return {label: central_sites(slab, widths[label]) for label, slab in slabs.items()
            if widths.get(label) is not None}


def _check_frozen_regions(structures: Mapping[str, Any], fixed_sites: Mapping[str, Sequence[int]], bulk: Any) -> None:
    """A frozen slab centre is counted as bulk in the vibrations, so it must have the bulk stoichiometry."""

    from math import isclose

    reference = bulk.composition.get_el_amt_dict()
    wrong = []
    for label, sites in fixed_sites.items():
        if not sites:
            continue
        counts: dict[str, int] = {}
        for index in sites:
            symbol = structures[label][index].specie.symbol
            counts[symbol] = counts.get(symbol, 0) + 1
        ratios = [counts.get(element, 0) / amount for element, amount in reference.items()]
        if set(counts) != set(reference) or any(not isclose(ratio, ratios[0]) for ratio in ratios):
            wrong.append(f"{label} ({''.join(f'{e}{n}' for e, n in sorted(counts.items()))})")
    if wrong:
        raise ValueError(
            f"vibrations: the frozen centre of {', '.join(wrong)} does not have the bulk stoichiometry "
            f"{bulk.composition.reduced_formula}, so it cannot be counted as bulk; change "
            "fixed_half_width_angstrom, for one slab with {label: half-width} (None frees every site)"
        )


def _free_energies(
    energies: Mapping[str, float],
    source: Any,
    temperature_k: float | None,
    vibrations: Mapping[str, Any] | None,
    molecules: Sequence[str],
) -> dict[str, float]:
    """``energies``, or with ``temperature_k`` their free energies with the vibrations of ``source`` or ``vibrations``."""

    if temperature_k is None:
        if vibrations is not None:
            raise ValueError("vibrations are used only with temperature_k")
        return dict(energies)
    if vibrations is None:
        if isinstance(source, Mapping):
            raise ValueError("with energies typed by hand, give vibrations={label: HarmonicVibrations} too")
        from psteros.vibrations_workflow import read_vibrations

        vibrations = read_vibrations(source, list(energies))
    from psteros.vibrations import free_energies_with_vibrations

    return free_energies_with_vibrations(energies, vibrations, temperature_k, bulk_label="bulk",
                                         molecule_labels=tuple(molecules))


def _dft_vibrations(
    relaxed: Mapping[str, Any],
    static: SurfaceWorkflowConfig,
    vibrations: VibrationsConfig | None,
    fixed_sites: Mapping[str, Sequence[int]],
    molecules: Sequence[str],
    submit: bool,
) -> Any:
    """:func:`psteros.build_vibrations_workgraph` (also with CP2K) on the model's geometries, fixed sites kept fixed."""

    from psteros.vibrations_workflow import _vibrations_workgraph

    overrides = dict(static.role_overrides)
    for label, sites in fixed_sites.items():
        if sites:
            overrides[label] = replace(overrides.get(label) or CalculationOverride(), fixed_sites=tuple(sites))
    return _vibrations_workgraph(relaxed, replace(static, role_overrides=overrides),
                                 _with_molecules(vibrations or VibrationsConfig(), molecules), submit=submit,
                                 allow_cp2k=True)


@dataclass
class MlipSurfaceStudy:
    """The MLIP route for your own slabs of a binary compound, from structures to the phase diagram.

    Args:
        bulk: Bulk (pymatgen ``Structure``) relaxed with the DFT recipe of
            the study; it must be a metal (the graph checks its static)
            unless :meth:`build` gets ``allow_non_metallic=True``, because
            the electron count of hand-cut slabs is right only for a metal.
            For semiconductors and insulators :class:`MlipRoute` cuts the
            slabs.
        references: ``{element: structure}`` of the reference of each element,
            relaxed with the same recipe: the elemental solid or a molecule in
            a box (O2). ``mu = E_ref / atoms + Delta mu``.
        slabs: ``{label: slab}`` cut at the cell of ``bulk``, e.g. the
            terminations of :func:`psteros.rutile_110_slab`.
        fixed_half_width_angstrom: Sites within this distance of the slab
            mid-plane stay in place in the MACE relaxations
            (:func:`psteros.central_sites`); ``None`` frees every site; a
            mapping ``{slab label: half-width or None}`` sets it per slab
            (the slabs it leaves out are free).
        rattled: Rattled copies of each slab added as training statics
            (:func:`psteros.rattled_copies`, ``rattle_sigma_angstrom``,
            ``rattle_seed``). 0, with one active-learning round, is the
            default recipe; 5 with ``active_learning_rounds=0`` replaces the
            round (more DFT up front, one fine-tune).

    The calculations are labelled ``bulk``, ``ref_<element>`` and the slab
    labels. :meth:`build` returns the graph of
    :func:`build_mlip_surface_workgraph`; :meth:`analyse` the phase diagram
    from its energies, all of one source, and with ``temperature_k`` from the
    free energies of the vibrations.
    """

    bulk: Any
    references: Mapping[str, Any]
    slabs: Mapping[str, Any]
    fixed_half_width_angstrom: float | Mapping[str, float | None] | None = 1.5
    rattled: int = 0
    rattle_sigma_angstrom: float = 0.05
    rattle_seed: int = 20261017
    structures: dict = field(init=False)

    def __post_init__(self) -> None:
        from psteros.workflow import _check_labels

        self.elements = sorted({site.specie.symbol for site in self.bulk})
        if len(self.elements) != 2:
            raise ValueError(
                f"MlipSurfaceStudy is for binary compounds, got {self.bulk.composition.reduced_formula}"
            )
        missing = sorted(set(self.elements).difference(self.references))
        if missing:
            raise ValueError(f"references lacks a reference structure for {missing}")
        if not self.slabs:
            raise ValueError("at least one slab is required")
        if isinstance(self.rattled, bool) or not isinstance(self.rattled, int) or self.rattled < 0:
            raise ValueError("rattled must be a non-negative integer")
        _check_half_width(self.fixed_half_width_angstrom, self.slabs)
        self.formula = self.bulk.composition.reduced_formula
        self.structures = {"bulk": self.bulk, **{f"ref_{e}": self.references[e] for e in self.elements}}
        for label, slab in self.slabs.items():
            if label in self.structures:
                raise ValueError(f"slab label {label!r} is taken by a reference")
            if {site.specie.symbol for site in slab} - set(self.elements):
                raise ValueError(f"{label}: the slab has elements that are not in the bulk")
            self.structures[label] = slab
        _check_labels(self.structures)

    @property
    def fixed_sites(self) -> dict[str, tuple[int, ...]]:
        """Sites of each slab that the MACE relaxations keep in place."""

        if self.fixed_half_width_angstrom is None:
            return {}
        return _central_sites(self.slabs, self.fixed_half_width_angstrom)

    @property
    def training_structures(self) -> dict[str, Any]:
        """The rattled slabs (only DFT statics for training), ``{"<slab>_r<i>": structure}``."""

        if not self.rattled:
            return {}
        from psteros.mlip import rattled_copies

        return rattled_copies(self.slabs, count=self.rattled, sigma_angstrom=self.rattle_sigma_angstrom,
                              seed=self.rattle_seed)

    @property
    def molecule_labels(self) -> tuple[str, ...]:
        """Labels of the gas references (molecules in a box)."""

        from psteros.surface_study import GAS_REFERENCES

        return tuple(f"ref_{element}" for element in self.elements if element in GAS_REFERENCES)

    def calculation_count(
        self, active_learning_rounds: int = 1, dft_check: bool = True, vibrations: bool = False
    ) -> dict[str, int]:
        """DFT statics and GPU jobs of :meth:`build` with these options."""

        structures = len(self.structures)
        statics = structures + len(self.slabs) * self.rattled
        statics += structures * active_learning_rounds + (structures if dft_check else 0)
        return {"dft_statics": statics, "gpu_jobs": 2 * (active_learning_rounds + 1) + (1 if vibrations else 0)}

    def summary(self) -> str:
        counts = self.calculation_count()
        lines = [f"{self.formula}: {len(self.slabs)} slabs, references {', '.join(self.elements)}; "
                 f"{counts['dft_statics']} DFT statics and {counts['gpu_jobs']} GPU jobs with the defaults"]
        fixed = self.fixed_sites
        for label, slab in self.slabs.items():
            lines.append(f"  {label}: {slab.composition.formula.replace(' ', '')}, {len(slab)} atoms, "
                         f"{len(fixed.get(label, ()))} fixed")
        return "\n".join(lines)

    __str__ = summary

    def build(
        self,
        static: SurfaceWorkflowConfig,
        fine_tune: MaceFineTuneConfig | UmaFineTuneConfig | NequipFineTuneConfig,
        relax: MlipRelaxConfig = MlipRelaxConfig(),
        *,
        job: MlipJobConfig,
        active_learning_rounds: int = 1,
        dft_check: bool = True,
        vibrations: VibrationsConfig | None = None,
        allow_non_metallic: bool = False,
        submit: bool = False,
    ) -> Any:
        """The graph of :func:`build_mlip_surface_workgraph` for these structures (bulk checked as a metal).

        ``vibrations`` adds the model's harmonic modes of every structure
        (the gas references are added to its ``molecules``); the frozen
        centre of each slab must then have the bulk stoichiometry.
        ``allow_non_metallic`` skips the metal check (see
        :func:`build_mlip_surface_workgraph`).
        """

        if vibrations is not None:
            _check_frozen_regions(self.slabs, self.fixed_sites, self.bulk)
        return build_mlip_surface_workgraph(
            self.structures, static, fine_tune, relax, job=job, metallic_bulk="bulk", fixed_sites=self.fixed_sites,
            training_structures=self.training_structures or None, active_learning_rounds=active_learning_rounds,
            dft_check=dft_check, vibrations=_with_molecules(vibrations, self.molecule_labels),
            allow_non_metallic=allow_non_metallic, submit=submit,
        )

    def dft_vibrations(
        self,
        source: Any,
        static: SurfaceWorkflowConfig,
        vibrations: VibrationsConfig | None = None,
        *,
        submit: bool = False,
    ) -> Any:
        """DFT finite differences at the model's geometries, to check the model's vibrations.

        ``source`` is the finished graph of :meth:`build`; the graph of
        :func:`psteros.build_vibrations_workgraph` is built on its relaxed
        structures with ``static`` (QE, VASP, or CP2K with ``RUN_TYPE
        ENERGY_FORCE``), the slab centres kept fixed as in the relaxations
        and the gas references as molecules. Its ``<label>_vibrations`` go to
        :meth:`analyse` as ``vibrations=psteros.read_vibrations(...)``.
        """

        relaxed = _relaxed_structures(source, list(self.structures))
        return _dft_vibrations(relaxed, static, vibrations, self.fixed_sites, self.molecule_labels, submit)

    def references_from(self, energies_ev: Mapping[str, float]):
        """:class:`psteros.BinaryReferences` from the energies of ``bulk`` and ``ref_<element>``."""

        from psteros.phase_diagram import BinaryReferences
        from psteros.surface_study import reference_energies_per_atom, reservoir_labels

        per_atom = reference_energies_per_atom(
            self.elements, self.references, {f"ref_{e}": energies_ev[f"ref_{e}"] for e in self.elements})
        return BinaryReferences(
            bulk_energy_ev=energies_ev["bulk"], bulk_composition=self.bulk.composition,
            reference_energies_per_atom_ev=per_atom, reservoir_labels=reservoir_labels(self.elements),
        )

    def analyse(
        self,
        source: Any,
        *,
        energy_source: str = "mlip_dft",
        temperature_k: float | None = None,
        vibrations: Mapping[str, Any] | None = None,
        **diagram_options,
    ) -> "MlipStudyResult":
        """Surface phase diagram of the slabs.

        ``source`` is the finished graph (``orm.load_node(pk)``), whose
        energies are read with :func:`read_mlip_energies` from
        ``energy_source`` (``"mlip_dft"``, the default, or ``"mlip"``), or a
        mapping ``{label: energy in eV}`` typed by hand, which must come from
        one source too. Composition and area are those of the slabs (the
        relaxations keep the cell). ``diagram_options`` go to
        :func:`psteros.surface_phase_diagram`.

        With ``temperature_k`` every energy becomes a free energy
        (:func:`psteros.free_energies_with_vibrations`): solids
        ``E + F_vib(T)``, the slab centres counted as bulk, gas references
        ``E + ZPE``. The modes are those of the graph (built with
        ``vibrations``) or ``vibrations={label: HarmonicVibrations}``, for
        example ``psteros.read_vibrations`` of :meth:`dft_vibrations`.
        """

        from psteros.phase_diagram import SlabTermination, surface_phase_diagram

        if isinstance(source, Mapping):
            energies = {label: float(value) for label, value in source.items()}
            missing = sorted(set(self.structures).difference(energies))
            if missing:
                raise ValueError(f"no energy for {missing}")
        else:
            energies = read_mlip_energies(source, list(self.structures), energy_source=energy_source)
        energies = _free_energies({label: energies[label] for label in self.structures}, source, temperature_k,
                                  vibrations, self.molecule_labels)
        references = self.references_from(energies)
        terminations = tuple(
            SlabTermination.from_structure(label, energies[label], slab) for label, slab in self.slabs.items())
        diagram = surface_phase_diagram(terminations, references, **diagram_options)
        return MlipStudyResult(self.formula, references, diagram, terminations, energy_source=energy_source,
                               temperature_k=temperature_k)


@dataclass(frozen=True)
class MlipStudyResult(ChargeNeutralStudyResult):
    """References, phase diagram and terminations of an :class:`MlipSurfaceStudy`, and the energy source used.

    ``temperature_k`` is the temperature of the vibrational free energies, or
    ``None`` for total energies.
    """

    energy_source: str = "mlip_dft"
    temperature_k: float | None = None


@dataclass
class MlipRoute:
    """The MLIP route for the slabs that psteros cuts with the right electron count.

    Args:
        study: A :class:`psteros.ChargeNeutralSurfaceStudy` (charge-neutral
            symmetric slabs of an element, a binary or a ternary compound,
            one or more orientations) or a :class:`psteros.PolarSurfaceStudy`
            (polar faces on a pseudo-hydrogen passivated bottom, VASP only).
            Its bulk and references must already be relaxed with the DFT
            recipe of the route; its slabs are not relaxed with DFT.
        fixed_half_width_angstrom: Sites within this distance of the
            mid-plane of each charge-neutral slab stay in place in the MACE
            relaxations (:func:`psteros.central_sites`); ``None`` frees every
            site; a mapping ``{slab label: half-width or None}`` sets it per
            slab. Polar slabs relax every atom, as in the polar method, and
            their bottoms are checked by :meth:`analyse`.
        rattled, rattle_sigma_angstrom, rattle_seed: Rattled copies of each
            slab as extra training statics, as in :class:`MlipSurfaceStudy`.

    The slabs are charge-neutral or passivated by construction, so the bulk
    is not checked to be a metal: the route serves metals, semiconductors and
    insulators. :meth:`build` returns the graph of
    :func:`build_mlip_surface_workgraph`; :meth:`analyse` gives the study's
    own result (phase diagram, bottom checks, ...) from the graph's energies
    and relaxed structures.
    """

    study: Any
    fixed_half_width_angstrom: float | Mapping[str, float | None] | None = 1.5
    rattled: int = 0
    rattle_sigma_angstrom: float = 0.05
    rattle_seed: int = 20261017

    def __post_init__(self) -> None:
        from psteros.polar_workflow import PolarSurfaceStudy
        from psteros.surface_study import ChargeNeutralSurfaceStudy

        if not isinstance(self.study, (ChargeNeutralSurfaceStudy, PolarSurfaceStudy)):
            raise TypeError(
                "MlipRoute takes a ChargeNeutralSurfaceStudy or a PolarSurfaceStudy, whose slabs have the right "
                "electron count; for your own slabs use MlipSurfaceStudy"
            )
        if isinstance(self.rattled, bool) or not isinstance(self.rattled, int) or self.rattled < 0:
            raise ValueError("rattled must be a non-negative integer")
        _check_half_width(self.fixed_half_width_angstrom, self.slabs)
        self.polar = isinstance(self.study, PolarSurfaceStudy)
        self.formula = self.study.formula

    @property
    def structures(self) -> dict[str, Any]:
        """Every calculation of the study, by label."""

        return dict(self.study.structures)

    @property
    def slabs(self) -> dict[str, Any]:
        """The slabs of the diagram (role ``slab``), by label."""

        return {label: self.study.structures[label] for label, role in self.study.roles.items() if role == "slab"}

    @property
    def fixed_sites(self) -> dict[str, tuple[int, ...]]:
        """Sites of each charge-neutral slab that the MACE relaxations keep in place (none for polar slabs)."""

        if self.polar or self.fixed_half_width_angstrom is None:
            return {}
        return _central_sites(self.slabs, self.fixed_half_width_angstrom)

    @property
    def training_structures(self) -> dict[str, Any]:
        """The rattled slabs (only DFT statics for training), ``{"<slab>_r<i>": structure}``."""

        if not self.rattled:
            return {}
        from psteros.mlip import rattled_copies

        return rattled_copies(self.slabs, count=self.rattled, sigma_angstrom=self.rattle_sigma_angstrom,
                              seed=self.rattle_seed)

    @property
    def molecule_labels(self) -> tuple[str, ...]:
        """Labels of the gas references and, for a polar study, of the pseudo-molecules and clusters."""

        from psteros.surface_study import GAS_REFERENCES

        labels = []
        for label, role in self.study.roles.items():
            if role == "reference" and label[len("ref_"):] in GAS_REFERENCES:
                labels.append(label)
            elif role in ("pseudo_molecule", "cluster"):
                labels.append(label)
        return tuple(labels)

    def calculation_count(
        self, active_learning_rounds: int = 1, dft_check: bool = True, vibrations: bool = False
    ) -> dict[str, int]:
        """DFT statics and GPU jobs of :meth:`build` with these options."""

        structures = len(self.study.structures)
        statics = structures + len(self.slabs) * self.rattled
        statics += structures * active_learning_rounds + (structures if dft_check else 0)
        return {"dft_statics": statics, "gpu_jobs": 2 * (active_learning_rounds + 1) + (1 if vibrations else 0)}

    def summary(self) -> str:
        counts = self.calculation_count()
        kind = "polar, passivated bottom" if self.polar else "charge-neutral"
        lines = [f"{self.formula}: {len(self.slabs)} slabs ({kind}), {len(self.study.structures)} calculations; "
                 f"{counts['dft_statics']} DFT statics and {counts['gpu_jobs']} GPU jobs with the defaults"]
        fixed = self.fixed_sites
        for label, slab in self.slabs.items():
            lines.append(f"  {label}: {slab.composition.formula.replace(' ', '')}, {len(slab)} atoms, "
                         f"{len(fixed.get(label, ()))} fixed")
        return "\n".join(lines)

    __str__ = summary

    def build(
        self,
        static: SurfaceWorkflowConfig,
        fine_tune: MaceFineTuneConfig | UmaFineTuneConfig | NequipFineTuneConfig,
        relax: MlipRelaxConfig = MlipRelaxConfig(),
        *,
        job: MlipJobConfig,
        active_learning_rounds: int = 1,
        dft_check: bool = True,
        vibrations: VibrationsConfig | None = None,
        submit: bool = False,
    ) -> Any:
        """The graph of :func:`build_mlip_surface_workgraph` for the study (no metal check).

        ``static`` is the DFT static recipe of every calculation, with the
        study's overrides of the static stage (``study.vasp_overrides("static")``,
        ``study.qe_overrides("static")``, or your own for CP2K); the rattled
        copies of a slab use the slab's override. A polar study needs VASP
        (the pseudo-hydrogen potentials). ``vibrations`` adds the model's
        harmonic modes of every structure (gas references are added to its
        ``molecules``; the frozen centre of each slab must have the bulk
        stoichiometry); polar studies have none, because the pseudo chemical
        potentials are 0 K energies.
        """

        if self.polar and static.backend != "vasp":
            raise ValueError(
                f"{self.formula}: a polar study runs with VASP (pseudo-hydrogen potentials), got backend {static.backend!r}"
            )
        if self.polar and vibrations is not None:
            raise ValueError(
                f"{self.formula}: polar slabs have no vibrations yet (their pseudo chemical potentials are 0 K energies)"
            )
        extra = self.training_structures
        if extra:
            overrides = dict(static.role_overrides)
            for label in extra:
                slab = label.rsplit("_r", 1)[0]
                if slab in overrides and label not in overrides:
                    overrides[label] = overrides[slab]
            static = replace(static, role_overrides=overrides)
        if vibrations is not None:
            _check_frozen_regions(self.slabs, self.fixed_sites, self.study.bulk)
        return _mlip_workgraph(
            self.structures, static, fine_tune, relax, job=job, metallic_bulk="bulk", metal_check=False,
            fixed_sites=self.fixed_sites, training_structures=extra or None,
            active_learning_rounds=active_learning_rounds, dft_check=dft_check,
            vibrations=_with_molecules(vibrations, self.molecule_labels), submit=submit,
        )

    def read(self, source: Any, *, energy_source: str = "mlip_dft") -> tuple[dict[str, float], dict[str, Any]]:
        """Energies (eV, all from ``energy_source``) and relaxed structures of every calculation of a graph.

        The same pair as :func:`psteros.read_vasp_results`, so it can go to
        the study's own ``analyse``.
        """

        labels = list(self.study.structures)
        return read_mlip_energies(source, labels, energy_source=energy_source), _relaxed_structures(source, labels)

    def analyse(
        self,
        source: Any,
        *,
        energy_source: str = "mlip_dft",
        relaxed_structures: Mapping[str, Any] | None = None,
        temperature_k: float | None = None,
        vibrations: Mapping[str, Any] | None = None,
        **options,
    ) -> Any:
        """The study's result from the graph of :meth:`build`.

        ``source`` is the finished graph (energies from ``energy_source``,
        relaxed structures of the last round) or energies typed by hand
        ``{label: eV}``, with ``relaxed_structures`` if you have them (the
        built structures otherwise; for a polar study the bottom check then
        compares nothing). ``temperature_k`` and ``vibrations`` turn the
        energies into free energies as in :meth:`MlipSurfaceStudy.analyse`.
        ``options`` go to the study's ``analyse``. A charge-neutral study
        gives an :class:`MlipStudyResult`, a polar one a
        :class:`psteros.PolarStudyResult` (with the bottom checks of the
        MACE-relaxed slabs).
        """

        if self.polar and temperature_k is not None:
            raise ValueError(f"{self.formula}: polar slabs have no vibrations yet; leave temperature_k out")
        if isinstance(source, Mapping):
            energies = {label: float(value) for label, value in source.items()}
            missing = sorted(set(self.study.structures).difference(energies))
            if missing:
                raise ValueError(f"no energy for {missing}")
            relaxed = dict(relaxed_structures or self.study.structures)
        else:
            energies, relaxed = self.read(source, energy_source=energy_source)
            relaxed.update(relaxed_structures or {})
        energies = _free_energies({label: energies[label] for label in self.study.structures}, source,
                                  temperature_k, vibrations, self.molecule_labels)
        if self.polar:
            return self.study.analyse(energies, relaxed, **options)
        result = self.study.analyse(energies, relaxed, **options)
        return MlipStudyResult(result.formula, result.references, result.diagram, result.terminations,
                               result.surface_energies_j_per_m2, energy_source=energy_source,
                               temperature_k=temperature_k)

    def dft_vibrations(
        self,
        source: Any,
        static: SurfaceWorkflowConfig,
        vibrations: VibrationsConfig | None = None,
        *,
        submit: bool = False,
    ) -> Any:
        """DFT finite differences at the model's geometries (see :meth:`MlipSurfaceStudy.dft_vibrations`)."""

        if self.polar:
            raise ValueError(f"{self.formula}: polar slabs have no vibrations yet")
        relaxed = _relaxed_structures(source, list(self.study.structures))
        return _dft_vibrations(relaxed, static, vibrations, self.fixed_sites, self.molecule_labels, submit)
