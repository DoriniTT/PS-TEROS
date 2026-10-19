"""WorkGraph of the MLIP route to a surface phase diagram (metallic systems).

:func:`build_mlip_surface_workgraph` runs the route (``docs/MLIP_SURFACES.md``):

1. one DFT static calculation (energy and forces) of every structure, plus
   optional extra training structures such as :func:`psteros.rattled_copies`;
2. a replay fine-tune of a MACE foundation model on those labels (one model
   per seed of :class:`~psteros.MaceFineTuneConfig`), on a GPU computer;
3. the model of the first seed relaxes every structure (cell fixed);
4. ``active_learning_rounds`` times: one DFT static at each relaxed geometry
   is added to the labels, the model is fine-tuned again from the foundation
   model on all labels and relaxes the structures again from the start;
5. with ``dft_check``, one DFT static at each final geometry.

The energies of one diagram come from one source (see
:func:`psteros.check_one_energy_source`): ``"mlip_dft"``, the DFT statics of
step 5 (the default), or ``"mlip"``, the model at its own geometries.
:func:`read_mlip_energies` returns one of
them for every label; :func:`read_mlip_results` the relaxations.

:class:`MlipSurfaceStudy` is the route for the slabs of a metallic binary
compound from structures to the phase diagram: it names the calculations,
fixes the slab centres, builds the graph and turns its energies into a
:class:`psteros.SurfacePhaseDiagram`.
"""

from __future__ import annotations

import json
import tempfile
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride, ExecutionPolicy, SurfaceWorkflowConfig
from psteros.mlip import MaceFineTuneConfig, MlipRelaxation, MlipRelaxConfig
from psteros.surface_study import ChargeNeutralStudyResult

MLIP_ENERGY_SOURCES = ("mlip", "mlip_dft")


@dataclass(frozen=True)
class MlipJobConfig:
    """Where the MACE fine-tunes and relaxations run.

    ``code_label`` is an AiiDA code whose executable is the Python of an
    environment with ase, torch and mace-torch on the GPU computer that holds
    ``MaceFineTuneConfig.foundation_model``; ``execution`` gives the scheduler
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
    fine_tune: MaceFineTuneConfig,
    relax: MlipRelaxConfig = MlipRelaxConfig(),
    *,
    job: MlipJobConfig,
    metallic_bulk: str,
    fixed_sites: Mapping[str, Sequence[int]] | None = None,
    training_structures: Mapping[str, Any] | None = None,
    active_learning_rounds: int = 1,
    dft_check: bool = True,
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
        static); the first fine-tune stops otherwise. The route is for
        metallic systems only.
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
    submit
        Submit the graph when true; otherwise return it unsubmitted.

    Graph outputs (all new; no other builder changes): per label
    ``<label>_static_parameters`` and ``<label>_static_retrieved`` (step 1),
    ``<label>_mlip_relaxed_structure`` and ``<label>_mlip_energy`` (eV, the
    first seed's model of the last round), with ``dft_check``
    ``<label>_mlip_dft_energy`` (eV), ``<label>_mlip_dft_parameters`` and
    ``<label>_mlip_dft_retrieved``; per round ``r`` (0 is the first fine-tune)
    ``mlip_round<r>_finetune_summary``, ``mlip_round<r>_model_folder`` and
    ``mlip_round<r>_relax_results`` (every seed's relaxations).
    """

    from psteros.mlip import select_replay_structures
    from psteros.workflow import _check_labels

    _check_labels(structures)
    extra = dict(training_structures or {})
    if extra:
        _check_labels(extra)
    if not isinstance(job, MlipJobConfig):
        raise TypeError("job must be an MlipJobConfig")
    if not isinstance(fine_tune, MaceFineTuneConfig) or not isinstance(relax, MlipRelaxConfig):
        raise TypeError("fine_tune must be a MaceFineTuneConfig and relax an MlipRelaxConfig")
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
    for label, sites in fixed_sites.items():
        if any(site < 0 or site >= len(plain[label]) for site in sites):
            raise ValueError(f"{label}: fixed_sites must be within [0, {len(plain[label])})")

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

        elements = sorted({site.specie.symbol for item in (*plain.values(), *plain_extra.values()) for site in item})
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
    for round_ in range(active_learning_rounds + 1):
        name = f"{static.name}_round{round_}"
        settings = orm.Dict({
            "backend": static.backend, "name": name, "device": job.device, "metallic_bulk": metallic_bulk,
            "frames": list(labels_of_frames), "fine_tune": json.loads(json.dumps(asdict(fine_tune))),
            "replay_selection": replay_record,
        })
        folder = workgraph.add_task(
            mlip_training_folder, name=f"mlip_round{round_}_training", settings=settings,
            inputs={**frames, **replay_nodes},
        )
        finetune = workgraph.add_task(
            shelljob, name=f"mlip_round{round_}_finetune", command=code, arguments=["-u", "finetune.py"],
            nodes={"job": folder.outputs.result}, outputs=["summary.json"],
            metadata={"label": f"{name}_finetune", "options": dict(execution)},
        )
        models = {f"seed{seed}": (f"finetune/models/{name}_seed{seed}.model", relax.head) for seed in fine_tune.seeds}
        with tempfile.TemporaryDirectory() as scratch:
            relax_folder = orm.FolderData(tree=str(prepare_relax_job(
                Path(scratch) / "relax", plain, models, relax, fixed_sites=fixed_sites, device=job.device)))
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
    if submit:
        workgraph.submit()
    return workgraph


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


@dataclass
class MlipSurfaceStudy:
    """The MLIP route for the slabs of a metallic binary compound, from structures to the phase diagram.

    Args:
        bulk: Bulk (pymatgen ``Structure``) relaxed with the DFT recipe of
            the study; it must be a metal (the graph checks its static).
        references: ``{element: structure}`` of the reference of each element,
            relaxed with the same recipe: the elemental solid or a molecule in
            a box (O2). ``mu = E_ref / atoms + Delta mu``.
        slabs: ``{label: slab}`` cut at the cell of ``bulk``, e.g. the
            terminations of :func:`psteros.rutile_110_slab`.
        fixed_half_width_angstrom: Sites within this distance of the slab
            mid-plane stay in place in the MACE relaxations
            (:func:`psteros.central_sites`); ``None`` frees every site.
        rattled: Rattled copies of each slab added as training statics
            (:func:`psteros.rattled_copies`, ``rattle_sigma_angstrom``,
            ``rattle_seed``). 0, with one active-learning round, is the
            default recipe; 5 with ``active_learning_rounds=0`` replaces the
            round (more DFT up front, one fine-tune).

    The calculations are labelled ``bulk``, ``ref_<element>`` and the slab
    labels. :meth:`build` returns the graph of
    :func:`build_mlip_surface_workgraph`; :meth:`analyse` the phase diagram
    from its energies, all of one source.
    """

    bulk: Any
    references: Mapping[str, Any]
    slabs: Mapping[str, Any]
    fixed_half_width_angstrom: float | None = 1.5
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
        if self.fixed_half_width_angstrom is not None and self.fixed_half_width_angstrom <= 0:
            raise ValueError("fixed_half_width_angstrom must be positive (or None)")
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
        from psteros.structures import central_sites

        return {label: central_sites(slab, self.fixed_half_width_angstrom) for label, slab in self.slabs.items()}

    @property
    def training_structures(self) -> dict[str, Any]:
        """The rattled slabs (only DFT statics for training), ``{"<slab>_r<i>": structure}``."""

        if not self.rattled:
            return {}
        from psteros.mlip import rattled_copies

        return rattled_copies(self.slabs, count=self.rattled, sigma_angstrom=self.rattle_sigma_angstrom,
                              seed=self.rattle_seed)

    def calculation_count(self, active_learning_rounds: int = 1, dft_check: bool = True) -> dict[str, int]:
        """DFT statics and GPU jobs of :meth:`build` with these options."""

        structures = len(self.structures)
        statics = structures + len(self.slabs) * self.rattled
        statics += structures * active_learning_rounds + (structures if dft_check else 0)
        return {"dft_statics": statics, "gpu_jobs": 2 * (active_learning_rounds + 1)}

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
        fine_tune: MaceFineTuneConfig,
        relax: MlipRelaxConfig = MlipRelaxConfig(),
        *,
        job: MlipJobConfig,
        active_learning_rounds: int = 1,
        dft_check: bool = True,
        submit: bool = False,
    ) -> Any:
        """The graph of :func:`build_mlip_surface_workgraph` for these structures (bulk checked as a metal)."""

        return build_mlip_surface_workgraph(
            self.structures, static, fine_tune, relax, job=job, metallic_bulk="bulk", fixed_sites=self.fixed_sites,
            training_structures=self.training_structures or None, active_learning_rounds=active_learning_rounds,
            dft_check=dft_check, submit=submit,
        )

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

    def analyse(self, source: Any, *, energy_source: str = "mlip_dft", **diagram_options) -> "MlipStudyResult":
        """Surface phase diagram of the slabs.

        ``source`` is the finished graph (``orm.load_node(pk)``), whose
        energies are read with :func:`read_mlip_energies` from
        ``energy_source`` (``"mlip_dft"``, the default, or ``"mlip"``), or a
        mapping ``{label: energy in eV}`` typed by hand, which must come from
        one source too. Composition and area are those of the slabs (the
        relaxations keep the cell). ``diagram_options`` go to
        :func:`psteros.surface_phase_diagram`.
        """

        from psteros.phase_diagram import SlabTermination, surface_phase_diagram

        if isinstance(source, Mapping):
            energies = {label: float(value) for label, value in source.items()}
            missing = sorted(set(self.structures).difference(energies))
            if missing:
                raise ValueError(f"no energy for {missing}")
        else:
            energies = read_mlip_energies(source, list(self.structures), energy_source=energy_source)
        references = self.references_from(energies)
        terminations = tuple(
            SlabTermination.from_structure(label, energies[label], slab) for label, slab in self.slabs.items())
        diagram = surface_phase_diagram(terminations, references, **diagram_options)
        return MlipStudyResult(self.formula, references, diagram, terminations, energy_source=energy_source)


@dataclass(frozen=True)
class MlipStudyResult(ChargeNeutralStudyResult):
    """References, phase diagram and terminations of an :class:`MlipSurfaceStudy`, and the energy source used."""

    energy_source: str = "mlip_dft"
