"""Machine-learned potentials for surface energies of metallic systems (no AiiDA).

The route: DFT static calculations of the bulk, the references and every
termination label a fine-tune of a MACE foundation model; the fine-tuned
model relaxes the terminations; the surface phase diagram uses its energies
(or one DFT static on each relaxed geometry).

This module holds the parts that need neither AiiDA nor a GPU: the labelled
structures and their extended-XYZ file, the per-element energy shifts of the
fine-tuned head, the MACE training options, a relaxation with any ASE
calculator, the metallic-bulk check and the one-energy-source rule.

Why one energy source: with γ = [E_slab − Σ_i N_i μ_i] / (n A) and
μ_i = E_i(reference) + Δμ_i, adding a constant δ_Z per atom of element Z to
every energy changes γ by Σ_Z N_Z δ_Z − Σ_Z N_Z δ_Z = 0 only if the slab and
all the references carry the same δ_Z, that is, come from the same model or
the same DFT. Mixing sources brings the offsets straight into γ.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

import numpy as np

TRAINING_ROLES = ("train", "valid", "test")
OPTIMIZERS = ("FIRE", "BFGS", "LBFGS")
ENERGY_SOURCES = ("dft", "mlip", "mlip_dft")


@dataclass(frozen=True)
class LabelledStructure:
    """One DFT static calculation as a training example.

    ``structure`` is a pymatgen ``Structure``; ``energy_ev`` is the energy
    whose derivatives are ``forces_ev_per_angstrom`` (with smearing, the
    electronic free energy); one force row per site, in eV/A.
    """

    label: str
    structure: Any
    energy_ev: float
    forces_ev_per_angstrom: Sequence[Sequence[float]]
    role: str = "train"

    def __post_init__(self) -> None:
        if not self.label or any(character.isspace() for character in self.label):
            raise ValueError(f"label must be a non-empty word, got {self.label!r}")
        if self.role not in TRAINING_ROLES:
            raise ValueError(f"{self.label}: role must be one of {TRAINING_ROLES}, got {self.role!r}")
        forces = np.asarray(self.forces_ev_per_angstrom, dtype=float)
        if forces.shape != (len(self.structure), 3):
            raise ValueError(
                f"{self.label}: forces must have one row of 3 per site "
                f"({len(self.structure)} sites), got shape {forces.shape}"
            )
        if not math.isfinite(self.energy_ev) or not np.isfinite(forces).all():
            raise ValueError(f"{self.label}: energy and forces must be finite")
        object.__setattr__(self, "forces_ev_per_angstrom", tuple(tuple(row) for row in forces.tolist()))

    @property
    def composition(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for site in self.structure:
            counts[site.specie.symbol] = counts.get(site.specie.symbol, 0) + 1
        return counts


def write_extxyz(frames: Sequence[LabelledStructure], destination: str | Path | TextIO) -> None:
    """Write frames as extended XYZ with ``REF_energy`` (eV) and ``REF_forces`` (eV/A), as MACE reads them.

    Each frame also records ``label`` and ``role``; every frame is periodic.
    """

    if not frames:
        raise ValueError("no frames to write")
    labels = [frame.label for frame in frames]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"duplicate frame labels: {duplicates}")
    lines = []
    for frame in frames:
        lattice = " ".join(f"{value:.10f}" for value in np.asarray(frame.structure.lattice.matrix).ravel())
        lines.append(str(len(frame.structure)))
        lines.append(
            f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:REF_forces:R:3 '
            f"REF_energy={frame.energy_ev:.10f} label={frame.label} role={frame.role} "
            'pbc="T T T"'
        )
        for site, force in zip(frame.structure, frame.forces_ev_per_angstrom):
            position = " ".join(f"{value:.10f}" for value in site.coords)
            values = " ".join(f"{value:.10f}" for value in force)
            lines.append(f"{site.specie.symbol} {position} {values}")
    text = "\n".join(lines) + "\n"
    if hasattr(destination, "write"):
        destination.write(text)
    else:
        Path(destination).write_text(text)


def fit_element_shifts(
    compositions: Sequence[Mapping[str, int]], residuals_ev: Sequence[float]
) -> dict[str, float]:
    """Per-element energy shifts δ_Z (eV) that best explain the residuals in least squares.

    ``residuals_ev[i]`` is E_DFT − E_base of frame ``i`` (DFT label minus the
    foundation model). Added to the foundation model's E0s, the shifts make the
    fine-tuned head start close to the DFT energies. Several compositions are
    needed to separate the elements (e.g. the bulk oxide, the metal and O2);
    otherwise a ``ValueError`` names the elements. By the one-source rule these
    shifts do not change surface energies; they only help the fit.
    """

    if len(compositions) != len(residuals_ev) or not compositions:
        raise ValueError("give one residual per composition, at least one")
    elements = sorted({element for composition in compositions for element in composition})
    matrix = np.array([[composition.get(element, 0) for element in elements] for composition in compositions], float)
    if np.linalg.matrix_rank(matrix) < len(elements):
        raise ValueError(
            f"the compositions do not separate the elements {elements}; add a reference that changes their ratio "
            "(the metal, the gas molecule)"
        )
    shifts, *_ = np.linalg.lstsq(matrix, np.asarray(residuals_ev, float), rcond=None)
    return {element: float(value) for element, value in zip(elements, shifts)}


@dataclass(frozen=True)
class MaceFineTuneConfig:
    """Recipe of a MACE fine-tune (by default the Ag3PO4 r4/goal-07 recipe with a replay head).

    ``foundation_model`` is the path of the model file on the GPU computer
    (checked against ``foundation_model_sha256`` when given) and
    ``foundation_head`` its head (MACE-MH-1: ``"omat_pbe"``). With ``replay``,
    the fine-tune is multi-head: :func:`select_replay_structures` picks from the
    local ``replay_file`` (the foundation model's replay data, e.g.
    ``replay-data-mh-1-omat-pbe.xyz``) every distinct structure made only of
    the elements of the training set, with ``replay_min_atoms`` to
    ``replay_max_atoms`` atoms; the job labels them with the foundation model
    and trains the ``pt_head`` on them; every ``replay_validation_every``-th
    one is kept for validation. The fine-tuned head is ``"Default"``.

    The training runs ``epochs`` passes at the constant ``learning_rate`` with
    an EMA of ``ema_decay``, the universal loss (energy and force weights,
    Huber ``huber_delta``) and batch size 1; the final EMA model is kept, the
    validation is reported but never picks the model. One model per seed.
    """

    foundation_model: str
    foundation_model_sha256: str | None = None
    foundation_head: str = "omat_pbe"
    replay: bool = True
    replay_file: str | None = None
    replay_min_atoms: int = 2
    replay_max_atoms: int = 128
    replay_validation_every: int = 5
    epochs: int = 100
    learning_rate: float = 1.0e-3
    ema_decay: float = 0.999
    energy_weight: float = 10.0
    forces_weight: float = 10.0
    huber_delta: float = 0.01
    clip_grad: float = 1.0
    seeds: tuple[int, ...] = (20260920, 20260921, 20260922)
    default_dtype: str = "float64"

    def __post_init__(self) -> None:
        if not self.foundation_model or not self.foundation_head:
            raise ValueError("MaceFineTuneConfig needs foundation_model and foundation_head")
        if self.replay and not self.replay_file:
            raise ValueError("MaceFineTuneConfig: replay=True needs replay_file")
        if not 1 <= self.replay_min_atoms <= self.replay_max_atoms:
            raise ValueError("MaceFineTuneConfig: need 1 <= replay_min_atoms <= replay_max_atoms")
        if self.replay_validation_every < 2:
            raise ValueError("MaceFineTuneConfig: replay_validation_every must be at least 2")
        if isinstance(self.epochs, bool) or not isinstance(self.epochs, int) or self.epochs < 1:
            raise ValueError("MaceFineTuneConfig: epochs must be a positive integer")
        for name in ("learning_rate", "energy_weight", "forces_weight", "huber_delta", "clip_grad"):
            if getattr(self, name) <= 0:
                raise ValueError(f"MaceFineTuneConfig: {name} must be positive")
        if not 0 < self.ema_decay < 1:
            raise ValueError("MaceFineTuneConfig: ema_decay must be between 0 and 1")
        seeds = tuple(self.seeds)
        if not seeds or len(set(seeds)) != len(seeds) or any(isinstance(s, bool) or not isinstance(s, int) for s in seeds):
            raise ValueError("MaceFineTuneConfig: seeds must be distinct integers, at least one")
        object.__setattr__(self, "seeds", seeds)
        if self.default_dtype not in ("float32", "float64"):
            raise ValueError("MaceFineTuneConfig: default_dtype must be 'float32' or 'float64'")


def select_replay_structures(
    path: str | Path,
    *,
    elements: Sequence[str],
    min_atoms: int = 2,
    max_atoms: int = 128,
    validation_every: int = 5,
) -> tuple[list[Any], list[Any], dict[str, Any]]:
    """Replay geometries for a multi-head fine-tune (the goal-07 rule of the Ag3PO4 project).

    Reads the extended-XYZ ``path`` frame by frame, in file order, and keeps
    every structure of ``min_atoms`` to ``max_atoms`` atoms made only of
    ``elements``, once per distinct geometry (numbers, positions, cell, pbc).
    Every ``validation_every``-th kept structure goes to validation. Returns
    the training and validation ASE ``Atoms`` (geometry only; the job labels
    them with the foundation model) and a record of the selection.
    """

    from ase.io import read

    allowed = set(elements)
    kept, seen, total = [], set(), 0
    digest = hashlib.sha256()
    with Path(path).open("rb") as raw:
        for chunk in iter(lambda: raw.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    with Path(path).open() as stream:
        while True:
            first = stream.readline()
            if not first:
                break
            if not first.strip():
                continue
            count = int(first)
            header = stream.readline()
            rows = [stream.readline() for _ in range(count)]
            index, total = total, total + 1
            if not min_atoms <= count <= max_atoms or not {row.split()[0] for row in rows} <= allowed:
                continue
            atoms = read(io.StringIO(first + header + "".join(rows)), format="extxyz")
            key = hashlib.sha256(
                atoms.numbers.tobytes() + atoms.positions.tobytes() + atoms.cell.array.tobytes() + atoms.pbc.tobytes()
            ).hexdigest()
            if key in seen:
                continue
            seen.add(key)
            clean = atoms.copy()
            clean.calc = None
            clean.info = {"replay_source_frame": index, "replay_geometry_sha256": key, "head": "pt_head"}
            kept.append(clean)
    if not kept:
        raise ValueError(f"no replay structure of {min_atoms}-{max_atoms} atoms made only of {sorted(allowed)}")
    valid = [atoms for i, atoms in enumerate(kept) if i % validation_every == validation_every - 1]
    train = [atoms for i, atoms in enumerate(kept) if i % validation_every != validation_every - 1]
    record = {
        "source": str(path), "source_sha256": digest.hexdigest(), "source_frames": total, "kept": len(kept),
        "train": len(train), "valid": len(valid), "elements": sorted(allowed),
        "atoms": [min_atoms, max_atoms], "validation_every": validation_every,
    }
    return train, valid, record


def mace_training_options(
    config: MaceFineTuneConfig,
    *,
    name: str,
    seed: int,
    atomic_numbers: Sequence[int],
    train_file: str = "train.extxyz",
    valid_file: str = "train.extxyz",
    e0s_file: str = "E0s.json",
) -> dict[str, Any]:
    """The ``mace_run_train`` options (``--config`` file) of one fine-tune.

    The job adds ``mean`` and ``std`` (the scale of the foundation head) and
    writes ``e0s_file``; the replay files are ``replay_labels/train.extxyz``
    and ``replay_labels/valid.extxyz``, written by the job.
    """

    if seed not in config.seeds:
        raise ValueError(f"seed {seed} is not one of the recipe's seeds {config.seeds}")
    options: dict[str, Any] = {
        "name": name,
        "foundation_model": config.foundation_model,
        "foundation_head": config.foundation_head,
        "multiheads_finetuning": config.replay,
        "train_file": train_file,
        "valid_file": valid_file,
        "energy_key": "REF_energy",
        "forces_key": "REF_forces",
        "E0s": e0s_file,
        "atomic_numbers": "[" + ",".join(str(number) for number in sorted(set(atomic_numbers))) + "]",
        "default_dtype": config.default_dtype,
        "device": "cuda",
        "batch_size": 1,
        "valid_batch_size": 1,
        "loss": "universal",
        "energy_weight": config.energy_weight,
        "forces_weight": config.forces_weight,
        "stress_weight": 0.0,
        "virials_weight": 0.0,
        "huber_delta": config.huber_delta,
        "lr": config.learning_rate,
        "weight_decay": 0.0,
        "ema": True,
        "ema_decay": config.ema_decay,
        "clip_grad": config.clip_grad,
        "max_num_epochs": config.epochs,
        "patience": config.epochs + 1,
        "eval_interval": 1,
        "seed": seed,
        "compute_avg_num_neighbors": False,
        "save_cpu": True,
        "log_dir": "results/logs",
        "model_dir": "results/models",
        "checkpoints_dir": "results/checkpoints",
        "results_dir": "results/metrics",
        "num_workers": 0,
        "error_table": "PerAtomRMSE",
        "scheduler": "ExponentialLR",
        "lr_scheduler_gamma": 1.0,
        "save_all_checkpoints": True,
    }
    if config.replay:
        options.update(
            pt_train_file="replay_labels/train.extxyz",
            pt_valid_file="replay_labels/valid.extxyz",
            force_mh_ft_lr=True,
            real_pt_data_ratio_threshold=0.0,
        )
    return options


@dataclass(frozen=True)
class MlipRelaxConfig:
    """How a fine-tuned model relaxes a structure: ASE ``optimizer`` until every free atom feels less than
    ``fmax_ev_per_angstrom`` or ``max_steps`` steps, with the model's ``head`` in ``default_dtype``."""

    fmax_ev_per_angstrom: float = 0.01
    max_steps: int = 2000
    optimizer: str = "FIRE"
    head: str = "Default"
    default_dtype: str = "float64"

    def __post_init__(self) -> None:
        if self.fmax_ev_per_angstrom <= 0:
            raise ValueError("MlipRelaxConfig: fmax_ev_per_angstrom must be positive")
        if isinstance(self.max_steps, bool) or not isinstance(self.max_steps, int) or self.max_steps < 1:
            raise ValueError("MlipRelaxConfig: max_steps must be a positive integer")
        if self.optimizer not in OPTIMIZERS:
            raise ValueError(f"MlipRelaxConfig: optimizer must be one of {OPTIMIZERS}, got {self.optimizer!r}")
        if not self.head:
            raise ValueError("MlipRelaxConfig: head must not be empty")
        if self.default_dtype not in ("float32", "float64"):
            raise ValueError("MlipRelaxConfig: default_dtype must be 'float32' or 'float64'")


@dataclass(frozen=True)
class MlipRelaxation:
    """Result of one relaxation: the final ``structure`` (pymatgen), its ``energy_ev``, the largest force on a
    free atom, the number of optimizer steps and whether ``fmax`` was reached."""

    structure: Any
    energy_ev: float
    max_force_ev_per_angstrom: float
    steps: int
    converged: bool


def relax_with_calculator(
    structure: Any,
    calculator: Any,
    *,
    config: MlipRelaxConfig = MlipRelaxConfig(),
    fixed_sites: Sequence[int] = (),
) -> MlipRelaxation:
    """Relax a pymatgen structure (positions only, cell fixed) with an ASE calculator.

    ``fixed_sites`` (zero-based) stay in place, as ``CalculationOverride.fixed_sites``
    in the DFT relaxations. ``config.head`` and ``config.default_dtype`` are
    for building the calculator; this function uses the calculator it gets.
    """

    from ase.constraints import FixAtoms
    from ase.optimize import BFGS, FIRE, LBFGS
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

    plain = Structure(structure.lattice, structure.species, structure.frac_coords)
    fixed = sorted(set(int(index) for index in fixed_sites))
    if any(index < 0 or index >= len(plain) for index in fixed):
        raise ValueError(f"fixed_sites must be within [0, {len(plain)})")
    atoms = AseAtomsAdaptor.get_atoms(plain)
    atoms.pbc = True
    if fixed:
        atoms.set_constraint(FixAtoms(indices=fixed))
    atoms.calc = calculator
    optimizer = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}[config.optimizer](atoms, logfile=None)
    optimizer.run(fmax=config.fmax_ev_per_angstrom, steps=config.max_steps)
    forces = atoms.get_forces()  # constrained atoms report zero force
    largest = float(np.linalg.norm(forces, axis=1).max()) if len(atoms) else 0.0
    relaxed = AseAtomsAdaptor.get_structure(atoms)
    return MlipRelaxation(
        structure=Structure(relaxed.lattice, relaxed.species, relaxed.frac_coords),
        energy_ev=float(atoms.get_potential_energy()),
        max_force_ev_per_angstrom=largest,
        steps=int(optimizer.get_number_of_steps()),
        converged=largest < config.fmax_ev_per_angstrom,
    )


def partially_occupied_count(occupations: Any, *, tolerance: float = 0.01) -> int:
    """Number of states filled between ``tolerance`` and ``1 - tolerance`` of the fullest state.

    ``occupations`` is any array of occupation numbers (per spin, k-point and
    band) of a converged smeared SCF. More than zero means partly filled
    states at the Fermi level: a metal. For CP2K, use
    ``Cp2kDos.partially_occupied_bins``, which applies the same rule to the
    occupations per energy bin.
    """

    if not 0 < tolerance < 0.5:
        raise ValueError("tolerance must be between 0 and 0.5")
    values = np.asarray(occupations, dtype=float).ravel()
    full = float(values.max()) if values.size else 0.0
    if full <= 0:
        raise ValueError("no occupied states")
    filling = values / full
    return int(((filling > tolerance) & (filling < 1 - tolerance)).sum())


def require_metallic(label: str, *, partially_occupied: int) -> None:
    """Raise ``ValueError`` unless the bulk ``label`` has partly filled states (a metal).

    The MLIP route is validated for metallic systems only; semiconductors and
    insulators need a separate treatment.
    """

    if partially_occupied <= 0:
        raise ValueError(
            f"{label}: the bulk has no partly filled states at the Fermi level, so it is not a metal; "
            "the MLIP route of psteros is for metallic systems only"
        )


def check_one_energy_source(sources: Mapping[str, str]) -> str:
    """Return the single energy source of a set of energies, or raise ``ValueError`` naming the labels.

    ``sources`` maps each label (slabs and references alike) to ``"dft"``,
    ``"mlip"`` (the model at its own geometry) or ``"mlip_dft"`` (a DFT static
    at the model's geometry).
    """

    if not sources:
        raise ValueError("no energies given")
    unknown = sorted(label for label, source in sources.items() if source not in ENERGY_SOURCES)
    if unknown:
        raise ValueError(f"unknown energy source for {unknown}; use one of {ENERGY_SOURCES}")
    groups: dict[str, list[str]] = {}
    for label, source in sources.items():
        groups.setdefault(source, []).append(label)
    if len(groups) > 1:
        detail = "; ".join(f"{source}: {sorted(labels)}" for source, labels in sorted(groups.items()))
        raise ValueError(
            "every energy of one surface phase diagram must come from one source "
            f"(slabs and references alike), got {detail}"
        )
    return next(iter(groups))


def prepare_finetune_job(
    directory: str | Path,
    frames: Sequence[LabelledStructure],
    config: MaceFineTuneConfig,
    *,
    name: str,
    replay_train: Sequence[Any] = (),
    replay_valid: Sequence[Any] = (),
    device: str = "cuda",
    seeds: Sequence[int] | None = None,
) -> Path:
    """Write the folder of one fine-tune job (``job.json``, ``train.extxyz``, replay geometries, ``finetune.py``).

    One run per seed (``config.seeds`` unless ``seeds`` picks some), named
    ``<name>_seed<seed>``. ``replay_train`` and ``replay_valid`` are the ASE
    structures of :func:`select_replay_structures`, needed when
    ``config.replay``. The training frames are the ``"train"`` frames; the
    validation file is the training file (validation is only reported).
    Returns the folder, which must not exist yet.
    """

    from psteros.mlip_jobs import script_path

    train = [frame for frame in frames if frame.role == "train"]
    if not train:
        raise ValueError("no training frames (role 'train')")
    if config.replay and (not replay_train or not replay_valid):
        raise ValueError("a replay fine-tune needs replay_train and replay_valid structures")
    if not config.replay and (replay_train or replay_valid):
        raise ValueError("replay structures were given to a fine-tune without replay")
    if device not in ("cuda", "cpu"):
        raise ValueError("device must be 'cuda' or 'cpu'")
    chosen = tuple(config.seeds if seeds is None else seeds)
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=False)
    write_extxyz(train, folder / "train.extxyz")
    if config.replay:
        from ase.io import write

        write(folder / "replay_train_geometries.extxyz", list(replay_train), format="extxyz")
        write(folder / "replay_valid_geometries.extxyz", list(replay_valid), format="extxyz")
    from pymatgen.core import Element

    numbers = sorted({Element(element).Z for frame in train for element in frame.composition})
    runs = [
        {"name": f"{name}_seed{seed}",
         "options": mace_training_options(config, name=f"{name}_seed{seed}", seed=seed, atomic_numbers=numbers)}
        for seed in chosen
    ]
    job = {
        "foundation_model": config.foundation_model, "foundation_model_sha256": config.foundation_model_sha256,
        "foundation_head": config.foundation_head, "replay": config.replay, "device": device, "runs": runs,
    }
    (folder / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    (folder / "finetune.py").write_text(script_path("finetune.py").read_text())
    return folder


def prepare_relax_job(
    directory: str | Path,
    structures: Mapping[str, Any],
    models: Mapping[str, tuple[str, str]],
    config: MlipRelaxConfig = MlipRelaxConfig(),
    *,
    fixed_sites: Mapping[str, Sequence[int]] | None = None,
    device: str = "cuda",
) -> Path:
    """Write the folder of one relaxation job: every model of ``models`` (name -> (path, head)) relaxes every
    pymatgen structure of ``structures`` (label -> structure), keeping ``fixed_sites[label]`` in place."""

    from psteros.mlip_jobs import script_path

    if not structures or not models:
        raise ValueError("give at least one structure and one model")
    fixed_sites = dict(fixed_sites or {})
    unknown = sorted(set(fixed_sites).difference(structures))
    if unknown:
        raise ValueError(f"fixed_sites names labels that are not in structures: {unknown}")
    if device not in ("cuda", "cpu"):
        raise ValueError("device must be 'cuda' or 'cpu'")
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=False)
    items = []
    for label, structure in structures.items():
        frame = LabelledStructure(label, structure, 0.0, np.zeros((len(structure), 3)))
        write_extxyz([frame], folder / f"{label}.extxyz")
        fixed = sorted(set(int(index) for index in fixed_sites.get(label, ())))
        if any(index < 0 or index >= len(structure) for index in fixed):
            raise ValueError(f"{label}: fixed_sites must be within [0, {len(structure)})")
        items.append({"label": label, "file": f"{label}.extxyz", "fixed": fixed})
    job = {
        "device": device, "default_dtype": config.default_dtype, "optimizer": config.optimizer,
        "fmax_ev_per_angstrom": config.fmax_ev_per_angstrom, "max_steps": config.max_steps,
        "models": {name: {"path": str(path), "head": head} for name, (path, head) in models.items()},
        "structures": items,
    }
    (folder / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    (folder / "relax.py").write_text(script_path("relax.py").read_text())
    return folder
