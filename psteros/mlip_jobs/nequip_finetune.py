"""psteros NequIP/Allegro fine-tune job: runs on the GPU computer and needs only ase, numpy, torch and nequip
(and nequip-allegro for Allegro models).

Files in the working directory: ``job.json`` (written by ``psteros.nequip.prepare_nequip_finetune_job``),
``train.extxyz`` (DFT labels, REF_energy in eV and REF_forces in eV/A) and, with replay,
``replay_train_geometries.extxyz`` and ``replay_valid_geometries.extxyz`` (geometries only). One run per entry of
``job.json["runs"]``, in ``runs/<name>/``; the model of each run is packaged to ``models/<name>.nequip.zip``.
``summary.json`` records every step.

1. The foundation package is checked (size and SHA-256 when given).
2. The foundation model predicts the training frames and labels the replay geometries.
3. Per-type energy shifts of the new model: the foundation's plus per-element shifts δ_Z fitted in least squares to
   E_DFT − E_foundation over the training frames; the replay labels get the same shifts, so the DFT and the replay
   labels share one reference.
4. Datasets: extended-XYZ training set (DFT frames and replay structures) and validation set (replay validation
   structures, or the DFT frames without replay).
5. Per run: ``nequip-train`` from the package with the shifted per-type energies (nequip's ``modify_PerTypeScaleShift``).
   A first pass with a zero learning rate checks that this model reproduces the foundation model shifted by
   Σ_Z n_Z δ_Z (1e-5 eV/atom, 1e-4 eV/A) on the first training frame and the first replay structure, through the
   checkpoint the training writes; then the training proper, whose last checkpoint is packaged.
6. Evaluation of each model: energy and force errors on the training frames, force RMSE against the foundation on
   the replay validation structures.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import yaml
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
DEVICE = JOB.get("device", "cuda")
SUMMARY: dict = {"status": "starting", "runs": {}}
INITIAL_ENERGY_TOLERANCE = 1e-5  # eV/atom
INITIAL_FORCE_TOLERANCE = 1e-4  # eV/A


def dump() -> None:
    (ROOT / "summary.json").write_text(json.dumps(SUMMARY, indent=2, allow_nan=False) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def calculator(path: Path):
    """A checkpoint or package of the NequIP framework as an ASE calculator, run eagerly (nothing compiled)."""

    from nequip.integrations.ase import NequIPCalculator

    return NequIPCalculator._from_saved_model(str(path), device=DEVICE, chemical_species_to_atom_type_map=True)


def predict(calc, atoms):
    atoms = atoms.copy()
    atoms.calc = calc
    return float(atoms.get_potential_energy()), np.asarray(atoms.get_forces(apply_constraint=False))


def counts(atoms, numbers):
    return np.array([(atoms.numbers == z).sum() for z in numbers], dtype=float)


def write_dataset(path: Path, items) -> None:
    labelled = []
    for atoms, energy, forces in items:
        copy = atoms.copy()
        copy.info = {}
        copy.arrays = {key: value for key, value in copy.arrays.items() if key in ("numbers", "positions")}
        copy.calc = SinglePointCalculator(copy, energy=float(energy), forces=np.asarray(forces))
        labelled.append(copy)
    write(path, labelled, format="extxyz")


def foundation_shifts(path: Path):
    """The per-type energy shifts (eV per atom) and type names of the foundation model."""

    import torch
    from nequip.model.saved_models.load_utils import load_saved_model

    model = load_saved_model(str(path))
    modules = [module for module in model.modules()
               if type(module).__name__ == "PerTypeScaleShift" and getattr(module, "has_shifts", False)]
    assert len(modules) == 1, f"expected one per-type energy shift in the foundation model, found {len(modules)}"
    (module,) = modules
    shifts = module.shifts.detach().to(torch.float64).reshape(-1).cpu().numpy()
    if shifts.size == 1:
        shifts = np.full(len(module.type_names), shifts[0])
    return dict(zip(module.type_names, shifts.tolist()))


def run_config(folder: Path, seed: int, learning_rate: float, epochs: int, shifts: dict, cutoff: float,
               type_names: list, train_file: Path, valid_file: Path) -> Path:
    training = JOB["training"]
    loader = {"_target_": "torch.utils.data.DataLoader", "batch_size": training["batch_size"], "num_workers": 0}
    config = {
        "run": ["train"],
        "data": {
            "_target_": "nequip.data.datamodule.ASEDataModule", "seed": seed,
            "train_file_path": str(train_file), "val_file_path": str(valid_file),
            "transforms": [
                {"_target_": "nequip.data.transforms.ChemicalSpeciesToAtomTypeMapper", "model_type_names": type_names},
                {"_target_": "nequip.data.transforms.NeighborListTransform", "r_max": cutoff},
            ],
            "train_dataloader": {**loader, "shuffle": True},
            "val_dataloader": loader,
        },
        "trainer": {
            "_target_": "lightning.Trainer", "accelerator": "gpu" if DEVICE == "cuda" else "cpu", "devices": 1,
            "max_epochs": epochs, "logger": False, "enable_progress_bar": False, "num_sanity_val_steps": 0,
            "callbacks": [{"_target_": "lightning.pytorch.callbacks.ModelCheckpoint", "dirpath": str(folder),
                           "save_last": True, "save_top_k": 0}],
        },
        "training_module": {
            "_target_": "nequip.train.EMALightningModule", "ema_decay": training["ema_decay"],
            "loss": {"_target_": "nequip.train.EnergyForceLoss", "per_atom_energy": True,
                     "coeffs": {"total_energy": training["energy_weight"], "forces": training["forces_weight"]}},
            "val_metrics": {"_target_": "nequip.train.EnergyForceMetrics",
                            "coeffs": {"total_energy_mae": 1.0, "forces_mae": 1.0}},
            "optimizer": {"_target_": "torch.optim.Adam", "lr": learning_rate},
            "model": {
                "_target_": "nequip.model.modify",
                "modifiers": [{"modifier": "modify_PerTypeScaleShift", "shifts": shifts,
                               "shifts_trainable": False, "scales_trainable": False}],
                "model": {"_target_": "nequip.model.ModelFromPackage", "package_path": JOB["foundation_model"]},
            },
        },
    }
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    return folder


def train(folder: Path) -> Path:
    """``nequip-train`` in ``folder``; returns its last checkpoint."""

    with (folder / "train.log").open("w") as log:
        # the modules of nequip-train and nequip-package with this Python: the job works wherever the environment is
        done = subprocess.run([sys.executable, "-m", "nequip.scripts.train", "-cp", str(folder), "-cn", "config.yaml",
                               f"hydra.run.dir={folder}"], cwd=folder, stdout=log, stderr=subprocess.STDOUT)
    if done.returncode:
        raise RuntimeError(f"nequip-train failed in {folder}: see {folder / 'train.log'}")
    last = folder / "last.ckpt"
    assert last.is_file(), f"no last checkpoint in {folder}"
    return last


try:
    # ---- 1. foundation package
    base_path = Path(JOB["foundation_model"])
    assert base_path.is_file(), f"no foundation package at {base_path}"
    SUMMARY["foundation_model"] = {"path": str(base_path), "bytes": base_path.stat().st_size}
    if JOB.get("foundation_model_sha256"):
        SUMMARY["foundation_model"]["sha256"] = sha256(base_path)
        assert SUMMARY["foundation_model"]["sha256"] == JOB["foundation_model_sha256"], "foundation model changed"
    dump()

    frames = read(ROOT / "train.extxyz", ":")
    numbers = sorted({int(z) for atoms in frames for z in atoms.numbers})
    replay = bool(JOB["replay"])

    # ---- 2. foundation predictions and replay labels
    base_calc = calculator(base_path)
    base_pred = [predict(base_calc, atoms) for atoms in frames]
    replay_sets = {}
    if replay:
        for split in ("train", "valid"):
            geometries = read(ROOT / f"replay_{split}_geometries.extxyz", ":")
            replay_sets[split] = [(atoms, *predict(base_calc, atoms)) for atoms in geometries]
        SUMMARY["replay"] = {split: len(items) for split, items in replay_sets.items()}
    del base_calc

    # ---- 3. per-type energy shifts of the new model
    from ase.data import chemical_symbols

    matrix = np.array([counts(atoms, numbers) for atoms in frames])
    residuals = np.array([atoms.info["REF_energy"] - energy for atoms, (energy, _) in zip(frames, base_pred)])
    assert np.linalg.matrix_rank(matrix) == len(numbers), "the training compositions do not separate the elements"
    deltas = np.linalg.lstsq(matrix, residuals, rcond=None)[0]
    base_shifts = foundation_shifts(base_path)
    type_names = list(base_shifts)
    symbols = [chemical_symbols[z] for z in numbers]
    missing = sorted(set(symbols).difference(type_names))
    assert not missing, f"the foundation model has no atom type for {missing} (its types are chemical symbols)"
    shifts = {symbol: base_shifts[symbol] + float(delta) for symbol, delta in zip(symbols, deltas)}
    from nequip.model.saved_models.load_utils import load_saved_model
    from nequip.nn import graph_model

    cutoff = float(load_saved_model(str(base_path)).metadata[graph_model.R_MAX_KEY])
    SUMMARY["references"] = {
        "numbers": numbers, "foundation_shifts_eV": {s: base_shifts[s] for s in symbols},
        "shifts_eV": deltas.tolist(), "per_type_shifts_eV": shifts, "cutoff_angstrom": cutoff,
        "residual_rms_after_shift_eV": float(np.sqrt(np.mean((residuals - matrix @ deltas) ** 2))),
        "per_frame_DFT_minus_foundation_eV": {atoms.info.get("label", str(i)): float(r)
                                              for i, (atoms, r) in enumerate(zip(frames, residuals))},
    }
    dump()

    # ---- 4. datasets
    data = ROOT / "data"
    data.mkdir()
    dft = [(atoms, atoms.info["REF_energy"], atoms.arrays["REF_forces"]) for atoms in frames]
    shifted = {split: [(atoms, energy + float(counts(atoms, numbers) @ deltas), forces)
                       for atoms, energy, forces in items] for split, items in replay_sets.items()}
    assert all(set(atoms.numbers) <= set(numbers) for items in shifted.values() for atoms, _, _ in items), \
        "a replay structure has an element without a DFT label"
    write_dataset(data / "train.extxyz", dft + shifted.get("train", []))
    write_dataset(data / "val.extxyz", shifted["valid"] if replay else dft)

    # ---- 5. one fine-tune per run
    training = JOB["training"]
    (ROOT / "models").mkdir(exist_ok=True)
    SUMMARY["status"] = "training"
    for spec in JOB["runs"]:
        name, seed = spec["name"], int(spec["seed"])
        record = SUMMARY["runs"].setdefault(name, {"seed": seed, "status": "starting"})
        folder = ROOT / "runs" / name
        common = dict(seed=seed, shifts=shifts, cutoff=cutoff, type_names=type_names,
                      train_file=data / "train.extxyz", valid_file=data / "val.extxyz")
        # the check: one pass at a zero learning rate leaves the initialised model, saved like the trained one
        start = time.monotonic()
        initial = calculator(train(run_config(folder / "check", learning_rate=0.0, epochs=1, **common)))
        energy, forces = predict(initial, frames[0])
        expected = base_pred[0][0] + float(counts(frames[0], numbers) @ deltas)
        record.update(initial_energy_error_eV_atom=abs(energy - expected) / len(frames[0]),
                      initial_force_error_eV_A=float(np.max(np.abs(forces - base_pred[0][1]))))
        if replay:
            atoms, energy_ref, forces_ref = shifted["train"][0]
            energy, forces = predict(initial, atoms)
            record.update(replay_initial_energy_error_eV_atom=abs(energy - energy_ref) / len(atoms),
                          replay_initial_force_error_eV_A=float(np.max(np.abs(forces - forces_ref))))
        del initial
        shutil.rmtree(folder / "check")
        dump()
        assert record["initial_energy_error_eV_atom"] < INITIAL_ENERGY_TOLERANCE, record
        assert record["initial_force_error_eV_A"] < INITIAL_FORCE_TOLERANCE, record
        if replay:
            assert record["replay_initial_energy_error_eV_atom"] < INITIAL_ENERGY_TOLERANCE, record
            assert record["replay_initial_force_error_eV_A"] < INITIAL_FORCE_TOLERANCE, record
        record["initialization_passed"] = True
        record["status"] = "training"
        dump()
        last = train(run_config(folder / "train", learning_rate=training["learning_rate"], epochs=training["epochs"],
                                **common))
        final = ROOT / "models" / f"{name}.nequip.zip"
        with (folder / "package.log").open("w") as log:
            done = subprocess.run([sys.executable, "-m", "nequip.scripts.package", "build", str(last), str(final)],
                                  cwd=folder, stdout=log, stderr=subprocess.STDOUT)
        assert done.returncode == 0 and final.is_file(), f"nequip-package failed: see {folder / 'package.log'}"
        record.update(model=str(final.relative_to(ROOT)), model_sha256=sha256(final),
                      elapsed_seconds=time.monotonic() - start, status="trained")
        dump()
        # ---- 6. evaluation
        model = calculator(final)
        rows = {}
        for atoms in frames:
            energy, forces = predict(model, atoms)
            rows[atoms.info.get("label", str(len(rows)))] = {
                "natoms": len(atoms), "energy_error_eV": energy - atoms.info["REF_energy"],
                "force_rmse_meV_A": float(1000 * np.sqrt(((forces - atoms.arrays["REF_forces"]) ** 2).mean()))}
        evaluation = {"train": rows}
        if replay:
            evaluation["replay_valid_vs_foundation"] = {
                str(i): {"natoms": len(atoms), "force_rmse_meV_A": float(
                    1000 * np.sqrt(((predict(model, atoms)[1] - forces_ref) ** 2).mean()))}
                for i, (atoms, _, forces_ref) in enumerate(shifted["valid"])}
        (folder / "evaluation.json").write_text(json.dumps(evaluation, indent=2) + "\n")
        for checkpoint in (folder / "train").glob("*.ckpt"):  # the package in models/ is the model
            checkpoint.unlink()
        record["status"] = "completed"
        dump()
    SUMMARY["status"] = "completed"
except Exception:
    SUMMARY.update(status="failed", traceback=traceback.format_exc())
    raise
finally:
    dump()
