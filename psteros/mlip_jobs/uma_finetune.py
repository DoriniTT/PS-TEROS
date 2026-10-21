"""psteros UMA fine-tune job: runs on the GPU computer and needs only ase, numpy, torch and fairchem-core.

Files in the working directory: ``job.json`` (written by ``psteros.uma.prepare_uma_finetune_job``), ``train.extxyz``
(DFT labels, REF_energy in eV and REF_forces in eV/A) and, with replay, ``replay_train_geometries.extxyz`` and
``replay_valid_geometries.extxyz`` (geometries only). One run per entry of ``job.json["runs"]``, in ``runs/<name>/``;
the model of each run is copied to ``models/<name>.model`` (a fairchem inference checkpoint). ``summary.json``
records every step.

1. The foundation checkpoint is checked (size and SHA-256 when given).
2. The foundation model (its ``task``) predicts the training frames and labels the replay geometries.
3. Element references of the new head: the foundation's for ``task`` plus per-element shifts δ_Z fitted in least
   squares to E_DFT − E_foundation over the training frames; the replay labels get the same shifts, so the DFT and
   the replay labels share one reference.
4. Datasets: ASE-LMDB training set (DFT frames and replay structures) and validation set (replay validation
   structures, or the DFT frames without replay).
5. Per run: fairchem's trainer with a new energy head that starts as a copy of the foundation's, and the
   foundation's normalizer. A first pass with a zero learning rate checks that the initialised model reproduces
   the foundation model shifted by Σ_Z n_Z δ_Z (1e-5 eV/atom, 1e-4 eV/A in float32) on the first training frame and
   the first replay structure, through the same checkpoint the training writes; then the training proper.
6. Evaluation of each model: energy and force errors on the training frames, force RMSE against the foundation on
   the replay validation structures.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import torch
import yaml
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.io import read

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
DEVICE = JOB.get("device", "cuda")
TASK = JOB["task"]
SUMMARY: dict = {"status": "starting", "runs": {}}
INITIAL_ENERGY_TOLERANCE = 1e-5  # eV/atom; UMA runs in float32
INITIAL_FORCE_TOLERANCE = 1e-4  # eV/A
if DEVICE == "cuda":
    assert torch.cuda.is_available(), "job.json asks for cuda but no GPU is visible"


def dump() -> None:
    (ROOT / "summary.json").write_text(json.dumps(SUMMARY, indent=2, allow_nan=False) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def calculator(path: Path):
    """UMA at the training precision (float32), MoLE merged (exact for one composition), without torch.compile."""

    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings

    settings = InferenceSettings(tf32=False, activation_checkpointing=False, merge_mole=True, compile=False,
                                 base_precision_dtype=torch.float32)
    return FAIRChemCalculator(load_predict_unit(str(path), inference_settings=settings, device=DEVICE), task_name=TASK)


def predict(calc, atoms):
    atoms = atoms.copy()
    atoms.calc = calc
    return float(atoms.get_potential_energy()), np.asarray(atoms.get_forces(apply_constraint=False))


def counts(atoms, numbers):
    return np.array([(atoms.numbers == z).sum() for z in numbers], dtype=float)


def write_dataset(folder: Path, items) -> None:
    """An ASE-LMDB dataset of fairchem: ``data.0000.aselmdb`` and ``metadata.npz`` (atoms per structure)."""

    import fairchem.core.datasets  # noqa: F401  (registers the aselmdb backend of ase.db)

    folder.mkdir(parents=True)
    with connect(str(folder / "data.0000.aselmdb")) as db:
        for atoms, energy, forces in items:
            labelled = atoms.copy()
            labelled.calc = SinglePointCalculator(labelled, energy=float(energy), forces=np.asarray(forces))
            labelled.info = {}
            db.write(labelled)
    np.savez_compressed(folder / "metadata.npz", natoms=np.array([len(atoms) for atoms, _, _ in items]))


def task_config(name: str, prop: str, level: str, coefficient: float, loss: str, rmsd: float, references=None) -> dict:
    config = {
        "_target_": "fairchem.core.units.mlip_unit.mlip_unit.Task",
        "name": name, "level": level, "property": prop,
        "loss_fn": {"_target_": "fairchem.core.modules.loss.DDPMTLoss",
                    "loss_fn": {"_target_": loss}, "coefficient": coefficient},
        "out_spec": {"dim": [1] if prop == "energy" else [3], "dtype": "float32"},
        "normalizer": {"_target_": "fairchem.core.modules.normalization.normalizer.Normalizer",
                       "mean": 0.0, "rmsd": float(rmsd)},
        "datasets": [TASK],
        "metrics": ["mae", "per_atom_mae"] if prop == "energy" else ["mae"],
    }
    if prop == "forces":
        config.update(train_on_free_atoms=True, eval_on_free_atoms=True)
        config["loss_fn"]["reduction"] = "mean"
    if references is not None:
        config["element_references"] = {
            "_target_": "fairchem.core.modules.normalization.element_references.ElementReferences",
            "element_references": {"_target_": "torch.DoubleTensor", "_args_": [list(map(float, references))]},
        }
    return config


def foundation_task(checkpoint, prop: str) -> dict:
    """The task of the foundation checkpoint for ``TASK`` and ``prop`` (energy or forces), as a plain dict."""

    from omegaconf import OmegaConf

    for item in checkpoint.tasks_config:
        item = OmegaConf.to_container(item, resolve=True) if not isinstance(item, dict) else item
        if item.get("property") == prop and TASK in item.get("datasets", []):
            return item
    raise ValueError(f"the foundation checkpoint has no {prop} task for {TASK!r}")


def run_config(run_dir: Path, run_id: str, seed: int, learning_rate: float, epochs: int, energy_refs, rmsd_energy,
               rmsd_forces, rmsd_stress, train_dir: Path, valid_dir: Path, backbone_overrides: dict) -> Path:
    training = JOB["training"]
    stress = task_config("stress", "stress", "system", 1.0, "fairchem.core.modules.loss.L2NormLoss", rmsd_stress)
    stress.update(out_spec={"dim": [1, 9], "dtype": "float32"}, inference_only=True)
    del stress["loss_fn"], stress["metrics"]
    tasks = [  # as fairchem's energy+forces fine-tune; the stress is predicted, not trained
        task_config("energy", "energy", "system", training["energy_weight"], "fairchem.core.modules.loss.PerAtomMAELoss",
                    rmsd_energy, energy_refs),
        task_config("forces", "forces", "atom", training["forces_weight"], "fairchem.core.modules.loss.L2NormLoss",
                    rmsd_forces),
        stress,
    ]
    dataset = {"splits": {}, "format": "ase_db", "transforms": {"common_transform": {"dataset_name": TASK}}}

    def concat(split: str, src: Path) -> dict:
        return {"_target_": "fairchem.core.datasets.mt_concat_dataset.create_concat_dataset",
                "dataset_configs": {TASK: {**dataset, "splits": {split: {"src": str(src)}}}},
                "combined_dataset_config": {"sampling": {"type": "temperature", "temperature": 1.0}}}

    def loader(dataset_key: str, shuffle: bool) -> dict:
        return {"_target_": "fairchem.core.components.common.dataloader_builder.get_dataloader",
                "dataset": f"${{{dataset_key}}}",
                "batch_sampler_fn": {"_target_": "fairchem.core.common.data_parallel.BalancedBatchSampler",
                                     "_partial_": True, "batch_size": training["batch_size"], "shuffle": shuffle,
                                     "seed": seed},
                "num_workers": 0,
                "collate_fn": {"_target_": "fairchem.core.units.mlip_unit.mlip_unit.mt_collater_adapter",
                               "tasks": "${tasks}"}}

    steps_per_epoch = None  # epochs drive the schedule
    config = {
        "job": {"device_type": "CUDA" if DEVICE == "cuda" else "CPU",
                "scheduler": {"mode": "LOCAL", "ranks_per_node": 1, "num_nodes": 1},
                "debug": True, "run_dir": str(run_dir), "run_name": run_id, "timestamp_id": run_id,
                "seed": seed, "logger": None},
        "tasks": tasks,
        "train_dataset": concat("train", train_dir),
        "val_dataset": concat("val", valid_dir),
        "train_dataloader": loader("train_dataset", True),
        "eval_dataloader": loader("val_dataset", False),
        "runner": {
            "_target_": "fairchem.core.components.train.train_runner.TrainEvalRunner",
            "train_dataloader": "${train_dataloader}", "eval_dataloader": "${eval_dataloader}",
            "train_eval_unit": {
                "_target_": "fairchem.core.units.mlip_unit.mlip_unit.MLIPTrainEvalUnit",
                "job_config": "${job}", "tasks": "${tasks}",
                "model": {"_target_": "psteros_uma_init.initialize_from_foundation", "task": TASK,
                          "checkpoint_location": JOB["foundation_model"], "overrides": {"backbone": backbone_overrides},
                          "heads": {"efs": {"module": "fairchem.core.models.uma.escn_md.MLP_EFS_Head"}}},
                "optimizer_fn": {"_target_": "torch.optim.AdamW", "_partial_": True, "lr": learning_rate,
                                 "weight_decay": training["weight_decay"]},
                "cosine_lr_scheduler_fn": {"_target_": "fairchem.core.units.mlip_unit.mlip_unit._get_consine_lr_scheduler",
                                           "_partial_": True, "warmup_factor": 0.2, "warmup_epochs": 0.01,
                                           "lr_min_factor": 0.01, "epochs": epochs, "steps": steps_per_epoch},
                "ema_decay": training["ema_decay"], "print_every": 10, "clip_grad_norm": 100},
            "max_epochs": epochs, "max_steps": None, "evaluate_every_n_steps": 1000000,
            "callbacks": [{"_target_": "fairchem.core.components.callbacks.TrainCheckpointCallback",
                           "checkpoint_every_n_steps": 1000000, "max_saved_checkpoints": 1}],
        },
    }
    path = run_dir / f"{run_id}.yaml"
    run_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


INIT_MODULE = '''"""Initialise a fairchem fine-tune as the foundation model that fairchem's inference runs (written by psteros).

fairchem's predictions load the exponential moving average (EMA) of a checkpoint's weights, while its fine-tune starts
from the raw weights and a new random energy head. Here the backbone and the energy head of the new model are the EMA
weights of the foundation; a head with one expert per task (UMA's DatasetSpecificMoEWrapper) is first merged on the
task of the fine-tune, as fairchem's inference merges it, so the untrained model is the one the foundation's
predictions for that task come from.
"""
from types import SimpleNamespace

import torch


def initialize_from_foundation(checkpoint_location, overrides=None, heads=None, strict=True, task=None):
    from fairchem.core.units.mlip_unit.mlip_unit import initialize_finetuning_model
    from fairchem.core.units.mlip_unit.utils import load_inference_model

    averaged, _ = load_inference_model(checkpoint_location, overrides, use_ema=True, strict=strict)
    pretrained = averaged.module
    (source,) = pretrained.output_heads.values()
    if hasattr(source, "merge_MOLE_model") and hasattr(source, "dataset_name_to_exp"):  # one expert per task
        if task not in source.dataset_name_to_exp:
            raise ValueError(f"the foundation head has no expert for task {task!r}: {sorted(source.dataset_name_to_exp)}")
        source.merge_MOLE_model(SimpleNamespace(dataset=[task], pos=torch.zeros(1, 3)))
    source = getattr(source, "head", source)  # the MLP_EFS_Head inside the wrapper
    if not isinstance(getattr(source, "energy_block", None), torch.nn.Sequential):
        raise TypeError(f"the foundation head {type(source).__name__} has no plain energy block to copy")
    model = initialize_finetuning_model(checkpoint_location, overrides, heads, strict)
    model.backbone.load_state_dict(pretrained.backbone.state_dict())
    (target,) = model.output_heads.values()
    target.energy_block.load_state_dict(source.energy_block.state_dict())
    return model
'''


def train(config_path: Path) -> None:
    """fairchem's CLI in this process (one rank), as its own tests launch it: a clean distributed and Hydra state."""

    import hydra
    from fairchem.core._cli import main
    from fairchem.core.common import distutils, gp_utils

    def clean() -> None:
        if gp_utils.initialized():
            gp_utils.cleanup_gp()
        distutils.cleanup()
        hydra.core.global_hydra.GlobalHydra.instance().clear()

    clean()
    try:
        main(type("Args", (), {"config": str(config_path)})(), [])
    finally:
        clean()


try:
    # ---- 1. foundation checkpoint
    base_path = Path(JOB["foundation_model"])
    assert base_path.is_file(), f"no foundation checkpoint at {base_path}"
    SUMMARY["foundation_model"] = {"path": str(base_path), "bytes": base_path.stat().st_size, "task": TASK}
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

    # ---- 3. element references of the new head
    matrix = np.array([counts(atoms, numbers) for atoms in frames])
    residuals = np.array([atoms.info["REF_energy"] - energy for atoms, (energy, _) in zip(frames, base_pred)])
    assert np.linalg.matrix_rank(matrix) == len(numbers), "the training compositions do not separate the elements"
    shifts = np.linalg.lstsq(matrix, residuals, rcond=None)[0]
    checkpoint = torch.load(base_path, map_location="cpu", weights_only=False)
    energy_task, forces_task = foundation_task(checkpoint, "energy"), foundation_task(checkpoint, "forces")
    stress_task = foundation_task(checkpoint, "stress")
    references = energy_task.get("element_references")  # UMA's tasks have them; a model without has zeros
    base_refs = (np.array(references["element_references"]["_args_"][0], dtype=float) if references
                 else np.zeros(max(100, numbers[-1] + 1)))
    refs = base_refs.copy()
    refs[numbers] += shifts
    backbone = dict(checkpoint.model_config["backbone"])
    # fairchem's fine-tune settings (graphs built on the fly, stress regressed), with the foundation's neighbours
    backbone_overrides = {"otf_graph": True, "always_use_pbc": False, "regress_stress": True}
    if "max_neighbors" in backbone:
        backbone_overrides["max_neighbors"] = backbone["max_neighbors"]
    del checkpoint
    SUMMARY["references"] = {
        "numbers": numbers, "foundation_refs_eV": base_refs[numbers].tolist(), "shifts_eV": shifts.tolist(),
        "residual_rms_after_shift_eV": float(np.sqrt(np.mean((residuals - matrix @ shifts) ** 2))),
        "per_frame_DFT_minus_foundation_eV": {atoms.info.get("label", str(i)): float(r)
                                              for i, (atoms, r) in enumerate(zip(frames, residuals))},
        "normalizer_rmsd": {"energy": energy_task["normalizer"]["rmsd"], "forces": forces_task["normalizer"]["rmsd"],
                            "stress": stress_task["normalizer"]["rmsd"]},
        "backbone_overrides": backbone_overrides,
    }
    dump()

    # ---- 4. datasets
    data = ROOT / "data"
    dft = [(atoms, atoms.info["REF_energy"], atoms.arrays["REF_forces"]) for atoms in frames]
    shifted = {split: [(atoms, energy + float(counts(atoms, numbers) @ shifts) if set(atoms.numbers) <= set(numbers)
                        else np.nan, forces) for atoms, energy, forces in items]
               for split, items in replay_sets.items()}
    assert all(np.isfinite(energy) for items in shifted.values() for _, energy, _ in items), \
        "a replay structure has an element without a DFT label"
    write_dataset(data / "train", dft + shifted.get("train", []))
    write_dataset(data / "val", shifted["valid"] if replay else dft)
    (ROOT / "psteros_uma_init.py").write_text(INIT_MODULE)
    sys.path.insert(0, str(ROOT))
    del base_calc
    if DEVICE == "cuda":
        torch.cuda.empty_cache()

    # ---- 5. one fine-tune per run
    training = JOB["training"]
    (ROOT / "models").mkdir(exist_ok=True)
    SUMMARY["status"] = "training"
    for spec in JOB["runs"]:
        name, seed = spec["name"], int(spec["seed"])
        record = SUMMARY["runs"].setdefault(name, {"seed": seed, "status": "starting"})
        folder = ROOT / "runs" / name
        common = dict(seed=seed, energy_refs=refs, rmsd_energy=energy_task["normalizer"]["rmsd"],
                      rmsd_forces=forces_task["normalizer"]["rmsd"], rmsd_stress=stress_task["normalizer"]["rmsd"],
                      train_dir=data / "train", valid_dir=data / "val",
                      backbone_overrides=backbone_overrides)
        # the check: one pass at a zero learning rate leaves the initialised model, saved like the trained one
        start = time.monotonic()
        train(run_config(folder / "check", "check", learning_rate=0.0, epochs=1, **common))
        initial = calculator(folder / "check" / "check" / "checkpoints" / "final" / "inference_ckpt.pt")
        energy, forces = predict(initial, frames[0])
        expected = base_pred[0][0] + float(counts(frames[0], numbers) @ shifts)
        record.update(initial_energy_error_eV_atom=abs(energy - expected) / len(frames[0]),
                      initial_force_error_eV_A=float(np.max(np.abs(forces - base_pred[0][1]))))
        if replay:
            atoms, energy_ref, forces_ref = shifted["train"][0]
            energy, forces = predict(initial, atoms)
            record.update(replay_initial_energy_error_eV_atom=abs(energy - energy_ref) / len(atoms),
                          replay_initial_force_error_eV_A=float(np.max(np.abs(forces - forces_ref))))
        del initial
        shutil.rmtree(folder / "check")  # its checkpoints are the size of the foundation's
        dump()
        assert record["initial_energy_error_eV_atom"] < INITIAL_ENERGY_TOLERANCE, record
        assert record["initial_force_error_eV_A"] < INITIAL_FORCE_TOLERANCE, record
        if replay:
            assert record["replay_initial_energy_error_eV_atom"] < INITIAL_ENERGY_TOLERANCE, record
            assert record["replay_initial_force_error_eV_A"] < INITIAL_FORCE_TOLERANCE, record
        record["initialization_passed"] = True
        record["status"] = "training"
        dump()
        train(run_config(folder, "train", learning_rate=training["learning_rate"], epochs=training["epochs"], **common))
        exported = folder / "train" / "checkpoints" / "final" / "inference_ckpt.pt"
        assert exported.is_file(), "no final checkpoint"
        final = ROOT / "models" / f"{name}.model"
        shutil.copy(exported, final)
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
        shutil.rmtree(folder / "train" / "checkpoints")  # the model is in models/; the logs and config stay
        record["status"] = "completed"
        dump()
    SUMMARY["status"] = "completed"
except Exception:
    SUMMARY.update(status="failed", traceback=traceback.format_exc())
    raise
finally:
    dump()
