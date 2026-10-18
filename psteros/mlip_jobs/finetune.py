"""psteros MACE fine-tune job: runs on the GPU computer and needs only ase, torch and mace-torch.

Files in the working directory: ``job.json`` (written by psteros), ``train.extxyz`` (DFT labels, REF_energy in eV and
REF_forces in eV/A) and, for a multi-head replay fine-tune, ``replay_train_geometries.extxyz`` and
``replay_valid_geometries.extxyz`` (geometries only). For every run of ``job.json`` (one per seed) it works in
``runs/<name>/`` and copies the final model to ``models/<name>.model``. ``summary.json`` records every step.

1. The foundation model is checked (size and SHA-256 when given).
2. The foundation model (its ``foundation_head``, float64) predicts the training frames and labels the replay
   geometries for the ``pt_head``.
3. E0s of the fine-tuned head: the foundation E0s plus per-element shifts δ_Z fitted in least squares to
   E_DFT − E_foundation over the training frames (as ``psteros.mlip.fit_element_shifts``).
4. Before training, the initialised model must reproduce the foundation model: the ``Default`` head on the first
   training frame, shifted by Σ_Z n_Z δ_Z, and the ``pt_head`` on the first replay structure (1e-8 eV/atom,
   1e-6 eV/A).
5. ``passes`` epochs; only the checkpoint of the last epoch is written, and the exported model must be that
   epoch's EMA model (validation is reported, never used to pick).
6. Evaluation: energy and force errors of the ``Default`` head on the training frames; force RMSE of both heads
   against the foundation model on the replay validation structures (what replay keeps).

Short queues (optional): with ``"persistent_root"`` in ``job.json`` the runs and models live in that folder instead
of the job folder, and a run with ``"stop_epoch"`` below its number of passes stops after that epoch and saves its
whole training state (raw weights, EMA, optimizer, scheduler, the shuffling generator, the random-number states, the
accumulated gradients of trainable parameters outside the optimizer, which enter MACE's gradient clipping, and the
epoch) to
``runs/<name>/psteros_state.pt``. The next job with the same run continues from it, so a fine-tune split into several
jobs follows the same path as one uninterrupted run. Light copies of the run's results are kept in the job folder.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import re
import shutil
import time
import traceback
from pathlib import Path

import numpy as np
import torch
from ase.io import read, write
from mace import data, tools
from mace.cli import run_train
from mace.tools import torch_geometric
from mace.tools.utils import AtomicNumberTable

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
DEVICE = JOB.get("device", "cuda")
torch.set_default_dtype(torch.float64)
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))
if DEVICE == "cuda":
    assert torch.cuda.is_available(), "job.json asks for cuda but no GPU is visible"
SUMMARY: dict = {"status": "starting", "runs": {}}


def dump() -> None:
    (ROOT / "summary.json").write_text(json.dumps(SUMMARY, indent=2, allow_nan=False) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def calculator(path, head):
    from mace.calculators import MACECalculator

    return MACECalculator(model_paths=str(path), device=DEVICE, default_dtype="float64", head=head)


def predict(calc, atoms):
    atoms = atoms.copy()
    atoms.calc = calc
    return float(atoms.get_potential_energy()), np.asarray(atoms.get_forces(apply_constraint=False))


def counts(atoms, numbers):
    return np.array([(atoms.numbers == z).sum() for z in numbers], dtype=float)


# ---- 1. foundation model
base_path = Path(JOB["foundation_model"])
assert base_path.is_file(), f"no foundation model at {base_path}"
SUMMARY["foundation_model"] = {"path": str(base_path), "bytes": base_path.stat().st_size, "head": JOB["foundation_head"]}
if JOB.get("foundation_model_sha256"):
    SUMMARY["foundation_model"]["sha256"] = sha256(base_path)
    assert SUMMARY["foundation_model"]["sha256"] == JOB["foundation_model_sha256"], "foundation model changed"
dump()

frames = read(ROOT / "train.extxyz", ":")
numbers = sorted({int(z) for atoms in frames for z in atoms.numbers})
replay = bool(JOB["replay"])

# ---- 2. foundation predictions and replay labels
base_calc = calculator(base_path, JOB["foundation_head"])
base_pred = [predict(base_calc, atoms) for atoms in frames]
replay_labelled = {}
if replay:
    labels_dir = ROOT / "replay_labels"
    labels_dir.mkdir(exist_ok=False)
    for split in ("train", "valid"):
        labelled = []
        for atoms in read(ROOT / f"replay_{split}_geometries.extxyz", ":"):
            energy, forces = predict(base_calc, atoms)
            item = atoms.copy()
            item.calc = None
            item.info.update(REF_energy=energy, head="pt_head")
            item.arrays["REF_forces"] = forces
            labelled.append(item)
        write(labels_dir / f"{split}.extxyz", labelled, format="extxyz")
        replay_labelled[split] = read(labels_dir / f"{split}.extxyz", ":")
    SUMMARY["replay"] = {split: len(items) for split, items in replay_labelled.items()}

# ---- 3. E0s of the fine-tuned head
matrix = np.array([counts(atoms, numbers) for atoms in frames])
residuals = np.array([atoms.info["REF_energy"] - energy for atoms, (energy, _) in zip(frames, base_pred)])
assert np.linalg.matrix_rank(matrix) == len(numbers), "the training compositions do not separate the elements"
shifts = np.linalg.lstsq(matrix, residuals, rcond=None)[0]
base_model = torch.load(base_path, map_location="cpu", weights_only=False)
heads = list(getattr(base_model, "heads", ["Default"]))
head_index = heads.index(JOB["foundation_head"]) if JOB["foundation_head"] in heads else 0
atomic = base_model.atomic_energies_fn.atomic_energies
z_table = base_model.atomic_numbers.tolist()
base_e0 = np.array([float(atomic[head_index, z_table.index(z)] if atomic.dim() == 2 else atomic[z_table.index(z)])
                    for z in numbers])
scale = base_model.scale_shift.scale
shift = base_model.scale_shift.shift
std = float(scale[head_index] if scale.dim() and scale.numel() > 1 else scale.reshape(-1)[0])
mean = float(shift[head_index] if shift.dim() and shift.numel() > 1 else shift.reshape(-1)[0])
e0s = base_e0 + shifts
SUMMARY["e0s"] = {
    "numbers": numbers, "foundation_E0s_eV": base_e0.tolist(), "shifts_eV": shifts.tolist(), "E0s_eV": e0s.tolist(),
    "residual_rms_after_shift_eV": float(np.sqrt(np.mean((residuals - matrix @ shifts) ** 2))),
    "per_frame_DFT_minus_foundation_eV": {atoms.info.get("label", str(i)): float(r)
                                          for i, (atoms, r) in enumerate(zip(frames, residuals))},
    "std": std, "mean": mean,
}
del base_calc, base_model
dump()


# ---- 4.-6. one fine-tune per run
def check_configure(original, expected, expected_shift, first_frame, first_replay, record):
    def checked(*positional, **keyword):
        model, outputs = original(*positional, **keyword)
        model.to(DEVICE)
        model.eval()

        def run(atoms, head):
            datum = data.AtomicData.from_config(
                data.config_from_atoms(atoms, head_name=head), z_table=AtomicNumberTable(model.atomic_numbers.tolist()),
                cutoff=float(model.r_max), heads=model.heads)
            batch = next(iter(torch_geometric.dataloader.DataLoader(dataset=[datum], batch_size=1))).to(DEVICE)
            out = model(batch.to_dict(), training=False, compute_force=True, compute_stress=False)
            return float(out["energy"].detach().cpu().item()), out["forces"].detach().cpu().numpy()

        target = next(name for name in model.heads if name != "pt_head")
        energy, forces = run(first_frame, target)
        record.update(heads=list(model.heads), target_head=target,
                      initial_energy_error_eV_atom=abs(energy - expected[0] - expected_shift) / len(first_frame),
                      initial_force_error_eV_A=float(np.max(np.abs(forces - expected[1]))))
        dump()
        assert record["initial_energy_error_eV_atom"] < 1e-8 and record["initial_force_error_eV_A"] < 1e-6
        if first_replay is not None:
            energy, forces = run(first_replay, "pt_head")
            record.update(
                replay_initial_energy_error_eV_atom=abs(energy - first_replay.info["REF_energy"]) / len(first_replay),
                replay_initial_force_error_eV_A=float(np.max(np.abs(forces - first_replay.arrays["REF_forces"]))))
            dump()
            assert record["replay_initial_energy_error_eV_atom"] < 1e-8
            assert record["replay_initial_force_error_eV_A"] < 1e-6
        record["initialization_passed"] = True
        dump()
        return model, outputs

    return checked


STATE_FILE = "psteros_state.pt"


class ChunkDone(Exception):
    """Raised after the last epoch of a chunk (the state is saved) to skip MACE's export and evaluation."""


def check_train(original, record, passes, folder=None, stop=None):
    def checked(**keyword):
        loaders = keyword["valid_loaders"]
        assert set(loaders) == ({"Default", "pt_head"} if replay else {"Default"}), sorted(loaders)
        assert keyword["max_num_epochs"] == passes and keyword["patience"] == passes + 1
        record.update(validation_heads=list(loaders), training_examples_per_epoch=len(keyword["train_loader"].dataset))
        last = passes if stop is None else stop
        state_path = None if folder is None else folder / STATE_FILE
        start = 0
        if state_path is not None and state_path.exists():
            # on the CPU: the generator states must stay CPU byte tensors; the model, optimizer and EMA move what
            # they load to the device of their parameters
            state = torch.load(state_path, map_location="cpu", weights_only=False)
            keyword["model"].load_state_dict(state["model"])
            keyword["optimizer"].load_state_dict(state["optimizer"])
            keyword["lr_scheduler"].lr_scheduler.load_state_dict(state["lr_scheduler"])
            if keyword.get("ema") is not None:
                keyword["ema"].load_state_dict(state["ema"])
            keyword["train_loader"].generator.set_state(state["shuffle"])
            torch.set_rng_state(state["torch_rng"])
            if DEVICE == "cuda" and state.get("cuda_rng") is not None:
                torch.cuda.set_rng_state(state["cuda_rng"])
            np.random.set_state(state["numpy_rng"])
            random.setstate(state["python_rng"])
            for name, parameter in keyword["model"].named_parameters():
                if name in state["outside_gradients"]:
                    parameter.grad = state["outside_gradients"][name].to(parameter.device)
            start = int(state["epoch"])
            # MACE steps the scheduler before every pass but the first of a run (epoch > start_epoch): take the
            # step of pass ``start`` that an uninterrupted run takes (a no-op for a constant learning rate)
            keyword["lr_scheduler"].step()
        assert start < last <= passes, f"nothing to train: epochs {start} to {last} of {passes}"
        keyword.update(start_epoch=start, max_num_epochs=last)
        record.update(start_epoch=start, stop_epoch=last)
        dump()
        result = original(**keyword)
        if last < passes:
            saved = {
                "epoch": last, "model": keyword["model"].state_dict(), "optimizer": keyword["optimizer"].state_dict(),
                "lr_scheduler": keyword["lr_scheduler"].lr_scheduler.state_dict(),
                "ema": keyword["ema"].state_dict() if keyword.get("ema") is not None else None,
                "shuffle": keyword["train_loader"].generator.get_state(),
                "torch_rng": torch.get_rng_state(),
                "cuda_rng": torch.cuda.get_rng_state() if DEVICE == "cuda" else None,
                "numpy_rng": np.random.get_state(), "python_rng": random.getstate(),
                # MACE clips the gradient norm over every model parameter, while optimizer.zero_grad clears only the
                # optimised ones: the gradient of a trainable parameter outside the optimizer (the Bessel weights of
                # the radial embedding) keeps accumulating and enters the clipping of every later step
                "outside_gradients": {
                    name: parameter.grad.detach().clone()
                    for name, parameter in keyword["model"].named_parameters()
                    if parameter.grad is not None and id(parameter) not in {
                        id(member) for group in keyword["optimizer"].param_groups for member in group["params"]}
                },
            }
            temporary = state_path.with_suffix(".tmp")
            torch.save(saved, temporary)
            os.replace(temporary, state_path)
            record["state_saved_after_epoch"] = last
            record["outside_gradient_norms"] = {name: float(grad.norm()) for name, grad in saved["outside_gradients"].items()}
            dump()
            raise ChunkDone()
        return result

    return checked


def last_only_save(original, record, last):
    def save(self, state, epochs, keep_last=False):
        record.setdefault("checkpoint_save_calls", []).append(int(epochs))
        if epochs == last:
            original(self, state=state, epochs=epochs, keep_last=False)
            record["checkpoint_written"] = int(epochs)

    return save


def force_rmse(calc, items):
    rows = {}
    for i, atoms in enumerate(items):
        energy, forces = predict(calc, atoms)
        delta = forces - atoms.arrays["REF_forces"]
        rows[atoms.info.get("label", str(atoms.info.get("replay_source_frame", i)))] = {
            "natoms": len(atoms), "energy_error_eV": energy - atoms.info["REF_energy"],
            "force_rmse_meV_A": float(1000 * np.sqrt((delta ** 2).mean()))}
    return rows


originals = (run_train.configure_model, run_train.tools.train, tools.CheckpointHandler.save)
BASE = Path(JOB["persistent_root"]) if JOB.get("persistent_root") else ROOT
(BASE / "models").mkdir(parents=True, exist_ok=True)
SUMMARY["status"] = "training"


def keep_light_copies(folder: Path, name: str) -> None:
    """With a persistent root, copy the small results of a run into the job folder (what the job returns)."""

    if BASE == ROOT:
        return
    target = ROOT / "runs" / name
    for pattern in ("evaluation.json", "config.json", "E0s.json", "results/logs/*.log"):
        for path in folder.glob(pattern):
            destination = target / path.relative_to(folder)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(path, destination)


try:
    for spec in JOB["runs"]:
        name, passes = spec["name"], int(spec["options"]["max_num_epochs"])
        stop = int(spec.get("stop_epoch", passes))
        record = SUMMARY["runs"].setdefault(name, {"seed": spec["options"]["seed"], "status": "starting"})
        folder = BASE / "runs" / name
        e0s_text = json.dumps(dict(zip(map(str, numbers), e0s.tolist())), indent=2) + "\n"
        options = dict(spec["options"], device=DEVICE, mean=mean, std=std)
        config_text = json.dumps(options, indent=2) + "\n"
        if (folder / STATE_FILE).exists():
            # a later chunk: the recipe and the data must be those the run started with
            assert (folder / "config.json").read_text() == config_text, f"{name}: the MACE options changed"
            # the E0s come from GPU predictions, whose last bits can change between jobs; the run keeps its own
            previous = json.loads((folder / "E0s.json").read_text())
            assert all(abs(previous[key] - value) < 1e-6 for key, value in json.loads(e0s_text).items()), \
                f"{name}: the E0s changed"
            assert (folder / "train.extxyz").read_bytes() == (ROOT / "train.extxyz").read_bytes(), f"{name}: data"
        else:
            if folder.exists() and BASE != ROOT:
                # a first chunk that ended before saving a state (a walltime): keep it aside and start again
                folder.rename(folder.with_name(f"{name}.incomplete.{int(time.time())}"))
            folder.mkdir(parents=True, exist_ok=False)
            shutil.copy(ROOT / "train.extxyz", folder / "train.extxyz")
            if replay:
                shutil.copytree(ROOT / "replay_labels", folder / "replay_labels")
            (folder / "E0s.json").write_text(e0s_text)
            (folder / "config.json").write_text(config_text)
        os.chdir(folder)
        expected_shift = float(counts(frames[0], numbers) @ shifts)
        run_train.configure_model = check_configure(
            originals[0], base_pred[0], expected_shift, frames[0],
            replay_labelled["train"][0] if replay else None, record)
        run_train.tools.train = check_train(originals[1], record, passes, folder, stop)
        tools.CheckpointHandler.save = last_only_save(originals[2], record, passes - 1)
        start = time.monotonic()
        if DEVICE == "cuda":
            torch.cuda.reset_peak_memory_stats()
        chunk_done = False
        try:
            run_train.run(tools.build_default_arg_parser().parse_args(["--config", "config.json"]))
        except ChunkDone:
            chunk_done = True
        finally:
            run_train.configure_model, run_train.tools.train, tools.CheckpointHandler.save = originals
            os.chdir(ROOT)
            record["elapsed_seconds"] = time.monotonic() - start
            if DEVICE == "cuda":
                record["peak_allocated_GiB"] = torch.cuda.max_memory_allocated() / 2**30
            dump()
        assert record.get("initialization_passed"), "the initialisation check did not run"
        assert record["checkpoint_save_calls"] == list(range(record["start_epoch"], record["stop_epoch"])), \
            record["checkpoint_save_calls"]
        if chunk_done:
            record["status"] = f"chunk done ({record['stop_epoch']} of {passes} passes)"
            keep_light_copies(folder, name)
            dump()
            continue
        logs = [p for p in (folder / "results/logs").glob("*.log") if not p.name.endswith("_debug.log")]
        loaded = re.findall(r"Loaded Stage one model from epoch (\d+)", logs[0].read_text()) if len(logs) == 1 else []
        assert loaded == [str(passes - 1)], f"the exported model is not the last epoch: {loaded}"
        exported = folder / "results/models" / f"{name}.model"
        if not exported.is_file():
            exported = folder / "results/checkpoints" / f"{name}.model"
        assert exported.is_file(), "no exported model"
        final = BASE / "models" / f"{name}.model"
        shutil.copy(exported, final)
        record.update(model=str(final.relative_to(BASE)), model_path=str(final), model_sha256=sha256(final),
                      status="trained")
        dump()
        evaluation = {"train": force_rmse(calculator(final, "Default"), frames)}
        if replay:
            for head in ("Default", "pt_head"):
                evaluation[f"replay_valid_{head}_vs_foundation"] = force_rmse(calculator(final, head),
                                                                             replay_labelled["valid"])
        (folder / "evaluation.json").write_text(json.dumps(evaluation, indent=2) + "\n")
        keep_light_copies(folder, name)
        record["status"] = "completed"
        dump()
    SUMMARY["status"] = "completed"
except Exception:
    SUMMARY.update(status="failed", traceback=traceback.format_exc())
    raise
finally:
    dump()
