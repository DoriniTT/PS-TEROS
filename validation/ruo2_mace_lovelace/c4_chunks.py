"""C4/C5 in chunks on Lovelace ``testegpu`` (user, 2026-10-18: "Split into testegpu chunks").

The Bohr A100 was held by a 168 h job with another queued ahead. ``testegpu`` gives one A100 job of at most 30 min
per user (and two queued). So every fine-tune of the main study (M1, M2, M3 x 3 seeds, 100 passes each, the same
recipe as ``c4_gpu.py main``) runs as a sequence of jobs, each one chunk of one run, with
``psteros/mlip_jobs/finetune.py``'s ``stop_epoch`` and ``persistent_root``: the training state and the models stay
in ``/work/dorinitt/psteros_mace_ruo2`` on Lovelace between jobs, and a chunked run follows the same path as one
uninterrupted run (checked on CPU to machine precision). Chunks are sized from the measured time per pass to about
20 min of training. When every run is done, one relaxation job per model (M0 and the nine fine-tunes) relaxes the
nine structures.

Usage (PYTHONPATH=<worktree>):
    python c4_chunks.py prepare       # base folders and plan (local/jobs/chunks)
    python c4_chunks.py add-variant   # M2L, the extra variant (below), into the plan
    python c4_chunks.py run           # the orchestrator loop; resumable from local/jobs/chunks/plan.json
    python c4_chunks.py status

Extra variant M2L (user, 2026-10-18, after the first M2 model was not converged at 100 passes; fixed before any of
its results; reported, not gated): the M2 recipe with 500 passes and a learning rate decaying from 1e-3 to 1e-5
(MACE ExponentialLR, gamma = 0.01 ** (1 / 500) = 0.990832 per pass), deployment seed 20260920 only. Priority: the
main-study fine-tunes first, then the relaxation of every model whose training is done, then M2L.
"""

from __future__ import annotations

import json
import math
import shutil
import sys
import time
from pathlib import Path

import psteros
from psteros import mlip

import c4_gpu
import cluster

HERE = Path(__file__).resolve().parent
JOBS = HERE / "local" / "jobs" / "chunks"
PLAN = JOBS / "plan.json"
LOG = JOBS / "orchestrator.log"
REMOTE_ROOT = "/work/dorinitt/psteros_mace_ruo2"
COMPUTER = c4_gpu.COMPUTERS["lovelace"]
WALLTIME = 1800
TRAIN_BUDGET_SECONDS = 1200
FIRST_GUESS_SECONDS_PER_PASS = {"m1": 6.0, "m2": 52.0, "m3": 64.0}
FAMILIES = ("m2", "m1", "m3")
EXTRA = {"m2l": {"epochs": 500, "lr_scheduler_gamma": 0.990832, "seeds": (20260920,), "replay": True,
                 "frames": "plain", "rate_from": "m2"}}
MAX_IN_FLIGHT = 2
MAX_FAILURES = 3


def log(message: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
    print(line, flush=True)
    with LOG.open("a") as stream:
        stream.write(line + "\n")


def prepare() -> None:
    if JOBS.exists():
        raise FileExistsError(f"{JOBS} exists")
    plain = cluster.with_retry(c4_gpu.frames, False)
    rattled = cluster.with_retry(c4_gpu.frames, True)
    structures = cluster.with_retry(c4_gpu.relaxation_structures)
    train, valid, record = mlip.select_replay_structures(c4_gpu.REPLAY_FILE, elements=["Ru", "O"])
    base = dict(foundation_model=COMPUTER["foundation"], foundation_model_sha256=c4_gpu.FOUNDATION_SHA256,
                replay_file=c4_gpu.REPLAY_FILE)
    configs = {"m1": (mlip.MaceFineTuneConfig(**{**base, "replay_file": None}, replay=False), plain),
               "m2": (mlip.MaceFineTuneConfig(**base), plain), "m3": (mlip.MaceFineTuneConfig(**base), rattled)}
    runs = {}
    for family, (config, items) in configs.items():
        replay = dict(replay_train=train, replay_valid=valid) if config.replay else {}
        mlip.prepare_finetune_job(JOBS / f"ft_{family}", items, config, name=f"ruo2_{family}", **replay)
        for seed in config.seeds:
            runs[f"ruo2_{family}_seed{seed}"] = {"family": family, "seed": seed, "passes": config.epochs,
                                                  "epochs_done": 0, "failures": 0, "chunks": []}
    models = {"m0": (COMPUTER["foundation"], "omat_pbe")}
    for name, run in runs.items():
        models[f"{run['family']}_seed{run['seed']}"] = (f"{REMOTE_ROOT}/models/{name}.model", "Default")
    fixed = {label: psteros.central_sites(structures[label]) for label in psteros.RUTILE_110_TERMINATIONS}
    for model, (path, head) in models.items():
        mlip.prepare_relax_job(JOBS / f"relax_{model}", structures, {model: (path, head)}, fixed_sites=fixed)
    plan = {"remote_root": REMOTE_ROOT, "seconds_per_pass": dict(FIRST_GUESS_SECONDS_PER_PASS), "runs": runs,
            "relax": {model: {"state": "waiting", "pk": None, "failures": 0} for model in models},
            "replay_selection": record}
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")
    log(f"prepared {len(runs)} runs and {len(models)} relaxation jobs")


def submit_folder(folder: Path, label: str) -> int:
    from aiida import engine, orm
    from aiida_shell import ShellJob

    plan = json.loads((folder / "driver.json").read_text())
    nodes = {"driver": orm.SinglefileData(file=str(folder / "gpu_driver.py")),
             "driver_json": orm.SinglefileData(file=str(folder / "driver.json"))}
    filenames = {"driver": "gpu_driver.py", "driver_json": "driver.json"}
    for name in list(plan["finetune"]) + ([plan["relax"]] if plan["relax"] else []):
        nodes[name] = orm.FolderData(tree=str(folder / name))
        filenames[name] = name
    code = orm.load_code(COMPUTER["code"])
    options = dict(c4_gpu.QUEUES["testegpu"], max_wallclock_seconds=WALLTIME, withmpi=False, redirect_stderr=True,
                   import_sys_environment=False,
                   environment_variables={"OMP_NUM_THREADS": "8", "MKL_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1",
                                          "PYTHONNOUSERSITE": "1", "PYTORCH_CUDA_ALLOC_CONF": "backend:cudaMallocAsync"})
    node = engine.submit(ShellJob, code=code, nodes=nodes, filenames=filenames, arguments=["-u", "{driver}"],
                         outputs=["report"], metadata={"label": label, "options": options})
    return node.pk


def chunk_folder(name: str, run: dict, start: int, stop: int) -> Path:
    folder = JOBS / "jobs" / f"{name}_{start:03d}-{stop:03d}_{int(time.time())}"
    shutil.copytree(JOBS / f"ft_{run['family']}", folder / "ft")
    job = json.loads((folder / "ft" / "job.json").read_text())
    spec = next(item for item in job["runs"] if item["name"] == name)
    if stop < run["passes"]:
        spec["stop_epoch"] = stop
    job["runs"] = [spec]
    job["persistent_root"] = REMOTE_ROOT
    (folder / "ft" / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    (folder / "driver.json").write_text(json.dumps({"gpus": 1, "finetune": ["ft"], "relax": None}) + "\n")
    shutil.copy(HERE / "gpu_driver.py", folder / "gpu_driver.py")
    return folder


def relax_folder(model: str) -> Path:
    folder = JOBS / "jobs" / f"relax_{model}_{int(time.time())}"
    shutil.copytree(JOBS / f"relax_{model}", folder / "relax")
    (folder / "driver.json").write_text(json.dumps({"gpus": 1, "finetune": [], "relax": "relax"}) + "\n")
    shutil.copy(HERE / "gpu_driver.py", folder / "gpu_driver.py")
    return folder


def report_json(pk: int, path: str):
    from aiida import orm

    node = orm.load_node(pk)
    if "report" not in node.outputs:
        return None
    try:
        return json.loads(node.outputs.report.base.repository.get_object_content(path))
    except (FileNotFoundError, OSError):
        return None


def job_state(pk: int) -> str:
    from aiida import orm

    node = orm.load_node(pk)
    if not node.is_terminated:
        if node.paused:
            log(f"PROBLEM job {pk} is paused (refused qsub or lost connection): play it with verdi process play")
        return "active"
    return "ok" if node.exit_status == 0 else f"failed[{node.exit_status}]"


def step(plan: dict) -> bool:
    """One pass of the loop; returns True when everything is done."""

    runs, relax = plan["runs"], plan["relax"]
    in_flight = 0
    for name, run in runs.items():
        chunk = run["chunks"][-1] if run["chunks"] else None
        if chunk and chunk["state"] == "active":
            state = job_state(chunk["pk"])
            if state == "active":
                in_flight += 1
                continue
            summary = report_json(chunk["pk"], "ft/summary.json")
            record = (summary or {}).get("runs", {}).get(name, {})
            done = record.get("status", "").startswith(("chunk done", "completed"))
            if done:
                chunk.update(state="done", job=state, elapsed_seconds=record.get("elapsed_seconds"),
                             outside_gradient_norms=record.get("outside_gradient_norms"))
                run["epochs_done"] = chunk["stop"]
                passes = chunk["stop"] - chunk["start"]
                if passes > 0 and record.get("elapsed_seconds"):
                    family = run["family"]
                    measured = record["elapsed_seconds"] / passes
                    plan["seconds_per_pass"][family] = max(plan["seconds_per_pass"][family] * 0.5, measured)
                log(f"{name}: passes {chunk['start']}-{chunk['stop']} done in {record.get('elapsed_seconds', 0):.0f} s"
                    f" (job {chunk['pk']}, {state})")
            else:
                chunk["state"] = "failed"
                run["failures"] += 1
                plan["seconds_per_pass"][run["family"]] *= 1.5
                log(f"PROBLEM {name}: chunk {chunk['start']}-{chunk['stop']} failed (job {chunk['pk']}, {state}); "
                    f"failure {run['failures']} of {MAX_FAILURES}")
                if run["failures"] >= MAX_FAILURES:
                    raise RuntimeError(f"{name} failed {MAX_FAILURES} times")
    for model, entry in relax.items():
        if entry["state"] == "active":
            state = job_state(entry["pk"])
            if state == "active":
                in_flight += 1
                continue
            results = report_json(entry["pk"], "relax/results.json")
            if results and results.get("status") == "completed":
                entry["state"] = "done"
                log(f"relax {model} done (job {entry['pk']})")
            else:
                entry.update(state="waiting", failures=entry["failures"] + 1)
                log(f"PROBLEM relax {model} failed (job {entry['pk']}, {state})")
                if entry["failures"] >= MAX_FAILURES:
                    raise RuntimeError(f"relax {model} failed {MAX_FAILURES} times")
    busy = {name for name, run in runs.items() if run["chunks"] and run["chunks"][-1]["state"] == "active"}

    def trained(model: str) -> bool:
        if model == "m0":
            return True
        family, seed = model.split("_seed")
        run = runs[f"ruo2_{family}_seed{seed}"]
        return run["epochs_done"] >= run["passes"]

    def next_task():
        """Main fine-tunes, then relaxations of trained models (main first), then the extra variants."""

        families = list(FAMILIES) + [family for family in plan.get("extra_families", []) if family not in FAMILIES]
        main = set(FAMILIES)
        for group in ("main", "relax", "extra"):
            if group == "relax":
                ready = [model for model, entry in relax.items() if entry["state"] == "waiting" and trained(model)]
                ready.sort(key=lambda model: model.split("_seed")[0] not in main and model != "m0")
                if ready:
                    return ("relax", ready[0])
                continue
            for family in families:
                if (family in main) != (group == "main"):
                    continue
                for name, run in runs.items():
                    if run["family"] == family and name not in busy and run["epochs_done"] < run["passes"]:
                        return ("chunk", name)
        return None

    while in_flight < MAX_IN_FLIGHT:
        task = next_task()
        if task is None:
            break
        kind, key = task
        if kind == "chunk":
            run = runs[key]
            per_pass = plan["seconds_per_pass"][run["family"]]
            size = max(1, math.floor(TRAIN_BUDGET_SECONDS / per_pass))
            start, stop = run["epochs_done"], min(run["passes"], run["epochs_done"] + size)
            folder = chunk_folder(key, run, start, stop)
            pk = cluster.with_retry(submit_folder, folder, f"ruo2_{key}_{start}-{stop}")
            run["chunks"].append({"start": start, "stop": stop, "pk": pk, "state": "active", "folder": folder.name})
            busy.add(key)
            log(f"submitted {key} passes {start}-{stop} as job {pk}")
        else:
            pk = cluster.with_retry(submit_folder, relax_folder(key), f"ruo2_relax_{key}")
            relax[key].update(state="active", pk=pk)
            log(f"submitted relaxations with {key} as job {pk}")
        in_flight += 1
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")
    return all(run["epochs_done"] >= run["passes"] for run in runs.values()) and all(
        entry["state"] == "done" for entry in relax.values())


def add_variant() -> None:
    plan = json.loads(PLAN.read_text())
    structures = cluster.with_retry(c4_gpu.relaxation_structures)
    plain = cluster.with_retry(c4_gpu.frames, False)
    train, valid, _ = mlip.select_replay_structures(c4_gpu.REPLAY_FILE, elements=["Ru", "O"])
    fixed = {label: psteros.central_sites(structures[label]) for label in psteros.RUTILE_110_TERMINATIONS}
    for family, spec in EXTRA.items():
        if family in plan.get("extra_families", []):
            raise ValueError(f"{family} is already in the plan")
        config = mlip.MaceFineTuneConfig(
            foundation_model=COMPUTER["foundation"], foundation_model_sha256=c4_gpu.FOUNDATION_SHA256,
            replay_file=c4_gpu.REPLAY_FILE, epochs=spec["epochs"], seeds=spec["seeds"],
            lr_scheduler_gamma=spec["lr_scheduler_gamma"])
        mlip.prepare_finetune_job(JOBS / f"ft_{family}", plain, config, name=f"ruo2_{family}",
                                  replay_train=train, replay_valid=valid)
        for seed in config.seeds:
            name = f"ruo2_{family}_seed{seed}"
            plan["runs"][name] = {"family": family, "seed": seed, "passes": config.epochs, "epochs_done": 0,
                                  "failures": 0, "chunks": []}
            model = f"{family}_seed{seed}"
            mlip.prepare_relax_job(JOBS / f"relax_{model}", structures,
                                   {model: (f"{REMOTE_ROOT}/models/{name}.model", "Default")}, fixed_sites=fixed)
            plan["relax"][model] = {"state": "waiting", "pk": None, "failures": 0}
        plan["seconds_per_pass"][family] = plan["seconds_per_pass"][spec["rate_from"]]
        plan.setdefault("extra_families", []).append(family)
        log(f"added the extra variant {family}: {spec}")
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")


def run_loop() -> None:
    while True:
        plan = json.loads(PLAN.read_text())
        if cluster.with_retry(step, plan):
            log("ALL DONE")
            return
        time.sleep(120)


def status() -> None:
    plan = json.loads(PLAN.read_text())
    for name, run in plan["runs"].items():
        print(f"{name:24s} {run['epochs_done']:3d}/{run['passes']} chunks={len(run['chunks'])} failures={run['failures']}")
    print({model: entry["state"] for model, entry in plan["relax"].items()})
    print("seconds per pass:", plan["seconds_per_pass"])


if __name__ == "__main__":
    from aiida import load_profile

    load_profile("psteros_mace_ruo2")
    {"prepare": prepare, "add-variant": add_variant, "run": run_loop, "status": status}[sys.argv[1]]()


# ---------------------------------------------------------------- readers for c6_route_b.py and c7_analysis.py

def relax_results() -> dict:
    """``{model: {"path", "head", "structures": {label: row}}}`` merged from the relaxation job of every model."""

    plan = json.loads(PLAN.read_text())
    merged = {}
    for model, entry in plan["relax"].items():
        if entry["state"] == "done":
            merged.update(report_json(entry["pk"], "relax/results.json")["models"])
    return merged


def relaxed_text(model: str, label: str) -> str:
    """The extended XYZ of ``label`` relaxed by ``model``."""

    from aiida import orm

    plan = json.loads(PLAN.read_text())
    node = orm.load_node(plan["relax"][model]["pk"])
    return node.outputs.report.base.repository.get_object_content(f"relax/relaxed/{model}/{label}.extxyz")
