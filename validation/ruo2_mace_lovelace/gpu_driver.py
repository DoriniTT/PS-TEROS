"""Run on the GPU node: the fine-tune folders on the visible GPUs in parallel, then the relaxation folder.

``driver.json``: ``{"gpus": 2, "finetune": ["ft_m1", "ft_m2"], "relax": "relax"}``. Fine-tune folder ``i`` runs on
GPU ``i % gpus`` (one after the other on the same GPU), each with ``python finetune.py`` (psteros.mlip_jobs); the
relaxation folder then runs ``python relax.py``, whose model paths point into the fine-tune folders. Needs only the
Python of the MACE environment. Exit status: the largest exit status of the steps; ``driver_status.json`` records each.

At the end, also after a failure, ``report/`` collects what the job retrieves: the status, every ``summary.json``,
log, ``evaluation.json``, ``config.json`` and ``E0s.json`` of the fine-tunes, and the relaxation results. Models and
checkpoints stay in the remote folder (their SHA-256 are in the summaries).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path.cwd()
PLAN = json.loads((ROOT / "driver.json").read_text())
STATUS: dict = {"started": time.time(), "finetune": {}, "relax": None}


def dump() -> None:
    (ROOT / "driver_status.json").write_text(json.dumps(STATUS, indent=2) + "\n")


def report() -> None:
    target = ROOT / "report"
    target.mkdir(exist_ok=True)
    patterns = ["driver_status.json", "finetune_exit_codes.txt"]
    for folder in PLAN.get("finetune", []):
        patterns += [f"{folder}/summary.json", f"{folder}/finetune.log", f"{folder}/runs/*/evaluation.json",
                     f"{folder}/runs/*/config.json", f"{folder}/runs/*/E0s.json", f"{folder}/runs/*/results/logs/*.log"]
    if PLAN.get("relax"):
        patterns += [f"{PLAN['relax']}/results.json", f"{PLAN['relax']}/relax.log", f"{PLAN['relax']}/relaxed/*/*.extxyz"]
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            if path.is_file() and not path.name.endswith("_debug.log"):
                destination = target / path.relative_to(ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(path, destination)


codes: dict[str, int] = {}
try:
    gpus = max(1, int(PLAN.get("gpus", 1)))
    lanes: list[list[str]] = [[] for _ in range(gpus)]
    for index, folder in enumerate(PLAN.get("finetune", [])):
        lanes[index % gpus].append(folder)
    processes = []
    for gpu, folders in enumerate(lanes):
        if not folders:
            continue
        steps = " ; ".join(
            f'(cd "{ROOT / folder}" && "{sys.executable}" -u finetune.py > finetune.log 2>&1; '
            f'echo "{folder} $?" >> "{ROOT}/finetune_exit_codes.txt")'
            for folder in folders
        )
        environment = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu))
        processes.append(subprocess.Popen(["bash", "-c", steps], env=environment))
    for process in processes:
        process.wait()
    if (ROOT / "finetune_exit_codes.txt").exists():
        for line in (ROOT / "finetune_exit_codes.txt").read_text().split("\n"):
            if line.strip():
                folder, code = line.split()
                codes[folder] = int(code)
    STATUS["finetune"] = codes
    dump()
    if PLAN.get("relax"):
        folder = ROOT / PLAN["relax"]
        with (folder / "relax.log").open("w") as log:
            STATUS["relax"] = subprocess.call([sys.executable, "-u", "relax.py"], cwd=folder, stdout=log,
                                              stderr=subprocess.STDOUT, env=dict(os.environ, CUDA_VISIBLE_DEVICES="0"))
    STATUS["finished"] = time.time()
    dump()
finally:
    dump()
    report()
missing = [folder for folder in PLAN.get("finetune", []) if folder not in codes]
sys.exit(max([*codes.values(), STATUS["relax"] or 0, 1 if missing else 0]))
