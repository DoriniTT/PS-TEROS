"""psteros MACE vibrations job: runs on the GPU computer and needs only ase, numpy, torch and mace-torch.

``job.json`` (written by ``psteros.mlip.prepare_vibrations_job``) names the model (``{"path": ..., "head": ...}``), the
settings (``device``, ``default_dtype``, ``displacement_angstrom``) and the structures
(``[{"label": ..., "file": "<label>.extxyz", "sites": [zero-based indices to displace]}]``). Every listed site is moved
by ``+d`` and ``-d`` along x, y and z, as the displaced DFT calculations of ``psteros.build_vibrations_workgraph`` are,
and the model's forces on every atom are kept (no atom is held fixed). Results: ``vibrations/<label>.npz`` with
``forces_plus`` and ``forces_minus`` (eV/A, shape ``(3 * sites, atoms, 3)``, row ``3 k + a`` for the ``k``-th site moved
along axis ``a``) and ``results.json`` (largest force at the undisplaced geometry in eV/A, number of force evaluations).
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

import numpy as np
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
RESULTS: dict = {"status": "running", "settings": {k: v for k, v in JOB.items() if k != "structures"}, "structures": {}}


def dump() -> None:
    (ROOT / "results.json").write_text(json.dumps(RESULTS, indent=2, allow_nan=False) + "\n")


try:
    calc = MACECalculator(model_paths=str(JOB["model"]["path"]), device=JOB.get("device", "cuda"),
                          default_dtype=JOB.get("default_dtype", "float64"), head=JOB["model"]["head"])
    step = float(JOB["displacement_angstrom"])
    out = ROOT / "vibrations"
    out.mkdir(exist_ok=False)
    for item in JOB["structures"]:
        start = time.monotonic()
        atoms = read(ROOT / item["file"])
        atoms.pbc = True
        atoms.calc = calc
        reference = atoms.get_positions()
        residual = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        plus, minus = [], []
        for site in (int(index) for index in item["sites"]):
            for axis in range(3):
                for sign, forces in ((1.0, plus), (-1.0, minus)):
                    positions = reference.copy()
                    positions[site, axis] += sign * step
                    atoms.set_positions(positions)
                    forces.append(np.array(atoms.get_forces()))
        atoms.set_positions(reference)
        np.savez_compressed(out / f"{item['label']}.npz", forces_plus=np.array(plus), forces_minus=np.array(minus))
        RESULTS["structures"][item["label"]] = {
            "natoms": len(atoms), "displaced_sites": len(item["sites"]), "force_evaluations": 2 * len(plus) + 1,
            "max_force_at_start_ev_per_angstrom": residual, "seconds": time.monotonic() - start,
        }
        dump()
    RESULTS["status"] = "completed"
except Exception:
    RESULTS.update(status="failed", traceback=traceback.format_exc())
    raise
finally:
    dump()
