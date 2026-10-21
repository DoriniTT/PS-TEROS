"""psteros vibrations job: runs on the GPU computer and needs only ase, numpy, torch and mace-torch (or fairchem-core, or nequip).

``job.json`` (written by ``psteros.mlip.prepare_vibrations_job``) names the model (``{"path": ..., "head": ...}``,
with ``"family": "uma"`` for a UMA model, whose ``head`` is its task), the
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

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
RESULTS: dict = {"status": "running", "settings": {k: v for k, v in JOB.items() if k != "structures"}, "structures": {}}


def dump() -> None:
    (ROOT / "results.json").write_text(json.dumps(RESULTS, indent=2, allow_nan=False) + "\n")


def calculator(model: dict):
    """The ASE calculator of ``model``: MACE (``family`` left out or ``"mace"``), UMA (``"uma"``, ``head`` = task) or
    the NequIP framework (``"nequip"``: NequIP and Allegro packages, ``head`` unused)."""

    family = model.get("family", "mace")
    if family == "mace":
        from mace.calculators import MACECalculator

        return MACECalculator(model_paths=str(model["path"]), device=JOB.get("device", "cuda"),
                              default_dtype=JOB.get("default_dtype", "float64"), head=model["head"])
    if family == "uma":
        from fairchem.core import FAIRChemCalculator
        from fairchem.core.units.mlip_unit import load_predict_unit
        from fairchem.core.units.mlip_unit.api.inference import InferenceSettings

        # MoLE merged for the fixed composition (exact), no torch.compile (it needs a recent C++ compiler)
        settings = InferenceSettings(tf32=False, activation_checkpointing=False, merge_mole=True, compile=False,
                                     base_precision_dtype=JOB.get("default_dtype", "float64"))
        unit = load_predict_unit(str(model["path"]), inference_settings=settings, device=JOB.get("device", "cuda"))
        return FAIRChemCalculator(unit, task_name=model["head"])
    if family == "nequip":
        from nequip.integrations.ase import NequIPCalculator

        # a NequIP or Allegro package, run eagerly (nothing compiled); its type names are chemical symbols
        return NequIPCalculator._from_saved_model(str(model["path"]), device=JOB.get("device", "cuda"),
                                                  chemical_species_to_atom_type_map=True)
    raise ValueError(f"unknown model family {family!r}")


try:
    calc = calculator(JOB["model"])
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
