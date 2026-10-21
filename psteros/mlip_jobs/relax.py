"""psteros relaxation job: runs on the GPU computer and needs only ase, torch and mace-torch (or fairchem-core, or nequip).

``job.json`` (written by psteros) names the models (``{name: {"path": ..., "head": ...}}``, with ``"family": "uma"``
for a UMA model, whose ``head`` is its task), the structures
(``[{"label": ..., "file": "<label>.extxyz", "fixed": [zero-based indices]}]``) and the settings (``device``,
``default_dtype``, ``optimizer``, ``fmax_ev_per_angstrom``, ``max_steps``). Every model relaxes every structure
(positions only, cell fixed, fixed atoms kept in place), as ``psteros.mlip.relax_with_calculator`` does. Results:
``relaxed/<model>/<label>.extxyz`` (final geometry, energy and forces) and ``results.json`` (energy at the start and
at the end in eV, largest force on a free atom in eV/A, steps, converged).
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS, FIRE, LBFGS

ROOT = Path.cwd()
JOB = json.loads((ROOT / "job.json").read_text())
OPTIMIZERS = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}
RESULTS: dict = {"status": "running", "settings": {k: v for k, v in JOB.items() if k not in ("models", "structures")},
                 "models": {}}


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
    for model_name, model in JOB["models"].items():
        calc = calculator(model)
        out = ROOT / "relaxed" / model_name
        out.mkdir(parents=True, exist_ok=False)
        rows = RESULTS["models"].setdefault(model_name, {"path": str(model["path"]), "head": model["head"],
                                                         "structures": {}})["structures"]
        for item in JOB["structures"]:
            atoms = read(ROOT / item["file"])
            atoms.pbc = True
            fixed = sorted(int(index) for index in item.get("fixed", []))
            if fixed:
                atoms.set_constraint(FixAtoms(indices=fixed))
            atoms.calc = calc
            start = time.monotonic()
            initial = float(atoms.get_potential_energy())
            optimizer = OPTIMIZERS[JOB.get("optimizer", "FIRE")](atoms, logfile=None)
            optimizer.run(fmax=float(JOB["fmax_ev_per_angstrom"]), steps=int(JOB["max_steps"]))
            forces = atoms.get_forces()
            largest = float(np.linalg.norm(forces, axis=1).max())
            energy = float(atoms.get_potential_energy())
            final = atoms.copy()
            final.calc = None
            final.info.update(label=item["label"], model=model_name, energy_ev=energy)
            final.arrays["model_forces"] = atoms.get_forces(apply_constraint=False)
            write(out / f"{item['label']}.extxyz", final, format="extxyz")
            rows[item["label"]] = {
                "natoms": len(atoms), "initial_energy_ev": initial, "energy_ev": energy,
                "max_force_ev_per_angstrom": largest, "steps": int(optimizer.get_number_of_steps()),
                "converged": largest < float(JOB["fmax_ev_per_angstrom"]), "fixed": fixed,
                "seconds": time.monotonic() - start,
            }
            dump()
    RESULTS["status"] = "completed"
except Exception:
    RESULTS.update(status="failed", traceback=traceback.format_exc())
    raise
finally:
    dump()
