"""Step C6 of docs/plans/MACE_SURFACE_PLAN.md: one CP2K static on each model-relaxed slab (route b).

Models: M0 (MH-1 zero-shot) and the deployment seed (20260920) of M1, M2 and M3, read from the relaxation jobs of
the chunked study (``c4_chunks.py``). Only relaxations that reached 0.01 eV/A are computed; the others are listed.
The CP2K recipe is the one of the labels (``c2_cp2k.recipe``, ``settings.json``).

Usage (PYTHONPATH=<worktree>):
    python c6_route_b.py submit [--dry]                                  # M0 and M1-M3 (c6_route_b.json)
    python c6_route_b.py submit --models m2r_seed20260920 --record c6_route_b_m2r.json
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import psteros

import c2_cp2k

HERE = Path(__file__).resolve().parent
MODELS = ("m0", "m1_seed20260920", "m2_seed20260920", "m3_seed20260920")


def relaxed(model: str, label: str):
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor

    import c4_chunks

    atoms = read(io.StringIO(c4_chunks.relaxed_text(model, label)), format="extxyz")
    return AseAtomsAdaptor.get_structure(atoms)


def structures(models=MODELS) -> tuple[dict, dict]:
    import c4_chunks

    results = c4_chunks.relax_results()  # the main study ran in chunks on testegpu (c4_chunks.py)
    chosen, skipped = {}, {}
    for model in models:
        for label in psteros.RUTILE_110_TERMINATIONS:
            row = results[model]["structures"][label]
            if row["converged"]:
                structure = relaxed(model, label)
                chosen[f"{model}_{label}"] = psteros.structures.Structure(
                    structure.lattice, structure.species, structure.frac_coords)
            else:
                skipped[f"{model}_{label}"] = row["max_force_ev_per_angstrom"]
    return chosen, skipped


def submit(dry: bool, models=MODELS, record: str = "c6_route_b.json") -> None:
    chosen, skipped = structures(models)
    graph = psteros.build_cp2k_surface_workgraph(
        chosen, c2_cp2k.config("ENERGY_FORCE", "ruo2_c6_route_b", 2), submit=not dry)  # cluster.py: <= 2 jobs
    print(len(chosen), "statics; not converged (left out):", skipped)
    if not dry:
        (HERE / record).write_text(json.dumps(
            {"workgraph_pk": graph.pk, "labels": list(chosen), "not_converged": skipped}, indent=2) + "\n")
        print("WorkGraph", graph.pk)


if __name__ == "__main__":
    from aiida import load_profile

    import cluster

    load_profile(c2_cp2k.PROFILE)
    models = tuple(sys.argv[sys.argv.index("--models") + 1].split(",")) if "--models" in sys.argv else MODELS
    record = sys.argv[sys.argv.index("--record") + 1] if "--record" in sys.argv else "c6_route_b.json"
    cluster.with_retry(submit, "--dry" in sys.argv, models, record)
