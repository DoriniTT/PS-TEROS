"""Step C0 of docs/plans/MACE_SURFACE_PLAN.md: smoke test of the CP2K code (``cluster.py``; first on Lovelace).

Three ENERGY_FORCE calculations of bulk RuO2 (6 atoms) with one O moved 0.05 A along x (``base``), and the same
O moved a further +-0.01 A (``plus``, ``minus``). ``check`` reports, per calculation, the MPI ranks and CP2K version
from the output, the energy, the parsed forces and the DOS; then:

* one CP2K run per job (the output lists the requested number of MPI processes, not 1 per copy);
* the energy, the force file and the DOS file are read by the psteros readers;
* bulk RuO2 is a metal for require_metallic (partly filled DOS bins);
* the energy is the one the forces derive from: -(E+ - E-) / (2 delta) against the printed force on that O
  (agreement within 0.01 eV/A is taken as consistent).

Usage (PYTHONPATH=<worktree>):
    python c0_smoke.py submit [--cutoff 800]   # builds and submits the graph, writes c0_graph_<cluster>_<cutoff>.json
    python c0_smoke.py check [--cutoff 800]    # after the graph has finished

The first run (500 Ry) agreed to 0.011 eV/A (1.2%); CP2K's total energy is the sum of its terms including the
electronic entropic energy, so the gap is expected to come from the real-space grid. The 800 Ry run checks that.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

import psteros

import cluster

HERE = Path(__file__).resolve().parent
PROFILE = "psteros_mace_ruo2"
CODE = cluster.CP2K_CODE
CUTOFF = int(sys.argv[sys.argv.index("--cutoff") + 1]) if "--cutoff" in sys.argv else 500
SUFFIX = f"_{cluster.TAG}_{CUTOFF}"
MOVED, SHIFT, DELTA = 2, 0.05, 0.01  # index of the moved O, its first shift and the finite-difference step (A)

DFT = {
    "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
    "MGRID": {"NGRIDS": 4, "CUTOFF": CUTOFF, "REL_CUTOFF": 60},
    "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
    "POISSON": {"PERIODIC": "XYZ", "PSOLVER": "PERIODIC"},
    "SCF": {"EPS_SCF": 1.0e-8, "MAX_SCF": 300, "ADDED_MOS": 20,
            "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
            "MIXING": {"METHOD": "BROYDEN_MIXING", "ALPHA": 0.2, "NBUFFER": 8},
            "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": 300}},
}


def structures():
    bulk = psteros.rutile_bulk("Ru", a=4.492, c=3.107)
    assert bulk[MOVED].specie.symbol == "O"
    result = {}
    for label, extra in (("base", 0.0), ("plus", DELTA), ("minus", -DELTA)):
        moved = bulk.copy()
        moved.translate_sites([MOVED], [SHIFT + extra, 0.0, 0.0], frac_coords=False, to_unit_cell=False)
        result[label] = moved
    return result


def config():
    return psteros.SurfaceWorkflowConfig(
        backend="cp2k", name=f"ruo2_c0_smoke{SUFFIX}",
        execution=cluster.execution("debug", 2, 1800),
        calculation=psteros.Cp2kCalculationConfig(
            code_label=CODE, parameters={"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": DFT}},
            basis_sets={"Ru": "DZVP-MOLOPT-PBE-GTH-q16", "O": "DZVP-MOLOPT-PBE-GTH-q6"},
            potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
            basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH",
            kpoints_spacing=0.30, max_iterations=1),
    )


def submit() -> None:
    graph = psteros.build_cp2k_surface_workgraph(structures(), config(), submit=True)
    record = {"workgraph_pk": graph.pk, "labels": list(structures()), "code": CODE, "moved_site": MOVED,
              "shift_angstrom": SHIFT, "delta_angstrom": DELTA, "cutoff_ry": CUTOFF}
    (HERE / f"c0_graph{SUFFIX}.json").write_text(json.dumps(record, indent=2) + "\n")
    print(record)


def check() -> None:
    from aiida import orm

    record = json.loads((HERE / f"c0_graph{SUFFIX}.json").read_text())
    node = orm.load_node(record["workgraph_pk"])
    labels = record["labels"]
    energies, _ = psteros.read_cp2k_results(node, labels)
    forces = psteros.read_cp2k_forces(node, labels)
    dos = psteros.read_cp2k_dos(node, ["base"])["base"]
    report = {"cutoff_ry": CUTOFF, "energies_ev": energies}
    for label in labels:
        text = node.outputs[f"{label}_retrieved"].base.repository.get_object_content("aiida.out")
        ranks = re.findall(r"Total number of message passing processes\s+(\d+)", text)
        version = re.findall(r"CP2K\| version string:\s+(.*)", text)
        entropy = re.findall(r"Electronic entropic energy:\s+(\S+)", text)
        report[label] = {"mpi_ranks": ranks, "cp2k": version[:1], "last_entropic_energy_hartree": entropy[-1:],
                         "max_force_ev_per_angstrom": float(np.abs(forces[label]).max())}
    finite = -(energies["plus"] - energies["minus"]) / (2 * record["delta_angstrom"])
    analytic = float(forces["base"][record["moved_site"], 0])
    report["finite_difference_force_x_ev_per_angstrom"] = finite
    report["printed_force_x_ev_per_angstrom"] = analytic
    report["difference_ev_per_angstrom"] = finite - analytic
    report["dos_partially_occupied_bins"] = dos.partially_occupied_bins()
    psteros.require_metallic("base", partially_occupied=report["dos_partially_occupied_bins"])
    report["energy_consistent_with_forces"] = abs(finite - analytic) < 0.01
    print(json.dumps(report, indent=2))
    (HERE / f"c0_check{SUFFIX}.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    from aiida import load_profile

    load_profile(PROFILE)
    {"submit": submit, "check": check}[sys.argv[1]]()
