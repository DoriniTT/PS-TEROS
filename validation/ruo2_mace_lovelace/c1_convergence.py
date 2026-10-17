"""Step C1 of docs/plans/MACE_SURFACE_PLAN.md: CP2K convergence for RuO2(110), before any training label.

Structures (experimental lattices; only differences between settings matter here): bulk RuO2, hcp Ru, triplet O2
and the unrelaxed (1x1) RuO2(110) slabs o_bridge, o_cus and metal (3 trilayers, symmetric).

Settings: the baseline (DZVP-MOLOPT-PBE-GTH, 600 Ry, REL_CUTOFF 60, k-spacing 0.20 A^-1, Fermi-Dirac 300 K,
22 A vacuum) and one change at a time: 400 and 800 Ry, TZV2P, k-spacing 0.15, 150 K (at k-spacing 0.15),
30 A vacuum (slabs only). Quantities: ΔH_f(RuO2) = E(RuO2) - E(Ru) - E(O2) per formula unit, and, at Δμ_O = 0,
γ(o_bridge) - γ(o_cus) and γ(o_bridge) - γ(metal) of the unrelaxed slabs (meV/A^2).

Converged (fixed in the plan): the next level changes ΔH_f by <= 0.05 eV/f.u. and both γ differences by
<= 2 meV/A^2.

Reported with them (not a criterion of the plan, added after C0 showed a 9 % change of one force between 500 and
800 Ry): the largest change of a force component between the two levels, over the slabs and the bulk.

Also two CELL_OPT runs of bulk RuO2 (analytical stress, k-points) at 600 and 800 Ry: they check that the stress
works with k-points and give the lattice of the labels (a and c should move by less than 0.01 A between them).

Usage (PYTHONPATH=<worktree>):
    python c1_convergence.py submit [--dry]
    python c1_convergence.py check
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import psteros
from psteros.backends import add_cp2k_task

import cluster

HERE = Path(__file__).resolve().parent
PROFILE = "psteros_mace_ruo2"
CODE = cluster.CP2K_CODE
A, C = 4.492, 3.107
SLABS = ("o_bridge", "o_cus", "metal")
BASELINE = {"cutoff": 600, "basis": "DZVP", "kspacing": 0.20, "temperature": 300, "vacuum": 22.0}
LEVELS = {
    "base": {},
    "c400": {"cutoff": 400},
    "c800": {"cutoff": 800},
    "tzv2p": {"basis": "TZV2P"},
    "k015": {"kspacing": 0.15},
    "t150": {"kspacing": 0.15, "temperature": 150},
    "vac30": {"vacuum": 30.0},
}
# what each level recomputes; the rest is taken from the level it is compared with
RECOMPUTE = {
    "base": ("bulk", "ru", "o2") + SLABS, "c400": ("bulk", "ru", "o2") + SLABS, "c800": ("bulk", "ru", "o2") + SLABS,
    "tzv2p": ("bulk", "ru", "o2") + SLABS, "k015": ("bulk", "ru") + SLABS, "t150": ("bulk", "ru") + SLABS,
    "vac30": SLABS,
}
COMPARE = {"c400": "base", "c800": "base", "tzv2p": "base", "k015": "base", "t150": "k015", "vac30": "base"}
EXECUTION = cluster.execution("node", 2, 4 * 3600)
CELL = psteros.CalculationOverride(parameters={
    "GLOBAL": {"RUN_TYPE": "CELL_OPT"},
    "FORCE_EVAL": {"STRESS_TENSOR": "ANALYTICAL"},
    "MOTION": {"CELL_OPT": {"OPTIMIZER": "BFGS", "MAX_FORCE": 1.94e-4, "KEEP_SYMMETRY": True, "KEEP_ANGLES": True,
                            "PRESSURE_TOLERANCE": 50.0, "MAX_ITER": 200}},
})
CELL_LEVELS = ("base", "c800")
O2 = psteros.CalculationOverride(
    kpoints_distance=100.0,
    parameters={"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3, "SCF": {"SMEAR": {"_": "OFF"}}}}})


def settings(level: str) -> dict:
    return {**BASELINE, **LEVELS[level]}


def recipe(level: str) -> psteros.Cp2kCalculationConfig:
    s = settings(level)
    dft = {
        "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
        "MGRID": {"NGRIDS": 4, "CUTOFF": s["cutoff"], "REL_CUTOFF": 60},
        "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
        "POISSON": {"PERIODIC": "XYZ", "PSOLVER": "PERIODIC"},
        "SCF": {"EPS_SCF": 1.0e-7, "MAX_SCF": 300, "ADDED_MOS": 40,
                "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
                "MIXING": {"METHOD": "BROYDEN_MIXING", "ALPHA": 0.2, "NBUFFER": 8},
                "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": s["temperature"]}},
    }
    basis = s["basis"]
    return psteros.Cp2kCalculationConfig(
        code_label=CODE, parameters={"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": dft}},
        basis_sets={"Ru": f"{basis}-MOLOPT-PBE-GTH-q16", "O": f"{basis}-MOLOPT-PBE-GTH-q6"},
        potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
        basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH",
        kpoints_spacing=s["kspacing"], max_iterations=2)


def structures(level: str) -> dict:
    from pymatgen.core import Lattice, Structure

    s = settings(level)
    result = {
        "bulk": psteros.rutile_bulk("Ru", a=A, c=C),
        "ru": Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.706, 4.282), ["Ru"], [[1 / 3, 2 / 3, 0.25]]),
        "o2": psteros.triplet_o2_cell(cell_length=15.0),
    }
    for termination in SLABS:
        slab, _ = psteros.rutile_110_slab(metal="Ru", termination=termination, a=A, c=C,
                                          vacuum_angstrom=s["vacuum"])
        result[termination] = Structure(slab.lattice, slab.species, slab.frac_coords)
    return {label: structure for label, structure in result.items() if label in RECOMPUTE[level]}


def build():
    from aiida_workgraph import WorkGraph

    graph = WorkGraph(name="ruo2_c1_convergence")
    graph.max_number_jobs = EXECUTION.max_concurrent_jobs
    for level in LEVELS:
        config = recipe(level)
        for name, structure in structures(level).items():
            label = f"{level}_{name}"
            task = add_cp2k_task(graph, label=label, structure=structure, config=config, execution=EXECUTION,
                                 override=O2 if name == "o2" else None)
            graph.outputs.__setattr__(f"{label}_parameters", task.outputs.output_parameters)
            graph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
    for level in CELL_LEVELS:
        label = f"{level}_bulk_cell"
        task = add_cp2k_task(graph, label=label, structure=psteros.rutile_bulk("Ru", a=A, c=C),
                             config=recipe(level), execution=EXECUTION, override=CELL)
        graph.outputs.__setattr__(f"{label}_parameters", task.outputs.output_parameters)
        graph.outputs.__setattr__(f"{label}_structure", task.outputs.output_structure)
    return graph


def submit(dry: bool) -> None:
    graph = build()
    print(len(graph.tasks), "CP2K calculations")
    if dry:
        return
    graph.submit()
    (HERE / "c1_graph.json").write_text(json.dumps({"workgraph_pk": graph.pk, "levels": LEVELS,
                                                    "baseline": BASELINE}, indent=2) + "\n")
    print("WorkGraph", graph.pk)


def energy(node, level: str, name: str) -> float:
    while name not in RECOMPUTE[level]:
        level = COMPARE[level]
    return psteros.read_cp2k_results(node, [f"{level}_{name}"])[0][f"{level}_{name}"]


def quantities(node, level: str) -> dict:
    bulk = energy(node, level, "bulk") / 2.0  # two formula units in the conventional cell
    ru = energy(node, level, "ru") / 2.0
    o2 = energy(node, level, "o2")
    result = {"delta_h_f_ev": bulk - ru - o2}
    gamma = {}
    for name in SLABS:
        slab = structures("vac30" if level == "vac30" else "base")[name]
        lattice = slab.lattice.matrix
        area = float(np.linalg.norm(np.cross(lattice[0], lattice[1])))
        n_ru = sum(site.specie.symbol == "Ru" for site in slab)
        n_o = len(slab) - n_ru
        gamma[name] = 1000 * (energy(node, level, name) - n_ru * bulk - (n_o - 2 * n_ru) * o2 / 2) / (2 * area)
    result["gamma_mev_per_angstrom2"] = gamma
    result["obridge_minus_ocus"] = gamma["o_bridge"] - gamma["o_cus"]
    result["obridge_minus_metal"] = gamma["o_bridge"] - gamma["metal"]
    return result


def check() -> None:
    from aiida import orm

    node = orm.load_node(json.loads((HERE / "c1_graph.json").read_text())["workgraph_pk"])
    table = {level: quantities(node, level) for level in LEVELS}
    for level, reference in COMPARE.items():
        changes = {key: table[level][key] - table[reference][key]
                   for key in ("delta_h_f_ev", "obridge_minus_ocus", "obridge_minus_metal")}
        same_geometry = [name for name in ("bulk",) + SLABS if name in RECOMPUTE[level] and level != "vac30"]
        if same_geometry:
            current = psteros.read_cp2k_forces(node, [f"{level}_{name}" for name in same_geometry])
            before = psteros.read_cp2k_forces(node, [f"{reference}_{name}" for name in same_geometry])
            changes["max_force_change_ev_per_angstrom"] = max(
                float(np.abs(current[f"{level}_{name}"] - before[f"{reference}_{name}"]).max())
                for name in same_geometry)
        table[level]["change_against_" + reference] = changes
        table[level]["within_criteria"] = (abs(changes["delta_h_f_ev"]) <= 0.05
                                           and abs(changes["obridge_minus_ocus"]) <= 2.0
                                           and abs(changes["obridge_minus_metal"]) <= 2.0)
    for level in CELL_LEVELS:
        lattice = node.outputs[f"{level}_bulk_cell_structure"].get_pymatgen().lattice
        table[level]["relaxed_lattice_angstrom"] = {"a": lattice.a, "b": lattice.b, "c": lattice.c}
    print(json.dumps(table, indent=2))
    (HERE / "c1_check.json").write_text(json.dumps(table, indent=2) + "\n")


if __name__ == "__main__":
    from aiida import load_profile

    load_profile(PROFILE)
    if sys.argv[1] == "submit":
        submit(dry="--dry" in sys.argv)
    else:
        check()
