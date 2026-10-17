"""Steps C2 and C3 of docs/plans/MACE_SURFACE_PLAN.md: CP2K references, training labels and ground truth.

The CP2K settings are those chosen in C1 and written to ``settings.json`` (basis, cutoff, k-spacing, electronic
temperature, vacuum), so every CP2K energy of the validation comes from one recipe.

Phase ``refs``    bulk RuO2 and hcp Ru: CELL_OPT -> static; triplet O2: GEO_OPT -> static.
Phase ``labels``  ENERGY_FORCE of the six RuO2(110) terminations as built on the CP2K lattice (R0, also the training
                  labels), of bulk RuO2 in the (110)-oriented cell of the slabs (``bulk110``), and of 5 rattled copies
                  of every termination (σ = 0.05 A, fixed seed; the M3 training set).
Phase ``truth``   GEO_OPT -> ENERGY_FORCE of the six terminations, central trilayer fixed (R1, the ground truth).

Usage (PYTHONPATH=<worktree>):
    python c2_cp2k.py refs   [--dry]
    python c2_cp2k.py labels [--dry]     # after refs; reads c2_refs.json
    python c2_cp2k.py truth  [--dry]     # after refs
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import psteros

import cluster

HERE = Path(__file__).resolve().parent
PROFILE = "psteros_mace_ruo2"
CODE = cluster.CP2K_CODE
RATTLES, SIGMA, SEED = 5, 0.05, 20261017
GEO = {"GEO_OPT": {"OPTIMIZER": "BFGS", "MAX_FORCE": 1.94e-4, "MAX_ITER": 300}}  # 0.01 eV/A
CELL = psteros.CalculationOverride(parameters={
    "GLOBAL": {"RUN_TYPE": "CELL_OPT"},
    "FORCE_EVAL": {"STRESS_TENSOR": "ANALYTICAL"},
    "MOTION": {"CELL_OPT": {"OPTIMIZER": "BFGS", "MAX_FORCE": 1.94e-4, "KEEP_SYMMETRY": True, "KEEP_ANGLES": True,
                            "PRESSURE_TOLERANCE": 50.0, "MAX_ITER": 200}},
})
O2_PARAMETERS = {"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3, "SCF": {"SMEAR": {"_": "OFF"}}}}}


def settings() -> dict:
    return json.loads((HERE / "settings.json").read_text())


def execution(concurrent: int) -> psteros.ExecutionPolicy:
    return cluster.execution("node", concurrent, 8 * 3600)


def recipe(run_type: str) -> psteros.Cp2kCalculationConfig:
    s = settings()
    parameters = {
        "GLOBAL": {"RUN_TYPE": run_type},
        "FORCE_EVAL": {"METHOD": "Quickstep", "DFT": {
            "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
            "MGRID": {"NGRIDS": 4, "CUTOFF": s["cutoff"], "REL_CUTOFF": 60},
            "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
            "POISSON": {"PERIODIC": "XYZ", "PSOLVER": "PERIODIC"},
            "SCF": {"EPS_SCF": 1.0e-7, "MAX_SCF": 300, "ADDED_MOS": 40,
                    "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
                    "MIXING": {"METHOD": "BROYDEN_MIXING", "ALPHA": 0.2, "NBUFFER": 8},
                    "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": s["temperature"]}},
        }},
    }
    if run_type == "GEO_OPT":
        parameters["MOTION"] = GEO
    basis = s["basis"]
    return psteros.Cp2kCalculationConfig(
        code_label=CODE, parameters=parameters,
        basis_sets={"Ru": f"{basis}-MOLOPT-PBE-GTH-q16", "O": f"{basis}-MOLOPT-PBE-GTH-q6"},
        potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
        basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH",
        kpoints_spacing=s["kspacing"], max_iterations=3)


def config(run_type, name, concurrent, overrides=None):
    return psteros.SurfaceWorkflowConfig(backend="cp2k", name=name, calculation=recipe(run_type),
                                         execution=execution(concurrent), role_overrides=overrides or {})


def refs(dry: bool) -> None:
    from pymatgen.core import Lattice, Structure

    structures = {
        "bulk": psteros.rutile_bulk("Ru", a=4.492, c=3.107),
        "ru": Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.706, 4.282), ["Ru"], [[1 / 3, 2 / 3, 0.25]]),
        "o2": psteros.triplet_o2_cell(cell_length=15.0),
    }
    o2 = psteros.CalculationOverride(parameters=O2_PARAMETERS)
    o2_relax = psteros.CalculationOverride(parameters={**O2_PARAMETERS, "MOTION": GEO})
    graph = psteros.build_cp2k_relax_static_workgraph(
        structures, config("GEO_OPT", "ruo2_c2_refs", 2, {"bulk": CELL, "ru": CELL, "o2": o2_relax}),
        config("ENERGY_FORCE", "ruo2_c2_refs_static", 2, {"o2": o2}), submit=not dry)
    print(len(graph.tasks), "tasks")
    if not dry:
        (HERE / "c2_refs.json").write_text(json.dumps({"workgraph_pk": graph.pk, "settings": settings()}, indent=2))
        print("WorkGraph", graph.pk)


def relaxed_bulk():
    from aiida import orm

    node = orm.load_node(json.loads((HERE / "c2_refs.json").read_text())["workgraph_pk"])
    from psteros.backends.cp2k import _stage_output

    bulk = _stage_output(node, "bulk", "structure").get_pymatgen()
    a, c = (bulk.lattice.a + bulk.lattice.b) / 2, bulk.lattice.c
    # psteros.rutile_bulk puts M at (0, 0, 0) first and O at (u, u, 0) third; CP2K keeps the order
    assert bulk[0].specie.symbol == "Ru" and bulk[2].specie.symbol == "O"
    u = float(np.mod(bulk[2].frac_coords[0] - bulk[0].frac_coords[0], 1.0))
    assert 0.28 < u < 0.33, u
    return a, c, u


def terminations():
    from pymatgen.core import Structure

    a, c, u = relaxed_bulk()
    s = settings()
    result, oriented = {}, None
    for termination in psteros.RUTILE_110_TERMINATIONS:
        slab, _ = psteros.rutile_110_slab(metal="Ru", termination=termination, a=a, c=c, u=u,
                                          vacuum_angstrom=s["vacuum"])
        if termination == "o_bridge":
            oriented = slab.oriented_unit_cell
        result[termination] = Structure(slab.lattice, slab.species, slab.frac_coords)
    return result, Structure(oriented.lattice, oriented.species, oriented.frac_coords), (a, c, u)


def labels(dry: bool) -> None:
    slabs, bulk110, lattice = terminations()
    structures = {**slabs, "bulk110": bulk110}
    rng = np.random.default_rng(SEED)
    for termination, slab in slabs.items():
        for index in range(RATTLES):
            rattled = slab.copy()
            for site_index in range(len(rattled)):
                rattled.translate_sites([site_index], rng.normal(0.0, SIGMA, 3), frac_coords=False, to_unit_cell=True)
            structures[f"{termination}_r{index}"] = rattled
    # labels and truth may run together: one CP2K job each, two in all
    graph = psteros.build_cp2k_surface_workgraph(structures, config("ENERGY_FORCE", "ruo2_c2_labels", 1),
                                                 submit=not dry)
    print(len(graph.tasks), "tasks; lattice a, c, u =", lattice)
    if not dry:
        (HERE / "c2_labels.json").write_text(json.dumps(
            {"workgraph_pk": graph.pk, "lattice_a_c_u": lattice, "labels": list(structures),
             "rattle": {"count": RATTLES, "sigma_angstrom": SIGMA, "seed": SEED}}, indent=2))
        print("WorkGraph", graph.pk)


def truth(dry: bool) -> None:
    slabs, _, lattice = terminations()
    fixed = {label: psteros.CalculationOverride(fixed_sites=psteros.central_sites(slab))
             for label, slab in slabs.items()}
    graph = psteros.build_cp2k_relax_static_workgraph(
        slabs, config("GEO_OPT", "ruo2_c3_truth", 1, fixed), config("ENERGY_FORCE", "ruo2_c3_truth_static", 1),
        submit=not dry)
    print(len(graph.tasks), "tasks; fixed sites:", {label: len(o.fixed_sites) for label, o in fixed.items()})
    if not dry:
        (HERE / "c3_truth.json").write_text(json.dumps(
            {"workgraph_pk": graph.pk, "lattice_a_c_u": lattice,
             "fixed_sites": {label: list(o.fixed_sites) for label, o in fixed.items()}}, indent=2))
        print("WorkGraph", graph.pk)


if __name__ == "__main__":
    from aiida import load_profile

    load_profile(PROFILE)
    {"refs": refs, "labels": labels, "truth": truth}[sys.argv[1]]("--dry" in sys.argv)
