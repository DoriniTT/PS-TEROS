#!/usr/bin/env python
"""End-to-end validation of the wurtzite polar-surface workflow on ZnO (VASP, Bohr par40).

Benchmark: ideal 1x1 ZnO(0001) (Zn-terminated) and ZnO(000-1) (O-terminated) at the O-rich limit,
147.7 and 63.1 meV/A^2 (Zhang et al., arXiv:1510.08961, Table 3, GGA, cluster method).

Phases (nothing runs without --submit), same as gaas_polar.py:

  refs      relax + static of ZnO bulk and Zn (ISIF=3) and of the O2 molecule (triplet, Gamma).
  checks    pseudo-hydrogen molecules, the Eq. 7 slab and the (10-10) non-polar pair, cut from
            the relaxed bulk.
  slabs     ideal (0001) and (000-1) slabs on their shared passivated bottoms.
  clusters  tetrahedral clusters (zinc-blende analogue) for the cluster method.
  analyse   energies -> bottom checks, consistency checks, phase diagram.

    PYTHONPATH=$PWD python validation/zno_polar.py refs --submit
"""

from __future__ import annotations

import argparse
import json
import math
import os

from pymatgen.core import Lattice, Structure

import psteros

PROFILE = "psteros_test"
CODE = "vasp-std-cpu-h5old@bohr"
COMPUTER = "bohr"
QUEUE = "par40"
RANKS = 40
POTENTIAL_FAMILY = "PBE"
POTCARS = {"Zn": "Zn", "O": "O"}

HERE = os.path.dirname(os.path.abspath(__file__))
# Convergence study: ZNO_TAG=conv ZNO_ENCUT=600 ZNO_KSPACING=0.10 writes to validation/zno_conv with its own runs.json.
TAG = os.environ.get("ZNO_TAG", "")
OUT = os.path.join(HERE, "zno" + (f"_{TAG}" if TAG else ""))
RUNS = os.path.join(OUT, "runs.json")

# O has ENMAX 400 eV, so 500 eV is >= 1.25 ENMAX for the cell relaxations.
ELECTRONIC = {
    "ENCUT": int(os.environ.get("ZNO_ENCUT", 500)), "PREC": "Accurate", "EDIFF": 1.0e-6, "ISMEAR": 0, "SIGMA": 0.05,
    "LREAL": False, "LWAVE": False, "LCHARG": False, "NCORE": 8,
}
RELAX = {"IBRION": 2, "NSW": 200, "ISIF": 2, "EDIFFG": -0.005}
STATIC = {"IBRION": -1, "NSW": 0}
# psteros states kpoints_spacing in 1/A with the 2*pi of VASP's KSPACING (10-15 k-points for a 3.25 A cell).
KPOINTS_SPACING = float(os.environ.get("ZNO_KSPACING", 0.15))
BILAYERS = 9
CLUSTER_SIZES = (2, 3, 4, 5)  # the four unknowns of Eq. 9; n=6 (175 atoms) is too slow here
FACES = [(0, 0, 1), (0, 0, -1)]


def starting_structures():
    """Experimental starting points; the cells are relaxed in phase 'refs'."""
    bulk = Structure.from_spacegroup("P6_3mc", Lattice.hexagonal(3.25, 5.21), ["Zn", "O"],
                                     [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])
    zinc = Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.66, 4.95), ["Zn"], [[1 / 3, 2 / 3, 0.25]])
    return bulk, zinc, psteros.triplet_o2_cell(cell_length=12.0)


def execution(walltime_h: float) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=COMPUTER, queue=QUEUE, max_concurrent_jobs=1, with_mpi=True,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": RANKS, "num_cores_per_machine": RANKS},
        max_wallclock_seconds=int(walltime_h * 3600),
    )


def recipes(name, mapping, relax_overrides, static_overrides, walltime_h):
    def vasp(ionic):
        return psteros.VaspCalculationConfig(
            code_label=CODE, incar={**ELECTRONIC, **ionic}, potential_family=POTENTIAL_FAMILY,
            potential_mapping=mapping, kpoints_spacing=KPOINTS_SPACING, max_iterations=3,
        )

    relax = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp(RELAX), execution=execution(walltime_h),
                                          name=name, role_overrides=relax_overrides)
    static = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp(STATIC), execution=execution(walltime_h),
                                           name=f"{name}_static", role_overrides=static_overrides)
    return relax, static


def load_runs():
    return json.load(open(RUNS)) if os.path.exists(RUNS) else {}


def save_run(name, graph, labels):
    runs = load_runs()
    runs[name] = {"pk": graph.pk, "labels": list(labels)}
    os.makedirs(OUT, exist_ok=True)
    json.dump(runs, open(RUNS, "w"), indent=2)


# Ideal polar slabs: charge reorganises slowly across the slab on a cold start, so the SCF energy drifts
# steadily and 60 (even 170) iterations are not enough. Only the iteration limit changes, not the physics.
# Tried and rejected: AMIX=0.2/BMIX=1e-4 made the SCF diverge (E = -1073 eV by iteration 7).
SLAB_SCF = {"NELM": 500}


def slab_scf_overrides(labels):
    return {label: psteros.CalculationOverride(parameters={"INCAR": SLAB_SCF}) for label in labels}


def launch(name, structures, study, walltime_h, submit, overrides_extra=None, static_extra=None):
    relax_ov = study.vasp_overrides("relax", extra=overrides_extra)
    static_ov = study.vasp_overrides("static", extra=static_extra)
    relax, static = recipes(
        name, study.potential_mapping(POTCARS),
        {k: v for k, v in relax_ov.items() if k in structures},
        {k: v for k, v in static_ov.items() if k in structures}, walltime_h,
    )
    graph = psteros.build_relax_static_workgraph(structures, relax, static, submit=submit)
    atoms = ", ".join(f"{k}:{len(v)}" for k, v in structures.items())
    print(f"graph {name!r}: {len(structures)} structures ({atoms})" + (f" -> PK {graph.pk}" if submit else " [not submitted]"))
    if submit:
        save_run(name, graph, structures)
    return graph


def relaxed_inputs(require_all=True):
    from aiida import orm

    runs = load_runs()
    out = {}
    for label in ("bulk", "ref_Zn", "ref_O"):
        try:
            run = next(r for r in runs.values() if label in r["labels"] and not label.startswith("cluster"))
            graph = orm.load_node(run["pk"])
            energies, relaxed = psteros.read_vasp_results(graph, [label])
            out[label] = (energies[label], relaxed[label].get_pymatgen_structure())
        except Exception:
            if require_all or label == "bulk":
                raise
    return out


def build_study(bulk, zinc, oxygen, method="molecules"):
    return psteros.PolarSurfaceStudy(
        bulk, faces=FACES, references={"Zn": zinc, "O": oxygen}, bilayers=BILAYERS,
        pseudo_hydrogen_method=method, cluster_sizes=CLUSTER_SIZES, nonpolar_check=(1, 0, 0),
        electron_counting=False, oxidation_states={"Zn": 2, "O": -2},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("refs", "checks", "slabs", "clusters", "analyse"))
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--method", choices=("molecules", "clusters"), default="molecules")
    parser.add_argument("--scf", choices=("both", "static", "none"), default="both",
                        help="slabs: SCF mixing help (NELM/AMIX/BMIX) for both stages, the static one only "
                             "(keeps a finished relaxation in the AiiDA cache), or none")
    parser.add_argument("--algo-all-static", action="store_true", help="slabs: ALGO=All for the static step only")
    parser.add_argument("--only", help="refs: one label (bulk, ref_Zn, ref_O); slabs: one face tag (001, 00m1)")
    parser.add_argument("--unrelaxed", action="store_true", help="dry run: build everything from the starting structures")
    args = parser.parse_args()

    from aiida import load_profile

    load_profile(PROFILE)
    start_bulk, start_zn, start_o2 = starting_structures()

    if args.phase == "refs":
        study = build_study(start_bulk, start_zn, start_o2)
        # ZNO_QN=1: quasi-Newton (IBRION=1) for the cell relaxations; the conjugate-gradient line search
        # can fail with "ZBRENT: fatal error in bracketing" on a flat energy surface (ZnO at 600 eV did).
        cell = {"ISIF": 3, **({"IBRION": 1} if os.environ.get("ZNO_QN") else {})}
        extra = {label: psteros.CalculationOverride(parameters={"INCAR": cell}) for label in ("bulk", "ref_Zn")}
        # One graph = one job at a time (the limit is 3 concurrent jobs for the whole campaign).
        structures = {label: structure for label, structure in (("bulk", start_bulk), ("ref_Zn", start_zn), ("ref_O", start_o2))
                      if not args.only or label == args.only}
        launch("zno_refs", structures, study, walltime_h=12, submit=args.submit, overrides_extra=extra)
        return

    if args.unrelaxed:
        relaxed = {"bulk": (0.0, start_bulk), "ref_Zn": (0.0, start_zn), "ref_O": (0.0, start_o2)}
    else:
        relaxed = relaxed_inputs(require_all=args.phase == "analyse")
    bulk = relaxed["bulk"][1]
    method = "clusters" if args.phase == "clusters" else args.method
    study = build_study(bulk, relaxed["ref_Zn"][1] if "ref_Zn" in relaxed else start_zn,
                        relaxed["ref_O"][1] if "ref_O" in relaxed else start_o2, method)
    print(study)
    print(f"bulk a = {bulk.lattice.a:.4f} A, c = {bulk.lattice.c:.4f} A")

    if args.phase == "checks":
        labels = [l for l, r in study.roles.items()
                  if r in ("pseudo_molecule", "eq7_check", "nonpolar_symmetric", "nonpolar_passivated")]
        launch("zno_checks", {l: study.structures[l] for l in labels}, study, walltime_h=24, submit=args.submit)
    elif args.phase == "slabs":
        for prefix, tag in (("ZnO_001", "001"), ("ZnO_00m1", "00m1")):
            if args.only and tag != args.only:
                continue
            labels = [l for l, r in study.roles.items() if r == "slab" and l.startswith(prefix + "_")]
            scf = slab_scf_overrides(labels) if args.scf != "none" else None
            static_scf = scf
            if args.algo_all_static:
                # ALGO=All for the static step only: the default solver can diverge on a cold start of a polar slab
                # (seen at 600 eV with a denser mesh). The relaxation INCAR stays as before, so its result is cached.
                static_scf = {l: psteros.CalculationOverride(parameters={"INCAR": {**SLAB_SCF, "ALGO": "All"}}) for l in labels}
            launch(f"zno_slabs_{tag}", {l: study.structures[l] for l in labels}, study, walltime_h=48, submit=args.submit,
                   overrides_extra=scf if args.scf == "both" else None, static_extra=static_scf)
    elif args.phase == "clusters" and args.only == "small":
        # every cluster below the largest size, in one graph (one job at a time)
        biggest = max(CLUSTER_SIZES)
        labels = [l for l in study.roles if study.roles[l] == "cluster" and not l.endswith(f"_n{biggest}")]
        launch("zno_clusters_small", {l: study.structures[l] for l in labels}, study, walltime_h=48, submit=args.submit)
    elif args.phase == "clusters" and args.only:
        launch(f"zno_{args.only}", {args.only: study.structures[args.only]}, study, walltime_h=48, submit=args.submit)
    elif args.phase == "clusters":
        biggest = max(CLUSTER_SIZES)
        groups = {
            "zno_clusters_Zn_big": [f"cluster_Zn_n{biggest}"],
            "zno_clusters_O_big": [f"cluster_O_n{biggest}"],
            "zno_clusters_small": [l for l in study.roles if study.roles[l] == "cluster" and not l.endswith(f"_n{biggest}")],
        }
        for name, labels in groups.items():
            launch(name, {l: study.structures[l] for l in labels}, study, walltime_h=48, submit=args.submit)
    else:
        analyse(study, relaxed, method)


def read_from_tasks(graph, label):
    """Energy and relaxed structure of one label straight from the finished work chains of ``graph``.

    A graph publishes its outputs only when it finishes normally; one that was killed or failed after
    some of its tasks finished still holds their results.
    """
    from aiida.common.links import LinkType

    calls = {link.link_label: link.node for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all()}
    static, relax = calls[f"{label}_static_vasp"], calls[f"{label}_relax_vasp"]
    if static.exit_status != 0 or relax.exit_status != 0:
        raise ValueError(f"{label}: the relaxation or the static calculation did not finish")
    return static.outputs.misc.get_dict()["total_energies"]["energy_extrapolated"], relax.outputs.structure


def analyse(study, relaxed, method="molecules"):
    from aiida import orm

    energies = {label: e for label, (e, _) in relaxed.items()}
    structures = {label: s for label, (_, s) in relaxed.items()}
    for name, run in load_runs().items():
        if name in ("zno_refs", "gaas_bulk", "gaas_ref_Ga", "gaas_ref_As"):
            continue
        for label in run["labels"]:
            if label.startswith("cluster") and method != "clusters":
                continue  # cluster calculations only enter the cluster method
            graph = orm.load_node(run["pk"])
            try:  # a label without a finished static calculation (still running, or failed and rerun elsewhere)
                e, s = psteros.read_vasp_results(graph, [label])
            except ValueError:
                try:  # a killed graph publishes no outputs, but its finished tasks hold the results
                    energy, structure = read_from_tasks(graph, label)
                except (KeyError, ValueError, AttributeError):
                    continue
                e, s = {label: energy}, {label: structure}
            energies.update(e)
            structures.update(s)
    result = study.analyse(energies, structures)
    os.makedirs(OUT, exist_ok=True)
    suffix = "" if method == "molecules" else f"_{method}"
    text = [str(result)]
    for face, report in result.bottom_checks.items():
        text += [f"\n{face}", str(report)]
    open(os.path.join(OUT, f"report{suffix}.txt"), "w").write("\n".join(text) + "\n")
    print("\n".join(text))
    print(result.diagram.plot(os.path.join(OUT, f"zno_polar_phase_diagram{suffix}.png"), title=f"ZnO polar faces ({method})"))
    print(result.diagram.to_csv(os.path.join(OUT, f"zno_polar_phase_diagram{suffix}.csv")))
    json.dump({k: float(v) for k, v in energies.items()}, open(os.path.join(OUT, f"energies{suffix}.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
