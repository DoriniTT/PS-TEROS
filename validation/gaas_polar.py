#!/usr/bin/env python
"""End-to-end validation of the polar-surface workflow on GaAs (VASP, Bohr par40).

Phases (each builds the WorkGraphs; nothing runs without --submit):

  refs     relax + static of the GaAs bulk (ISIF=3) and of the Ga and As solids
           (ISIF=3). One WorkGraph per structure, so three jobs can run at once.
  checks   pseudo-hydrogen molecules, the Eq. 7 slab and the (110) pair, cut from
           the *relaxed* bulk. These must pass before the polar numbers are trusted.
  slabs    the (111) and (-1-1-1) terminations on their shared passivated bottoms.
  analyse  gather every energy, run the bottom check and the consistency checks,
           and write the phase diagram (png + csv) and a text report.

Run from the repository root with the branch on PYTHONPATH:

    PYTHONPATH=$PWD python validation/gaas_polar.py refs --submit
    PYTHONPATH=$PWD python validation/gaas_polar.py checks --submit
    PYTHONPATH=$PWD python validation/gaas_polar.py slabs --submit
    PYTHONPATH=$PWD python validation/gaas_polar.py analyse

The PKs of the submitted graphs are kept in validation/gaas/runs.json.
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
POTCARS = {"Ga": "Ga_d", "As": "As"}

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gaas")
RUNS = os.path.join(OUT, "runs.json")

# Settings of Zhang et al., Sci. Rep. 6, 20055 (2016): PBE, 400 eV, forces 0.005 eV/A.
ELECTRONIC = {
    "ENCUT": 400, "PREC": "Accurate", "EDIFF": 1.0e-6, "ISMEAR": 0, "SIGMA": 0.05,
    "LREAL": False, "LWAVE": False, "LCHARG": False, "NCORE": 8,
}
RELAX = {"IBRION": 2, "NSW": 200, "ISIF": 2, "EDIFFG": -0.005}
STATIC = {"IBRION": -1, "NSW": 0}
# psteros states kpoints_spacing in 1/A with the 2*pi of VASP's KSPACING and converts it for aiida-vasp.
KPOINTS_SPACING = 0.15

BILAYERS = 9


def starting_structures():
    """Experimental starting points; the cells are relaxed in phase 'refs'."""
    bulk = Structure.from_spacegroup("F-43m", Lattice.cubic(5.75), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])
    gallium = Structure.from_spacegroup("Cmce", Lattice.orthorhombic(4.52, 7.66, 4.53), ["Ga"], [[0, 0.1549, 0.081]])
    arsenic = Structure.from_spacegroup("R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]])
    return bulk, gallium, arsenic


def execution(walltime_h: float) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=COMPUTER, queue=QUEUE, max_concurrent_jobs=1, with_mpi=True,
        # par40 requires ncpus=40 explicitly (queue resources_min.ncpus = 40).
        resources={"num_machines": 1, "num_mpiprocs_per_machine": RANKS, "num_cores_per_machine": RANKS},
        max_wallclock_seconds=int(walltime_h * 3600),
    )


def recipes(name, mapping, relax_overrides, static_overrides, walltime_h):
    def vasp(ionic):
        return psteros.VaspCalculationConfig(
            code_label=CODE, incar={**ELECTRONIC, **ionic}, potential_family=POTENTIAL_FAMILY,
            potential_mapping=mapping, kpoints_spacing=KPOINTS_SPACING, max_iterations=3,
        )

    relax = psteros.SurfaceWorkflowConfig(
        backend="vasp", calculation=vasp(RELAX), execution=execution(walltime_h),
        name=name, role_overrides=relax_overrides,
    )
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", calculation=vasp(STATIC), execution=execution(walltime_h),
        name=f"{name}_static", role_overrides=static_overrides,
    )
    return relax, static


def load_runs():
    return json.load(open(RUNS)) if os.path.exists(RUNS) else {}


def save_run(name, graph, labels):
    runs = load_runs()
    runs[name] = {"pk": graph.pk, "labels": list(labels)}
    os.makedirs(OUT, exist_ok=True)
    json.dump(runs, open(RUNS, "w"), indent=2)


def launch(name, structures, study, walltime_h, submit, overrides_extra=None):
    """One relax -> static WorkGraph holding ``structures`` (label -> pymatgen Structure)."""
    relax_ov = study.vasp_overrides("relax", extra=overrides_extra)
    static_ov = study.vasp_overrides("static")
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
    """Relaxed bulk and reference solids from the finished 'refs' graphs.

    The slab phases only need the relaxed bulk (the references enter the analysis, not the
    structures), so ``require_all=False`` skips references that have not finished yet.
    """
    from aiida import orm

    runs = load_runs()
    out = {}
    for label in ("bulk", "ref_Ga", "ref_As"):
        try:
            graph = orm.load_node(runs[f"gaas_{label}"]["pk"])
            energies, relaxed = psteros.read_vasp_results(graph, [label])
            out[label] = (energies[label], relaxed[label].get_pymatgen_structure())
        except Exception:
            if require_all or label == "bulk":
                raise
    return out


# Four sizes = the four unknowns of Eq. 9 (as in the paper; no residual). n=6 (175 atoms) needs about 3-4 h per
# ionic step and was dropped.
CLUSTER_SIZES = (2, 3, 4, 5)


def build_study(bulk, gallium, arsenic, method="molecules"):
    return psteros.PolarSurfaceStudy(
        bulk, faces=[(1, 1, 1), (-1, -1, -1)], references={"Ga": gallium, "As": arsenic},
        bilayers=BILAYERS, pseudo_hydrogen_method=method, cluster_sizes=CLUSTER_SIZES,
        nonpolar_check=(1, 1, 0),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("refs", "checks", "eq7", "slabs", "clusters", "analyse"))
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--method", choices=("molecules", "clusters"), default="molecules",
                        help="pseudo-hydrogen method of the study (analyse, clusters)")
    parser.add_argument("--quasi-newton", action="store_true", help="clusters --only: relax with IBRION=1")
    parser.add_argument("--only", help="refs: one label (bulk, ref_Ga, ref_As); slabs: one face tag (111, m1m1m1)")
    args = parser.parse_args()

    from aiida import load_profile

    load_profile(PROFILE)

    if args.phase == "refs":
        bulk, gallium, arsenic = starting_structures()
        study = build_study(bulk, gallium, arsenic)
        # Cell relaxation for all three solids (the study fixes the cell of its bulk).
        cell = {"ISIF": 3}
        extra = {label: psteros.CalculationOverride(parameters={"INCAR": cell}) for label in ("bulk", "ref_Ga", "ref_As")}
        # alpha-Ga has a very flat energy surface: the CG line search (IBRION=2) fails with
        # "ZBRENT: fatal error in bracketing" near convergence; the quasi-Newton optimiser does not.
        extra["ref_Ga"] = psteros.CalculationOverride(parameters={"INCAR": {**cell, "IBRION": 1}})
        for label, structure in (("bulk", bulk), ("ref_Ga", gallium), ("ref_As", arsenic)):
            if args.only and label != args.only:
                continue
            launch(f"gaas_{label}", {label: structure}, study, walltime_h=6, submit=args.submit,
                   overrides_extra=extra)
        return

    relaxed = relaxed_inputs(require_all=args.phase == "analyse")
    bulk = relaxed["bulk"][1]
    _, start_ga, start_as = starting_structures()
    method = "clusters" if args.phase == "clusters" else args.method
    study = build_study(bulk, relaxed["ref_Ga"][1] if "ref_Ga" in relaxed else start_ga,
                        relaxed["ref_As"][1] if "ref_As" in relaxed else start_as, method)
    print(study)
    print(f"relaxed bulk a = {bulk.lattice.a:.4f} A")

    if args.phase == "checks":
        labels = [l for l, r in study.roles.items() if r in ("pseudo_molecule", "eq7_check", "nonpolar_symmetric", "nonpolar_passivated")]
        launch("gaas_checks", {l: study.structures[l] for l in labels}, study, walltime_h=24, submit=args.submit)
    elif args.phase == "clusters" and args.only:
        # One named cluster in its own graph, e.g. --only cluster_Ga_n5
        # --quasi-newton: IBRION=1, for a cluster whose CG line search fails ("ZBRENT: fatal error in bracketing")
        extra = {args.only: psteros.CalculationOverride(parameters={"INCAR": {"IBRION": 1}})} if args.quasi_newton else None
        name = f"gaas_{args.only}" + ("_qn" if args.quasi_newton else "")
        launch(name, {args.only: study.structures[args.only]}, study, walltime_h=48, submit=args.submit, overrides_extra=extra)
    elif args.phase == "clusters":
        # Largest clusters get a graph each; the small ones share one. Three graphs = three jobs.
        biggest = max(CLUSTER_SIZES)
        groups = {
            "gaas_clusters_As_big": [f"cluster_As_n{biggest}"],
            "gaas_clusters_Ga_big": [f"cluster_Ga_n{biggest}"],
            "gaas_clusters_small": [l for l in study.roles if study.roles[l] == "cluster" and not l.endswith(f"_n{biggest}")],
        }
        for name, labels in groups.items():
            launch(name, {l: study.structures[l] for l in labels}, study, walltime_h=48, submit=args.submit)
    elif args.phase == "eq7":
        label = "GaAs_111_both_passivated"
        launch("gaas_eq7", {label: study.structures[label]}, study, walltime_h=24, submit=args.submit)
    elif args.phase == "slabs":
        for prefix, tag in (("GaAs_111", "111"), ("GaAs_m1m1m1", "m1m1m1")):
            if args.only and tag != args.only:
                continue
            labels = [l for l, r in study.roles.items() if r == "slab" and l.startswith(prefix + "_")]
            launch(f"gaas_slabs_{tag}", {l: study.structures[l] for l in labels}, study, walltime_h=48,
                   submit=args.submit)
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
    text = [str(result)]
    for face, report in result.bottom_checks.items():
        text += [f"\n{face}", str(report)]
    suffix = "" if method == "molecules" else f"_{method}"
    report_path = os.path.join(OUT, f"report{suffix}.txt")
    open(report_path, "w").write("\n".join(text) + "\n")
    print("\n".join(text))
    print(result.diagram.plot(os.path.join(OUT, f"gaas_polar_phase_diagram{suffix}.png"), title=f"GaAs polar faces ({method})"))
    print(result.diagram.to_csv(os.path.join(OUT, f"gaas_polar_phase_diagram{suffix}.csv")))
    json.dump({k: float(v) for k, v in energies.items()}, open(os.path.join(OUT, f"energies{suffix}.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
