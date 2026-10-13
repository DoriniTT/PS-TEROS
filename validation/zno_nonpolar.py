#!/usr/bin/env python
"""Why does the non-polar consistency check of ZnO(10-10) not close with clusters?

The check compares gamma of a non-polar face from a symmetric slab (two clean surfaces) and from a slab
with one surface and a pseudo-hydrogen passivated bottom. This driver runs the pair at other thicknesses
and for the second non-polar face, (11-20), with the settings of the main ZnO run (500 eV, 0.15 1/A), to
tell a finite-thickness effect from a limit of the pseudo-hydrogen potentials.

    PYTHONPATH=$PWD python validation/zno_nonpolar.py submit 100 14 --submit
    PYTHONPATH=$PWD python validation/zno_nonpolar.py analyse
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zno_polar as z  # noqa: E402

import psteros  # noqa: E402

OUT = os.path.join(HERE, "zno_nonpolar")
RUNS = os.path.join(OUT, "runs.json")
OX = {"Zn": 2, "O": -2}
FACES = {"100": (1, 0, 0), "110": (1, 1, 0)}
ASYMMETRIC = {"ISIF": 2, "LDIPOL": True, "IDIPOL": 3, "DIPOL": [0.5, 0.5, 0.5]}


def build_pair(bulk, miller, layers):
    """Passivated one-sided slab and the stoichiometric symmetric slab of the same thickness."""
    from psteros.core.terminations import find_charge_neutral_terminations
    from psteros.polar import find_polar_terminations

    kwargs = {} if layers is None else {"layers": layers}
    one = find_polar_terminations(bulk, miller, vacuum=15.0, oxidation_states=OX, electron_counting=False, **kwargs)
    if one.polar:
        raise ValueError(f"{miller} is polar")
    reference = one[0]
    structure = reference.structure
    heavy = [site.coords for site, p in zip(structure, structure.site_properties["pseudo_hydrogen"]) if p is None]
    normal = structure.lattice.matrix[2] / structure.lattice.c
    heights = [float(x @ normal) for x in heavy]
    symmetric = [t for t in find_charge_neutral_terminations(bulk, miller, max(heights) - min(heights) - 0.1, 15.0,
                                                             oxidation_states=OX) if t.is_stoichiometric]
    if not symmetric:
        raise ValueError(f"no stoichiometric symmetric slab for {miller}")
    return symmetric[0], reference


def tag(face, layers):
    return f"ZnO_{face}_L{layers}"


def load_runs():
    return json.load(open(RUNS)) if os.path.exists(RUNS) else {}


def submit(face, layers, do_submit, only=None):
    from aiida import orm  # noqa: F401

    r = z.relaxed_inputs(True)
    study = z.build_study(r["bulk"][1], r["ref_Zn"][1], r["ref_O"][1])
    symmetric, passivated = build_pair(r["bulk"][1], FACES[face], layers)
    name = tag(face, layers)
    structures = {f"{name}_symmetric": symmetric.structure, f"{name}_passivated": passivated.structure}
    if only:  # one slab of the pair in its own graph, so that both run at the same time
        structures = {k: v for k, v in structures.items() if k.endswith("_" + only)}
    scf = {"NELM": 500}
    relax_extra = {
        f"{name}_symmetric": psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 2, **scf}}),
        f"{name}_passivated": psteros.CalculationOverride(parameters={"INCAR": {**ASYMMETRIC, **scf}}),
    }
    static_extra = {
        f"{name}_symmetric": psteros.CalculationOverride(parameters={"INCAR": dict(scf)}),
        f"{name}_passivated": psteros.CalculationOverride(
            parameters={"INCAR": {k: v for k, v in ASYMMETRIC.items() if k != "ISIF"} | scf}),
    }
    relax_extra = {k: v for k, v in relax_extra.items() if k in structures}
    static_extra = {k: v for k, v in static_extra.items() if k in structures}
    relax, static = z.recipes(name + (f"_{only}" if only else ""), study.potential_mapping(z.POTCARS),
                              relax_extra, static_extra, 24)
    graph = psteros.build_relax_static_workgraph(structures, relax, static, submit=do_submit)
    print(f"graph {name!r}: " + ", ".join(f"{k}:{len(v)}" for k, v in structures.items())
          + (f" -> PK {graph.pk}" if do_submit else " [not submitted]"))
    if do_submit:
        os.makedirs(OUT, exist_ok=True)
        runs = load_runs()
        runs[name + (f"_{only}" if only else "")] = {"pk": graph.pk, "labels": list(structures), "face": face,
                                                        "layers": layers}
        json.dump(runs, open(RUNS, "w"), indent=2)


def analyse():
    from aiida import orm
    from psteros.phase_diagram import SlabTermination
    from psteros.polar import nonpolar_check

    r = z.relaxed_inputs(True)
    bulk = r["bulk"][1]
    energies = {k: v[0] for k, v in r.items()}
    energies.update(json.load(open(os.path.join(HERE, "zno", "energies.json"))))
    cluster_energies = dict(energies)
    cluster_energies.update(json.load(open(os.path.join(HERE, "zno", "energies_clusters.json"))))
    studies = {m: z.build_study(bulk, r["ref_Zn"][1], r["ref_O"][1], m) for m in ("molecules", "clusters")}
    refs = studies["molecules"].binary_references(energies)
    hydrogen = {"molecules": studies["molecules"].pseudo_hydrogen_references(energies, refs),
                "clusters": studies["clusters"].pseudo_hydrogen_references(cluster_energies, refs)}

    rows = []
    # the pair of the main run: (10-10), 12 planes
    sym, pas = build_pair(bulk, FACES["100"], None)
    pairs = {("100", 12): (sym, pas, energies["ZnO_100_symmetric"], energies["ZnO_100_passivated"])}
    found = {}  # (face, layers) -> {label: energy}, from every graph that has the label finished
    for name, run in load_runs().items():
        graph = orm.load_node(run["pk"])
        for label in run["labels"]:
            try:
                energy = psteros.read_vasp_results(graph, [label])[0][label]
            except ValueError:
                try:
                    energy = z_read_from_tasks(graph, label)
                except Exception:
                    continue
            found.setdefault((run["face"], run["layers"]), {})[label] = energy
    for (face, layers), got in found.items():
        name = tag(face, layers)
        if f"{name}_symmetric" in got and f"{name}_passivated" in got:
            sym, pas = build_pair(bulk, FACES[face], layers)
            pairs[(face, layers)] = (sym, pas, got[f"{name}_symmetric"], got[f"{name}_passivated"])
    for (face, layers), (sym, pas, e_sym, e_pas) in sorted(pairs.items()):
        symmetric = SlabTermination.from_structure("sym", e_sym, sym.structure)
        passivated = SlabTermination.from_polar(pas, e_pas)
        thickness = max(s.coords @ (sym.structure.lattice.matrix[2] / sym.structure.lattice.c) for s in sym.structure) \
            - min(s.coords @ (sym.structure.lattice.matrix[2] / sym.structure.lattice.c) for s in sym.structure)
        out = {m: nonpolar_check(symmetric, passivated, refs, hydrogen[m]) for m in hydrogen}
        rows.append((face, layers, thickness, symmetric.surface_area_angstrom2, out))
    print(f"{'face':>5s} {'planes':>6s} {'thickness':>10s} {'area':>7s} | {'molecules':>12s} {'clusters':>12s}   (meV/A^2; per cell: x area)")
    for face, layers, thickness, area, out in rows:
        m, c = out["molecules"], out["clusters"]
        print(f"{face:>5s} {layers:6d} {float(thickness):9.2f}A {area:7.2f} | {m.difference_mev_per_angstrom2:+12.2f} "
              f"{c.difference_mev_per_angstrom2:+12.2f}   ({m.difference_mev_per_angstrom2 * area:+.1f} / "
              f"{c.difference_mev_per_angstrom2 * area:+.1f} meV)")
    return rows


def z_read_from_tasks(graph, label):
    from aiida.common.links import LinkType

    calls = {link.link_label: link.node for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all()}
    static, relax = calls[f"{label}_static_vasp"], calls[f"{label}_relax_vasp"]
    if static.exit_status != 0 or relax.exit_status != 0:
        raise ValueError(label)
    return static.outputs.misc.get_dict()["total_energies"]["energy_extrapolated"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("submit")
    p.add_argument("face", choices=FACES)
    p.add_argument("layers", type=int)
    p.add_argument("--submit", action="store_true")
    p.add_argument("--only", choices=("symmetric", "passivated"), help="just this slab of the pair, in its own graph")
    sub.add_parser("analyse")
    args = parser.parse_args()

    from aiida import load_profile

    load_profile(z.PROFILE)
    if args.cmd == "submit":
        submit(args.face, args.layers, args.submit, args.only)
    else:
        analyse()


if __name__ == "__main__":
    main()
