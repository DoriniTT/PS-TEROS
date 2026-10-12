"""Analysis of the SrTiO3(001) ternary smoke test: free-energy ternary diagrams, tables, frequencies and checks T1-T8.

    python analyse_ternary.py --refs-pk R --slabs-pk S --vib-pk V [--rerun-pk ...] \\
        [--o2-refs-pk 977 --o2-vib-pk 3400] [--legacy-src <checkout of 8b74134>]

Writes to ``results/``: ``ternary_total_energy.{png,csv}`` (static energies, as before the feature),
``ternary_T<T>K.{png,csv}`` (E + F_vib(T), T = 0, 300, 600, 900 K), ``ternary_comparison.png``, ``ternary_table.csv``,
``frequencies_<label>.csv``, ``mode_summary.csv`` and ``checks.json``.

O2 is not part of this test: its energy and zero-point energy come from the SnO2 test (same code and settings).
The calculated numbers go through the public ``psteros`` API; the checks reuse the label-independent readers of
``validation/vibrations_lovelace/analyse.py``.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import psteros

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE.parent / "vibrations_lovelace"))
import analyse as sno2_analysis  # noqa: E402  (ModeGraph, check_displacements, check_independent_recomputation)
import run_test as sno2  # noqa: E402

BULKS = ("srtio3_bulk", "sr_metal", "ti_metal", "sro_bulk", "tio2_bulk")
SLABS = ("slab_tio2", "slab_sro")
LABELS = BULKS + SLABS
TEMPERATURES = (0.0, 300.0, 600.0, 900.0)
EXPECTED_MODES = {
    "srtio3_bulk": 27, "tio2_bulk": 15, "sro_bulk": 21, "sr_metal": 9, "ti_metal": 9, "slab_tio2": 24, "slab_sro": 21,
}
EXPECTED_ZERO_MODES = {"srtio3_bulk": 3, "tio2_bulk": 3, "sro_bulk": 3, "sr_metal": 3, "ti_metal": 3, "slab_tio2": 0, "slab_sro": 0}


def read_energies(refs_graph, slabs_graph, o2_refs_graph):
    energies, relaxed = psteros.read_vasp_results(refs_graph, BULKS)
    slab_energies, slab_structures = psteros.read_vasp_results(slabs_graph, SLABS)
    o2_energy, _ = psteros.read_vasp_results(o2_refs_graph, ("o2",))
    energies.update(slab_energies)
    energies.update(o2_energy)
    relaxed.update(slab_structures)
    return energies, relaxed


def references_and_slabs(energies, relaxed, modes=None, temperature_k=None):
    """Ternary references and terminations; with ``modes`` every energy becomes a free energy (O2: E + ZPE)."""

    def free(label):
        if modes is None:
            return energies[label]
        return psteros.solid_free_energy_ev(energies[label], modes[label], temperature_k)

    structure = {label: relaxed[label].get_pymatgen_structure() for label in LABELS}
    count = lambda label, element: structure[label].composition[element]
    o2 = energies["o2"] if modes is None else psteros.molecule_reference_energy_ev(energies["o2"], modes["o2"])
    references = psteros.TernaryOxideReferences(
        bulk_energy_ev=free("srtio3_bulk"),
        bulk_composition=structure["srtio3_bulk"].composition,
        oxygen_molecule_energy_ev=o2,
        element_energies_per_atom_ev={
            "Sr": free("sr_metal") / count("sr_metal", "Sr"), "Ti": free("ti_metal") / count("ti_metal", "Ti"),
        },
        competing_phases=(
            psteros.CompetingPhase("SrO", free("sro_bulk"), structure["sro_bulk"].composition),
            psteros.CompetingPhase("TiO2", free("tio2_bulk"), structure["tio2_bulk"].composition),
        ),
    )
    bulk = None if modes is None else modes["srtio3_bulk"]
    terminations = [
        psteros.SlabTermination.from_structure(
            label,
            energies[label] if modes is None else psteros.solid_free_energy_ev(
                energies[label], modes[label], temperature_k, bulk=bulk),
            relaxed[label],
        )
        for label in SLABS
    ]
    return references, terminations


def polygon_area(polygon) -> float:
    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(polygon, polygon[1:] + polygon[:1]))) if polygon else 0.0


def diagrams(energies, relaxed, modes):
    RESULTS.mkdir(exist_ok=True)
    out, rows = {}, []
    for stem, temperature in [("total_energy", None)] + [(f"T{int(t)}K", t) for t in TEMPERATURES]:
        if temperature is not None and modes is None:
            continue
        references, terminations = references_and_slabs(energies, relaxed, modes if temperature is not None else None, temperature)
        diagram = psteros.ternary_surface_phase_diagram(terminations, references)
        out[stem] = diagram
        title = "SrTiO$_3$(001), static total energies" if temperature is None else (
            f"SrTiO$_3$(001), E + F$_{{vib}}$ at {temperature:g} K" + (" (ZPE only)" if temperature == 0 else ""))
        diagram.plot(RESULTS / f"ternary_{stem}.png", title=title)
        diagram.to_csv(RESULTS / f"ternary_{stem}.csv")
        total = polygon_area(list(references.stability_region))
        row = {
            "case": stem, "temperature_K": "" if temperature is None else temperature,
            "formation_free_energy_eV": references.formation_enthalpy_ev,
            "stability_region_area_eV2": total,
        }
        for label in SLABS:
            row[f"region_fraction_{label}"] = polygon_area(list(diagram.regions[label])) / total if total else float("nan")
        # surface energies at the centroid of the stability region
        xs = [p[0] for p in references.stability_region]
        ys = [p[1] for p in references.stability_region]
        point = (float(np.mean(xs)), float(np.mean(ys)))
        row["centroid_dmu_Sr_eV"], row["centroid_dmu_O_eV"] = point
        for label in SLABS:
            row[f"gamma_{label}_at_centroid_J_per_m2"] = (
                diagram.gamma_ev_per_angstrom2(label, *point) * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2)
        row["stable_at_centroid"] = diagram.stable_termination(*point)
        rows.append(row)
    sno2_analysis.write_csv(RESULTS / "ternary_table.csv", list(rows[0]), [list(r.values()) for r in rows])
    return out, rows


def comparison_figure(temperature="T900K"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    for axis, (stem, title) in zip(axes, (("total_energy", "Without vibrations: static total energies"),
                                          (temperature, f"With vibrations: E + F$_{{vib}}$ at {temperature[1:-1]}"))):
        axis.imshow(mpimg.imread(RESULTS / f"ternary_{stem}.png"))
        axis.set_title(title, fontsize=10, loc="left")
        axis.axis("off")
    fig.savefig(RESULTS / "ternary_comparison.png", dpi=150)
    plt.close(fig)


def frequency_tables(mode_graph):
    summary = []
    for label in LABELS:
        data = mode_graph.outputs[f"{label}_vibrations"].get_dict()
        frequencies = sorted(data["frequencies_cm1"])
        sno2_analysis.write_csv(RESULTS / f"frequencies_{label}.csv", ["mode", "frequency_cm1"],
                                [[i, f"{v:.4f}"] for i, v in enumerate(frequencies)])
        modes = psteros.HarmonicVibrations.from_dict(data, imaginary_modes="drop")
        summary.append({
            "label": label, "modes": len(frequencies), "zero_modes_removed": data["zero_modes"],
            "displaced_sites": len(data["displaced_sites"]), "supercell_size": data["supercell_size"],
            "lowest_cm1": frequencies[0], "highest_cm1": frequencies[-1],
            "imaginary_modes": len(modes.imaginary_frequencies_cm1),
            "lowest_imaginary_cm1": -min(frequencies) if frequencies[0] < 0 else 0.0,
            "zpe_eV": modes.zero_point_energy_ev, "f_vib_300K_eV": modes.free_energy_ev(300.0),
            "f_vib_900K_eV": modes.free_energy_ev(900.0),
        })
    sno2_analysis.write_csv(RESULTS / "mode_summary.csv", list(summary[0]), [list(r.values()) for r in summary])
    return summary


# ------------------------------------------------------------------------------------------ checks


def check_mode_counts(summary):
    counts = {r["label"]: r["modes"] for r in summary}
    zero = {r["label"]: r["zero_modes_removed"] for r in summary}
    return counts == EXPECTED_MODES and zero == EXPECTED_ZERO_MODES, {"modes": counts, "zero_modes": zero, "expected": EXPECTED_MODES}


def check_frozen_centre(mode_graph):
    """The frozen SrTiO3 unit of each slab is counted as one bulk cell; a non-stoichiometric region is refused."""

    modes = psteros.read_vibrations(mode_graph, LABELS, imaginary_modes="drop")
    bulk = modes["srtio3_bulk"]
    detail, passed = {}, True
    for label in SLABS:
        vib = modes[label]
        added = psteros.solid_free_energy_ev(0.0, vib, 600.0, bulk=bulk) - vib.free_energy_ev(600.0)
        cells = added / bulk.free_energy_ev(600.0)
        detail[label] = {"frozen": dict(vib.frozen_composition), "bulk_cells_counted": round(cells, 9)}
        passed &= abs(cells - 1.0) < 1e-9 and dict(vib.frozen_composition) == {"Sr": 1, "Ti": 1, "O": 3}
    broken = psteros.HarmonicVibrations(
        modes["slab_tio2"].frequencies_cm1, composition={"Sr": 2, "Ti": 3, "O": 8}, frozen_composition={"Sr": 1, "Ti": 2, "O": 3})
    try:
        psteros.solid_free_energy_ev(0.0, broken, 600.0, bulk=bulk)
        detail["non_stoichiometric_frozen_region"] = "NOT refused"
        passed = False
    except ValueError as error:
        detail["non_stoichiometric_frozen_region"] = "refused: " + str(error)[:70]
    return passed, detail


def check_imaginary(mode_graph, summary):
    detail, passed = {}, True
    for row in summary:
        label = row["label"]
        entry = {"imaginary_modes": row["imaginary_modes"], "lowest_imaginary_cm1": row["lowest_imaginary_cm1"]}
        try:
            psteros.read_vibrations(mode_graph, [label])
            raised = False
        except ValueError:
            raised = True
        dropped = psteros.read_vibrations(mode_graph, [label], imaginary_modes="drop")[label]
        entry["read_default"] = "raises ValueError" if raised else "ok"
        entry["modes_after_drop"] = len(dropped.mode_energies_ev)
        passed &= raised == bool(row["imaginary_modes"])
        if row["imaginary_modes"]:
            passed &= row["lowest_imaginary_cm1"] < 200.0
        detail[label] = entry
    return passed, detail


def check_backward_compatibility(energies, relaxed, legacy_src):
    payload = {"energies": energies, "bulk": {}, "slabs": {}}
    for label in BULKS:
        composition = relaxed[label].get_pymatgen_structure().composition.get_el_amt_dict()
        payload["bulk"][label] = {"composition": dict(composition)}
    for label in SLABS:
        structure = relaxed[label].get_pymatgen_structure()
        payload["slabs"][label] = {
            "composition": dict(structure.composition.get_el_amt_dict()),
            "area": float(np.linalg.norm(np.cross(*structure.lattice.matrix[:2]))),
        }
    outputs = {}
    for name, source in (("legacy", Path(legacy_src)), ("current", REPOSITORY)):
        target = RESULTS / f"_{name}_ternary.csv"
        subprocess.run([sys.executable, str(HERE / "legacy_ternary.py"), str(target)], input=json.dumps(payload).encode(),
                       env={**os.environ, "PYTHONPATH": str(source)}, check=True, capture_output=True)
        outputs[name] = target.read_bytes()
        target.unlink()
    ours = (RESULTS / "ternary_total_energy.csv").read_bytes()
    detail = {"legacy_src": str(legacy_src), "identical_to_pre_feature_code": outputs["legacy"] == outputs["current"],
              "legacy_equals_analysis_csv": outputs["legacy"] == ours}
    return all(v for k, v in detail.items() if k != "legacy_src"), detail


def check_sanity(rows):
    by_case = {r["case"]: r for r in rows}
    static = by_case["total_energy"]
    detail = {"formation_free_energy_eV": {c: round(r["formation_free_energy_eV"], 4) for c, r in by_case.items()},
              "stability_region_area_eV2": {c: round(r["stability_region_area_eV2"], 4) for c, r in by_case.items()},
              "stable_at_centroid": {c: r["stable_at_centroid"] for c, r in by_case.items()}}
    shift = max(abs(r["formation_free_energy_eV"] - static["formation_free_energy_eV"]) for r in rows)
    detail["max_shift_of_formation_energy_eV"] = round(shift, 4)
    passed = -19.0 < static["formation_free_energy_eV"] < -15.0 and shift < 0.3 and all(
        r["stability_region_area_eV2"] > 0 for r in rows)
    return passed, detail


def expected_incar(label):
    return {**sno2.ELECTRONIC, **sno2.STATIC, "EDIFF": sno2.VIBRATIONS_EDIFF}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default=sno2.PROFILE)
    parser.add_argument("--refs-pk", type=int, required=True)
    parser.add_argument("--slabs-pk", type=int, required=True)
    parser.add_argument("--vib-pk", type=int)
    parser.add_argument("--rerun-pk", type=int, nargs="*", default=[])
    parser.add_argument("--o2-refs-pk", type=int, default=977, help="graph with the O2 energy (SnO2 test)")
    parser.add_argument("--o2-vib-pk", type=int, default=3400, help="graph with the O2 vibrations (SnO2 test)")
    parser.add_argument("--legacy-src")
    args = parser.parse_args(argv)

    from aiida import load_profile, orm

    load_profile(args.profile)
    RESULTS.mkdir(exist_ok=True)
    energies, relaxed = read_energies(orm.load_node(args.refs_pk), orm.load_node(args.slabs_pk), orm.load_node(args.o2_refs_pk))
    print("static energies (eV):", {k: round(v, 5) for k, v in energies.items()})
    modes = mode_graph = None
    checks = {}
    if args.vib_pk:
        vib_graph = orm.load_node(args.vib_pk)
        mode_graph = sno2_analysis.ModeGraph(vib_graph, [orm.load_node(pk) for pk in args.rerun_pk])
        print("mode results from:", mode_graph.source)
        modes = psteros.read_vibrations(mode_graph, LABELS, imaginary_modes="drop")
        o2_modes = psteros.read_vibrations(sno2_analysis.ModeGraph(orm.load_node(args.o2_vib_pk)), ["o2"])
        modes.update(o2_modes)
    cases, rows = diagrams(energies, relaxed, modes)
    if args.legacy_src:
        checks["T7_backward_compatibility"] = check_backward_compatibility(energies, relaxed, args.legacy_src)
    if modes is not None:
        comparison_figure()
        summary = frequency_tables(mode_graph)
        checks["T2_displacements"] = sno2_analysis.check_displacements(mode_graph, relaxed, expected_incar)
        checks["T3_mode_counts"] = check_mode_counts(summary)
        checks["T4_frozen_centre"] = check_frozen_centre(mode_graph)
        checks["T5_imaginary_modes"] = check_imaginary(mode_graph, summary)
        checks["T6_independent_recomputation"] = sno2_analysis.check_independent_recomputation(mode_graph, "slab_tio2")
        checks["T8_physical_sanity"] = check_sanity(rows)
    for name, (passed, detail) in checks.items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}: {json.dumps(detail, default=float)[:700]}")
    (RESULTS / "checks.json").write_text(json.dumps(
        {n: {"passed": bool(p), "detail": d} for n, (p, d) in checks.items()}, indent=2, default=float))
    for row in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
    return 0 if all(p for p, _ in checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
