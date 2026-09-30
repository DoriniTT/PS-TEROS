"""SrTiO3(001) validation campaign: the analysis.

Builds the psteros ternary surface phase diagram from the three finished
graphs of ``campaign.py`` and compares every derived quantity with experiment
and literature.  Writes ``<out>/phase_diagram.png``, ``<out>/phase_diagram.csv``,
``<out>/validation.md`` and ``<out>/validation.json``.

    python analysis.py --profile P --refs-pk R --slabs-pk S --unrelaxed-pk U --out results
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import psteros

EV_PER_KJ_MOL = 1.0 / 96.4853
EV_A2_TO_J_M2 = psteros.EV_PER_ANGSTROM2_TO_J_PER_M2

# Experimental standard enthalpies of formation at 298 K (kJ/mol).
EXPERIMENT = {
    "SrTiO3": (-1668.99, 9.2, "bomb calorimetry, J. Chem. Thermodyn. (strontium titanates)"),
    "SrO": (-592.04, None, "NIST-JANAF (Chase 1998)"),
    "TiO2": (-944.0, 0.8, "rutile, CODATA (Cox, Wagman et al. 1984)"),
}
# Eglitis & Vanderbilt, PRB 77, 195408 (2008): B3PW, 7-layer slabs, eV per 1x1 surface cell.
EGLITIS_ENERGIES = {"E_cleav": 1.39, "E_rel": {"SrO": -0.24, "TiO2": -0.16}, "E_surf": {"SrO": 1.15, "TiO2": 1.23}}
# Same paper, Table III: rumpling s, Delta d12, Delta d23 in % of the lattice constant.
RELAXATION_REFERENCES = {
    "SrO": {"B3PW": (5.66, -6.58, 1.75), "LEED": (4.1, -5.0, 2.0), "RHEED": (4.1, 2.6, 1.3)},
    "TiO2": {"B3PW": (2.12, -5.79, 3.55), "LEED": (2.1, 1.0, -1.0), "RHEED": (2.6, 1.8, 1.3)},
}
SLABS = {"SrO": "slab_sro", "TiO2": "slab_tio2"}


def _child(graph, task_name):
    from aiida.common.links import LinkType

    for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        if link.link_label == task_name:
            return link.node
    return None


def _output(graph, name, task_name, port):
    """Graph-level output, or the output of the task's work chain when the graph attached none."""
    if name in graph.outputs:
        return graph.outputs[name]
    child = _child(graph, task_name)
    if child is None or port not in child.outputs:
        raise LookupError(f"{name!r} not found in graph {graph.pk}")
    return child.outputs[port]


def static_energy(graph, label) -> float:
    return float(_output(graph, f"{label}_static_parameters", f"{label}_static_qe", "output_parameters")["energy"])


def relaxed(graph, label):
    return _output(graph, f"{label}_relaxed_structure", f"{label}_relax_qe", "output_structure").get_pymatgen_structure()


def collect(refs, slabs, unrelaxed) -> dict:
    labels = ("srtio3", "sr_metal", "ti_metal", "sro", "tio2_rutile", "tio2_anatase", "o2")
    data = {label: {"energy": static_energy(refs, label), "structure": relaxed(refs, label)} for label in labels}
    for termination, label in SLABS.items():
        initial = _child(slabs, f"{label}_relax_qe").inputs.pw.structure.get_pymatgen_structure()
        data[label] = {
            "energy": static_energy(slabs, label),
            "structure": relaxed(slabs, label),
            "initial": initial,
            "unrelaxed_energy": float(_output(unrelaxed, f"{label}_parameters", f"{label}_qe", "output_parameters")["energy"]),
        }
    return data


def per_formula(data, label):
    structure = data[label]["structure"]
    units = structure.composition.get_reduced_composition_and_factor()[1]
    return data[label]["energy"] / units, structure.composition.reduced_formula


def formation_enthalpies(data) -> dict:
    e_sr = data["sr_metal"]["energy"] / len(data["sr_metal"]["structure"])
    e_ti = data["ti_metal"]["energy"] / len(data["ti_metal"]["structure"])
    e_o = data["o2"]["energy"] / 2.0
    refs = {"Sr": e_sr, "Ti": e_ti, "O": e_o}
    result = {}
    for label in ("srtio3", "sro", "tio2_rutile", "tio2_anatase"):
        energy, formula = per_formula(data, label)
        composition = data[label]["structure"].composition.reduced_composition
        result[label] = energy - sum(float(n) * refs[str(e)] for e, n in composition.items())
    return result


def relaxation_geometry(slab_data, a) -> tuple[float, float, float]:
    """Rumpling s and Delta d12, Delta d23 of the upper surface, in % of a (outward positive)."""
    initial, final = slab_data["initial"], slab_data["structure"]
    heights = sorted({round(site.coords[2], 3) for site in initial}, reverse=True)[:3]
    layers = []
    for height in heights:
        metal, oxygen = [], []
        for site_i, site_f in zip(initial, final):
            if abs(site_i.coords[2] - height) < 1e-2:
                (oxygen if site_i.specie.symbol == "O" else metal).append(site_f.coords[2] - site_i.coords[2])
        layers.append((sum(metal) / len(metal), sum(oxygen) / len(oxygen)))
    (m1, o1), (m2, _), (m3, _) = layers
    return 100 * (o1 - m1) / a, 100 * (m1 - m2) / a, 100 * (m2 - m3) / a


def analyze(data) -> tuple[dict, psteros.TernarySurfacePhaseDiagram]:
    bulk = data["srtio3"]["structure"]
    a = sum(bulk.lattice.abc) / 3.0
    e_bulk, _ = per_formula(data, "srtio3")
    enthalpies = formation_enthalpies(data)
    tio2 = min(("tio2_rutile", "tio2_anatase"), key=lambda label: per_formula(data, label)[0])
    reaction = e_bulk - per_formula(data, "sro")[0] - per_formula(data, tio2)[0]

    references = psteros.TernaryOxideReferences(
        bulk_energy_ev=data["srtio3"]["energy"],
        bulk_composition=bulk.composition,
        oxygen_molecule_energy_ev=data["o2"]["energy"],
        element_energies_per_atom_ev={
            "Sr": data["sr_metal"]["energy"] / len(data["sr_metal"]["structure"]),
            "Ti": data["ti_metal"]["energy"] / len(data["ti_metal"]["structure"]),
        },
        competing_phases=tuple(
            psteros.CompetingPhase(name, data[label]["energy"], data[label]["structure"].composition)
            for name, label in (("SrO", "sro"), ("TiO2 (rutile)", "tio2_rutile"), ("TiO2 (anatase)", "tio2_anatase"))
        ),
    )
    terminations = [
        psteros.SlabTermination.from_structure(termination, data[label]["energy"], data[label]["structure"])
        for termination, label in SLABS.items()
    ]
    diagram = psteros.ternary_surface_phase_diagram(terminations, references)

    # Stability strip in Delta mu_SrO = Delta mu_Sr + Delta mu_O; its width should be -Delta H(SrO + TiO2 -> SrTiO3).
    strip = [x + y for x, y in references.stability_region]
    # gamma_SrO - gamma_TiO2 depends only on Delta mu_SrO; find where they cross.
    (g1, s1, _), (g2, s2, _) = diagram.planes["SrO"], diagram.planes["TiO2"]
    crossing = (g2 - g1) / (s1 - s2) if abs(s1 - s2) > 1e-12 else None

    # Eglitis & Vanderbilt definitions (eV per 1x1 surface cell).
    sro, tio2_slab = data["slab_sro"], data["slab_tio2"]
    e_cleav = 0.25 * (sro["unrelaxed_energy"] + tio2_slab["unrelaxed_energy"] - 7 * e_bulk)
    e_rel = {t: 0.5 * (data[SLABS[t]]["energy"] - data[SLABS[t]]["unrelaxed_energy"]) for t in SLABS}
    e_surf = {t: e_cleav + e_rel[t] for t in SLABS}
    area = terminations[0].surface_area_angstrom2
    average_psteros = [
        0.5 * (diagram.gamma_ev_per_angstrom2("SrO", x, y) + diagram.gamma_ev_per_angstrom2("TiO2", x, y)) * area
        for x, y in references.stability_region
    ]

    results = {
        "lattice_constant_A": a,
        "formation_enthalpy_eV": enthalpies,
        "reaction_SrO_TiO2_to_SrTiO3_eV": reaction,
        "lowest_TiO2": tio2,
        "stability_region_vertices": references.stability_region,
        "stability_boundaries": references.stability_boundaries,
        "strip_delta_mu_SrO_eV": [min(strip), max(strip)],
        "termination_crossing_delta_mu_SrO_eV": crossing,
        "E_cleav_eV_per_cell": e_cleav,
        "E_rel_eV_per_cell": e_rel,
        "E_surf_eV_per_cell": e_surf,
        "average_surface_energy_J_m2": 0.5 * (e_surf["SrO"] + e_surf["TiO2"]) / area * EV_A2_TO_J_M2,
        "average_from_psteros_eV_per_cell": average_psteros,
        "relaxation_percent_of_a": {t: relaxation_geometry(data[SLABS[t]], a) for t in SLABS},
        "o2_bond_A": data["o2"]["structure"].get_distance(0, 1),
        "area_A2": area,
    }
    return results, diagram


def report(results) -> str:
    kj = {k: (v * EV_PER_KJ_MOL, None if u is None else u * EV_PER_KJ_MOL, s) for k, (v, u, s) in EXPERIMENT.items()}
    dh = results["formation_enthalpy_eV"]
    exp_reaction = kj["SrTiO3"][0] - kj["SrO"][0] - kj["TiO2"][0]
    rows = [
        "| Quantity | psteros (QE, PBE) | Reference | Source |",
        "|---|---|---|---|",
        f"| a(SrTiO3), Å | {results['lattice_constant_A']:.4f} | 3.89 | experiment extrapolated to 0 K, quoted by Eglitis & Vanderbilt |",
        f"| O2 bond, Å | {results['o2_bond_A']:.4f} | 1.23 | plane-wave PBE, Alexandrov et al. (arXiv:1005.4833) |",
    ]
    for label, name in (("srtio3", "SrTiO3"), ("sro", "SrO"), ("tio2_rutile", "TiO2")):
        value, unc, source = kj[name]
        rows.append(
            f"| ΔH_f({name}), eV/f.u. | {dh[label]:.3f} | {value:.3f}" + (f" ± {unc:.2f}" if unc else "") + f" | {source} |"
        )
    rows += [
        f"| ΔH_f(TiO2 anatase), eV/f.u. | {dh['tio2_anatase']:.3f} | — | lowest PBE polymorph: {results['lowest_TiO2']} |",
        f"| ΔH(SrO + TiO2 → SrTiO3), eV | {results['reaction_SrO_TiO2_to_SrTiO3_eV']:.3f} | {exp_reaction:.3f} ± 0.10 | from the three experimental ΔH_f above |",
        f"| Strip width in Δμ_SrO, eV | {results['strip_delta_mu_SrO_eV'][1] - results['strip_delta_mu_SrO_eV'][0]:.3f} | {-results['reaction_SrO_TiO2_to_SrTiO3_eV']:.3f} | must equal −ΔH above (internal) |",
        f"| E_cleav, eV/cell | {results['E_cleav_eV_per_cell']:.3f} | {EGLITIS_ENERGIES['E_cleav']:.2f} | B3PW, Eglitis & Vanderbilt PRB 77, 195408 |",
    ]
    for t in SLABS:
        rows.append(
            f"| E_rel({t}), eV/cell | {results['E_rel_eV_per_cell'][t]:.3f} | {EGLITIS_ENERGIES['E_rel'][t]:.2f} | B3PW, same |"
        )
        rows.append(
            f"| E_surf({t}), eV/cell | {results['E_surf_eV_per_cell'][t]:.3f} | {EGLITIS_ENERGIES['E_surf'][t]:.2f} | B3PW, same |"
        )
    average = 0.5 * sum(results["E_surf_eV_per_cell"].values())
    spread = max(abs(v - average) for v in results["average_from_psteros_eV_per_cell"])
    rows.append(
        f"| Average γ from psteros planes, eV/cell | {average:.4f} (max deviation {spread:.1e}) | {average:.4f} | must be μ-independent (internal) |"
    )
    for t in SLABS:
        s, d12, d23 = results["relaxation_percent_of_a"][t]
        ref = RELAXATION_REFERENCES[t]
        rows.append(
            f"| {t}: s, Δd12, Δd23 (% of a) | {s:.2f}, {d12:.2f}, {d23:.2f} | B3PW {ref['B3PW']}; LEED {ref['LEED']}; RHEED {ref['RHEED']} | Eglitis & Vanderbilt Table III |"
        )
    low, high = results["strip_delta_mu_SrO_eV"]
    crossing = results["termination_crossing_delta_mu_SrO_eV"]
    lines = ["# SrTiO3(001) validation", "", *rows, "",
             f"Stability region bounded by: {', '.join(b for b in results['stability_boundaries'] if b)}.",
             f"SrO/TiO2 terminations cross at Δμ_SrO = {crossing:.3f} eV; the strip spans {low:.3f} to {high:.3f} eV "
             f"({100 * (crossing - low) / (high - low):.0f}% from its TiO2-rich edge)." if crossing is not None else ""]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--refs-pk", type=int, required=True)
    parser.add_argument("--slabs-pk", type=int, required=True)
    parser.add_argument("--unrelaxed-pk", type=int, required=True)
    parser.add_argument("--out", default="results")
    args = parser.parse_args(argv)

    from aiida import load_profile, orm

    load_profile(args.profile)
    data = collect(*(orm.load_node(pk) for pk in (args.refs_pk, args.slabs_pk, args.unrelaxed_pk)))
    results, diagram = analyze(data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    diagram.plot(out / "phase_diagram.png", title="SrTiO$_3$(001), QE PBE: stable 1×1 terminations")
    diagram.to_csv(out / "phase_diagram.csv")
    (out / "validation.json").write_text(json.dumps(results, indent=2, default=list))
    text = report(results)
    (out / "validation.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
