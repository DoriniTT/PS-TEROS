"""Step C7 of docs/plans/MACE_SURFACE_PLAN.md: surface phase diagrams of every energy source and gates G0-G5.

Energy sources (each diagram takes all its energies from one source, psteros.check_one_energy_source):
    R0      CP2K statics of the terminations as built (c2_labels), CP2K references
    R1      CP2K relax -> static (c3_truth), CP2K references: the ground truth
    <M>     route (a): model M at its own geometries, for slabs, bulk RuO2, hcp Ru and O2 (c4_main report)
    <M>_dft route (b): CP2K static at the geometry of model M (c6_route_b), CP2K references
Models: m0 (MH-1 zero-shot) and m1, m2, m3 with seeds 20260920 (deployment seed; the gates use it, decided
before any result), 20260921 and 20260922 (the committee spread).

CP2K references of the diagrams: bulk110 (the bulk in the slabs' (110) cell, same in-plane k-mesh), hcp Ru and O2;
ΔH_f of G0 and G4 from the conventional bulk, as the model's.

Gates (PLAN.md): G0 CP2K against Reuter and Scheffler (o_bridge, o_cus, metal only); G1 same stable termination
except within 0.05 eV of an R1 transition, every R1 transition within ±0.05 eV; G2 |γ_M − γ_R1| <= 5 meV/A^2
(route a), <= 2 (route b); G3 E(x0) − E(x_M) >= 0.9 [E(x0) − E(x_R1)] where the R1 relaxation lowers γ by more
than 2 meV/A^2; G4 |ΔH_f(M) − ΔH_f(CP2K)| <= 0.05 eV/f.u.; G5 is R0 already passing G1 and G2?

Usage (PYTHONPATH=<worktree>):
    python c7_analysis.py              # writes c7_analysis.json and prints the gates
    python c7_analysis.py --dft-only   # CP2K only (G0, G5, R0 and R1 diagrams): c7_dft.json
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np

import psteros

HERE = Path(__file__).resolve().parent
TERMINATIONS = psteros.RUTILE_110_TERMINATIONS
LITERATURE = ("o_bridge", "o_cus", "metal")
DEPLOY = 20260920
SEEDS = (20260920, 20260921, 20260922)
REUTER = {"transition_ev": -0.93, "transition_tol": 0.15, "dgamma_mev": 49.0, "dgamma_tol": 10.0,
          "dhf_ev": -3.16, "dhf_tol": 0.20}


# ---------------------------------------------------------------- pure part (no AiiDA)

def diagram(slab_energies, slabs, bulk_energy, bulk_composition, o2_energy, ru_energy_per_atom, labels, window):
    """psteros surface phase diagram of ``labels`` on the Δμ_O ``window`` (both ends included)."""

    terminations = []
    for label in labels:
        structure = slabs[label]
        composition = {}
        for site in structure:
            composition[site.specie.symbol] = composition.get(site.specie.symbol, 0) + 1
        lattice = structure.lattice.matrix
        area = float(np.linalg.norm(np.cross(lattice[0], lattice[1])))
        terminations.append(psteros.SlabTermination(label, slab_energies[label], composition, area, surfaces=2))
    references = psteros.BinaryOxideReferences(
        bulk_energy_ev=bulk_energy, bulk_composition=bulk_composition, oxygen_molecule_energy_ev=o2_energy,
        metal_energy_per_atom_ev=ru_energy_per_atom)
    return psteros.surface_phase_diagram(terminations, references, delta_mu_range=window, points=601)


def gammas(result) -> dict[str, np.ndarray]:
    return {label: 1000 * np.array([p.gamma_ev_per_angstrom2 for p in points]) for label, points in result.curves.items()}


def gate_g1(model, truth, window_ev=0.05) -> dict:
    grid = np.array(truth.delta_mu_ev)
    near = np.zeros(len(grid), dtype=bool)
    for mu, *_ in truth.transitions:
        near |= np.abs(grid - mu) <= window_ev
    disagree = [float(mu) for mu, a, b, n in zip(grid, model.stable, truth.stable, near) if a != b and not n]
    matched = {}
    for mu, below, above in truth.transitions:
        partners = [m for m, b, a in model.transitions if (b, a) == (below, above)]
        matched[f"{below}->{above}"] = {"truth_ev": mu, "model_ev": min(partners, key=lambda m: abs(m - mu))
                                        if partners else None}
    ok_transitions = all(v["model_ev"] is not None and abs(v["model_ev"] - v["truth_ev"]) <= window_ev
                         for v in matched.values())
    return {"pass": not disagree and ok_transitions, "disagreeing_delta_mu_ev": disagree[:10],
            "n_disagreeing": len(disagree), "transitions": matched,
            "model_transitions": [list(t) for t in model.transitions]}


def gate_g2(model, truth, limit) -> dict:
    a, b = gammas(model), gammas(truth)
    worst = {label: float(np.abs(a[label] - b[label]).max()) for label in b}
    return {"pass": max(worst.values()) <= limit, "limit_mev": limit, "max_abs_error_mev": worst}


def gate_g3(e_x0, e_model, e_truth, areas, limit=0.9, threshold_mev=2.0) -> dict:
    rows, ok = {}, True
    for label in e_truth:
        relaxation = e_x0[label] - e_truth[label]
        lowering = 1000 * relaxation / (2 * areas[label])
        row = {"dft_relaxation_mev_per_angstrom2": lowering}
        if label in e_model:
            row["recovered_fraction"] = (e_x0[label] - e_model[label]) / relaxation if relaxation else None
        if lowering > threshold_mev:
            row["gated"] = True
            row["pass"] = label in e_model and row["recovered_fraction"] >= limit
            ok &= row["pass"]
        rows[label] = row
    return {"pass": ok, "terminations": rows}


def formation(bulk_cell, bulk_fu, ru_atom, o2) -> float:
    return bulk_cell / bulk_fu - ru_atom - o2


def gate_g0(truth) -> dict:
    transitions = [t for t in truth.transitions if (t[1], t[2]) == ("o_bridge", "o_cus")]
    g = gammas(truth)
    grid = np.array(truth.delta_mu_ev)
    rich = int(np.argmin(np.abs(grid - 0.0)))
    dgamma = float(g["o_bridge"][rich] - g["o_cus"][rich])
    dhf = truth.references.formation_enthalpy_ev
    window = [i for i, mu in enumerate(grid) if truth.in_stability_window(mu)]
    metal_stable = any(truth.stable[i] == "metal" for i in window)
    result = {
        "transition_ev": transitions[0][0] if transitions else None,
        "obridge_minus_ocus_at_0_mev": dgamma, "delta_h_f_ev": dhf, "metal_stable_somewhere": metal_stable,
    }
    result["pass"] = bool(
        transitions and abs(transitions[0][0] - REUTER["transition_ev"]) <= REUTER["transition_tol"]
        and abs(dgamma - REUTER["dgamma_mev"]) <= REUTER["dgamma_tol"]
        and abs(dhf - REUTER["dhf_ev"]) <= REUTER["dhf_tol"] and not metal_stable)
    return result


# ---------------------------------------------------------------- data (AiiDA)

def load(name):
    from aiida import orm

    return orm.load_node(json.loads((HERE / name).read_text())[
        "shelljob_pk" if name.startswith("c4") else "workgraph_pk"])


def energies(node, labels):
    return psteros.read_cp2k_results(node, list(labels))[0]


def report(path):
    return load("c4_main.json").outputs.report.base.repository.get_object_content(path)


def main(dft_only: bool = False) -> None:
    refs, labels, truth_graph = load("c2_refs.json"), load("c2_labels.json"), load("c3_truth.json")
    from psteros.backends.cp2k import _stage_output

    slabs = {t: _stage_output(labels, t, "retrieved").creator.inputs.structure.get_pymatgen() for t in TERMINATIONS}
    ref = energies(refs, ("bulk", "ru", "o2"))
    bulk110 = energies(labels, ("bulk110",))["bulk110"]
    dft_refs = dict(bulk_energy=bulk110, bulk_composition={"Ru": 2, "O": 4}, o2_energy=ref["o2"],
                    ru_energy_per_atom=ref["ru"] / 2)
    dhf_dft = formation(ref["bulk"], 2, ref["ru"] / 2, ref["o2"])
    window = (dhf_dft / 2, 0.0)
    e_x0 = energies(labels, TERMINATIONS)
    e_r1 = energies(truth_graph, TERMINATIONS)
    out = {"window_ev": window, "delta_h_f_cp2k_ev": dhf_dft,
           "delta_h_f_cp2k_bulk110_ev": formation(bulk110, 2, ref["ru"] / 2, ref["o2"])}
    r1 = diagram(e_r1, slabs, labels=TERMINATIONS, window=window, **dft_refs)
    r0 = diagram(e_x0, slabs, labels=TERMINATIONS, window=window, **dft_refs)
    out["R1"] = {"transitions": [list(t) for t in r1.transitions], "gamma_at_0_mev": {k: float(v[-1]) for k, v in gammas(r1).items()}}
    out["G0"] = gate_g0(diagram(e_r1, slabs, labels=LITERATURE, window=(-1.8, 0.0), **dft_refs))
    out["G5"] = {"R0_G1": gate_g1(r0, r1), "R0_G2": gate_g2(r0, r1, 5.0)}
    out["G5"]["discriminating"] = not (out["G5"]["R0_G1"]["pass"] and out["G5"]["R0_G2"]["pass"])
    out["R0"] = {"transitions": [list(t) for t in r0.transitions],
                 "gamma_at_0_mev": {k: float(v[-1]) for k, v in gammas(r0).items()}}
    areas = {t: float(np.linalg.norm(np.cross(*slabs[t].lattice.matrix[:2]))) for t in TERMINATIONS}
    out["relaxation_lowers_gamma_mev"] = {t: 1000 * (e_x0[t] - e_r1[t]) / (2 * areas[t]) for t in TERMINATIONS}
    if dft_only:
        (HERE / "c7_dft.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
        print(json.dumps({k: out[k] for k in ("window_ev", "delta_h_f_cp2k_ev", "delta_h_f_cp2k_bulk110_ev", "G0", "R1",
                                                 "R0", "relaxation_lowers_gamma_mev")}, indent=1, default=float))
        print("G5 discriminating:", out["G5"]["discriminating"], "| R0 vs R1 G1:", out["G5"]["R0_G1"]["pass"],
              "G2:", out["G5"]["R0_G2"]["pass"], out["G5"]["R0_G2"]["max_abs_error_mev"])
        return

    results = json.loads(report("relax/results.json"))["models"]
    route_b = json.loads((HERE / "c6_route_b.json").read_text()) if (HERE / "c6_route_b.json").exists() else None
    b_node = load("c6_route_b.json") if route_b else None
    out["models"] = {}
    for model, data in results.items():
        rows = data["structures"]
        converged = all(rows[t]["converged"] for t in TERMINATIONS)
        e = {label: rows[label]["energy_ev"] for label in rows}
        entry = {"all_slabs_converged": converged,
                 "delta_h_f_ev": formation(e["bulk"], 2, e["ru"] / 2, e["o2"])}
        entry["G4"] = {"pass": abs(entry["delta_h_f_ev"] - dhf_dft) <= 0.05,
                       "difference_ev": entry["delta_h_f_ev"] - dhf_dft}
        route_a = diagram({t: e[t] for t in TERMINATIONS}, slabs, bulk_energy=e["bulk"],
                          bulk_composition={"Ru": 2, "O": 4}, o2_energy=e["o2"], ru_energy_per_atom=e["ru"] / 2,
                          labels=TERMINATIONS, window=window)
        entry["route_a"] = {"G1": gate_g1(route_a, r1), "G2": gate_g2(route_a, r1, 5.0),
                            "gamma_at_0_mev": {k: float(v[-1]) for k, v in gammas(route_a).items()}}
        if route_b and all(f"{model}_{t}" in route_b["labels"] for t in TERMINATIONS):
            e_b = energies(b_node, [f"{model}_{t}" for t in TERMINATIONS])
            e_b = {t: e_b[f"{model}_{t}"] for t in TERMINATIONS}
            diagram_b = diagram(e_b, slabs, labels=TERMINATIONS, window=window, **dft_refs)
            entry["route_b"] = {"G1": gate_g1(diagram_b, r1), "G2": gate_g2(diagram_b, r1, 2.0),
                                "G3": gate_g3(e_x0, e_b, e_r1, areas)}
        out["models"][model] = entry
    for family in ("m1", "m2", "m3"):
        members = [f"{family}_seed{seed}" for seed in SEEDS if f"{family}_seed{seed}" in out["models"]]
        if members:
            curves = np.array([[out["models"][m]["route_a"]["gamma_at_0_mev"][t] for t in TERMINATIONS] for m in members])
            out[f"{family}_committee_spread_at_0_mev"] = dict(zip(TERMINATIONS, curves.std(axis=0).tolist()))
    (HERE / "c7_analysis.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
    summary = {"G0": out["G0"]["pass"], "G5_discriminating": out["G5"]["discriminating"]}
    for model, entry in out["models"].items():
        summary[model] = {"G4": entry["G4"]["pass"], "a_G1": entry["route_a"]["G1"]["pass"],
                          "a_G2": entry["route_a"]["G2"]["pass"]}
        if "route_b" in entry:
            summary[model].update(b_G1=entry["route_b"]["G1"]["pass"], b_G2=entry["route_b"]["G2"]["pass"],
                                  b_G3=entry["route_b"]["G3"]["pass"])
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    from aiida import load_profile

    import sys

    import cluster

    load_profile("psteros_mace_ruo2")
    cluster.with_retry(main, "--dft-only" in sys.argv)
