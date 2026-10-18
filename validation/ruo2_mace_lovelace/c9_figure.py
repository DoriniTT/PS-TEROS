"""Figure of the RuO2(110) validation for the docs: ``docs/images/mlip_ruo2_validation.png`` (and .svg).

(a) γ(Δμ_O) of the six terminations: CP2K relaxations (R1) and the MLIP route (M2R: MACE geometry, one CP2K static).
(b) The MLIP route against Reuter and Scheffler, PRB 65, 035406 (2002), Fig. 3 (FP-LAPW, GGA, relaxed slabs), whose
    four lines are read from the vector graphics of arXiv:cond-mat/0107229 (``REUTER``; straight lines to 0.05
    meV/A^2). Their (2x1) half-covered O-cus line is our o_cus_half. The 2007 erratum (PRB 75, 049901) corrects a
    statement on the temperature dependence of the formation energy, not these lines.
(c) Largest γ error of each termination for M2R: the model's own energies (route a) against the CP2K static at the
    model geometry (route b).
(d) Share of the CP2K relaxation energy that the MACE relaxation recovers, before (M2) and after (M2R) the
    active-learning round.

Usage (PYTHONPATH=<worktree>):
    python c9_figure.py data     # reads the AiiDA graphs and c7_analysis.json, writes c9_figure.json
    python c9_figure.py plot     # draws the figure from c9_figure.json (no AiiDA)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "c9_figure.json"
OUT = HERE.parents[1] / "docs" / "images" / "mlip_ruo2_validation"
TERMINATIONS = ("o_bridge", "o_cus", "metal", "o_cus_half", "o_bridge_half", "o_cus_no_bridge")
NAMES = {"o_bridge": "O$^\\mathrm{br}$ (stoich.)", "o_cus": "O$^\\mathrm{br}$ + O$^\\mathrm{cus}$",
         "metal": "no O$^\\mathrm{br}$", "o_cus_half": "O$^\\mathrm{br}$ + ½ O$^\\mathrm{cus}$",
         "o_bridge_half": "½ O$^\\mathrm{br}$", "o_cus_no_bridge": "O$^\\mathrm{cus}$, no O$^\\mathrm{br}$"}
# Reuter and Scheffler (2002), Fig. 3: gamma = intercept + slope * Delta mu_O (meV/A^2, eV), read from the
# vector paths of page 7 of arXiv:cond-mat/0107229 (axis ticks: x -2.0 ... 0.5 eV, y 0 ... 300 meV/A^2).
REUTER = {"o_bridge": (79.88, 0.0), "o_cus": (35.79, -46.776), "metal": (200.85, 46.775), "o_cus_half": (58.24, -23.387)}
COLORS = {"o_bridge": "#2f6db5", "o_cus": "#c9721c", "metal": "#7b4fa0", "o_cus_half": "#2e8b57",
          "o_bridge_half": "#b83232", "o_cus_no_bridge": "#6b7280"}


def data() -> None:
    import c7_analysis as c7
    from psteros.backends.cp2k import _stage_output

    refs, labels, truth = c7.load("c2_refs.json"), c7.load("c2_labels.json"), c7.load("c3_truth.json")
    slabs = {t: _stage_output(labels, t, "retrieved").creator.inputs.structure.get_pymatgen() for t in TERMINATIONS}
    ref = c7.energies(refs, ("bulk", "ru", "o2"))
    bulk110 = c7.energies(labels, ("bulk110",))["bulk110"]
    dft_refs = dict(bulk_energy=bulk110, bulk_composition={"Ru": 2, "O": 4}, o2_energy=ref["o2"],
                    ru_energy_per_atom=ref["ru"] / 2)
    window = (c7.formation(ref["bulk"], 2, ref["ru"] / 2, ref["o2"]) / 2, 0.0)
    r1 = c7.diagram(c7.energies(truth, TERMINATIONS), slabs, labels=TERMINATIONS, window=window, **dft_refs)
    m2r = c7.load("c6_route_b_m2r.json")
    e_b = c7.energies(m2r, [f"m2r_seed20260920_{t}" for t in TERMINATIONS])
    route_b = c7.diagram({t: e_b[f"m2r_seed20260920_{t}"] for t in TERMINATIONS}, slabs, labels=TERMINATIONS,
                         window=window, **dft_refs)
    analysis = json.loads((HERE / "c7_analysis.json").read_text())
    models = analysis["models"]
    out = {
        "window_ev": list(window),
        "delta_mu_ev": list(map(float, r1.delta_mu_ev)),
        "r1_mev": {t: v.tolist() for t, v in c7.gammas(r1).items()},
        "m2r_b_mev": {t: c7.on_grid(route_b, np.array(r1.delta_mu_ev))[t].tolist() for t in TERMINATIONS},
        "r1_transitions": [list(t) for t in r1.transitions],
        "m2r_b_transitions": [list(t) for t in route_b.transitions],
        "error_route_a_mev": models["m2r_seed20260920"]["route_a"]["G2"]["max_abs_error_mev"],
        "error_route_b_mev": models["m2r_seed20260920"]["route_b"]["G2"]["max_abs_error_mev"],
        "recovered_m2": {t: v["recovered_fraction"] for t, v in models["m2_seed20260920"]["route_b"]["G3"]["terminations"].items()},
        "recovered_m2r": {t: v["recovered_fraction"] for t, v in models["m2r_seed20260920"]["route_b"]["G3"]["terminations"].items()},
    }
    DATA.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", DATA.name, "R1 transitions", out["r1_transitions"], "M2R", out["m2r_b_transitions"])


def _strip(ax, mu, curves, labels, row, height, text):
    """Stable termination along mu as a coloured band in axes coordinates ``row`` (bottom) to ``row + height``."""

    table = np.array([curves[t] for t in labels])
    stable = table.argmin(axis=0)
    edges = [0, *(i for i in range(1, len(mu)) if stable[i] != stable[i - 1]), len(mu) - 1]
    for lo, hi in zip(edges, edges[1:]):
        label = labels[stable[lo]]
        ax.axvspan(mu[lo], mu[hi], ymin=row, ymax=row + height, color=COLORS[label], alpha=0.85, linewidth=0)
    ax.text(mu[0] + 0.02, row + height / 2, text, transform=ax.get_xaxis_transform(), fontsize=8, color="white",
            fontweight="bold", va="center")
    return [(float(mu[i]), labels[stable[i - 1]], labels[stable[i]]) for i in edges[1:-1]]


def plot() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    d = json.loads(DATA.read_text())
    mu = np.array(d["delta_mu_ev"])
    fig = plt.figure(figsize=(13.6, 11.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.3, 1.0], hspace=0.5, wspace=0.27)
    dashed = (0, (3, 2.5))

    # (a) CP2K relaxations against the MLIP route
    ax = fig.add_subplot(grid[0, 0])
    for t in TERMINATIONS:
        ax.plot(mu, d["r1_mev"][t], color=COLORS[t], linewidth=2.6, alpha=0.85, label=NAMES[t])
        ax.plot(mu, d["m2r_b_mev"][t], color="black", linewidth=1.0, linestyle=dashed)
    for value, *_ in d["r1_transitions"]:
        ax.axvline(value, color="#9aa5b1", linewidth=0.8, zorder=0)
    ax.set_xlim(mu[0], mu[-1])
    ax.set_ylim(20, 252)
    ax.set_xlabel(r"$\Delta\mu_\mathrm{O}$ (eV)")
    ax.set_ylabel(r"$\gamma$ (meV/Å$^2$)")
    ax.text(0.02, 0.97, "every termination within 0.5 meV/Å$^2$\nphase boundaries within 0.016 eV",
            transform=ax.transAxes, ha="left", va="top", fontsize=9, color="#3e4c59")
    top = make_axes_locatable(ax).append_axes("top", size="10%", pad=0.06, sharex=ax)
    _strip(top, mu, d["r1_mev"], TERMINATIONS, 0.5, 0.5, "CP2K relaxations")
    _strip(top, mu, d["m2r_b_mev"], TERMINATIONS, 0.0, 0.48, "MLIP route")
    top.set_yticks([])
    top.tick_params(axis="x", labelbottom=False, length=0)
    top.set_title("(a) CP2K relaxations (solid) and the MLIP route (dashed)", loc="left",
                  fontsize=11, fontweight="bold", pad=8)

    # (b) the MLIP route against Reuter and Scheffler
    bx = fig.add_subplot(grid[0, 1], sharey=ax)
    reuter = {t: intercept + slope * mu for t, (intercept, slope) in REUTER.items()}
    for t in TERMINATIONS:
        bx.plot(mu, d["m2r_b_mev"][t], color=COLORS[t], linewidth=2.4, alpha=0.85 if t in REUTER else 0.35)
    for t, values in reuter.items():
        bx.plot(mu, values, color=COLORS[t], linewidth=1.6, linestyle=(0, (1, 1.6)))
        bx.plot(mu[::60], values[::60], "o", color=COLORS[t], markersize=4.2, markerfacecolor="white", markeredgewidth=1.2)
    bx.set_xlim(mu[0], mu[-1])
    bx.set_xlabel(r"$\Delta\mu_\mathrm{O}$ (eV)")
    bx.tick_params(axis="y", labelleft=False)
    ours = d["m2r_b_transitions"]
    offsets = [REUTER[t][0] - d["m2r_b_mev"][t][-1] for t in REUTER]  # at Delta mu_O = 0
    theirs = [(t, a, b) for t, a, b in _strip_transitions(mu, reuter)]
    bx.text(0.02, 0.975, "Reuter and Scheffler: O$^\\mathrm{br}$ → O$^\\mathrm{br}$ + O$^\\mathrm{cus}$ at "
            f"{theirs[0][0]:.2f} eV\nMLIP route: {ours[0][0]:.2f} and {ours[-1][0]:.2f} eV (narrow ½ O$^\\mathrm{{cus}}$ window)\n"
            "no O$^\\mathrm{br}$: never stable in either\n"
            f"γ of the MLIP route {min(offsets):.0f}–{max(offsets):.0f} meV/Å$^2$ lower\n"
            "(their stated accuracy: 10 meV/Å$^2$)",
            transform=bx.transAxes, ha="left", va="top", fontsize=8.8, color="#3e4c59",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    top = make_axes_locatable(bx).append_axes("top", size="10%", pad=0.06, sharex=bx)
    _strip(top, mu, reuter, tuple(REUTER), 0.5, 0.5, "Reuter and Scheffler (2002)")
    _strip(top, mu, d["m2r_b_mev"], TERMINATIONS, 0.0, 0.48, "MLIP route")
    top.set_yticks([])
    top.tick_params(axis="x", labelbottom=False, length=0)
    top.set_title("(b) The MLIP route (solid) and Reuter and Scheffler 2002 (dotted)", loc="left",
                  fontsize=11, fontweight="bold", pad=8)
    handles = [plt.Line2D([], [], color=COLORS[t], linewidth=2.6) for t in TERMINATIONS]
    names = [NAMES[t] for t in TERMINATIONS]
    handles += [plt.Line2D([], [], color="black", linewidth=1.0, linestyle=dashed),
                plt.Line2D([], [], color="black", linewidth=1.6, linestyle=(0, (1, 1.6)), marker="o",
                           markerfacecolor="white", markersize=4.2)]
    names += ["MLIP route (in a)", "Reuter and Scheffler (in b)"]
    below = ax.get_position().y0
    fig.legend(handles, names, fontsize=8.8, loc="upper center", bbox_to_anchor=(0.5, below - 0.035), ncol=4,
               frameon=False, columnspacing=1.4, handlelength=2.2)

    labels = [NAMES[t] for t in TERMINATIONS]
    y = np.arange(len(TERMINATIONS))
    cx = fig.add_subplot(grid[1, 0])
    a = [d["error_route_a_mev"][t] for t in TERMINATIONS]
    b = [d["error_route_b_mev"][t] for t in TERMINATIONS]
    cx.barh(y - 0.19, a, height=0.36, color="#c9721c", label="model's own energy")
    cx.barh(y + 0.19, b, height=0.36, color="#2f6db5", label="CP2K static at the MACE geometry")
    cx.axvline(2.0, color="#2f6db5", linewidth=0.9, linestyle="--")
    cx.text(2.08, -0.95, "gate: 2 meV/Å$^2$", fontsize=7.8, color="#2f6db5", va="top")
    cx.set_ylim(len(TERMINATIONS) - 0.4, -1.1)
    cx.set_yticks(y, labels, fontsize=8.4)
    cx.set_xlabel(r"largest $|\gamma - \gamma_\mathrm{CP2K}|$ (meV/Å$^2$)", fontsize=9.5)
    cx.set_title("(c) MACE for the geometry, DFT for the energy", loc="left", fontsize=11, fontweight="bold")
    cx.legend(fontsize=8.4, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)

    ex = fig.add_subplot(grid[1, 1])
    m2 = [100 * d["recovered_m2"][t] for t in TERMINATIONS]
    m2r = [100 * d["recovered_m2r"][t] for t in TERMINATIONS]
    ex.barh(y - 0.19, m2, height=0.36, color="#c7cdd4", label="first model")
    ex.barh(y + 0.19, m2r, height=0.36, color="#2e8b57", label="after one active-learning round")
    ex.axvline(90, color="#3e4c59", linewidth=0.9, linestyle="--")
    ex.text(89.6, -0.95, "gate: 90 %", fontsize=7.8, color="#3e4c59", va="top", ha="right")
    ex.set_xlim(60, 101)
    ex.set_ylim(len(TERMINATIONS) - 0.4, -1.1)
    ex.set_yticks(y, labels, fontsize=8.4)
    ex.set_xlabel("share of the CP2K relaxation energy recovered (%)", fontsize=9.5)
    ex.set_title("(d) Why the active-learning round", loc="left", fontsize=11, fontweight="bold")
    ex.legend(fontsize=8.4, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)
    for suffix in ("png", "svg"):
        fig.savefig(f"{OUT}.{suffix}", dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote", OUT, "Reuter transitions", theirs)


def _strip_transitions(mu, curves):
    labels = list(curves)
    table = np.array([curves[t] for t in labels])
    stable = table.argmin(axis=0)
    out = []
    for i in range(1, len(mu)):
        if stable[i] != stable[i - 1]:
            a, b = labels[stable[i - 1]], labels[stable[i]]
            # exact crossing of the two lines
            (ia, sa), (ib, sb) = REUTER[a], REUTER[b]
            out.append(((ib - ia) / (sa - sb), a, b))
    return out


if __name__ == "__main__":
    if sys.argv[1] == "data":
        from aiida import load_profile

        import cluster

        load_profile("psteros_mace_ruo2")
        cluster.with_retry(data)
    else:
        plot()
