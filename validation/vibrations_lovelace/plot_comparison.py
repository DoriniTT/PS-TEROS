"""Phase diagram without and with the vibrational contributions, side by side.

Reads the CSV files that ``analyse.py`` writes (no AiiDA profile needed):

    python plot_comparison.py [--temperature 900] [--out results/phase_diagram_comparison.png]

Left: static total energies (the diagram as it was before the feature). Right: E + F_vib(T). Both panels share the axes
and the colour of each termination; the stability window of the bulk oxide ends at the O-poor limit of each case, and the
vertical line marks where the stable termination changes.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
LABELS = {"slab_o": "o  (Sn$_6$O$_{12}$, stoichiometric)", "slab_sn2o": "sn2o  (Sn$_6$O$_8$, reduced)"}
COLOR = {"slab_o": "#2a78d6", "slab_sn2o": "#eb6834"}  # categorical slots 1 and 2
INK, INK_SECONDARY, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"


def read(path: Path) -> dict:
    rows = list(csv.DictReader(path.open()))
    data = {"x": [float(r["delta_mu_O_eV"]) for r in rows], "stable": [r["stable_termination"] for r in rows]}
    for label in LABELS:
        data[label] = [float(r[f"gamma_{label}_Jm2"]) for r in rows]
    data["poor_limit"] = data["x"][0]
    switch = next(i for i in range(1, len(rows)) if data["stable"][i] != data["stable"][i - 1])
    # linear interpolation of the crossing between the two grid points around the switch
    i, j = switch - 1, switch
    d0 = data["slab_o"][i] - data["slab_sn2o"][i]
    d1 = data["slab_o"][j] - data["slab_sn2o"][j]
    data["transition"] = data["x"][i] + (data["x"][j] - data["x"][i]) * d0 / (d0 - d1)
    return data


def panel(axis, data: dict, title: str) -> None:
    for label in LABELS:
        axis.plot(data["x"], data[label], color=COLOR[label], linewidth=2.0, label=LABELS[label])
    axis.axvline(data["transition"], color=MUTED, linewidth=1.0, linestyle=(0, (4, 3)))
    axis.annotate(f"{data['transition']:.2f} eV", (data["transition"], 3.95), xytext=(4, 0), textcoords="offset points",
                  fontsize=8, color=INK_SECONDARY, va="top")
    axis.annotate("O-poor limit", (data["poor_limit"], 3.95), xytext=(3, 0), textcoords="offset points",
                  fontsize=8, color=MUTED, va="top")
    # thin strip under the curves: which termination is stable
    for label in LABELS:
        low = data["poor_limit"] if label == "slab_sn2o" else data["transition"]
        high = data["transition"] if label == "slab_sn2o" else 0.0
        axis.fill_between([low, high], -0.45, -0.3, color=COLOR[label], linewidth=0)
    axis.set_title(title, color=INK, fontsize=10, loc="left")
    axis.set_xlabel(r"$\Delta\mu_\mathrm{O} = \mu_\mathrm{O} - \frac{1}{2}E(\mathrm{O_2})$  (eV)")
    axis.set_xlim(-2.6, 0.0)
    axis.set_ylim(-0.45, 4.0)
    axis.grid(color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)


def main(argv=None) -> Path:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--temperature", type=int, default=900, choices=(0, 300, 600, 900))
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    out = args.out or HERE / "results" / f"phase_diagram_comparison_T{args.temperature}K.png"
    static = read(HERE / "results" / "phase_diagram_total_energy.csv")
    vib = read(HERE / "results" / f"phase_diagram_T{args.temperature}K.csv")

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK_SECONDARY,
                         "xtick.color": INK_SECONDARY, "ytick.color": INK_SECONDARY})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True, constrained_layout=True)
    panel(axes[0], static, "Without vibrations: static total energies")
    zpe = " (zero-point energy only)" if args.temperature == 0 else ""
    panel(axes[1], vib, f"With vibrations: E + F$_{{vib}}$ at {args.temperature} K{zpe}")
    axes[0].set_ylabel(r"Surface energy $\gamma$  (J m$^{-2}$)")
    axes[0].legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(0.0, 0.62))
    gamma_o_change = vib["slab_o"][-1] - static["slab_o"][-1]
    shift = vib["transition"] - static["transition"]
    fig.suptitle(
        f"SnO$_2$(110): $\\gamma_o$ at $\\Delta\\mu_O=0$ changes by {gamma_o_change:+.2f} J m$^{{-2}}$, "
        f"the transition moves by {shift:+.2f} eV",
        fontsize=9.5, color=INK_SECONDARY, y=1.04,
    )
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"static: transition {static['transition']:.3f} eV, gamma_o(0) {static['slab_o'][-1]:.3f}, gamma_sn2o(0) {static['slab_sn2o'][-1]:.3f} J/m2")
    print(f"{args.temperature} K:  transition {vib['transition']:.3f} eV, gamma_o(0) {vib['slab_o'][-1]:.3f}, gamma_sn2o(0) {vib['slab_sn2o'][-1]:.3f} J/m2")
    print("figure:", out)
    return out


if __name__ == "__main__":
    main()
