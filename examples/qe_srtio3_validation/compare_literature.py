"""Compare the psteros SrTiO3(001) results with published calculations.

All surface grand potentials are put on one convention: gamma (J/m2) of the
1x1 SrO- and TiO2-terminated surfaces across the window where SrTiO3 is stable
against SrO and TiO2 (rutile), from the SrO-rich edge (0) to the TiO2-rich edge (1).
Literature lines are straight lines through the end points read from the
published figures, so they carry the reading uncertainty noted below.

    python compare_literature.py results/validation.json --out results
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

EV_A2_TO_J_M2 = 16.021766

# Values quoted in the text of each paper (exact), or read from a figure (approximate).
LITERATURE = {
    "Padilla & Vanderbilt 1998 (LDA)": {
        "reference": "J. Padilla and D. Vanderbilt, Surf. Sci. 418, 64 (1998), arXiv:cond-mat/9802207",
        "average_J_m2": 1.358,  # text: 1.26 eV per surface cell (1358 erg/cm2)
        "relaxation_eV_per_cell": 0.18,  # text: E_unrel 1.44 eV -> 1.26 eV
        "geometry_A": {"SrO": (0.22, -0.26, 0.10), "TiO2": (0.07, -0.27, 0.12)},  # Table 2: s, dd12, dd23
        # Fig. 1, F in eV per surface cell vs mu_TiO2; the axis spans about -3.2 to 0 eV, the
        # window implied by their LDA formation energy (their BaTiO3 paper reports E_f = 3.23 eV).
        # Read: SrO 0.65 -> 2.25, TiO2 1.9 -> 0.3 eV; a = 3.86 A (A = 14.90 A2). Reading error ~0.1 eV.
        "window_eV": 3.2,
        "lines_J_m2": {"SrO": (0.65, 2.25), "TiO2": (1.9, 0.3)},
        "line_units": "eV per cell",
        "area_A2": 3.86**2,
    },
    "Johnston et al. 2004 (LDA)": {
        "reference": "K. Johnston, M. R. Castell, A. T. Paxton, M. W. Finnis, Phys. Rev. B 70, 085415 (2004)",
        "average_J_m2": 1.46,  # text
        "formation_from_oxides_eV": -0.1094 * 13.605693,  # text: Delta G_f = -0.1094 Ry
        # Fig. 4, sigma in J/m2 vs mu_TiO2 (rutile): SrO 0.8 -> 1.62, TiO2 2.1 -> 1.3. Reading error ~0.05 J/m2.
        "window_eV": 0.1094 * 13.605693,
        "lines_J_m2": {"SrO": (0.8, 1.62), "TiO2": (2.1, 1.3)},
        "line_units": "J/m2",
    },
    "Eglitis & Vanderbilt 2008 (B3PW)": {
        "reference": "R. I. Eglitis and D. Vanderbilt, Phys. Rev. B 77, 195408 (2008), Table VII",
        "average_J_m2": 0.5 * (1.15 + 1.23) / 3.904**2 * EV_A2_TO_J_M2,  # eV per cell on a = 3.904 A
        "relaxation_eV_per_cell": 0.5 * (0.24 + 0.16),
        "geometry_A": {"SrO": tuple(x * 3.904 / 100 for x in (5.66, -6.58, 1.75)),
                       "TiO2": tuple(x * 3.904 / 100 for x in (2.12, -5.79, 3.55))},
    },
}
PBE_LATTICE = ("3.94 Å (VASP PBE)", "K. V. Sopiha et al., arXiv:1705.05250")


def ours(results: dict) -> dict:
    """Our gamma lines on the rutile-referenced window, from the psteros planes' end points."""
    area = results["area_A2"]
    dh = results["formation_enthalpy_eV"]
    window = -(dh["srtio3"] - dh["sro"] - dh["tio2_rutile"])
    e_surf = results["E_surf_eV_per_cell"]
    average = 0.5 * (e_surf["SrO"] + e_surf["TiO2"])
    # gamma_SrO - gamma_TiO2 is linear in mu_TiO2 with slope 1 eV/cell per eV; it vanishes at the crossing.
    low, high = results["strip_delta_mu_SrO_eV"]
    crossing_srO = results["termination_crossing_delta_mu_SrO_eV"]
    # Delta mu_SrO is high at the SrO-rich edge; with rutile the TiO2-rich edge sits at high - window.
    distance_from_tio2_rich = crossing_srO - (high - window)
    diff_tio2_rich = distance_from_tio2_rich  # gamma_SrO - gamma_TiO2 at the TiO2-rich edge, eV/cell
    to_j = EV_A2_TO_J_M2 / area
    lines = {
        "SrO": ((average - (window - diff_tio2_rich) / 2) * to_j, (average + diff_tio2_rich / 2) * to_j),
        "TiO2": ((average + (window - diff_tio2_rich) / 2) * to_j, (average - diff_tio2_rich / 2) * to_j),
    }
    return {"window_eV": window, "lines_J_m2": lines, "average_J_m2": average * to_j, "area_A2": area}


def crossing_fraction(lines) -> float:
    """Position (0 = SrO-rich, 1 = TiO2-rich) where the two straight lines cross."""
    (s0, s1), (t0, t1) = lines["SrO"], lines["TiO2"]
    return (t0 - s0) / ((s1 - s0) - (t1 - t0))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("validation_json")
    parser.add_argument("--out", default="results")
    args = parser.parse_args(argv)
    results = json.loads(Path(args.validation_json).read_text())
    mine = ours(results)

    rows = ["| Source | Window ΔG(SrO+TiO₂→SrTiO₃), eV | Average γ, J/m² | TiO₂ termination stable over | Crossing from TiO₂-rich edge |",
            "|---|---|---|---|---|"]
    series = {"This work (QE PBE, psteros)": mine}
    for name, data in LITERATURE.items():
        if "lines_J_m2" not in data:
            continue
        lines = data["lines_J_m2"]
        if data.get("line_units") == "eV per cell":
            factor = EV_A2_TO_J_M2 / data["area_A2"]
            lines = {k: (v[0] * factor, v[1] * factor) for k, v in lines.items()}
        series[name] = {"window_eV": data["window_eV"], "lines_J_m2": lines, "average_J_m2": data["average_J_m2"]}
    for name, data in series.items():
        f = crossing_fraction(data["lines_J_m2"])
        rows.append(
            f"| {name} | {-data['window_eV']:.2f} | {data['average_J_m2']:.2f} | {100 * (1 - f):.0f}% of the window "
            f"| {(1 - f) * data['window_eV']:.2f} eV |"
        )
    rows.append("| Experiment | −1.38 ± 0.10 | — | — | — |")

    geo = ["| Source | SrO: s, Δd₁₂, Δd₂₃ (Å) | TiO₂: s, Δd₁₂, Δd₂₃ (Å) | Relaxation energy, eV/cell |", "|---|---|---|---|"]
    a = results["lattice_constant_A"]
    mine_geo = {t: tuple(v * a / 100 for v in results["relaxation_percent_of_a"][t]) for t in ("SrO", "TiO2")}
    relax = -0.5 * sum(results["E_rel_eV_per_cell"].values())
    fmt = lambda t: ", ".join(f"{x:+.2f}" for x in t)  # noqa: E731
    geo.append(f"| This work (QE PBE) | {fmt(mine_geo['SrO'])} | {fmt(mine_geo['TiO2'])} | {relax:.2f} |")
    for name, data in LITERATURE.items():
        if "geometry_A" in data:
            geo.append(f"| {name} | {fmt(data['geometry_A']['SrO'])} | {fmt(data['geometry_A']['TiO2'])} | {data['relaxation_eV_per_cell']:.2f} |")

    text = "\n".join([
        "# SrTiO3(001): comparison with the literature", "",
        "## Phase diagram of the 1×1 terminations (window against SrO and rutile TiO₂)", "", *rows, "",
        "## Relaxation geometry and energy", "", *geo, "",
        f"PBE lattice constant: this work {a:.3f} Å; literature {PBE_LATTICE[0]} ({PBE_LATTICE[1]}).", "",
        "Literature lines are read from published figures (Padilla & Vanderbilt Fig. 1, ~0.1 eV; "
        "Johnston et al. Fig. 4, ~0.05 J/m²). References:",
        *[f"- {data['reference']}" for data in LITERATURE.values()],
    ]) + "\n"
    out = Path(args.out)
    (out / "literature_comparison.md").write_text(text)
    print(text)
    plot(series, out / "literature_comparison.png")


def plot(series, path):
    from matplotlib.figure import Figure

    colors = {"SrO": "#2a78d6", "TiO2": "#eb6834"}
    styles = ["-", (0, (5, 3)), (0, (1, 2))]
    fig = Figure(figsize=(7.4, 4.8), facecolor="#fcfcfb")
    ax = fig.subplots()
    ax.set_facecolor("#fcfcfb")
    for style, (name, data) in zip(styles, series.items()):
        for termination, (start, end) in data["lines_J_m2"].items():
            ax.plot([0, 1], [start, end], color=colors[termination], lw=2, linestyle=style)
        f = crossing_fraction(data["lines_J_m2"])
        s0, s1 = data["lines_J_m2"]["SrO"]
        ax.plot([f], [s0 + f * (s1 - s0)], marker="o", ms=6, color="#0b0b0b", zorder=5)
        ax.plot([], [], color="#52514e", lw=2, linestyle=style, label=name)
    ax.plot([], [], color=colors["SrO"], lw=6, label="1×1 SrO termination")
    ax.plot([], [], color=colors["TiO2"], lw=6, label="1×1 TiO$_2$ termination")
    ax.set_xticks([0, 0.5, 1], ["SrO-rich\n(SrTiO$_3$ + SrO)", "middle of the window", "TiO$_2$-rich\n(SrTiO$_3$ + rutile)"])
    ax.set_ylabel("Surface free energy γ (J/m$^2$)", color="#52514e")
    ax.set_title("SrTiO$_3$(001): this work (PBE) vs. literature (LDA); ● = crossing", loc="left", fontsize=11, color="#0b0b0b")
    ax.grid(axis="y", color="#e1e0d9", lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.tick_params(colors="#898781", length=0)
    ax.legend(frameon=False, fontsize=8.5, labelcolor="#52514e", loc="upper center", ncol=2)
    ax.set_ylim(0, 3.2)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
