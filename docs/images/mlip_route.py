"""Draw the flowchart of the MLIP route of psteros: ``python docs/images/mlip_route.py`` writes mlip_route.png/.svg."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

HERE = Path(__file__).resolve().parent
INK = "#1f2933"
MUTED = "#52606d"
STYLE = {
    "input": ("#ffffff", "#3e4c59"),
    "dft": ("#e3edf8", "#2f6db5"),
    "gpu": ("#fdebd7", "#c9721c"),
    "check": ("#f0f0f0", "#6b7280"),
    "result": ("#e1f3e8", "#2e8b57"),
    "stop": ("#fbe4e4", "#b83232"),
}


def box(ax, x, y, w, h, kind, title, body="", number=None, title_size=11.5):
    face, edge = STYLE[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.16",
                                facecolor=face, edgecolor=edge, linewidth=1.6, zorder=2))
    left = x + 0.22
    if number is not None:
        ax.add_patch(plt.Circle((x + 0.42, y + h / 2), 0.25, facecolor=edge, edgecolor="none", zorder=3))
        ax.text(x + 0.42, y + h / 2, str(number), ha="center", va="center", color="white",
                fontsize=11, fontweight="bold", zorder=4)
        left = x + 0.88
    if body:
        ax.text(left, y + h - 0.2, title, ha="left", va="top", fontsize=title_size, fontweight="bold",
                color=INK, zorder=4)
        ax.text(left, y + h - 0.56, body, ha="left", va="top", fontsize=9.4, color=MUTED, zorder=4,
                linespacing=1.35)
    else:
        ax.text(left, y + h / 2, title, ha="left", va="center", fontsize=title_size, fontweight="bold",
                color=INK, zorder=4)


def tag(ax, right, top, text, kind):
    """A small label on the top border of a box, near its right end."""

    face, edge = STYLE[kind]
    ax.text(right - 0.25, top, text, ha="right", va="center", fontsize=8.2, fontweight="bold", color=edge, zorder=5,
            bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.12", facecolor="white", edgecolor=edge, linewidth=1.2))


def arrow(ax, start, end, color=INK, style="-|>", connection="arc3", width=1.6, dashed=False):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=14, color=color, linewidth=width,
                                 connectionstyle=connection, linestyle=(0, (4, 3)) if dashed else "-", zorder=1,
                                 shrinkA=0, shrinkB=0))


def main() -> None:
    width, height = 11.0, 13.0
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")
    x, w = 2.25, 6.0  # main column
    cx, right = x + w / 2, x + w
    side = right + 0.45  # side notes

    ax.text(0.3, 12.75, "The MLIP route of PS-TEROS", fontsize=18, fontweight="bold", color=INK, va="top")
    ax.text(0.3, 12.25, "MACE finds the relaxed geometry, DFT gives its energy  ·  for metallic systems",
            fontsize=10.8, color=MUTED, va="top")
    legend_x = 0.3
    for kind, label in (("dft", "DFT, on CPUs"), ("gpu", "MACE, on a GPU"), ("check", "check"), ("result", "result")):
        face, edge = STYLE[kind]
        ax.add_patch(FancyBboxPatch((legend_x, 11.6), 0.32, 0.22, boxstyle="round,pad=0,rounding_size=0.05",
                                    facecolor=face, edgecolor=edge, linewidth=1.3))
        ax.text(legend_x + 0.45, 11.71, label, fontsize=9.4, color=MUTED, va="center")
        legend_x += 0.45 + 0.105 * len(label) + 0.55

    # 0 inputs
    box(ax, x, 10.15, w, 1.15, "input", "Structures",
        "bulk and references (metal, O$_2$ in a box) relaxed with DFT\nslabs as built, not relaxed")
    arrow(ax, (cx, 10.15), (cx, 9.8))
    # 1 first statics
    box(ax, x, 8.8, w, 1.0, "dft", "One DFT static of every structure",
        "energies and forces: the training labels", number=1)
    tag(ax, right, 9.8, "DFT", "dft")
    arrow(ax, (cx, 8.8), (cx, 8.47))
    # 2 metal check (diamond)
    dy, dh, dw = 7.83, 0.66, 2.6
    ax.add_patch(Polygon([(cx, dy + dh), (cx + dw, dy), (cx, dy - dh), (cx - dw, dy)], closed=True,
                         facecolor=STYLE["check"][0], edgecolor=STYLE["check"][1], linewidth=1.6, zorder=2))
    ax.add_patch(plt.Circle((cx - dw - 0.4, dy), 0.25, facecolor=STYLE["check"][1], edgecolor="none", zorder=3))
    ax.text(cx - dw - 0.4, dy, "2", ha="center", va="center", color="white", fontsize=11, fontweight="bold", zorder=4)
    ax.text(cx, dy + 0.12, "Is the bulk a metal?", ha="center", va="center", fontsize=11, fontweight="bold",
            color=INK, zorder=4)
    ax.text(cx, dy - 0.2, "partly filled states at the Fermi level", ha="center", va="center", fontsize=8.6,
            color=MUTED, zorder=4)
    box(ax, side, dy - 0.5, 2.0, 1.0, "stop", "Stop", "validated for\nmetals only", title_size=10.5)
    arrow(ax, (cx + dw, dy), (side, dy), color=STYLE["stop"][1])
    ax.text(cx + dw + 0.1, dy + 0.14, "no", fontsize=9, color=STYLE["stop"][1], fontweight="bold")
    arrow(ax, (cx, dy - dh), (cx, 6.85))
    ax.text(cx + 0.12, dy - dh - 0.17, "yes", fontsize=9, color=STYLE["result"][1], fontweight="bold")
    # 3 fine-tune
    box(ax, x, 5.5, w, 1.35, "gpu", "Replay fine-tune of MACE-MH-1",
        "a new head learns your DFT labels; the replay head keeps\n"
        "training on ~150 small structures labelled by MH-1,\nso the model does not forget general chemistry", number=3)
    tag(ax, right, 6.85, "GPU", "gpu")
    arrow(ax, (cx, 5.5), (cx, 5.15))
    # 4 relax
    box(ax, x, 4.15, w, 1.0, "gpu", "MACE relaxes every structure",
        "cell fixed, slab centre held in place  ·  minutes", number=4)
    tag(ax, right, 5.15, "GPU", "gpu")
    arrow(ax, (cx, 4.15), (cx, 3.8))
    # 5 active learning
    box(ax, x, 2.8, w, 1.0, "dft", "Active-learning round (once)",
        "one DFT static at each MACE geometry, added to the labels", number=5)
    tag(ax, right, 3.8, "DFT", "dft")
    # loop back to the fine-tune, in the left margin
    loop = STYLE["gpu"][1]
    lx = x - 0.55
    ax.plot([x, lx, lx], [3.3, 3.3, 6.17], color=loop, linewidth=1.9, solid_capstyle="round", zorder=1)
    arrow(ax, (lx, 6.17), (x, 6.17), color=loop, width=1.9)
    ax.text(lx - 0.15, 4.72, "fine-tune again\nfrom MH-1 on\nall labels, and\nrelax again from\nthe start", fontsize=8.6,
            color=loop, ha="right", va="center", fontweight="bold", linespacing=1.3)
    arrow(ax, (cx, 2.8), (cx, 2.47), dashed=True)
    ax.text(cx - 0.14, 2.635, "after the second relaxation", fontsize=8.6, color=MUTED, va="center", ha="right")
    # 6 final statics
    box(ax, x, 1.47, w, 1.0, "dft", "Final DFT static at each relaxed geometry",
        "these energies, and no model energy, enter the diagram", number=6)
    tag(ax, right, 2.47, "DFT", "dft")
    arrow(ax, (cx, 1.47), (cx, 1.12))
    # 7 result
    box(ax, x, 0.12, w, 1.0, "result", "Surface phase diagram  γ(Δμ$_\\mathrm{O}$)",
        "slabs and references from one energy source", number=7)

    # side notes
    def note(y, head, text, color):
        ax.text(side, y, head, fontsize=10.2, fontweight="bold", color=color, va="top")
        ax.text(side, y - 0.32, text, fontsize=9, color=MUTED, va="top", linespacing=1.38)

    note(6.8, "Why the extra round", "the first model has only seen\nunrelaxed slabs, so it stops\nshort of the minimum on\n"
         "some terminations", loop)
    note(4.75, "Validated on RuO$_2$(110)", "γ within 0.5 meV/Å$^2$\nphase boundaries within 0.02 eV\n"
         "96–99 % of the DFT relaxation", STYLE["result"][1])
    note(2.45, "DFT cost", "3 statics per termination:\n34 % of relaxing with DFT\n(RuO$_2$; more saving for\n"
         "longer relaxations)", STYLE["dft"][1])
    note(0.95, "Do not", "use the model's own energies:\n6–9 meV/Å$^2$ off", STYLE["stop"][1])

    for suffix in ("png", "svg"):
        fig.savefig(HERE / f"mlip_route.{suffix}", dpi=170, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
