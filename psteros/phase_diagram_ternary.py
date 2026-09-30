"""Surface phase diagrams of ternary oxides A_xB_yO_z over (Delta mu_A, Delta mu_O).

Bulk equilibrium, ``x dmu_A + y dmu_B + z dmu_O = Delta H_f``, eliminates the
chemical potential of B, so every termination has a surface energy that is a
plane over the two remaining chemical potentials.  The bulk is stable inside a
convex polygon: no element precipitates (every ``dmu <= 0``) and no competing
phase forms.  Both that polygon and the region where each termination is the
most stable are computed exactly.

As for binary oxides, the diagram is written directly as a figure
(:meth:`TernarySurfacePhaseDiagram.plot`) or exported as a CSV table
(:meth:`TernarySurfacePhaseDiagram.to_csv`).
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from functools import reduce
from math import gcd
from pathlib import Path
from typing import Any, Iterable, Mapping

from psteros.phase_diagram import (
    _AXIS,
    _INK,
    _INK_SECONDARY,
    _MUTED,
    _OUTSIDE,
    _SERIES_COLORS,
    _SURFACE,
    SlabTermination,
    Units,
    _integer_composition,
    _unit_conversion,
)
from psteros.thermodynamics import surface_energy_elemental

Point = tuple[float, float]
_TOLERANCE = 1e-9


@dataclass(frozen=True)
class CompetingPhase:
    """A bulk phase that must not be favoured over the oxide, e.g. SrO for SrTiO3.

    ``energy_ev`` is the energy of the calculated cell whose element counts are
    ``composition``; elemental phases are already covered by the references.
    """

    label: str
    energy_ev: float
    composition: Mapping[str, int]

    def __post_init__(self) -> None:
        if not self.label or not str(self.label).strip():
            raise ValueError("competing phase label must not be empty")
        object.__setattr__(self, "composition", _integer_composition(self.composition, "composition"))


@dataclass(frozen=True)
class TernaryOxideReferences:
    """Bulk, elemental and competing-phase references of a ternary oxide A_xB_yO_z.

    ``element_energies_per_atom_ev`` gives the energy per atom of the two
    non-oxygen elements in their reference phases, where ``Delta mu = 0``.
    ``independent`` is the element on the horizontal axis of the diagram (the
    alphabetically first by default); the other one is eliminated through bulk
    equilibrium.
    """

    bulk_energy_ev: float
    bulk_composition: Mapping[str, int]
    oxygen_molecule_energy_ev: float
    element_energies_per_atom_ev: Mapping[str, float]
    competing_phases: tuple[CompetingPhase, ...] = ()
    independent: str | None = None
    eliminated: str = field(init=False)
    formula_unit: tuple[int, int, int] = field(init=False)
    stability_region: tuple[Point, ...] = field(init=False)
    stability_boundaries: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        composition = _integer_composition(self.bulk_composition, "bulk_composition")
        cations = sorted(element for element in composition if element != "O")
        if "O" not in composition or len(cations) != 2:
            raise ValueError(f"bulk_composition must describe a ternary oxide A_xB_yO_z, got {composition}")
        independent = self.independent or cations[0]
        if independent not in cations:
            raise ValueError(f"independent must be one of {cations}, got {independent!r}")
        eliminated = next(element for element in cations if element != independent)
        missing = sorted(set(cations).difference(self.element_energies_per_atom_ev))
        if missing:
            raise ValueError(f"element_energies_per_atom_ev lacks {missing}")
        divisor = reduce(gcd, (composition[independent], composition[eliminated], composition["O"]))
        for name, value in (
            ("bulk_composition", composition),
            ("independent", independent),
            ("eliminated", eliminated),
            ("formula_unit", tuple(composition[e] // divisor for e in (independent, eliminated, "O"))),
            ("competing_phases", tuple(self.competing_phases)),
            ("element_energies_per_atom_ev", {e: float(self.element_energies_per_atom_ev[e]) for e in cations}),
        ):
            object.__setattr__(self, name, value)
        for phase in self.competing_phases:
            foreign = sorted(set(phase.composition).difference(composition))
            if foreign:
                raise ValueError(f"competing phase {phase.label!r} contains elements {foreign} absent from the bulk")
        if self.formation_enthalpy_ev >= 0:
            raise ValueError(
                f"formation enthalpy {self.formation_enthalpy_ev:.4f} eV is not negative: "
                "the oxide is unstable against its elements"
            )
        region, boundaries = self._stability_region()
        object.__setattr__(self, "stability_region", region)
        object.__setattr__(self, "stability_boundaries", boundaries)

    @property
    def formula(self) -> str:
        counts = zip((self.independent, self.eliminated, "O"), self.formula_unit)
        return "".join(f"{element}{count if count > 1 else ''}" for element, count in counts)

    @property
    def bulk_energy_per_formula_unit_ev(self) -> float:
        return self.bulk_energy_ev * self.formula_unit[0] / self.bulk_composition[self.independent]

    def _reference_energy(self, element: str) -> float:
        if element == "O":
            return self.oxygen_molecule_energy_ev / 2.0
        return self.element_energies_per_atom_ev[element]

    def _formation_energy(self, energy_ev: float, composition: Mapping[str, int]) -> float:
        return energy_ev - sum(count * self._reference_energy(e) for e, count in composition.items())

    @property
    def formation_enthalpy_ev(self) -> float:
        """Delta H_f per A_xB_yO_z formula unit relative to the elements, at 0 K."""

        x, y, z = self.formula_unit
        composition = {self.independent: x, self.eliminated: y, "O": z}
        return self._formation_energy(self.bulk_energy_per_formula_unit_ev, composition)

    def delta_mu_eliminated_ev(self, delta_mu_independent_ev: float, delta_mu_oxygen_ev: float) -> float:
        """Delta mu of the eliminated element that keeps the bulk in equilibrium."""

        x, y, z = self.formula_unit
        return (self.formation_enthalpy_ev - x * delta_mu_independent_ev - z * delta_mu_oxygen_ev) / y

    def chemical_potentials_ev(self, delta_mu_independent_ev: float, delta_mu_oxygen_ev: float) -> dict[str, float]:
        """Absolute chemical potentials of all three elements at one point of the diagram."""

        return {
            self.independent: self._reference_energy(self.independent) + delta_mu_independent_ev,
            self.eliminated: self._reference_energy(self.eliminated)
            + self.delta_mu_eliminated_ev(delta_mu_independent_ev, delta_mu_oxygen_ev),
            "O": self._reference_energy("O") + delta_mu_oxygen_ev,
        }

    def in_stability_region(self, delta_mu_independent_ev: float, delta_mu_oxygen_ev: float) -> bool:
        return all(
            a * delta_mu_independent_ev + b * delta_mu_oxygen_ev <= c + _TOLERANCE
            for a, b, c, _ in self._constraints()
        )

    def _constraints(self) -> list[tuple[float, float, float, str]]:
        """Half-planes ``a dmu_A + b dmu_O <= c`` of the bulk stability region."""

        x, y, z = self.formula_unit
        enthalpy = self.formation_enthalpy_ev
        constraints = [
            (1.0, 0.0, 0.0, self.independent),
            (0.0, 1.0, 0.0, "O2"),
            # dmu_B <= 0 with dmu_B from bulk equilibrium
            (-x / y, -z / y, -enthalpy / y, self.eliminated),
        ]
        for phase in self.competing_phases:
            counts = phase.composition
            n_a, n_b, n_o = (counts.get(e, 0) for e in (self.independent, self.eliminated, "O"))
            constraints.append((
                n_a - n_b * x / y,
                n_o - n_b * z / y,
                self._formation_energy(phase.energy_ev, counts) - n_b * enthalpy / y,
                phase.label,
            ))
        return constraints

    def _stability_region(self) -> tuple[tuple[Point, ...], tuple[str, ...]]:
        constraints = self._constraints()
        x, _, z = self.formula_unit
        # Elements alone bound a triangle; competing phases cut it further.
        polygon = [(0.0, 0.0), (self.formation_enthalpy_ev / x, 0.0), (0.0, self.formation_enthalpy_ev / z)]
        for a, b, c, label in constraints[3:]:
            polygon = _clip(polygon, a, b, c)
            if _area(polygon) <= _TOLERANCE:
                raise ValueError(
                    f"{self.formula} has no stability region: it is unstable against {label!r} "
                    "(together with the phases listed before it)"
                )
        boundaries = []
        for start, end in zip(polygon, polygon[1:] + polygon[:1]):
            labels = [
                label for a, b, c, label in constraints
                if abs(a * start[0] + b * start[1] - c) <= 1e-7 and abs(a * end[0] + b * end[1] - c) <= 1e-7
            ]
            boundaries.append(labels[-1] if labels else "")
        return tuple(polygon), tuple(boundaries)


@dataclass(frozen=True)
class TernarySurfacePhaseDiagram:
    """Surface energies over (Delta mu_A, Delta mu_O) and the stable termination regions.

    ``planes`` maps each label to ``(gamma_0, d gamma / d dmu_A, d gamma / d dmu_O)``
    in eV/A^2, the exact plane ``gamma = gamma_0 + s_A dmu_A + s_O dmu_O``.
    ``regions`` maps each label to the polygon, inside the stability region,
    where it is the most stable termination (empty if it never is).
    """

    references: TernaryOxideReferences
    terminations: tuple[SlabTermination, ...]
    delta_mu_independent_ev: tuple[float, ...]
    delta_mu_oxygen_ev: tuple[float, ...]
    planes: Mapping[str, tuple[float, float, float]]
    regions: Mapping[str, tuple[Point, ...]]

    def gamma_ev_per_angstrom2(self, label: str, delta_mu_independent_ev: float, delta_mu_oxygen_ev: float) -> float:
        gamma0, slope_a, slope_o = self.planes[label]
        return gamma0 + slope_a * delta_mu_independent_ev + slope_o * delta_mu_oxygen_ev

    def stable_termination(self, delta_mu_independent_ev: float, delta_mu_oxygen_ev: float) -> str:
        return min(
            self.planes,
            key=lambda label: self.gamma_ev_per_angstrom2(label, delta_mu_independent_ev, delta_mu_oxygen_ev),
        )

    def to_csv(self, path: str | Path, *, units: Units = "J/m2") -> Path:
        """Write one CSV row per grid point and return the path.

        Columns: ``delta_mu_<A>_eV``, ``delta_mu_O_eV``, ``delta_mu_<B>_eV`` (from
        bulk equilibrium), one ``gamma_<label>_Jm2`` (or ``_eVA2``) column per
        termination, ``stable_termination`` and ``in_stability_region``.
        """

        factor, suffix = _unit_conversion(units)
        references = self.references
        path = Path(path)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [f"delta_mu_{references.independent}_eV", "delta_mu_O_eV", f"delta_mu_{references.eliminated}_eV"]
                + [f"gamma_{label}_{suffix}" for label in self.planes]
                + ["stable_termination", "in_stability_region"]
            )
            for delta_mu_o in self.delta_mu_oxygen_ev:
                for delta_mu_a in self.delta_mu_independent_ev:
                    writer.writerow(
                        [repr(delta_mu_a), repr(delta_mu_o), repr(references.delta_mu_eliminated_ev(delta_mu_a, delta_mu_o))]
                        + [repr(self.gamma_ev_per_angstrom2(label, delta_mu_a, delta_mu_o) * factor) for label in self.planes]
                        + [
                            self.stable_termination(delta_mu_a, delta_mu_o),
                            references.in_stability_region(delta_mu_a, delta_mu_o),
                        ]
                    )
        return path

    def figure(self, *, title: str | None = None) -> Any:
        """Return a matplotlib ``Figure`` of the stable-termination regions."""

        from matplotlib.figure import Figure
        from matplotlib.patches import Patch

        references = self.references
        region = references.stability_region
        xs, ys = [p[0] for p in region], [p[1] for p in region]
        pad_x, pad_y = 0.12 * (max(xs) - min(xs)), 0.12 * (max(ys) - min(ys))

        fig = Figure(figsize=(7.2, 6.2), facecolor=_SURFACE)
        ax = fig.subplots()
        ax.set_facecolor(_OUTSIDE)  # everything outside the stability region
        colors = {label: _SERIES_COLORS[i % len(_SERIES_COLORS)] for i, label in enumerate(self.planes)}
        total = _area(list(region))
        for label, polygon in self.regions.items():
            if not polygon:
                continue
            ax.fill(
                [p[0] for p in polygon], [p[1] for p in polygon], facecolor=colors[label],
                edgecolor=_SURFACE, linewidth=1.5, zorder=2,
            )
            # Label only regions wide enough to hold the text; thin ones rely on the legend.
            if _area(list(polygon)) > 0.15 * total:
                cx, cy = _centroid(polygon)
                ax.text(
                    cx, cy, _formula_text(label), ha="center", va="center", fontsize=9, color=_INK, zorder=4,
                    bbox={"boxstyle": "round,pad=0.25", "fc": _SURFACE, "ec": "none", "alpha": 0.85},
                )
        ax.plot(xs + xs[:1], ys + ys[:1], color=_INK_SECONDARY, lw=1, zorder=3)

        # Name the phase that forms beyond each edge of the stability region, set
        # parallel to the edge and offset along its outward normal on screen.
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)
        to_screen = ax.transData.transform
        cx, cy = to_screen(_centroid(region))
        for (start, end), label in zip(zip(region, region[1:] + region[:1]), references.stability_boundaries):
            if not label:
                continue
            (x0, y0), (x1, y1) = to_screen(start), to_screen(end)
            length = math.hypot(x1 - x0, y1 - y0) or 1.0
            nx, ny = (y1 - y0) / length, -(x1 - x0) / length
            if nx * ((x0 + x1) / 2 - cx) + ny * ((y0 + y1) / 2 - cy) < 0:
                nx, ny = -nx, -ny
            angle = math.degrees(math.atan2(y1 - y0, x1 - x0))
            if angle > 90:
                angle -= 180
            elif angle < -90:
                angle += 180
            ax.annotate(
                _formula_text(label), xy=((start[0] + end[0]) / 2, (start[1] + end[1]) / 2),
                xytext=(7 * nx, 7 * ny), textcoords="offset points", rotation=angle,
                rotation_mode="anchor", ha="center", va="center", fontsize=8, color=_MUTED, zorder=5,
            )

        present = [label for label, polygon in self.regions.items() if polygon]
        if len(present) > 1:
            ax.legend(
                handles=[Patch(color=colors[label], label=_formula_text(label)) for label in present],
                loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=min(len(present), 4),
                frameon=False, fontsize=9, labelcolor=_INK_SECONDARY, borderaxespad=0.2,
            )
        ax.set_xlabel(rf"$\Delta\mu_\mathrm{{{references.independent}}}$  (eV)", color=_INK_SECONDARY)
        ax.set_ylabel(r"$\Delta\mu_\mathrm{O}$  (eV)", color=_INK_SECONDARY)
        if title:
            ax.set_title(title, loc="left", fontsize=11.5, color=_INK, pad=30 if len(present) > 1 else 10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("bottom", "left"):
            ax.spines[side].set_color(_AXIS)
        ax.tick_params(colors=_MUTED, length=0)
        return fig

    def plot(self, path: str | Path, *, title: str | None = None, dpi: int = 200) -> Path:
        """Draw the stable-termination map and save it; the format follows the suffix."""

        path = Path(path)
        self.figure(title=title).savefig(path, dpi=dpi, bbox_inches="tight", facecolor=_SURFACE)
        return path


def ternary_surface_phase_diagram(
    terminations: Iterable[SlabTermination],
    references: TernaryOxideReferences,
    *,
    points: int = 101,
) -> TernarySurfacePhaseDiagram:
    """Evaluate every termination over the stability region of a ternary oxide.

    ``points`` is the number of grid values per axis used by the CSV export;
    the regions and planes are exact.
    """

    terminations = tuple(terminations)
    if not terminations:
        raise ValueError("at least one termination is required")
    labels = [termination.label for termination in terminations]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"termination labels must be unique: {duplicates}")
    allowed = {references.independent, references.eliminated, "O"}
    for termination in terminations:
        foreign = sorted(set(termination.composition).difference(allowed))
        if foreign:
            raise ValueError(
                f"termination {termination.label!r} contains {foreign}; expected only {sorted(allowed)}"
            )
    if points < 2:
        raise ValueError("points must be at least 2")

    def gamma(termination: SlabTermination, delta_mu_a: float, delta_mu_o: float) -> float:
        return surface_energy_elemental(
            slab_energy_ev=termination.slab_energy_ev,
            stoichiometry=termination.composition,
            chemical_potentials_ev=references.chemical_potentials_ev(delta_mu_a, delta_mu_o),
            surface_area_angstrom2=termination.surface_area_angstrom2,
            surfaces=termination.surfaces,
        )

    planes = {}
    for termination in terminations:
        gamma0 = gamma(termination, 0.0, 0.0)
        planes[termination.label] = (
            gamma0, gamma(termination, 1.0, 0.0) - gamma0, gamma(termination, 0.0, 1.0) - gamma0,
        )

    regions = {}
    for index, (label, (g_i, a_i, o_i)) in enumerate(planes.items()):
        polygon = list(references.stability_region)
        for other_index, (g_j, a_j, o_j) in enumerate(planes.values()):
            if other_index == index or not polygon:
                continue
            a, b, c = a_i - a_j, o_i - o_j, g_j - g_i  # gamma_i <= gamma_j
            if abs(a) <= _TOLERANCE and abs(b) <= _TOLERANCE:
                # Parallel planes: lower one wins; identical ones go to the first label.
                if c < -_TOLERANCE or (abs(c) <= _TOLERANCE and other_index < index):
                    polygon = []
                continue
            polygon = _clip(polygon, a, b, c)
            if _area(polygon) <= _TOLERANCE:
                polygon = []
        regions[label] = tuple(polygon)

    xs = [p[0] for p in references.stability_region]
    ys = [p[1] for p in references.stability_region]
    return TernarySurfacePhaseDiagram(
        references=references,
        terminations=terminations,
        delta_mu_independent_ev=_linspace(min(xs), max(xs), points),
        delta_mu_oxygen_ev=_linspace(min(ys), max(ys), points),
        planes=planes,
        regions=regions,
    )


def _formula_text(label: str) -> str:
    """Subscript the counts of formula-like words for matplotlib (TiO2 (anatase) -> TiO$_2$ (anatase))."""

    def subscript(match: re.Match) -> str:
        return re.sub(r"(\d+)", r"$_{\1}$", match.group(0))

    return re.sub(r"(?<![\w$])(?:[A-Z][a-z]?\d*)+(?![\w$])", subscript, label)


def _linspace(low: float, high: float, points: int) -> tuple[float, ...]:
    step = (high - low) / (points - 1)
    return tuple(low + index * step for index in range(points - 1)) + (high,)


def _clip(polygon: list[Point], a: float, b: float, c: float) -> list[Point]:
    """Keep the part of a convex polygon where ``a x + b y <= c`` (Sutherland-Hodgman)."""

    result: list[Point] = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        f_start = a * start[0] + b * start[1] - c
        f_end = a * end[0] + b * end[1] - c
        if f_start <= _TOLERANCE:
            result.append(start)
        if (f_start < -_TOLERANCE and f_end > _TOLERANCE) or (f_start > _TOLERANCE and f_end < -_TOLERANCE):
            t = f_start / (f_start - f_end)
            result.append((start[0] + t * (end[0] - start[0]), start[1] + t * (end[1] - start[1])))
    deduplicated: list[Point] = []
    for point in result:
        if not deduplicated or max(abs(point[0] - deduplicated[-1][0]), abs(point[1] - deduplicated[-1][1])) > 1e-12:
            deduplicated.append(point)
    if len(deduplicated) > 1 and max(
        abs(deduplicated[0][0] - deduplicated[-1][0]), abs(deduplicated[0][1] - deduplicated[-1][1])
    ) <= 1e-12:
        deduplicated.pop()
    return deduplicated


def _area(polygon: list[Point] | tuple[Point, ...]) -> float:
    if len(polygon) < 3:
        return 0.0
    return abs(sum(
        x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(polygon, list(polygon[1:]) + [polygon[0]])
    )) / 2.0


def _centroid(polygon: tuple[Point, ...] | list[Point]) -> Point:
    signed = sum(
        x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(polygon, list(polygon[1:]) + [polygon[0]])
    ) / 2.0
    if abs(signed) <= 1e-15:
        return sum(p[0] for p in polygon) / len(polygon), sum(p[1] for p in polygon) / len(polygon)
    cx = sum(
        (x0 + x1) * (x0 * y1 - x1 * y0) for (x0, y0), (x1, y1) in zip(polygon, list(polygon[1:]) + [polygon[0]])
    ) / (6.0 * signed)
    cy = sum(
        (y0 + y1) * (x0 * y1 - x1 * y0) for (x0, y0), (x1, y1) in zip(polygon, list(polygon[1:]) + [polygon[0]])
    ) / (6.0 * signed)
    return cx, cy
