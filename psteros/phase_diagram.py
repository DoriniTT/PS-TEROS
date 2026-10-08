"""Surface phase diagrams of binary compounds as a function of Delta mu.

For an oxide the axis is Delta mu_O (:class:`BinaryOxideReferences`); any
other binary compound A_xB_y uses :class:`BinaryReferences` with the axis on
Delta mu_B (for example Delta mu_As for GaAs).

Energies in, a :class:`SurfacePhaseDiagram` out.  The diagram can be written
directly as a figure (:meth:`SurfacePhaseDiagram.plot`) or exported as a CSV
table (:meth:`SurfacePhaseDiagram.to_csv`) for plotting elsewhere.  Nothing
here needs an AiiDA profile; matplotlib is imported only when a figure is
requested.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from math import gcd, isclose
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping

from psteros.thermodynamics import (
    EV_PER_ANGSTROM2_TO_J_PER_M2,
    SurfaceEnergyPoint,
    stable_termination,
    surface_energy_binary_equilibrium,
)

Units = Literal["J/m2", "eV/A2"]

# Categorical slots in fixed order (validated for colour-vision deficiency on
# a light surface).  Terminations beyond eight reuse the hues with a dashed
# line, so identity never rests on colour alone.
_SERIES_COLORS = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
)
_INK, _INK_SECONDARY, _MUTED = "#0b0b0b", "#52514e", "#898781"
_GRID, _AXIS, _SURFACE, _OUTSIDE = "#e1e0d9", "#c3c2b7", "#fcfcfb", "#f0efec"


def _integer_composition(composition: Any, name: str) -> dict[str, int]:
    """Return ``{element: count}`` from a mapping or a pymatgen ``Composition``."""

    if hasattr(composition, "get_el_amt_dict"):
        composition = composition.get_el_amt_dict()
    result: dict[str, int] = {}
    for element, amount in dict(composition).items():
        if not isclose(float(amount), round(float(amount)), abs_tol=1e-8) or amount < 0:
            raise ValueError(f"{name} must contain non-negative integer counts, got {element}={amount!r}")
        if round(float(amount)):
            result[str(element)] = int(round(float(amount)))
    if not result:
        raise ValueError(f"{name} must not be empty")
    return result


@dataclass(frozen=True)
class SlabTermination:
    """One slab termination entering a surface phase diagram.

    ``surface_area_angstrom2`` is the area of one exposed face and
    ``surfaces`` the number of equivalent faces the slab represents (2 for a
    symmetric slab).

    A polar slab with a pseudo-hydrogen passivated bottom exposes one face
    (``surfaces=1``); ``pseudo_hydrogen`` counts its pseudo-hydrogen by the
    element they passivate, ``composition`` excludes them, and
    ``bottom_fingerprint`` and ``face`` identify its bottom and its face. Use
    :meth:`from_polar` to fill these from :func:`psteros.find_polar_terminations`.
    """

    label: str
    slab_energy_ev: float
    composition: Mapping[str, int]
    surface_area_angstrom2: float
    surfaces: int = 2
    pseudo_hydrogen: Mapping[str, int] = field(default_factory=dict)
    bottom_fingerprint: str | None = None
    face: str | None = None

    def __post_init__(self) -> None:
        if not self.label or not str(self.label).strip():
            raise ValueError("termination label must not be empty")
        if self.surface_area_angstrom2 <= 0:
            raise ValueError("surface_area_angstrom2 must be positive")
        if self.surfaces <= 0:
            raise ValueError("surfaces must be positive")
        object.__setattr__(self, "composition", _integer_composition(self.composition, "composition"))
        hydrogen = {str(k): int(v) for k, v in dict(self.pseudo_hydrogen).items() if int(v)}
        if any(v < 0 for v in hydrogen.values()):
            raise ValueError("pseudo_hydrogen counts must be non-negative")
        object.__setattr__(self, "pseudo_hydrogen", hydrogen)
        if hydrogen and self.surfaces != 1:
            raise ValueError("a slab with a passivated bottom exposes one face: use surfaces=1")

    @property
    def is_passivated(self) -> bool:
        return bool(self.pseudo_hydrogen)

    @classmethod
    def from_polar(
        cls, termination: Any, slab_energy_ev: float, *, label: str | None = None
    ) -> "SlabTermination":
        """From a :class:`psteros.PolarTermination` and the energy of its relaxed slab.

        The counts, area, bottom fingerprint and face come from the built
        slab; a relaxation at fixed cell does not change them.
        """

        from psteros.core.terminations import hkl_label

        face = f"{termination.bulk_formula}{hkl_label(termination.miller_index)}"
        return cls(
            label=label or termination.label,
            slab_energy_ev=slab_energy_ev,
            composition=termination.composition,
            surface_area_angstrom2=termination.area,
            surfaces=1,
            pseudo_hydrogen=termination.pseudo_hydrogen_counts,
            bottom_fingerprint=termination.bottom_fingerprint,
            face=face,
        )

    @classmethod
    def from_structure(
        cls, label: str, slab_energy_ev: float, structure: Any, *, surfaces: int = 2
    ) -> "SlabTermination":
        """Read composition and area from a pymatgen structure or AiiDA ``StructureData``.

        The exposed face is taken as the plane of the first two lattice
        vectors, the convention of pymatgen and psteros slabs.
        """

        if hasattr(structure, "get_pymatgen_structure"):
            structure = structure.get_pymatgen_structure()
        import numpy as np

        a, b = structure.lattice.matrix[0], structure.lattice.matrix[1]
        return cls(
            label=label,
            slab_energy_ev=slab_energy_ev,
            composition=structure.composition,
            surface_area_angstrom2=float(np.linalg.norm(np.cross(a, b))),
            surfaces=surfaces,
        )


@dataclass(frozen=True)
class BinaryOxideReferences:
    """Bulk oxide and gas references that fix the oxygen chemical potential.

    ``bulk_energy_ev`` is the energy of the calculated bulk cell whose
    composition is ``bulk_composition`` (for example ``{"Sn": 2, "O": 4}``).
    ``oxygen_molecule_energy_ev`` is E(O2) from a triplet calculation.  With
    the optional ``metal_energy_per_atom_ev`` the O-poor limit
    ``Delta mu_O = Delta H_f / y`` is known; without it, callers must choose
    the Delta mu_O range themselves.
    """

    bulk_energy_ev: float
    bulk_composition: Mapping[str, int]
    oxygen_molecule_energy_ev: float
    metal_energy_per_atom_ev: float | None = None
    metal: str = field(init=False)
    formula_unit: tuple[int, int] = field(init=False)

    def __post_init__(self) -> None:
        composition = _integer_composition(self.bulk_composition, "bulk_composition")
        metals = sorted(element for element in composition if element != "O")
        if "O" not in composition or len(metals) != 1:
            raise ValueError(
                "bulk_composition must describe a binary oxide M_xO_y, "
                f"got {composition}"
            )
        divisor = gcd(composition[metals[0]], composition["O"])
        object.__setattr__(self, "bulk_composition", composition)
        object.__setattr__(self, "metal", metals[0])
        object.__setattr__(
            self, "formula_unit", (composition[metals[0]] // divisor, composition["O"] // divisor)
        )
        enthalpy = self.formation_enthalpy_ev
        if enthalpy is not None and enthalpy >= 0:
            raise ValueError(
                f"formation enthalpy {enthalpy:.4f} eV is not negative: the oxide is "
                "unstable against metal + O2 and no Delta mu_O window exists"
            )

    @property
    def formula(self) -> str:
        x, y = self.formula_unit
        return f"{self.metal}{x if x > 1 else ''}O{y if y > 1 else ''}"

    @property
    def bulk_energy_per_formula_unit_ev(self) -> float:
        return self.bulk_energy_ev * self.formula_unit[0] / self.bulk_composition[self.metal]

    @property
    def formation_enthalpy_ev(self) -> float | None:
        """Delta H_f per M_xO_y formula unit at 0 K, or ``None`` without a metal reference."""

        if self.metal_energy_per_atom_ev is None:
            return None
        x, y = self.formula_unit
        return (
            self.bulk_energy_per_formula_unit_ev
            - x * self.metal_energy_per_atom_ev
            - y * self.oxygen_molecule_energy_ev / 2.0
        )

    @property
    def oxygen_poor_limit_ev(self) -> float | None:
        """Delta mu_O below which the oxide decomposes into the metal."""

        enthalpy = self.formation_enthalpy_ev
        return None if enthalpy is None else enthalpy / self.formula_unit[1]

    oxygen_rich_limit_ev = 0.0

    # Element-neutral interface shared with BinaryReferences.
    variable = "O"
    rich_limit_ev = 0.0
    axis_label = (
        r"Oxygen chemical potential  $\Delta\mu_\mathrm{O} = \mu_\mathrm{O} - \frac{1}{2}E(\mathrm{O_2})$  (eV)"
    )

    @property
    def other(self) -> str:
        return self.metal

    @property
    def variable_reference_energy_ev(self) -> float:
        return self.oxygen_molecule_energy_ev / 2.0

    @property
    def poor_limit_ev(self) -> float | None:
        return self.oxygen_poor_limit_ev

    def chemical_potentials_ev(self, delta_mu_ev: float) -> dict[str, float]:
        """mu_M and mu_O at ``Delta mu_O``, with x mu_M + y mu_O = E_bulk."""

        x, y = self.formula_unit
        mu_oxygen = self.variable_reference_energy_ev + delta_mu_ev
        return {"O": mu_oxygen, self.metal: (self.bulk_energy_per_formula_unit_ev - y * mu_oxygen) / x}

    @property
    def poor_limit_label(self) -> str:
        return f" O-poor limit\n ({self.metal} metal)"

    @property
    def rich_limit_label(self) -> str:
        return "O-rich limit \n(O$_2$ gas) "


@dataclass(frozen=True)
class BinaryReferences:
    """Bulk and elemental references of any binary compound A_xB_y.

    ``reference_energies_per_atom_ev`` gives the reference energy per atom
    of each element, which fixes ``mu = reference + Delta mu``: the energy
    per atom of the elemental solid (Ga, As, Zn, ...) or half the energy of
    the molecule (O2, N2).  The axis element ``variable`` (B) must have a
    reference; with a reference for the other element too, the B-poor limit
    ``Delta mu_B = Delta H_f / y`` is known.  By default ``variable`` is the
    more electronegative element (As in GaAs, N in GaN, O in ZnO).

    ``reservoir_labels`` names the references in the figure, for example
    ``{"As": "As bulk", "Ga": "Ga bulk"}``.

    For oxides with an O2 reference, :class:`BinaryOxideReferences` gives the
    same diagram with the oxide-specific labels.
    """

    bulk_energy_ev: float
    bulk_composition: Mapping[str, int]
    reference_energies_per_atom_ev: Mapping[str, float]
    variable: str | None = None
    reservoir_labels: Mapping[str, str] = field(default_factory=dict)
    other: str = field(init=False)
    formula_unit: tuple[int, int] = field(init=False)

    def __post_init__(self) -> None:
        composition = _integer_composition(self.bulk_composition, "bulk_composition")
        if len(composition) != 2:
            raise ValueError(f"bulk_composition must describe a binary compound A_xB_y, got {composition}")
        references = {str(element): float(value) for element, value in dict(self.reference_energies_per_atom_ev).items()}
        foreign = sorted(set(references).difference(composition))
        if foreign:
            raise ValueError(f"reference energies given for elements not in the compound: {foreign}")
        variable = self.variable
        if variable is None:
            from pymatgen.core import Element

            variable = max(composition, key=lambda element: Element(element).X)
        if variable not in composition:
            raise ValueError(f"variable element {variable!r} is not in {sorted(composition)}")
        if variable not in references:
            raise ValueError(
                f"reference_energies_per_atom_ev must include the axis element {variable!r}; "
                "it defines Delta mu = 0"
            )
        other = next(element for element in composition if element != variable)
        divisor = gcd(composition[other], composition[variable])
        object.__setattr__(self, "bulk_composition", composition)
        object.__setattr__(self, "reference_energies_per_atom_ev", references)
        object.__setattr__(self, "reservoir_labels", dict(self.reservoir_labels))
        object.__setattr__(self, "variable", variable)
        object.__setattr__(self, "other", other)
        object.__setattr__(
            self, "formula_unit", (composition[other] // divisor, composition[variable] // divisor)
        )
        enthalpy = self.formation_enthalpy_ev
        if enthalpy is not None and enthalpy >= 0:
            raise ValueError(
                f"formation enthalpy {enthalpy:.4f} eV is not negative: {self.formula} is "
                f"unstable against its elements and no Delta mu_{variable} window exists"
            )

    @property
    def formula(self) -> str:
        x, y = self.formula_unit
        return f"{self.other}{x if x > 1 else ''}{self.variable}{y if y > 1 else ''}"

    @property
    def bulk_energy_per_formula_unit_ev(self) -> float:
        return self.bulk_energy_ev * self.formula_unit[0] / self.bulk_composition[self.other]

    @property
    def variable_reference_energy_ev(self) -> float:
        return self.reference_energies_per_atom_ev[self.variable]

    @property
    def formation_enthalpy_ev(self) -> float | None:
        """Delta H_f per A_xB_y formula unit, or ``None`` without a reference for A."""

        if self.other not in self.reference_energies_per_atom_ev:
            return None
        x, y = self.formula_unit
        return (
            self.bulk_energy_per_formula_unit_ev
            - x * self.reference_energies_per_atom_ev[self.other]
            - y * self.variable_reference_energy_ev
        )

    @property
    def poor_limit_ev(self) -> float | None:
        """Delta mu_B below which the compound decomposes into A."""

        enthalpy = self.formation_enthalpy_ev
        return None if enthalpy is None else enthalpy / self.formula_unit[1]

    rich_limit_ev = 0.0

    def chemical_potentials_ev(self, delta_mu_ev: float) -> dict[str, float]:
        """mu_A and mu_B at ``Delta mu_B``, with x mu_A + y mu_B = E_bulk."""

        x, y = self.formula_unit
        mu_variable = self.variable_reference_energy_ev + delta_mu_ev
        return {
            self.variable: mu_variable,
            self.other: (self.bulk_energy_per_formula_unit_ev - y * mu_variable) / x,
        }

    @property
    def axis_label(self) -> str:
        b = self.variable
        return (
            rf"Chemical potential  $\Delta\mu_\mathrm{{{b}}} = \mu_\mathrm{{{b}}}"
            rf" - \mu_\mathrm{{{b}}}^\mathrm{{ref}}$  (eV)"
        )

    @property
    def poor_limit_label(self) -> str:
        reservoir = self.reservoir_labels.get(self.other, f"{self.other} reference")
        return f" {self.variable}-poor limit\n ({reservoir})"

    @property
    def rich_limit_label(self) -> str:
        reservoir = self.reservoir_labels.get(self.variable, f"{self.variable} reference")
        return f"{self.variable}-rich limit \n({reservoir}) "


@dataclass(frozen=True)
class SurfacePhaseDiagram:
    """gamma(Delta mu) of every termination and the stable one at each point.

    Delta mu is that of ``references.variable`` (O for an oxide).
    ``transitions`` holds the exact Delta mu values, inside the sampled
    range, where the lowest-energy termination changes, as
    ``(delta_mu_ev, stable_below, stable_above)``.
    """

    references: BinaryOxideReferences | BinaryReferences
    terminations: tuple[SlabTermination, ...]
    curves: Mapping[str, tuple[SurfaceEnergyPoint, ...]]
    stable: tuple[str, ...]
    transitions: tuple[tuple[float, str, str], ...]

    @property
    def delta_mu_ev(self) -> tuple[float, ...]:
        return tuple(point.delta_mu_ev for point in next(iter(self.curves.values())))

    @property
    def delta_mu_oxygen_ev(self) -> tuple[float, ...]:
        """Historical name of :attr:`delta_mu_ev`."""

        return self.delta_mu_ev

    def in_stability_window(self, delta_mu_ev: float) -> bool:
        """Whether the bulk compound is stable at this Delta mu."""

        lower = self.references.poor_limit_ev
        upper = self.references.rich_limit_ev
        return delta_mu_ev <= upper + 1e-12 and (lower is None or delta_mu_ev >= lower - 1e-12)

    def to_csv(self, path: str | Path, *, units: Units = "J/m2") -> Path:
        """Write the diagram as one CSV row per Delta mu_O point and return the path.

        Columns: ``delta_mu_<B>_eV`` (``delta_mu_O_eV`` for an oxide), one
        ``gamma_<label>_Jm2`` (or ``_eVA2``) column per termination,
        ``stable_termination`` and ``in_stability_window`` (False where the
        bulk compound would decompose or where Delta mu > 0).
        """

        factor, suffix = _unit_conversion(units)
        path = Path(path)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [f"delta_mu_{self.references.variable}_eV"]
                + [f"gamma_{label}_{suffix}" for label in self.curves]
                + ["stable_termination", "in_stability_window"]
            )
            for index, delta_mu in enumerate(self.delta_mu_ev):
                writer.writerow(
                    [repr(delta_mu)]
                    + [repr(points[index].gamma_ev_per_angstrom2 * factor) for points in self.curves.values()]
                    + [self.stable[index], self.in_stability_window(delta_mu)]
                )
        return path

    def figure(self, *, units: Units = "J/m2", title: str | None = None) -> Any:
        """Return a matplotlib ``Figure`` of gamma(Delta mu) with a stability strip."""

        from matplotlib.figure import Figure

        factor, _ = _unit_conversion(units)
        grid = self.delta_mu_ev
        low, high = grid[0], grid[-1]
        poor = self.references.poor_limit_ev
        rich = self.references.rich_limit_ev

        fig = Figure(figsize=(7.6, 5.6), facecolor=_SURFACE)
        ax, strip = fig.subplots(
            2, 1, sharex=True, gridspec_kw={"height_ratios": [8, 0.6], "hspace": 0.08}
        )
        for axis in (ax, strip):
            axis.set_facecolor(_SURFACE)
            if poor is not None and poor > low:
                axis.axvspan(low, poor, color=_OUTSIDE, lw=0, zorder=0)
            if rich < high:
                axis.axvspan(rich, high, color=_OUTSIDE, lw=0, zorder=0)

        styles = {}
        for index, label in enumerate(self.curves):
            styles[label] = {
                "color": _SERIES_COLORS[index % len(_SERIES_COLORS)],
                "linestyle": "-" if index < len(_SERIES_COLORS) else (0, (5, 2)),
            }
            gamma = [point.gamma_ev_per_angstrom2 * factor for point in self.curves[label]]
            ax.plot(grid, gamma, lw=2, zorder=3, label=label, solid_capstyle="round", **styles[label])

        # Headroom above the curves keeps the limit labels clear of the data.
        ymin, ymax = ax.get_ylim()
        ymax += 0.14 * (ymax - ymin)
        for limit, text, align in (
            (poor, self.references.poor_limit_label, "left"),
            (rich, self.references.rich_limit_label, "right"),
        ):
            if limit is not None and low <= limit <= high:
                if low < limit < high:
                    ax.axvline(limit, color=_MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
                ax.text(limit, ymax, text, va="top", ha=align, fontsize=8.5, color=_MUTED)
        for delta_mu, _, _ in self.transitions:
            ax.axvline(delta_mu, color=_AXIS, lw=1, zorder=1)
            ax.annotate(
                f"{delta_mu:.2f} eV".replace("-", "\N{MINUS SIGN}"), xy=(delta_mu, ymin), xytext=(4, 5),
                textcoords="offset points", va="bottom", ha="left", fontsize=8.5, color=_INK_SECONDARY,
            )
        ax.set_ylim(ymin, ymax)

        # Direct labels at the right edge for a few series, nudged apart.
        if len(self.curves) <= 4:
            ends = sorted(
                (points[-1].gamma_ev_per_angstrom2 * factor, label) for label, points in self.curves.items()
            )
            gap, previous = 0.06 * (ymax - ymin), None
            for value, label in ends:
                if previous is not None and value - previous < gap:
                    value = previous + gap
                previous = value
                ax.annotate(
                    label, xy=(high, value), xytext=(6, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=9, color=_INK_SECONDARY, annotation_clip=False,
                )
        if len(self.curves) > 1:
            ax.legend(
                loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=min(len(self.curves), 4),
                frameon=False, fontsize=9, labelcolor=_INK_SECONDARY, handlelength=1.6, borderaxespad=0.2,
            )

        # Stability strip: exact segments of the lower envelope.
        edges = [low] + [delta_mu for delta_mu, _, _ in self.transitions] + [high]
        first = self.transitions[0][1] if self.transitions else self.stable[0]
        owners = [first] + [above for _, _, above in self.transitions]
        for start, end, label in zip(edges, edges[1:], owners):
            for part_start, part_end, inside in _split_by_window(start, end, poor, rich):
                strip.axvspan(
                    part_start, part_end, color=styles[label]["color"], lw=0, zorder=2,
                    alpha=1.0 if inside else 0.35,
                )

        unit_label = "J/m$^2$" if units == "J/m2" else "eV/Å$^2$"
        ax.set_ylabel(f"Surface free energy γ ({unit_label})", color=_INK_SECONDARY)
        if title:
            ax.set_title(title, loc="left", fontsize=11.5, color=_INK, pad=30 if len(self.curves) > 1 else 10)
        ax.grid(axis="y", color=_GRID, lw=0.8)
        ax.set_axisbelow(True)
        strip.set_yticks([])
        strip.set_ylabel("stable", rotation=0, ha="right", va="center", color=_MUTED, fontsize=9)
        strip.set_xlabel(self.references.axis_label, color=_INK_SECONDARY)
        strip.set_xlim(low, high)
        for axis in (ax, strip):
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            axis.spines["bottom"].set_color(_AXIS)
            axis.spines["left"].set_color(_AXIS)
            axis.tick_params(colors=_MUTED, length=0)
        strip.spines["left"].set_visible(False)
        return fig

    def plot(
        self, path: str | Path, *, units: Units = "J/m2", title: str | None = None, dpi: int = 200
    ) -> Path:
        """Draw the phase diagram and save it; the format follows the suffix (png, pdf, svg)."""

        path = Path(path)
        self.figure(units=units, title=title).savefig(
            path, dpi=dpi, bbox_inches="tight", facecolor=_SURFACE
        )
        return path


def surface_phase_diagram(
    terminations: Iterable[SlabTermination],
    references: BinaryOxideReferences | BinaryReferences,
    *,
    delta_mu_range: tuple[float, float] | None = None,
    points: int = 201,
    pseudo_hydrogen: Any = None,
) -> SurfacePhaseDiagram:
    """Evaluate gamma(Delta mu) for every termination on a common grid.

    Delta mu is that of ``references.variable`` (O for an oxide). The default
    range is the stability window of the bulk compound, from the poor limit
    (requires the reference of the other element) to ``Delta mu = 0``. A
    wider ``delta_mu_range`` is allowed; the window limits are then added to
    the grid so that they appear exactly in the CSV export.

    Polar slabs with a passivated bottom (``surfaces=1`` and pseudo-hydrogen
    counts) need ``pseudo_hydrogen``, a
    :class:`psteros.PseudoHydrogenReferences`; their gamma is the absolute
    energy of the top face,
    ``[E_slab - sum_i n_i mu_i - sum_k n_k muhat_k] / A``, on the same scale
    as symmetric slabs. Passivated slabs of one face must share one bottom
    (equal ``bottom_fingerprint``).
    """

    terminations = tuple(terminations)
    if not terminations:
        raise ValueError("at least one termination is required")
    labels = [termination.label for termination in terminations]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"termination labels must be unique: {duplicates}")
    other, variable = references.other, references.variable
    _check_passivated(terminations, pseudo_hydrogen, {other, variable})
    allowed = {other, variable}
    for termination in terminations:
        foreign = sorted(set(termination.composition).difference(allowed))
        if foreign or other not in termination.composition or variable not in termination.composition:
            raise ValueError(
                f"termination {termination.label!r} has composition {termination.composition}; "
                f"expected only {other} and {variable} with both present"
            )
    if points < 2:
        raise ValueError("points must be at least 2")
    if delta_mu_range is None:
        if references.poor_limit_ev is None:
            raise ValueError(
                f"delta_mu_range is required when no reference energy of {other} is given"
            )
        delta_mu_range = (references.poor_limit_ev, references.rich_limit_ev)
    low, high = (float(value) for value in delta_mu_range)
    if not low < high:
        raise ValueError(f"delta_mu_range must be increasing, got {delta_mu_range!r}")

    step = (high - low) / (points - 1)
    grid = [low + index * step for index in range(points - 1)] + [high]
    for limit in (references.poor_limit_ev, references.rich_limit_ev):
        if limit is not None and low < limit < high and not any(isclose(limit, x, abs_tol=1e-12) for x in grid):
            grid.append(limit)
    grid.sort()

    bulk_energy = references.bulk_energy_per_formula_unit_ev

    def gamma(termination: SlabTermination, delta_mu: float) -> SurfaceEnergyPoint:
        correction = 0.0
        if termination.is_passivated:
            correction = pseudo_hydrogen.reservoir_energy_ev(
                termination.pseudo_hydrogen, references.chemical_potentials_ev(delta_mu)
            )
        return surface_energy_binary_equilibrium(
            slab_energy_ev=termination.slab_energy_ev,
            n_other=termination.composition[other],
            n_variable=termination.composition[variable],
            bulk_formula_energy_ev=bulk_energy,
            variable_reference_energy_ev=references.variable_reference_energy_ev,
            delta_mu_ev=delta_mu,
            surface_area_angstrom2=termination.surface_area_angstrom2,
            surfaces=termination.surfaces,
            formula_unit=references.formula_unit,
            reservoir_correction_ev=correction,
        )

    curves = {
        termination.label: tuple(gamma(termination, delta_mu) for delta_mu in grid)
        for termination in terminations
    }
    stable = tuple(label for _, label, _ in stable_termination(curves))
    # gamma is linear in Delta mu: intercept at 0 and slope from the endpoints.
    lines = [
        (
            termination.label,
            gamma(termination, 0.0).gamma_ev_per_angstrom2,
            (gamma(termination, high).gamma_ev_per_angstrom2 - gamma(termination, low).gamma_ev_per_angstrom2)
            / (high - low),
        )
        for termination in terminations
    ]
    return SurfacePhaseDiagram(
        references=references,
        terminations=terminations,
        curves=curves,
        stable=stable,
        transitions=tuple(_lower_envelope_transitions(lines, low, high)),
    )


def _check_passivated(terminations, pseudo_hydrogen, elements: set[str]) -> None:
    """Passivated slabs need pseudo-H references, known elements and, per face, one shared bottom."""

    passivated = [termination for termination in terminations if termination.is_passivated]
    if not passivated:
        return
    if pseudo_hydrogen is None:
        raise ValueError(
            "terminations with pseudo-hydrogen need pseudo_hydrogen=PseudoHydrogenReferences(...): "
            + ", ".join(t.label for t in passivated)
        )
    faces: dict[str, dict[str, list[str]]] = {}
    for termination in passivated:
        foreign = sorted(set(termination.pseudo_hydrogen).difference(elements))
        if foreign:
            raise ValueError(f"termination {termination.label!r} passivates {foreign}, not in the compound")
        if not termination.bottom_fingerprint:
            raise ValueError(
                f"termination {termination.label!r} has pseudo-hydrogen but no bottom_fingerprint; build it "
                "with psteros.find_polar_terminations and SlabTermination.from_polar"
            )
        key = termination.face or "?"
        faces.setdefault(key, {}).setdefault(termination.bottom_fingerprint, []).append(termination.label)
    for face, bottoms in faces.items():
        if len(bottoms) > 1:
            listing = "; ".join(f"{fingerprint}: {', '.join(labels)}" for fingerprint, labels in bottoms.items())
            raise ValueError(
                f"passivated slabs of {face} must share one bottom, but they have {len(bottoms)} ({listing}). "
                "Build all of them in one find_polar_terminations call."
            )


def _lower_envelope_transitions(
    lines: list[tuple[str, float, float]], low: float, high: float, tolerance: float = 1e-12
) -> list[tuple[float, str, str]]:
    """Exact crossings of the lower envelope of ``gamma = intercept + slope * x`` on (low, high)."""

    def value(line: tuple[str, float, float], x: float) -> float:
        return line[1] + line[2] * x

    lowest = min(value(line, low) for line in lines)
    # Among lines tied at the left edge, the smallest slope stays lowest.
    current = min(
        (line for line in lines if value(line, low) <= lowest + tolerance), key=lambda line: line[2]
    )
    x, transitions = low, []
    while True:
        best = None
        for line in lines:
            if line[2] >= current[2] - tolerance:
                continue
            crossing = (line[1] - current[1]) / (current[2] - line[2])
            if crossing <= x + tolerance or crossing >= high - tolerance:
                continue
            if best is None or crossing < best[0] - tolerance or (
                abs(crossing - best[0]) <= tolerance and line[2] < best[1][2]
            ):
                best = (crossing, line)
        if best is None:
            return transitions
        transitions.append((best[0], current[0], best[1][0]))
        x, current = best


def _split_by_window(start: float, end: float, poor: float | None, rich: float):
    """Split ``[start, end]`` at the stability-window limits, flagging the inside parts."""

    cuts = sorted({start, end, *(c for c in (poor, rich) if c is not None and start < c < end)})
    for part_start, part_end in zip(cuts, cuts[1:]):
        middle = (part_start + part_end) / 2
        inside = middle <= rich and (poor is None or middle >= poor)
        yield part_start, part_end, inside


def _unit_conversion(units: str) -> tuple[float, str]:
    if units == "J/m2":
        return EV_PER_ANGSTROM2_TO_J_PER_M2, "Jm2"
    if units == "eV/A2":
        return 1.0, "eVA2"
    raise ValueError(f"units must be 'J/m2' or 'eV/A2', got {units!r}")
