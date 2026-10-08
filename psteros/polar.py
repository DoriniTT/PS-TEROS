"""Absolute surface energies of polar surfaces with pseudo-hydrogen.

A polar slab of a tetrahedrally bonded compound (zinc blende (111), wurtzite
(0001)) has two different faces. Its bottom is passivated with
pseudo-hydrogen: hydrogen with a fractional nuclear charge that completes the
electron count of each broken bond, so the bottom is closed-shell and
electronically separate from the top. The energy of the passivated bottom is
carried by the pseudo chemical potential of the pseudo-hydrogen, and the top
face gets an absolute surface energy (Zhang et al., Sci. Rep. 6, 20055
(2016); Zhang et al., arXiv:1510.08961):

    gamma_top = [E_slab - sum_i n_i mu_i - sum_k n_k muhat_k] / A

This module provides the pseudo-hydrogen model, the slab and reference
builders, and the data needed by :func:`psteros.surface_phase_diagram`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

#: Coordination of every atom in the bulk of the supported structures.
TETRAHEDRAL_COORDINATION = 4


def valence_electrons(element: str) -> int:
    """Valence electrons used in bond counting (group number for s and p elements).

    Group 12 elements (Zn, Cd, Hg) count 2. Transition metals of groups 3-11
    have no unique count and are rejected.
    """

    from pymatgen.core import Element

    group = Element(element).group
    if group in (1, 2):
        return group
    if group == 12:
        return 2
    if 13 <= group <= 18:
        return group - 10
    raise ValueError(
        f"{element} (group {group}) has no well-defined number of bonding electrons; "
        "pseudo-hydrogen passivation supports s- and p-block elements and group 12"
    )


def pseudo_hydrogen_charge(element: str, coordination: int = TETRAHEDRAL_COORDINATION) -> float:
    """Nuclear charge of the pseudo-hydrogen that saturates one broken bond of ``element``.

    Each of the ``coordination`` bonds of an atom with Z valence electrons
    holds Z/N electrons from it; the pseudo-hydrogen supplies the rest of the
    pair, 2 - Z/N. For four-fold coordination: 0.5 on O and S, 0.75 on N,
    P and As, 1.0 (real hydrogen) on Si, 1.25 on Ga, 1.5 on Zn.
    """

    if coordination <= 0:
        raise ValueError("coordination must be positive")
    charge = 2.0 - valence_electrons(element) / coordination
    if not 0.0 < charge < 2.0:
        raise ValueError(f"no pseudo-hydrogen for {element}: charge {charge:g} is outside (0, 2)")
    return charge


def _charge_text(charge: float) -> str:
    """0.75 -> '.75', 1.25 -> '1.25', 1.0 -> '' (VASP POTCAR naming)."""

    if abs(charge - 1.0) < 1e-9:
        return ""
    text = f"{charge:.2f}".rstrip("0").rstrip(".")
    return text[1:] if text.startswith("0") else text


@dataclass(frozen=True)
class PseudoHydrogen:
    """The pseudo-hydrogen bonded to one element of the compound.

    Attributes:
        bonded_to: Element it passivates (the bottom-surface atom).
        charge: Nuclear charge, 2 - Z/4.
        formal_charge: Its share of the electron count, minus the oxidation
            state of ``bonded_to`` divided by the coordination (+0.75 on As,
            -0.75 on Ga). A slab whose formal charges, pseudo-hydrogen
            included, add up to zero satisfies the electron-counting rule.
    """

    bonded_to: str
    charge: float
    formal_charge: float

    @property
    def kind_name(self) -> str:
        """AiiDA kind name, e.g. ``H0p75`` (``H`` for real hydrogen)."""

        if abs(self.charge - 1.0) < 1e-9:
            return "H"
        return "H" + f"{self.charge:.2f}".rstrip("0").rstrip(".").replace(".", "p")

    @property
    def vasp_potential(self) -> str:
        """Name of the VASP POTCAR, e.g. ``H.75``, ``H1.25`` (``H`` for real hydrogen)."""

        return "H" + _charge_text(self.charge)

    @property
    def label(self) -> str:
        """Short text label, e.g. ``H0.75(As)``."""

        return f"H{self.charge:g}({self.bonded_to})"

    def site_properties(self) -> dict[str, Any]:
        """pymatgen site properties for this pseudo-hydrogen (``kind_name`` is read by AiiDA)."""

        return {"kind_name": self.kind_name, "pseudo_hydrogen": self.bonded_to}


def pseudo_hydrogens(
    bulk: Any,
    oxidation_states: Mapping[str, float] | None = None,
    *,
    check_coordination: bool = True,
) -> dict[str, PseudoHydrogen]:
    """The pseudo-hydrogen for each element of a tetrahedrally bonded compound.

    Args:
        bulk: pymatgen ``Structure`` of the bulk (zinc blende, wurtzite,
            diamond, ...).
        oxidation_states: Element -> oxidation state; guessed if omitted
            (Ga +3, As -3).
        check_coordination: Reject structures in which an atom does not have
            four nearest neighbours.

    Returns:
        ``{element: PseudoHydrogen}``.
    """

    from psteros.core.terminations import resolve_oxidation_states

    states = resolve_oxidation_states(bulk, oxidation_states)
    if check_coordination:
        coordination = tetrahedral_coordination(bulk)
        wrong = {element: count for element, count in coordination.items() if count != TETRAHEDRAL_COORDINATION}
        if wrong:
            raise ValueError(
                f"pseudo-hydrogen passivation needs four-fold coordinated atoms; {bulk.composition.reduced_formula} "
                f"has coordination {wrong}"
            )
    return {
        element: PseudoHydrogen(
            bonded_to=element,
            charge=pseudo_hydrogen_charge(element),
            formal_charge=-states[element] / TETRAHEDRAL_COORDINATION + 0.0,
        )
        for element in states
    }


def default_hydrogen_bond_length(element: str) -> float:
    """Initial X-H distance (A): sum of the atomic radii of X and H; relaxed afterwards."""

    from pymatgen.core import Element

    return float(Element(element).atomic_radius or 1.0) + float(Element("H").atomic_radius or 0.25)


def tetrahedral_coordination(bulk: Any, tolerance: float = 0.15) -> dict[str, int]:
    """Number of nearest neighbours of each element (within ``tolerance`` of the shortest bond)."""

    neighbours = bulk.get_all_neighbors(4.0)
    shortest = min(neighbour.nn_distance for site in neighbours for neighbour in site)
    counts: dict[str, set[int]] = {}
    for site, found in zip(bulk, neighbours):
        count = sum(1 for neighbour in found if neighbour.nn_distance <= shortest * (1 + tolerance))
        counts.setdefault(site.specie.symbol, set()).add(count)
    return {element: (values.pop() if len(values) == 1 else -1) for element, values in counts.items()}


# =============================================================================
# POLAR SLABS WITH A PASSIVATED BOTTOM
# =============================================================================

#: Number of atomic planes of the default slab: 9 bilayers for a polar
#: direction (the Sci. Rep. setting), 12 layers for a non-polar one.
DEFAULT_BILAYERS = 9
DEFAULT_NONPOLAR_LAYERS = 12


@dataclass(frozen=True)
class PolarTermination:
    """One slab with a fixed, pseudo-hydrogen passivated bottom.

    Attributes:
        structure: pymatgen ``Structure`` with c normal to the surface. Site
            properties: ``kind_name`` (AiiDA kind; ``H0p75`` etc. for
            pseudo-hydrogen), ``pseudo_hydrogen`` (element it passivates, or
            ``None``) and ``bottom`` (True for the bottom reference region).
        label: ``term_0``, ``term_1``, ...
        miller_index: Orientation; the top face is this (hkl) face.
        supercell: Surface cell, the same for every slab of one set.
        top_element: Species of the ideal top plane.
        bottom_element: Species of the bottom plane, carrying the pseudo-H.
        electron_counting: Whether the formal charges, pseudo-hydrogen
            included, add up to zero (the electron-counting rule).
        removed_per_cell: Formulas of the top atoms removed from the ideal top.
        thickness: Distance between the outermost nuclei, pseudo-H included.
        bottom_fingerprint: Hash of the bottom region and the surface cell;
            equal for slabs that share the same bottom.
        removed_sites: Species and position of the removed atoms (for plots).
    """

    structure: Any
    label: str
    miller_index: tuple[int, int, int]
    supercell: tuple[int, int]
    top_element: str
    bottom_element: str
    electron_counting: bool
    removed_per_cell: tuple[str, ...]
    thickness: float
    bottom_fingerprint: str
    removed_sites: tuple = ()
    bulk_reduced_composition: Mapping[str, float] | None = None
    bulk_formula: str = ""

    @property
    def composition(self) -> dict[str, int]:
        """Element counts without the pseudo-hydrogen."""

        counts: dict[str, int] = {}
        for site, passivates in zip(self.structure, self.structure.site_properties["pseudo_hydrogen"]):
            if passivates is None:
                counts[site.specie.symbol] = counts.get(site.specie.symbol, 0) + 1
        return counts

    @property
    def pseudo_hydrogen_counts(self) -> dict[str, int]:
        """Number of pseudo-hydrogen atoms, by the element they passivate."""

        counts: dict[str, int] = {}
        for passivates in self.structure.site_properties["pseudo_hydrogen"]:
            if passivates is not None:
                counts[passivates] = counts.get(passivates, 0) + 1
        return counts

    @property
    def bottom_indices(self) -> tuple[int, ...]:
        """Indices of the sites in the bottom reference region (pseudo-H included)."""

        return tuple(i for i, flag in enumerate(self.structure.site_properties["bottom"]) if flag)

    @property
    def formula(self) -> str:
        return "".join(f"{element}{count}" for element, count in self.composition.items())

    @property
    def is_stoichiometric(self) -> bool:
        if not self.bulk_reduced_composition:
            return False
        from pymatgen.core import Composition

        return Composition(self.composition).reduced_composition == Composition(self.bulk_reduced_composition)

    @property
    def area(self) -> float:
        import numpy as np

        matrix = self.structure.lattice.matrix
        return float(np.linalg.norm(np.cross(matrix[0], matrix[1])))

    @property
    def origin(self) -> str:
        cell = "x".join(str(value) for value in self.supercell)
        if not self.removed_per_cell:
            return f"ideal {self.top_element}-terminated"
        from psteros.core.terminations import _count_formulas

        removed = " + ".join(_count_formulas(self.removed_per_cell))
        return f"{self.top_element}-terminated minus {removed} per {cell} cell"

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "formula": self.formula,
            "pseudo_hydrogen": self.pseudo_hydrogen_counts,
            "n_atoms": len(self.structure),
            "miller_index": list(self.miller_index),
            "supercell": list(self.supercell),
            "top_element": self.top_element,
            "bottom_element": self.bottom_element,
            "electron_counting": self.electron_counting,
            "removed_per_cell": list(self.removed_per_cell),
            "thickness": self.thickness,
            "area": self.area,
            "bottom_fingerprint": self.bottom_fingerprint,
            "bottom_indices": list(self.bottom_indices),
        }


class PolarTerminationSet(list):
    """Slabs of one face that share one passivated bottom.

    Behaves like a list of :class:`PolarTermination`. Printing it shows a
    table; :meth:`write` saves the slabs, :meth:`plot` draws side views.
    """

    COLUMNS = ("label", "formula", "pseudo-H", "atoms", "thickness (Å)", "e-count", "origin")

    def __init__(self, terminations=(), *, bulk_formula: str = "", miller_index=(0, 0, 0), supercell=(1, 1),
                 oxidation_states=None, hydrogens=None, polar: bool = True, message: str = ""):
        super().__init__(terminations)
        self.bulk_formula = bulk_formula
        self.miller_index = tuple(miller_index)
        self.supercell = tuple(supercell)
        self.oxidation_states = dict(oxidation_states or {})
        self.hydrogens = dict(hydrogens or {})
        self.polar = polar
        self.message = message

    @property
    def bottom_fingerprint(self) -> str:
        fingerprints = {termination.bottom_fingerprint for termination in self}
        if len(fingerprints) != 1:
            raise ValueError(f"terminations do not share one bottom: {sorted(fingerprints)}")
        return fingerprints.pop()

    def header(self) -> str:
        from psteros.core.terminations import hkl_label

        first = self[0] if self else None
        cell = "x".join(str(value) for value in self.supercell)
        parts = [f"{self.bulk_formula}{hkl_label(self.miller_index)}", "polar" if self.polar else "non-polar",
                 f"{cell} surface cell"]
        if first is not None:
            hydrogen = self.hydrogens.get(first.bottom_element)
            passivation = f" + {hydrogen.label}" if hydrogen else ""
            parts.append(f"top {first.top_element}, bottom {first.bottom_element}{passivation} (shared)")
        return " | ".join(parts)

    def rows(self) -> list[list[str]]:
        rows = []
        for termination in self:
            hydrogen = ", ".join(
                f"{count} {self.hydrogens[element].label}" if element in self.hydrogens else f"{count} H"
                for element, count in termination.pseudo_hydrogen_counts.items()
            )
            rows.append([
                termination.label, termination.formula, hydrogen, str(len(termination.structure)),
                f"{termination.thickness:.2f}", "yes" if termination.electron_counting else "no",
                termination.origin,
            ])
        return rows

    def summary(self) -> str:
        rows = self.rows()
        widths = [max(len(cell) for cell in column) for column in zip(self.COLUMNS, *rows)]

        def line(cells):
            return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip()

        text = [self.header(), "", line(self.COLUMNS), line(["-" * width for width in widths])]
        text += [line(row) for row in rows]
        noun = "slab" if len(self) == 1 else "slabs"
        text += ["", f"{len(self)} {noun} on one bottom (fingerprint {self[0].bottom_fingerprint})."
                 if self else "No slabs."]
        if self.message:
            text.append(self.message)
        return "\n".join(text)

    __str__ = summary
    __repr__ = summary

    def _repr_html_(self) -> str:
        import html

        head = "".join(f"<th style='text-align:left'>{html.escape(c)}</th>" for c in self.COLUMNS)
        body = "".join("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in row) + "</tr>" for row in self.rows())
        return f"<p><b>{html.escape(self.header())}</b></p><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"

    def to_dicts(self) -> list[dict]:
        return [termination.to_dict() for termination in self]

    def write(self, directory: str) -> list[str]:
        """Write each slab as a VASP POSCAR and a ``terminations.json`` summary.

        Pseudo-hydrogen of different charges are written as separate species
        groups; ``terminations.json`` lists the POTCAR of each group in order.
        """

        import json
        import os

        os.makedirs(directory, exist_ok=True)
        paths, potcars = [], {}
        for termination in self:
            path = os.path.join(directory, f"{termination.label}_{termination.formula}.vasp")
            potcars[termination.label] = write_poscar(termination.structure, path, comment=termination.origin)
            paths.append(path)
        summary = os.path.join(directory, "terminations.json")
        with open(summary, "w") as handle:
            json.dump({"header": self.header(), "terminations": self.to_dicts(), "potcar_order": potcars},
                      handle, indent=2)
        return paths + [summary]

    def plot(self, filename: str | None = None, repeat: int = 2, **kwargs):
        from psteros.core.terminations import plot_terminations

        return plot_terminations(self, filename=filename, repeat=repeat, **kwargs)


def write_poscar(structure: Any, path: str, comment: str = "") -> list[str]:
    """Write a POSCAR with one species group per AiiDA kind and return the group order.

    Each entry of the returned list is the VASP POTCAR of one group, in order
    (``H.75``, ``H1.25`` for pseudo-hydrogen, the element symbol otherwise).
    """

    kinds = structure.site_properties.get("kind_name") or [site.specie.symbol for site in structure]
    order: list[str] = []
    for kind in kinds:
        if kind not in order:
            order.append(kind)
    lattice = structure.lattice.matrix
    lines = [comment or structure.composition.formula, "1.0"]
    lines += ["  " + " ".join(f"{value:.10f}" for value in row) for row in lattice]
    groups = [[i for i, kind in enumerate(kinds) if kind == name] for name in order]
    lines.append("  " + " ".join(structure[group[0]].specie.symbol for group in groups))
    lines.append("  " + " ".join(str(len(group)) for group in groups))
    lines.append("Cartesian")
    for group in groups:
        for index in group:
            lines.append("  " + " ".join(f"{value:.10f}" for value in structure[index].coords) + f"  {kinds[index]}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    potcars = []
    for group, name in zip(groups, order):
        passivates = structure.site_properties.get("pseudo_hydrogen", [None] * len(structure))[group[0]]
        potcars.append(_potcar_for_kind(name) if passivates is not None else structure[group[0]].specie.symbol)
    return potcars


def _potcar_for_kind(kind_name: str) -> str:
    """'H0p75' -> 'H.75', 'H1p25' -> 'H1.25', 'H' -> 'H'."""

    if kind_name == "H":
        return "H"
    return "H" + _charge_text(float(kind_name[1:].replace("p", ".")))


def find_polar_terminations(
    bulk: Any,
    miller_index,
    *,
    bilayers: int | None = None,
    layers: int | None = None,
    vacuum: float = 15.0,
    oxidation_states: Mapping[str, float] | None = None,
    electron_counting: bool = True,
    include_ideal: bool = True,
    supercell: tuple[int, int] | None = None,
    hydrogen_bond_lengths: Mapping[str, float] | None = None,
    max_variants: int = 20,
    plane_tolerance: float = 0.3,
    symprec: float = 0.1,
    lll_reduce: bool = True,
    max_normal_search: int | None = None,
) -> PolarTerminationSet:
    """Slabs of one (hkl) face built on one shared, pseudo-hydrogen passivated bottom.

    The top face is the (hkl) face; use (-h-k-l) for the opposite polar
    face. Both cuts break the fewest bonds (one per surface atom for zinc
    blende (111) and wurtzite (0001)), and every broken bond of the bottom
    is saturated by a pseudo-hydrogen along the bond.

    Args:
        bulk: Relaxed bulk (pymatgen ``Structure``) of a tetrahedrally
            bonded compound.
        miller_index: Surface orientation of the top face.
        bilayers: Slab thickness in bilayers (two atomic planes each);
            default 9 for polar directions.
        layers: Slab thickness in atomic planes; default 12 for non-polar
            directions (useful to validate the method against a symmetric
            slab).
        vacuum: Vacuum in A.
        oxidation_states: Element -> oxidation state; guessed if omitted.
        electron_counting: Also build tops that satisfy the electron-counting
            rule by removing top atoms (e.g. the Ga vacancy (111)-2x2).
        include_ideal: Include the ideal, unreconstructed top.
        supercell: Surface cell. By default the smallest of 1x1, 2x1, 1x2
            and 2x2 in which the electron-counting tops exist (1x1 when none
            are requested). All slabs of the set use the same cell.
        hydrogen_bond_lengths: Element -> X-H distance (A) for the initial
            pseudo-H positions; default: sum of atomic radii. Positions are
            relaxed afterwards.
        max_variants: Largest number of symmetry-distinct electron-counting
            tops kept.

    Returns:
        A :class:`PolarTerminationSet` whose slabs share one bottom
        fingerprint, labelled ``term_0``, ``term_1``, ...
    """

    import warnings

    from psteros.core.terminations import (
        AUTO_SUPERCELLS,
        has_face_reversing_operation,
        normal_repeat,
        resolve_oxidation_states,
    )

    miller = tuple(int(value) for value in miller_index)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")
        states = resolve_oxidation_states(bulk, oxidation_states)
        hydrogens = pseudo_hydrogens(bulk, states)
        polar = not has_face_reversing_operation(bulk, miller, symprec)
        if layers is None:
            layers = 2 * (bilayers or DEFAULT_BILAYERS) if (polar or bilayers) else DEFAULT_NONPOLAR_LAYERS
        repeat = normal_repeat(bulk, miller, symprec)
        cells = [(1, 1)] if not electron_counting and supercell is None else (
            AUTO_SUPERCELLS if supercell is None else [(int(supercell[0]), int(supercell[1]))]
        )
        message = ""
        for cell in cells:
            ideal = _ideal_slab(bulk, miller, layers, cell, vacuum, states, hydrogens, repeat,
                                hydrogen_bond_lengths, plane_tolerance, symprec, lll_reduce, max_normal_search)
            variants = []
            if electron_counting and abs(ideal.charge) > 1e-6:
                variants = _top_removals(ideal, repeat, max_variants, symprec)
                if not variants and supercell is None and cell != cells[-1]:
                    continue
                if not variants:
                    message = (f"No electron-counting top in the {'x'.join(map(str, cell))} cell: the ideal top "
                               f"has formal charge {ideal.charge:+g}.")
            break

        built = []
        if include_ideal or abs(ideal.charge) <= 1e-6:
            built.append(ideal.termination(frozenset(), abs(ideal.charge) <= 1e-6))
        for removed in variants:
            built.append(ideal.termination(removed, True))
    labelled = [_relabel(termination, f"term_{index}") for index, termination in enumerate(built)]
    return PolarTerminationSet(
        labelled, bulk_formula=bulk.composition.reduced_formula, miller_index=miller, supercell=cell,
        oxidation_states=states, hydrogens=hydrogens, polar=polar, message=message,
    )


def _relabel(termination: PolarTermination, label: str) -> PolarTermination:
    from dataclasses import replace

    return replace(termination, label=label)


class _IdealSlab:
    """The ideal slab of one cell: atoms, pseudo-H and the bookkeeping for top removals."""

    def __init__(self, *, lattice, species, coords, kinds, passivates, plane, bottom, charges, miller, cell,
                 top_element, bottom_element, thickness, bulk_reduced):
        self.lattice, self.species, self.coords = lattice, species, coords
        self.kinds, self.passivates, self.plane, self.bottom = kinds, passivates, plane, bottom
        self.charges, self.miller, self.cell = charges, miller, cell
        self.top_element, self.bottom_element, self.thickness = top_element, bottom_element, thickness
        self.bulk_reduced = bulk_reduced
        from pymatgen.core import Composition

        self.bulk_formula = Composition(bulk_reduced).reduced_formula
        self.charge = float(sum(charges))
        self.fingerprint = _bottom_fingerprint(lattice, species, coords, kinds, bottom)

    def structure(self, keep):
        from pymatgen.core import Structure

        indices = [i for i in range(len(self.species)) if i in keep]
        return Structure(
            self.lattice, [self.species[i] for i in indices], [self.coords[i] for i in indices],
            coords_are_cartesian=True,
            site_properties={
                "kind_name": [self.kinds[i] for i in indices],
                "pseudo_hydrogen": [self.passivates[i] for i in indices],
                "bottom": [bool(self.bottom[i]) for i in indices],
            },
        )

    def termination(self, removed: frozenset, electron_counting: bool) -> PolarTermination:
        keep = set(range(len(self.species))).difference(removed)
        structure = self.structure(keep)
        import numpy as np

        normal = np.asarray(self.lattice.matrix[2]) / np.linalg.norm(self.lattice.matrix[2])
        heights = structure.cart_coords @ normal
        removed_formulas = tuple(sorted(self.species[i] for i in removed))
        ghosts = tuple((self.species[i], tuple(float(x) for x in self.coords[i])) for i in sorted(removed))
        return PolarTermination(
            structure=structure, label="", miller_index=self.miller, supercell=self.cell,
            top_element=self.top_element, bottom_element=self.bottom_element,
            electron_counting=electron_counting, removed_per_cell=removed_formulas,
            thickness=float(heights.max() - heights.min()), bottom_fingerprint=self.fingerprint,
            removed_sites=ghosts, bulk_reduced_composition=self.bulk_reduced,
            bulk_formula=self.bulk_formula,
        )


def _bottom_fingerprint(lattice, species, coords, kinds, bottom) -> str:
    """Hash of the bottom region (kinds and positions relative to its lowest atom) and the cell."""

    import hashlib

    import numpy as np

    normal = np.asarray(lattice.matrix[2]) / np.linalg.norm(lattice.matrix[2])
    indices = [i for i, flag in enumerate(bottom) if flag]
    relative = np.array([coords[i] for i in indices])
    relative = relative - relative[np.argmin(relative @ normal)]
    # In-plane position wrapped into the cell, and height above the lowest
    # bottom site: independent of the slab thickness and of the vacuum.
    frac = lattice.get_fractional_coords(relative)[:, :2]
    frac -= np.floor(frac + 1e-6)
    in_plane = frac @ np.asarray(lattice.matrix[:2])
    heights = relative @ normal
    rows = sorted(
        (kinds[i], *(round(float(v), 3) + 0.0 for v in (*xy, height)))
        for i, xy, height in zip(indices, in_plane, heights)
    )
    cell = [round(float(v), 4) + 0.0 for v in np.asarray(lattice.matrix[:2]).ravel()]
    return hashlib.sha256(repr((cell, rows)).encode()).hexdigest()[:16]


def _ideal_slab(bulk, miller, layers, cell, vacuum, states, hydrogens, repeat, hydrogen_bond_lengths,
                plane_tolerance, symprec, lll_reduce, max_normal_search) -> _IdealSlab:
    import numpy as np
    from pymatgen.core import Element, Lattice

    from psteros.core.terminations import _oriented_parent

    neighbours = bulk.get_all_neighbors(4.0)
    bond = min(neighbour.nn_distance for site in neighbours for neighbour in site)
    thickness = (layers + 8) * bond
    while True:
        parent = _oriented_parent(bulk, miller, thickness, vacuum, lll_reduce, max_normal_search, cell,
                                  repeat, symprec)
        normal = parent.lattice.matrix[2] / np.linalg.norm(parent.lattice.matrix[2])
        heights = parent.cart_coords @ normal
        order = np.argsort(heights)
        planes: list[list[int]] = []
        for index in order:
            if planes and heights[index] - heights[planes[-1][-1]] < plane_tolerance:
                planes[-1].append(int(index))
            else:
                planes.append([int(index)])
        if len(planes) >= layers + 6:
            break
        thickness *= 1.5

    plane_of = np.empty(len(parent), dtype=int)
    for number, members in enumerate(planes):
        plane_of[members] = number
    bonds = []  # (i, j, vector from i to j) with plane(i) < plane(j)
    centres, points, images, distances = parent.get_neighbor_list(bond * 1.15)
    for i, j, image, distance in zip(centres, points, images, distances):
        if plane_of[i] < plane_of[j]:
            vector = parent.lattice.get_cartesian_coords(parent[j].frac_coords + image - parent[i].frac_coords)
            bonds.append((int(i), int(j), vector))
    crossing = [0] * (len(planes) - 1)
    for i, j, _ in bonds:
        for gap in range(plane_of[i], plane_of[j]):
            crossing[gap] += 1
    usable = range(2, len(planes) - layers - 2)
    fewest = min(crossing[gap] for gap in usable)
    bottom_gap = next(gap for gap in usable if crossing[gap] == fewest)
    top_gap = bottom_gap + layers
    if crossing[top_gap] != fewest:
        raise ValueError(
            f"{layers} planes do not end on a cut with the fewest broken bonds; use an even number of "
            "planes (whole bilayers) for this orientation"
        )
    window = [index for plane in planes[bottom_gap + 1:top_gap + 1] for index in plane]
    window_set = set(window)

    # SlabGenerator may return the stack upside down relative to (hkl). The
    # (hkl) face is the one whose cut bonds point along +n in the bulk; when
    # that is the parent's lower side, passivate the upper side instead and
    # turn the slab over.
    cut = [(parent[i].specie.symbol, parent[j].specie.symbol, float(vector @ normal))
           for i, j, vector in bonds if plane_of[i] <= bottom_gap < plane_of[j]]
    flip = _parent_is_reversed(bulk, miller, cut)
    # (surface atom inside the slab, its missing neighbour, vector from the first to the second)
    if flip:
        passivated = [(i, j, vector) for i, j, vector in bonds
                      if i in window_set and j not in window_set and plane_of[j] > top_gap]
    else:
        passivated = [(j, i, -vector) for i, j, vector in bonds
                      if j in window_set and i not in window_set and plane_of[i] <= bottom_gap]

    lengths = dict(hydrogen_bond_lengths or {})
    species, coords, kinds, passivates, charges = [], [], [], [], []
    for index in window:
        symbol = parent[index].specie.symbol
        species.append(symbol)
        coords.append(parent[index].coords)
        kinds.append(symbol)
        passivates.append(None)
        charges.append(states[symbol])
    for inside, _, outward in passivated:
        # outward: from the surface atom towards its missing neighbour
        symbol = parent[inside].specie.symbol
        hydrogen = hydrogens[symbol]
        length = float(lengths.get(symbol, default_hydrogen_bond_length(symbol)))
        species.append("H")
        coords.append(parent[inside].coords + outward / np.linalg.norm(outward) * length)
        kinds.append(hydrogen.kind_name)
        passivates.append(symbol)
        charges.append(hydrogen.formal_charge)
    coords = np.array(coords)
    in_plane = np.array(parent.lattice.matrix[:2], dtype=float)
    if flip:
        axis = in_plane[0] / np.linalg.norm(in_plane[0])
        rotation = 2.0 * np.outer(axis, axis) - np.eye(3)  # 180 degrees about a: proper rotation
        # Rotate the atoms and the in-plane vectors; c keeps pointing up, so
        # the passivated side ends at the bottom.
        coords = coords @ rotation.T
        in_plane = in_plane @ rotation.T
    heights = coords @ normal
    slab_thickness = float(heights.max() - heights.min())
    length = slab_thickness + vacuum
    lattice = Lattice(np.vstack([in_plane, normal * length]))
    coords = coords + normal * (0.5 * length - 0.5 * (heights.max() + heights.min()))
    frac = lattice.get_fractional_coords(coords)
    frac[:, :2] -= np.floor(frac[:, :2])
    coords = lattice.get_cartesian_coords(frac)
    heights = coords @ normal

    # Planes of the slab, numbered upwards; pseudo-H get -1.
    atoms = [i for i, p in enumerate(passivates) if p is None]
    order = sorted(atoms, key=lambda i: heights[i])
    plane = [-1] * len(species)
    number = -1
    previous = None
    for index in order:
        if previous is None or heights[index] - previous >= plane_tolerance:
            number += 1
        plane[index] = number
        previous = heights[index]
    top_species = {species[i] for i in atoms if plane[i] == number}
    bottom_species = {species[i] for i in atoms if plane[i] == 0}
    top_element = "/".join(sorted(top_species))
    bottom_element = "/".join(sorted(bottom_species))
    lowest = min(heights[i] for i in atoms)
    bottom_depth = repeat + plane_tolerance
    bottom = [p is not None or heights[i] - lowest <= bottom_depth for i, p in enumerate(passivates)]
    return _IdealSlab(
        lattice=lattice, species=species, coords=list(coords), kinds=kinds, passivates=passivates,
        plane=plane, bottom=bottom, charges=charges, miller=miller, cell=cell, top_element=top_element,
        bottom_element=bottom_element, thickness=slab_thickness,
        bulk_reduced=dict(bulk.composition.reduced_composition.as_dict()),
    )


def _parent_is_reversed(bulk, miller, cut, tolerance: float = 0.05) -> bool:
    """Whether the parent slab's normal points along -n_(hkl).

    ``cut`` holds ``(lower species, upper species, normal component)`` of the
    bonds crossing the chosen gap of the parent. In the bulk, bonds between
    the same species with the same component along +n_(hkl) mean the parent
    is oriented like (hkl); only along -n_(hkl) means it is reversed. A
    non-polar cut matches both, and is left as it is.
    """

    import numpy as np

    normal = bulk.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(miller)
    normal = normal / np.linalg.norm(normal)
    lower, upper, component = cut[0]
    neighbours_of = bulk.get_all_neighbors(4.0)
    bond = min(neighbour.nn_distance for found in neighbours_of for neighbour in found)
    forward = backward = False
    for site, neighbours in zip(bulk, neighbours_of):
        if site.specie.symbol != lower:
            continue
        for neighbour in neighbours:
            if neighbour.specie.symbol != upper or neighbour.nn_distance > bond * 1.15:
                continue
            value = float((neighbour.coords - site.coords) @ normal)
            forward |= abs(value - component) <= tolerance
            backward |= abs(value + component) <= tolerance
    return backward and not forward


def _top_removals(ideal: _IdealSlab, repeat: float, max_variants: int, symprec: float) -> list[frozenset]:
    """Sets of top atoms whose removal brings the formal charge to zero.

    Whole outer planes go first, then part of the next plane, with the fewest
    atoms; sets related by a symmetry of the ideal slab are kept once.
    """

    import itertools
    import warnings

    import numpy as np

    from psteros.core.terminations import _slab_symmetry

    target = ideal.charge
    top_planes = sorted({p for p in ideal.plane if p >= 0}, reverse=True)
    lowest_top = max(top_planes) - max(1, int(round(repeat / 1.0)))
    candidates = [p for p in top_planes if p >= lowest_top][:4]
    members = {p: [i for i, q in enumerate(ideal.plane) if q == p] for p in candidates}
    solutions: dict[int, list[tuple[int, ...]]] = {}
    for depth, p in enumerate(candidates):
        stripped = tuple(i for q in candidates[:depth] for i in members[q])
        stripped_charge = sum(ideal.charges[i] for i in stripped)
        for count in range(1, len(members[p]) + 1):
            if len(stripped) + count > 4 * ideal.cell[0] * ideal.cell[1]:
                break
            for subset in itertools.combinations(members[p], count):
                total = stripped_charge + sum(ideal.charges[i] for i in subset)
                if abs(total - target) <= 1e-6:
                    solutions.setdefault(len(stripped) + count, []).append(stripped + subset)
    if not solutions:
        return []
    structure = ideal.structure(set(range(len(ideal.species))))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")
        operations = [permutation for _, permutation, reverses in _slab_symmetry(structure, symprec) if not reverses]
    if not operations:
        operations = [np.arange(len(structure))]
    seen, result = set(), []
    for subset in solutions[min(solutions)]:
        key = min(tuple(sorted(int(permutation[i]) for i in subset)) for permutation in operations)
        if key in seen:
            continue
        seen.add(key)
        result.append(frozenset(subset))
        if len(result) >= max_variants:
            break
    return result


# =============================================================================
# PSEUDO CHEMICAL POTENTIALS OF PSEUDO-HYDROGEN
# =============================================================================

@dataclass(frozen=True)
class PseudoHydrogenReference:
    """Pseudo chemical potential of the pseudo-hydrogen bonded to one element.

    ``muhat = constant_ev - mu_X / 4`` for any chemical potential ``mu_X`` of
    the passivated element X (Sci. Rep. Eqs. 8 and 9 both give
    d muhat / d mu_X = -1/4). ``constant_ev`` is E(molecule)/4 for the
    pseudo-molecule method.
    """

    bonded_to: str
    constant_ev: float
    method: str = "pseudo-molecule"

    def value_ev(self, mu_bonded_ev: float) -> float:
        """muhat at the chemical potential ``mu_bonded_ev`` of the passivated element."""

        return self.constant_ev - mu_bonded_ev / TETRAHEDRAL_COORDINATION


class PseudoHydrogenReferences(dict):
    """``{element: PseudoHydrogenReference}`` for the pseudo-hydrogen of one compound."""

    @classmethod
    def from_pseudo_molecules(cls, energies_ev: Mapping[str, float]) -> "PseudoHydrogenReferences":
        """From the total energies of the relaxed pseudo-molecules X(H_X)4, by element X."""

        return cls({
            element: PseudoHydrogenReference(element, float(energy) / TETRAHEDRAL_COORDINATION)
            for element, energy in energies_ev.items()
        })

    def reservoir_energy_ev(self, counts: Mapping[str, int], chemical_potentials_ev: Mapping[str, float]) -> float:
        """sum_k n_k muhat_k for pseudo-hydrogen counts by passivated element."""

        missing = sorted(set(counts).difference(self))
        if missing:
            raise ValueError(f"no pseudo chemical potential for pseudo-hydrogen on {missing}")
        return sum(
            count * self[element].value_ev(chemical_potentials_ev[element]) for element, count in counts.items()
        )


def pseudo_molecule(
    bulk: Any,
    element: str,
    *,
    oxidation_states: Mapping[str, float] | None = None,
    bond_length: float | None = None,
    box: float | None = None,
) -> Any:
    """The CH4-like pseudo-molecule X(H_X)4 used for the pseudo chemical potential of H_X.

    The four pseudo-hydrogen sit at the tetrahedral directions; the
    molecule has eight valence electrons and is relaxed at Gamma in a cubic
    box with at least 15 A between images (``box`` sets the edge).
    """

    import numpy as np
    from pymatgen.core import Element, Lattice, Structure

    hydrogen = pseudo_hydrogens(bulk, oxidation_states)[element]
    if bond_length is None:
        bond_length = default_hydrogen_bond_length(element)
    edge = box or (2 * bond_length + 15.0)
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3.0)
    centre = np.full(3, edge / 2)
    coords = [centre] + [centre + bond_length * direction for direction in directions]
    return Structure(
        Lattice.cubic(edge), [element] + ["H"] * 4, coords, coords_are_cartesian=True,
        site_properties={
            "kind_name": [element] + [hydrogen.kind_name] * 4,
            "pseudo_hydrogen": [None] + [element] * 4,
            "bottom": [False] * 5,
        },
    )


def cluster_counts(size: int) -> dict[str, int]:
    """Atoms and pseudo-hydrogen of the tetrahedral cluster with ``size`` outer atoms per edge (Eq. 9)."""

    if size < 2:
        raise ValueError("cluster size must be at least 2")
    n = size
    return {
        "outer": n * (n + 1) * (n + 2) // 6,
        "inner": (n - 1) * n * (n + 1) // 6,
        "face": 2 * (n - 2) * (n - 3),
        "edge": 12 * (n - 2),
        "corner": 12,
    }


def zinc_blende_analogue(bulk: Any) -> Any:
    """Zinc-blende structure of the same compound with the same nearest-neighbour distance.

    The cluster method uses zinc-blende clusters for wurtzite compounds as
    well (arXiv:1510.08961).
    """

    import numpy as np
    from pymatgen.core import Lattice, Structure

    symbols = [site.specie.symbol for site in bulk]
    elements = sorted(set(symbols), key=symbols.index)
    if len(elements) != 2:
        raise ValueError("the cluster method needs a binary compound")
    neighbours = bulk.get_all_neighbors(4.0)
    bond = min(neighbour.nn_distance for found in neighbours for neighbour in found)
    a = 4.0 * bond / np.sqrt(3.0)
    return Structure.from_spacegroup("F-43m", Lattice.cubic(a), elements, [[0, 0, 0], [0.25, 0.25, 0.25]])


def tetrahedral_cluster(
    bulk: Any,
    outer: str,
    size: int,
    *,
    oxidation_states: Mapping[str, float] | None = None,
    bond_length: float | None = None,
    vacuum: float = 15.0,
) -> Any:
    """Zinc-blende tetrahedral cluster with four (111) facets of ``outer`` atoms, all passivated.

    ``size`` is the number of outer atoms on an edge. Corner atoms carry
    three pseudo-hydrogen, edge atoms two and face atoms one (like a (111)
    surface atom); the counts are those of :func:`cluster_counts`. For a
    wurtzite bulk the zinc-blende analogue is used.
    """

    import numpy as np
    from pymatgen.core import Element, Lattice, Structure

    counts = cluster_counts(size)
    zinc_blende = zinc_blende_analogue(bulk)
    elements = sorted({site.specie.symbol for site in zinc_blende})
    if outer not in elements:
        raise ValueError(f"outer must be one of {elements}, got {outer!r}")
    inner = next(element for element in elements if element != outer)
    hydrogen = pseudo_hydrogens(bulk, oxidation_states)[outer]
    a = zinc_blende.lattice.a
    if bond_length is None:
        bond_length = default_hydrogen_bond_length(outer)
    steps = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) * a / 2
    bonds = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * a / 4
    outer_sites = {}
    for i in range(size):
        for j in range(size - i):
            for k in range(size - i - j):
                outer_sites[(i, j, k)] = i * steps[0] + j * steps[1] + k * steps[2]
    inner_sites = {key: position + bonds[0] for key, position in outer_sites.items() if sum(key) <= size - 2}
    species, coords, kinds, passivates = [], [], [], []
    for key, position in outer_sites.items():
        species.append(outer)
        coords.append(position)
        kinds.append(outer)
        passivates.append(None)
    for position in inner_sites.values():
        species.append(inner)
        coords.append(position)
        kinds.append(inner)
        passivates.append(None)
    inner_positions = np.array(list(inner_sites.values()))
    for position in outer_sites.values():
        for bond in bonds:
            target = position + bond
            if len(inner_positions) and np.min(np.linalg.norm(inner_positions - target, axis=1)) < 1e-6:
                continue
            species.append("H")
            coords.append(position + bond / np.linalg.norm(bond) * bond_length)
            kinds.append(hydrogen.kind_name)
            passivates.append(outer)
    coords = np.array(coords)
    n_hydrogen = sum(1 for p in passivates if p)
    expected = counts["face"] + counts["edge"] + counts["corner"]
    if (len(outer_sites), len(inner_sites), n_hydrogen) != (counts["outer"], counts["inner"], expected):
        raise RuntimeError("cluster construction does not match the counts of Eq. 9")
    span = coords.max(axis=0) - coords.min(axis=0)
    edge = float(span.max() + vacuum)
    coords = coords - coords.mean(axis=0) + edge / 2
    return Structure(
        Lattice.cubic(edge), species, coords, coords_are_cartesian=True,
        site_properties={"kind_name": kinds, "pseudo_hydrogen": passivates, "bottom": [False] * len(species)},
    )


@dataclass(frozen=True)
class ClusterFit:
    """Result of fitting Eq. 9 to tetrahedral-cluster energies.

    ``face_ev``, ``edge_ev`` and ``corner_ev`` are the pseudo chemical
    potentials of face, edge and corner pseudo-hydrogen at ``mu_outer_ev``;
    the face value is the one for (111) surfaces. ``bulk_energy_ev`` is the
    fitted energy per formula unit and ``residual_ev`` the largest misfit.
    """

    outer: str
    mu_outer_ev: float
    bulk_energy_ev: float
    face_ev: float
    edge_ev: float
    corner_ev: float
    residual_ev: float
    sizes: tuple[int, ...]

    @property
    def reference(self) -> PseudoHydrogenReference:
        return PseudoHydrogenReference(
            self.outer, self.face_ev + self.mu_outer_ev / TETRAHEDRAL_COORDINATION, method="cluster",
        )


def fit_cluster_pseudo_chemical_potentials(
    outer: str, energies_ev: Mapping[int, float], mu_outer_ev: float,
) -> ClusterFit:
    """Solve Eq. 9 of the Sci. Rep. paper for E_AB and the face, edge and corner muhat.

    Args:
        outer: Element at the cluster surface (the passivated one).
        energies_ev: Total energy of each cluster, by size (at least four
            sizes, e.g. 2, 3, 8, 9 as in the paper).
        mu_outer_ev: Chemical potential of ``outer`` used in Eq. 9 (any
            value: the returned reference holds for every mu).
    """

    import numpy as np

    sizes = tuple(sorted(int(size) for size in energies_ev))
    if len(sizes) < 4:
        raise ValueError("at least four cluster sizes are needed for the four unknowns of Eq. 9")
    rows, right = [], []
    for size in sizes:
        counts = cluster_counts(size)
        rows.append([counts["inner"], counts["face"], counts["edge"], counts["corner"]])
        right.append(energies_ev[size] - (counts["outer"] - counts["inner"]) * mu_outer_ev)
    matrix, right = np.array(rows, dtype=float), np.array(right)
    solution, *_ = np.linalg.lstsq(matrix, right, rcond=None)
    if np.linalg.matrix_rank(matrix) < 4:
        raise ValueError(f"cluster sizes {sizes} do not determine the four unknowns; add larger clusters")
    residual = float(np.max(np.abs(matrix @ solution - right)))
    bulk, face, edge, corner = (float(value) for value in solution)
    return ClusterFit(outer, float(mu_outer_ev), bulk, face, edge, corner, residual, sizes)


def doubly_passivated_slab(bulk: Any, miller_index, *, bilayers: int | None = None, **kwargs) -> Any:
    """The ideal slab with both faces passivated, for the Eq. 7 self-consistency check.

    Its energy gives n_A muhat_A + n_B muhat_B = E - n_A mu_A - n_B mu_B
    directly (Sci. Rep. Eq. 7); compare with the values from the
    pseudo-molecules or clusters.
    """

    import numpy as np
    from pymatgen.core import Element, Structure

    slabs = find_polar_terminations(bulk, miller_index, bilayers=bilayers, electron_counting=False, **kwargs)
    ideal = slabs[0]
    structure = ideal.structure
    neighbours = bulk.get_all_neighbors(4.0)
    bond = min(neighbour.nn_distance for found in neighbours for neighbour in found)
    normal = structure.lattice.matrix[2] / structure.lattice.c
    heights = structure.cart_coords @ normal
    atoms = [i for i, p in enumerate(structure.site_properties["pseudo_hydrogen"]) if p is None]
    top = max(heights[i] for i in atoms)
    hydrogens = pseudo_hydrogens(bulk, slabs.oxidation_states)
    species, coords = [site.specie.symbol for site in structure], [site.coords for site in structure]
    properties = {key: list(values) for key, values in structure.site_properties.items()}
    for i in atoms:
        if top - heights[i] > 0.3:
            continue
        found = [n for n in structure.get_neighbors(structure[i], bond * 1.15)]
        missing = TETRAHEDRAL_COORDINATION - len(found)
        if missing != 1:
            raise ValueError("double passivation supports tops with one broken bond per atom")
        symbol = structure[i].specie.symbol
        species.append("H")
        coords.append(structure[i].coords + normal * default_hydrogen_bond_length(symbol))
        properties["kind_name"].append(hydrogens[symbol].kind_name)
        properties["pseudo_hydrogen"].append(symbol)
        properties["bottom"].append(False)
    return Structure(structure.lattice, species, coords, coords_are_cartesian=True, site_properties=properties)


# =============================================================================
# BOTTOM CHECK AFTER RELAXATION
# =============================================================================

#: Default tolerances of the post-relaxation bottom check (A).
BOTTOM_RMSD_TOLERANCE = 0.02
BOTTOM_MAX_TOLERANCE = 0.05


@dataclass(frozen=True)
class BottomCheck:
    """How the relaxed bottom of one slab compares with the reference slab.

    ``rmsd`` and ``max_deviation`` compare the relaxed bottom (pseudo-H
    included, rigid shift removed) with the relaxed bottom of the reference
    slab; ``relaxation`` is the RMS displacement of the bottom from its
    built positions, for information.
    """

    label: str
    rmsd: float
    max_deviation: float
    relaxation: float
    cell_unchanged: bool
    passed: bool


class BottomCheckReport(list):
    """Results of :func:`check_bottoms`, one :class:`BottomCheck` per slab."""

    def __init__(self, checks=(), *, reference: str = "", rmsd_tolerance: float = 0.0, max_tolerance: float = 0.0):
        super().__init__(checks)
        self.reference = reference
        self.rmsd_tolerance = rmsd_tolerance
        self.max_tolerance = max_tolerance

    @property
    def passed(self) -> list[str]:
        return [check.label for check in self if check.passed]

    @property
    def failed(self) -> list[str]:
        return [check.label for check in self if not check.passed]

    @property
    def all_passed(self) -> bool:
        return not self.failed

    def summary(self) -> str:
        columns = ("label", "RMSD (Å)", "max (Å)", "relaxed by (Å)", "cell", "result")
        rows = [
            [c.label, f"{c.rmsd:.3f}", f"{c.max_deviation:.3f}", f"{c.relaxation:.3f}",
             "same" if c.cell_unchanged else "changed", "ok" if c.passed else "FAILED"]
            for c in self
        ]
        widths = [max(len(cell) for cell in column) for column in zip(columns, *rows)]

        def line(cells):
            return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip()

        text = [
            f"Bottom after relaxation, compared with {self.reference} "
            f"(tolerance: RMSD {self.rmsd_tolerance} Å, max {self.max_tolerance} Å)",
            "", line(columns), line(["-" * w for w in widths]),
        ]
        text += [line(row) for row in rows]
        text.append("")
        text.append("All bottoms agree." if self.all_passed else
                    f"Left out of comparisons: {', '.join(self.failed)}.")
        return "\n".join(text)

    __str__ = summary
    __repr__ = summary


def _as_pymatgen(structure: Any) -> Any:
    if hasattr(structure, "get_pymatgen_structure"):
        return structure.get_pymatgen_structure()
    return structure


def _bottom_displacements(termination: PolarTermination, relaxed: Any):
    """Displacement of each bottom site from its built position (nearest same-element site)."""

    import numpy as np

    built = termination.structure
    lattice = relaxed.lattice
    displacements = []
    taken: set[int] = set()
    symbols = [site.specie.symbol for site in relaxed]
    frac_relaxed = relaxed.frac_coords
    for index in termination.bottom_indices:
        symbol = built[index].specie.symbol
        frac = lattice.get_fractional_coords(built[index].coords)
        delta = frac_relaxed - frac
        delta[:, :2] -= np.round(delta[:, :2])
        distances = np.linalg.norm(lattice.get_cartesian_coords(delta), axis=1)
        order = [j for j in np.argsort(distances) if symbols[j] == symbol and j not in taken]
        if not order:
            raise ValueError(f"relaxed structure of {termination.label} lacks a {symbol} site for the bottom")
        taken.add(int(order[0]))
        displacements.append(lattice.get_cartesian_coords(delta[order[0]]))
    return np.array(displacements)


def check_bottoms(
    terminations,
    relaxed_structures: Mapping[str, Any],
    *,
    reference: str | None = None,
    rmsd_tolerance: float = BOTTOM_RMSD_TOLERANCE,
    max_tolerance: float = BOTTOM_MAX_TOLERANCE,
) -> BottomCheckReport:
    """Check that the relaxed bottom is the same in every slab of one face.

    Args:
        terminations: The slabs as built (a :class:`PolarTerminationSet`),
            all on one bottom.
        relaxed_structures: Relaxed structure of each slab by label
            (pymatgen ``Structure`` or AiiDA ``StructureData``). Atoms are
            matched by element and position, so the order may differ.
        reference: Label of the slab whose relaxed bottom is the reference;
            default the first one.
        rmsd_tolerance, max_tolerance: Largest RMS and largest single
            deviation (A) from the reference bottom, after removing a rigid
            shift of each slab.

    Returns:
        A :class:`BottomCheckReport`; ``failed`` lists the slabs to leave out.
    """

    import numpy as np

    terminations = list(terminations)
    fingerprints = {termination.bottom_fingerprint for termination in terminations}
    if len(fingerprints) != 1:
        raise ValueError("the slabs were not built on one bottom; check each face separately")
    missing = [t.label for t in terminations if t.label not in relaxed_structures]
    if missing:
        raise ValueError(f"no relaxed structure for {missing}")
    reference = reference or terminations[0].label
    relative = {}
    relaxation = {}
    cell_unchanged = {}
    for termination in terminations:
        relaxed = _as_pymatgen(relaxed_structures[termination.label])
        cell_unchanged[termination.label] = bool(np.allclose(
            relaxed.lattice.matrix[:2], termination.structure.lattice.matrix[:2], atol=1e-3))
        displacements = _bottom_displacements(termination, relaxed)
        shift = displacements.mean(axis=0)
        relative[termination.label] = displacements - shift
        relaxation[termination.label] = float(np.sqrt(np.mean(np.sum(displacements ** 2, axis=1))))
    if reference not in relative:
        raise ValueError(f"unknown reference slab {reference!r}")
    checks = []
    for termination in terminations:
        difference = relative[termination.label] - relative[reference]
        norms = np.linalg.norm(difference, axis=1)
        rmsd = float(np.sqrt(np.mean(norms ** 2)))
        largest = float(norms.max())
        ok = cell_unchanged[termination.label] and rmsd <= rmsd_tolerance and largest <= max_tolerance
        checks.append(BottomCheck(termination.label, rmsd, largest, relaxation[termination.label],
                                  cell_unchanged[termination.label], ok))
    return BottomCheckReport(checks, reference=reference, rmsd_tolerance=rmsd_tolerance, max_tolerance=max_tolerance)


def polar_slab_terminations(
    terminations,
    energies_ev: Mapping[str, float],
    relaxed_structures: Mapping[str, Any],
    **check_options,
):
    """``SlabTermination`` objects for the slabs whose relaxed bottom passes :func:`check_bottoms`.

    Returns ``(slab_terminations, report)``; slabs that fail the check are
    left out, and the report says why.
    """

    import warnings

    from psteros.phase_diagram import SlabTermination

    report = check_bottoms(terminations, relaxed_structures, **check_options)
    kept = [
        SlabTermination.from_polar(termination, energies_ev[termination.label])
        for termination in terminations if termination.label in report.passed
    ]
    if report.failed:
        warnings.warn(
            f"bottom changed after relaxation in {', '.join(report.failed)}; left out of the phase diagram",
            stacklevel=2,
        )
    return kept, report
