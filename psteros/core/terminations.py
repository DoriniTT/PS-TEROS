"""
Charge-Neutral Slab Terminations

This module finds slab terminations suitable for semiconductors and
insulators. Pymatgen's ``SlabGenerator`` returns every distinct cut, and
``symmetrize=True`` makes the two faces equivalent by deleting sites without
regard to valence. For an ionic semiconductor most of those slabs carry a net
formal charge, so DFT places the surplus electrons (or holes) in the
conduction (or valence) band and the reference surface becomes a metal.

The criterion used here is the electron-counting rule in its simplest form:
with nominal oxidation states, a usable slab has zero net formal charge and
two symmetry-equivalent faces. Stoichiometry is not required. A neutral,
off-stoichiometric slab (for example an extra Ag2O unit) is closed-shell and
is handled by the chemical-potential dependent surface thermodynamics.

Algorithm:
    1. Assign oxidation states to the bulk and check that it is neutral.
       Stop if no bulk operation reverses the (hkl) normal (polar direction).
    2. Build a thick oriented slab and split it into units: single atoms, or
       finite clusters joined by ``unit_bonds`` (e.g. P-O for PO4).
       Units are never split.
    3. Group units into planes along the surface normal and take every
       contiguous window of planes.
    4. Keep windows with a symmetry operation that maps one face onto the
       other.
    5. If a symmetric window has formal charge Q != 0, remove surface units
       in symmetry-related pairs until Q = 0: whole outer planes first, then
       part of the next plane, with the fewest units per face.
    6. Keep slabs whose thickness lies in
       ``[min_slab_thickness, min_slab_thickness + d)``, d being the bulk
       repeat along the normal, and drop slabs whose environment near the
       face matches one already kept.

Usage:
    >>> terminations = find_charge_neutral_terminations(
    ...     bulk, (1, 1, 0), 10.0, oxidation_states={'Ag': 1, 'P': 5, 'O': -2},
    ...     unit_bonds={('P', 'O'): 1.9})
    >>> print(terminations)          # summary table
    >>> terminations.write('slabs/')  # POSCAR files + terminations.json
    >>> terminations.plot('slabs.png')

    or from a shell: ``psteros-terminations bulk.cif 110 --keep P-O:1.9``.
    See docs/CHARGE_NEUTRAL_TERMINATIONS.md.

Formal charge is a necessary, not a sufficient, condition: the band gap of
each candidate should still be confirmed with a single-point calculation.
Covalent semiconductors (Si, GaAs) need dangling-bond counting and surface
reconstructions, which are outside the scope of this module.
"""

from __future__ import annotations

import argparse
import itertools
import math
import re
import sys
import typing as t
import warnings
from dataclasses import dataclass, replace

import numpy as np
from pymatgen.core import Composition, Element, Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

CHARGE_TOLERANCE = 1e-6


class NoChargeNeutralTerminationError(ValueError):
    """Raised when no symmetric, charge-neutral termination can be built."""


@dataclass(frozen=True)
class Termination:
    """A symmetric, charge-neutral slab and how it was obtained.

    Attributes:
        structure: Slab with the c axis normal to the surface and the slab
            centred in the cell.
        formal_charge: Net formal charge of the slab (zero within tolerance).
        parent_formula: Formula of the symmetric whole-plane window it was
            built from.
        parent_formal_charge: Formal charge of that window before repair.
        removed_per_face: Formulas of the units removed from each face.
        thickness: Distance between the outermost nuclei, in Angstrom.
        is_stoichiometric: Whether the slab has the bulk composition ratio.
        label: Name used in workflows and file names (``term_0``, ...).
        miller_index: Surface orientation.
        supercell: In-plane repetition of the surface cell.
        removed_sites: Species and Cartesian position, in the frame of
            ``structure``, of every removed atom (both faces). Used to draw
            what the repair took away.
    """

    structure: Structure
    formal_charge: float
    parent_formula: str
    parent_formal_charge: float
    removed_per_face: tuple[str, ...]
    thickness: float
    is_stoichiometric: bool
    label: str = ''
    miller_index: tuple[int, int, int] = (0, 0, 0)
    supercell: tuple[int, int] = (1, 1)
    removed_sites: tuple[tuple[str, tuple[float, float, float]], ...] = ()

    @property
    def formula(self) -> str:
        return self.structure.composition.formula.replace(' ', '')

    @property
    def area(self) -> float:
        """Surface area of one face, in A^2."""
        matrix = self.structure.lattice.matrix
        return float(np.linalg.norm(np.cross(matrix[0], matrix[1])))

    @property
    def origin(self) -> str:
        """How the slab was obtained, in words."""
        if not self.removed_per_face:
            return 'neutral as cut'
        removed = ' + '.join(_count_formulas(self.removed_per_face))
        return f"{self.parent_formula} ({self.parent_formal_charge:+g}) minus {removed} per face"

    def to_dict(self) -> dict:
        return {
            'label': self.label,
            'formula': self.formula,
            'n_atoms': len(self.structure),
            'miller_index': list(self.miller_index),
            'supercell': list(self.supercell),
            'thickness': self.thickness,
            'area': self.area,
            'formal_charge': self.formal_charge,
            'is_stoichiometric': self.is_stoichiometric,
            'parent_formula': self.parent_formula,
            'parent_formal_charge': self.parent_formal_charge,
            'removed_per_face': list(self.removed_per_face),
        }


def hkl_label(miller_index: t.Sequence[int]) -> str:
    """(1, -1, 0) -> '(1-10)'; indices above 9 are comma-separated."""
    values = [int(value) for value in miller_index]
    separator = ',' if any(abs(value) > 9 for value in values) else ''
    return '(' + separator.join(str(value) for value in values) + ')'


def _count_formulas(formulas: t.Sequence[str]) -> list[str]:
    """('Ag1', 'Ag1', 'P1O4') -> ['2 Ag', 'PO4']."""
    counts: dict[str, int] = {}
    for formula in formulas:
        name = re.sub(r'(?<=[A-Za-z])1(?![0-9])', '', formula.replace(' ', ''))
        counts[name] = counts.get(name, 0) + 1
    return [f"{count} {name}" if count > 1 else name for name, count in counts.items()]


class TerminationSet(list):
    """The terminations of one surface, with a readable summary.

    Behaves like a list of :class:`Termination`. Printing it (or displaying
    it in Jupyter) shows a table; :meth:`write` saves the slabs and
    :meth:`plot` draws side views.
    """

    def __init__(
        self,
        terminations: t.Iterable[Termination] = (),
        *,
        bulk_formula: str = '',
        miller_index: tuple[int, int, int] = (0, 0, 0),
        supercell: tuple[int, int] = (1, 1),
        oxidation_states: t.Mapping[str, float] | None = None,
        unit_bonds: t.Mapping[tuple[str, str], float] | None = None,
        repeat: float = 0.0,
    ):
        super().__init__(terminations)
        self.bulk_formula = bulk_formula
        self.miller_index = tuple(miller_index)
        self.supercell = tuple(supercell)
        self.oxidation_states = dict(oxidation_states or {})
        self.unit_bonds = dict(unit_bonds or {})
        self.repeat = repeat

    @property
    def is_metallic(self) -> bool:
        """True when every formal charge is zero (metals, alloys, Si...)."""
        return all(value == 0 for value in self.oxidation_states.values())

    def header(self) -> str:
        cell = 'x'.join(str(value) for value in self.supercell)
        parts = [f"{self.bulk_formula}{hkl_label(self.miller_index)}", f"{cell} surface cell"]
        if self.is_metallic:
            parts.append('formal charges all zero: only face symmetry is enforced')
        else:
            ordered = sorted(self.oxidation_states, key=lambda element: Element(element).X)
            parts.append(' '.join(f"{element}{self.oxidation_states[element]:+g}" for element in ordered))
        if self.unit_bonds:
            pairs = (sorted(pair, key=lambda symbol: Element(symbol).X) for pair in self.unit_bonds)
            parts.append('kept whole: ' + ', '.join(f"{a}-{b}" for a, b in pairs))
        return ' | '.join(parts)

    def rows(self) -> list[list[str]]:
        return [
            [
                termination.label,
                termination.formula,
                str(len(termination.structure)),
                f"{termination.thickness:.2f}",
                'yes' if termination.is_stoichiometric else 'no',
                termination.origin,
            ]
            for termination in self
        ]

    COLUMNS = ('label', 'formula', 'atoms', 'thickness (Å)', 'stoich.', 'origin')

    def summary(self) -> str:
        """A plain-text table of the terminations."""
        rows = self.rows()
        widths = [max(len(cell) for cell in column) for column in zip(self.COLUMNS, *rows)]
        line = lambda cells: '  '.join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip()  # noqa: E731
        text = [self.header(), '', line(self.COLUMNS), line(['-' * width for width in widths])]
        text += [line(row) for row in rows]
        noun = 'termination' if len(self) == 1 else 'terminations'
        text += ['', f"{len(self)} {noun}, each symmetric with zero formal charge. "
                     "Confirm the band gap with a single-point calculation."]
        return '\n'.join(text)

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        return self.summary()

    def _repr_html_(self) -> str:
        import html

        head = ''.join(f"<th style='text-align:left'>{html.escape(column)}</th>" for column in self.COLUMNS)
        body = ''.join(
            '<tr>' + ''.join(f"<td>{html.escape(cell)}</td>" for cell in row) + '</tr>'
            for row in self.rows()
        )
        return (
            f"<p><b>{html.escape(self.header())}</b></p>"
            f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
            "<p><i>Each slab is symmetric with zero formal charge. "
            "Confirm the band gap with a single-point calculation.</i></p>"
        )

    def to_dicts(self) -> list[dict]:
        return [termination.to_dict() for termination in self]

    def write(self, directory: str, fmt: str = 'poscar') -> list[str]:
        """Write every slab and a ``terminations.json`` summary.

        Args:
            directory: Output directory (created if needed).
            fmt: ``'poscar'`` (``term_0_Ag18P6O24.vasp``), ``'cif'`` or
                ``'extxyz'``.

        Returns:
            The paths written.
        """
        import json
        import os

        extensions = {'poscar': 'vasp', 'cif': 'cif', 'extxyz': 'extxyz'}
        if fmt not in extensions:
            raise ValueError(f"fmt must be one of {sorted(extensions)}, got {fmt!r}")
        os.makedirs(directory, exist_ok=True)
        paths = []
        for termination in self:
            path = os.path.join(directory, f"{termination.label}_{termination.formula}.{extensions[fmt]}")
            if fmt == 'extxyz':
                from pymatgen.io.ase import AseAtomsAdaptor
                from ase.io import write as ase_write
                ase_write(path, AseAtomsAdaptor.get_atoms(termination.structure), format='extxyz')
            else:
                termination.structure.to(filename=path, fmt=fmt)
            paths.append(path)
        summary = os.path.join(directory, 'terminations.json')
        with open(summary, 'w') as handle:
            json.dump({'header': self.header(), 'terminations': self.to_dicts()}, handle, indent=2)
        paths.append(summary)
        return paths

    def plot(self, filename: str | None = None, repeat: int = 2, **kwargs):
        """Side views of every termination; see :func:`plot_terminations`."""
        return plot_terminations(self, filename=filename, repeat=repeat, **kwargs)


@dataclass
class _Unit:
    """A finite group of sites that is kept or removed as a whole."""

    indices: tuple[int, ...]
    coords: np.ndarray  # unwrapped Cartesian coordinates, one row per site
    charge: float
    composition: Composition
    height: float = 0.0


# =============================================================================
# OXIDATION STATES AND FORMAL CHARGE
# =============================================================================

def resolve_oxidation_states(
    bulk: Structure,
    oxidation_states: t.Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Return one oxidation state per element and check the bulk is neutral.

    Args:
        bulk: Bulk crystal.
        oxidation_states: Element -> oxidation state, e.g.
            ``{'Ag': 1, 'P': 5, 'O': -2}``. If omitted, pymatgen's most
            probable guess for the bulk composition is used.

    Raises:
        ValueError: If an element is missing, no guess is found, or the bulk
            is not neutral with the chosen states.
    """
    elements = sorted({site.specie.symbol for site in bulk})
    if oxidation_states is None:
        guesses = bulk.composition.oxi_state_guesses()
        if guesses:
            oxidation_states = guesses[0]
        elif all(Element(element).is_metal for element in elements):
            # Intermetallics have no ionic model: every formal charge is zero
            # and only face symmetry constrains the slab.
            oxidation_states = {element: 0 for element in elements}
        else:
            raise ValueError(
                f"Could not guess oxidation states for {bulk.composition.reduced_formula}; "
                "pass oxidation_states explicitly, e.g. {'Ag': 1, 'P': 5, 'O': -2}."
            )
    states = {str(element): float(value) for element, value in oxidation_states.items()}
    missing = [element for element in elements if element not in states]
    if missing:
        raise ValueError(f"oxidation_states is missing elements: {missing}")
    charge = formal_charge(bulk, states)
    if abs(charge) > CHARGE_TOLERANCE:
        raise ValueError(
            f"Bulk {bulk.composition.reduced_formula} is not neutral with oxidation "
            f"states {states} (net charge {charge:+g} per cell)."
        )
    return {element: states[element] for element in elements}


def formal_charge(structure: Structure, oxidation_states: t.Mapping[str, float]) -> float:
    """Sum of nominal oxidation states over all sites."""
    return float(sum(oxidation_states[site.specie.symbol] for site in structure))


def is_stoichiometric(structure: Structure, bulk: Structure) -> bool:
    """Whether ``structure`` has the same composition ratios as ``bulk``."""
    return (
        structure.composition.reduced_composition
        == bulk.composition.reduced_composition
    )


# =============================================================================
# SYMMETRY BETWEEN THE TWO FACES
# =============================================================================

def _analyzer(structure: Structure, symprec: float) -> SpacegroupAnalyzer:
    """SpacegroupAnalyzer without spglib's deprecation chatter."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        return SpacegroupAnalyzer(structure, symprec=symprec)


def _quiet(call: t.Callable, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        return call(*args, **kwargs)


def _site_permutation(structure: Structure, operation, tolerance: float) -> np.ndarray | None:
    """Index of the site each site is mapped onto, or None if it fails."""
    frac = structure.frac_coords
    mapped = operation.operate_multi(frac)
    species = [site.specie.symbol for site in structure]
    lattice = structure.lattice
    permutation = np.full(len(structure), -1, dtype=int)
    for index, point in enumerate(mapped):
        delta = frac - point
        delta -= np.round(delta)
        distances = np.linalg.norm(lattice.get_cartesian_coords(delta), axis=1)
        candidates = [
            other for other in np.argsort(distances)[:4]
            if distances[other] < tolerance and species[other] == species[index]
        ]
        if not candidates:
            return None
        permutation[index] = candidates[0]
    if len(set(permutation.tolist())) != len(structure):
        return None
    return permutation


def _slab_symmetry(
    structure: Structure,
    symprec: float = 0.1,
) -> list[tuple[t.Any, np.ndarray, bool]]:
    """Symmetry operations of a slab whose c axis is normal to the surface.

    Returns:
        ``(operation, permutation, reverses_normal)`` for every operation that
        keeps the surface plane, where ``permutation[i]`` is the site that site
        ``i`` is mapped onto. Operations that tilt the normal are dropped.
    """
    try:
        analyzer = _analyzer(structure, symprec)
        operations = _quiet(analyzer.get_symmetry_operations)
    except Exception:  # spglib can fail on degenerate inputs
        return []
    result = []
    for operation in operations:
        rotation = np.round(operation.rotation_matrix).astype(int)
        if rotation[0, 2] != 0 or rotation[1, 2] != 0 or rotation[2, 0] != 0 or rotation[2, 1] != 0:
            continue
        permutation = _site_permutation(structure, operation, tolerance=max(2 * symprec, 0.05))
        if permutation is not None:
            result.append((operation, permutation, bool(rotation[2, 2] == -1)))
    return result


def find_face_operation(
    structure: Structure,
    symprec: float = 0.1,
) -> tuple[t.Any, np.ndarray] | None:
    """Find a symmetry operation that maps one face of a slab onto the other.

    The slab must have its c axis normal to the surface. Accepted operations
    reverse c (inversion, a mirror parallel to the surface, a two-fold axis or
    glide in the surface plane, ...).

    Returns:
        ``(operation, permutation)`` with the fractional-coordinate operation
        and the site each site is mapped onto, or None if the faces are not
        equivalent.
    """
    for operation, permutation, reverses in _slab_symmetry(structure, symprec):
        if reverses:
            return operation, permutation
    return None


def has_face_reversing_operation(
    bulk: Structure,
    miller_index: t.Sequence[int],
    symprec: float = 0.1,
) -> bool:
    """Whether a bulk point-group operation maps the (hkl) normal n onto -n.

    Without one, (hkl) is a polar direction and no slab of it has two
    equivalent faces (e.g. zinc blende (111)).
    """
    normal = bulk.lattice.reciprocal_lattice.get_cartesian_coords(miller_index)
    normal = normal / np.linalg.norm(normal)
    analyzer = _analyzer(bulk, symprec)
    operations = _quiet(analyzer.get_point_group_operations, cartesian=True)
    return any(
        np.allclose(operation.rotation_matrix @ normal, -normal, atol=1e-3)
        for operation in operations
    )


def classify_slab(
    slab: Structure,
    bulk: Structure,
    oxidation_states: t.Mapping[str, float] | None = None,
    symprec: float = 0.1,
) -> dict:
    """Report formal charge, face symmetry and stoichiometry of a slab.

    Useful for checking slabs from other sources (``input_slabs``) before
    spending DFT time on them. The slab's c axis must be normal to the
    surface.
    """
    states = resolve_oxidation_states(bulk, oxidation_states)
    charge = formal_charge(slab, states)
    return {
        'formula': slab.composition.formula.replace(' ', ''),
        'formal_charge': charge,
        'is_charge_neutral': abs(charge) <= CHARGE_TOLERANCE,
        'is_symmetric': find_face_operation(slab, symprec) is not None,
        'is_stoichiometric': is_stoichiometric(slab, bulk),
    }


# =============================================================================
# UNITS (ATOMS AND FINITE CLUSTERS)
# =============================================================================

def _normalise_bonds(unit_bonds) -> dict[tuple[str, str], float]:
    """Accept {(A, B): r}, {'A-B': r} or [[A, B, r], ...]."""
    if not unit_bonds:
        return {}
    items = unit_bonds.items() if isinstance(unit_bonds, t.Mapping) else (
        ((entry[0], entry[1]), entry[2]) for entry in unit_bonds
    )
    bonds = {}
    for key, cutoff in items:
        first, second = key.split('-') if isinstance(key, str) else key
        bonds[tuple(sorted((str(first), str(second))))] = float(cutoff)
    return bonds


def _find_units(
    structure: Structure,
    bonds: dict[tuple[str, str], float],
    oxidation_states: t.Mapping[str, float],
) -> list[_Unit]:
    """Split a structure into finite bonded clusters, unwrapping images.

    Raises:
        ValueError: If the bonds form an extended network (a site reaches
            itself through a lattice translation).
    """
    neighbours: list[list[tuple[int, np.ndarray]]] = [[] for _ in structure]
    if bonds:
        cutoff = max(bonds.values())
        centres, points, images, _ = structure.get_neighbor_list(cutoff)
        for centre, point, image in zip(centres, points, images):
            pair = tuple(sorted((structure[centre].specie.symbol, structure[point].specie.symbol)))
            limit = bonds.get(pair)
            if limit is None:
                continue
            if structure[centre].distance(structure[point], jimage=image) <= limit:
                neighbours[centre].append((int(point), np.asarray(image, dtype=int)))

    frac = structure.frac_coords
    seen: dict[int, np.ndarray] = {}
    units = []
    for start in range(len(structure)):
        if start in seen:
            continue
        seen[start] = np.zeros(3, dtype=int)
        queue, members = [start], [start]
        while queue:
            current = queue.pop()
            for other, image in neighbours[current]:
                offset = seen[current] + image
                if other in seen:
                    if not np.array_equal(seen[other], offset):
                        raise ValueError(
                            "unit_bonds connect sites into an extended network around "
                            f"{structure[other].specie.symbol}; only finite units such as "
                            "PO4 or SO4 can be kept whole."
                        )
                    continue
                seen[other] = offset
                members.append(other)
                queue.append(other)
        members.sort()
        coords = structure.lattice.get_cartesian_coords(
            np.array([frac[index] + seen[index] for index in members])
        )
        composition = Composition(
            ' '.join(structure[index].specie.symbol for index in members)
        )
        units.append(_Unit(
            indices=tuple(members),
            coords=coords,
            charge=float(sum(oxidation_states[structure[index].specie.symbol] for index in members)),
            composition=composition,
        ))
    return units


# =============================================================================
# SLAB CONSTRUCTION
# =============================================================================

def _assemble(
    units: t.Sequence[_Unit],
    species_of: t.Callable[[int], str],
    in_plane: np.ndarray,
    normal: np.ndarray,
    vacuum: float,
    ghosts: t.Sequence[_Unit] = (),
) -> tuple[Structure, float, np.ndarray, tuple]:
    """Assemble units into a slab with c along ``normal``, centred in the cell.

    Returns the slab, its thickness, for every site the position of its unit
    in ``units``, and the species and Cartesian positions of the ``ghosts``
    (removed units) in the same frame.
    """
    coords = np.vstack([unit.coords for unit in units])
    species = [species_of(index) for unit in units for index in unit.indices]
    site_unit = np.repeat(np.arange(len(units)), [len(unit.indices) for unit in units])
    heights = coords @ normal
    thickness = float(heights.max() - heights.min())
    length = thickness + vacuum
    lattice = Lattice(np.vstack([in_plane, normal * length]))
    shift = normal * (0.5 * length - 0.5 * (heights.max() + heights.min()))
    slab = Structure(lattice, species, coords + shift, coords_are_cartesian=True, to_unit_cell=True)
    ghost_sites = []
    for unit in ghosts:
        for index, position in zip(unit.indices, unit.coords + shift):
            frac = lattice.get_fractional_coords(position)
            frac[:2] -= np.floor(frac[:2])
            ghost_sites.append((species_of(index), tuple(float(x) for x in lattice.get_cartesian_coords(frac))))
    return slab, thickness, site_unit, tuple(ghost_sites)


def normal_repeat(bulk: Structure, miller_index: t.Sequence[int], symprec: float = 0.1) -> float:
    """Shortest bulk translation along the (hkl) normal, in Angstrom.

    Includes centring translations, so rock-salt (100) gives a/2, not a.
    """
    miller = np.array(miller_index, dtype=int)
    miller = miller // (math.gcd(*(abs(int(value)) for value in miller)) or 1)
    vector = bulk.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(miller)
    spacing = 1.0 / float(np.linalg.norm(vector))
    normal = vector * spacing
    repeat = spacing
    analyzer = _analyzer(bulk, symprec)
    for operation in _quiet(analyzer.get_symmetry_operations, cartesian=True):
        if not np.allclose(operation.rotation_matrix, np.eye(3), atol=1e-3):
            continue
        shift = float(operation.translation_vector @ normal) % spacing
        if 1e-3 < shift < spacing - 1e-3:
            repeat = min(repeat, shift)
    return repeat


def _oriented_parent(
    bulk: Structure,
    miller_index: tuple[int, int, int],
    thickness: float,
    vacuum: float,
    lll_reduce: bool,
    max_normal_search: int | None,
    supercell: tuple[int, int],
    repeat: float,
    symprec: float,
) -> Structure:
    """A thick slab with c normal to the surface and a primitive surface cell."""
    generator = SlabGenerator(
        bulk,
        miller_index,
        min_slab_size=thickness + 2.0 * repeat,
        min_vacuum_size=vacuum,
        center_slab=True,
        in_unit_planes=False,
        primitive=True,
        max_normal_search=max_normal_search,
        lll_reduce=lll_reduce,
    )
    parent = Structure.from_sites(generator.get_slab(shift=0.0).get_orthogonal_c_slab())
    # SlabGenerator may cut through an atomic plane, leaving partial planes at
    # the edges, and its surface cell can be a supercell of the primitive one
    # (4x for cubic (111)). Trimming by height keeps every in-plane symmetry,
    # after which the primitive surface cell can be found with c held normal.
    normal = parent.lattice.matrix[2] / np.linalg.norm(parent.lattice.matrix[2])
    heights = parent.cart_coords @ normal
    levels = np.sort(heights)
    gaps = [(low + high) / 2 for low, high in zip(levels, levels[1:]) if high - low > 0.05]
    # Cut in the gaps between atomic planes, never through one.
    low_cut = min(gap for gap in gaps if gap >= levels[0] + repeat)
    high_cut = max(gap for gap in gaps if gap <= levels[-1] - repeat)
    outside = np.where((heights < low_cut) | (heights > high_cut))[0]
    parent.remove_sites(outside.tolist())
    reduced = parent.get_primitive_structure(
        tolerance=symprec,
        constrain_latt={'c': parent.lattice.c, 'alpha': 90, 'beta': 90},
    )
    if len(reduced) < len(parent):
        parent = Structure.from_sites(reduced)
    if tuple(supercell) != (1, 1):
        parent.make_supercell([supercell[0], supercell[1], 1])
    return parent


def _group_planes(units: t.Sequence[_Unit], tolerance: float) -> list[list[_Unit]]:
    """Group units sorted by height into planes."""
    planes: list[list[_Unit]] = []
    for unit in units:
        if planes and unit.height - planes[-1][-1].height < tolerance:
            planes[-1].append(unit)
        else:
            planes.append([unit])
    return planes


def find_charge_neutral_terminations(
    bulk: Structure,
    miller_index: t.Sequence[int],
    min_slab_thickness: float,
    min_vacuum_thickness: float = 15.0,
    *,
    oxidation_states: t.Mapping[str, float] | None = None,
    unit_bonds=None,
    supercell: tuple[int, int] | None = None,
    surface_depth: float | None = None,
    max_removed_per_face: int | None = None,
    max_combinations: int = 20000,
    max_variants_per_window: int = 20,
    plane_tolerance: float = 0.3,
    symprec: float = 0.1,
    lll_reduce: bool = True,
    max_normal_search: int | None = None,
) -> TerminationSet:
    """Find the symmetric, charge-neutral terminations of a surface.

    Args:
        bulk: Bulk crystal.
        miller_index: Surface orientation, e.g. ``(1, 1, 0)``.
        min_slab_thickness: Minimum distance between outermost nuclei (A).
            Slabs are returned within one bulk repeat above this value, so
            each termination appears once.
        min_vacuum_thickness: Vacuum added on top of the slab thickness (A).
        oxidation_states: Element -> nominal oxidation state. Guessed from the
            bulk composition when omitted.
        unit_bonds: Bonds that must never be broken, so their clusters are
            kept or removed as whole units, e.g. ``{('P', 'O'): 1.9}``,
            ``{'P-O': 1.9}`` or ``[['P', 'O', 1.9]]``. Cation-anion bonds of
            the framework itself must not be listed (they form an extended
            network).
        supercell: In-plane repetition (n1, n2) of the surface cell. Larger
            cells allow partial surface coverages, e.g. removing one cation in
            two. By default 1x1, 2x1, 1x2 and 2x2 are tried in turn and the
            first cell with a solution is used.
        surface_depth: How far below each face (A) the partially emptied plane
            may lie. Defaults to half the bulk repeat along the normal.
        max_removed_per_face: Largest number of units removed per face.
            Defaults to 4 per in-plane cell of the supercell.
        max_combinations: Skip a plane whose number of removal subsets of one
            size exceeds this.
        max_variants_per_window: Largest number of symmetry-distinct removal
            patterns kept per window. Large supercells with partial coverages
            can have thousands; a warning is issued when the list is cut.
        plane_tolerance: Units whose heights differ by less than this (A)
            belong to the same plane.
        symprec: Symmetry tolerance (A).
        lll_reduce: Passed to ``SlabGenerator``.
        max_normal_search: Passed to ``SlabGenerator``.

    Returns:
        A :class:`TerminationSet` (a list of :class:`Termination`), sorted by
        number of removed units and then by thickness, labelled ``term_0``,
        ``term_1``, ... Print it for a summary table. Each slab has c normal
        to the surface.

    Raises:
        NoChargeNeutralTerminationError: If (hkl) is polar or no termination
            is found.
    """
    with warnings.catch_warnings():
        # spglib (called by pymatgen) warns about its error-handling API on
        # every call; that is noise for the user.
        warnings.filterwarnings('ignore', category=DeprecationWarning, module='spglib')
        return _find_terminations(
            bulk, miller_index, min_slab_thickness, min_vacuum_thickness, oxidation_states,
            unit_bonds, supercell, surface_depth, max_removed_per_face, max_combinations,
            max_variants_per_window, plane_tolerance, symprec, lll_reduce, max_normal_search,
        )


def _find_terminations(
    bulk, miller_index, min_slab_thickness, min_vacuum_thickness, oxidation_states,
    unit_bonds, supercell, surface_depth, max_removed_per_face, max_combinations,
    max_variants_per_window, plane_tolerance, symprec, lll_reduce, max_normal_search,
) -> TerminationSet:
    miller_index = tuple(int(value) for value in miller_index)
    states = resolve_oxidation_states(bulk, oxidation_states)
    bonds = _normalise_bonds(unit_bonds)
    formula = bulk.composition.reduced_formula

    if not has_face_reversing_operation(bulk, miller_index, symprec):
        raise NoChargeNeutralTerminationError(
            f"{formula}{hkl_label(miller_index)} is a polar direction: no bulk symmetry operation "
            "reverses the surface normal, so the two faces of a slab can never be "
            "equivalent. Use an asymmetric slab with a dipole correction or a "
            "passivated back face instead."
        )

    repeat = normal_repeat(bulk, miller_index, symprec)
    cells = AUTO_SUPERCELLS if supercell is None else [(int(supercell[0]), int(supercell[1]))]
    failures = []
    for cell in cells:
        try:
            found = _search(
                bulk, miller_index, min_slab_thickness, min_vacuum_thickness, states, bonds,
                cell, repeat, surface_depth,
                4 * cell[0] * cell[1] if max_removed_per_face is None else max_removed_per_face,
                max_combinations, max_variants_per_window, plane_tolerance, symprec,
                lll_reduce, max_normal_search,
            )
        except NoChargeNeutralTerminationError as error:
            failures.append((cell, error))
            continue
        labelled = [
            replace(termination, label=f"term_{index}", miller_index=miller_index, supercell=cell)
            for index, termination in enumerate(found)
        ]
        return TerminationSet(
            labelled, bulk_formula=formula, miller_index=miller_index, supercell=cell,
            oxidation_states=states, unit_bonds=bonds, repeat=repeat,
        )

    tried = ', '.join('x'.join(map(str, cell)) for cell, _ in failures)
    details = '\n'.join(f"  {'x'.join(map(str, cell))}: {error}" for cell, error in failures)
    raise NoChargeNeutralTerminationError(
        f"No symmetric, charge-neutral termination of {formula}{hkl_label(miller_index)} "
        f"in the surface cells tried ({tried}).\n{details}\n"
        "  Try a larger supercell, surface_depth or max_removed_per_face, or check "
        "oxidation_states and unit_bonds."
    )


#: Surface cells tried, in order, when ``supercell`` is not given.
AUTO_SUPERCELLS = [(1, 1), (2, 1), (1, 2), (2, 2)]


def _search(
    bulk: Structure,
    miller_index: tuple[int, int, int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    states: dict[str, float],
    bonds: dict[tuple[str, str], float],
    supercell: tuple[int, int],
    repeat: float,
    surface_depth: float | None,
    max_removed_per_face: int,
    max_combinations: int,
    max_variants_per_window: int,
    plane_tolerance: float,
    symprec: float,
    lll_reduce: bool,
    max_normal_search: int | None,
) -> list[Termination]:
    """Search one surface cell; see :func:`find_charge_neutral_terminations`."""
    # Bulk units define what a complete unit looks like.
    unit_types = {unit.composition for unit in _find_units(bulk, bonds, states)}

    depth = 0.5 * repeat if surface_depth is None else float(surface_depth)
    # Removing units thins a window by up to 2 * depth, so windows are taken up
    # to that much thicker. Final slabs must lie within one repeat above the
    # minimum. The parent adds a few repeats of margin for broken edge units.
    lower, upper = min_slab_thickness, min_slab_thickness + repeat
    parent = _oriented_parent(
        bulk, miller_index, upper + 2.0 * depth + 4.0 * repeat, min_vacuum_thickness,
        lll_reduce, max_normal_search, supercell, repeat, symprec,
    )

    normal = parent.lattice.matrix[2] / np.linalg.norm(parent.lattice.matrix[2])
    in_plane = parent.lattice.matrix[:2]
    units = _find_units(parent, bonds, states)
    for unit in units:
        unit.height = float(unit.coords.mean(axis=0) @ normal)
    units.sort(key=lambda unit: unit.height)

    # Keep the longest run of planes made only of complete units; incomplete
    # units are the ones cut at the parent's edges.
    runs, current = [], []
    for plane in _group_planes(units, plane_tolerance):
        if all(unit.composition in unit_types for unit in plane):
            current.append(plane)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    planes = max(runs, key=len) if runs else []

    species_of = lambda index: parent[index].specie.symbol  # noqa: E731
    candidates: list[tuple[int, float, Termination]] = []
    windows_tried = symmetric_windows = 0
    charges_seen: set[float] = set()

    for start in range(len(planes)):
        for stop in range(start, len(planes)):
            window = [unit for plane in planes[start:stop + 1] for unit in plane]
            heights = np.concatenate([unit.coords @ normal for unit in window])
            thickness = float(heights.max() - heights.min())
            if thickness < lower:
                continue
            if thickness >= upper + 2.0 * depth:
                break
            windows_tried += 1
            slab, thickness, site_unit, _ = _assemble(
                window, species_of, in_plane, normal, min_vacuum_thickness,
            )
            symmetry = _slab_symmetry(slab, symprec)
            if not any(reverses for _, _, reverses in symmetry):
                continue
            symmetric_windows += 1
            charge = formal_charge(slab, states)
            charges_seen.add(round(charge, 6))
            parent_formula = slab.composition.formula.replace(' ', '')
            if abs(charge) <= CHARGE_TOLERANCE:
                if thickness < upper:
                    candidates.append((0, thickness, Termination(
                        slab.get_sorted_structure(), 0.0, parent_formula, 0.0, (),
                        thickness, is_stoichiometric(slab, bulk),
                    )))
                continue
            for removed in _neutralising_removals(
                window, charge, symmetry, site_unit, plane_tolerance, depth,
                max_removed_per_face, max_combinations, max_variants_per_window,
            ):
                kept = [unit for index, unit in enumerate(window) if index not in removed]
                gone = [window[index] for index in sorted(removed)]
                repaired, repaired_thickness, _, ghosts = _assemble(
                    kept, species_of, in_plane, normal, min_vacuum_thickness, ghosts=gone,
                )
                if not lower <= repaired_thickness < upper:
                    continue
                residual = formal_charge(repaired, states)
                if abs(residual) > CHARGE_TOLERANCE:
                    continue
                top_half = sorted(removed, key=lambda index: -window[index].height)[:len(removed) // 2]
                formulas = tuple(sorted(
                    window[index].composition.formula.replace(' ', '') for index in top_half
                ))
                candidates.append((len(top_half), repaired_thickness, Termination(
                    repaired.get_sorted_structure(), residual, parent_formula, charge,
                    formulas, repaired_thickness, is_stoichiometric(repaired, bulk),
                    removed_sites=ghosts,
                )))

    if not candidates:
        charges = ', '.join(f"{charge:+g}" for charge in sorted(charges_seen)) or 'none'
        raise NoChargeNeutralTerminationError(
            f"{windows_tried} slabs tried, {symmetric_windows} symmetric, formal charges "
            f"of the symmetric ones: {charges}; none could be neutralised."
        )

    # Remove duplicates. Two slabs are the same termination when the
    # environment within ``probe`` of the face matches; this also catches
    # identical slabs. Candidates are visited thinnest first, then with the
    # fewest removed units, so the thinner, simpler slab is kept.
    candidates.sort(key=lambda item: (round(item[1], 2), item[0]))
    probe = max(3.0, depth + 1.5)
    unique: list[tuple[int, dict, Termination]] = []
    for removed_count, _, termination in candidates:
        face = _distance_fingerprint(termination.structure, probe=probe)
        if not any(_same_fingerprint(face, kept_face) for _, kept_face, _ in unique):
            unique.append((removed_count, face, termination))
    unique.sort(key=lambda item: (item[0], item[2].thickness))
    return [termination for *_, termination in unique]


def _distance_fingerprint(
    structure: Structure,
    cutoff: float = 4.0,
    probe: float | None = None,
) -> dict:
    """Sorted heights and interatomic distances, unchanged by any isometry.

    With ``probe``, only sites within ``probe`` of the top face (and their
    distances to any neighbour) are included, which characterises the
    termination rather than the whole slab.
    """
    normal = structure.lattice.matrix[2] / np.linalg.norm(structure.lattice.matrix[2])
    heights = structure.cart_coords @ normal
    depths = heights.max() - heights
    species = [site.specie.symbol for site in structure]
    selected = np.ones(len(structure), dtype=bool) if probe is None else depths <= probe
    fingerprint: dict[tuple, list[float]] = {}
    for index in np.where(selected)[0]:
        fingerprint.setdefault(('depth', species[index]), []).append(depths[index])
    centres, points, _, distances = structure.get_neighbor_list(cutoff)
    for i, j, distance in zip(centres, points, distances):
        if selected[i] and (probe is not None or i < j):
            fingerprint.setdefault(('pair', species[i], species[j]), []).append(distance)
    return {key: np.sort(values) for key, values in fingerprint.items()}


def _same_fingerprint(first: dict, second: dict, tolerance: float = 0.05) -> bool:
    return first.keys() == second.keys() and all(
        len(first[key]) == len(second[key])
        and np.allclose(first[key], second[key], atol=tolerance)
        for key in first
    )


def _neutralising_removals(
    window: list[_Unit],
    charge: float,
    symmetry: list[tuple[t.Any, np.ndarray, bool]],
    site_unit: np.ndarray,
    plane_tolerance: float,
    depth: float,
    max_removed: int,
    max_combinations: int,
    max_variants: int,
) -> t.Iterator[frozenset[int]]:
    """Yield sets of window units whose removal makes the slab neutral.

    Units are removed from the top face layer by layer: all units of the
    outermost planes, then a subset of the next plane, so no vacancy is left
    below a kept surface unit. Each unit removed from the top face is paired
    with its image on the bottom face under a face-reversing operation, so the
    slab stays symmetric. Only the smallest number of units per face that
    reaches zero charge is used, and removal sets related by a symmetry of
    the window are yielded once.
    """
    # Unit permutations for every symmetry operation of the window.
    unit_permutations = []
    for _, permutation, reverses in symmetry:
        mapping = np.full(len(window), -1, dtype=int)
        for site, image in enumerate(permutation):
            unit, image_unit = site_unit[site], site_unit[image]
            if mapping[unit] not in (-1, image_unit):
                mapping = None
                break
            mapping[unit] = image_unit
        if mapping is not None and (mapping >= 0).all():
            unit_permutations.append((mapping, reverses))
    faces = []
    for mapping, reverses in unit_permutations:
        if reverses and not any(np.array_equal(mapping, other) for other in faces):
            faces.append(mapping)
    if not faces:
        return
    face = faces[0]

    heights = np.array([unit.height for unit in window])
    top_units = sorted(
        (index for index in range(len(window)) if heights[face[index]] < heights[index] - 1e-6),
        key=lambda index: -heights[index],
    )
    planes: list[list[int]] = []
    for index in top_units:
        if planes and heights[planes[-1][-1]] - heights[index] < plane_tolerance:
            planes[-1].append(index)
        else:
            planes.append([index])

    target = charge / 2.0
    solutions: dict[int, list[tuple[int, ...]]] = {}
    for depth_index, plane in enumerate(planes):
        if heights[plane[0]] < heights.max() - depth - 1e-6:
            break
        stripped = tuple(index for upper_plane in planes[:depth_index] for index in upper_plane)
        stripped_charge = sum(window[index].charge for index in stripped)
        for count in range(1, len(plane) + 1):
            size = len(stripped) + count
            if size > max_removed or math.comb(len(plane), count) > max_combinations:
                break
            for subset in itertools.combinations(plane, count):
                total = stripped_charge + sum(window[index].charge for index in subset)
                if abs(total - target) <= CHARGE_TOLERANCE:
                    solutions.setdefault(size, []).append(stripped + subset)
    if not solutions:
        return

    seen: set[tuple[int, ...]] = set()
    for top_set in solutions[min(solutions)]:
        # Different face-reversing operations pair the top units with
        # different bottom units, giving different symmetric slabs.
        for face in faces:
            removed = frozenset(top_set) | frozenset(int(face[index]) for index in top_set)
            if len(removed) != 2 * len(top_set):
                continue
            key = min(tuple(sorted(int(mapping[index]) for index in removed))
                      for mapping, _ in unit_permutations)
            if key in seen:
                continue
            if len(seen) >= max_variants:
                warnings.warn(
                    f"More than {max_variants} symmetry-distinct ways to neutralise a "
                    f"{len(window)}-unit window; only the first {max_variants} are kept. "
                    "Raise max_variants_per_window to see more.",
                    stacklevel=3,
                )
                return
            seen.add(key)
            yield removed


# =============================================================================
# PLOTTING
# =============================================================================

def plot_terminations(
    terminations: t.Sequence[Termination],
    filename: str | None = None,
    repeat: int = 2,
    columns: int | None = None,
    show_removed: bool = True,
):
    """Draw a side view of each termination.

    Atoms are drawn to scale (covalent radii, Jmol colours) looking along the
    second surface vector; the cell is repeated ``repeat`` times along the
    first. Atoms removed to neutralise the slab are drawn as dashed outlines.

    Args:
        terminations: A :class:`TerminationSet` or list of terminations.
        filename: Save the figure here (PNG, PDF, SVG, ...) if given.
        repeat: In-plane repetitions shown.
        columns: Panels per row (default: up to 4).
        show_removed: Draw the removed atoms.

    Returns:
        The matplotlib ``Figure``.
    """
    from ase.data import atomic_numbers, covalent_radii
    from ase.data.colors import jmol_colors
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle

    terminations = list(terminations)
    if not terminations:
        raise ValueError('nothing to plot')
    columns = columns or min(len(terminations), 4)
    rows = math.ceil(len(terminations) / columns)

    def frame(termination):
        matrix = termination.structure.lattice.matrix
        along = matrix[0] / np.linalg.norm(matrix[0])
        normal = matrix[2] / np.linalg.norm(matrix[2])
        depth = np.cross(normal, along)
        return matrix, along, depth, normal

    extents = []
    for termination in terminations:
        matrix, along, _, normal = frame(termination)
        heights = termination.structure.cart_coords @ normal
        extents.append((np.linalg.norm(matrix[0]) * repeat, heights.max() - heights.min()))
    width = max(extent[0] for extent in extents) + 2.0
    height = max(extent[1] for extent in extents) + 5.0
    scale = 3.2 / max(width, height)
    figure = Figure(figsize=(columns * max(2.6, width * scale), rows * (height * scale + 0.9)),
                    constrained_layout=True)
    axes = figure.subplots(rows, columns, squeeze=False)
    species_seen: dict[str, tuple] = {}

    for panel, termination in enumerate(terminations):
        axis = axes[panel // columns][panel % columns]
        matrix, along, depth, normal = frame(termination)
        sites = [(site.specie.symbol, site.coords, False) for site in termination.structure]
        if show_removed:
            sites += [(symbol, np.array(position), True) for symbol, position in termination.removed_sites]
        bottom = (termination.structure.cart_coords @ normal).min()
        drawn = []
        for symbol, position, removed in sites:
            for copy in range(repeat):
                point = position + copy * matrix[0]
                drawn.append((float(point @ depth), float(point @ along), float(point @ normal - bottom),
                              symbol, removed))
        for _, x, z, symbol, removed in sorted(drawn, key=lambda item: -item[0]):
            number = atomic_numbers[symbol]
            colour = tuple(jmol_colors[number])
            radius = 0.55 * covalent_radii[number]
            species_seen.setdefault(symbol, colour)
            if removed:
                axis.add_patch(Circle((x, z), radius, facecolor='none', edgecolor=colour,
                                      linestyle='--', linewidth=1.2, zorder=3))
            else:
                axis.add_patch(Circle((x, z), radius, facecolor=colour, edgecolor='0.25',
                                      linewidth=0.5, zorder=2))
        axis.set_xlim(-1.0, width - 1.0)
        axis.set_ylim(-2.5, height - 2.5)
        axis.set_aspect('equal')
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_color('0.8')
        stoich = 'stoichiometric' if termination.is_stoichiometric else 'non-stoichiometric'
        axis.set_title(f"{termination.label}  {termination.formula}\n{stoich}, {termination.thickness:.1f} Å",
                       fontsize=9)
        axis.text(0.5, 0.02, termination.origin, transform=axis.transAxes, ha='center',
                  va='bottom', fontsize=7, color='0.35', wrap=True)

    for panel in range(len(terminations), rows * columns):
        axes[panel // columns][panel % columns].set_visible(False)

    handles = [Line2D([], [], marker='o', linestyle='', markerfacecolor=colour, markeredgecolor='0.25',
                      markersize=8, label=symbol) for symbol, colour in species_seen.items()]
    if show_removed and any(termination.removed_sites for termination in terminations):
        handles.append(Line2D([], [], marker='o', linestyle='', markerfacecolor='none',
                              markeredgecolor='0.4', markersize=8, label='removed'))
    figure.legend(handles=handles, loc='outside lower center', ncol=len(handles), frameon=False, fontsize=8)
    if filename:
        figure.savefig(filename, dpi=150)
    return figure


# =============================================================================
# COMMAND LINE
# =============================================================================

def _parse_miller(values: t.Sequence[str]) -> tuple[int, int, int]:
    text = ' '.join(values).replace(',', ' ').strip('()[] ')
    parts = text.split()
    if len(parts) == 1 and re.fullmatch(r'-?\d-?\d-?\d', parts[0]):
        parts = re.findall(r'-?\d', parts[0])
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"Miller index must have three integers, got {text!r}")
    return tuple(int(part) for part in parts)


def _parse_states(text: str) -> dict[str, float]:
    states = {}
    for item in text.replace(' ', '').split(','):
        element, _, value = item.partition('=')
        if not value:
            raise argparse.ArgumentTypeError(f"expected Element=state, got {item!r}")
        states[element] = float(value)
    return states


def _parse_bond(text: str) -> tuple[tuple[str, str], float]:
    match = re.fullmatch(r'([A-Z][a-z]?)-([A-Z][a-z]?):([0-9.]+)', text.strip())
    if not match:
        raise argparse.ArgumentTypeError(f"expected A-B:cutoff such as P-O:1.9, got {text!r}")
    return (match.group(1), match.group(2)), float(match.group(3))


def _parse_cell(text: str) -> tuple[int, int]:
    match = re.fullmatch(r'(\d+)[x,](\d+)', text.strip())
    if not match:
        raise argparse.ArgumentTypeError(f"expected n1xn2 such as 2x1, got {text!r}")
    return int(match.group(1)), int(match.group(2))


def main(argv: t.Sequence[str] | None = None) -> int:
    """Command-line entry point: ``psteros-terminations``."""
    parser = argparse.ArgumentParser(
        prog='psteros-terminations',
        description='List the symmetric, charge-neutral slab terminations of a crystal surface.',
        epilog='Example: psteros-terminations Ag3PO4.cif 110 --oxidation Ag=1,P=5,O=-2 '
               '--keep P-O:1.9 --write slabs/ --plot slabs.png',
    )
    parser.add_argument('bulk', help='bulk structure file (CIF, POSCAR, ...)')
    parser.add_argument('miller', nargs='+', help='Miller index: 110, "1 1 0" or 1,1,0')
    parser.add_argument('--thickness', type=float, default=10.0,
                        help='minimum slab thickness in Å (default: 10)')
    parser.add_argument('--vacuum', type=float, default=15.0, help='vacuum in Å (default: 15)')
    parser.add_argument('--oxidation', type=_parse_states,
                        help='oxidation states, e.g. Ag=1,P=5,O=-2 (default: guessed)')
    parser.add_argument('--keep', type=_parse_bond, action='append', default=[], metavar='A-B:CUTOFF',
                        help='bond never to break, e.g. P-O:1.9 (repeatable)')
    parser.add_argument('--supercell', type=_parse_cell,
                        help='surface cell, e.g. 2x1 (default: smallest of 1x1, 2x1, 1x2, 2x2 that works)')
    parser.add_argument('--write', metavar='DIR', help='write the slabs and terminations.json here')
    parser.add_argument('--format', default='poscar', choices=['poscar', 'cif', 'extxyz'],
                        help='file format for --write (default: poscar)')
    parser.add_argument('--plot', metavar='FILE', help='save side views (PNG, PDF, SVG)')
    args = parser.parse_args(argv)

    try:
        miller = _parse_miller(args.miller)
    except argparse.ArgumentTypeError as error:
        parser.error(str(error))
    bulk = Structure.from_file(args.bulk)
    try:
        terminations = find_charge_neutral_terminations(
            bulk, miller, args.thickness, args.vacuum,
            oxidation_states=args.oxidation, unit_bonds=dict(args.keep) or None,
            supercell=args.supercell,
        )
    except (NoChargeNeutralTerminationError, ValueError) as error:
        print(f"No termination found.\n{error}", file=sys.stderr)
        return 1
    print(terminations.summary())
    if args.write:
        paths = terminations.write(args.write, fmt=args.format)
        print(f"\nWrote {len(paths) - 1} slabs and terminations.json to {args.write}")
    if args.plot:
        terminations.plot(args.plot)
        print(f"Saved side views to {args.plot}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
