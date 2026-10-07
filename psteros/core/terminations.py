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

Formal charge is a necessary, not a sufficient, condition: the band gap of
each candidate should still be confirmed with a single-point calculation.
Covalent semiconductors (Si, GaAs) need dangling-bond counting and surface
reconstructions, which are outside the scope of this module.
"""

from __future__ import annotations

import itertools
import math
import typing as t
import warnings
from dataclasses import dataclass

import numpy as np
from pymatgen.core import Composition, Lattice, Structure
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
    """

    structure: Structure
    formal_charge: float
    parent_formula: str
    parent_formal_charge: float
    removed_per_face: tuple[str, ...]
    thickness: float
    is_stoichiometric: bool

    @property
    def formula(self) -> str:
        return self.structure.composition.formula.replace(' ', '')

    def to_dict(self) -> dict:
        return {
            'formula': self.formula,
            'formal_charge': self.formal_charge,
            'parent_formula': self.parent_formula,
            'parent_formal_charge': self.parent_formal_charge,
            'removed_per_face': list(self.removed_per_face),
            'thickness': self.thickness,
            'is_stoichiometric': self.is_stoichiometric,
        }


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
        if not guesses:
            raise ValueError(
                f"Could not guess oxidation states for {bulk.composition.reduced_formula}; "
                "pass oxidation_states explicitly."
            )
        oxidation_states = guesses[0]
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
        operations = SpacegroupAnalyzer(structure, symprec=symprec).get_symmetry_operations()
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
    operations = SpacegroupAnalyzer(bulk, symprec=symprec).get_point_group_operations(cartesian=True)
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
) -> tuple[Structure, float, np.ndarray]:
    """Assemble units into a slab with c along ``normal``, centred in the cell.

    Returns the slab, its thickness and, for every site, the position of its
    unit in ``units``.
    """
    coords = np.vstack([unit.coords for unit in units])
    species = [species_of(index) for unit in units for index in unit.indices]
    site_unit = np.repeat(np.arange(len(units)), [len(unit.indices) for unit in units])
    heights = coords @ normal
    thickness = float(heights.max() - heights.min())
    length = thickness + vacuum
    lattice = Lattice(np.vstack([in_plane, normal * length]))
    shifted = coords + normal * (0.5 * length - 0.5 * (heights.max() + heights.min()))
    slab = Structure(lattice, species, shifted, coords_are_cartesian=True, to_unit_cell=True)
    return slab, thickness, site_unit


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
    for operation in SpacegroupAnalyzer(bulk, symprec=symprec).get_symmetry_operations(cartesian=True):
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
) -> Structure:
    """A thick slab with c normal to the surface."""
    generator = SlabGenerator(
        bulk,
        miller_index,
        min_slab_size=thickness,
        min_vacuum_size=vacuum,
        center_slab=True,
        in_unit_planes=False,
        primitive=True,
        max_normal_search=max_normal_search,
        lll_reduce=lll_reduce,
    )
    # Any cut will do: units broken at the parent's edges are discarded later.
    parent = Structure.from_sites(generator.get_slab(shift=0.0).get_orthogonal_c_slab())
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
    supercell: tuple[int, int] = (1, 1),
    surface_depth: float | None = None,
    max_removed_per_face: int | None = None,
    max_combinations: int = 20000,
    max_variants_per_window: int = 20,
    plane_tolerance: float = 0.3,
    symprec: float = 0.1,
    lll_reduce: bool = True,
    max_normal_search: int | None = None,
) -> list[Termination]:
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
        supercell: In-plane repetition (n1, n2) applied before searching.
            Larger cells allow partial surface coverages, e.g. removing one
            cation in four.
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
        Distinct terminations, sorted by number of removed units and then by
        thickness. Each slab has c normal to the surface.

    Raises:
        NoChargeNeutralTerminationError: If (hkl) is polar or no termination
            is found.
    """
    miller_index = tuple(int(value) for value in miller_index)
    states = resolve_oxidation_states(bulk, oxidation_states)
    bonds = _normalise_bonds(unit_bonds)
    supercell = (int(supercell[0]), int(supercell[1]))
    if max_removed_per_face is None:
        max_removed_per_face = 4 * supercell[0] * supercell[1]

    if not has_face_reversing_operation(bulk, miller_index, symprec):
        raise NoChargeNeutralTerminationError(
            f"{miller_index} is a polar direction of {bulk.composition.reduced_formula}: "
            "no bulk symmetry operation reverses the surface normal, so the two faces "
            "of a slab can never be equivalent. Use an asymmetric slab with a dipole "
            "correction or a passivated back face instead."
        )

    # Bulk units define what a complete unit looks like.
    unit_types = {unit.composition for unit in _find_units(bulk, bonds, states)}

    repeat = normal_repeat(bulk, miller_index)
    depth = 0.5 * repeat if surface_depth is None else float(surface_depth)
    # Removing units thins a window by up to 2 * depth, so windows are taken up
    # to that much thicker. Final slabs must lie within one repeat above the
    # minimum. The parent adds a few repeats of margin for broken edge units.
    lower, upper = min_slab_thickness, min_slab_thickness + repeat
    parent = _oriented_parent(
        bulk, miller_index, upper + 2.0 * depth + 4.0 * repeat, min_vacuum_thickness,
        lll_reduce, max_normal_search, supercell,
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
            slab, thickness, site_unit = _assemble(
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
                repaired, repaired_thickness, _ = _assemble(
                    kept, species_of, in_plane, normal, min_vacuum_thickness,
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
                )))

    if not candidates:
        raise NoChargeNeutralTerminationError(
            f"No symmetric, charge-neutral termination found for {miller_index} "
            f"with oxidation states {states}.\n"
            f"  Windows tried: {windows_tried}; symmetric: {symmetric_windows}; "
            f"formal charges of symmetric windows: {sorted(charges_seen)}.\n"
            "  Try a larger supercell (partial coverages), a larger surface_depth "
            "or max_removed_per_face, or check oxidation_states and unit_bonds."
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
