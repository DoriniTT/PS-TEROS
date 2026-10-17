"""Reproducible starting structures for the maintained SnO2(110) campaign."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from typing import Literal

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import Slab, SlabGenerator


Termination = Literal["o", "sno", "sn2o"]


@dataclass(frozen=True)
class SlabIdentity:
    """Published structural identity carried beside every generated slab."""

    miller_index: tuple[int, int, int]
    termination: Termination
    triple_layers: int
    vacuum_angstrom: float


def rutile_sno2_bulk(*, a: float = 4.737, c: float = 3.186) -> Structure:
    """Return conventional rutile SnO2 using experimentally close lattice values."""

    lattice = Lattice.tetragonal(a, c)
    return Structure.from_spacegroup(
        "P4_2/mnm",
        lattice,
        ["Sn", "O"],
        [[0.0, 0.0, 0.0], [0.3056, 0.3056, 0.0]],
    )


def alpha_sn_bulk(*, a: float = 6.489) -> Structure:
    """Return conventional diamond-cubic alpha-Sn for a zero-K reference."""

    return Structure.from_spacegroup(
        "Fd-3m",
        Lattice.cubic(a),
        ["Sn"],
        [[0.0, 0.0, 0.0]],
    )


def litharge_sno_bulk(*, a: float = 3.802, c: float = 4.838, z_sn: float = 0.2367) -> Structure:
    """Return a tetragonal litharge SnO starting cell.

    This is a reproducible experimental-neighbourhood starting geometry. The
    campaign relaxes both positions and cell before using its reference energy.
    """

    return Structure(
        Lattice.tetragonal(a, c),
        ["Sn", "Sn", "O", "O"],
        [
            [0.25, 0.25, z_sn],
            [0.75, 0.75, 1.0 - z_sn],
            [0.75, 0.25, 0.0],
            [0.25, 0.75, 0.0],
        ],
    )


def triplet_o2_cell(*, cell_length: float = 18.0, bond_length: float = 1.208) -> Structure:
    """Return an isolated O2 molecule in a cubic periodic cell.

    The geometry is suitable for a spin-polarised triplet calculation; the
    triplet occupation is an explicit QE input, not a structure property.
    """

    if cell_length <= bond_length:
        raise ValueError("cell_length must exceed bond_length")
    half_fractional_bond = bond_length / (2.0 * cell_length)
    return Structure(
        Lattice.cubic(cell_length),
        ["O", "O"],
        [
            [0.5, 0.5, 0.5 - half_fractional_bond],
            [0.5, 0.5, 0.5 + half_fractional_bond],
        ],
    )


def _select_symmetric_slab(
    slabs: list[Slab], *, vacuum_angstrom: float
) -> Slab:
    if not slabs:
        raise RuntimeError("pymatgen did not generate a SnO2(110) slab")
    # The published O-terminated reference is stoichiometric.  Pymatgen also
    # returns oxygen-rich and oxygen-poor symmetric cuts; choosing one merely
    # because it has fewer sites would silently change the thermodynamic model.
    stoichiometric = [
        slab
        for slab in slabs
        if slab.composition["O"] == 2 * slab.composition["Sn"]
    ]
    if not stoichiometric:
        formulas = [slab.composition.formula for slab in slabs]
        raise RuntimeError(
            "no stoichiometric symmetric SnO2(110) termination was generated; "
            f"candidates were {formulas}"
        )
    return min(stoichiometric, key=len)


def sno2_110_slab(
    *,
    termination: Termination = "o",
    triple_layers: int = 9,
    vacuum_angstrom: float = 20.0,
    a: float = 4.737,
    c: float = 3.186,
) -> tuple[Slab, SlabIdentity]:
    """Build a symmetric 1x1 rutile SnO2(110) starting slab.

    ``sno`` and ``sn2o`` use paired oxygen removal from the two equivalent
    outermost surfaces.  The pairing keeps the cell stoichiometrically and
    electrostatically symmetric; relaxation remains the authority for the
    final structure.
    """

    if termination not in ("o", "sno", "sn2o"):
        raise ValueError(f"termination must be 'o', 'sno' or 'sn2o', got {termination!r}")
    if isinstance(triple_layers, bool) or not isinstance(triple_layers, int):
        raise TypeError(f"triple_layers must be an integer, got {triple_layers!r}")
    if triple_layers < 3 or triple_layers % 2 == 0:
        raise ValueError("triple_layers must be an odd integer of at least 3")
    if vacuum_angstrom <= 0:
        raise ValueError("vacuum_angstrom must be positive")
    bulk = rutile_sno2_bulk(a=a, c=c)
    # A rutile (110) triple layer is approximately a/√2 thick.  This maps the
    # scientist-facing layer count to SlabGenerator's Angstrom requirement.
    minimum_slab = triple_layers * a / (2.0**0.5)
    generator = SlabGenerator(
        initial_structure=bulk,
        miller_index=(1, 1, 0),
        min_slab_size=minimum_slab,
        min_vacuum_size=vacuum_angstrom,
        center_slab=True,
        in_unit_planes=False,
        primitive=False,
        lll_reduce=True,
    )
    slab = _select_symmetric_slab(generator.get_slabs(symmetrize=True), vacuum_angstrom=vacuum_angstrom)
    if termination != "o":
        slab = _remove_outer_oxygen_pairs(slab, pairs=1 if termination == "sno" else 2)
    identity = SlabIdentity((1, 1, 0), termination, triple_layers, vacuum_angstrom)
    slab.add_site_property("psteros_termination", [termination] * len(slab))
    return slab, identity


def _remove_outer_oxygen_pairs(slab: Slab, *, pairs: int) -> Slab:
    """Remove equivalent oxygen atoms at the two exterior surfaces.

    This is intentionally geometric rather than index-based so it remains
    stable across harmless pymatgen ordering changes.
    """

    result = slab.copy()
    normal = result.normal
    projections = [site.coords.dot(normal) for site in result]
    oxygen = [index for index, site in enumerate(result) if site.specie.symbol == "O"]
    if len(oxygen) < 2 * pairs:
        raise RuntimeError("not enough oxygen atoms for requested reduction")
    low = sorted(oxygen, key=lambda index: projections[index])[:pairs]
    high = sorted(oxygen, key=lambda index: projections[index], reverse=True)[:pairs]
    result.remove_sites(sorted(set(low + high), reverse=True))
    return result


def central_sites(structure: Structure, half_width: float = 1.5) -> tuple[int, ...]:
    """Indices of the sites within ``half_width`` Å of the slab mid-plane.

    The surface normal is taken normal to the first two lattice vectors (the
    psteros and pymatgen slab convention). Pass the result as
    ``CalculationOverride(fixed_sites=...)`` to keep the centre of a slab at
    its bulk positions during the relaxation.
    """

    import numpy as np

    if half_width <= 0:
        raise ValueError("half_width must be positive")
    a, b = structure.lattice.matrix[0], structure.lattice.matrix[1]
    normal = np.cross(a, b)
    normal /= np.linalg.norm(normal)
    heights = [float(np.dot(site.coords, normal)) for site in structure]
    middle = (max(heights) + min(heights)) / 2.0
    return tuple(index for index, height in enumerate(heights) if abs(height - middle) < half_width)


RutileTermination = Literal[
    "o_bridge", "o_cus", "metal", "o_cus_half", "o_bridge_half", "o_cus_no_bridge"
]
RUTILE_110_TERMINATIONS: tuple[str, ...] = (
    "o_bridge", "o_cus", "metal", "o_cus_half", "o_bridge_half", "o_cus_no_bridge"
)


def rutile_bulk(metal: str, *, a: float, c: float, u: float = 0.3056) -> Structure:
    """Return a conventional rutile MO2 cell (P4_2/mnm, M at 2a, O at 4f (u, u, 0))."""

    if a <= 0 or c <= 0 or not 0 < u < 0.5:
        raise ValueError("rutile lattice needs a > 0, c > 0 and 0 < u < 0.5")
    return Structure.from_spacegroup("P4_2/mnm", Lattice.tetragonal(a, c), [metal, "O"], [[0, 0, 0], [u, u, 0]])


def _heights(structure: Structure) -> tuple[np.ndarray, np.ndarray]:
    a, b = structure.lattice.matrix[0], structure.lattice.matrix[1]
    normal = np.cross(a, b)
    return structure.cart_coords @ (normal / np.linalg.norm(normal)), normal / np.linalg.norm(normal)


def _rutile_110_sites(slab: Structure, metal: str) -> dict[str, list[list[int]]]:
    """Per surface (bottom, top): indices of the bridging O and of the five-fold (cus) metal atoms."""

    heights, _ = _heights(slab)
    metals = [index for index, site in enumerate(slab) if site.specie.symbol == metal]
    oxygens = [index for index, site in enumerate(slab) if site.specie.symbol == "O"]
    found: dict[str, list[list[int]]] = {"bridge": [], "cus": []}
    for sign in (-1.0, 1.0):  # bottom, top
        outer = max(sign * heights[index] for index in metals)
        plane = [index for index in metals if outer - sign * heights[index] < 0.3]
        bridge = [index for index in oxygens if sign * heights[index] > outer + 0.5]
        cus = [
            index for index in plane
            if not any(slab.get_distance(index, other) < 2.3 for other in bridge)
        ]
        if not bridge or len(cus) * 2 != len(plane):
            raise RuntimeError("could not identify the bridging O and cus metal rows of the rutile (110) slab")
        found["bridge"].append(sorted(bridge))
        found["cus"].append(sorted(cus))
    return found


def _inversion_partner(slab: Structure, index: int, candidates: list[int]) -> int:
    """The candidate site closest to the image of ``index`` through the centre of the slab."""

    centre = slab.frac_coords.mean(axis=0)
    centre[2] = (slab.frac_coords[:, 2].max() + slab.frac_coords[:, 2].min()) / 2.0
    image = 2.0 * centre - slab.frac_coords[index]
    distances = [slab.lattice.get_distance_and_image(image, slab.frac_coords[other])[0] for other in candidates]
    return candidates[int(np.argmin(distances))]


def rutile_110_slab(
    *,
    metal: str,
    termination: str = "o_bridge",
    triple_layers: int = 3,
    vacuum_angstrom: float = 15.0,
    a: float,
    c: float,
    u: float = 0.3056,
    repeat_c: int = 1,
    cus_bond_angstrom: float = 1.70,
) -> tuple[Slab, SlabIdentity]:
    """Build a symmetric rutile MO2(110) slab with one of six terminations.

    Starting from the stoichiometric O-bridge slab of pymatgen (outermost
    plane: bridging O rows above six-fold M, five-fold "cus" M exposed),
    on both surfaces:

    ``"o_bridge"``        stoichiometric (Reuter and Scheffler's RuO2(110)-O\\ :sup:`bridge`)
    ``"o_cus"``           one O atop every cus M, ``cus_bond_angstrom`` above it (O-rich)
    ``"metal"``           the bridging O removed (O-poor; their RuO2(110)-Ru)
    ``"o_cus_half"``      O atop every second cus M (cell doubled along [001])
    ``"o_bridge_half"``   every second bridging O removed (cell doubled along [001])
    ``"o_cus_no_bridge"`` O-cus added and bridging O removed (stoichiometric, non-bulk-like)

    The cell is ``c`` (along [001], repeated ``repeat_c`` times) by ``a*sqrt(2)``
    (along [1-10]). The bottom surface gets the inversion image of every
    change at the top, so the slab stays symmetric (checked). The atoms that
    were added sit at a guessed height; a relaxation finds their real place.
    """

    if termination not in RUTILE_110_TERMINATIONS:
        raise ValueError(f"termination must be one of {RUTILE_110_TERMINATIONS}, got {termination!r}")
    if isinstance(triple_layers, bool) or not isinstance(triple_layers, int):
        raise TypeError(f"triple_layers must be an integer, got {triple_layers!r}")
    if triple_layers < 3 or triple_layers % 2 == 0:
        raise ValueError("triple_layers must be an odd integer of at least 3")
    if vacuum_angstrom <= 0 or cus_bond_angstrom <= 0:
        raise ValueError("vacuum_angstrom and cus_bond_angstrom must be positive")
    if isinstance(repeat_c, bool) or not isinstance(repeat_c, int) or repeat_c < 1:
        raise ValueError("repeat_c must be a positive integer")
    half = termination in ("o_cus_half", "o_bridge_half")
    if half and repeat_c % 2:
        repeat_c *= 2

    bulk = rutile_bulk(metal, a=a, c=c, u=u)
    generator = SlabGenerator(
        initial_structure=bulk, miller_index=(1, 1, 0), min_slab_size=triple_layers * a / (2.0**0.5),
        min_vacuum_size=vacuum_angstrom, center_slab=True, in_unit_planes=False, primitive=False,
        lll_reduce=True,
    )
    candidates = [
        slab for slab in generator.get_slabs(symmetrize=True)
        if slab.composition["O"] == 2 * slab.composition[metal]
    ]
    if not candidates:
        raise RuntimeError(f"no stoichiometric symmetric {metal}O2(110) slab was generated")
    slab = min(candidates, key=len)
    lengths = slab.lattice.abc
    along_c = int(np.argmin([abs(length - c) for length in lengths[:2]]))
    scaling = [1, 1, 1]
    scaling[along_c] = repeat_c
    slab.make_supercell(scaling)

    sites = _rutile_110_sites(slab, metal)
    heights, normal = _heights(slab)
    remove: list[int] = []
    add: list[np.ndarray] = []
    bottom_bridge, top_bridge = sites["bridge"]
    bottom_cus, top_cus = sites["cus"]
    if termination in ("metal", "o_cus_no_bridge"):
        remove += top_bridge + bottom_bridge
    if termination == "o_bridge_half":
        chosen = top_bridge[::2]
        remove += chosen + [_inversion_partner(slab, index, bottom_bridge) for index in chosen]
    if termination in ("o_cus", "o_cus_no_bridge", "o_cus_half"):
        chosen = top_cus[::2] if termination == "o_cus_half" else top_cus
        for index in chosen:
            partner = _inversion_partner(slab, index, bottom_cus)
            add.append(slab.cart_coords[index] + cus_bond_angstrom * normal)
            add.append(slab.cart_coords[partner] - cus_bond_angstrom * normal)
    for position in add:
        slab.append("O", position, coords_are_cartesian=True)
    if remove:
        slab.remove_sites(sorted(set(remove), reverse=True))
    if not slab.is_symmetric():
        raise RuntimeError(f"the {metal}O2(110) {termination} slab lost its symmetry")
    identity = SlabIdentity((1, 1, 0), termination, triple_layers, vacuum_angstrom)
    slab.add_site_property("psteros_termination", [termination] * len(slab))
    return slab, identity
