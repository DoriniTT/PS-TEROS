"""Reproducible starting structures for the maintained SnO2(110) campaign."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from typing import Literal

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
