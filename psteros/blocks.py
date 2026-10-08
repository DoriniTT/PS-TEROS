"""Calculation blocks: the steps of a reference calculation, linked by name.

A block is one calculation applied to every reference of a graph, such as a
relaxation, a static calculation or a vibrational analysis.  Blocks run in the
order given; each takes the structure of an earlier block
(``structure_from``, by default the block just before it) so a protocol reads
like a list of steps::

    blocks = (
        psteros.Relax(incar={"ibrion": 2, "nsw": 100, "ediffg": -0.005}),
        psteros.Static(),
        psteros.Vibrations(),
    )

Blocks hold no AiiDA nodes.  The INCAR of one block is, in increasing
priority: the recipe INCAR, the block defaults, ``Block.incar``, the
reference's ``override`` and its ``block_overrides[block.name]``.  Tags a
block needs to be what it is (``NSW = 0`` for :class:`Static`, ``IBRION``,
``POTIM``, ``NFREE``, ``NSW`` for :class:`Vibrations`) are applied last.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, ClassVar, Literal, Mapping, Sequence, Union

Phase = Literal["gas", "solid"]

_NAME = re.compile(r"\w+", re.ASCII)


def _check_block(block: Any) -> None:
    if not isinstance(block.name, str) or not _NAME.fullmatch(block.name):
        raise ValueError(f"block name must be letters, digits and underscores, got {block.name!r}")
    if block.structure_from is not None and not isinstance(block.structure_from, str):
        raise ValueError(f"{block.name}: structure_from must be a block name or None")
    if not isinstance(block.incar, Mapping):
        raise TypeError(f"{block.name}: incar must be a mapping of INCAR tags")
    object.__setattr__(block, "incar", {str(key).lower(): value for key, value in block.incar.items()})


@dataclass(frozen=True)
class Relax:
    """Relax the structure; later blocks take the relaxed structure.

    ``incar`` must lead, together with the recipe, to ``NSW > 0``.  A cell
    relaxation of a bulk (``ISIF = 3``) belongs in that reference's
    ``block_overrides``, so that molecules keep a fixed box.
    """

    name: str = "relax"
    incar: Mapping[str, Any] = field(default_factory=dict)
    structure_from: str | None = None

    kind: ClassVar[str] = "relax"
    moves_ions: ClassVar[bool] = True

    def __post_init__(self) -> None:
        _check_block(self)

    def defaults(self, phase: Phase) -> dict[str, Any]:
        return {}

    def required(self, phase: Phase) -> dict[str, Any]:
        return {}


@dataclass(frozen=True)
class Static:
    """Single-point energy on the structure of an earlier block (``NSW = 0``)."""

    name: str = "static"
    incar: Mapping[str, Any] = field(default_factory=dict)
    structure_from: str | None = None

    kind: ClassVar[str] = "static"
    moves_ions: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _check_block(self)

    def defaults(self, phase: Phase) -> dict[str, Any]:
        return {"ibrion": -1}

    def required(self, phase: Phase) -> dict[str, Any]:
        return {"nsw": 0}


@dataclass(frozen=True)
class Vibrations:
    """Harmonic frequencies by VASP finite differences on a relaxed structure.

    Runs ``IBRION = 5`` for a gas (all atoms displaced) and ``IBRION = 6`` for
    a solid (displacements reduced by symmetry) unless ``ibrion`` is given;
    ``potim`` is the displacement (A) and ``nfree`` the number of
    displacements per degree of freedom.  A solid is displaced in the
    supercell given by its :class:`~psteros.references.ReferenceSystem`.
    The relaxation before it must be tight (forces of a few meV/A) and the
    electronic convergence strict (``EDIFF`` of 1e-7 or less), otherwise
    spurious imaginary modes appear.

    Displacing atoms lowers the symmetry, so VASP has to change its k-point
    set during the run, and it refuses to do that with band parallelisation
    (``NCORE > 1`` or ``NPAR``): it stops with "VASP internal routines have
    requested a change of the k-point set ... remove the tag NPAR".  The block
    therefore defaults to ``ISIF = 2`` and, unless ``incar`` says otherwise,

    * for a **gas** ``ISYM = 0`` (no symmetry, so the k-point set cannot
      change; a molecule in a box is cheap and has a single k-point), keeping
      the recipe's ``NCORE``, because ``NCORE = 1`` on many ranks would pad the
      bands to one per rank (128 bands for 6 occupied ones) and VASP's
      diagonalisation then fails (``EDDDAV: Call to ZHEGV failed``);
    * for a **solid** ``NCORE = 1``, which keeps the symmetry-reduced set of
      displacements (``ISYM = 0`` would displace every atom in every
      direction) and is harmless for a supercell with hundreds of bands.

    ``NPAR`` must not be in the recipe INCAR.
    """

    name: str = "vibrations"
    incar: Mapping[str, Any] = field(default_factory=dict)
    structure_from: str | None = None
    ibrion: int | None = None
    potim: float = 0.015
    nfree: int = 2

    kind: ClassVar[str] = "vibrations"
    moves_ions: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _check_block(self)
        if self.ibrion is not None and self.ibrion not in (5, 6):
            raise ValueError(f"{self.name}: ibrion must be 5 or 6 (finite differences), got {self.ibrion!r}")
        if self.potim <= 0:
            raise ValueError(f"{self.name}: potim must be positive, got {self.potim!r}")
        if self.nfree not in (1, 2, 4):
            raise ValueError(f"{self.name}: nfree must be 1, 2 or 4, got {self.nfree!r}")
        reserved = {"ibrion", "potim", "nfree", "nsw"}.intersection(self.incar)
        if reserved:
            raise ValueError(
                f"{self.name}: set {sorted(reserved)} through the Vibrations fields, not incar"
            )

    def defaults(self, phase: Phase) -> dict[str, Any]:
        # ISIF >= 3 would make IBRION = 6 also strain the cell (elastic constants).
        # VASP cannot change its k-point set under band parallelisation, which
        # displaced (lower-symmetry) cells need: no symmetry for a gas, NCORE = 1 for a solid.
        return {"isif": 2, "isym": 0} if phase == "gas" else {"isif": 2, "ncore": 1}

    def required(self, phase: Phase) -> dict[str, Any]:
        ibrion = self.ibrion if self.ibrion is not None else (5 if phase == "gas" else 6)
        return {"ibrion": ibrion, "potim": self.potim, "nfree": self.nfree, "nsw": 1}


Block = Union[Relax, Static, Vibrations]

DEFAULT_BLOCKS: tuple[Block, ...] = (Relax(), Static(), Vibrations())


def check_blocks(blocks: Sequence[Block]) -> tuple[Block, ...]:
    """Validate a block sequence and return it as a tuple.

    Names must be unique and ``structure_from`` must name an earlier block.
    """

    blocks = tuple(blocks)
    if not blocks:
        raise ValueError("at least one block is required")
    seen: list[str] = []
    for block in blocks:
        if not isinstance(block, (Relax, Static, Vibrations)):
            raise TypeError(f"blocks must be Relax, Static or Vibrations, got {type(block).__name__}")
        if block.name in seen:
            raise ValueError(f"duplicate block name {block.name!r}")
        if block.structure_from is not None and block.structure_from not in seen:
            raise ValueError(
                f"{block.name}: structure_from must name an earlier block, got {block.structure_from!r}"
            )
        seen.append(block.name)
    return blocks
