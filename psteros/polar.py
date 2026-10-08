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


def tetrahedral_coordination(bulk: Any, tolerance: float = 0.15) -> dict[str, int]:
    """Number of nearest neighbours of each element (within ``tolerance`` of the shortest bond)."""

    neighbours = bulk.get_all_neighbors(4.0)
    shortest = min(neighbour.nn_distance for site in neighbours for neighbour in site)
    counts: dict[str, set[int]] = {}
    for site, found in zip(bulk, neighbours):
        count = sum(1 for neighbour in found if neighbour.nn_distance <= shortest * (1 + tolerance))
        counts.setdefault(site.specie.symbol, set()).add(count)
    return {element: (values.pop() if len(values) == 1 else -1) for element, values in counts.items()}
