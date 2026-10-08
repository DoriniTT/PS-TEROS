"""Surface free-energy analysis independent of a workflow engine."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping


EV_PER_ANGSTROM2_TO_J_PER_M2 = 16.021_766_34


@dataclass(frozen=True)
class SurfaceEnergyPoint:
    """A surface energy at one chemical-potential offset.

    The offset is that of the element on the phase-diagram axis (oxygen for
    an oxide); ``delta_mu_oxygen_ev`` keeps its historical name and
    ``delta_mu_ev`` is the element-neutral alias.
    """

    delta_mu_oxygen_ev: float
    gamma_ev_per_angstrom2: float

    @property
    def delta_mu_ev(self) -> float:
        return self.delta_mu_oxygen_ev

    @property
    def gamma_j_per_m2(self) -> float:
        return self.gamma_ev_per_angstrom2 * EV_PER_ANGSTROM2_TO_J_PER_M2


def surface_energy_binary_equilibrium(
    *,
    slab_energy_ev: float,
    n_other: int,
    n_variable: int,
    bulk_formula_energy_ev: float,
    variable_reference_energy_ev: float,
    delta_mu_ev: float,
    surface_area_angstrom2: float,
    surfaces: int = 2,
    formula_unit: tuple[int, int] = (1, 1),
    reservoir_correction_ev: float = 0.0,
) -> SurfaceEnergyPoint:
    """Return gamma for an A_xB_y slab in equilibrium with its bulk compound.

    B is the element on the phase-diagram axis, with
    ``mu_B = variable_reference_energy_ev + delta_mu_ev`` (for an oxide the
    reference is E(O2)/2; for GaAs it can be the energy per atom of bulk As).
    A is the other element, fixed by the bulk: ``x mu_A + y mu_B = E_bulk``.
    ``formula_unit`` is ``(x, y)``.  Then

    ``gamma = [Eslab - (N_A/x) Ebulk - (N_B - (y/x) N_A) mu_B - C] / (n_surface A)``

    where ``C = reservoir_correction_ev`` collects any further reservoir
    terms evaluated at the same chemical potentials (for example the energy
    of pseudo-hydrogen on a passivated slab bottom).
    """

    if n_other <= 0 or n_variable < 0:
        raise ValueError("slab atom counts must be positive")
    if surface_area_angstrom2 <= 0:
        raise ValueError("surface_area_angstrom2 must be positive")
    if surfaces <= 0:
        raise ValueError("surfaces must be positive")
    other_per_formula, variable_per_formula = formula_unit
    if other_per_formula <= 0 or variable_per_formula <= 0:
        raise ValueError(f"formula_unit must contain positive atom counts, got {formula_unit!r}")
    mu_variable = variable_reference_energy_ev + delta_mu_ev
    excess_variable = n_variable - variable_per_formula * n_other / other_per_formula
    gamma = (
        slab_energy_ev
        - n_other * bulk_formula_energy_ev / other_per_formula
        - excess_variable * mu_variable
        - reservoir_correction_ev
    ) / (surfaces * surface_area_angstrom2)
    return SurfaceEnergyPoint(delta_mu_ev, gamma)


def surface_energy_oxide_equilibrium(
    *,
    slab_energy_ev: float,
    n_metal: int,
    n_oxygen: int,
    bulk_formula_energy_ev: float,
    oxygen_reference_energy_ev: float,
    delta_mu_oxygen_ev: float,
    surface_area_angstrom2: float,
    surfaces: int = 2,
    formula_unit: tuple[int, int] = (1, 2),
) -> SurfaceEnergyPoint:
    """Return gamma for an M_xO_y slab in equilibrium with its bulk oxide.

    ``formula_unit`` is ``(x, y)``, the metal and oxygen atoms in the formula
    unit whose energy is ``bulk_formula_energy_ev``: ``(1, 2)`` for MO2 (the
    default), ``(1, 1)`` for MO, ``(2, 3)`` for M2O3.  The convention is
    ``mu_O = E(O2)/2 + Delta mu_O``.  The expression is valid for a symmetric
    slab and exposes the oxygen-excess term explicitly:

    ``gamma = [Eslab - (N_M/x) Ebulk - (N_O - (y/x) N_M) mu_O] / (n_surface A)``.
    """

    if n_metal <= 0 or n_oxygen <= 0:
        raise ValueError("slab atom counts must be positive")
    return surface_energy_binary_equilibrium(
        slab_energy_ev=slab_energy_ev,
        n_other=n_metal,
        n_variable=n_oxygen,
        bulk_formula_energy_ev=bulk_formula_energy_ev,
        variable_reference_energy_ev=oxygen_reference_energy_ev / 2.0,
        delta_mu_ev=delta_mu_oxygen_ev,
        surface_area_angstrom2=surface_area_angstrom2,
        surfaces=surfaces,
        formula_unit=formula_unit,
    )


def surface_energy_elemental(
    *,
    slab_energy_ev: float,
    stoichiometry: Mapping[str, int],
    chemical_potentials_ev: Mapping[str, float],
    surface_area_angstrom2: float,
    surfaces: int = 2,
) -> float:
    """Return gamma using explicit elemental chemical potentials in eV/A2."""

    if surface_area_angstrom2 <= 0 or surfaces <= 0:
        raise ValueError("surface area and number of surfaces must be positive")
    missing = set(stoichiometry).difference(chemical_potentials_ev)
    if missing:
        raise ValueError(f"missing chemical potentials for {sorted(missing)}")
    reservoir = sum(
        count * chemical_potentials_ev[element]
        for element, count in stoichiometry.items()
    )
    return (slab_energy_ev - reservoir) / (surfaces * surface_area_angstrom2)


def stable_termination(
    points_by_label: Mapping[str, Iterable[SurfaceEnergyPoint]],
) -> list[tuple[float, str, float]]:
    """Choose the lowest-energy termination at each common chemical potential.

    This small deterministic helper deliberately does not interpolate: callers
    must supply a common grid, making numerical comparisons auditable.
    """

    grids = {label: list(points) for label, points in points_by_label.items()}
    if not grids:
        return []
    reference_grid = [point.delta_mu_ev for point in next(iter(grids.values()))]
    for label, points in grids.items():
        if [point.delta_mu_ev for point in points] != reference_grid:
            raise ValueError(f"chemical-potential grid for {label!r} differs from reference")
    return [
        (
            delta_mu,
            min(
                ((label, points[index].gamma_ev_per_angstrom2) for label, points in grids.items()),
                key=lambda item: item[1],
            )[0],
            min(points[index].gamma_ev_per_angstrom2 for points in grids.values()),
        )
        for index, delta_mu in enumerate(reference_grid)
    ]
