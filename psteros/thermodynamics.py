"""Surface free-energy analysis independent of a workflow engine."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping


EV_PER_ANGSTROM2_TO_J_PER_M2 = 16.021_766_34


@dataclass(frozen=True)
class SurfaceEnergyPoint:
    """A surface energy at one oxygen chemical-potential offset."""

    delta_mu_oxygen_ev: float
    gamma_ev_per_angstrom2: float

    @property
    def gamma_j_per_m2(self) -> float:
        return self.gamma_ev_per_angstrom2 * EV_PER_ANGSTROM2_TO_J_PER_M2


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
) -> SurfaceEnergyPoint:
    """Return gamma for an MO2 slab in equilibrium with its bulk oxide.

    The convention is ``mu_O = E(O2)/2 + Delta mu_O``.  The expression is
    valid for a symmetric slab and exposes the oxygen-excess term explicitly:

    ``gamma = [Eslab - N_M Ebulk - (N_O - 2 N_M) mu_O] / (n_surface A)``.
    """

    if n_metal <= 0 or n_oxygen <= 0:
        raise ValueError("slab atom counts must be positive")
    if surface_area_angstrom2 <= 0:
        raise ValueError("surface_area_angstrom2 must be positive")
    if surfaces <= 0:
        raise ValueError("surfaces must be positive")
    mu_oxygen = oxygen_reference_energy_ev / 2.0 + delta_mu_oxygen_ev
    excess_oxygen = n_oxygen - 2 * n_metal
    gamma = (
        slab_energy_ev - n_metal * bulk_formula_energy_ev - excess_oxygen * mu_oxygen
    ) / (surfaces * surface_area_angstrom2)
    return SurfaceEnergyPoint(delta_mu_oxygen_ev, gamma)


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
    reference_grid = [point.delta_mu_oxygen_ev for point in next(iter(grids.values()))]
    for label, points in grids.items():
        if [point.delta_mu_oxygen_ev for point in points] != reference_grid:
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

