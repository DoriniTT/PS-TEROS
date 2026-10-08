"""Harmonic vibrational free energies for surface thermodynamics.

The DFT total energy E of a slab or bulk is the 0 K energy of fixed nuclei.
Its free energy adds the vibrational term of the harmonic approximation

``G(T) ~ F(T) = E + F_vib(T)``,  ``F_vib = sum_i [h nu_i / 2 + k_B T ln(1 - exp(-h nu_i / k_B T))]``

(``pV`` is negligible for solids).  The free energies of
:func:`solid_free_energy_ev` and :func:`molecule_reference_energy_ev` are
passed to the existing references and terminations in place of the total
energies; the phase-diagram code is unchanged and, at fixed temperature,
gamma stays linear in Delta mu.

Slabs whose central layers are frozen are treated with a partial Hessian:
only the free sites are displaced, and the frozen region is counted as bulk,
``F_vib(slab) = F_vib(free sites) + r F_vib(bulk cell)`` with ``r`` the number
of bulk cells in the frozen region (it must have the bulk stoichiometry).

Frequencies come from Gamma-point finite differences (see
:func:`psteros.build_vibrations_workgraph`); a bulk is calculated in a
Gamma-only supercell.  A future implementation may add the bulk vibrations
from phonopy on a q-point mesh.

Nothing here needs an AiiDA profile.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from math import isclose
from typing import Any, Literal, Mapping, Sequence

#: Boltzmann constant in eV/K.
KB_EV_PER_K = 8.617_333_262e-5
#: Energy of one wavenumber, h c (eV per cm^-1).
EV_PER_CM1 = 1.239_841_984e-4

# sqrt(eV / (A^2 amu)) in rad/s, and the matching wavenumber.
_ANGULAR_UNIT = math.sqrt(1.602_176_634e-19 / (1.0e-20 * 1.660_539_066_60e-27))
_CM1_PER_SQRT_EIGENVALUE = _ANGULAR_UNIT / (2.0 * math.pi * 2.997_924_58e10)

ImaginaryModes = Literal["raise", "drop"]


def _composition(value: Any, name: str) -> dict[str, int]:
    if hasattr(value, "get_el_amt_dict"):
        value = value.get_el_amt_dict()
    result: dict[str, int] = {}
    for element, amount in dict(value).items():
        if not isclose(float(amount), round(float(amount)), abs_tol=1e-8) or float(amount) < 0:
            raise ValueError(f"{name} must contain non-negative integer counts, got {element}={amount!r}")
        if round(float(amount)):
            result[str(element)] = int(round(float(amount)))
    return result


def _mode_free_energy_ev(energy_ev: float, temperature_k: float) -> float:
    if temperature_k == 0:
        return energy_ev / 2.0
    return energy_ev / 2.0 + KB_EV_PER_K * temperature_k * math.log1p(-math.exp(-energy_ev / (KB_EV_PER_K * temperature_k)))


def _mode_internal_energy_ev(energy_ev: float, temperature_k: float) -> float:
    if temperature_k == 0:
        return energy_ev / 2.0
    return energy_ev / 2.0 + energy_ev / math.expm1(energy_ev / (KB_EV_PER_K * temperature_k))


@dataclass(frozen=True)
class HarmonicVibrations:
    """Harmonic vibrational modes of one calculated structure.

    ``frequencies_cm1`` are the vibrational wavenumbers in cm^-1 without the
    zero modes (translations, and rotations of a molecule); an imaginary mode
    is given as a negative number.  ``composition`` is that of the structure
    whose energy is corrected (the input cell, not the supercell) and
    ``supercell_size`` the number of such cells in the calculated supercell:
    every energy below is per input cell.  ``frozen_composition`` counts the
    sites that were not displaced (the frozen centre of a slab).

    Imaginary modes raise an error (``imaginary_modes="raise"``): they mean
    the structure is not at a minimum, so its harmonic free energy is not
    defined.  ``imaginary_modes="drop"`` leaves them out instead.
    ``low_frequency_cutoff_cm1`` raises the real modes below it to the cutoff,
    which limits the entropy of very soft (and poorly converged) modes.
    """

    frequencies_cm1: tuple[float, ...]
    composition: Mapping[str, int] = field(default_factory=dict)
    frozen_composition: Mapping[str, int] = field(default_factory=dict)
    supercell_size: int = 1
    imaginary_modes: ImaginaryModes = "raise"
    low_frequency_cutoff_cm1: float | None = None

    def __post_init__(self) -> None:
        frequencies = tuple(float(value) for value in self.frequencies_cm1)
        if any(not math.isfinite(value) for value in frequencies):
            raise ValueError("frequencies_cm1 must be finite numbers")
        object.__setattr__(self, "frequencies_cm1", frequencies)
        object.__setattr__(self, "composition", _composition(self.composition, "composition"))
        object.__setattr__(self, "frozen_composition", _composition(self.frozen_composition, "frozen_composition"))
        if isinstance(self.supercell_size, bool) or not isinstance(self.supercell_size, int) or self.supercell_size <= 0:
            raise ValueError(f"supercell_size must be a positive integer, got {self.supercell_size!r}")
        if self.imaginary_modes not in ("raise", "drop"):
            raise ValueError(f"imaginary_modes must be 'raise' or 'drop', got {self.imaginary_modes!r}")
        if self.low_frequency_cutoff_cm1 is not None and self.low_frequency_cutoff_cm1 <= 0:
            raise ValueError("low_frequency_cutoff_cm1 must be positive")
        for element, count in self.frozen_composition.items():
            if self.composition and count > self.composition.get(element, 0):
                raise ValueError(f"frozen_composition has more {element} than composition")
        imaginary = self.imaginary_frequencies_cm1
        if imaginary and self.imaginary_modes == "raise":
            listing = ", ".join(f"{value:.1f}i" for value in imaginary)
            raise ValueError(
                f"{len(imaginary)} imaginary mode(s) ({listing} cm^-1): the structure is not at a minimum. "
                "Relax it further, or pass imaginary_modes='drop' to leave these modes out"
            )

    @property
    def imaginary_frequencies_cm1(self) -> tuple[float, ...]:
        """Magnitudes of the imaginary modes."""

        return tuple(-value for value in self.frequencies_cm1 if value < 0)

    @property
    def mode_energies_ev(self) -> tuple[float, ...]:
        """h nu of the modes entering the thermodynamics (real modes, after the cutoff)."""

        cutoff = self.low_frequency_cutoff_cm1 or 0.0
        return tuple(max(value, cutoff) * EV_PER_CM1 for value in self.frequencies_cm1 if value > 0)

    @property
    def zero_point_energy_ev(self) -> float:
        return sum(energy / 2.0 for energy in self.mode_energies_ev) / self.supercell_size

    def free_energy_ev(self, temperature_k: float) -> float:
        """Vibrational Helmholtz free energy F_vib(T) per input cell, zero-point energy included."""

        _check_temperature(temperature_k)
        return sum(_mode_free_energy_ev(e, temperature_k) for e in self.mode_energies_ev) / self.supercell_size

    def internal_energy_ev(self, temperature_k: float) -> float:
        """Vibrational internal energy U_vib(T) per input cell, zero-point energy included."""

        _check_temperature(temperature_k)
        return sum(_mode_internal_energy_ev(e, temperature_k) for e in self.mode_energies_ev) / self.supercell_size

    def entropy_ev_per_k(self, temperature_k: float) -> float:
        """Vibrational entropy S_vib(T) per input cell, in eV/K."""

        _check_temperature(temperature_k)
        if temperature_k == 0:
            return 0.0
        return (self.internal_energy_ev(temperature_k) - self.free_energy_ev(temperature_k)) / temperature_k

    def to_dict(self) -> dict[str, Any]:
        """Plain mapping (stored in AiiDA as a ``Dict``); :meth:`from_dict` reads it back."""

        return {
            "frequencies_cm1": list(self.frequencies_cm1),
            "composition": dict(self.composition),
            "frozen_composition": dict(self.frozen_composition),
            "supercell_size": self.supercell_size,
        }

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        *,
        imaginary_modes: ImaginaryModes = "raise",
        low_frequency_cutoff_cm1: float | None = None,
    ) -> "HarmonicVibrations":
        """From :meth:`to_dict`, an AiiDA ``Dict`` or the ``<label>_vibrations`` output of a graph."""

        if hasattr(data, "get_dict"):
            data = data.get_dict()
        return cls(
            frequencies_cm1=tuple(data["frequencies_cm1"]),
            composition=data.get("composition", {}),
            frozen_composition=data.get("frozen_composition", {}),
            supercell_size=int(data.get("supercell_size", 1)),
            imaginary_modes=imaginary_modes,
            low_frequency_cutoff_cm1=low_frequency_cutoff_cm1,
        )


def _check_temperature(temperature_k: float) -> None:
    if not math.isfinite(temperature_k) or temperature_k < 0:
        raise ValueError(f"temperature_k must be a non-negative number of kelvin, got {temperature_k!r}")


def solid_free_energy_ev(
    energy_ev: float,
    vibrations: HarmonicVibrations,
    temperature_k: float,
    *,
    bulk: HarmonicVibrations | None = None,
) -> float:
    """Free energy ``E + F_vib(T)`` of a slab or bulk, to use in place of its total energy.

    ``energy_ev`` is the DFT energy of the structure whose vibrations are
    ``vibrations``.  When sites were frozen (``vibrations.frozen_composition``),
    the frozen region is counted as bulk: ``bulk`` gives the bulk vibrations,
    and the frozen region must contain the bulk elements in the bulk ratio.
    """

    free_energy = energy_ev + vibrations.free_energy_ev(temperature_k)
    frozen = vibrations.frozen_composition
    if not frozen:
        return free_energy
    if bulk is None:
        raise ValueError(
            f"the structure has frozen sites {frozen}: pass bulk=<bulk HarmonicVibrations> "
            "to count the frozen region as bulk, or displace every site"
        )
    if bulk.frozen_composition:
        raise ValueError("the bulk reference must have every site displaced (no frozen_composition)")
    if not bulk.composition:
        raise ValueError("the bulk reference needs its composition")
    if set(frozen) != set(bulk.composition):
        raise ValueError(
            f"the frozen region {frozen} does not have the elements of the bulk {bulk.composition}; "
            "it cannot be counted as bulk"
        )
    ratios = {element: frozen[element] / bulk.composition[element] for element in frozen}
    cells = next(iter(ratios.values()))
    if any(not isclose(ratio, cells, rel_tol=1e-9) for ratio in ratios.values()):
        raise ValueError(
            f"the frozen region {frozen} does not have the bulk stoichiometry {bulk.composition}; "
            "freeze whole bulk formula units"
        )
    return free_energy + cells * bulk.free_energy_ev(temperature_k)


def molecule_reference_energy_ev(energy_ev: float, vibrations: HarmonicVibrations) -> float:
    """Energy of a gas molecule with its zero-point energy, ``E + ZPE``.

    This is the reference that fixes ``Delta mu = 0`` (for example
    ``oxygen_molecule_energy_ev`` for E(O2)): the thermal part of the gas
    belongs to ``Delta mu(T, p)`` itself, so only the zero-point energy is
    added here.
    """

    return energy_ev + vibrations.zero_point_energy_ev


def harmonic_vibrations_from_forces(
    *,
    masses_amu: Sequence[float],
    displaced_sites: Sequence[int],
    displacement_angstrom: float,
    forces_plus: Sequence[Any],
    forces_minus: Sequence[Any],
    composition: Mapping[str, int] | None = None,
    frozen_composition: Mapping[str, int] | None = None,
    zero_modes: int = 0,
    supercell_size: int = 1,
    imaginary_modes: ImaginaryModes = "raise",
    low_frequency_cutoff_cm1: float | None = None,
) -> HarmonicVibrations:
    """Gamma-point modes from central finite differences of the forces.

    For every displaced site ``i`` (in ``displaced_sites`` order) and every
    Cartesian direction ``a``, ``forces_plus[3 k + a]`` and
    ``forces_minus[3 k + a]`` are the forces (eV/A, one row per site of the
    structure) after moving site ``i`` by ``+d`` and ``-d`` along ``a``.  The
    Hessian of the displaced sites is
    ``H[(i,a),(j,b)] = -[F_jb(+d_ia) - F_jb(-d_ia)] / 2d``, symmetrised and
    mass-weighted.  ``zero_modes`` modes of smallest magnitude are removed: 3
    for a periodic cell with every site displaced, 5 for a linear molecule,
    6 for another molecule, 0 for a slab with frozen sites.
    """

    import numpy as np

    sites = list(displaced_sites)
    if not sites:
        raise ValueError("displaced_sites must not be empty")
    if len(set(sites)) != len(sites):
        raise ValueError("displaced_sites must not repeat a site")
    if displacement_angstrom <= 0:
        raise ValueError("displacement_angstrom must be positive")
    masses = np.asarray(masses_amu, dtype=float)
    plus = np.asarray(forces_plus, dtype=float)
    minus = np.asarray(forces_minus, dtype=float)
    expected = (3 * len(sites), len(masses), 3)
    if plus.shape != expected or minus.shape != expected:
        raise ValueError(f"forces_plus and forces_minus must have shape {expected}, got {plus.shape} and {minus.shape}")
    if any(site < 0 or site >= len(masses) for site in sites):
        raise ValueError(f"displaced_sites must be site indices in [0, {len(masses)})")
    if np.any(masses <= 0):
        raise ValueError("masses_amu must be positive")
    if not 0 <= zero_modes < 3 * len(sites):
        raise ValueError(f"zero_modes must be in [0, {3 * len(sites)})")

    hessian = -(plus - minus)[:, sites, :].reshape(3 * len(sites), 3 * len(sites)) / (2.0 * displacement_angstrom)
    hessian = (hessian + hessian.T) / 2.0
    weights = 1.0 / np.sqrt(np.repeat(masses[sites], 3))
    eigenvalues = np.linalg.eigvalsh(hessian * np.outer(weights, weights))
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * _CM1_PER_SQRT_EIGENVALUE
    keep = sorted(np.argsort(np.abs(frequencies))[zero_modes:])
    return HarmonicVibrations(
        frequencies_cm1=tuple(float(frequencies[index]) for index in keep),
        composition=composition or {},
        frozen_composition=frozen_composition or {},
        supercell_size=supercell_size,
        imaginary_modes=imaginary_modes,
        low_frequency_cutoff_cm1=low_frequency_cutoff_cm1,
    )
