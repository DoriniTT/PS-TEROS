"""Free energies of reference phases from DFT energies and vibrational frequencies.

Two models, both independent of AiiDA:

* :class:`IdealGasMolecule`: ideal gas, rigid rotor and harmonic oscillator
  for molecules such as O2, H2 or H2O.
* :class:`HarmonicSolid`: harmonic lattice vibrations of a bulk reference,
  from the modes of a supercell.

Both return a :class:`FreeEnergy` that keeps every term apart (DFT energy,
zero-point energy, thermal enthalpy, entropy, explicit correction), so each
number entering a phase diagram can be audited.  Frequencies are in cm^-1;
an imaginary mode is given as a negative number, as
:func:`parse_vasp_frequencies_cm1` reads it from a VASP OUTCAR.

:func:`delta_mu_oxygen_ev` places a temperature and O2 pressure on the psteros
oxygen axis ``mu_O = E(O2)/2 + Delta mu_O``, where ``E(O2)`` is the bare DFT
energy.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from math import exp, log, pi, sqrt
from typing import Any, Literal, Mapping, Sequence

# CODATA 2018.
BOLTZMANN_EV_PER_K = 8.617_333_262e-5
PLANCK_EV_S = 4.135_667_696e-15
PLANCK_J_S = 6.626_070_15e-34
BOLTZMANN_J_PER_K = 1.380_649e-23
SPEED_OF_LIGHT_CM_PER_S = 2.997_924_58e10
AMU_KG = 1.660_539_066_60e-27
ANGSTROM_M = 1.0e-10
CM1_TO_EV = PLANCK_EV_S * SPEED_OF_LIGHT_CM_PER_S
BAR_PA = 1.0e5
STANDARD_PRESSURE_BAR = 1.0

Geometry = Literal["monatomic", "linear", "nonlinear"]

# A principal moment below this (amu A^2) is zero: the molecule is linear.
_LINEAR_MOMENT_TOLERANCE = 1.0e-3


@dataclass(frozen=True)
class FreeEnergy:
    """Terms of a free energy at one temperature and pressure, in eV.

    ``free_energy_ev = electronic_energy_ev + zero_point_energy_ev
    + thermal_enthalpy_ev - temperature_k * entropy_ev_per_k + correction_ev``.

    For a gas, ``thermal_enthalpy_ev`` is ``H(T) - H(0)`` (translation,
    rotation, vibration and ``pV = kT``) and ``entropy_ev_per_k`` includes the
    pressure term.  For a solid, it is the vibrational ``U(T) - U(0)`` and the
    ``pV`` term is neglected.
    """

    temperature_k: float
    pressure_bar: float | None
    electronic_energy_ev: float
    zero_point_energy_ev: float
    thermal_enthalpy_ev: float
    entropy_ev_per_k: float
    correction_ev: float = 0.0

    @property
    def entropy_term_ev(self) -> float:
        """``-T S`` in eV."""

        return -self.temperature_k * self.entropy_ev_per_k

    @property
    def enthalpy_ev(self) -> float:
        return (
            self.electronic_energy_ev
            + self.zero_point_energy_ev
            + self.thermal_enthalpy_ev
            + self.correction_ev
        )

    @property
    def free_energy_ev(self) -> float:
        return self.enthalpy_ev + self.entropy_term_ev

    @property
    def vibrational_free_energy_ev(self) -> float:
        """Everything except the DFT energy and the correction: ``G - E - correction``."""

        return self.zero_point_energy_ev + self.thermal_enthalpy_ev + self.entropy_term_ev

    def as_dict(self) -> dict[str, float | None]:
        """All terms and totals, for tables and CSV export."""

        return {
            "temperature_K": self.temperature_k,
            "pressure_bar": self.pressure_bar,
            "electronic_energy_eV": self.electronic_energy_ev,
            "zero_point_energy_eV": self.zero_point_energy_ev,
            "thermal_enthalpy_eV": self.thermal_enthalpy_ev,
            "entropy_eV_per_K": self.entropy_ev_per_k,
            "minus_TS_eV": self.entropy_term_ev,
            "correction_eV": self.correction_ev,
            "enthalpy_eV": self.enthalpy_ev,
            "free_energy_eV": self.free_energy_ev,
        }


def _check_temperature(temperature_k: float) -> None:
    if temperature_k < 0:
        raise ValueError(f"temperature_k must not be negative, got {temperature_k!r}")


def _real_modes(
    frequencies_cm1: Sequence[float],
    *,
    expected: int,
    discard: int,
    imaginary_tolerance_cm1: float,
    owner: str,
) -> tuple[float, ...]:
    """Drop the ``discard`` softest modes and return the ``expected`` remaining ones.

    A full list of ``expected + discard`` modes (as VASP prints it) loses its
    ``discard`` modes of smallest magnitude, the translations and rotations
    (molecule) or acoustic modes at Gamma (solid).  A list of exactly
    ``expected`` modes is taken as already filtered.  An imaginary mode left
    after that is an error unless its magnitude is below
    ``imaginary_tolerance_cm1``, in which case it is dropped.
    """

    modes = [float(value) for value in frequencies_cm1]
    if len(modes) == expected + discard:
        modes = sorted(modes, key=abs)[discard:]
    elif len(modes) != expected:
        raise ValueError(
            f"{owner}: expected {expected + discard} modes (all, as VASP prints them) or "
            f"{expected} vibrational modes, got {len(modes)}"
        )
    imaginary = [value for value in modes if value < 0]
    too_large = [value for value in imaginary if -value >= imaginary_tolerance_cm1]
    if too_large:
        raise ValueError(
            f"{owner}: imaginary modes {sorted(too_large)} cm^-1 remain; the structure is "
            "not at a minimum (tighten the relaxation) or, for numerical noise, raise "
            "imaginary_tolerance_cm1 to drop them"
        )
    return tuple(value for value in modes if value > 0)


def _harmonic_terms(modes_cm1: Sequence[float], temperature_k: float) -> tuple[float, float, float]:
    """Zero-point energy (eV), U(T) - U(0) (eV) and entropy (eV/K) of harmonic modes."""

    energies = [CM1_TO_EV * value for value in modes_cm1]
    zero_point = 0.5 * sum(energies)
    if temperature_k == 0:
        return zero_point, 0.0, 0.0
    kt = BOLTZMANN_EV_PER_K * temperature_k
    thermal = 0.0
    entropy = 0.0
    for energy in energies:
        x = energy / kt
        if x > 700:  # frozen mode: no thermal energy or entropy, avoids overflow
            continue
        thermal += energy / (exp(x) - 1.0)
        entropy += BOLTZMANN_EV_PER_K * (x / (exp(x) - 1.0) - log(1.0 - exp(-x)))
    return zero_point, thermal, entropy


def _cartesian_coordinates(structure: Any) -> tuple[list[str], list[tuple[float, float, float]]]:
    """Symbols and unwrapped Cartesian coordinates (A) of a molecule.

    A periodic cell is unwrapped to the minimum image of each atom relative to
    the first one, so a molecule split by the cell boundary stays whole.
    """

    if hasattr(structure, "get_pymatgen_structure"):
        structure = structure.get_pymatgen_structure()
    symbols = [site.specie.symbol for site in structure]
    if not hasattr(structure, "lattice"):
        return symbols, [tuple(float(x) for x in site.coords) for site in structure]
    origin = structure[0].frac_coords
    coordinates = []
    for site in structure:
        shift = [value - round(value) for value in site.frac_coords - origin]
        cartesian = structure.lattice.get_cartesian_coords([o + s for o, s in zip(origin, shift)])
        coordinates.append(tuple(float(x) for x in cartesian))
    return symbols, coordinates


def _principal_moments(
    masses_amu: Sequence[float], positions_angstrom: Sequence[Sequence[float]]
) -> tuple[float, float, float]:
    """Principal moments of inertia in amu A^2, smallest first."""

    import numpy as np

    masses = np.asarray(masses_amu, dtype=float)
    positions = np.asarray(positions_angstrom, dtype=float)
    centred = positions - (masses[:, None] * positions).sum(axis=0) / masses.sum()
    tensor = np.zeros((3, 3))
    for mass, (x, y, z) in zip(masses, centred):
        tensor += mass * np.array(
            [
                [y * y + z * z, -x * y, -x * z],
                [-x * y, x * x + z * z, -y * z],
                [-x * z, -y * z, x * x + y * y],
            ]
        )
    return tuple(float(value) for value in sorted(np.linalg.eigvalsh(tensor)))


@dataclass(frozen=True)
class IdealGasMolecule:
    """A gas-phase reference: ideal gas, rigid rotor, harmonic oscillator.

    ``frequencies_cm1`` are either all ``3N`` modes as printed by VASP
    (translations and rotations are dropped as the softest ``3N - n_vib``) or
    exactly the ``n_vib`` vibrational modes (``3N - 5`` linear, ``3N - 6``
    nonlinear).  ``symmetry_number`` is the rotational symmetry number (2 for
    O2, H2 and H2O, 1 for CO) and ``spin`` the total electron spin S (1 for
    triplet O2, 0 for a closed shell); the electronic entropy is
    ``k ln(2S + 1)``.  ``correction_ev`` is an explicit, opt-in correction of
    the energy (for example a fitted O2 binding correction); it is reported
    as its own term.  Use :meth:`from_structure` to read masses and geometry
    from a structure.
    """

    electronic_energy_ev: float
    frequencies_cm1: tuple[float, ...]
    masses_amu: tuple[float, ...]
    positions_angstrom: tuple[tuple[float, float, float], ...]
    symmetry_number: int
    spin: float
    geometry: Geometry | None = None
    correction_ev: float = 0.0
    imaginary_tolerance_cm1: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "frequencies_cm1", tuple(float(v) for v in self.frequencies_cm1))
        object.__setattr__(self, "masses_amu", tuple(float(v) for v in self.masses_amu))
        object.__setattr__(
            self,
            "positions_angstrom",
            tuple(tuple(float(x) for x in position) for position in self.positions_angstrom),
        )
        if not self.masses_amu or len(self.masses_amu) != len(self.positions_angstrom):
            raise ValueError("masses_amu and positions_angstrom must describe the same atoms")
        if any(mass <= 0 for mass in self.masses_amu):
            raise ValueError("masses_amu must be positive")
        if isinstance(self.symmetry_number, bool) or not isinstance(self.symmetry_number, int) or self.symmetry_number <= 0:
            raise ValueError(f"symmetry_number must be a positive integer, got {self.symmetry_number!r}")
        if self.spin < 0 or not float(2 * self.spin).is_integer():
            raise ValueError(f"spin must be a non-negative multiple of 1/2, got {self.spin!r}")
        if self.imaginary_tolerance_cm1 < 0:
            raise ValueError("imaginary_tolerance_cm1 must not be negative")
        detected = self._detected_geometry()
        if self.geometry is None:
            object.__setattr__(self, "geometry", detected)
        elif self.geometry not in ("monatomic", "linear", "nonlinear"):
            raise ValueError(
                f"geometry must be 'monatomic', 'linear' or 'nonlinear', got {self.geometry!r}"
            )
        elif self.geometry != detected:
            raise ValueError(f"geometry {self.geometry!r} contradicts the positions, which are {detected}")
        # Validate the modes now, so a bad molecule fails when it is defined.
        self.vibrational_frequencies_cm1

    @classmethod
    def from_structure(
        cls,
        structure: Any,
        *,
        electronic_energy_ev: float,
        frequencies_cm1: Sequence[float],
        symmetry_number: int,
        spin: float,
        geometry: Geometry | None = None,
        correction_ev: float = 0.0,
        imaginary_tolerance_cm1: float = 0.0,
    ) -> "IdealGasMolecule":
        """Read masses and positions from a pymatgen ``Structure``/``Molecule`` or ``StructureData``."""

        from pymatgen.core import Element

        symbols, positions = _cartesian_coordinates(structure)
        return cls(
            electronic_energy_ev=electronic_energy_ev,
            frequencies_cm1=tuple(frequencies_cm1),
            masses_amu=tuple(float(Element(symbol).atomic_mass) for symbol in symbols),
            positions_angstrom=tuple(positions),
            symmetry_number=symmetry_number,
            spin=spin,
            geometry=geometry,
            correction_ev=correction_ev,
            imaginary_tolerance_cm1=imaginary_tolerance_cm1,
        )

    def _detected_geometry(self) -> Geometry:
        if len(self.masses_amu) == 1:
            return "monatomic"
        smallest = _principal_moments(self.masses_amu, self.positions_angstrom)[0]
        return "linear" if smallest < _LINEAR_MOMENT_TOLERANCE else "nonlinear"

    @property
    def number_of_atoms(self) -> int:
        return len(self.masses_amu)

    @property
    def vibrational_frequencies_cm1(self) -> tuple[float, ...]:
        """The real vibrational modes used in the thermochemistry."""

        n_atoms = self.number_of_atoms
        discard = {"monatomic": 3, "linear": 5, "nonlinear": 6}[self.geometry]
        return _real_modes(
            self.frequencies_cm1,
            expected=3 * n_atoms - discard,
            discard=discard,
            imaginary_tolerance_cm1=self.imaginary_tolerance_cm1,
            owner=f"{self.geometry} molecule of {n_atoms} atoms",
        )

    def free_energy(
        self, temperature_k: float, pressure_bar: float = STANDARD_PRESSURE_BAR
    ) -> FreeEnergy:
        """Gibbs free energy terms of one molecule at ``temperature_k`` and ``pressure_bar``."""

        _check_temperature(temperature_k)
        if pressure_bar <= 0:
            raise ValueError(f"pressure_bar must be positive, got {pressure_bar!r}")
        zero_point, vib_thermal, vib_entropy = _harmonic_terms(
            self.vibrational_frequencies_cm1, temperature_k
        )
        if temperature_k == 0:
            return FreeEnergy(
                temperature_k, pressure_bar, self.electronic_energy_ev, zero_point, 0.0, 0.0,
                self.correction_ev,
            )
        kt = BOLTZMANN_EV_PER_K * temperature_k
        rotational_dof = {"monatomic": 0, "linear": 2, "nonlinear": 3}[self.geometry]
        thermal = 2.5 * kt + 0.5 * rotational_dof * kt + vib_thermal
        entropy = (
            self._translational_entropy(temperature_k, pressure_bar)
            + self._rotational_entropy(temperature_k)
            + vib_entropy
            + BOLTZMANN_EV_PER_K * log(2 * self.spin + 1)
        )
        return FreeEnergy(
            temperature_k, pressure_bar, self.electronic_energy_ev, zero_point, thermal, entropy,
            self.correction_ev,
        )

    def _translational_entropy(self, temperature_k: float, pressure_bar: float) -> float:
        mass = sum(self.masses_amu) * AMU_KG
        kt_j = BOLTZMANN_J_PER_K * temperature_k
        volume_per_molecule = kt_j / (pressure_bar * BAR_PA)
        quantum = (2 * pi * mass * kt_j / PLANCK_J_S**2) ** 1.5
        return BOLTZMANN_EV_PER_K * (log(quantum * volume_per_molecule) + 2.5)

    def _rotational_entropy(self, temperature_k: float) -> float:
        if self.geometry == "monatomic":
            return 0.0
        to_si = AMU_KG * ANGSTROM_M**2
        moments = _principal_moments(self.masses_amu, self.positions_angstrom)
        kt_j = BOLTZMANN_J_PER_K * temperature_k
        if self.geometry == "linear":
            inertia = moments[2] * to_si
            partition = 8 * pi**2 * inertia * kt_j / (self.symmetry_number * PLANCK_J_S**2)
            return BOLTZMANN_EV_PER_K * (log(partition) + 1.0)
        product = moments[0] * moments[1] * moments[2] * to_si**3
        partition = sqrt(pi * product) / self.symmetry_number * (8 * pi**2 * kt_j / PLANCK_J_S**2) ** 1.5
        return BOLTZMANN_EV_PER_K * (log(partition) + 1.5)


@dataclass(frozen=True)
class HarmonicSolid:
    """A crystalline reference with harmonic lattice vibrations.

    ``electronic_energy_ev`` is the DFT energy of the reference cell, which
    holds ``atoms_in_cell`` atoms (the cell whose energy enters the phase
    diagram).  ``frequencies_cm1`` are the modes of a supercell of the same
    crystal with ``atoms_in_supercell`` atoms (a supercell, because the Gamma
    point of a small cell misses most of the phonon spectrum; for an
    elemental metal it has only the zero acoustic modes).  Either all
    ``3 * atoms_in_supercell`` modes are given, and the three acoustic modes
    (softest) are dropped, or the ``3 * atoms_in_supercell - 3`` others.  The
    vibrational terms are scaled to the reference cell by
    ``atoms_in_cell / atoms_in_supercell``.  The ``pV`` term is neglected.
    """

    electronic_energy_ev: float
    atoms_in_cell: int
    frequencies_cm1: tuple[float, ...]
    atoms_in_supercell: int
    correction_ev: float = 0.0
    imaginary_tolerance_cm1: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "frequencies_cm1", tuple(float(v) for v in self.frequencies_cm1))
        for name in ("atoms_in_cell", "atoms_in_supercell"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer, got {value!r}")
        if self.imaginary_tolerance_cm1 < 0:
            raise ValueError("imaginary_tolerance_cm1 must not be negative")
        self.vibrational_frequencies_cm1

    @property
    def vibrational_frequencies_cm1(self) -> tuple[float, ...]:
        """The real supercell modes used in the thermochemistry."""

        return _real_modes(
            self.frequencies_cm1,
            expected=3 * self.atoms_in_supercell - 3,
            discard=3,
            imaginary_tolerance_cm1=self.imaginary_tolerance_cm1,
            owner=f"solid supercell of {self.atoms_in_supercell} atoms",
        )

    def free_energy(self, temperature_k: float) -> FreeEnergy:
        """Harmonic Helmholtz free energy terms of the reference cell at ``temperature_k``."""

        _check_temperature(temperature_k)
        zero_point, thermal, entropy = _harmonic_terms(self.vibrational_frequencies_cm1, temperature_k)
        scale = self.atoms_in_cell / self.atoms_in_supercell
        return FreeEnergy(
            temperature_k,
            None,
            self.electronic_energy_ev,
            zero_point * scale,
            thermal * scale,
            entropy * scale,
            self.correction_ev,
        )


def free_energies(
    systems: Mapping[str, "IdealGasMolecule | HarmonicSolid"],
    temperature_k: float,
    pressure_bar: float = STANDARD_PRESSURE_BAR,
) -> dict[str, FreeEnergy]:
    """Free energy terms of labelled references at one temperature (and gas pressure).

    ``pressure_bar`` applies to every gas; solids do not depend on it.
    """

    return {
        label: system.free_energy(temperature_k, pressure_bar)
        if isinstance(system, IdealGasMolecule)
        else system.free_energy(temperature_k)
        for label, system in systems.items()
    }


def delta_mu_oxygen_ev(
    oxygen: IdealGasMolecule,
    temperature_k: float,
    pressure_bar: float = STANDARD_PRESSURE_BAR,
    *,
    include_zero_point: bool = True,
) -> float:
    """Delta mu_O(T, p) of O2 gas on the psteros axis ``mu_O = E(O2)/2 + Delta mu_O``.

    ``Delta mu_O = [G_O2(T, p) - E_DFT(O2)] / 2``: half of the zero-point,
    thermal, entropy, pressure and correction terms of ``oxygen``.  With
    ``include_zero_point=False`` the zero-point energy and the correction are
    left out, which is the tabulated ``Delta mu_O`` of Reuter and Scheffler,
    Phys. Rev. B 65, 035406 (2001).
    """

    terms = oxygen.free_energy(temperature_k, pressure_bar)
    shift = terms.thermal_enthalpy_ev + terms.entropy_term_ev
    if include_zero_point:
        shift += terms.zero_point_energy_ev + terms.correction_ev
    return shift / 2.0


def oxygen_pressure_bar(
    oxygen: IdealGasMolecule,
    temperature_k: float,
    delta_mu_oxygen: float,
    *,
    include_zero_point: bool = True,
) -> float:
    """O2 pressure (bar) at which O2 gas at ``temperature_k`` has ``delta_mu_oxygen`` (eV).

    The inverse of :func:`delta_mu_oxygen_ev`, using
    ``Delta mu_O(T, p) = Delta mu_O(T, 1 bar) + kT ln(p / 1 bar) / 2``.
    """

    if temperature_k <= 0:
        raise ValueError(f"temperature_k must be positive, got {temperature_k!r}")
    standard = delta_mu_oxygen_ev(
        oxygen, temperature_k, STANDARD_PRESSURE_BAR, include_zero_point=include_zero_point
    )
    return STANDARD_PRESSURE_BAR * exp(
        2.0 * (delta_mu_oxygen - standard) / (BOLTZMANN_EV_PER_K * temperature_k)
    )


_OUTCAR_MODE = re.compile(
    r"^\s*\d+\s+f(/i)?\s*=\s*[-\d.]+\s+THz\s+[-\d.]+\s+2PiTHz\s+([-\d.]+)\s+cm-1",
    re.MULTILINE,
)


def parse_vasp_frequencies_cm1(outcar: str) -> tuple[float, ...]:
    """Frequencies (cm^-1) of the last dynamical matrix in a VASP OUTCAR.

    Imaginary modes (``f/i=``) are returned as negative numbers.  Works for
    finite differences (``IBRION = 5, 6``) and DFPT (``IBRION = 7, 8``).
    """

    marker = "Eigenvectors and eigenvalues of the dynamical matrix"
    position = outcar.rfind(marker)
    if position < 0:
        raise ValueError("no dynamical matrix in the OUTCAR; was it run with IBRION = 5-8?")
    block = outcar[position:]
    # VASP then prints the same modes again, divided by sqrt(mass); keep one set.
    end = block.find("Eigenvectors after division by SQRT(mass)")
    if end > 0:
        block = block[:end]
    modes = [(-1.0 if imaginary else 1.0) * float(value) for imaginary, value in _OUTCAR_MODE.findall(block)]
    if not modes:
        raise ValueError("the OUTCAR dynamical matrix block contains no frequencies")
    return tuple(modes)
