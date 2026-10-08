"""Tier-1 tests for the reference thermochemistry (ideal-gas molecules, harmonic solids)."""

from __future__ import annotations

from math import cos, radians, sin

import pytest

import psteros
from psteros.thermochemistry import (
    BOLTZMANN_EV_PER_K,
    CM1_TO_EV,
    HarmonicSolid,
    IdealGasMolecule,
    delta_mu_oxygen_ev,
    oxygen_pressure_bar,
    parse_vasp_frequencies_cm1,
)

# eV per molecule and K -> J/(mol K)
EV_PER_K_TO_J_PER_MOL_K = 96_485.332_12
ATM_BAR = 1.013_25

# Experimental O2: r = 1.2075 A, omega = 1580 cm^-1, triplet, sigma = 2.
O2_BOND = 1.2075
O2_MODE = 1580.0


def oxygen(frequencies=(O2_MODE,), **kwargs) -> IdealGasMolecule:
    return IdealGasMolecule(
        electronic_energy_ev=-9.86,
        frequencies_cm1=frequencies,
        masses_amu=(15.999, 15.999),
        positions_angstrom=((0.0, 0.0, 0.0), (0.0, 0.0, O2_BOND)),
        symmetry_number=2,
        spin=1.0,
        **kwargs,
    )


def water() -> IdealGasMolecule:
    # Experimental geometry (r = 0.9572 A, HOH = 104.52 deg) and fundamentals.
    half = radians(104.52) / 2
    return IdealGasMolecule(
        electronic_energy_ev=-14.22,
        frequencies_cm1=(1595.0, 3657.0, 3756.0),
        masses_amu=(15.999, 1.008, 1.008),
        positions_angstrom=(
            (0.0, 0.0, 0.0),
            (0.9572 * sin(half), 0.0, 0.9572 * cos(half)),
            (-0.9572 * sin(half), 0.0, 0.9572 * cos(half)),
        ),
        symmetry_number=2,
        spin=0.0,
    )


def test_free_energy_terms_add_up() -> None:
    terms = oxygen(correction_ev=-0.1).free_energy(500.0, 0.2)
    assert terms.free_energy_ev == pytest.approx(
        terms.electronic_energy_ev
        + terms.zero_point_energy_ev
        + terms.thermal_enthalpy_ev
        - 500.0 * terms.entropy_ev_per_k
        + terms.correction_ev
    )
    assert terms.vibrational_free_energy_ev == pytest.approx(
        terms.free_energy_ev - terms.electronic_energy_ev - terms.correction_ev
    )
    assert terms.as_dict()["free_energy_eV"] == pytest.approx(terms.free_energy_ev)
    assert terms.as_dict()["minus_TS_eV"] == pytest.approx(-500.0 * terms.entropy_ev_per_k)


def test_zero_kelvin_keeps_only_zero_point_energy() -> None:
    terms = oxygen().free_energy(0.0)
    assert terms.zero_point_energy_ev == pytest.approx(0.5 * O2_MODE * CM1_TO_EV)
    assert terms.thermal_enthalpy_ev == 0.0
    assert terms.entropy_ev_per_k == 0.0
    assert terms.free_energy_ev == pytest.approx(-9.86 + 0.5 * O2_MODE * CM1_TO_EV)


def test_standard_entropies_match_nist() -> None:
    # NIST-JANAF S(298.15 K, 1 bar): O2 205.15, H2O(g) 188.84 J/(mol K).
    s_o2 = oxygen().free_energy(298.15, 1.0).entropy_ev_per_k * EV_PER_K_TO_J_PER_MOL_K
    s_h2o = water().free_energy(298.15, 1.0).entropy_ev_per_k * EV_PER_K_TO_J_PER_MOL_K
    assert s_o2 == pytest.approx(205.15, abs=0.5)
    assert s_h2o == pytest.approx(188.84, abs=0.5)


@pytest.mark.parametrize(
    "temperature, expected",
    # Reuter and Scheffler, Phys. Rev. B 65, 035406 (2001), Table I, p = 1 atm.
    [(300.0, -0.27), (600.0, -0.61), (1000.0, -1.10)],
)
def test_delta_mu_oxygen_matches_reuter_scheffler(temperature, expected) -> None:
    value = delta_mu_oxygen_ev(oxygen(), temperature, ATM_BAR, include_zero_point=False)
    assert value == pytest.approx(expected, abs=0.02)


def test_delta_mu_oxygen_includes_half_the_zero_point_and_pressure() -> None:
    o2 = oxygen(correction_ev=-0.2)
    bare = delta_mu_oxygen_ev(o2, 700.0, include_zero_point=False)
    full = delta_mu_oxygen_ev(o2, 700.0)
    assert full - bare == pytest.approx(0.5 * (0.5 * O2_MODE * CM1_TO_EV - 0.2))
    shifted = delta_mu_oxygen_ev(o2, 700.0, 1.0e-10)
    assert shifted - full == pytest.approx(0.5 * BOLTZMANN_EV_PER_K * 700.0 * -23.025850929940457)


def test_oxygen_pressure_inverts_delta_mu() -> None:
    o2 = oxygen()
    value = delta_mu_oxygen_ev(o2, 800.0, 3.0e-6)
    assert oxygen_pressure_bar(o2, 800.0, value) == pytest.approx(3.0e-6)


def test_vasp_mode_list_drops_translations_and_rotations() -> None:
    all_modes = (O2_MODE, 12.0, -8.0, 3.0, -1.5, 0.4)
    assert oxygen(all_modes).vibrational_frequencies_cm1 == (O2_MODE,)
    assert oxygen(all_modes).free_energy(300.0).free_energy_ev == pytest.approx(
        oxygen().free_energy(300.0).free_energy_ev
    )


def test_imaginary_vibrational_mode_is_rejected_unless_tolerated() -> None:
    modes = (1595.0, -40.0, 3756.0)
    with pytest.raises(ValueError, match="imaginary modes"):
        IdealGasMolecule(
            -14.22, modes, water().masses_amu, water().positions_angstrom, 2, 0.0
        )
    tolerated = IdealGasMolecule(
        -14.22, modes, water().masses_amu, water().positions_angstrom, 2, 0.0,
        imaginary_tolerance_cm1=50.0,
    )
    assert tolerated.vibrational_frequencies_cm1 == (1595.0, 3756.0)


def test_molecule_validation_names_the_problem() -> None:
    with pytest.raises(ValueError, match="expected 6 modes"):
        oxygen((1580.0, 10.0))
    with pytest.raises(ValueError, match="symmetry_number"):
        IdealGasMolecule(-9.86, (O2_MODE,), (16.0, 16.0), ((0, 0, 0), (0, 0, 1.2)), 0, 1.0)
    with pytest.raises(ValueError, match="spin"):
        IdealGasMolecule(-9.86, (O2_MODE,), (16.0, 16.0), ((0, 0, 0), (0, 0, 1.2)), 2, 0.3)
    with pytest.raises(ValueError, match="contradicts"):
        IdealGasMolecule(
            -9.86, (O2_MODE,), (16.0, 16.0), ((0, 0, 0), (0, 0, 1.2)), 2, 1.0, geometry="nonlinear"
        )
    with pytest.raises(ValueError, match="pressure_bar"):
        oxygen().free_energy(300.0, 0.0)
    with pytest.raises(ValueError, match="temperature_k"):
        oxygen().free_energy(-1.0)


def test_from_structure_detects_geometry_across_the_cell_boundary() -> None:
    from pymatgen.core import Lattice, Structure

    # O2 split by the periodic boundary of a 10 A box.
    split = Structure(Lattice.cubic(10.0), ["O", "O"], [[0.5, 0.5, 0.97], [0.5, 0.5, 0.0907]])
    molecule = IdealGasMolecule.from_structure(
        split, electronic_energy_ev=-9.86, frequencies_cm1=(O2_MODE,), symmetry_number=2, spin=1.0
    )
    assert molecule.geometry == "linear"
    assert molecule.free_energy(298.15).entropy_ev_per_k == pytest.approx(
        oxygen().free_energy(298.15).entropy_ev_per_k, rel=1e-3
    )
    triplet = IdealGasMolecule.from_structure(
        psteros.triplet_o2_cell(), electronic_energy_ev=-9.86, frequencies_cm1=(O2_MODE,),
        symmetry_number=2, spin=1.0,
    )
    assert triplet.geometry == "linear"


def test_harmonic_solid_scales_supercell_modes_to_the_reference_cell() -> None:
    # Two-atom cell, 16-atom supercell: 3 acoustic + 45 optical modes.
    modes = (0.5, -0.3, 0.2) + tuple(100.0 + 10.0 * i for i in range(45))
    supercell = HarmonicSolid(-50.0, 16, modes, 16)
    cell = HarmonicSolid(-50.0 / 8, 2, modes, 16)
    assert cell.vibrational_frequencies_cm1 == supercell.vibrational_frequencies_cm1
    big, small = supercell.free_energy(600.0), cell.free_energy(600.0)
    assert small.vibrational_free_energy_ev == pytest.approx(big.vibrational_free_energy_ev / 8)
    assert small.pressure_bar is None


def test_harmonic_solid_high_temperature_limit() -> None:
    # Classical limit of one mode: F -> kT ln(h nu / kT).
    modes = (0.0, 0.0, 0.0, 50.0, 50.0, 50.0)
    solid = HarmonicSolid(0.0, 2, modes, 2)
    temperature = 5000.0
    kt = BOLTZMANN_EV_PER_K * temperature
    expected = 3 * kt * __import__("math").log(50.0 * CM1_TO_EV / kt)
    assert solid.free_energy(temperature).vibrational_free_energy_ev == pytest.approx(expected, rel=1e-3)


def test_harmonic_solid_validation() -> None:
    with pytest.raises(ValueError, match="atoms_in_supercell"):
        HarmonicSolid(-1.0, 2, (100.0,) * 3, 0)
    with pytest.raises(ValueError, match="expected 6 modes"):
        HarmonicSolid(-1.0, 2, (100.0,) * 4, 2)
    with pytest.raises(ValueError, match="imaginary"):
        HarmonicSolid(-1.0, 2, (0.0, 0.0, 0.0, -80.0, 100.0, 100.0), 2)


def test_matches_ase_thermochemistry() -> None:
    thermochemistry = pytest.importorskip("ase.thermochemistry")
    from ase import Atoms

    h2o = water()
    atoms = Atoms("OH2", positions=h2o.positions_angstrom, masses=h2o.masses_amu)
    ase_gas = thermochemistry.IdealGasThermo(
        vib_energies=[f * CM1_TO_EV for f in h2o.vibrational_frequencies_cm1],
        geometry="nonlinear",
        potentialenergy=h2o.electronic_energy_ev,
        atoms=atoms,
        symmetrynumber=2,
        spin=0,
    )
    for temperature, pressure in ((298.15, 1.0), (900.0, 1.0e-4)):
        expected = ase_gas.get_gibbs_energy(temperature, pressure * 1.0e5, verbose=False)
        assert h2o.free_energy(temperature, pressure).free_energy_ev == pytest.approx(expected, abs=2e-4)

    modes = tuple(80.0 + 25.0 * i for i in range(9))
    ase_solid = thermochemistry.HarmonicThermo(
        vib_energies=[f * CM1_TO_EV for f in modes], potentialenergy=-12.0
    )
    solid = HarmonicSolid(-12.0, 4, modes, 4)
    assert solid.free_energy(700.0).free_energy_ev == pytest.approx(
        ase_solid.get_helmholtz_energy(700.0, verbose=False), abs=1e-6
    )


OUTCAR_SNIPPET = """\
 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------


   1 f  =   47.316210 THz   297.297474 2PiTHz 1578.301557 cm-1   195.683925 meV
             X         Y         Z           dx          dy          dz
      9.000000  9.000000  8.396000            0.000000    0.000000   -0.707107
      9.000000  9.000000  9.604000            0.000000    0.000000    0.707107

   2 f  =    0.912345 THz     5.732458 2PiTHz   30.432812 cm-1     3.773184 meV
             X         Y         Z           dx          dy          dz
      9.000000  9.000000  8.396000            0.707107    0.000000    0.000000
      9.000000  9.000000  9.604000           -0.707107    0.000000    0.000000

   3 f/i=    0.212345 THz     1.334200 2PiTHz    7.083031 cm-1     0.878184 meV
             X         Y         Z           dx          dy          dz
      9.000000  9.000000  8.396000            0.000000    0.707107    0.000000
      9.000000  9.000000  9.604000            0.000000    0.707107    0.000000

 Eigenvectors after division by SQRT(mass)

 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------


   1 f  =   99.000000 THz   622.035345 2PiTHz 3302.293000 cm-1   409.430000 meV
"""


def test_parse_vasp_frequencies_reads_the_first_set_with_imaginary_sign() -> None:
    outcar = "header\n" + OUTCAR_SNIPPET.split(" Eigenvectors after division")[0] + (
        " Eigenvectors after division by SQRT(mass)\n\n"
        "   1 f  =   47.316210 THz   297.297474 2PiTHz 9999.000000 cm-1   195.683925 meV\n"
    )
    assert parse_vasp_frequencies_cm1(outcar) == pytest.approx((1578.301557, 30.432812, -7.083031))


def test_parse_vasp_frequencies_uses_the_last_dynamical_matrix() -> None:
    assert parse_vasp_frequencies_cm1(OUTCAR_SNIPPET) == pytest.approx((3302.293,))
    with pytest.raises(ValueError, match="IBRION"):
        parse_vasp_frequencies_cm1("no modes here")
