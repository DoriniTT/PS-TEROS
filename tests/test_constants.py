"""
Tests for Physical Constants and Unit Conversions

This module verifies that all physical constants in teros.core.constants
match CODATA 2018 recommended values and that unit conversions are
internally consistent.

References:
    - CODATA 2018: https://physics.nist.gov/cuu/Constants/
    - NIST Reference on Constants, Units, and Uncertainty
"""

import pytest
import math
from teros.core import constants


@pytest.mark.tier1
class TestEnergyUnitConversions:
    """Test energy unit conversion constants."""

    def test_ev_per_angstrom2_to_j_per_m2(self):
        """Test eV/Å² to J/m² conversion factor.

        Derivation: 1 eV = 1.602176634e-19 J, 1 Å = 1e-10 m
        1 eV/Å² = 1.602176634e-19 J / (1e-10 m)² = 16.02176634 J/m²
        """
        expected = 16.02176634
        assert math.isclose(
            constants.EV_PER_ANGSTROM2_TO_J_PER_M2,
            expected,
            rel_tol=1e-6,
        ), f"Expected {expected}, got {constants.EV_PER_ANGSTROM2_TO_J_PER_M2}"

    def test_ev_to_kj_per_mol(self):
        """Test eV to kJ/mol conversion factor.

        Derivation: 1 eV × N_A = 1.602176634e-19 × 6.02214076e23
                    = 96485.33212 J/mol = 96.48533212 kJ/mol
        """
        expected = 96.48533212
        assert math.isclose(
            constants.EV_TO_KJ_PER_MOL,
            expected,
            rel_tol=1e-3,  # Stored value is rounded to 96.485
        ), f"Expected ~{expected}, got {constants.EV_TO_KJ_PER_MOL}"

    def test_ev_to_joule(self):
        """Test eV to Joule conversion (CODATA 2018 exact value)."""
        expected = 1.602176634e-19  # Exact as of 2019 SI redefinition
        assert (
            constants.EV_TO_JOULE == expected
        ), f"Expected {expected}, got {constants.EV_TO_JOULE}"

    def test_ev_to_mev(self):
        """Test eV to meV conversion."""
        expected = 1000.0
        assert (
            constants.EV_TO_MEV == expected
        ), f"Expected {expected}, got {constants.EV_TO_MEV}"

    def test_hartree_to_ev(self):
        """Test Hartree to eV conversion (CODATA 2018).

        CODATA 2018: 27.211386245988(53) eV
        """
        expected = 27.211386245988
        assert math.isclose(
            constants.HARTREE_TO_EV,
            expected,
            rel_tol=1e-10,
        ), f"Expected {expected}, got {constants.HARTREE_TO_EV}"

    def test_rydberg_to_ev(self):
        """Test Rydberg to eV conversion (CODATA 2018).

        CODATA 2018: 13.605693122994(26) eV (Ry = 0.5 Hartree)
        """
        expected = 13.605693122994
        assert math.isclose(
            constants.RYDBERG_TO_EV,
            expected,
            rel_tol=1e-10,
        ), f"Expected {expected}, got {constants.RYDBERG_TO_EV}"

    def test_rydberg_hartree_relationship(self):
        """Test that Rydberg = 0.5 × Hartree."""
        expected_ratio = 0.5
        actual_ratio = constants.RYDBERG_TO_EV / constants.HARTREE_TO_EV
        assert math.isclose(
            actual_ratio, expected_ratio, rel_tol=1e-10
        ), f"Rydberg should be half of Hartree, got ratio {actual_ratio}"


@pytest.mark.tier1
class TestLengthUnitConversions:
    """Test length unit conversion constants."""

    def test_angstrom_to_meter(self):
        """Test Ångström to meter conversion."""
        expected = 1e-10
        assert (
            constants.ANGSTROM_TO_METER == expected
        ), f"Expected {expected}, got {constants.ANGSTROM_TO_METER}"

    def test_bohr_to_angstrom(self):
        """Test Bohr radius to Ångström conversion (CODATA 2018).

        CODATA 2018: 0.529177210903(80) Å
        """
        expected = 0.529177210903
        assert math.isclose(
            constants.BOHR_TO_ANGSTROM,
            expected,
            rel_tol=1e-10,
        ), f"Expected {expected}, got {constants.BOHR_TO_ANGSTROM}"

    def test_angstrom_to_bohr(self):
        """Test Ångström to Bohr conversion (derived)."""
        expected = 1.0 / 0.529177210903
        assert math.isclose(
            constants.ANGSTROM_TO_BOHR,
            expected,
            rel_tol=1e-10,
        ), f"Expected {expected}, got {constants.ANGSTROM_TO_BOHR}"

    def test_bohr_angstrom_reciprocal(self):
        """Test that BOHR_TO_ANGSTROM × ANGSTROM_TO_BOHR = 1."""
        product = constants.BOHR_TO_ANGSTROM * constants.ANGSTROM_TO_BOHR
        assert math.isclose(
            product, 1.0, rel_tol=1e-10
        ), f"Product should be 1.0, got {product}"


@pytest.mark.tier1
class TestFundamentalConstants:
    """Test fundamental physical constants."""

    def test_avogadro(self):
        """Test Avogadro constant (CODATA 2018 exact value).

        CODATA 2018: 6.02214076e23 mol⁻¹ (exact as of 2019 SI redefinition)
        """
        expected = 6.02214076e23
        assert (
            constants.AVOGADRO == expected
        ), f"Expected {expected}, got {constants.AVOGADRO}"

    def test_boltzmann_ev(self):
        """Test Boltzmann constant in eV/K (CODATA 2018).

        CODATA 2018: 8.617333262...e-5 eV/K
        """
        expected = 8.617333262e-5
        assert math.isclose(
            constants.BOLTZMANN_EV,
            expected,
            rel_tol=1e-8,
        ), f"Expected {expected}, got {constants.BOLTZMANN_EV}"

    def test_boltzmann_j(self):
        """Test Boltzmann constant in J/K (CODATA 2018 exact value).

        CODATA 2018: 1.380649e-23 J/K (exact as of 2019 SI redefinition)
        """
        expected = 1.380649e-23
        assert (
            constants.BOLTZMANN_J == expected
        ), f"Expected {expected}, got {constants.BOLTZMANN_J}"

    def test_gas_constant(self):
        """Test gas constant (CODATA 2018).

        CODATA 2018: 8.314462618... J/(mol·K)
        Derived: R = N_A × k_B
        """
        expected = 8.314462618
        assert math.isclose(
            constants.GAS_CONSTANT,
            expected,
            rel_tol=1e-8,
        ), f"Expected {expected}, got {constants.GAS_CONSTANT}"

    def test_gas_constant_derived(self):
        """Test that R = N_A × k_B."""
        derived_r = constants.AVOGADRO * constants.BOLTZMANN_J
        assert math.isclose(
            constants.GAS_CONSTANT, derived_r, rel_tol=1e-6
        ), f"R should equal N_A × k_B, got {derived_r} vs {constants.GAS_CONSTANT}"

    def test_planck_ev(self):
        """Test Planck constant in eV·s (CODATA 2018).

        CODATA 2018: h = 6.62607015e-34 J·s (exact)
        h (eV·s) = 6.62607015e-34 / 1.602176634e-19 = 4.135667696e-15 eV·s
        """
        expected = 4.135667696e-15
        assert math.isclose(
            constants.PLANCK_EV,
            expected,
            rel_tol=1e-8,
        ), f"Expected {expected}, got {constants.PLANCK_EV}"


@pytest.mark.tier1
class TestDerivedRelationships:
    """Test consistency between related constants."""

    def test_ev_to_kj_per_mol_derivation(self):
        """Test that EV_TO_KJ_PER_MOL = EV_TO_JOULE × AVOGADRO / 1000."""
        derived = constants.EV_TO_JOULE * constants.AVOGADRO / 1000.0
        assert math.isclose(
            constants.EV_TO_KJ_PER_MOL, derived, rel_tol=1e-3
        ), f"Derived {derived}, stored {constants.EV_TO_KJ_PER_MOL}"

    def test_ev_per_angstrom2_to_j_per_m2_derivation(self):
        """Test eV/Å² to J/m² conversion derivation."""
        # 1 eV/Å² = (eV in J) / (Å² in m²)
        # = 1.602176634e-19 / (1e-10)² = 1.602176634e-19 / 1e-20
        derived = constants.EV_TO_JOULE / (constants.ANGSTROM_TO_METER**2)
        assert math.isclose(
            constants.EV_PER_ANGSTROM2_TO_J_PER_M2, derived, rel_tol=1e-6
        ), f"Derived {derived}, stored {constants.EV_PER_ANGSTROM2_TO_J_PER_M2}"

    def test_boltzmann_ev_j_conversion(self):
        """Test consistency between Boltzmann constants in different units."""
        # k_B (eV/K) = k_B (J/K) / (eV to J conversion)
        derived_boltzmann_ev = constants.BOLTZMANN_J / constants.EV_TO_JOULE
        assert math.isclose(
            constants.BOLTZMANN_EV, derived_boltzmann_ev, rel_tol=1e-6
        ), f"Derived {derived_boltzmann_ev}, stored {constants.BOLTZMANN_EV}"


@pytest.mark.tier1
class TestNumericalTolerances:
    """Test numerical tolerance constants."""

    def test_stoichiometry_rtol(self):
        """Test stoichiometry relative tolerance."""
        assert constants.STOICHIOMETRY_RTOL == 1e-3
        assert isinstance(constants.STOICHIOMETRY_RTOL, float)

    def test_float_atol(self):
        """Test floating point absolute tolerance."""
        assert constants.FLOAT_ATOL == 1e-10
        assert isinstance(constants.FLOAT_ATOL, float)

    def test_tolerances_positive(self):
        """Test that all tolerances are positive."""
        assert constants.STOICHIOMETRY_RTOL > 0
        assert constants.FLOAT_ATOL > 0


@pytest.mark.tier1
class TestVASPConstants:
    """Test VASP-specific constants."""

    def test_default_kpoints_spacing(self):
        """Test default k-point spacing."""
        assert constants.DEFAULT_KPOINTS_SPACING == 0.03
        assert isinstance(constants.DEFAULT_KPOINTS_SPACING, float)

    def test_encut_values(self):
        """Test ENCUT preset values."""
        assert constants.ENCUT_LOW == 400
        assert constants.ENCUT_STANDARD == 520
        assert constants.ENCUT_HIGH == 600

    def test_encut_ordering(self):
        """Test that ENCUT values are properly ordered."""
        assert constants.ENCUT_LOW < constants.ENCUT_STANDARD < constants.ENCUT_HIGH

    def test_encut_types(self):
        """Test that ENCUT values are integers."""
        assert isinstance(constants.ENCUT_LOW, int)
        assert isinstance(constants.ENCUT_STANDARD, int)
        assert isinstance(constants.ENCUT_HIGH, int)


@pytest.mark.tier1
class TestTemperatureConstants:
    """Test temperature constants."""

    def test_room_temperature(self):
        """Test room temperature value."""
        assert constants.ROOM_TEMPERATURE == 298.15
        assert isinstance(constants.ROOM_TEMPERATURE, float)

    def test_standard_temperature(self):
        """Test standard state temperature."""
        assert constants.STANDARD_TEMPERATURE == 298.15
        assert isinstance(constants.STANDARD_TEMPERATURE, float)

    def test_absolute_zero(self):
        """Test absolute zero."""
        assert constants.ABSOLUTE_ZERO == 0.0
        assert isinstance(constants.ABSOLUTE_ZERO, float)

    def test_room_standard_temperature_match(self):
        """Test that room temperature equals standard temperature."""
        assert constants.ROOM_TEMPERATURE == constants.STANDARD_TEMPERATURE

    def test_temperatures_physical(self):
        """Test that temperatures are physically reasonable."""
        assert constants.ABSOLUTE_ZERO >= 0.0
        assert constants.ROOM_TEMPERATURE > constants.ABSOLUTE_ZERO
        assert constants.STANDARD_TEMPERATURE > constants.ABSOLUTE_ZERO


@pytest.mark.tier1
class TestModuleExports:
    """Test that all constants are accessible from the module."""

    def test_energy_conversions_exported(self):
        """Test energy conversion constants are accessible."""
        assert hasattr(constants, 'EV_PER_ANGSTROM2_TO_J_PER_M2')
        assert hasattr(constants, 'EV_TO_KJ_PER_MOL')
        assert hasattr(constants, 'EV_TO_JOULE')
        assert hasattr(constants, 'EV_TO_MEV')
        assert hasattr(constants, 'HARTREE_TO_EV')
        assert hasattr(constants, 'RYDBERG_TO_EV')

    def test_length_conversions_exported(self):
        """Test length conversion constants are accessible."""
        assert hasattr(constants, 'ANGSTROM_TO_METER')
        assert hasattr(constants, 'BOHR_TO_ANGSTROM')
        assert hasattr(constants, 'ANGSTROM_TO_BOHR')

    def test_fundamental_constants_exported(self):
        """Test fundamental constants are accessible."""
        assert hasattr(constants, 'AVOGADRO')
        assert hasattr(constants, 'BOLTZMANN_EV')
        assert hasattr(constants, 'BOLTZMANN_J')
        assert hasattr(constants, 'GAS_CONSTANT')
        assert hasattr(constants, 'PLANCK_EV')

    def test_tolerances_exported(self):
        """Test tolerance constants are accessible."""
        assert hasattr(constants, 'STOICHIOMETRY_RTOL')
        assert hasattr(constants, 'FLOAT_ATOL')

    def test_vasp_constants_exported(self):
        """Test VASP constants are accessible."""
        assert hasattr(constants, 'DEFAULT_KPOINTS_SPACING')
        assert hasattr(constants, 'ENCUT_LOW')
        assert hasattr(constants, 'ENCUT_STANDARD')
        assert hasattr(constants, 'ENCUT_HIGH')

    def test_temperature_constants_exported(self):
        """Test temperature constants are accessible."""
        assert hasattr(constants, 'ROOM_TEMPERATURE')
        assert hasattr(constants, 'STANDARD_TEMPERATURE')
        assert hasattr(constants, 'ABSOLUTE_ZERO')


@pytest.mark.tier1
class TestConstantTypes:
    """Test that all constants have correct types."""

    def test_energy_conversion_types(self):
        """Test energy conversion constants are floats."""
        assert isinstance(constants.EV_PER_ANGSTROM2_TO_J_PER_M2, float)
        assert isinstance(constants.EV_TO_KJ_PER_MOL, float)
        assert isinstance(constants.EV_TO_JOULE, float)
        assert isinstance(constants.EV_TO_MEV, float)
        assert isinstance(constants.HARTREE_TO_EV, float)
        assert isinstance(constants.RYDBERG_TO_EV, float)

    def test_length_conversion_types(self):
        """Test length conversion constants are floats."""
        assert isinstance(constants.ANGSTROM_TO_METER, float)
        assert isinstance(constants.BOHR_TO_ANGSTROM, float)
        assert isinstance(constants.ANGSTROM_TO_BOHR, float)

    def test_fundamental_constant_types(self):
        """Test fundamental constants are floats."""
        assert isinstance(constants.AVOGADRO, float)
        assert isinstance(constants.BOLTZMANN_EV, float)
        assert isinstance(constants.BOLTZMANN_J, float)
        assert isinstance(constants.GAS_CONSTANT, float)
        assert isinstance(constants.PLANCK_EV, float)


@pytest.mark.tier1
class TestPhysicalReasonableness:
    """Test that constants have physically reasonable values."""

    def test_positive_conversion_factors(self):
        """Test that all conversion factors are positive."""
        assert constants.EV_PER_ANGSTROM2_TO_J_PER_M2 > 0
        assert constants.EV_TO_KJ_PER_MOL > 0
        assert constants.EV_TO_JOULE > 0
        assert constants.EV_TO_MEV > 0
        assert constants.HARTREE_TO_EV > 0
        assert constants.RYDBERG_TO_EV > 0
        assert constants.BOHR_TO_ANGSTROM > 0
        assert constants.ANGSTROM_TO_BOHR > 0

    def test_positive_fundamental_constants(self):
        """Test that fundamental constants are positive."""
        assert constants.AVOGADRO > 0
        assert constants.BOLTZMANN_EV > 0
        assert constants.BOLTZMANN_J > 0
        assert constants.GAS_CONSTANT > 0
        assert constants.PLANCK_EV > 0

    def test_reasonable_magnitudes(self):
        """Test that constants have reasonable orders of magnitude."""
        # Surface energy conversion should be ~10-20 J/m² per eV/Å²
        assert 10 < constants.EV_PER_ANGSTROM2_TO_J_PER_M2 < 20

        # Energy conversion should be ~90-100 kJ/mol per eV
        assert 90 < constants.EV_TO_KJ_PER_MOL < 100

        # Avogadro should be ~6e23
        assert 6e23 < constants.AVOGADRO < 7e23

        # Room temperature should be ~298 K
        assert 290 < constants.ROOM_TEMPERATURE < 310
