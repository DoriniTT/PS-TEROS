"""
Fukui Module Tests

These tests verify the Fukui module's pure Python functions that don't require AiiDA.
They run quickly and are executed on every push/PR.
"""

import pytest
import warnings


# =============================================================================
# LABEL GENERATION TESTS
# =============================================================================

class TestMakeDeltaLabel:
    """Test the make_delta_label utility function."""

    def test_zero_value(self):
        """Test label generation for delta_n = 0.0."""
        from teros.core.fukui.utils import make_delta_label

        result = make_delta_label(0.0)
        assert result == 'delta_0_00'

    def test_standard_values(self):
        """Test label generation for standard delta_n values."""
        from teros.core.fukui.utils import make_delta_label

        test_cases = [
            (0.05, 'delta_0_05'),
            (0.10, 'delta_0_10'),
            (0.15, 'delta_0_15'),
            (0.20, 'delta_0_20'),
        ]

        for delta_n, expected in test_cases:
            result = make_delta_label(delta_n)
            assert result == expected, f"Failed for delta_n={delta_n}"

    def test_unusual_values(self):
        """Test label generation for unusual delta_n values."""
        from teros.core.fukui.utils import make_delta_label

        test_cases = [
            (0.01, 'delta_0_01'),
            (0.001, 'delta_0_00'),  # Rounds to 2 decimal places
            (0.999, 'delta_1_00'),  # Rounds up
            (1.0, 'delta_1_00'),
        ]

        for delta_n, expected in test_cases:
            result = make_delta_label(delta_n)
            assert result == expected, f"Failed for delta_n={delta_n}"

    def test_returns_string(self):
        """Test that make_delta_label returns a string."""
        from teros.core.fukui.utils import make_delta_label

        result = make_delta_label(0.05)
        assert isinstance(result, str)

    def test_valid_python_identifier(self):
        """Test that generated labels are valid Python identifiers."""
        from teros.core.fukui.utils import make_delta_label

        test_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.5, 1.0]

        for delta_n in test_values:
            result = make_delta_label(delta_n)
            assert result.isidentifier(), f"'{result}' is not a valid Python identifier"


# =============================================================================
# INPUT VALIDATION TESTS
# =============================================================================

class TestValidateFukuiInputs:
    """Test the validate_fukui_inputs function."""

    def test_valid_inputs_plus(self):
        """Test validation passes for valid Fukui+ inputs."""
        from teros.core.fukui.utils import validate_fukui_inputs

        # Should not raise
        validate_fukui_inputs(
            nelect_neutral=100,
            delta_n_values=[0.0, 0.05, 0.10, 0.15],
            fukui_type='plus',
        )

    def test_valid_inputs_minus(self):
        """Test validation passes for valid Fukui- inputs."""
        from teros.core.fukui.utils import validate_fukui_inputs

        # Should not raise
        validate_fukui_inputs(
            nelect_neutral=100,
            delta_n_values=[0.0, 0.05, 0.10, 0.15],
            fukui_type='minus',
        )

    def test_invalid_nelect_neutral_zero(self):
        """Test validation fails for nelect_neutral = 0."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with pytest.raises(ValueError, match="nelect_neutral must be positive"):
            validate_fukui_inputs(
                nelect_neutral=0,
                delta_n_values=[0.0, 0.05],
                fukui_type='plus',
            )

    def test_invalid_nelect_neutral_negative(self):
        """Test validation fails for negative nelect_neutral."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with pytest.raises(ValueError, match="nelect_neutral must be positive"):
            validate_fukui_inputs(
                nelect_neutral=-10,
                delta_n_values=[0.0, 0.05],
                fukui_type='plus',
            )

    def test_invalid_negative_delta_n(self):
        """Test validation fails for negative delta_n values."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with pytest.raises(ValueError, match="non-negative"):
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.0, -0.05, 0.10],
                fukui_type='plus',
            )

    def test_invalid_fukui_type(self):
        """Test validation fails for invalid fukui_type."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with pytest.raises(ValueError, match="must be 'plus' or 'minus'"):
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.0, 0.05],
                fukui_type='invalid',
            )

    def test_invalid_delta_n_too_large(self):
        """Test validation fails when delta_n would result in NELECT < 1."""
        from teros.core.fukui.utils import validate_fukui_inputs

        # With nelect_neutral=5 and delta_n=10, Fukui+ would give NELECT=-5
        with pytest.raises(ValueError, match="delta_n values too large"):
            validate_fukui_inputs(
                nelect_neutral=5,
                delta_n_values=[0.0, 10.0],
                fukui_type='plus',
            )

    def test_too_many_delta_n_values(self):
        """Test validation fails for more than 4 delta_n values."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with pytest.raises(ValueError, match="Maximum of 4 delta_n values"):
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.0, 0.05, 0.10, 0.15, 0.20],
                fukui_type='plus',
            )

    def test_warning_missing_neutral_reference(self):
        """Test warning when delta_n=0.0 is not included."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.05, 0.10, 0.15],  # Missing 0.0
                fukui_type='plus',
            )
            assert len(w) >= 1
            assert any("delta_n=0.0" in str(warning.message) for warning in w)

    def test_warning_small_delta_n(self):
        """Test warning for very small delta_n values."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.0, 0.005, 0.01],  # 0.005 is very small
                fukui_type='plus',
            )
            assert len(w) >= 1
            assert any("Very small delta_n" in str(warning.message) for warning in w)

    def test_warning_large_delta_n(self):
        """Test warning for large delta_n values."""
        from teros.core.fukui.utils import validate_fukui_inputs

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_fukui_inputs(
                nelect_neutral=100,
                delta_n_values=[0.0, 0.05, 0.6],  # 0.6 is large
                fukui_type='plus',
            )
            assert len(w) >= 1
            assert any("Large delta_n" in str(warning.message) for warning in w)


# =============================================================================
# DEFAULT VALUES TESTS
# =============================================================================

class TestDefaultValues:
    """Test default constant values."""

    def test_default_delta_n_values_exists(self):
        """Test that DEFAULT_DELTA_N_VALUES is defined."""
        from teros.core.fukui.utils import DEFAULT_DELTA_N_VALUES

        assert DEFAULT_DELTA_N_VALUES is not None
        assert isinstance(DEFAULT_DELTA_N_VALUES, list)

    def test_default_delta_n_values_content(self):
        """Test that DEFAULT_DELTA_N_VALUES has expected content."""
        from teros.core.fukui.utils import DEFAULT_DELTA_N_VALUES

        # Check it includes neutral reference
        assert 0.0 in DEFAULT_DELTA_N_VALUES

        # Check it has reasonable number of points
        assert len(DEFAULT_DELTA_N_VALUES) >= 3
        assert len(DEFAULT_DELTA_N_VALUES) <= 4

        # Check all values are non-negative
        assert all(dn >= 0.0 for dn in DEFAULT_DELTA_N_VALUES)

    def test_default_delta_n_values_valid(self):
        """Test that DEFAULT_DELTA_N_VALUES passes validation."""
        from teros.core.fukui.utils import (
            DEFAULT_DELTA_N_VALUES,
            validate_fukui_inputs,
        )

        # Should not raise
        validate_fukui_inputs(
            nelect_neutral=100,
            delta_n_values=DEFAULT_DELTA_N_VALUES,
            fukui_type='plus',
        )


# =============================================================================
# MODULE IMPORT TESTS
# =============================================================================

class TestFukuiModuleImports:
    """Test that the fukui module can be imported correctly."""

    def test_import_main_module(self):
        """Test importing the main fukui module."""
        from teros.core import fukui
        assert fukui is not None

    def test_import_public_functions(self):
        """Test importing public functions from fukui module."""
        from teros.core.fukui import (
            build_fukui_workgraph,
            get_fukui_results,
            print_fukui_summary,
            make_delta_label,
            validate_fukui_inputs,
            DEFAULT_DELTA_N_VALUES,
        )

        # Verify they are callable/accessible
        assert callable(build_fukui_workgraph)
        assert callable(get_fukui_results)
        assert callable(print_fukui_summary)
        assert callable(make_delta_label)
        assert callable(validate_fukui_inputs)
        assert isinstance(DEFAULT_DELTA_N_VALUES, list)

    def test_import_calcfunctions(self):
        """Test importing calcfunctions from fukui module."""
        from teros.core.fukui import (
            collect_chgcar_files,
            generate_fukui_summary,
            extract_total_energy,
        )

        # These should be task-decorated functions
        assert callable(collect_chgcar_files)
        assert callable(generate_fukui_summary)
        assert callable(extract_total_energy)

    def test_import_phase2_functions(self):
        """Test importing Phase 2 (electrodes) functions."""
        from teros.core.fukui import (
            extract_dielectric_constant,
            run_fukui_electrodes_calcfunc,
        )

        assert callable(extract_dielectric_constant)
        assert callable(run_fukui_electrodes_calcfunc)

    def test_import_phase4_functions(self):
        """Test importing Phase 4 (perturbative expansion) functions."""
        from teros.core.fukui import (
            extract_locpot_from_retrieved,
            run_perturbative_expansion_calcfunc,
        )

        assert callable(extract_locpot_from_retrieved)
        assert callable(run_perturbative_expansion_calcfunc)
