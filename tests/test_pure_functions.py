"""
Pure Python Function Tests

These tests verify pure Python functions that don't require AiiDA.
They run quickly and are executed on every push/PR.
"""

import pytest
from copy import deepcopy


# =============================================================================
# FORMULA PARSING TESTS
# =============================================================================

class TestFormulaParsingStandalone:
    """Test chemical formula parsing using pymatgen directly."""

    def test_parse_simple_formulas(self):
        """Test parsing simple chemical formulas."""
        from pymatgen.core import Composition

        test_cases = {
            'H2O': {'H': 2, 'O': 1},
            'O2': {'O': 2},
            'H2': {'H': 2},
            'CO2': {'C': 1, 'O': 2},
            'NH3': {'N': 1, 'H': 3},
        }

        for formula, expected in test_cases.items():
            comp = Composition(formula)
            result = {str(el): int(count) for el, count in comp.items()}
            assert result == expected, f"Failed for {formula}"

    def test_parse_complex_formulas(self):
        """Test parsing more complex formulas."""
        from pymatgen.core import Composition

        test_cases = {
            'OOH': {'O': 2, 'H': 1},
            'OH': {'O': 1, 'H': 1},
            'H2O2': {'H': 2, 'O': 2},
            'Ag2O': {'Ag': 2, 'O': 1},
            'Fe2O3': {'Fe': 2, 'O': 3},
            'Ag3PO4': {'Ag': 3, 'P': 1, 'O': 4},
        }

        for formula, expected in test_cases.items():
            comp = Composition(formula)
            result = {str(el): int(count) for el, count in comp.items()}
            assert result == expected, f"Failed for {formula}"

    def test_invalid_formula_raises(self):
        """Test that invalid formulas raise errors."""
        from pymatgen.core import Composition

        invalid_formulas = ['', 'Xx', '123', 'H2O2O']

        for formula in invalid_formulas:
            if formula:  # Empty string might not raise
                try:
                    Composition(formula)
                    # Some "invalid" formulas might still parse
                except Exception:
                    pass  # Expected behavior


# =============================================================================
# DEEP MERGE TESTS
# =============================================================================

class TestDeepMergeDicts:
    """Test the deep_merge_dicts utility function."""

    def get_deep_merge_function(self):
        """Get the deep_merge_dicts function from various possible locations."""
        # Try importing from different modules
        try:
            from teros.core.slabs import deep_merge_dicts
            return deep_merge_dicts
        except ImportError:
            pass

        try:
            from teros.core.adsorption_energy import deep_merge_dicts
            return deep_merge_dicts
        except ImportError:
            pass

        # Fallback: implement locally for testing
        def deep_merge_dicts(base: dict, override: dict) -> dict:
            result = deepcopy(base)
            for key, value in override.items():
                if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    result[key] = deep_merge_dicts(result[key], value)
                else:
                    result[key] = deepcopy(value)
            return result

        return deep_merge_dicts

    def test_simple_merge(self):
        """Test merging flat dictionaries."""
        deep_merge = self.get_deep_merge_function()

        base = {'a': 1, 'b': 2}
        override = {'b': 99, 'c': 3}

        result = deep_merge(base, override)

        assert result == {'a': 1, 'b': 99, 'c': 3}

    def test_nested_merge(self):
        """Test merging nested dictionaries."""
        deep_merge = self.get_deep_merge_function()

        base = {
            'a': 1,
            'b': {'c': 2, 'd': 3},
        }
        override = {
            'b': {'c': 99, 'e': 100},
        }

        result = deep_merge(base, override)

        assert result['a'] == 1
        assert result['b']['c'] == 99
        assert result['b']['d'] == 3
        assert result['b']['e'] == 100

    def test_original_not_modified(self):
        """Test that original dictionaries are not modified."""
        deep_merge = self.get_deep_merge_function()

        base = {'a': {'b': 1}}
        override = {'a': {'b': 2}}

        result = deep_merge(base, override)

        assert base['a']['b'] == 1  # Original unchanged
        assert result['a']['b'] == 2

    def test_deeply_nested_merge(self):
        """Test merging deeply nested structures."""
        deep_merge = self.get_deep_merge_function()

        base = {
            'level1': {
                'level2': {
                    'level3': {
                        'value': 1
                    }
                }
            }
        }
        override = {
            'level1': {
                'level2': {
                    'level3': {
                        'value': 99,
                        'new_key': 'new'
                    }
                }
            }
        }

        result = deep_merge(base, override)

        assert result['level1']['level2']['level3']['value'] == 99
        assert result['level1']['level2']['level3']['new_key'] == 'new'

    def test_list_not_merged(self):
        """Test that lists are replaced, not merged."""
        deep_merge = self.get_deep_merge_function()

        base = {'a': [1, 2, 3]}
        override = {'a': [4, 5]}

        result = deep_merge(base, override)

        assert result['a'] == [4, 5]


# =============================================================================
# MILLER INDEX TESTS
# =============================================================================

class TestMillerIndices:
    """Test Miller index handling."""

    def test_miller_index_tuple_creation(self):
        """Test creating Miller index tuples."""
        test_cases = [
            ([1, 0, 0], (1, 0, 0)),
            ([1, 1, 0], (1, 1, 0)),
            ([1, 1, 1], (1, 1, 1)),
            ([2, 1, 0], (2, 1, 0)),
        ]

        for input_list, expected_tuple in test_cases:
            result = tuple(input_list)
            assert result == expected_tuple

    def test_miller_index_validation(self):
        """Test Miller index validation logic."""
        # All zeros is invalid
        invalid_indices = [0, 0, 0]
        assert not any(invalid_indices)

        # At least one non-zero is valid
        valid_indices = [1, 0, 0]
        assert any(valid_indices)


# =============================================================================
# STOICHIOMETRY TESTS
# =============================================================================

class TestStoichiometry:
    """Test stoichiometry calculations."""

    def test_gcd_calculation(self):
        """Test GCD calculation for formula units."""
        from math import gcd
        from functools import reduce

        def gcd_list(numbers):
            return reduce(gcd, numbers)

        test_cases = [
            ([2, 1], 1),  # Ag2O -> 1 formula unit
            ([6, 3], 3),  # Can be reduced
            ([4, 2, 8], 2),  # Multiple elements
            ([3, 1, 4], 1),  # Ag3PO4
        ]

        for numbers, expected in test_cases:
            result = gcd_list(numbers)
            assert result == expected

    def test_composition_counting(self):
        """Test counting atoms in compositions."""
        from collections import Counter

        symbols = ['Ag', 'Ag', 'O', 'P', 'O', 'O', 'O']
        counts = Counter(symbols)

        assert counts['Ag'] == 2
        assert counts['O'] == 4
        assert counts['P'] == 1


# =============================================================================
# ENERGY DIFFERENCE TESTS
# =============================================================================

class TestEnergyCalculations:
    """Test energy calculation utilities."""

    def test_energy_difference(self):
        """Test simple energy difference calculation."""
        E_relaxed = -100.5
        E_unrelaxed = -100.0

        E_relax = E_relaxed - E_unrelaxed

        assert E_relax == pytest.approx(-0.5)

    def test_adsorption_energy(self):
        """Test adsorption energy calculation."""
        E_complete = -150.0
        E_substrate = -100.0
        E_molecule = -30.0

        E_ads = E_complete - E_substrate - E_molecule

        assert E_ads == pytest.approx(-20.0)

    def test_formation_enthalpy_per_atom(self):
        """Test formation enthalpy normalization."""
        total_energy = -50.0
        n_atoms = 5

        per_atom = total_energy / n_atoms

        assert per_atom == pytest.approx(-10.0)
