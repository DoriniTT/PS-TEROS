"""
DOS Calculation Module Tests

These tests verify that the DOS calculation WorkGraph functions work correctly.
"""

import pytest


def check_aiida():
    """Check if AiiDA is available and configured."""
    try:
        from aiida import load_profile
        load_profile()
        return True
    except Exception:
        return False


AIIDA_AVAILABLE = check_aiida()

if not AIIDA_AVAILABLE:
    pytest.skip("AiiDA not configured", allow_module_level=True)


# =============================================================================
# DOS FUNCTION TESTS
# =============================================================================

class TestDOSFunctions:
    """Test DOS calculation functions without AiiDA requirement."""

    def test_build_dos_calculation_workgraph_import(self):
        """Test that build_dos_calculation_workgraph can be imported."""
        from teros.core.custom_calculation import build_dos_calculation_workgraph
        assert callable(build_dos_calculation_workgraph)

    def test_get_dos_results_import(self):
        """Test that get_dos_results can be imported."""
        from teros.core.custom_calculation import get_dos_results
        assert callable(get_dos_results)

    def test_dos_functions_in_init(self):
        """Test that DOS functions are properly exported in __init__.py."""
        from teros.core.custom_calculation import (
            build_dos_calculation_workgraph,
            get_dos_results,
        )
        assert build_dos_calculation_workgraph is not None
        assert get_dos_results is not None


class TestDOSWorkGraphConstruction:
    """Test DOS WorkGraph construction.

    Note: These tests require a configured AiiDA code for VASP.
    They are skipped if no code is available.
    """

    @pytest.fixture
    def mock_structure(self):
        """Create a mock structure for testing."""
        from aiida import orm
        from ase import Atoms

        # Create a simple test structure
        atoms = Atoms(
            'Ag4O2',
            positions=[
                [0, 0, 0], [2, 2, 0], [2, 0, 2], [0, 2, 2],
                [1, 1, 1], [3, 3, 1]
            ],
            cell=[4, 4, 4],
            pbc=True
        )
        structure = orm.StructureData(ase=atoms)
        structure.label = 'test_structure'
        return structure

    @pytest.fixture
    def mock_code_label(self):
        """Return a mock code label for testing."""
        from aiida import orm
        try:
            qb = orm.QueryBuilder()
            qb.append(orm.Code, project=['label', 'uuid'])
            codes = qb.all()
            if codes:
                label_counts = {}
                for label, uuid in codes:
                    label_counts[label] = label_counts.get(label, 0) + 1

                for label, count in label_counts.items():
                    if count == 1:
                        return label

                if codes:
                    label, uuid = codes[0]
                    return f"{label}@{uuid[:8]}"
        except Exception:
            pass
        pytest.skip("No AiiDA code configured for testing")

    @pytest.fixture
    def scf_inputs(self):
        """Return SCF inputs for testing."""
        return {
            'parameters': {
                'incar': {
                    'PREC': 'Accurate',
                    'ENCUT': 500,
                    'EDIFF': 1e-5,
                    'ALGO': 'Fast',
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'LWAVE': True,
                    'LCHARG': True,
                    'LORBIT': 11,
                }
            },
            'kpoints_spacing': 0.3,
        }

    @pytest.fixture
    def dos_inputs(self):
        """Return DOS inputs for testing."""
        return {
            'parameters': {
                'incar': {
                    'PREC': 'Accurate',
                    'ENCUT': 500,
                    'EDIFF': 1e-5,
                    'ALGO': 'Normal',
                    'ISTART': 1,
                    'ICHARG': 11,
                    'ISMEAR': -5,
                    'LORBIT': 11,
                    'NEDOS': 2001,
                }
            },
        }

    @pytest.mark.skip(reason="Requires unique configured VASP code")
    def test_build_dos_workgraph_single_structure(
        self, mock_structure, mock_code_label, scf_inputs, dos_inputs
    ):
        """Test building a DOS workgraph for a single structure."""
        from teros.core.custom_calculation import build_dos_calculation_workgraph

        wg = build_dos_calculation_workgraph(
            structure=mock_structure,
            code_label=mock_code_label,
            scf_inputs=scf_inputs,
            dos_inputs=dos_inputs,
            dos_kpoints_distance=0.2,
            potential_family='PBE',
            potential_mapping={'Ag': 'Ag', 'O': 'O'},
            name='test_dos_single',
        )

        assert wg is not None
        assert wg.name == 'test_dos_single'

        # Check that expected tasks exist
        task_names = list(wg.tasks.keys())
        assert 'dos_calc' in task_names

    @pytest.mark.skip(reason="Requires unique configured VASP code")
    def test_build_dos_workgraph_multiple_structures(
        self, mock_structure, mock_code_label, scf_inputs, dos_inputs
    ):
        """Test building a DOS workgraph for multiple structures."""
        from teros.core.custom_calculation import build_dos_calculation_workgraph
        from aiida import orm
        from ase import Atoms

        # Create second structure
        atoms2 = Atoms(
            'Ag4O2',
            positions=[
                [0, 0, 0], [2, 2, 0], [2, 0, 2], [0, 2, 2],
                [1, 1, 1], [3, 3, 1]
            ],
            cell=[4.1, 4.1, 4.1],
            pbc=True
        )
        structure2 = orm.StructureData(ase=atoms2)
        structure2.label = 'test_structure_2'

        structures = [mock_structure, structure2]

        wg = build_dos_calculation_workgraph(
            structure=structures,
            code_label=mock_code_label,
            scf_inputs=scf_inputs,
            dos_inputs=dos_inputs,
            dos_kpoints_distance=0.2,
            potential_family='PBE',
            potential_mapping={'Ag': 'Ag', 'O': 'O'},
            name='test_dos_multiple',
            max_concurrent_jobs=2,
        )

        assert wg is not None
        assert wg.name == 'test_dos_multiple'

        # Check that multiple DOS tasks exist
        task_names = list(wg.tasks.keys())
        assert len([t for t in task_names if t.startswith('dos_calc')]) == 2


class TestDOSResultsExtraction:
    """Test get_dos_results function."""

    def test_get_dos_results_returns_dict_with_expected_keys(self):
        """Test that get_dos_results returns a dict with expected keys."""
        from teros.core.custom_calculation import get_dos_results

        # Create a mock workgraph with no tasks
        class MockWorkGraph:
            tasks = {}

        mock_wg = MockWorkGraph()
        results = get_dos_results(mock_wg)

        assert isinstance(results, dict)
        assert 'dos' in results
        assert 'projectors' in results
        assert 'primitive_structure' in results
        assert isinstance(results['dos'], dict)
        assert isinstance(results['projectors'], dict)
        assert isinstance(results['primitive_structure'], dict)


class TestDOSInputPreparation:
    """Test internal DOS input preparation functions."""

    def test_sanitize_key_function(self):
        """Test the _sanitize_key function."""
        from teros.core.custom_calculation.workgraph import _sanitize_key

        # Test normal label
        assert _sanitize_key('my_structure', 0) == 'my_structure'

        # Test label starting with digit
        assert _sanitize_key('123_structure', 0) == 's_123_structure'

        # Test empty label
        assert _sanitize_key('', 0) == 'structure_0'

        # Test label with special characters
        assert _sanitize_key('my-structure.cif', 0) == 'my_structure_cif'

        # Test label with spaces
        assert _sanitize_key('my structure', 0) == 'my_structure'
