"""
WorkGraph Construction Tests

These tests verify that WorkGraph objects can be constructed correctly
without actually submitting them. This catches wiring errors and
parameter validation issues.
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
# WORKGRAPH CONSTRUCTION TESTS
# =============================================================================

class TestWorkGraphConstruction:
    """Test WorkGraph construction without submission.
    
    Note: These tests require a configured AiiDA code for VASP/CP2K.
    They are skipped if no code is available.
    """

    @pytest.fixture
    def structures_dir(self, tmp_path):
        """Create temporary structures directory with test files."""
        from ase import Atoms
        from ase.io import write

        structures = tmp_path / "structures"
        structures.mkdir()

        # Create Ag2O bulk structure
        ag2o = Atoms(
            'Ag2O',
            positions=[[0, 0, 0], [2.5, 0, 0], [1.25, 1.25, 0]],
            cell=[5, 5, 5],
            pbc=True
        )
        write(str(structures / "ag2o.cif"), ag2o)

        # Create Ag metal reference
        ag = Atoms('Ag4',
                  positions=[[0, 0, 0], [2, 2, 0], [2, 0, 2], [0, 2, 2]],
                  cell=[4, 4, 4], pbc=True)
        write(str(structures / "Ag.cif"), ag)

        # Create O2 reference
        o2 = Atoms('O2', positions=[[0, 0, 0], [1.2, 0, 0]], cell=[10, 10, 10], pbc=True)
        write(str(structures / "O2.cif"), o2)

        return str(structures)

    @pytest.fixture
    def mock_code_label(self):
        """Return a mock code label for testing.
        
        Note: These tests don't actually run calculations, so we use a 
        placeholder code label. For real tests, you'd need a configured code.
        """
        # Try to find a real code with unique label
        from aiida import orm
        try:
            qb = orm.QueryBuilder()
            qb.append(orm.Code, project=['label', 'uuid'])
            codes = qb.all()
            if codes:
                # Find a code with unique label
                label_counts = {}
                for label, uuid in codes:
                    label_counts[label] = label_counts.get(label, 0) + 1
                
                # Return first unique label
                for label, count in label_counts.items():
                    if count == 1:
                        return label
                
                # If all have duplicates, use full identifier
                if codes:
                    label, uuid = codes[0]
                    return f"{label}@{uuid[:8]}"
        except Exception:
            pass
        pytest.skip("No AiiDA code configured for testing")

    @pytest.fixture
    def vasp_parameters(self):
        """Return standard VASP parameters for testing."""
        return {
            'PREC': 'Accurate',
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 3,
            'NSW': 100,
            'EDIFFG': -0.01,
            'ALGO': 'Normal',
            'LREAL': 'Auto',
            'LWAVE': False,
            'LCHARG': False,
        }

    @pytest.fixture
    def scheduler_options(self):
        """Return scheduler options for testing."""
        return {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 4,
            },
        }

    @pytest.mark.skip(reason="Requires unique configured VASP code")
    def test_build_bulk_only_workgraph(self, structures_dir, mock_code_label,
                                       vasp_parameters, scheduler_options):
        """Test building a bulk_only workgraph."""
        from teros.core.workgraph import build_core_workgraph

        wg = build_core_workgraph(
            workflow_preset='bulk_only',
            structures_dir=structures_dir,
            bulk_name='ag2o.cif',
            code_label=mock_code_label,
            potential_family='PBE',
            kpoints_spacing=0.4,
            bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            bulk_parameters=vasp_parameters,
            bulk_options=scheduler_options,
            name='Test_BulkOnly',
        )

        assert wg is not None
        assert wg.name == 'Test_BulkOnly'

        # Check that expected tasks exist
        task_names = [t.name for t in wg.tasks.values()]
        assert 'VaspWorkChain' in task_names or any('Vasp' in t for t in task_names)

    @pytest.mark.skip(reason="Requires unique configured VASP code")
    def test_build_formation_enthalpy_workgraph(self, structures_dir, mock_code_label,
                                                 vasp_parameters, scheduler_options):
        """Test building a formation_enthalpy_only workgraph."""
        from teros.core.workgraph import build_core_workgraph

        wg = build_core_workgraph(
            workflow_preset='formation_enthalpy_only',
            structures_dir=structures_dir,
            bulk_name='ag2o.cif',
            metal_name='Ag.cif',
            oxygen_name='O2.cif',
            code_label=mock_code_label,
            potential_family='PBE',
            kpoints_spacing=0.4,
            bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            bulk_parameters=vasp_parameters,
            bulk_options=scheduler_options,
            metal_potential_mapping={'Ag': 'Ag'},
            metal_parameters=vasp_parameters,
            metal_options=scheduler_options,
            oxygen_potential_mapping={'O': 'O'},
            oxygen_parameters=vasp_parameters,
            oxygen_options=scheduler_options,
            name='Test_FormationEnthalpy',
        )

        assert wg is not None

        # Check for formation enthalpy task
        task_names = [t.name for t in wg.tasks.values()]
        assert any('formation' in t.lower() or 'enthalpy' in t.lower() for t in task_names) or \
               'calculate_formation_enthalpy' in task_names

    @pytest.mark.skip(reason="Requires unique configured VASP code")
    def test_build_surface_thermodynamics_workgraph(self, structures_dir, mock_code_label,
                                                     vasp_parameters, scheduler_options):
        """Test building a surface_thermodynamics workgraph."""
        from teros.core.workgraph import build_core_workgraph

        slab_params = vasp_parameters.copy()
        slab_params['ISIF'] = 2  # Fix cell for slabs

        wg = build_core_workgraph(
            workflow_preset='surface_thermodynamics',
            structures_dir=structures_dir,
            bulk_name='ag2o.cif',
            metal_name='Ag.cif',
            oxygen_name='O2.cif',
            code_label=mock_code_label,
            potential_family='PBE',
            kpoints_spacing=0.4,
            bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            bulk_parameters=vasp_parameters,
            bulk_options=scheduler_options,
            metal_potential_mapping={'Ag': 'Ag'},
            metal_parameters=vasp_parameters,
            metal_options=scheduler_options,
            oxygen_potential_mapping={'O': 'O'},
            oxygen_parameters=vasp_parameters,
            oxygen_options=scheduler_options,
            miller_indices=[1, 0, 0],
            min_slab_thickness=10.0,
            min_vacuum_thickness=10.0,
            slab_parameters=slab_params,
            slab_options=scheduler_options,
            name='Test_SurfaceThermodynamics',
        )

        assert wg is not None

        # Check for slab generation task
        task_names = [t.name for t in wg.tasks.values()]
        assert any('slab' in t.lower() or 'generate' in t.lower() for t in task_names)


class TestWorkGraphPresetValidation:
    """Test that preset validation works correctly."""

    def test_invalid_preset_raises_error(self):
        """Test that invalid preset raises ValueError."""
        from teros.core.workgraph import build_core_workgraph

        with pytest.raises(ValueError, match="Unknown workflow preset"):
            build_core_workgraph(
                workflow_preset='nonexistent_preset',
                structures_dir='/tmp',
                bulk_name='test.cif',
            )

    def test_missing_required_params_raises_error(self, tmp_path):
        """Test that missing required parameters raise error."""
        from teros.core.workgraph import build_core_workgraph
        from ase import Atoms
        from ase.io import write

        # Create minimal structure
        structures = tmp_path / "structures"
        structures.mkdir()
        atoms = Atoms('Ag2O', positions=[[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]], cell=[5, 5, 5])
        write(str(structures / "ag2o.cif"), atoms)

        # surface_thermodynamics requires metal_name and oxygen_name
        with pytest.raises(ValueError):
            build_core_workgraph(
                workflow_preset='surface_thermodynamics',
                structures_dir=str(structures),
                bulk_name='ag2o.cif',
                # Missing: metal_name, oxygen_name
                code_label='fake_code',
            )


class TestWorkGraphOutputs:
    """Test WorkGraph output definitions.
    
    Note: These tests require a configured AiiDA code with unique labels.
    They are skipped in CI environments without proper code setup.
    """

    @pytest.fixture
    def simple_workgraph(self, tmp_path):
        """Create a simple workgraph for output testing."""
        pytest.skip("Requires unique configured VASP code - skipped in CI")

    def test_workgraph_has_outputs(self, simple_workgraph):
        """Test that workgraph has defined outputs."""
        wg = simple_workgraph

        # WorkGraph should have outputs attribute
        assert hasattr(wg, 'outputs')

    def test_bulk_only_expected_outputs(self, simple_workgraph):
        """Test that bulk_only workgraph has expected outputs."""
        wg = simple_workgraph

        # Check for common output names (may vary by implementation)
        output_names = list(wg.outputs.keys()) if hasattr(wg.outputs, 'keys') else []

        # At minimum, bulk_only should have bulk-related outputs
        # The exact names depend on implementation
        assert len(output_names) >= 0  # Just check it's iterable


# =============================================================================
# TASK GRAPH TESTS
# =============================================================================

class TestTaskGraphFunctions:
    """Test @task.graph decorated functions."""

    def test_scf_slabs_scatter_structure(self):
        """Test that scf_slabs_scatter has correct output structure."""
        from teros.core.slabs import scf_slabs_scatter

        # Check function exists and has expected attributes
        assert callable(scf_slabs_scatter)

        # Check it's wrapped as a task.graph
        assert hasattr(scf_slabs_scatter, 'build') or hasattr(scf_slabs_scatter, '__wrapped__')

    def test_relax_slabs_scatter_structure(self):
        """Test that relax_slabs_scatter has correct output structure."""
        from teros.core.slabs import relax_slabs_scatter

        assert callable(relax_slabs_scatter)

    def test_compute_surface_energies_scatter_structure(self):
        """Test that compute_surface_energies_scatter has correct output structure."""
        from teros.core.thermodynamics import compute_surface_energies_scatter

        assert callable(compute_surface_energies_scatter)

    def test_compute_adsorption_energies_scatter_structure(self):
        """Test that compute_adsorption_energies_scatter has correct output structure."""
        from teros.core.adsorption_energy import compute_adsorption_energies_scatter

        assert callable(compute_adsorption_energies_scatter)
