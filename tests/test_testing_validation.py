"""Tests for teros.core.testing.validation module.

These are Tier 1 tests - pure Python, no AiiDA required.
"""
import pytest


# =============================================================================
# Test ValidationResult
# =============================================================================

class TestValidationResult:
    """Test the ValidationResult dataclass."""

    def test_is_valid_no_errors(self):
        """Test is_valid returns True when no errors."""
        from teros.core.testing import ValidationResult

        result = ValidationResult()
        assert result.is_valid is True

    def test_is_valid_with_errors(self):
        """Test is_valid returns False when errors present."""
        from teros.core.testing import ValidationResult

        result = ValidationResult(errors=['Something wrong'])
        assert result.is_valid is False

    def test_is_valid_with_warnings_only(self):
        """Test is_valid returns True with warnings but no errors."""
        from teros.core.testing import ValidationResult

        result = ValidationResult(warnings=['A warning'])
        assert result.is_valid is True

    def test_str_representation(self):
        """Test string representation."""
        from teros.core.testing import ValidationResult

        result = ValidationResult(
            errors=['Error 1'],
            warnings=['Warning 1'],
        )
        s = str(result)
        assert 'ERRORS' in s
        assert 'Error 1' in s
        assert 'WARNINGS' in s
        assert 'Warning 1' in s


# =============================================================================
# Test validate_incar
# =============================================================================

class TestValidateIncar:
    """Test INCAR validation rules."""

    def test_valid_relaxation(self):
        """Test valid relaxation settings pass."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'encut': 520,
            'ibrion': 2,
            'nsw': 100,
            'isif': 2,
        })
        # Should not have the IBRION/NSW warning
        assert not any('IBRION' in w and 'NSW' in w for w in result.warnings)

    def test_ibrion_without_nsw(self):
        """Test warning for IBRION without NSW."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'encut': 520,
            'ibrion': 2,
            'nsw': 0,  # Problem!
        })
        assert any('IBRION' in w and 'NSW' in w for w in result.warnings)

    def test_isif_3_without_ibrion(self):
        """Test warning for ISIF=3 without relaxation."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'encut': 520,
            'isif': 3,
            'ibrion': -1,  # SCF only
        })
        assert any('ISIF=3' in w for w in result.warnings)

    def test_low_encut_warning(self):
        """Test warning for very low ENCUT."""
        from teros.core.testing import validate_incar

        result = validate_incar({'encut': 100})
        assert any('ENCUT=100' in w and 'low' in w for w in result.warnings)

    def test_high_encut_warning(self):
        """Test warning for very high ENCUT."""
        from teros.core.testing import validate_incar

        result = validate_incar({'encut': 900})
        assert any('ENCUT=900' in w and 'high' in w for w in result.warnings)

    def test_sigma_with_ismear_0(self):
        """Test warning for large SIGMA with ISMEAR=0."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'ismear': 0,
            'sigma': 0.2,  # Too large for insulators
        })
        assert any('SIGMA' in w for w in result.warnings)

    def test_ispin_magmom_consistency(self):
        """Test warning for MAGMOM without ISPIN=2."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'ispin': 1,
            'magmom': [0.6, 0.6],  # Will be ignored
        })
        assert any('MAGMOM' in w and 'ISPIN=1' in w for w in result.warnings)

    def test_case_insensitivity(self):
        """Test that uppercase keys are handled."""
        from teros.core.testing import validate_incar

        result = validate_incar({
            'ENCUT': 520,
            'IBRION': 2,
            'NSW': 0,
        })
        # Should still catch the IBRION/NSW issue
        assert any('IBRION' in w for w in result.warnings)


# =============================================================================
# Test validate_builder_inputs
# =============================================================================

class TestValidateBuilderInputs:
    """Test builder_inputs validation."""

    def test_valid_builder_inputs(self):
        """Test valid builder_inputs pass."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'parameters': {
                'incar': {
                    'encut': 520,
                    'ismear': 0,
                    'sigma': 0.05,
                }
            },
            'options': {
                'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4},
            },
            'kpoints_spacing': 0.03,
        })
        assert result.is_valid

    def test_missing_parameters(self):
        """Test error for missing parameters key."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'options': {'resources': {}},
        })
        assert not result.is_valid
        assert any('parameters' in e for e in result.errors)

    def test_missing_incar(self):
        """Test error for missing incar key."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'parameters': {},  # No incar
        })
        assert not result.is_valid
        assert any('incar' in e for e in result.errors)

    def test_coarse_kpoints_warning(self):
        """Test warning for coarse k-points spacing."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'parameters': {'incar': {'encut': 520}},
            'kpoints_spacing': 0.2,  # Very coarse
        })
        assert any('kpoints_spacing' in w and 'coarse' in w for w in result.warnings)

    def test_fine_kpoints_warning(self):
        """Test warning for very fine k-points spacing."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'parameters': {'incar': {'encut': 520}},
            'kpoints_spacing': 0.005,  # Very fine
        })
        assert any('kpoints_spacing' in w and 'fine' in w for w in result.warnings)

    def test_missing_resources_warning(self):
        """Test warning for missing resources."""
        from teros.core.testing import validate_builder_inputs

        result = validate_builder_inputs({
            'parameters': {'incar': {'encut': 520}},
            'options': {},  # No resources
        })
        assert any('resources' in w for w in result.warnings)


# =============================================================================
# Test estimate_kpoints_mesh
# =============================================================================

class TestEstimateKpointsMesh:
    """Test k-points mesh estimation."""

    def test_cubic_cell(self):
        """Test k-points for cubic cell."""
        from teros.core.testing import estimate_kpoints_mesh

        # Simple cubic with a=4 A
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        mesh = estimate_kpoints_mesh(lattice, spacing=0.5)

        # Should be roughly equal in all directions
        assert mesh[0] == mesh[1] == mesh[2]
        assert all(k >= 1 for k in mesh)

    def test_orthorhombic_cell(self):
        """Test k-points for orthorhombic cell."""
        from teros.core.testing import estimate_kpoints_mesh

        # Different lattice parameters
        lattice = [[10, 0, 0], [0, 5, 0], [0, 0, 2]]
        mesh = estimate_kpoints_mesh(lattice, spacing=0.1)

        # Shorter lattice vector -> more k-points
        assert mesh[2] > mesh[1] > mesh[0]

    def test_fine_spacing_more_kpoints(self):
        """Test finer spacing gives more k-points."""
        from teros.core.testing import estimate_kpoints_mesh

        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        mesh_coarse = estimate_kpoints_mesh(lattice, spacing=0.1)
        mesh_fine = estimate_kpoints_mesh(lattice, spacing=0.05)

        assert all(f >= c for f, c in zip(mesh_fine, mesh_coarse))


# =============================================================================
# Test INCAR file generation (pure Python)
# =============================================================================

class TestGenerateIncarFromDict:
    """Test INCAR file generation."""

    def test_basic_incar(self):
        """Test basic INCAR generation."""
        from teros.core.testing import generate_incar_from_dict

        content = generate_incar_from_dict({
            'ENCUT': 520,
            'ISMEAR': 0,
            'SIGMA': 0.05,
        })

        assert 'ENCUT = 520' in content
        assert 'ISMEAR = 0' in content
        assert 'SIGMA = 0.05' in content

    def test_boolean_values(self):
        """Test boolean conversion."""
        from teros.core.testing import generate_incar_from_dict

        content = generate_incar_from_dict({
            'LWAVE': True,
            'LCHARG': False,
        })

        assert '.TRUE.' in content
        assert '.FALSE.' in content

    def test_list_values(self):
        """Test list conversion."""
        from teros.core.testing import generate_incar_from_dict

        content = generate_incar_from_dict({
            'MAGMOM': [0.6, 0.6, -0.6, -0.6],
        })

        assert '0.6 0.6 -0.6 -0.6' in content


# =============================================================================
# Test KPOINTS file generation
# =============================================================================

class TestGenerateKpointsFromMesh:
    """Test KPOINTS file generation."""

    def test_gamma_centered(self):
        """Test Gamma-centered mesh generation."""
        from teros.core.testing import generate_kpoints_from_mesh

        content = generate_kpoints_from_mesh((4, 4, 4))

        assert 'Gamma' in content
        assert '4 4 4' in content

    def test_with_shift(self):
        """Test mesh with shift."""
        from teros.core.testing import generate_kpoints_from_mesh

        content = generate_kpoints_from_mesh((4, 4, 4), shift=(0.5, 0.5, 0.5))

        assert '0.5 0.5 0.5' in content


# =============================================================================
# Test module imports
# =============================================================================

class TestModuleImports:
    """Test that module imports work correctly."""

    def test_import_validation_result(self):
        """Test ValidationResult import."""
        from teros.core.testing import ValidationResult
        assert ValidationResult is not None

    def test_import_validate_incar(self):
        """Test validate_incar import."""
        from teros.core.testing import validate_incar
        assert callable(validate_incar)

    def test_import_validate_builder_inputs(self):
        """Test validate_builder_inputs import."""
        from teros.core.testing import validate_builder_inputs
        assert callable(validate_builder_inputs)

    def test_import_estimate_kpoints_mesh(self):
        """Test estimate_kpoints_mesh import."""
        from teros.core.testing import estimate_kpoints_mesh
        assert callable(estimate_kpoints_mesh)

    def test_import_incar_rules(self):
        """Test INCAR_RULES import."""
        from teros.core.testing import INCAR_RULES
        assert isinstance(INCAR_RULES, list)
        assert len(INCAR_RULES) > 0

    def test_import_generate_incar_from_dict(self):
        """Test generate_incar_from_dict import."""
        from teros.core.testing import generate_incar_from_dict
        assert callable(generate_incar_from_dict)

    def test_import_generate_kpoints_from_mesh(self):
        """Test generate_kpoints_from_mesh import."""
        from teros.core.testing import generate_kpoints_from_mesh
        assert callable(generate_kpoints_from_mesh)
