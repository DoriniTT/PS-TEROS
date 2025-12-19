"""
Workflow Preset Tests

These tests verify the workflow preset system works correctly.
"""

import pytest


class TestWorkflowPresetsModule:
    """Test the workflow_presets module."""

    def test_presets_module_imports(self):
        """Test that workflow_presets module imports successfully."""
        from teros.core.workflow_presets import (
            WORKFLOW_PRESETS,
            DEFAULT_PRESET,
            list_workflow_presets,
            get_preset_config,
            get_preset_summary,
            resolve_preset,
            validate_preset_inputs,
            validate_flag_dependencies,
        )

        assert WORKFLOW_PRESETS is not None
        assert DEFAULT_PRESET is not None

    def test_default_preset_exists(self):
        """Test that default preset is defined and valid."""
        from teros.core.workflow_presets import WORKFLOW_PRESETS, DEFAULT_PRESET

        assert DEFAULT_PRESET in WORKFLOW_PRESETS
        assert DEFAULT_PRESET == 'surface_thermodynamics'

    def test_all_expected_presets_exist(self, workflow_presets_list):
        """Test that all expected presets are defined."""
        from teros.core.workflow_presets import WORKFLOW_PRESETS

        for preset_name in workflow_presets_list:
            assert preset_name in WORKFLOW_PRESETS, f"Missing preset: {preset_name}"

    def test_preset_structure(self):
        """Test that each preset has required keys."""
        from teros.core.workflow_presets import WORKFLOW_PRESETS

        required_keys = ['name', 'description', 'flags', 'requires', 'dependencies', 'use_cases']
        flag_keys = [
            'relax_slabs',
            'compute_thermodynamics',
            'compute_cleavage',
            'compute_relaxation_energy',
            'compute_electronic_properties_bulk',
            'compute_electronic_properties_slabs',
            'run_aimd',
            'run_adsorption_energy',
        ]

        for preset_name, preset_config in WORKFLOW_PRESETS.items():
            # Check top-level keys
            for key in required_keys:
                assert key in preset_config, f"Preset '{preset_name}' missing key: {key}"

            # Check flags
            for flag in flag_keys:
                assert flag in preset_config['flags'], \
                    f"Preset '{preset_name}' missing flag: {flag}"

            # Check requires structure
            assert 'parameters' in preset_config['requires']
            assert 'optional' in preset_config['requires']

    def test_preset_flags_are_booleans(self):
        """Test that all flags are boolean values."""
        from teros.core.workflow_presets import WORKFLOW_PRESETS

        for preset_name, preset_config in WORKFLOW_PRESETS.items():
            for flag_name, flag_value in preset_config['flags'].items():
                assert isinstance(flag_value, bool), \
                    f"Preset '{preset_name}' flag '{flag_name}' is not boolean: {type(flag_value)}"


class TestResolvePreset:
    """Test preset resolution logic."""

    def test_resolve_default_preset(self):
        """Test resolving with no arguments uses default."""
        from teros.core.workflow_presets import resolve_preset, DEFAULT_PRESET

        preset_name, flags = resolve_preset()

        assert preset_name == DEFAULT_PRESET
        assert isinstance(flags, dict)

    def test_resolve_specific_preset(self):
        """Test resolving a specific preset."""
        from teros.core.workflow_presets import resolve_preset

        preset_name, flags = resolve_preset('bulk_only')

        assert preset_name == 'bulk_only'
        assert flags['relax_slabs'] is False
        assert flags['compute_thermodynamics'] is False

    def test_resolve_with_override(self):
        """Test that explicit flags override preset defaults."""
        from teros.core.workflow_presets import resolve_preset

        # surface_thermodynamics has compute_cleavage=False by default
        preset_name, flags = resolve_preset(
            'surface_thermodynamics',
            compute_cleavage=True
        )

        assert preset_name == 'surface_thermodynamics'
        assert flags['compute_cleavage'] is True

    def test_resolve_invalid_preset_raises(self):
        """Test that invalid preset name raises ValueError."""
        from teros.core.workflow_presets import resolve_preset

        with pytest.raises(ValueError, match="Unknown workflow preset"):
            resolve_preset('nonexistent_preset')

    def test_resolve_none_overrides_ignored(self):
        """Test that None values don't override preset defaults."""
        from teros.core.workflow_presets import resolve_preset, WORKFLOW_PRESETS

        preset_name, flags = resolve_preset(
            'surface_thermodynamics',
            relax_slabs=None,  # Should not override
            compute_cleavage=None,  # Should not override
        )

        # Should match preset defaults
        expected_flags = WORKFLOW_PRESETS['surface_thermodynamics']['flags']
        assert flags['relax_slabs'] == expected_flags['relax_slabs']
        assert flags['compute_cleavage'] == expected_flags['compute_cleavage']


class TestGetPresetConfig:
    """Test get_preset_config function."""

    def test_get_valid_preset(self):
        """Test getting a valid preset config."""
        from teros.core.workflow_presets import get_preset_config

        config = get_preset_config('bulk_only')

        assert config['name'] == 'bulk_only'
        assert 'flags' in config
        assert 'description' in config

    def test_get_invalid_preset_raises(self):
        """Test that invalid preset raises ValueError."""
        from teros.core.workflow_presets import get_preset_config

        with pytest.raises(ValueError, match="Unknown workflow preset"):
            get_preset_config('nonexistent_preset')

    def test_config_is_copy(self):
        """Test that returned config is a copy, not original."""
        from teros.core.workflow_presets import get_preset_config, WORKFLOW_PRESETS

        config = get_preset_config('bulk_only')
        config['name'] = 'modified'

        # Original should be unchanged
        assert WORKFLOW_PRESETS['bulk_only']['name'] == 'bulk_only'


class TestValidatePresetInputs:
    """Test preset input validation."""

    def test_valid_inputs(self):
        """Test validation with all required inputs provided."""
        from teros.core.workflow_presets import validate_preset_inputs

        errors = validate_preset_inputs(
            'surface_thermodynamics',
            metal_name='Ag.cif',
            oxygen_name='O2.cif',
        )

        assert len(errors) == 0

    def test_missing_required_input(self):
        """Test validation with missing required input."""
        from teros.core.workflow_presets import validate_preset_inputs

        errors = validate_preset_inputs(
            'surface_thermodynamics',
            metal_name='Ag.cif',
            oxygen_name=None,  # Missing!
        )

        assert len(errors) > 0
        assert any('oxygen_name' in e for e in errors)

    def test_bulk_only_no_requirements(self):
        """Test that bulk_only has no required parameters."""
        from teros.core.workflow_presets import validate_preset_inputs

        errors = validate_preset_inputs('bulk_only')

        assert len(errors) == 0


class TestValidateFlagDependencies:
    """Test flag dependency validation."""

    def test_thermodynamics_without_references(self):
        """Test warning when thermodynamics enabled without references."""
        from teros.core.workflow_presets import validate_flag_dependencies

        flags = {'compute_thermodynamics': True}
        messages = validate_flag_dependencies(flags, metal_name=None, oxygen_name=None)

        assert len(messages) > 0
        assert any('metal_name' in m or 'oxygen_name' in m for m in messages)

    def test_cleavage_without_relax(self):
        """Test warning when cleavage enabled without slab relaxation."""
        from teros.core.workflow_presets import validate_flag_dependencies

        flags = {'compute_cleavage': True, 'relax_slabs': False}
        messages = validate_flag_dependencies(flags)

        assert len(messages) > 0
        assert any('relax_slabs' in m for m in messages)

    def test_relaxation_energy_without_relax(self):
        """Test warning when relaxation energy enabled without slab relaxation."""
        from teros.core.workflow_presets import validate_flag_dependencies

        flags = {'compute_relaxation_energy': True, 'relax_slabs': False}
        messages = validate_flag_dependencies(flags)

        assert len(messages) > 0
        assert any('relax_slabs' in m for m in messages)

    def test_valid_configuration(self):
        """Test that valid configuration produces no warnings."""
        from teros.core.workflow_presets import validate_flag_dependencies

        flags = {
            'relax_slabs': True,
            'compute_thermodynamics': True,
            'compute_cleavage': True,
            'compute_relaxation_energy': True,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        }
        messages = validate_flag_dependencies(
            flags,
            metal_name='Ag.cif',
            oxygen_name='O2.cif',
            miller_indices=[1, 0, 0],
        )

        # Filter out warnings (keep only errors)
        errors = [m for m in messages if m.startswith('ERROR')]
        assert len(errors) == 0


class TestGetPresetSummary:
    """Test preset summary generation."""

    def test_summary_contains_name(self):
        """Test that summary contains preset name."""
        from teros.core.workflow_presets import get_preset_summary

        summary = get_preset_summary('bulk_only')

        assert 'bulk_only' in summary

    def test_summary_contains_flags(self):
        """Test that summary contains flag information."""
        from teros.core.workflow_presets import get_preset_summary

        summary = get_preset_summary('surface_thermodynamics')

        assert 'relax_slabs' in summary
        assert 'compute_thermodynamics' in summary

    def test_invalid_preset_returns_error_message(self):
        """Test that invalid preset returns error message."""
        from teros.core.workflow_presets import get_preset_summary

        summary = get_preset_summary('nonexistent')

        assert 'Unknown preset' in summary
