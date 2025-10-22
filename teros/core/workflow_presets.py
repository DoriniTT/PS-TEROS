"""
PS-TEROS Workflow Preset System

This module provides a three-tier workflow system for PS-TEROS:
1. Named workflow presets (high-level convenience)
2. Individual component flags (fine-grained control)
3. Automatic dependency resolution (smart defaults)

The preset system simplifies common workflows while maintaining full flexibility
for advanced users.
"""

import warnings
from typing import Dict, List, Optional, Set, Tuple


# ============================================================================
# PRESET DEFINITIONS
# ============================================================================

WORKFLOW_PRESETS = {
    'surface_thermodynamics': {
        'name': 'surface_thermodynamics',
        'description': 'Complete surface thermodynamics workflow with relaxation (cleavage/relaxation optional)',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': True,
            'compute_cleavage': False,  # CHANGED: Now optional
            'compute_relaxation_energy': False,  # CHANGED: Now optional
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['metal_name', 'oxygen_name'],
            'optional': ['nonmetal_name', 'miller_indices'],
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy', 'slabs'],
        'use_cases': [
            'Surface energy calculations',
            'Pourbaix-like stability analysis',
            'Complete thermodynamic characterization',
        ],
    },
    
    'surface_thermodynamics_unrelaxed': {
        'name': 'surface_thermodynamics_unrelaxed',
        'description': 'Surface thermodynamics with unrelaxed slabs only',
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': True,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['metal_name', 'oxygen_name'],
            'optional': ['nonmetal_name', 'miller_indices'],
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy', 'slabs'],
        'use_cases': [
            'Quick surface energy screening',
            'Testing slab terminations',
            'Low-cost initial assessment',
        ],
    },
    
    'cleavage_only': {
        'name': 'cleavage_only',
        'description': 'Cleavage energy calculations only',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': True,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': [],
            'optional': ['miller_indices'],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Cleavage energy analysis',
            'Material brittleness assessment',
            'Surface creation energy',
        ],
    },
    
    'relaxation_energy_only': {
        'name': 'relaxation_energy_only',
        'description': 'Relaxation energy calculations only',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': True,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': [],
            'optional': ['miller_indices'],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Surface reconstruction analysis',
            'Relaxation energy comparison',
            'Convergence studies',
        ],
    },
    
    'bulk_only': {
        'name': 'bulk_only',
        'description': 'Bulk relaxation only (no surfaces)',
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': [],
            'optional': [],
        },
        'dependencies': ['bulk'],
        'use_cases': [
            'Bulk structure optimization',
            'Testing calculation parameters',
            'Initial structure validation',
        ],
    },
    
    'formation_enthalpy_only': {
        'name': 'formation_enthalpy_only',
        'description': 'Formation enthalpy calculation without surfaces',
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['metal_name', 'oxygen_name'],
            'optional': ['nonmetal_name'],
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy'],
        'use_cases': [
            'Stability analysis',
            'Phase diagram construction',
            'Thermodynamic validation',
        ],
    },
    
    'electronic_structure_bulk_only': {
        'name': 'electronic_structure_bulk_only',
        'description': 'Electronic structure (DOS/bands) for bulk only',
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': True,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['bands_parameters'],
            'optional': ['bands_options', 'band_settings'],
        },
        'dependencies': ['bulk', 'electronic_properties'],
        'use_cases': [
            'Band structure analysis',
            'Density of states calculation',
            'Electronic property characterization',
        ],
    },
    
    'electronic_structure_slabs_only': {
        'name': 'electronic_structure_slabs_only',
        'description': 'Electronic properties (DOS and bands) for slabs only',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': True,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['slab_bands_parameters', 'slab_band_settings'],
            'optional': ['slab_bands_options', 'slab_electronic_properties', 'miller_indices'],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Slab band structure calculation',
            'Slab density of states (DOS)',
            'Surface electronic property analysis',
        ],
    },
    
    'electronic_structure_bulk_and_slabs': {
        'name': 'electronic_structure_bulk_and_slabs',
        'description': 'Electronic properties (DOS and bands) for both bulk and slabs',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': True,
            'compute_electronic_properties_slabs': True,
            'run_aimd': False,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['bands_parameters', 'band_settings', 'slab_bands_parameters', 'slab_band_settings'],
            'optional': ['bands_options', 'slab_bands_options', 'slab_electronic_properties', 'miller_indices'],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Complete electronic structure analysis',
            'Bulk vs surface electronic comparison',
            'Comprehensive band structure study',
        ],
    },
    
    'aimd_only': {
        'name': 'aimd_only',
        'description': 'AIMD simulation on slabs only (no relax_slabs)',
        'flags': {
            'relax_slabs': False,  # CHANGED: Don't relax slabs, just run AIMD
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': True,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': ['aimd_sequence', 'aimd_parameters'],
            'optional': ['aimd_options', 'aimd_potential_mapping', 'aimd_kpoints_spacing', 'miller_indices'],
        },
        'dependencies': ['bulk', 'slabs', 'aimd'],
        'use_cases': [
            'Molecular dynamics simulation',
            'Temperature effects',
            'Dynamic property analysis',
        ],
    },

    'adsorption_energy': {
        'name': 'adsorption_energy',
        'description': 'Adsorption energy calculation with automatic structure separation',
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
            'run_adsorption_energy': True,
        },
        'requires': {
            'parameters': ['adsorption_structures', 'adsorption_formulas'],
            'optional': ['adsorption_parameters', 'adsorption_options', 'adsorption_potential_mapping', 'adsorption_kpoints_spacing'],
        },
        'dependencies': ['adsorption'],
        'use_cases': [
            'Adsorption energy calculations',
            'Adsorbate binding strength analysis',
            'Catalytic activity assessment',
        ],
    },

    'comprehensive': {
        'name': 'comprehensive',
        'description': 'Complete analysis: thermodynamics + electronic properties + AIMD',
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': True,
            'compute_cleavage': True,
            'compute_relaxation_energy': True,
            'compute_electronic_properties_bulk': True,
            'compute_electronic_properties_slabs': True,
            'run_aimd': True,
            'run_adsorption_energy': False,
        },
        'requires': {
            'parameters': [
                'metal_name', 'oxygen_name',
                'bands_parameters',
                'aimd_sequence', 'aimd_parameters',
                'slab_bands_parameters',
            ],
            'optional': [
                'nonmetal_name', 'miller_indices',
                'bands_options', 'band_settings',
                'aimd_options', 'aimd_potential_mapping', 'aimd_kpoints_spacing',
                'slab_bands_options', 'slab_band_settings', 'slab_electronic_properties',
            ],
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy', 'slabs', 'electronic_properties', 'aimd'],
        'use_cases': [
            'Complete material characterization',
            'Publication-ready comprehensive analysis',
            'Full property investigation',
        ],
    },
}


# Default preset when none specified
DEFAULT_PRESET = 'surface_thermodynamics'


# ============================================================================
# PRESET RESOLUTION AND VALIDATION
# ============================================================================

def resolve_preset(
    workflow_preset: Optional[str] = None,
    relax_slabs: Optional[bool] = None,
    compute_thermodynamics: Optional[bool] = None,
    compute_cleavage: Optional[bool] = None,
    compute_relaxation_energy: Optional[bool] = None,
    compute_electronic_properties_bulk: Optional[bool] = None,
    compute_electronic_properties_slabs: Optional[bool] = None,
    run_aimd: Optional[bool] = None,
    run_adsorption_energy: Optional[bool] = None,
) -> Tuple[str, Dict[str, bool]]:
    """
    Resolve workflow preset and apply user overrides.

    Resolution order:
    1. Load preset (or use default)
    2. Apply user overrides for individual flags
    3. Return resolved preset name and final flag configuration

    Args:
        workflow_preset: Name of preset to use. If None, uses DEFAULT_PRESET
        relax_slabs: Override preset default
        compute_thermodynamics: Override preset default
        compute_cleavage: Override preset default
        compute_relaxation_energy: Override preset default
        compute_electronic_properties_bulk: Override preset default
        compute_electronic_properties_slabs: Override preset default
        run_aimd: Override preset default
        run_adsorption_energy: Override preset default
        
    Returns:
        Tuple of (preset_name, resolved_flags_dict)
        
    Raises:
        ValueError: If preset name is invalid
        
    Example:
        >>> preset, flags = resolve_preset('surface_thermodynamics', compute_cleavage=False)
        >>> print(preset)
        'surface_thermodynamics'
        >>> print(flags['compute_cleavage'])
        False
    """
    # Use default if no preset specified
    if workflow_preset is None:
        workflow_preset = DEFAULT_PRESET
    
    # Validate preset name
    if workflow_preset not in WORKFLOW_PRESETS:
        available = ', '.join(WORKFLOW_PRESETS.keys())
        raise ValueError(
            f"Unknown workflow preset: '{workflow_preset}'. "
            f"Available presets: {available}"
        )
    
    # Load preset configuration
    preset_config = WORKFLOW_PRESETS[workflow_preset]
    resolved_flags = preset_config['flags'].copy()
    
    # Apply user overrides (only if explicitly set, not None)
    overrides = {
        'relax_slabs': relax_slabs,
        'compute_thermodynamics': compute_thermodynamics,
        'compute_cleavage': compute_cleavage,
        'compute_relaxation_energy': compute_relaxation_energy,
        'compute_electronic_properties_bulk': compute_electronic_properties_bulk,
        'compute_electronic_properties_slabs': compute_electronic_properties_slabs,
        'run_aimd': run_aimd,
        'run_adsorption_energy': run_adsorption_energy,
    }
    
    for flag_name, flag_value in overrides.items():
        if flag_value is not None:
            resolved_flags[flag_name] = flag_value
    
    return workflow_preset, resolved_flags


def validate_preset_inputs(
    preset_name: str,
    **kwargs
) -> List[str]:
    """
    Validate that required parameters for a preset are provided.
    
    Args:
        preset_name: Name of the preset to validate
        **kwargs: All parameters passed to build_core_workgraph
        
    Returns:
        List of error messages (empty if valid)
        
    Example:
        >>> errors = validate_preset_inputs('surface_thermodynamics',
        ...                                  metal_name='Ag.cif',
        ...                                  oxygen_name='O2.cif')
        >>> if errors:
        ...     print("\\n".join(errors))
    """
    if preset_name not in WORKFLOW_PRESETS:
        return [f"Unknown preset: {preset_name}"]
    
    preset_config = WORKFLOW_PRESETS[preset_name]
    required_params = preset_config['requires']['parameters']
    
    errors = []
    for param in required_params:
        if param not in kwargs or kwargs[param] is None:
            errors.append(
                f"Preset '{preset_name}' requires parameter '{param}' but it was not provided or is None"
            )
    
    return errors


def validate_flag_dependencies(flags: Dict[str, bool], **kwargs) -> List[str]:
    """
    Validate flag dependencies and parameter requirements.
    
    Args:
        flags: Resolved flag configuration
        **kwargs: All parameters passed to build_core_workgraph
        
    Returns:
        List of warning/error messages
        
    Example:
        >>> flags = {'compute_thermodynamics': True, 'relax_slabs': False}
        >>> warnings = validate_flag_dependencies(flags, metal_name=None)
        >>> for w in warnings:
        ...     print(w)
    """
    messages = []
    
    # Check thermodynamics dependencies
    if flags.get('compute_thermodynamics', False):
        if kwargs.get('metal_name') is None or kwargs.get('oxygen_name') is None:
            messages.append(
                "WARNING: compute_thermodynamics=True requires both metal_name and oxygen_name. "
                "Thermodynamics will be skipped."
            )
    
    # Check cleavage dependencies
    if flags.get('compute_cleavage', False):
        if not flags.get('relax_slabs', False):
            messages.append(
                "WARNING: compute_cleavage=True requires relax_slabs=True. "
                "Cleavage energies will not be computed."
            )
    
    # Check relaxation energy dependencies
    if flags.get('compute_relaxation_energy', False):
        if not flags.get('relax_slabs', False):
            messages.append(
                "WARNING: compute_relaxation_energy=True requires relax_slabs=True. "
                "Relaxation energies will not be computed."
            )
    
    # Check electronic properties dependencies
    if flags.get('compute_electronic_properties_bulk', False):
        if kwargs.get('bands_parameters') is None:
            messages.append(
                "WARNING: compute_electronic_properties_bulk=True requires bands_parameters. "
                "Electronic properties will not be computed."
            )
    
    # Check slab electronic properties dependencies
    if flags.get('compute_electronic_properties_slabs', False):
        if kwargs.get('slab_bands_parameters') is None:
            messages.append(
                "WARNING: compute_electronic_properties_slabs=True requires slab_bands_parameters. "
                "Slab electronic properties will not be computed."
            )
    
    # Check AIMD dependencies
    if flags.get('run_aimd', False):
        if kwargs.get('aimd_sequence') is None:
            messages.append(
                "ERROR: run_aimd=True requires aimd_sequence parameter"
            )
        if kwargs.get('aimd_parameters') is None:
            messages.append(
                "ERROR: run_aimd=True requires aimd_parameters parameter"
            )
    
    # Check slab generation
    if flags.get('relax_slabs', False) or flags.get('run_aimd', False):
        if kwargs.get('miller_indices') is None and kwargs.get('input_slabs') is None:
            messages.append(
                "WARNING: No slabs will be generated. Provide either miller_indices or input_slabs."
            )
    
    return messages


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def list_workflow_presets() -> None:
    """
    Print all available workflow presets with descriptions.
    
    Example:
        >>> from teros.core import list_workflow_presets
        >>> list_workflow_presets()
    """
    print("\n" + "="*80)
    print("PS-TEROS WORKFLOW PRESETS")
    print("="*80 + "\n")
    
    for preset_name, config in WORKFLOW_PRESETS.items():
        print(f"📋 {preset_name}")
        print(f"   {config['description']}")
        print(f"   Use cases: {', '.join(config['use_cases'][:2])}")
        
        # Show required parameters
        if config['requires']['parameters']:
            req_params = ', '.join(config['requires']['parameters'])
            print(f"   Requires: {req_params}")
        
        print()
    
    print("="*80)
    print(f"Default preset: {DEFAULT_PRESET}")
    print("="*80 + "\n")
    print("Usage:")
    print("  wg = build_core_workgraph(")
    print("      workflow_preset='surface_thermodynamics',")
    print("      structures_dir='...',")
    print("      bulk_name='...',")
    print("      ...)")
    print()


def get_preset_config(preset_name: str) -> Dict:
    """
    Get the full configuration for a preset.
    
    Args:
        preset_name: Name of preset
        
    Returns:
        Preset configuration dictionary
        
    Raises:
        ValueError: If preset name is invalid
        
    Example:
        >>> config = get_preset_config('surface_thermodynamics')
        >>> print(config['description'])
    """
    if preset_name not in WORKFLOW_PRESETS:
        available = ', '.join(WORKFLOW_PRESETS.keys())
        raise ValueError(
            f"Unknown workflow preset: '{preset_name}'. "
            f"Available presets: {available}"
        )
    
    return WORKFLOW_PRESETS[preset_name].copy()


def get_preset_summary(preset_name: str) -> str:
    """
    Get a formatted summary of a preset configuration.
    
    Args:
        preset_name: Name of preset
        
    Returns:
        Formatted string summary
        
    Example:
        >>> print(get_preset_summary('surface_thermodynamics'))
    """
    if preset_name not in WORKFLOW_PRESETS:
        return f"Unknown preset: {preset_name}"
    
    config = WORKFLOW_PRESETS[preset_name]
    
    lines = [
        f"\n{'='*60}",
        f"Preset: {preset_name}",
        f"{'='*60}",
        f"\n{config['description']}\n",
        "Flags:",
    ]
    
    for flag, value in config['flags'].items():
        status = "✓" if value else "✗"
        lines.append(f"  {status} {flag}: {value}")
    
    if config['requires']['parameters']:
        lines.append("\nRequired Parameters:")
        for param in config['requires']['parameters']:
            lines.append(f"  • {param}")
    
    if config['requires']['optional']:
        lines.append("\nOptional Parameters:")
        for param in config['requires']['optional']:
            lines.append(f"  • {param}")
    
    lines.append("\nUse Cases:")
    for use_case in config['use_cases']:
        lines.append(f"  • {use_case}")
    
    lines.append(f"\n{'='*60}\n")
    
    return "\n".join(lines)


# ============================================================================
# DEPRECATION HELPERS
# ============================================================================

def check_old_style_api(
    workflow_preset: Optional[str],
    relax_slabs: Optional[bool],
    compute_thermodynamics: Optional[bool],
    compute_cleavage: Optional[bool],
    compute_relaxation_energy: Optional[bool],
    compute_electronic_properties_bulk: Optional[bool],
    compute_electronic_properties_slabs: Optional[bool],
    run_aimd: Optional[bool],
    run_adsorption_energy: Optional[bool],
) -> None:
    """
    Check if user is using old-style API (explicit flags without preset).

    Emit deprecation warning if old style is detected.

    Args:
        workflow_preset: Preset parameter
        All flag parameters: Individual flag values
    """
    # If preset is explicitly set, user is using new API
    if workflow_preset is not None:
        return
    
    # Check if any flags are explicitly set (not None and not default)
    flags_set = [
        relax_slabs is not None,
        compute_thermodynamics is not None,
        compute_cleavage is not None,
        compute_relaxation_energy is not None,
        compute_electronic_properties_bulk is not None,
        compute_electronic_properties_slabs is not None,
        run_aimd is not None,
        run_adsorption_energy is not None,
    ]
    
    if any(flags_set):
        warnings.warn(
            "Calling build_core_workgraph() with explicit boolean flags without "
            "specifying workflow_preset is deprecated. Please use workflow_preset "
            "parameter for better clarity and maintainability.\n"
            "Example: workflow_preset='surface_thermodynamics'\n"
            "Use list_workflow_presets() to see available presets.",
            DeprecationWarning,
            stacklevel=3
        )


# ============================================================================
# TESTING
# ============================================================================

if __name__ == '__main__':
    print("Testing workflow_presets module...\n")
    
    # Test 1: List all presets
    print("Test 1: List all presets")
    list_workflow_presets()
    
    # Test 2: Get preset config
    print("\nTest 2: Get preset config")
    config = get_preset_config('surface_thermodynamics')
    print(f"Config keys: {list(config.keys())}")
    
    # Test 3: Resolve preset with defaults
    print("\nTest 3: Resolve preset with defaults")
    preset, flags = resolve_preset('surface_thermodynamics')
    print(f"Preset: {preset}")
    print(f"Flags: {flags}")
    
    # Test 4: Resolve preset with overrides
    print("\nTest 4: Resolve preset with overrides")
    preset, flags = resolve_preset('surface_thermodynamics', compute_cleavage=False)
    print(f"compute_cleavage overridden to: {flags['compute_cleavage']}")
    
    # Test 5: Validate preset inputs
    print("\nTest 5: Validate preset inputs")
    errors = validate_preset_inputs(
        'surface_thermodynamics',
        metal_name='Ag.cif',
        oxygen_name=None  # Missing required parameter
    )
    print(f"Validation errors: {errors}")
    
    # Test 6: Get preset summary
    print("\nTest 6: Get preset summary")
    print(get_preset_summary('aimd_only'))
    
    # Test 7: Validate flag dependencies
    print("\nTest 7: Validate flag dependencies")
    flags = {'compute_thermodynamics': True}
    warnings_list = validate_flag_dependencies(flags, metal_name=None)
    for w in warnings_list:
        print(f"  {w}")
    
    print("\n✓ All tests completed")
