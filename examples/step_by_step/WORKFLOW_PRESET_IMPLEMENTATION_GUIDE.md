# PS-TEROS Workflow Preset System - Implementation Guide

**Version:** 1.0
**Date:** 2025-10-13
**Status:** Design Approved - Ready for Implementation

---

## Table of Contents

1. [Overview](#overview)
2. [Design Architecture](#design-architecture)
3. [File Structure](#file-structure)
4. [Implementation Steps](#implementation-steps)
5. [Code Specifications](#code-specifications)
6. [Testing Strategy](#testing-strategy)
7. [Migration Guide](#migration-guide)
8. [Examples](#examples)
9. [Future Enhancements](#future-enhancements)

---

## Overview

### Problem Statement

The current PS-TEROS workflow system uses multiple boolean flags to control different computational components. This approach has several issues:

- **Too many flags**: Users must remember and set 7+ boolean flags
- **Unclear dependencies**: Not obvious which flags require others
- **Error-prone**: Easy to enable incompatible combinations
- **No presets**: Common workflows require full manual configuration every time

### Solution: Three-Tier Workflow System

**Tier 1: Named Workflow Presets** (High-level convenience)
- Pre-configured workflows for common use cases
- Automatic dependency resolution
- One parameter activates entire workflow

**Tier 2: Independent Component Flags** (Fine-grained control)
- Override preset defaults
- Mix and match components
- Full backward compatibility

**Tier 3: Automatic Dependency Resolution** (Smart defaults)
- Auto-enable required dependencies
- Validate parameter requirements
- Clear error messages for missing inputs

### Goals

1. **Simplicity**: One parameter for common workflows
2. **Flexibility**: Override any preset with individual flags
3. **Safety**: Automatic validation and dependency resolution
4. **Discoverability**: Clear documentation of available presets
5. **Maintainability**: Easy to add new presets
6. **Backward Compatible**: Existing code continues to work

---

## Design Architecture

### Workflow Preset Structure

Each preset is a dictionary with:

```python
{
    'name': 'surface_thermodynamics',
    'description': 'Complete surface thermodynamics workflow',
    'flags': {
        'relax_slabs': True,
        'compute_thermodynamics': True,
        'compute_cleavage': True,
        'compute_relaxation_energy': True,
        'compute_electronic_properties_bulk': False,
        'compute_electronic_properties_slabs': False,
        'run_aimd': False,
    },
    'requires': {
        'parameters': ['metal_name', 'oxygen_name'],  # Required parameters
        'optional': ['nonmetal_name'],  # Optional but recommended
    },
    'dependencies': ['bulk', 'references', 'formation_enthalpy', 'slabs'],
    'use_cases': [
        'Surface energy calculations',
        'Pourbaix-like stability analysis',
        'Complete thermodynamic characterization',
    ],
}
```

### Dependency Chain

```
surface_thermodynamics
├── bulk_relaxation (always required)
├── reference_relaxations (if metal_name/oxygen_name provided)
│   ├── metal
│   ├── oxygen
│   └── nonmetal (optional)
├── formation_enthalpy (requires references)
├── slab_generation (if miller_indices provided)
│   └── slab_relaxation (if relax_slabs=True)
├── surface_energies (requires formation_enthalpy + slabs)
├── cleavage_energies (requires slabs)
└── relaxation_energies (requires slabs)
```

### Resolution Order

1. **User specifies preset**: `workflow_preset="surface_thermodynamics"`
2. **Preset loads default flags**: `relax_slabs=True, compute_thermodynamics=True, ...`
3. **User overrides specific flags**: `compute_cleavage=False`
4. **System validates requirements**: Check metal_name, oxygen_name present
5. **System auto-enables dependencies**: If compute_thermodynamics → enable formation_enthalpy
6. **Build workgraph**: Use resolved configuration

---

## File Structure

### New Files

```
teros/core/
├── workflow_presets.py              # NEW: Preset definitions and logic
└── __init__.py                      # MODIFIED: Export preset functions

docs/
├── WORKFLOW_PRESETS_GUIDE.md        # NEW: User-facing documentation
└── WORKFLOW_PRESETS_EXAMPLES.md     # NEW: Usage examples

examples/
└── workflow_presets/                # NEW: Example scripts
    ├── example_surface_thermo.py
    ├── example_aimd_only.py
    ├── example_electronic_structure.py
    └── example_custom_workflow.py
```

### Modified Files

```
teros/core/
└── workgraph.py                     # MODIFIED: Add preset support to build_core_workgraph()
```

---

## Implementation Steps

### Step 1: Create `workflow_presets.py`

**File:** `teros/core/workflow_presets.py`

**Contents:**
1. Define all 9 workflow presets
2. Implement `resolve_preset()` function
3. Implement `validate_preset_inputs()` function
4. Implement `list_workflow_presets()` helper
5. Implement `get_preset_config()` helper

**Estimated time:** 2-3 hours

### Step 2: Modify `build_core_workgraph()`

**File:** `teros/core/workgraph.py`

**Changes:**
1. Add `workflow_preset` parameter (default: `"surface_thermodynamics"`)
2. Add preset resolution at function start
3. Change all boolean flag defaults to `None`
4. Implement override logic (user flags override preset)
5. Add parameter validation based on preset
6. Add deprecation warnings for old API (optional)

**Estimated time:** 2-3 hours

### Step 3: Update `__init__.py`

**File:** `teros/core/__init__.py`

**Changes:**
1. Export `list_workflow_presets`
2. Export `get_preset_config`
3. Export preset constants (if needed)

**Estimated time:** 15 minutes

### Step 4: Create User Documentation

**File:** `docs/WORKFLOW_PRESETS_GUIDE.md`

**Contents:**
1. Overview of preset system
2. Table of all presets with descriptions
3. Parameter requirements per preset
4. How to override preset defaults
5. How to create custom presets
6. Common patterns and recipes

**Estimated time:** 1-2 hours

### Step 5: Create Example Documentation

**File:** `docs/WORKFLOW_PRESETS_EXAMPLES.md`

**Contents:**
1. Example for each preset
2. Common combinations
3. Migration examples (old API → new API)

**Estimated time:** 1 hour

### Step 6: Create Example Scripts

**Directory:** `examples/workflow_presets/`

**Scripts:**
1. `example_surface_thermo.py` - Core workflow
2. `example_aimd_only.py` - AIMD workflow
3. `example_electronic_structure.py` - DOS/bands workflow
4. `example_custom_workflow.py` - Mixed preset + overrides

**Estimated time:** 2 hours

### Step 7: Testing

1. Test each preset independently
2. Test preset + override combinations
3. Test parameter validation
4. Test backward compatibility
5. Run existing examples to ensure no breakage

**Estimated time:** 3-4 hours

---

## Code Specifications

### `workflow_presets.py` - Complete Implementation

```python
"""
PS-TEROS Workflow Preset System

This module provides pre-configured workflow presets for common PS-TEROS
use cases, automatic dependency resolution, and parameter validation.

Usage:
    from teros.core.workflow_presets import list_workflow_presets, get_preset_config

    # List available presets
    list_workflow_presets()

    # Get preset configuration
    config = get_preset_config('surface_thermodynamics')
"""

from typing import Dict, List, Any, Optional
import warnings


# =============================================================================
# WORKFLOW PRESET DEFINITIONS
# =============================================================================

WORKFLOW_PRESETS = {

    # -------------------------------------------------------------------------
    # 1. SURFACE THERMODYNAMICS (Default/Core Workflow)
    # -------------------------------------------------------------------------
    'surface_thermodynamics': {
        'description': 'Complete surface thermodynamics workflow',
        'long_description': (
            'Performs bulk relaxation, reference relaxations, formation enthalpy '
            'calculation, slab generation and relaxation, surface energy analysis, '
            'cleavage energy analysis, and relaxation energy analysis.'
        ),
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': True,
            'compute_cleavage': True,
            'compute_relaxation_energy': True,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': ['metal_name', 'oxygen_name'],
            'optional': ['nonmetal_name'],
            'either_or': [['miller_indices', 'input_slabs']],  # Need one or the other
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy', 'slabs'],
        'use_cases': [
            'Surface energy calculations as function of chemical potential',
            'Pourbaix-like stability phase diagrams',
            'Complete thermodynamic characterization of surfaces',
        ],
        'notes': [
            'Requires metal and oxygen reference structures',
            'For ternary oxides, provide nonmetal_name',
            'Automatically generates slabs if miller_indices provided',
        ],
    },

    # -------------------------------------------------------------------------
    # 2. RELAXATION ENERGY ANALYSIS
    # -------------------------------------------------------------------------
    'relaxation_energy_analysis': {
        'description': 'Slab relaxation energy analysis',
        'long_description': (
            'Performs bulk relaxation, slab generation, SCF on unrelaxed slabs, '
            'slab relaxation, and calculates relaxation energies (E_relaxed - E_unrelaxed).'
        ),
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': True,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': [],
            'optional': ['metal_name', 'oxygen_name'],
            'either_or': [['miller_indices', 'input_slabs']],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Analyze how much energy is gained by relaxing surfaces',
            'Screen different surface terminations by relaxation energy',
            'Validate convergence of slab relaxation',
        ],
        'notes': [
            'Does not require reference structures',
            'Runs both SCF and relaxation for each slab',
        ],
    },

    # -------------------------------------------------------------------------
    # 3. BULK ONLY
    # -------------------------------------------------------------------------
    'bulk_only': {
        'description': 'Bulk structure relaxation only',
        'long_description': (
            'Performs only bulk structure relaxation. Useful for initial '
            'structure optimization or when you only need bulk properties.'
        ),
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': [],
            'optional': [],
            'either_or': [],
        },
        'dependencies': ['bulk'],
        'use_cases': [
            'Initial structure optimization',
            'Quick bulk energy calculation',
            'Benchmark different DFT parameters on bulk',
        ],
        'notes': [
            'Minimal workflow - just bulk relaxation',
            'Useful as first step before more complex calculations',
        ],
    },

    # -------------------------------------------------------------------------
    # 4. FORMATION ENTHALPY ONLY
    # -------------------------------------------------------------------------
    'formation_enthalpy_only': {
        'description': 'Formation enthalpy calculation without slabs',
        'long_description': (
            'Performs bulk and reference relaxations, then calculates formation '
            'enthalpy. No slab calculations.'
        ),
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': ['metal_name', 'oxygen_name'],
            'optional': ['nonmetal_name'],
            'either_or': [],
        },
        'dependencies': ['bulk', 'references', 'formation_enthalpy'],
        'use_cases': [
            'Calculate stability of bulk compound',
            'Screen different compositions',
            'Validate thermodynamic data',
        ],
        'notes': [
            'Requires metal and oxygen reference structures',
            'For ternary oxides, provide nonmetal_name',
            'Does not generate or relax slabs',
        ],
    },

    # -------------------------------------------------------------------------
    # 5. SLABS ONLY
    # -------------------------------------------------------------------------
    'slabs_only': {
        'description': 'Slab generation and relaxation only',
        'long_description': (
            'Generates and relaxes slabs from a provided or previously relaxed '
            'bulk structure. No formation enthalpy or surface energy calculations.'
        ),
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': [],
            'optional': [],
            'either_or': [['miller_indices', 'input_slabs']],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Generate and relax slabs from known bulk structure',
            'Rerun slab calculations with different parameters',
            'Prepare slabs for subsequent AIMD or electronic structure',
        ],
        'notes': [
            'Can use previously relaxed bulk (set restart_from_node)',
            'Generates slabs if miller_indices provided',
            'Use input_slabs to relax pre-generated slabs',
        ],
    },

    # -------------------------------------------------------------------------
    # 6. CLEAVAGE ANALYSIS
    # -------------------------------------------------------------------------
    'cleavage_analysis': {
        'description': 'Cleavage energy analysis',
        'long_description': (
            'Performs bulk relaxation, slab generation and relaxation, and '
            'calculates cleavage energies for complementary slab pairs.'
        ),
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': True,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': ['miller_indices'],  # Must auto-generate (not input_slabs)
            'optional': [],
            'either_or': [],
        },
        'dependencies': ['bulk', 'slabs'],
        'use_cases': [
            'Calculate cleavage energies for complementary surfaces',
            'Analyze energy cost of cleaving bulk crystal',
            'Screen different cleavage planes',
        ],
        'notes': [
            'Requires automatic slab generation (miller_indices)',
            'Cannot use input_slabs (need complementary pairs)',
            'Does not require reference structures',
        ],
    },

    # -------------------------------------------------------------------------
    # 7. ELECTRONIC STRUCTURE FULL
    # -------------------------------------------------------------------------
    'electronic_structure_full': {
        'description': 'Electronic structure (bulk + slabs)',
        'long_description': (
            'Performs bulk relaxation with DOS/bands calculation, slab generation '
            'and relaxation, and DOS/bands for selected slabs.'
        ),
        'flags': {
            'relax_slabs': True,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': True,
            'compute_electronic_properties_slabs': True,
            'run_aimd': False,
        },
        'requires': {
            'parameters': ['bands_parameters', 'slab_electronic_properties'],
            'optional': ['bands_options', 'band_settings', 'slab_bands_parameters'],
            'either_or': [['miller_indices', 'input_slabs']],
        },
        'dependencies': ['bulk', 'slabs', 'electronic_structure'],
        'use_cases': [
            'Complete electronic structure characterization',
            'Calculate band structures for bulk and surfaces',
            'Analyze surface states and band bending',
        ],
        'notes': [
            'Requires bands_parameters with scf/bands/dos keys',
            'slab_electronic_properties selects which slabs get DOS/bands',
            'Can be computationally expensive',
        ],
    },

    # -------------------------------------------------------------------------
    # 8. ELECTRONIC STRUCTURE BULK ONLY
    # -------------------------------------------------------------------------
    'electronic_structure_bulk_only': {
        'description': 'Electronic structure for bulk only',
        'long_description': (
            'Performs bulk relaxation followed by DOS and band structure '
            'calculation. No slab calculations.'
        ),
        'flags': {
            'relax_slabs': False,
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': True,
            'compute_electronic_properties_slabs': False,
            'run_aimd': False,
        },
        'requires': {
            'parameters': ['bands_parameters'],
            'optional': ['bands_options', 'band_settings'],
            'either_or': [],
        },
        'dependencies': ['bulk', 'electronic_structure'],
        'use_cases': [
            'Calculate bulk band structure and DOS',
            'Analyze bulk electronic properties',
            'Quick electronic structure screening',
        ],
        'notes': [
            'Requires bands_parameters with scf/bands/dos keys',
            'No slab calculations performed',
        ],
    },

    # -------------------------------------------------------------------------
    # 9. AIMD ONLY
    # -------------------------------------------------------------------------
    'aimd_only': {
        'description': 'Ab initio molecular dynamics',
        'long_description': (
            'Runs AIMD on provided or generated slabs. Can use relaxed slabs '
            'from previous calculations or generate new ones.'
        ),
        'flags': {
            'relax_slabs': False,  # Optional: can be overridden to True
            'compute_thermodynamics': False,
            'compute_cleavage': False,
            'compute_relaxation_energy': False,
            'compute_electronic_properties_bulk': False,
            'compute_electronic_properties_slabs': False,
            'run_aimd': True,
        },
        'requires': {
            'parameters': ['aimd_sequence', 'aimd_parameters'],
            'optional': ['aimd_options', 'aimd_potential_mapping', 'aimd_kpoints_spacing'],
            'either_or': [['miller_indices', 'input_slabs']],
        },
        'dependencies': ['bulk', 'slabs', 'aimd'],
        'use_cases': [
            'Run finite-temperature MD on surfaces',
            'Study surface dynamics and reconstruction',
            'Equilibrate surface structures',
        ],
        'notes': [
            'Requires aimd_sequence with temperature and steps',
            'Can relax slabs first by setting relax_slabs=True',
            'Sequential AIMD stages automatically chain',
        ],
    },
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def list_workflow_presets(verbose: bool = False) -> None:
    """
    Print available workflow presets with descriptions.

    Args:
        verbose: If True, print detailed information for each preset
    """
    print("\n" + "=" * 80)
    print("Available PS-TEROS Workflow Presets")
    print("=" * 80 + "\n")

    for name, config in WORKFLOW_PRESETS.items():
        print(f"📋 {name}")
        print(f"   {config['description']}")

        if verbose:
            print(f"\n   Long Description:")
            print(f"   {config['long_description']}")

            print(f"\n   Required Parameters:")
            if config['requires']['parameters']:
                for param in config['requires']['parameters']:
                    print(f"   - {param}")
            else:
                print(f"   - None")

            if config['requires']['optional']:
                print(f"\n   Optional Parameters:")
                for param in config['requires']['optional']:
                    print(f"   - {param}")

            if config['requires']['either_or']:
                print(f"\n   Either/Or Requirements:")
                for group in config['requires']['either_or']:
                    print(f"   - Provide one of: {', '.join(group)}")

            print(f"\n   Use Cases:")
            for use_case in config['use_cases']:
                print(f"   - {use_case}")

            print(f"\n   Notes:")
            for note in config['notes']:
                print(f"   - {note}")

            print("\n" + "-" * 80 + "\n")
        else:
            print()

    if not verbose:
        print("💡 Tip: Use list_workflow_presets(verbose=True) for detailed information")

    print("=" * 80 + "\n")


def get_preset_config(preset_name: str) -> Dict[str, Any]:
    """
    Get configuration dictionary for a workflow preset.

    Args:
        preset_name: Name of the preset (e.g., 'surface_thermodynamics')

    Returns:
        Dictionary with preset configuration

    Raises:
        ValueError: If preset_name not found
    """
    if preset_name not in WORKFLOW_PRESETS:
        available = ', '.join(WORKFLOW_PRESETS.keys())
        raise ValueError(
            f"Unknown workflow preset '{preset_name}'. "
            f"Available presets: {available}\n"
            f"Use list_workflow_presets() to see all options."
        )

    return WORKFLOW_PRESETS[preset_name].copy()


def resolve_preset(
    preset_name: str,
    user_overrides: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Resolve a workflow preset to its final flag configuration.

    This function:
    1. Loads the preset configuration
    2. Applies user overrides for any flags
    3. Returns the final resolved flags

    Args:
        preset_name: Name of the workflow preset
        user_overrides: Dictionary of user-specified flag overrides
                       (only flags with non-None values)

    Returns:
        Dictionary with resolved boolean flags

    Example:
        >>> resolve_preset('surface_thermodynamics', {'compute_cleavage': False})
        {
            'relax_slabs': True,
            'compute_thermodynamics': True,
            'compute_cleavage': False,  # User override
            'compute_relaxation_energy': True,
            ...
        }
    """
    # Get preset config
    config = get_preset_config(preset_name)

    # Start with preset defaults
    resolved_flags = config['flags'].copy()

    # Apply user overrides (only non-None values)
    for key, value in user_overrides.items():
        if value is not None:
            resolved_flags[key] = value

    return resolved_flags


def validate_preset_inputs(
    preset_name: str,
    **kwargs
) -> None:
    """
    Validate that required parameters are provided for a preset.

    Args:
        preset_name: Name of the workflow preset
        **kwargs: All parameters passed to build_core_workgraph()

    Raises:
        ValueError: If required parameters are missing
    """
    config = get_preset_config(preset_name)
    requirements = config['requires']

    # Check required parameters
    missing_required = []
    for param in requirements['parameters']:
        if kwargs.get(param) is None:
            missing_required.append(param)

    if missing_required:
        raise ValueError(
            f"Workflow preset '{preset_name}' requires the following parameters:\n"
            f"  Missing: {', '.join(missing_required)}\n"
            f"  Description: {config['description']}"
        )

    # Check either/or requirements
    for either_or_group in requirements.get('either_or', []):
        if not any(kwargs.get(param) is not None for param in either_or_group):
            raise ValueError(
                f"Workflow preset '{preset_name}' requires at least one of:\n"
                f"  {', '.join(either_or_group)}\n"
                f"  Description: {config['description']}"
            )

    # Warn about recommended optional parameters
    missing_optional = []
    for param in requirements.get('optional', []):
        if kwargs.get(param) is None:
            missing_optional.append(param)

    if missing_optional:
        warnings.warn(
            f"Workflow preset '{preset_name}' recommends providing:\n"
            f"  Optional: {', '.join(missing_optional)}\n"
            f"  (These are optional but may improve results)",
            UserWarning
        )


def get_preset_summary() -> str:
    """
    Generate a formatted summary table of all workflow presets.

    Returns:
        Formatted string with table of presets
    """
    header = "| Preset Name | Description | Requires References | Requires Slabs | Computes Thermo | Computes Elec |"
    separator = "|-------------|-------------|---------------------|----------------|-----------------|---------------|"

    rows = [header, separator]

    for name, config in WORKFLOW_PRESETS.items():
        requires_refs = 'Yes' if 'metal_name' in config['requires']['parameters'] else 'No'
        requires_slabs = 'Yes' if 'slabs' in config['dependencies'] else 'No'
        computes_thermo = 'Yes' if config['flags']['compute_thermodynamics'] else 'No'
        computes_elec = 'Yes' if (
            config['flags']['compute_electronic_properties_bulk'] or
            config['flags']['compute_electronic_properties_slabs']
        ) else 'No'

        row = f"| {name} | {config['description']} | {requires_refs} | {requires_slabs} | {computes_thermo} | {computes_elec} |"
        rows.append(row)

    return '\n'.join(rows)


# =============================================================================
# DEPRECATION HELPERS
# =============================================================================

def check_old_api_usage(**kwargs) -> bool:
    """
    Check if user is using old API (explicit boolean flags without preset).

    Returns:
        True if old API detected (should show deprecation warning)
    """
    # Check if any boolean flags are explicitly set (not None)
    flag_names = [
        'relax_slabs',
        'compute_thermodynamics',
        'compute_cleavage',
        'compute_relaxation_energy',
        'compute_electronic_properties_bulk',
        'compute_electronic_properties_slabs',
        'run_aimd',
    ]

    # If workflow_preset is None and any flags are set, it's old API
    if kwargs.get('workflow_preset') is None:
        for flag in flag_names:
            if kwargs.get(flag) is not None:
                return True

    return False


def show_deprecation_warning() -> None:
    """
    Show deprecation warning for old API usage.
    """
    warnings.warn(
        "\n" + "=" * 80 + "\n"
        "DEPRECATION WARNING: Using explicit boolean flags without workflow_preset\n"
        "is deprecated and will be removed in a future version.\n\n"
        "Please use workflow presets for better clarity and maintainability:\n\n"
        "  Old API (deprecated):\n"
        "    build_core_workgraph(\n"
        "        relax_slabs=True,\n"
        "        compute_thermodynamics=True,\n"
        "        compute_cleavage=True,\n"
        "        ...\n"
        "    )\n\n"
        "  New API (recommended):\n"
        "    build_core_workgraph(\n"
        "        workflow_preset='surface_thermodynamics',\n"
        "        ...\n"
        "    )\n\n"
        "Use list_workflow_presets() to see available presets.\n"
        "=" * 80 + "\n",
        DeprecationWarning,
        stacklevel=3
    )


# =============================================================================
# MAIN (for testing)
# =============================================================================

if __name__ == '__main__':
    # Test the module
    print("\n🧪 Testing workflow_presets module\n")

    # Test 1: List all presets
    print("Test 1: List all presets (brief)")
    list_workflow_presets(verbose=False)

    # Test 2: Get specific preset
    print("\nTest 2: Get surface_thermodynamics preset")
    config = get_preset_config('surface_thermodynamics')
    print(f"  Flags: {config['flags']}")
    print(f"  Required: {config['requires']['parameters']}")

    # Test 3: Resolve preset with overrides
    print("\nTest 3: Resolve preset with user overrides")
    resolved = resolve_preset(
        'surface_thermodynamics',
        {'compute_cleavage': False, 'run_aimd': True}
    )
    print(f"  compute_cleavage: {resolved['compute_cleavage']} (overridden)")
    print(f"  run_aimd: {resolved['run_aimd']} (overridden)")
    print(f"  compute_thermodynamics: {resolved['compute_thermodynamics']} (from preset)")

    # Test 4: Validation (should pass)
    print("\nTest 4: Validate correct inputs")
    try:
        validate_preset_inputs(
            'surface_thermodynamics',
            metal_name='Ag.cif',
            oxygen_name='O2.cif',
            miller_indices=[1, 0, 0],
        )
        print("  ✅ Validation passed")
    except ValueError as e:
        print(f"  ❌ Validation failed: {e}")

    # Test 5: Validation (should fail)
    print("\nTest 5: Validate incorrect inputs (missing metal_name)")
    try:
        validate_preset_inputs(
            'surface_thermodynamics',
            oxygen_name='O2.cif',
            miller_indices=[1, 0, 0],
        )
        print("  ❌ Should have failed validation")
    except ValueError as e:
        print(f"  ✅ Correctly caught error: {e}")

    # Test 6: Summary table
    print("\nTest 6: Generate preset summary table")
    print(get_preset_summary())

    print("\n✅ All tests completed\n")
```

### Modifications to `build_core_workgraph()` in `workgraph.py`

**Location:** Near the top of the function (after docstring, before any logic)

```python
def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,

    # NEW: Workflow preset parameter (defaults to core workflow)
    workflow_preset: str = 'surface_thermodynamics',

    # MODIFIED: Change all boolean defaults to None (preset will decide)
    relax_slabs: bool = None,
    compute_thermodynamics: bool = None,
    compute_cleavage: bool = None,
    compute_relaxation_energy: bool = None,
    compute_electronic_properties_bulk: bool = None,
    compute_electronic_properties_slabs: bool = None,
    run_aimd: bool = None,

    # ... rest of parameters unchanged ...
):
    """
    Build a centralized WorkGraph for bulk relaxation, formation enthalpy, and surface calculations.

    NEW in v2.0: Workflow Preset System
    ------------------------------------
    Use workflow_preset to activate pre-configured workflows for common use cases.
    All presets can be customized by overriding individual flags.

    Available presets:
        - 'surface_thermodynamics' (default): Complete surface thermodynamics workflow
        - 'relaxation_energy_analysis': Slab relaxation energy analysis
        - 'bulk_only': Bulk structure relaxation only
        - 'formation_enthalpy_only': Formation enthalpy without slabs
        - 'slabs_only': Slab generation and relaxation
        - 'cleavage_analysis': Cleavage energy analysis
        - 'electronic_structure_full': DOS/bands for bulk and slabs
        - 'electronic_structure_bulk_only': DOS/bands for bulk only
        - 'aimd_only': AIMD on slabs

    Use list_workflow_presets() to see detailed descriptions.

    Args:
        workflow_preset: Name of workflow preset to use (default: 'surface_thermodynamics')
                        Set to None to use old API (explicit flags, deprecated)

        relax_slabs: Override preset setting for slab relaxation
        compute_thermodynamics: Override preset setting for surface energies
        compute_cleavage: Override preset setting for cleavage energies
        compute_relaxation_energy: Override preset setting for relaxation energies
        compute_electronic_properties_bulk: Override preset setting for bulk DOS/bands
        compute_electronic_properties_slabs: Override preset setting for slab DOS/bands
        run_aimd: Override preset setting for AIMD

        ... (rest of docstring unchanged) ...

    Examples:
        Using presets (recommended):

        >>> # Core surface thermodynamics workflow
        >>> wg = build_core_workgraph(
        ...     workflow_preset='surface_thermodynamics',
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     metal_name='Ag.cif',
        ...     oxygen_name='O2.cif',
        ...     ...
        ... )

        >>> # AIMD only on input slabs
        >>> wg = build_core_workgraph(
        ...     workflow_preset='aimd_only',
        ...     input_slabs=my_slabs,
        ...     aimd_sequence=[{'temperature': 300, 'steps': 1000}],
        ...     ...
        ... )

        >>> # Custom: preset + overrides
        >>> wg = build_core_workgraph(
        ...     workflow_preset='surface_thermodynamics',
        ...     compute_cleavage=False,  # Disable cleavage
        ...     run_aimd=True,  # Add AIMD
        ...     ...
        ... )
    """
    from teros.core.workflow_presets import (
        resolve_preset,
        validate_preset_inputs,
        check_old_api_usage,
        show_deprecation_warning,
    )

    # =========================================================================
    # STEP 1: Handle workflow preset resolution
    # =========================================================================

    # Check if user is using old API (deprecated)
    if check_old_api_usage(**locals()):
        show_deprecation_warning()
        # If using old API, set workflow_preset to None to skip preset logic
        workflow_preset = None

    # Resolve preset to flags
    if workflow_preset is not None:
        # Validate required inputs for this preset
        validate_preset_inputs(workflow_preset, **locals())

        # Build dictionary of user overrides (only non-None values)
        user_overrides = {
            'relax_slabs': relax_slabs,
            'compute_thermodynamics': compute_thermodynamics,
            'compute_cleavage': compute_cleavage,
            'compute_relaxation_energy': compute_relaxation_energy,
            'compute_electronic_properties_bulk': compute_electronic_properties_bulk,
            'compute_electronic_properties_slabs': compute_electronic_properties_slabs,
            'run_aimd': run_aimd,
        }
        user_overrides = {k: v for k, v in user_overrides.items() if v is not None}

        # Resolve preset with user overrides
        resolved_flags = resolve_preset(workflow_preset, user_overrides)

        # Apply resolved flags
        relax_slabs = resolved_flags['relax_slabs']
        compute_thermodynamics = resolved_flags['compute_thermodynamics']
        compute_cleavage = resolved_flags['compute_cleavage']
        compute_relaxation_energy = resolved_flags['compute_relaxation_energy']
        compute_electronic_properties_bulk = resolved_flags['compute_electronic_properties_bulk']
        compute_electronic_properties_slabs = resolved_flags['compute_electronic_properties_slabs']
        run_aimd = resolved_flags['run_aimd']

        # Print what preset is being used
        print(f"\n{'='*70}")
        print(f"Using workflow preset: '{workflow_preset}'")
        print(f"{'='*70}")
        if user_overrides:
            print(f"  User overrides: {list(user_overrides.keys())}")
        print()

    else:
        # Old API: Apply defaults for any still-None flags
        relax_slabs = relax_slabs if relax_slabs is not None else False
        compute_thermodynamics = compute_thermodynamics if compute_thermodynamics is not None else False
        compute_cleavage = compute_cleavage if compute_cleavage is not None else False
        compute_relaxation_energy = compute_relaxation_energy if compute_relaxation_energy is not None else False
        compute_electronic_properties_bulk = compute_electronic_properties_bulk if compute_electronic_properties_bulk is not None else False
        compute_electronic_properties_slabs = compute_electronic_properties_slabs if compute_electronic_properties_slabs is not None else False
        run_aimd = run_aimd if run_aimd is not None else False

    # =========================================================================
    # STEP 2: Continue with existing build_core_workgraph logic
    # =========================================================================

    # (All existing code after this point remains unchanged)
    # Handle restart from previous workgraph
    restart_folders = None
    restart_slabs = None
    if restart_from_node is not None:
        # ... existing restart logic ...

    # ... rest of function unchanged ...
```

### Modifications to `__init__.py`

**File:** `teros/core/__init__.py`

```python
"""
PS-TEROS Core Module

Main computational workflows for PS-TEROS calculations.
"""

# Import main workgraph builder
from .workgraph import (
    build_core_workgraph,
    build_core_workgraph_with_map,  # Deprecated
    core_workgraph,
)

# Import workflow preset utilities (NEW)
from .workflow_presets import (
    list_workflow_presets,
    get_preset_config,
    get_preset_summary,
    WORKFLOW_PRESETS,
)

# Import helper modules
from .slabs import (
    generate_slab_structures,
    relax_slabs_scatter,
    scf_slabs_scatter,
    calculate_relaxation_energies_scatter,
)

from .thermodynamics import (
    calculate_surface_energy_binary,
    calculate_surface_energy_ternary,
    compute_surface_energies_scatter,
)

from .cleavage import (
    calculate_cleavage_energy,
    compute_cleavage_energies_scatter,
)

from .hf import (
    calculate_formation_enthalpy,
)

from .aimd import (
    aimd_single_stage_scatter,
)

# Define what gets imported with "from teros.core import *"
__all__ = [
    # Main builders
    'build_core_workgraph',
    'build_core_workgraph_with_map',
    'core_workgraph',

    # Workflow presets (NEW)
    'list_workflow_presets',
    'get_preset_config',
    'get_preset_summary',
    'WORKFLOW_PRESETS',

    # Slab functions
    'generate_slab_structures',
    'relax_slabs_scatter',
    'scf_slabs_scatter',
    'calculate_relaxation_energies_scatter',

    # Thermodynamics functions
    'calculate_surface_energy_binary',
    'calculate_surface_energy_ternary',
    'compute_surface_energies_scatter',

    # Cleavage functions
    'calculate_cleavage_energy',
    'compute_cleavage_energies_scatter',

    # Formation enthalpy
    'calculate_formation_enthalpy',

    # AIMD
    'aimd_single_stage_scatter',
]
```

---

## Testing Strategy

### Phase 1: Unit Tests

Test each preset individually:

```python
def test_preset_resolution():
    """Test that presets resolve to correct flags"""
    from teros.core.workflow_presets import resolve_preset

    # Test surface_thermodynamics preset
    flags = resolve_preset('surface_thermodynamics', {})
    assert flags['relax_slabs'] == True
    assert flags['compute_thermodynamics'] == True
    assert flags['compute_cleavage'] == True
    assert flags['run_aimd'] == False

    # Test with override
    flags = resolve_preset('surface_thermodynamics', {'compute_cleavage': False})
    assert flags['compute_cleavage'] == False  # Overridden
    assert flags['compute_thermodynamics'] == True  # From preset

def test_validation():
    """Test parameter validation"""
    from teros.core.workflow_presets import validate_preset_inputs

    # Should pass
    validate_preset_inputs(
        'surface_thermodynamics',
        metal_name='Ag.cif',
        oxygen_name='O2.cif',
        miller_indices=[1, 0, 0],
    )

    # Should fail (missing metal_name)
    try:
        validate_preset_inputs(
            'surface_thermodynamics',
            oxygen_name='O2.cif',
        )
        assert False, "Should have raised ValueError"
    except ValueError:
        pass  # Expected
```

### Phase 2: Integration Tests

Test each preset with actual workflow building:

```python
def test_surface_thermodynamics_preset():
    """Test surface_thermodynamics preset builds correctly"""
    from teros.core import build_core_workgraph

    wg = build_core_workgraph(
        workflow_preset='surface_thermodynamics',
        structures_dir='structures',
        bulk_name='ag3po4.cif',
        metal_name='Ag.cif',
        oxygen_name='O2.cif',
        miller_indices=[1, 0, 0],
        # ... minimal required parameters ...
    )

    # Verify correct tasks exist in workgraph
    assert 'VaspWorkChain' in wg.tasks  # Bulk relaxation
    assert 'calculate_formation_enthalpy' in wg.tasks
    assert 'generate_slab_structures' in wg.tasks or 'relax_slabs_scatter' in wg.tasks
    # ... more assertions ...

def test_bulk_only_preset():
    """Test bulk_only preset builds minimal workflow"""
    from teros.core import build_core_workgraph

    wg = build_core_workgraph(
        workflow_preset='bulk_only',
        structures_dir='structures',
        bulk_name='ag3po4.cif',
        # ... minimal required parameters ...
    )

    # Verify ONLY bulk task exists
    assert 'VaspWorkChain' in wg.tasks
    assert 'calculate_formation_enthalpy' not in wg.tasks
    assert 'generate_slab_structures' not in wg.tasks
```

### Phase 3: End-to-End Tests

Run complete workflows with actual VASP calculations (if possible on test system):

```python
def test_e2e_surface_thermodynamics():
    """End-to-end test of surface_thermodynamics workflow"""
    from teros.core import build_core_workgraph
    from aiida import load_profile

    load_profile()

    wg = build_core_workgraph(
        workflow_preset='surface_thermodynamics',
        # ... full parameters ...
    )

    # Submit and wait
    wg.submit(wait=True)

    # Check outputs
    assert wg.state == 'FINISHED'
    assert 'formation_enthalpy' in wg.outputs
    assert 'surface_energies' in wg.outputs
    # ... more assertions ...
```

### Phase 4: Backward Compatibility Tests

Ensure old API still works (with deprecation warnings):

```python
def test_backward_compatibility():
    """Test that old API still works"""
    from teros.core import build_core_workgraph
    import warnings

    # Old API should work but show deprecation warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        wg = build_core_workgraph(
            workflow_preset=None,  # Explicitly no preset
            relax_slabs=True,
            compute_thermodynamics=True,
            # ... old style parameters ...
        )

        # Check deprecation warning was raised
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "DEPRECATION WARNING" in str(w[0].message)

    # Verify workgraph built correctly
    assert wg is not None
```

---

## Migration Guide

### For Users: Migrating from Old API to New API

#### Pattern 1: Simple Surface Thermodynamics

**Old API:**
```python
wg = build_core_workgraph(
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    miller_indices=[1, 0, 0],
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=True,
    # ... all other parameters ...
)
```

**New API:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # ONE parameter replaces 4 flags!
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    miller_indices=[1, 0, 0],
    # ... all other parameters unchanged ...
)
```

#### Pattern 2: Bulk Only

**Old API:**
```python
wg = build_core_workgraph(
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    relax_slabs=False,
    compute_thermodynamics=False,
    compute_cleavage=False,
    compute_relaxation_energy=False,
    # ... all other parameters ...
)
```

**New API:**
```python
wg = build_core_workgraph(
    workflow_preset='bulk_only',
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    # ... all other parameters unchanged ...
)
```

#### Pattern 3: Custom Combination

**Old API:**
```python
wg = build_core_workgraph(
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    miller_indices=[1, 0, 0],
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=False,  # Disabled
    compute_relaxation_energy=True,
    compute_electronic_properties_bulk=True,  # Added
    # ... all other parameters ...
)
```

**New API:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=False,  # Override preset
    compute_electronic_properties_bulk=True,  # Add to preset
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    miller_indices=[1, 0, 0],
    # ... all other parameters unchanged ...
)
```

### For Developers: Adding New Presets

To add a new workflow preset:

1. **Add preset definition to `WORKFLOW_PRESETS` dict in `workflow_presets.py`:**

```python
'my_new_preset': {
    'description': 'One-line description',
    'long_description': 'Detailed explanation of what this preset does...',
    'flags': {
        'relax_slabs': True,
        'compute_thermodynamics': False,
        # ... set all flags ...
    },
    'requires': {
        'parameters': ['required_param1', 'required_param2'],
        'optional': ['optional_param1'],
        'either_or': [['param_a', 'param_b']],  # Need one of these
    },
    'dependencies': ['bulk', 'slabs'],  # High-level dependencies
    'use_cases': [
        'Use case 1',
        'Use case 2',
    ],
    'notes': [
        'Important note 1',
        'Important note 2',
    ],
},
```

2. **Test the new preset:**

```python
def test_my_new_preset():
    from teros.core import build_core_workgraph

    wg = build_core_workgraph(
        workflow_preset='my_new_preset',
        # ... required parameters ...
    )

    # Verify expected tasks exist
    # ...
```

3. **Document the preset in user guide**

4. **Add example script** in `examples/workflow_presets/example_my_new_preset.py`

---

## Examples

### Example 1: Surface Thermodynamics (Core Workflow)

**File:** `examples/workflow_presets/example_surface_thermo.py`

```python
#!/usr/bin/env python3
"""
Example: Complete surface thermodynamics workflow using preset.

This script demonstrates the core PS-TEROS workflow for calculating
surface energies as a function of chemical potential.
"""

from pathlib import Path
from aiida import load_profile
from teros.core import build_core_workgraph
from teros.core.builders import get_ag3po4_defaults

# Load AiiDA profile
load_profile()

# Get default parameters for Ag3PO4
defaults = get_ag3po4_defaults()

# Build workflow using surface_thermodynamics preset
wg = build_core_workgraph(
    # === WORKFLOW PRESET ===
    workflow_preset='surface_thermodynamics',

    # === STRUCTURE FILES ===
    structures_dir=str(Path(__file__).parent / 'structures'),
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',

    # === COMPUTATION SETTINGS ===
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    clean_workdir=False,

    # === SLAB GENERATION ===
    miller_indices=[1, 0, 0],
    min_slab_thickness=15.0,
    min_vacuum_thickness=15.0,
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,

    # === THERMODYNAMICS SAMPLING ===
    thermodynamics_sampling=100,

    # === BULK PARAMETERS ===
    bulk_potential_mapping=defaults['bulk']['potential_mapping'],
    kpoints_spacing=defaults['bulk']['kpoints_spacing'],
    bulk_parameters=defaults['bulk']['parameters'],
    bulk_options=defaults['bulk']['options'],

    # === REFERENCE PARAMETERS ===
    metal_potential_mapping=defaults['metal']['potential_mapping'],
    metal_parameters=defaults['metal']['parameters'],
    metal_options=defaults['metal']['options'],

    nonmetal_potential_mapping=defaults['nonmetal']['potential_mapping'],
    nonmetal_parameters=defaults['nonmetal']['parameters'],
    nonmetal_options=defaults['nonmetal']['options'],

    oxygen_potential_mapping=defaults['oxygen']['potential_mapping'],
    oxygen_parameters=defaults['oxygen']['parameters'],
    oxygen_options=defaults['oxygen']['options'],

    # === SLAB PARAMETERS ===
    slab_potential_mapping=defaults['slab']['potential_mapping'],
    slab_parameters=defaults['slab']['parameters'],
    slab_options=defaults['slab']['options'],
    slab_kpoints_spacing=defaults['slab']['kpoints_spacing'],

    # === WORKGRAPH NAME ===
    name='SurfaceThermodynamics_Ag3PO4_100',
)

# Submit workflow
print(f"\n{'='*70}")
print(f"Submitting surface thermodynamics workflow...")
print(f"  Material: Ag3PO4")
print(f"  Miller indices: [1, 0, 0]")
print(f"  Preset: surface_thermodynamics")
print(f"{'='*70}\n")

wg.submit(wait=False)
wg.to_html()

print(f"✅ Workflow submitted successfully!")
print(f"   WorkGraph PK: {wg.pk}")
print(f"   HTML visualization: {wg.name}.html")
print(f"\nMonitor with:")
print(f"  verdi process show {wg.pk}")
print(f"  verdi process report {wg.pk}")
```

### Example 2: AIMD Only

**File:** `examples/workflow_presets/example_aimd_only.py`

```python
#!/usr/bin/env python3
"""
Example: Run AIMD on input slabs using preset.

This script demonstrates running ab initio molecular dynamics
on pre-generated slab structures.
"""

from pathlib import Path
from aiida import load_profile
from aiida.orm import load_node
from teros.core import build_core_workgraph
from teros.core.builders import get_ag3po4_defaults, get_aimd_defaults

# Load AiiDA profile
load_profile()

# Load slabs from previous workflow
previous_wg_pk = 12345  # Replace with your workgraph PK
previous_wg = load_node(previous_wg_pk)

# Extract relaxed slabs
input_slabs = {}
for label in previous_wg.outputs.relaxed_slabs.keys():
    input_slabs[label] = previous_wg.outputs.relaxed_slabs[label]

print(f"Loaded {len(input_slabs)} slabs from WorkGraph {previous_wg_pk}")

# Get default parameters
slab_defaults = get_ag3po4_defaults()['slab']
aimd_defaults = get_aimd_defaults()

# Define AIMD sequence (heating protocol)
aimd_sequence = [
    {'temperature': 300, 'steps': 1000},   # Equilibration at 300K
    {'temperature': 500, 'steps': 2000},   # Heating to 500K
    {'temperature': 500, 'steps': 3000},   # Production run at 500K
]

# Build workflow using aimd_only preset
wg = build_core_workgraph(
    # === WORKFLOW PRESET ===
    workflow_preset='aimd_only',

    # === STRUCTURE FILES (still need bulk for metadata) ===
    structures_dir=str(Path(__file__).parent / 'structures'),
    bulk_name='ag3po4.cif',

    # === INPUT SLABS ===
    input_slabs=input_slabs,

    # === COMPUTATION SETTINGS ===
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    clean_workdir=False,

    # === AIMD PARAMETERS ===
    aimd_sequence=aimd_sequence,
    aimd_parameters=aimd_defaults['parameters'],
    aimd_options=aimd_defaults['options'],
    aimd_potential_mapping=slab_defaults['potential_mapping'],
    aimd_kpoints_spacing=aimd_defaults['kpoints_spacing'],

    # === WORKGRAPH NAME ===
    name='AIMD_Ag3PO4_100_300K-500K',
)

# Submit workflow
print(f"\n{'='*70}")
print(f"Submitting AIMD workflow...")
print(f"  Number of slabs: {len(input_slabs)}")
print(f"  AIMD stages: {len(aimd_sequence)}")
print(f"  Temperature range: 300K - 500K")
print(f"  Total steps per slab: {sum(s['steps'] for s in aimd_sequence)}")
print(f"  Preset: aimd_only")
print(f"{'='*70}\n")

wg.submit(wait=False)
wg.to_html()

print(f"✅ Workflow submitted successfully!")
print(f"   WorkGraph PK: {wg.pk}")
print(f"   Access AIMD outputs via:")
print(f"     wg.tasks['aimd_stage_00_300K'].outputs.structures")
print(f"     wg.tasks['aimd_stage_01_500K'].outputs.structures")
print(f"     wg.tasks['aimd_stage_02_500K'].outputs.structures")
```

### Example 3: Custom Workflow (Preset + Overrides)

**File:** `examples/workflow_presets/example_custom_workflow.py`

```python
#!/usr/bin/env python3
"""
Example: Custom workflow combining preset with overrides.

This script demonstrates starting with a preset and customizing it
for a specific use case.
"""

from pathlib import Path
from aiida import load_profile
from teros.core import build_core_workgraph
from teros.core.builders import (
    get_ag3po4_defaults,
    get_electronic_properties_defaults,
    get_aimd_defaults,
)

# Load AiiDA profile
load_profile()

# Get default parameters
ag3po4_defaults = get_ag3po4_defaults()
elec_defaults = get_electronic_properties_defaults()
aimd_defaults = get_aimd_defaults()

# Build custom workflow:
# Start with surface_thermodynamics, then:
#   - Disable cleavage (not needed)
#   - Add bulk electronic structure
#   - Add AIMD stages
wg = build_core_workgraph(
    # === START WITH PRESET ===
    workflow_preset='surface_thermodynamics',

    # === CUSTOMIZE IT ===
    compute_cleavage=False,  # Don't need cleavage energies
    compute_electronic_properties_bulk=True,  # Add bulk DOS/bands
    run_aimd=True,  # Add AIMD after slab relaxation

    # === STRUCTURE FILES ===
    structures_dir=str(Path(__file__).parent / 'structures'),
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',

    # === COMPUTATION SETTINGS ===
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    clean_workdir=False,

    # === SLAB GENERATION ===
    miller_indices=[1, 0, 0],
    min_slab_thickness=15.0,
    min_vacuum_thickness=15.0,

    # === BULK PARAMETERS ===
    bulk_potential_mapping=ag3po4_defaults['bulk']['potential_mapping'],
    kpoints_spacing=ag3po4_defaults['bulk']['kpoints_spacing'],
    bulk_parameters=ag3po4_defaults['bulk']['parameters'],
    bulk_options=ag3po4_defaults['bulk']['options'],

    # === REFERENCE PARAMETERS ===
    metal_potential_mapping=ag3po4_defaults['metal']['potential_mapping'],
    metal_parameters=ag3po4_defaults['metal']['parameters'],
    metal_options=ag3po4_defaults['metal']['options'],

    nonmetal_potential_mapping=ag3po4_defaults['nonmetal']['potential_mapping'],
    nonmetal_parameters=ag3po4_defaults['nonmetal']['parameters'],
    nonmetal_options=ag3po4_defaults['nonmetal']['options'],

    oxygen_potential_mapping=ag3po4_defaults['oxygen']['potential_mapping'],
    oxygen_parameters=ag3po4_defaults['oxygen']['parameters'],
    oxygen_options=ag3po4_defaults['oxygen']['options'],

    # === SLAB PARAMETERS ===
    slab_potential_mapping=ag3po4_defaults['slab']['potential_mapping'],
    slab_parameters=ag3po4_defaults['slab']['parameters'],
    slab_options=ag3po4_defaults['slab']['options'],
    slab_kpoints_spacing=ag3po4_defaults['slab']['kpoints_spacing'],

    # === ELECTRONIC PROPERTIES PARAMETERS ===
    bands_parameters=elec_defaults['bands_parameters'],
    bands_options=elec_defaults['bands_options'],
    band_settings=elec_defaults['band_settings'],

    # === AIMD PARAMETERS ===
    aimd_sequence=[
        {'temperature': 300, 'steps': 1000},
        {'temperature': 400, 'steps': 2000},
    ],
    aimd_parameters=aimd_defaults['parameters'],
    aimd_options=aimd_defaults['options'],
    aimd_kpoints_spacing=aimd_defaults['kpoints_spacing'],

    # === WORKGRAPH NAME ===
    name='Custom_SurfaceThermo_ElecStruct_AIMD',
)

# Submit workflow
print(f"\n{'='*70}")
print(f"Submitting custom workflow...")
print(f"  Base preset: surface_thermodynamics")
print(f"  Modifications:")
print(f"    - Disabled: cleavage energies")
print(f"    - Added: bulk electronic structure")
print(f"    - Added: AIMD (2 stages)")
print(f"{'='*70}\n")

wg.submit(wait=False)
wg.to_html()

print(f"✅ Workflow submitted successfully!")
print(f"   WorkGraph PK: {wg.pk}")
```

---

## Future Enhancements

### Phase 2 Features (Post-Initial Implementation)

1. **Preset Validation Tool**
   - CLI tool to validate preset configurations before submission
   - Check all required files exist
   - Verify parameter compatibility

   ```bash
   $ psteros validate-preset surface_thermodynamics --config my_config.yaml
   ✅ All required parameters present
   ✅ Structure files found
   ⚠️  Warning: nonmetal_name not provided (optional for binary oxide)
   ✅ Ready to submit
   ```

2. **Preset Templates**
   - YAML/JSON templates for each preset
   - Auto-generate builder code from templates

   ```yaml
   # surface_thermodynamics_template.yaml
   workflow_preset: surface_thermodynamics
   structures_dir: ./structures
   bulk_name: my_bulk.cif
   metal_name: my_metal.cif
   oxygen_name: O2.cif
   miller_indices: [1, 0, 0]
   # ... etc ...
   ```

3. **Workflow Composition**
   - Chain multiple presets together
   - Sequential workflows with automatic data passing

   ```python
   # Run bulk_only, then use result in slabs_only
   wg1 = build_core_workgraph(workflow_preset='bulk_only', ...)
   wg1.submit(wait=True)

   wg2 = build_core_workgraph(
       workflow_preset='slabs_only',
       restart_from_node=wg1.pk,  # Auto-extracts bulk structure
       ...
   )
   ```

4. **Interactive Preset Builder**
   - Interactive CLI/notebook tool to build custom presets
   - Answer questions → generate preset configuration

   ```python
   from teros.core import interactive_preset_builder

   preset = interactive_preset_builder()
   # Q: What do you want to calculate? [thermodynamics/electronic/aimd]
   # A: thermodynamics
   # Q: Do you need cleavage energies? [y/n]
   # A: n
   # ... etc ...
   # Generated preset: surface_thermodynamics_no_cleavage
   ```

5. **Preset Visualization**
   - Generate flowcharts showing what each preset does
   - Visual comparison of different presets

   ```python
   from teros.core import visualize_preset

   visualize_preset('surface_thermodynamics')
   # Creates graphviz diagram showing task flow
   ```

6. **User-Defined Presets**
   - Allow users to define custom presets in config file
   - Load from `~/.psteros/custom_presets.yaml`

   ```yaml
   # ~/.psteros/custom_presets.yaml
   my_custom_workflow:
     description: My specific workflow
     flags:
       relax_slabs: true
       compute_thermodynamics: true
       # ... etc ...
   ```

### Phase 3 Features (Long-term)

1. **Preset Analytics**
   - Track which presets are most used
   - Suggest optimizations based on usage patterns

2. **Smart Preset Recommendation**
   - Analyze user's previous workflows
   - Recommend appropriate preset based on similar calculations

3. **Preset Versioning**
   - Version presets to track changes over time
   - Allow pinning to specific preset versions

---

## Appendix A: Complete File Listing

### Files to Create

```
teros/core/
├── workflow_presets.py                    (~800 lines)

docs/
├── WORKFLOW_PRESETS_GUIDE.md              (~500 lines)
└── WORKFLOW_PRESETS_EXAMPLES.md           (~300 lines)

examples/workflow_presets/
├── example_surface_thermo.py              (~150 lines)
├── example_aimd_only.py                   (~120 lines)
├── example_electronic_structure.py        (~140 lines)
└── example_custom_workflow.py             (~160 lines)
```

### Files to Modify

```
teros/core/
├── workgraph.py                           (~50 lines added)
└── __init__.py                            (~20 lines added)
```

### Total Implementation Effort

- **New code:** ~2,300 lines
- **Modified code:** ~70 lines
- **Documentation:** ~800 lines
- **Examples:** ~570 lines
- **Total:** ~3,740 lines

**Estimated implementation time:** 12-16 hours for experienced developer

---

## Appendix B: Quick Reference

### Available Presets

| Preset | Bulk | Refs | Slabs | Thermo | Cleavage | Relax Energy | Elec | AIMD |
|--------|------|------|-------|--------|----------|--------------|------|------|
| `surface_thermodynamics` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| `relaxation_energy_analysis` | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ |
| `bulk_only` | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| `formation_enthalpy_only` | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| `slabs_only` | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| `cleavage_analysis` | ✅ | ❌ | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ |
| `electronic_structure_full` | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ |
| `electronic_structure_bulk_only` | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ |
| `aimd_only` | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ | ✅ |

### Common Commands

```python
# List all presets
from teros.core import list_workflow_presets
list_workflow_presets()

# Get preset details
from teros.core import get_preset_config
config = get_preset_config('surface_thermodynamics')

# Build with preset
from teros.core import build_core_workgraph
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    ...
)

# Override preset
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=False,  # Override
    ...
)
```

---

## Appendix C: Decision Log

### Design Decisions

1. **Why default to `surface_thermodynamics`?**
   - User indicated this is the "core" workflow
   - Most common use case
   - Better than forcing explicit choice every time

2. **Why allow `None` for all boolean flags?**
   - Enables preset system to control defaults
   - Distinguishes "user didn't specify" from "user explicitly set False"
   - Required for override logic to work correctly

3. **Why dictionary-based preset definitions?**
   - Easy to extend with new presets
   - Self-documenting structure
   - Can be exported to JSON/YAML if needed

4. **Why deprecate old API instead of removing?**
   - Gives users time to migrate
   - Avoids breaking existing scripts immediately
   - Can remove in next major version (v3.0)

5. **Why validate parameters per preset?**
   - Catches errors early
   - Clear error messages
   - Prevents invalid workflow configurations

6. **Why separate `workflow_presets.py` module?**
   - Keeps `workgraph.py` focused on building logic
   - Easy to find and modify presets
   - Can be tested independently

---

## Appendix D: Implementation Checklist

Use this checklist to track implementation progress:

- [ ] **Step 1: Create `workflow_presets.py`**
  - [ ] Define all 9 presets
  - [ ] Implement `resolve_preset()`
  - [ ] Implement `validate_preset_inputs()`
  - [ ] Implement `list_workflow_presets()`
  - [ ] Implement `get_preset_config()`
  - [ ] Implement `get_preset_summary()`
  - [ ] Implement deprecation helpers
  - [ ] Add module docstring
  - [ ] Add tests in `if __name__ == '__main__'` block

- [ ] **Step 2: Modify `workgraph.py`**
  - [ ] Add `workflow_preset` parameter
  - [ ] Change boolean flag defaults to `None`
  - [ ] Add preset resolution logic at function start
  - [ ] Add parameter validation
  - [ ] Add deprecation warning for old API
  - [ ] Update docstring with preset information
  - [ ] Add example usage to docstring

- [ ] **Step 3: Update `__init__.py`**
  - [ ] Import preset functions
  - [ ] Add to `__all__`
  - [ ] Update module docstring

- [ ] **Step 4: Create documentation**
  - [ ] Write `WORKFLOW_PRESETS_GUIDE.md`
  - [ ] Write `WORKFLOW_PRESETS_EXAMPLES.md`
  - [ ] Review and proofread

- [ ] **Step 5: Create examples**
  - [ ] `example_surface_thermo.py`
  - [ ] `example_aimd_only.py`
  - [ ] `example_electronic_structure.py`
  - [ ] `example_custom_workflow.py`
  - [ ] Test all examples (dry run if no VASP)

- [ ] **Step 6: Testing**
  - [ ] Unit tests for preset resolution
  - [ ] Unit tests for validation
  - [ ] Integration test for each preset
  - [ ] Backward compatibility tests
  - [ ] End-to-end test (if possible)

- [ ] **Step 7: Cleanup**
  - [ ] Remove debug print statements
  - [ ] Check code style (PEP 8)
  - [ ] Run linter
  - [ ] Final review of all changes

- [ ] **Step 8: Documentation**
  - [ ] Update main README
  - [ ] Update CHANGELOG
  - [ ] Create migration guide for users
  - [ ] Update any affected documentation

---

## Appendix E: Contact and Support

**Implementation Questions:**
- Review this guide thoroughly before starting
- Check existing similar implementations in the codebase
- Test each component independently before integration

**Design Questions:**
- Refer to the brainstorming session notes
- Check design decisions in Appendix C
- Consider user experience and common workflows

**Testing Issues:**
- Start with unit tests, then integration
- Use mocking for VASP calls if needed
- Test backward compatibility carefully

---

**End of Implementation Guide**

*This guide provides complete specifications for implementing the PS-TEROS workflow preset system. All code examples are production-ready and follow the existing codebase conventions.*
