#!/usr/bin/env python
"""
Example: Custom Workflow by Overriding Preset

This example demonstrates how to start with a preset and override
specific components to create a custom workflow.

Starting preset: surface_thermodynamics
Overrides:
- Disable cleavage energy calculations
- Disable relaxation energy calculations

Result: Surface thermodynamics without cleavage/relaxation analysis
"""

from teros.core.workgraph import build_core_workgraph
from teros.core import list_workflow_presets, get_preset_summary
from aiida.engine import submit

# ============================================================================
# CONFIGURATION
# ============================================================================

STRUCTURES_DIR = '/path/to/your/structures'  # CHANGE THIS
BULK_NAME = 'ag3po4.cif'
METAL_NAME = 'Ag.cif'
OXYGEN_NAME = 'O2.cif'
NONMETAL_NAME = 'P.cif'

CODE_LABEL = 'VASP-VTST-6.4.3@bohr'  # CHANGE THIS
POTENTIAL_FAMILY = 'PBE'

VASP_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'ALGO': 'Normal',
    'IBRION': 2,
    'NSW': 100,
    'ISIF': 3,
}

SCHEDULER_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 48,
    },
    'max_wallclock_seconds': 3600 * 10,
}

# ============================================================================
# BUILD AND SUBMIT WORKFLOW
# ============================================================================

def main():
    """Build and submit the workflow."""
    
    print("\n" + "="*70)
    print("CUSTOM WORKFLOW - OVERRIDING PRESET EXAMPLE")
    print("="*70)
    print(f"Material: Ag3PO4")
    print(f"Surface: (111)")
    print(f"Base preset: surface_thermodynamics")
    print(f"Overrides:")
    print(f"  - compute_cleavage = False")
    print(f"  - compute_relaxation_energy = False")
    print(f"\nResult: Surface thermodynamics only (no cleavage/relaxation)")
    print("="*70 + "\n")
    
    # Optional: Show what the base preset does
    print("Base preset configuration:")
    print(get_preset_summary('surface_thermodynamics'))
    
    input("Press Enter to submit workflow with overrides...")
    
    wg = build_core_workgraph(
        # ====== PRESET (base configuration) ======
        workflow_preset='surface_thermodynamics',
        
        # ====== OVERRIDES (custom modifications) ======
        compute_cleavage=False,  # Disable cleavage energies
        compute_relaxation_energy=False,  # Disable relaxation energies
        
        # ====== Structure files ======
        structures_dir=STRUCTURES_DIR,
        bulk_name=BULK_NAME,
        metal_name=METAL_NAME,
        oxygen_name=OXYGEN_NAME,
        nonmetal_name=NONMETAL_NAME,
        
        # ====== Code ======
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.4,
        
        # ====== Bulk ======
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        bulk_parameters=VASP_PARAMS,
        bulk_options=SCHEDULER_OPTIONS,
        
        # ====== References ======
        metal_potential_mapping={'Ag': 'Ag'},
        metal_parameters=VASP_PARAMS.copy(),
        metal_options=SCHEDULER_OPTIONS,
        
        oxygen_potential_mapping={'O': 'O'},
        oxygen_parameters=VASP_PARAMS.copy(),
        oxygen_options=SCHEDULER_OPTIONS,
        
        nonmetal_potential_mapping={'P': 'P'},
        nonmetal_parameters=VASP_PARAMS.copy(),
        nonmetal_options=SCHEDULER_OPTIONS,
        
        # ====== Slabs ======
        miller_indices=[1, 1, 1],  # (111) surface
        min_slab_thickness=18.0,
        min_vacuum_thickness=15.0,
        slab_parameters=VASP_PARAMS.copy(),
        slab_options=SCHEDULER_OPTIONS,
        
        # ====== Thermodynamics ======
        thermodynamics_sampling=100,
        
        # ====== Name ======
        name='Ag3PO4_111_Custom',
    )
    
    print("\nSubmitting custom workflow...")
    result = submit(wg)
    
    print(f"\n{'='*70}")
    print(f"✓ Custom workflow submitted!")
    print(f"{'='*70}")
    print(f"WorkGraph PK: {result.pk}")
    print(f"\nActive components:")
    print(f"  ✓ Bulk relaxation")
    print(f"  ✓ Reference relaxations")
    print(f"  ✓ Formation enthalpy")
    print(f"  ✓ Slab generation and relaxation")
    print(f"  ✓ Surface thermodynamics")
    print(f"  ✗ Cleavage energies (disabled)")
    print(f"  ✗ Relaxation energies (disabled)")
    print(f"\nMonitor: verdi process status {result.pk}")
    print(f"{'='*70}\n")
    
    return result


def explore_presets():
    """Interactive exploration of presets."""
    print("\n" + "="*70)
    print("EXPLORING WORKFLOW PRESETS")
    print("="*70 + "\n")
    
    # List all presets
    list_workflow_presets()
    
    # Show details of specific presets
    presets_to_show = ['surface_thermodynamics', 'cleavage_only', 'aimd_only']
    
    print("\nDetailed configurations:")
    for preset in presets_to_show:
        print(get_preset_summary(preset))
        input("Press Enter for next preset...")


if __name__ == '__main__':
    # Uncomment to explore presets first:
    # explore_presets()
    
    main()
