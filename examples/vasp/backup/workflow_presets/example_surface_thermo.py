#!/usr/bin/env python
"""
Example: Surface Thermodynamics Workflow using Preset

This example demonstrates the simplest way to run a complete surface
thermodynamics workflow using the 'surface_thermodynamics' preset.

What this calculates:
- Bulk relaxation
- Reference relaxations (Ag, P, O2)
- Formation enthalpy
- Slab generation and relaxation for (100) surface
- Surface energies with chemical potential sampling
- Cleavage energies
- Relaxation energies
"""

from teros.core.workgraph import build_core_workgraph
from aiida.engine import submit

# ============================================================================
# CONFIGURATION
# ============================================================================

# Structure files
STRUCTURES_DIR = '/path/to/your/structures'  # CHANGE THIS
BULK_NAME = 'ag3po4.cif'
METAL_NAME = 'Ag.cif'
OXYGEN_NAME = 'O2.cif'
NONMETAL_NAME = 'P.cif'

# AiiDA configuration
CODE_LABEL = 'VASP-VTST-6.4.3@bohr'  # CHANGE THIS to your code label
POTENTIAL_FAMILY = 'PBE'

# VASP parameters for bulk and references
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
    'LREAL': False,
    'LWAVE': False,
    'LCHARG': False,
}

# Scheduler options
SCHEDULER_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 48,
    },
    'max_wallclock_seconds': 3600 * 10,  # 10 hours
    'queue_name': 'regular',
}

# Slab generation parameters
MILLER_INDICES = [1, 0, 0]  # (100) surface
MIN_SLAB_THICKNESS = 18.0  # Angstroms
MIN_VACUUM_THICKNESS = 15.0  # Angstroms

# ============================================================================
# BUILD AND SUBMIT WORKFLOW
# ============================================================================

def main():
    """Build and submit the workflow."""
    
    print("\n" + "="*70)
    print("SURFACE THERMODYNAMICS WORKFLOW - PRESET EXAMPLE")
    print("="*70)
    print(f"Material: Ag3PO4")
    print(f"Surface: {MILLER_INDICES}")
    print(f"Preset: surface_thermodynamics")
    print("="*70 + "\n")
    
    # Build workflow using preset
    wg = build_core_workgraph(
        # ====== PRESET (this line activates full workflow!) ======
        workflow_preset='surface_thermodynamics',
        
        # ====== Structure files ======
        structures_dir=STRUCTURES_DIR,
        bulk_name=BULK_NAME,
        metal_name=METAL_NAME,
        oxygen_name=OXYGEN_NAME,
        nonmetal_name=NONMETAL_NAME,
        
        # ====== Code configuration ======
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.4,
        clean_workdir=False,
        
        # ====== Bulk calculation ======
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        bulk_parameters=VASP_PARAMS,
        bulk_options=SCHEDULER_OPTIONS,
        
        # ====== Metal reference (Ag) ======
        metal_potential_mapping={'Ag': 'Ag'},
        metal_parameters=VASP_PARAMS.copy(),
        metal_options=SCHEDULER_OPTIONS,
        
        # ====== Oxygen reference (O2) ======
        oxygen_potential_mapping={'O': 'O'},
        oxygen_parameters=VASP_PARAMS.copy(),
        oxygen_options=SCHEDULER_OPTIONS,
        
        # ====== Nonmetal reference (P) ======
        nonmetal_potential_mapping={'P': 'P'},
        nonmetal_parameters=VASP_PARAMS.copy(),
        nonmetal_options=SCHEDULER_OPTIONS,
        
        # ====== Slab configuration ======
        miller_indices=MILLER_INDICES,
        min_slab_thickness=MIN_SLAB_THICKNESS,
        min_vacuum_thickness=MIN_VACUUM_THICKNESS,
        slab_parameters=VASP_PARAMS.copy(),
        slab_options=SCHEDULER_OPTIONS,
        slab_kpoints_spacing=0.4,
        
        # ====== Slab generation options ======
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        
        # ====== Thermodynamics sampling ======
        thermodynamics_sampling=100,
        
        # ====== Workflow name ======
        name='Ag3PO4_100_SurfaceThermodynamics',
    )
    
    # Submit workflow
    print("Submitting workflow to AiiDA...")
    result = submit(wg)
    
    print(f"\n{'='*70}")
    print(f"✓ Workflow submitted successfully!")
    print(f"{'='*70}")
    print(f"WorkGraph PK: {result.pk}")
    print(f"\nMonitor progress with:")
    print(f"  verdi process status {result.pk}")
    print(f"  verdi process show {result.pk}")
    print(f"\nAfter completion, analyze results with:")
    print(f"  verdi process show {result.pk}")
    print(f"{'='*70}\n")
    
    return result


if __name__ == '__main__':
    main()
