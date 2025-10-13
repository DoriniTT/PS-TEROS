#!/usr/bin/env python
"""
Example: AIMD Simulation using Preset

This example demonstrates Ab Initio Molecular Dynamics simulation
on surface slabs using the 'aimd_only' preset.

What this calculates:
- Bulk relaxation
- Slab generation and relaxation
- AIMD simulation on all slabs
"""

from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_aimd_defaults
from aiida.engine import submit

# ============================================================================
# CONFIGURATION
# ============================================================================

STRUCTURES_DIR = '/path/to/your/structures'  # CHANGE THIS
BULK_NAME = 'ag3po4.cif'

CODE_LABEL = 'VASP-VTST-6.4.3@bohr'  # CHANGE THIS
POTENTIAL_FAMILY = 'PBE'

# Standard VASP parameters for relaxation
RELAX_PARAMS = {
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
        'num_machines': 2,  # AIMD may need more resources
        'num_mpiprocs_per_machine': 48,
    },
    'max_wallclock_seconds': 3600 * 48,  # 48 hours for AIMD
}

# AIMD sequence: equilibration + production
AIMD_SEQUENCE = [
    {'temperature': 300, 'steps': 1000},  # 1 ps equilibration
    {'temperature': 300, 'steps': 5000},  # 5 ps production
]

# ============================================================================
# BUILD AND SUBMIT WORKFLOW
# ============================================================================

def main():
    """Build and submit the workflow."""
    
    print("\n" + "="*70)
    print("AIMD SIMULATION WORKFLOW - PRESET EXAMPLE")
    print("="*70)
    print(f"Material: Ag3PO4")
    print(f"Surface: (100)")
    print(f"Preset: aimd_only")
    print(f"AIMD: {len(AIMD_SEQUENCE)} stages")
    for i, stage in enumerate(AIMD_SEQUENCE, 1):
        print(f"  Stage {i}: {stage['temperature']} K, {stage['steps']} steps")
    print("="*70 + "\n")
    
    # Get default AIMD parameters
    aimd_params = get_aimd_defaults()
    
    wg = build_core_workgraph(
        # ====== PRESET ======
        workflow_preset='aimd_only',
        
        # ====== Structure ======
        structures_dir=STRUCTURES_DIR,
        bulk_name=BULK_NAME,
        
        # ====== Code ======
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.4,
        
        # ====== Bulk ======
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        bulk_parameters=RELAX_PARAMS,
        bulk_options=SCHEDULER_OPTIONS,
        
        # ====== Slabs ======
        miller_indices=[1, 0, 0],
        min_slab_thickness=15.0,  # Smaller for AIMD
        min_vacuum_thickness=15.0,
        slab_parameters=RELAX_PARAMS.copy(),
        slab_options=SCHEDULER_OPTIONS,
        
        # ====== AIMD Configuration ======
        aimd_sequence=AIMD_SEQUENCE,
        aimd_parameters=aimd_params,
        aimd_options=SCHEDULER_OPTIONS,
        aimd_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        aimd_kpoints_spacing=0.5,  # Coarser k-points for AIMD
        
        # ====== Name ======
        name='Ag3PO4_100_AIMD',
    )
    
    print("Submitting AIMD workflow...")
    result = submit(wg)
    
    print(f"\n{'='*70}")
    print(f"✓ AIMD workflow submitted!")
    print(f"{'='*70}")
    print(f"WorkGraph PK: {result.pk}")
    print(f"\nThis will take a while (AIMD is expensive)...")
    print(f"Monitor: verdi process status {result.pk}")
    print(f"{'='*70}\n")
    
    return result


if __name__ == '__main__':
    main()
