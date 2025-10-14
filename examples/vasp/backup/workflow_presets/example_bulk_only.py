#!/usr/bin/env python
"""
Example: Bulk-Only Workflow using Preset

This example demonstrates a simple bulk structure optimization
using the 'bulk_only' preset. No surface calculations.

What this calculates:
- Bulk structure relaxation only
"""

from teros.core.workgraph import build_core_workgraph
from aiida.engine import submit

# ============================================================================
# CONFIGURATION
# ============================================================================

STRUCTURES_DIR = '/path/to/your/structures'  # CHANGE THIS
BULK_NAME = 'ag3po4.cif'

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
    'max_wallclock_seconds': 3600 * 5,
}

# ============================================================================
# BUILD AND SUBMIT WORKFLOW
# ============================================================================

def main():
    """Build and submit the workflow."""
    
    print("\n" + "="*70)
    print("BULK-ONLY WORKFLOW - PRESET EXAMPLE")
    print("="*70)
    print(f"Material: Ag3PO4")
    print(f"Preset: bulk_only")
    print("="*70 + "\n")
    
    wg = build_core_workgraph(
        # ====== PRESET ======
        workflow_preset='bulk_only',
        
        # ====== Structure ======
        structures_dir=STRUCTURES_DIR,
        bulk_name=BULK_NAME,
        
        # ====== Code ======
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.4,
        
        # ====== Bulk calculation ======
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        bulk_parameters=VASP_PARAMS,
        bulk_options=SCHEDULER_OPTIONS,
        
        # ====== Name ======
        name='Ag3PO4_Bulk',
    )
    
    print("Submitting workflow...")
    result = submit(wg)
    
    print(f"\n✓ Workflow submitted: PK = {result.pk}")
    print(f"Monitor: verdi process status {result.pk}\n")
    
    return result


if __name__ == '__main__':
    main()
