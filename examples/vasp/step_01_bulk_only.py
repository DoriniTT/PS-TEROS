#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 1: Bulk Relaxation Only

This script tests the most basic functionality of PS-TEROS:
relaxing a bulk structure using VASP through AiiDA.

What this tests:
- Loading AiiDA profile
- Setting up VASP code and potential
- Bulk structure relaxation
- Basic WorkGraph submission

Material: Ag2O (cuprite structure)

Usage:
    source ~/envs/psteros/bin/activate
    python step_01_bulk_only.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Step 1: Test bulk relaxation only."""
    
    print("\n" + "="*70)
    print("STEP 1: BULK RELAXATION ONLY")
    print("="*70)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    
    print(f"\n2. Structure location:")
    print(f"   {structures_dir}/ag2o.cif")
    
    # Define parameters
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'
    
    bulk_parameters = {
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
    
    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }
    
    print("\n3. Building workgraph...")
    print("   Using preset: 'bulk_only'")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='bulk_only',

        structures_dir=structures_dir,
        bulk_name='ag2o.cif',

        code_label=code_label,
        potential_family=potential_family,
        kpoints_spacing=0.4,
        clean_workdir=False,

        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,

        # Concurrency control (limits simultaneous VASP calculations)
        max_concurrent_jobs=4,  # Default: 4 concurrent calculations

        name='Step01_BulkOnly_Ag2O',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n4. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 1 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: Total energy of relaxed Ag2O")
    print(f"  - bulk_structure: Relaxed structure")
    print(f"{'='*70}\n")
    
    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
