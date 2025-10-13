#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 7: AIMD Simulation

This script tests:
- Bulk relaxation
- Slab generation and relaxation
- AIMD (Ab Initio Molecular Dynamics) simulation on slabs

This demonstrates dynamic properties at finite temperature.

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence (equilibration + production)

Usage:
    source ~/envs/psteros/bin/activate
    python step_07_aimd_simulation.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_aimd_defaults

def main():
    """Step 7: Test AIMD simulation."""
    
    print("\n" + "="*70)
    print("STEP 7: AIMD SIMULATION")
    print("="*70)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, '../complete/structures')
    
    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")
    
    # Code configuration
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'
    
    # VASP parameters
    bulk_params = {
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
    
    slab_params = bulk_params.copy()
    slab_params['ISIF'] = 2
    
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }
    
    # AIMD configuration
    print("\n3. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50}, 
        {'temperature': 300, 'steps': 100}, 
    ]
    
    for i, stage in enumerate(aimd_sequence):
        print(f"   Stage {i+1}: {stage['temperature']} K, {stage['steps']} steps")
    
    aimd_params = get_aimd_defaults()
    
    print("\n4. Building workgraph...")
    print("   Using preset: 'aimd_only'")
    print("   AIMD will run on ALL generated slabs")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        
        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        
        # Code
        code_label=code_label,
        potential_family=potential_family,
        kpoints_spacing=0.4,
        clean_workdir=False,
        
        # Bulk
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=bulk_params,
        bulk_options=common_options,
        
        # Slab generation
        miller_indices=[1, 1, 1],
        min_slab_thickness=15.0,  # Smaller for AIMD
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        
        # Slab relaxation
        slab_parameters=slab_params,
        slab_options=common_options,
        slab_kpoints_spacing=0.4,
        
        # AIMD
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=common_options,
        aimd_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        aimd_kpoints_spacing=0.5,  # Coarser for AIMD
        
        name='Step07_AIMD_Ag2O_111',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 7 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures")
    print(f"  - slab_energies (from relaxation)")
    print(f"  - aimd_results: Nested namespace with AIMD data")
    print(f"    - For each slab and stage:")
    print(f"      - trajectory, structure, energy, remote, retrieved")
    print(f"\nAIMD trajectories can be analyzed for:")
    print(f"  - Temperature evolution")
    print(f"  - Atomic diffusion")
    print(f"  - Structural stability")
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
