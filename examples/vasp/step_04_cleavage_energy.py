#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 4: Cleavage Energy Calculation

This script tests:
- Bulk relaxation
- Slab generation
- Slab relaxation
- Cleavage energy calculation for complementary slab pairs

Cleavage energy = Energy to cleave the bulk into two complementary surfaces.

Material: Ag2O
Surface: (111)

Usage:
    source ~/envs/psteros/bin/activate
    python step_04_cleavage_energy.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Step 4: Test cleavage energy calculation."""
    
    print("\n" + "="*70)
    print("STEP 4: CLEAVAGE ENERGY CALCULATION")
    print("="*70)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    
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
        'EDIFFG': -0.05,
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
    
    print("\n3. Building workgraph...")
    print("   Using preset: 'cleavage_only'")
    print("   Cleavage energy:")
    print("     E_cleave = E_slab1 + E_slab2 - E_bulk")
    print("     where slab1 and slab2 are complementary terminations")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='cleavage_only',
        
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
        miller_indices=[1, 0, 0],
        min_slab_thickness=18.0,
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        
        # Slab relaxation
        slab_parameters=slab_params,
        slab_options=common_options,
        slab_kpoints_spacing=0.4,
        
        name='Step04_CleavageEnergy_Ag2O_111',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n4. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 4 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(bulk)")
    print(f"  - slab_structures: Generated slabs")
    print(f"  - slab_energies: Relaxed slab energies")
    print(f"  - cleavage_energies: E_cleave for each pair")
    print(f"\nCleavage energies are in eV/Å²")
    print(f"Lower values indicate easier cleavage")
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
