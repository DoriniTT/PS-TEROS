#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 8: Electronic Structure for Slabs Only

This script tests:
- Bulk relaxation
- Slab generation and relaxation
- Electronic structure calculation (DOS and bands) for slabs only

This demonstrates the electronic properties workflow for surface structures.

Material: Ag2O
Surface: (111)

Usage:
    source ~/envs/psteros/bin/activate
    python step_08_electronic_structure_slabs.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_slab_electronic_properties_defaults

def main():
    """Step 8: Test electronic properties calculation for slabs only."""
    
    print("\n" + "="*70)
    print("STEP 8: ELECTRONIC STRUCTURE FOR SLABS ONLY")
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
    
    # Bulk parameters
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
    
    # Slab parameters
    slab_params = bulk_params.copy()
    slab_params['ISIF'] = 2
    
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }
    
    # Get slab electronic properties defaults
    print("\n3. Setting up slab electronic properties parameters...")
    slab_elec_defaults = get_slab_electronic_properties_defaults()
    
    print("   Band structure will be calculated along high-symmetry paths for slabs")
    print("   DOS will include total and projected components")
    
    print("\n4. Building workgraph...")
    print("   Using preset: 'electronic_structure_slabs_only'")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='electronic_structure_slabs_only',
        
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
        
        # Slab electronic properties
        slab_bands_parameters=slab_elec_defaults,
        slab_bands_options=common_options,
        slab_band_settings=slab_elec_defaults['band_settings'],
        # Concurrency control (limits simultaneous VASP calculations)
        max_concurrent_jobs=4,  # Default: 4 concurrent calculations

        
        name='Step08_ElectronicStructure_Slabs_Ag2O_111',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 8 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures: All terminations")
    print(f"  - slab_energies: Relaxed slab energies")
    print(f"  - slab_bands: Band structure for each slab")
    print(f"  - slab_dos: Density of states for each slab")
    print(f"  - slab_primitive_structures: Primitive cells used")
    print(f"  - slab_seekpath_parameters: Symmetry info")
    print(f"\nUse AiiDA tools to visualize:")
    print(f"  verdi data bands show <PK>")
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
