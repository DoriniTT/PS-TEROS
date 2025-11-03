#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 6: Electronic Properties (DOS and Band Structure)

This script tests:
- Bulk relaxation
- Electronic structure calculation (DOS and bands)

This demonstrates the electronic properties workflow for bulk materials.

Material: Ag2O

Usage:
    source ~/envs/psteros/bin/activate
    python step_06_electronic_properties.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_electronic_properties_defaults

def main():
    """Step 6: Test electronic properties calculation."""
    
    print("\n" + "="*70)
    print("STEP 6: ELECTRONIC PROPERTIES (DOS AND BANDS)")
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
    
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }
    
    # Get electronic properties defaults
    print("\n3. Setting up electronic properties parameters...")
    elec_defaults = get_electronic_properties_defaults()
    
    print("   Band structure will be calculated along high-symmetry paths")
    print("   DOS will include total and projected components")
    
    print("\n4. Building workgraph...")
    print("   Using preset: 'electronic_structure_bulk_only'")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='electronic_structure_bulk_only',
        
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
        
        # Electronic properties
        bands_parameters=elec_defaults,
        bands_options=common_options,
        band_settings=elec_defaults['band_settings'],
        # Concurrency control (limits simultaneous VASP calculations)
        max_concurrent_jobs=4,  # Default: 4 concurrent calculations

        
        name='Step06_ElectronicProperties_Ag2O',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 6 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(bulk)")
    print(f"  - bulk_structure: Relaxed structure")
    print(f"  - bulk_bands: Band structure data")
    print(f"  - bulk_dos: Density of states")
    print(f"  - bulk_primitive_structure: Primitive cell used")
    print(f"  - bulk_seekpath_parameters: Symmetry info")
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
