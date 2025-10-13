#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 5: Surface Thermodynamics (Complete Workflow)

This script tests the COMPLETE surface thermodynamics workflow:
- Bulk and reference relaxations
- Formation enthalpy
- Slab generation and relaxation
- Surface energies with chemical potential sampling
- Cleavage energies
- Relaxation energies

This is the default workflow preset 'surface_thermodynamics'.

Material: Ag2O
Surface: (111)
References: Ag, O2

Usage:
    source ~/envs/psteros/bin/activate
    python step_05_surface_thermodynamics.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Step 5: Test complete surface thermodynamics workflow."""
    
    print("\n" + "="*70)
    print("STEP 5: SURFACE THERMODYNAMICS (COMPLETE WORKFLOW)")
    print("="*70)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, '../complete/structures')
    
    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")
    
    # Code configuration
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'
    
    # Common VASP parameters
    vasp_params = {
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
    slab_params = vasp_params.copy()
    slab_params['ISIF'] = 2
    
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }
    
    print("\n3. Building workgraph...")
    print("   Using preset: 'surface_thermodynamics' (default)")
    print("   This calculates:")
    print("     - Formation enthalpy")
    print("     - Surface energies: γ(μ_O)")
    print("     - Cleavage energies")
    print("     - Relaxation energies")
    
    # Build workgraph using default preset
    wg = build_core_workgraph(
        workflow_preset='surface_thermodynamics',
        
        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        metal_name='Ag.cif',
        oxygen_name='O2.cif',
        
        # Code
        code_label=code_label,
        potential_family=potential_family,
        kpoints_spacing=0.4,
        clean_workdir=False,
        
        # Bulk
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=vasp_params.copy(),
        bulk_options=common_options,
        
        # Metal
        metal_potential_mapping={'Ag': 'Ag'},
        metal_parameters=vasp_params.copy(),
        metal_options=common_options,
        
        # Oxygen
        oxygen_potential_mapping={'O': 'O'},
        oxygen_parameters=vasp_params.copy(),
        oxygen_options=common_options,
        
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
        
        # Thermodynamics sampling
        thermodynamics_sampling=100,
        
        name='Step05_SurfaceThermodynamics_Ag2O_111',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n4. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 5 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  1. Energies:")
    print(f"     - bulk_energy, metal_energy, oxygen_energy")
    print(f"  2. Formation enthalpy:")
    print(f"     - formation_enthalpy")
    print(f"  3. Slab structures:")
    print(f"     - slab_structures (all terminations)")
    print(f"  4. Slab energies:")
    print(f"     - unrelaxed_slab_energies")
    print(f"     - slab_energies (relaxed)")
    print(f"     - relaxation_energies")
    print(f"  5. Surface properties:")
    print(f"     - cleavage_energies")
    print(f"     - surface_energies: γ(μ_O)")
    print(f"\nSurface energies show stability as function of chemical potential")
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
