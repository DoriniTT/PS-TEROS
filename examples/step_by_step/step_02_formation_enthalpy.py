#!/home/thiagotd/envs/psteros/bin/python
"""
STEP 2: Formation Enthalpy Calculation

This script tests:
- Bulk relaxation (Ag2O)
- Reference relaxations (Ag, O2)
- Formation enthalpy calculation

This is step 2 in understanding the PS-TEROS workflow chain.

Material: Ag2O
References: Ag (metal), O2 (oxygen)

Usage:
    source ~/envs/psteros/bin/activate
    python step_02_formation_enthalpy.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Step 2: Test formation enthalpy calculation."""
    
    print("\n" + "="*70)
    print("STEP 2: FORMATION ENTHALPY CALCULATION")
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
    
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'max_wallclock_seconds': 3600 * 10,
        'queue_name': 'qchem',
    }
    
    print("\n3. Building workgraph...")
    print("   Using preset: 'formation_enthalpy_only'")
    
    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='formation_enthalpy_only',
        
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
        
        # Bulk (Ag2O)
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=vasp_params.copy(),
        bulk_options=common_options,
        
        # Metal (Ag)
        metal_potential_mapping={'Ag': 'Ag'},
        metal_parameters=vasp_params.copy(),
        metal_options=common_options,
        
        # Oxygen (O2)
        oxygen_potential_mapping={'O': 'O'},
        oxygen_parameters=vasp_params.copy(),
        oxygen_options=common_options,
        
        name='Step02_FormationEnthalpy_Ag2O',
    )
    
    print("   ✓ WorkGraph built successfully")
    
    # Submit
    print("\n4. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("STEP 2 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(Ag2O)")
    print(f"  - metal_energy: E(Ag)")
    print(f"  - oxygen_energy: E(O2)")
    print(f"  - formation_enthalpy: ΔHf = E(Ag2O) - 2*E(Ag) - 0.5*E(O2)")
    print(f"\nFormation enthalpy is reported in eV/formula unit")
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
