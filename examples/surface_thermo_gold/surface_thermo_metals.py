#!/home/trevizam/envs/aiida/bin/python
"""
Metal Surface Energy Calculation for FCC Gold (Au)

This script demonstrates the surface_energy module for computing
surface energies of elemental metals using the simple formula:
γ = (E_slab - N·E_bulk/atom) / (2A)

Material: FCC Au
Surfaces: (111), (100), (110) - all in a single WorkGraph

Usage:
    source ~/envs/aiida/bin/activate
    python surface_thermo_metals.py
"""

import sys
import os
from aiida import load_profile
from teros.core.surface_energy import build_metal_surface_energy_workgraph


def main():
    """Run metal surface energy calculation for FCC Au."""
    
    print("\n" + "="*70)
    print("METAL SURFACE ENERGY CALCULATION")
    print("Material: FCC Au")
    print("="*70)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bulk_structure_path = os.path.join(script_dir, 'au.cif')
    
    print(f"\n2. Structure:")
    print(f"   Bulk:   {bulk_structure_path}")
    
    # Code configuration
    code_label = 'VASP-6.5.1-idefix-4-12@obelix'
    potential_family = 'PBE'
    
    # VASP parameters for metals
    # Note: ISMEAR=1 (Methfessel-Paxton) is recommended for metals
    bulk_parameters = {
        'prec': 'Accurate',
        'encut': 500,
        'ediff': 1e-6,
        'ismear': 1,      # Methfessel-Paxton for metals
        'sigma': 0.2,     # Smearing width for metals
        'ibrion': 2,
        'isif': 3,        # Full relaxation for bulk
        'nsw': 100,
        'ediffg': -0.01,
        'algo': 'Normal',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }
    
    # Slab parameters - only relax ionic positions
    slab_parameters = bulk_parameters.copy()
    slab_parameters['isif'] = 2  # Fix cell, relax ions
    
    # Scheduler options for IdeFix
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 4,  # PROCESS_MPI=4 (hybrid MPI+OpenMP)
        },
        'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=idefix-4-12:ppn=88:skylake
#PBS -j oe
#PBS -N Au_surf''',
    }
    
    # Surface orientations to study - all in ONE WorkGraph!
    miller_indices = [
        [1, 1, 1],  # Most stable for FCC
        [1, 0, 0],  # Second most stable
        [1, 1, 0],  # Third
    ]
    
    print(f"\n3. Miller indices: {miller_indices}")
    print("   All orientations in a single WorkGraph!")
    
    print("\n4. Building workgraph...")
    
    # Build ONE workflow containing all Miller indices
    wg = build_metal_surface_energy_workgraph(
        # Structure
        bulk_structure_path=bulk_structure_path,
        
        # Code
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping={'Au': 'Au'},
        kpoints_spacing=0.02,
        clean_workdir=False,
        
        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,
        
        # Slab generation - ALL Miller indices!
        miller_indices=miller_indices,
        min_slab_thickness=20.0,
        min_vacuum_thickness=20.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        
        # Slab relaxation
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.03,
        
        # Concurrency control
        max_concurrent_jobs=4,
        
        name='Au_surface_energy',
    )
    
    print("   ✓ WorkGraph built successfully")
    print(f"   Contains tasks: surface_hkl_111, surface_hkl_100, surface_hkl_110")
    
    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    
    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs per orientation:")
    print(f"  - bulk_energy, bulk_energy_per_atom, bulk_structure")
    print(f"  - slab_structures, relaxed_slabs, slab_energies")
    print(f"  - surface_energies (gamma_eV_A2, gamma_J_m2)")
    print(f"\nExpected surface energies for Au:")
    print(f"  (111): ~0.79 J/m² (most stable)")
    print(f"  (100): ~0.94 J/m²")
    print(f"  (110): ~1.33 J/m²")
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
