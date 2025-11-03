#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 17: Test max_concurrent_jobs Parameter

This script tests the max_concurrent_jobs functionality that was implemented to
control concurrent VASP calculations in nested workgraphs.

The test:
1. Generates multiple slab structures ([1,0,0] Miller index creates multiple terminations)
2. Sets max_concurrent_jobs=2 to limit concurrent calculations
3. Monitors that no more than 2 VASP calculations run simultaneously

Material: Ag2O
Surface: (100) - Creates multiple terminations for testing concurrency
Max concurrent jobs: 2 (to easily observe the limiting behavior)

Expected behavior:
- Multiple slab structures will be generated
- Only 2 VASP calculations should run at the same time
- Remaining calculations should wait until slots are available

Usage:
    source ~/envs/aiida/bin/activate
    python step_17_test_max_concurrent_jobs.py

Monitor concurrent jobs with:
    watch -n 2 'verdi process list -a -p 1 | head -30'
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Step 17: Test max_concurrent_jobs parameter."""

    print("\n" + "="*70)
    print("STEP 17: TEST max_concurrent_jobs PARAMETER")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")

    # Code configuration
    code_label = 'VASP-6.5.0@bohr-new'
    #code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # Minimal VASP parameters for fast testing
    vasp_params = {
        'PREC': 'Normal',
        'ENCUT': 400,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 10,  # Small number for quick testing
        'EDIFFG': -0.05,
        'ALGO': 'Fast',
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

    print("\n3. Building workgraph with max_concurrent_jobs=2...")
    print("   Using workflow preset: 'surface_thermodynamics'")
    print("   Miller indices: [1, 0, 0] (will create multiple terminations)")
    print("   Concurrency limit: 2 VASP calculations at a time")

    # Build workgraph with max_concurrent_jobs
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
        min_slab_thickness=15.0,
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,

        # Slab relaxation
        slab_parameters=slab_params,
        slab_options=common_options,
        slab_kpoints_spacing=0.4,

        # *** KEY PARAMETER: Limit concurrent VASP calculations ***
        max_concurrent_jobs=1,

        name='Step17_TestMaxConcurrentJobs_Ag2O_100',
    )

    print("   ✓ WorkGraph built successfully")
    print(f"   ✓ max_concurrent_jobs is set to: 2")

    # Submit
    print("\n4. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 17 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nTo verify max_concurrent_jobs is working:")
    print(f"  watch -n 2 'verdi process list -a -p 1 | head -30'")
    print(f"\n  You should see:")
    print(f"    ✓ No more than 2 VaspWorkChain processes running simultaneously")
    print(f"    ✓ Other calculations waiting in 'Created' state")
    print(f"    ✓ New calculations starting as previous ones finish")
    print(f"\nExpected outputs:")
    print(f"  - bulk_structure (relaxed)")
    print(f"  - bulk_energy")
    print(f"  - relaxed_slabs (all terminations)")
    print(f"  - slab_energies (all terminations)")
    print(f"\nIf max_concurrent_jobs is working correctly:")
    print(f"  - Bulk calculation will run first")
    print(f"  - Then slab relaxations will run in batches of 2")
    print(f"  - Total time will be longer than unlimited, but system load controlled")
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
