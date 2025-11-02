#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 14: Concurrency Control Example

Demonstrates the max_concurrent_jobs parameter for controlling
how many VASP calculations run simultaneously.

This example shows three modes:
1. Serial mode (max_concurrent_jobs=1)
2. Limited concurrency (max_concurrent_jobs=4, default)
3. Unlimited parallel (max_concurrent_jobs=None)

Material: Ag2O
Miller indices: (1,0,0), (1,1,0)

Usage:
    source ~/envs/aiida/bin/activate
    python step_14_concurrent_limit.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 14: Test max_concurrent_jobs parameter."""

    print("\n" + "="*70)
    print("STEP 14: CONCURRENCY CONTROL EXAMPLE")
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
    potential_family = 'PBE.64'

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

    slab_params = vasp_params.copy()
    slab_params['ISIF'] = 2

    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Choose mode
    print("\n3. Concurrency modes:")
    print("   A. Serial mode (max_concurrent_jobs=1)")
    print("   B. Limited mode (max_concurrent_jobs=4, default)")
    print("   C. Unlimited mode (max_concurrent_jobs=None)")

    mode = input("\nSelect mode (A/B/C) [B]: ").strip().upper() or 'B'

    mode_map = {
        'A': (1, 'Serial', 'Only 1 VASP at a time'),
        'B': (4, 'Limited', 'Max 4 VASP calculations at once (default)'),
        'C': (None, 'Unlimited', 'No limit (full parallel)'),
    }

    max_concurrent_jobs, mode_name, description = mode_map.get(mode, mode_map['B'])

    print(f"\n4. Building workgraph...")
    print(f"   Mode: {mode_name}")
    print(f"   max_concurrent_jobs: {max_concurrent_jobs}")
    print(f"   Description: {description}")

    # Build workgraph
    wg = build_core_workgraph(
        max_concurrent_jobs=max_concurrent_jobs,  # CONCURRENCY CONTROL

        # Enable slab generation + relaxation
        relax_slabs=True,
        miller_indices=[(1, 0, 0), (1, 1, 0)],
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,

        # Enable thermodynamics
        compute_thermodynamics=True,
        compute_cleavage=False,

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

        # Slabs
        slab_parameters=slab_params.copy(),
        slab_options=common_options,
        slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},

        name=f'Step14_ConcurrencyControl_{mode_name}_Ag2O',
    )

    print("   ✓ WorkGraph built successfully")
    print(f"   WorkGraph.max_number_jobs = {wg.max_number_jobs}")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 14 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"Mode: {mode_name} (max_concurrent_jobs={max_concurrent_jobs})")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"  watch -n 5 'verdi process list -a | grep VASP'")

    if max_concurrent_jobs == 1:
        print(f"\n✓ Serial mode: Only 1 VASP running at any time")
    elif max_concurrent_jobs is None:
        print(f"\n✓ Unlimited mode: All VASP jobs launch in parallel")
    else:
        print(f"\n✓ Limited mode: Max {max_concurrent_jobs} VASP jobs running at once")

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
