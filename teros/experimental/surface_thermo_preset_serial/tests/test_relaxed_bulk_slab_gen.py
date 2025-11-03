#!/home/thiagotd/envs/aiida/bin/python
"""
Test Script: Slab Generation from RELAXED Bulk

This script tests the NEW implementation that:
1. Relaxes the bulk structure
2. Generates slabs from the RELAXED bulk (not input bulk)
3. Runs SCF on generated slabs

Material: Ag2O
Miller indices: (100)
Concurrency limit: 2 jobs

Usage:
    source ~/envs/aiida/bin/activate
    python test_relaxed_bulk_slab_gen.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.experimental.surface_thermo_preset_serial import surface_thermodynamics_serial_workgraph


def main():
    """Test slab generation from relaxed bulk."""

    print("\n" + "="*70)
    print("TEST: SLAB GENERATION FROM RELAXED BULK")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structures directory: {structures_dir}")

    # Code configuration
    code_label = 'VASP-6.5.0@bohr-new'
    code = orm.load_code(code_label)
    potential_family = 'PBE'

    print(f"\n3. VASP Configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # Common options
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # VASP parameters - VERY LIGHT for quick testing
    bulk_params = {
        'PREC': 'Normal',
        'ENCUT': 300,
        'EDIFF': 1e-3,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,  # Full relaxation
        'NSW': 50,  # Few steps for speed
        'EDIFFG': -0.5,
        'ALGO': 'Fast',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    slab_scf_params = {
        'PREC': 'Normal',
        'ENCUT': 300,
        'EDIFF': 1e-3,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 0,  # SCF only
        'ALGO': 'Fast',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Only one Miller index for simplicity
    miller_indices = [
        (1, 0, 0),
    ]

    print(f"\n4. Workflow Configuration:")
    print(f"   Miller indices: {miller_indices}")
    print(f"   Bulk: NSW={bulk_params['NSW']}, ENCUT={bulk_params['ENCUT']}")
    print(f"   Slabs: NSW={slab_scf_params['NSW']}, ENCUT={slab_scf_params['ENCUT']}")
    print(f"   Slab generation: From RELAXED bulk (not input bulk)")

    print(f"\n5. Building workgraph...")

    try:
        wg = surface_thermodynamics_serial_workgraph(
            # Structures
            structures_dir=structures_dir,
            bulk_name='ag2o.cif',

            # Code
            code_label=code_label,
            potential_family=potential_family,
            kpoints_spacing=1.0,  # Coarse
            clean_workdir=False,

            # Bulk parameters
            bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            bulk_parameters=bulk_params,
            bulk_options=common_options,

            # Slab generation from RELAXED bulk
            miller_indices=miller_indices,
            min_slab_thickness=8,
            min_vacuum_thickness=10,
            lll_reduce=False,
            center_slab=True,
            primitive=True,
            in_unit_planes=False,
            max_normal_search=1,
            symmetrize=False,

            # Slab parameters (SCF only for quick test)
            slab_parameters=slab_scf_params,
            slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            slab_options=common_options,
            slab_kpoints_spacing=1.0,

            # Workflow flags - MINIMAL for testing
            relax_slabs=False,  # Only SCF, no relaxation
            compute_thermodynamics=False,  # No thermodynamics
            compute_relaxation_energy=False,
        )

        print("   ✓ WorkGraph built successfully")

    except Exception as e:
        print(f"   ✗ Error building workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None

    # Set concurrent job limit
    print(f"\n6. Configuring concurrency...")
    wg.max_number_jobs = 2
    print(f"   ✓ max_number_jobs = 2")

    # Submit
    print(f"\n7. Submitting to AiiDA daemon...")

    try:
        result = wg.submit()
        pk = result.pk if hasattr(result, 'pk') else result

        print(f"\n{'='*70}")
        print("TEST SUBMITTED SUCCESSFULLY")
        print(f"{'='*70}")
        print(f"\nWorkGraph PK: {pk}")
        print(f"\nMonitor with:")
        print(f"  verdi process show {pk}")
        print(f"  verdi process report {pk}")
        print(f"  watch -n 2 'verdi process list -p 1'")

        print(f"\nExpected workflow:")
        print(f"  1. Bulk relaxation (1 VASP job, NSW=50)")
        print(f"  2. Slab generation from RELAXED bulk (calcfunction)")
        print(f"  3. SCF calculations on generated slabs (N slabs, max 2 concurrent)")

        print(f"\nKey Test Points:")
        print(f"  ✓ Slabs generated from RELAXED bulk (not input)")
        print(f"  ✓ Slab generation waits for bulk relaxation")
        print(f"  ✓ Dynamic slab outputs connected to VASP tasks")
        print(f"  ✓ Flat-graph architecture maintained")
        print(f"  ✓ max_number_jobs=2 controls concurrency")

        print(f"\nTo verify success:")
        print(f"  1. Check bulk relaxation completes")
        print(f"  2. Check slab generation task runs (generates slabs from relaxed bulk)")
        print(f"  3. Check slab SCF tasks start after slabgeneration")
        print(f"  4. Verify only 2 VASP jobs run concurrently")
        print(f"{'='*70}\n")

        return pk

    except Exception as e:
        print(f"\n✗ Error submitting workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == '__main__':
    try:
        pk = main()
        if pk is None:
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
