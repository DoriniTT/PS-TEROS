#!/usr/bin/env python
"""
STEP 14: Concurrency Control Example

Demonstrates the max_concurrent_jobs parameter for controlling
how many VASP calculations run simultaneously.

Modes available (modify MAX_CONCURRENT_JOBS below):
- 1: Serial mode - only 1 VASP at a time
- 4: Limited mode - max 4 VASP calculations at once (default)
- None: Unlimited mode - full parallel execution

Material: Ag2O
Surface: (100)

Usage:
    # Activate your AiiDA environment first
    python step_14_concurrent_limit.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
try:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, SLAB_PARAMS, COMMON_OPTIONS,
        POTENTIAL_MAPPING_AG2O, POTENTIAL_MAPPING_AG, POTENTIAL_MAPPING_O2
    )
except ImportError:
    AIIDA_PROFILE = 'presto'
    VASP_CODE = 'VASP-6.4.1@cluster'
    POTENTIAL_FAMILY = 'PBE'
    BULK_PARAMS = {
        'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6, 'ISMEAR': 0,
        'SIGMA': 0.05, 'IBRION': 2, 'ISIF': 3, 'NSW': 100, 'EDIFFG': -0.01,
        'ALGO': 'Normal', 'LREAL': 'Auto', 'LWAVE': False, 'LCHARG': False,
    }
    SLAB_PARAMS = BULK_PARAMS.copy()
    SLAB_PARAMS['ISIF'] = 2
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal',
    }
    POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
    POTENTIAL_MAPPING_AG = {'Ag': 'Ag'}
    POTENTIAL_MAPPING_O2 = {'O': 'O'}

# Concurrency setting - modify this value
# Options: 1 (serial), 4 (limited, default), None (unlimited)
MAX_CONCURRENT_JOBS = 4
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 14: Test max_concurrent_jobs parameter."""

    print("\n" + "="*70)
    print("STEP 14: CONCURRENCY CONTROL EXAMPLE")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    # Miller indices
    miller_indices = [1, 0, 0]  # (100) surface

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")
    print(f"   Surface: ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]})")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    # Determine mode description
    if MAX_CONCURRENT_JOBS == 1:
        mode_name = 'Serial'
        description = 'Only 1 VASP calculation at a time'
    elif MAX_CONCURRENT_JOBS is None:
        mode_name = 'Unlimited'
        description = 'No limit (full parallel)'
    else:
        mode_name = 'Limited'
        description = f'Max {MAX_CONCURRENT_JOBS} VASP calculations at once'

    print(f"\n4. Concurrency configuration:")
    print(f"   Mode: {mode_name}")
    print(f"   max_concurrent_jobs: {MAX_CONCURRENT_JOBS}")
    print(f"   Description: {description}")

    print("\n5. Building workgraph...")

    # Oxygen parameters (ISIF=2 for molecule in box)
    oxygen_params = BULK_PARAMS.copy()
    oxygen_params['ISIF'] = 2

    # Build workgraph
    wg = build_core_workgraph(
        max_concurrent_jobs=MAX_CONCURRENT_JOBS,

        # Enable slab generation + relaxation
        relax_slabs=True,
        compute_thermodynamics=True,
        compute_cleavage=False,

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        metal_name='Ag.cif',
        oxygen_name='O2.cif',

        # VASP code
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.04,
        clean_workdir=False,

        # Bulk (Ag2O)
        bulk_potential_mapping=POTENTIAL_MAPPING_AG2O,
        bulk_parameters=BULK_PARAMS,
        bulk_options=COMMON_OPTIONS,

        # Metal (Ag)
        metal_potential_mapping=POTENTIAL_MAPPING_AG,
        metal_parameters=BULK_PARAMS,
        metal_options=COMMON_OPTIONS,

        # Oxygen (O2)
        oxygen_potential_mapping=POTENTIAL_MAPPING_O2,
        oxygen_parameters=oxygen_params,
        oxygen_options=COMMON_OPTIONS,

        # Slab generation
        miller_indices=miller_indices,
        min_slab_thickness=18.0,
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,

        # Slabs
        slab_parameters=SLAB_PARAMS,
        slab_options=COMMON_OPTIONS,

        name=f'Step14_ConcurrencyControl_{mode_name}_Ag2O',
    )

    print("   WorkGraph built successfully")
    print(f"   WorkGraph.max_number_jobs = {wg.max_number_jobs}")

    # Submit
    print("\n6. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 14 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"Mode: {mode_name} (max_concurrent_jobs={MAX_CONCURRENT_JOBS})")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"  watch -n 5 'verdi process list -a | grep VASP'")

    if MAX_CONCURRENT_JOBS == 1:
        print(f"\nSerial mode: Only 1 VASP running at any time")
    elif MAX_CONCURRENT_JOBS is None:
        print(f"\nUnlimited mode: All VASP jobs launch in parallel")
    else:
        print(f"\nLimited mode: Max {MAX_CONCURRENT_JOBS} VASP jobs running at once")

    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
