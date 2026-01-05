#!/usr/bin/env python
"""
STEP 7: AIMD Simulation

This script tests:
- Bulk relaxation
- Slab generation (no relaxation)
- AIMD (Ab Initio Molecular Dynamics) simulation on slabs

This demonstrates dynamic properties at finite temperature.
NOTE: With 'aimd_only' preset, slabs are NOT relaxed before AIMD.

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence (equilibration + production)

WARNING: AIMD is computationally expensive and will take many hours.

Usage:
    # Activate your AiiDA environment first
    python step_07_aimd_simulation.py
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
        POTENTIAL_MAPPING_AG2O, SLAB_DEFAULTS
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
    SLAB_DEFAULTS = {
        'min_slab_thickness': 15.0,
        'min_vacuum_thickness': 15.0,
        'lll_reduce': True,
        'center_slab': True,
        'symmetrize': True,
        'primitive': True,
    }
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_aimd_defaults


def main():
    """Step 7: Test AIMD simulation."""

    print("\n" + "="*70)
    print("STEP 7: AIMD SIMULATION")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    # Miller indices
    miller_indices = [1, 1, 1]  # (111) surface

    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")
    print(f"   Surface: ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]})")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    # AIMD configuration
    print("\n4. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    for i, stage in enumerate(aimd_sequence):
        print(f"   Stage {i+1}: {stage['temperature']} K, {stage['steps']} steps")

    aimd_params = get_aimd_defaults()

    print("\n5. Building workgraph...")
    print("   Preset: 'aimd_only'")
    print("   Note: Slabs are NOT relaxed before AIMD with this preset")

    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='aimd_only',

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',

        # VASP code
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.04,
        clean_workdir=False,

        # Bulk
        bulk_potential_mapping=POTENTIAL_MAPPING_AG2O,
        bulk_parameters=BULK_PARAMS,
        bulk_options=COMMON_OPTIONS,

        # Slab generation (smaller for AIMD)
        miller_indices=miller_indices,
        min_slab_thickness=15.0,  # Smaller for AIMD
        min_vacuum_thickness=15.0,
        lll_reduce=SLAB_DEFAULTS['lll_reduce'],
        center_slab=SLAB_DEFAULTS['center_slab'],
        symmetrize=SLAB_DEFAULTS['symmetrize'],
        primitive=SLAB_DEFAULTS['primitive'],

        # Slab relaxation (not used in aimd_only but needed for param)
        slab_parameters=SLAB_PARAMS,
        slab_options=COMMON_OPTIONS,
        slab_kpoints_spacing=0.04,

        # AIMD
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=COMMON_OPTIONS,
        aimd_potential_mapping=POTENTIAL_MAPPING_AG2O,
        aimd_kpoints_spacing=0.5,  # Coarser for AIMD

        # Concurrency control
        max_concurrent_jobs=4,

        name=f'Step07_AIMD_Ag2O_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n6. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 7 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nWARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures")
    print(f"  - aimd_results: Nested namespace with AIMD data")
    print(f"    - For each slab and stage:")
    print(f"      - trajectory, structure, energy, remote, retrieved")
    print(f"\nAIMD trajectories can be analyzed for:")
    print(f"  - Temperature evolution")
    print(f"  - Atomic diffusion")
    print(f"  - Structural stability")
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
