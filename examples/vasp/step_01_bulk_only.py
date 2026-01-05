#!/usr/bin/env python
"""
STEP 1: Bulk Relaxation Only

This script tests the most basic functionality of PS-TEROS:
relaxing a bulk structure using VASP through AiiDA.

What this tests:
- Loading AiiDA profile
- Setting up VASP code and potential
- Bulk structure relaxation (ISIF=3: cell + ions)
- Basic WorkGraph submission

Material: Ag2O (cuprite structure)

Usage:
    # Activate your AiiDA environment first
    python step_01_bulk_only.py

    # Or with explicit environment activation:
    source ~/envs/aiida/bin/activate && python step_01_bulk_only.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
# Option 1: Import from shared config (recommended)
try:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, COMMON_OPTIONS, POTENTIAL_MAPPING_AG2O
    )
except ImportError:
    # Option 2: Define locally if config.py not available
    AIIDA_PROFILE = 'presto'
    VASP_CODE = 'VASP-6.4.1@cluster'
    POTENTIAL_FAMILY = 'PBE'
    BULK_PARAMS = {
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
    COMMON_OPTIONS = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
        'queue_name': 'normal',
    }
    POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 1: Test bulk relaxation only."""

    print("\n" + "="*70)
    print("STEP 1: BULK RELAXATION ONLY")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structure:")
    print(f"   {structures_dir}/ag2o.cif")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")
    print(f"   ENCUT: {BULK_PARAMS['ENCUT']} eV")
    print(f"   ISIF: {BULK_PARAMS['ISIF']} (full cell + ions)")

    print("\n4. Building workgraph...")
    print("   Preset: 'bulk_only'")

    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='bulk_only',

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',

        # VASP code
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.04,
        clean_workdir=False,

        # Bulk calculation
        bulk_potential_mapping=POTENTIAL_MAPPING_AG2O,
        bulk_parameters=BULK_PARAMS,
        bulk_options=COMMON_OPTIONS,

        # Concurrency control
        max_concurrent_jobs=4,

        name='Step01_BulkOnly_Ag2O',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 1 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: Total energy of relaxed Ag2O")
    print(f"  - bulk_structure: Relaxed structure")
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
