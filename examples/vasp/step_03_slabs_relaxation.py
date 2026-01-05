#!/usr/bin/env python
"""
STEP 3: Slab Generation and Relaxation

This script tests:
- Bulk relaxation
- Slab generation from bulk (using pymatgen)
- Slab relaxation (ISIF=2: ions only)
- Relaxation energy calculation (E_relaxed - E_unrelaxed)

This demonstrates the slab workflow without surface thermodynamics.

Material: Ag2O
Surface: (100) - multiple terminations

Usage:
    # Activate your AiiDA environment first
    python step_03_slabs_relaxation.py
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
        'min_slab_thickness': 18.0,
        'min_vacuum_thickness': 15.0,
        'lll_reduce': True,
        'center_slab': True,
        'symmetrize': True,
        'primitive': True,
    }
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 3: Test slab generation and relaxation."""

    print("\n" + "="*70)
    print("STEP 3: SLAB GENERATION AND RELAXATION")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    # Miller indices for slab generation
    miller_indices = [1, 0, 0]  # (100) surface

    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")
    print(f"   Surface: ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]})")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    print("\n4. Building workgraph...")
    print("   Preset: 'relaxation_energy_only'")
    print("   Steps:")
    print("     1. Relax bulk Ag2O")
    print(f"     2. Generate ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]}) slabs")
    print("     3. Run SCF on unrelaxed slabs")
    print("     4. Relax slabs (ISIF=2)")
    print("     5. Calculate delta_E = E_relaxed - E_unrelaxed")

    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='relaxation_energy_only',

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

        # Slab generation
        miller_indices=miller_indices,
        min_slab_thickness=SLAB_DEFAULTS['min_slab_thickness'],
        min_vacuum_thickness=SLAB_DEFAULTS['min_vacuum_thickness'],
        lll_reduce=SLAB_DEFAULTS['lll_reduce'],
        center_slab=SLAB_DEFAULTS['center_slab'],
        symmetrize=SLAB_DEFAULTS['symmetrize'],
        primitive=SLAB_DEFAULTS['primitive'],

        # Slab relaxation
        slab_parameters=SLAB_PARAMS,
        slab_options=COMMON_OPTIONS,
        slab_kpoints_spacing=0.04,

        # Concurrency control
        max_concurrent_jobs=4,

        name=f'Step03_SlabRelaxation_Ag2O_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 3 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(bulk Ag2O)")
    print(f"  - bulk_structure: Relaxed bulk")
    print(f"  - slab_structures: Generated terminations")
    print(f"  - unrelaxed_slab_energies: SCF energies")
    print(f"  - slab_energies: Relaxed energies")
    print(f"  - relaxation_energies: delta_E = E_relaxed - E_unrelaxed")
    print(f"\nNote: Number of terminations depends on crystal symmetry")
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
