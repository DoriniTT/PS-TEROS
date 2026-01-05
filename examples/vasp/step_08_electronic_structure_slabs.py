#!/usr/bin/env python
"""
STEP 8: Electronic Structure for Slabs Only

This script tests:
- Bulk relaxation
- Slab generation and relaxation
- Electronic structure calculation (DOS and bands) for slabs only

This demonstrates the electronic properties workflow for surface structures.

Material: Ag2O
Surface: (111)

Usage:
    # Activate your AiiDA environment first
    python step_08_electronic_structure_slabs.py
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
from teros.core.builders import get_slab_electronic_properties_defaults


def main():
    """Step 8: Test electronic properties calculation for slabs only."""

    print("\n" + "="*70)
    print("STEP 8: ELECTRONIC STRUCTURE FOR SLABS ONLY")
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

    # Get slab electronic properties defaults
    print("\n4. Setting up slab electronic properties parameters...")
    slab_elec_defaults = get_slab_electronic_properties_defaults()

    print("   Band structure: high-symmetry paths for slabs")
    print("   DOS: total and projected components")

    print("\n5. Building workgraph...")
    print("   Preset: 'electronic_structure_slabs_only'")

    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='electronic_structure_slabs_only',

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

        # Slab electronic properties
        slab_bands_parameters=slab_elec_defaults,
        slab_bands_options=COMMON_OPTIONS,
        slab_band_settings=slab_elec_defaults['band_settings'],

        # Concurrency control
        max_concurrent_jobs=4,

        name=f'Step08_ElectronicStructure_Slabs_Ag2O_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n6. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 8 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures: All terminations")
    print(f"  - slab_energies: Relaxed slab energies")
    print(f"  - slab_bands: Band structure for each slab")
    print(f"  - slab_dos: Density of states for each slab")
    print(f"  - slab_primitive_structures: Primitive cells used")
    print(f"  - slab_seekpath_parameters: Symmetry info")
    print(f"\nVisualize with:")
    print(f"  verdi data bands show <PK>")
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
