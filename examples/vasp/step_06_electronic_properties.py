#!/usr/bin/env python
"""
STEP 6: Electronic Properties (DOS and Band Structure)

This script tests:
- Bulk relaxation
- Electronic structure calculation (DOS and band structure)

This demonstrates the electronic properties workflow for bulk materials,
using the vasp.v2.bands workchain from aiida-vasp.

Material: Ag2O (cuprite structure)

Usage:
    # Activate your AiiDA environment first
    python step_06_electronic_properties.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
try:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, COMMON_OPTIONS, POTENTIAL_MAPPING_AG2O
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
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal',
    }
    POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_electronic_properties_defaults


def main():
    """Step 6: Test electronic properties calculation."""

    print("\n" + "="*70)
    print("STEP 6: ELECTRONIC PROPERTIES (DOS AND BANDS)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    # Get electronic properties defaults
    print("\n4. Setting up electronic properties parameters...")
    elec_defaults = get_electronic_properties_defaults()

    print("   Band structure: high-symmetry path (auto-determined)")
    print("   DOS: total and projected components")

    print("\n5. Building workgraph...")
    print("   Preset: 'electronic_structure_bulk_only'")

    # Build workgraph using preset
    wg = build_core_workgraph(
        workflow_preset='electronic_structure_bulk_only',

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

        # Electronic properties
        bands_parameters=elec_defaults,
        bands_options=COMMON_OPTIONS,
        band_settings=elec_defaults['band_settings'],

        # Concurrency control
        max_concurrent_jobs=4,

        name='Step06_ElectronicProperties_Ag2O',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n6. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 6 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(bulk)")
    print(f"  - bulk_structure: Relaxed structure")
    print(f"  - bulk_bands: Band structure data")
    print(f"  - bulk_dos: Density of states")
    print(f"  - bulk_primitive_structure: Primitive cell used")
    print(f"  - bulk_seekpath_parameters: Symmetry info")
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
