#!/usr/bin/env python
"""
STEP 10: Custom Workflow Configuration

This script demonstrates how to build a custom workflow by setting individual flags
instead of using predefined presets. This gives maximum flexibility for combining
different workflow components.

Example: Compute surface thermodynamics + cleavage energies, but skip relaxation energies

Material: Ag2O
Surface: (100)

Usage:
    # Activate your AiiDA environment first
    python step_10_custom_workflow.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
try:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, METAL_PARAMS, SLAB_PARAMS, COMMON_OPTIONS,
        POTENTIAL_MAPPING_AG2O, POTENTIAL_MAPPING_AG, POTENTIAL_MAPPING_O2,
        SLAB_DEFAULTS
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
    METAL_PARAMS = BULK_PARAMS.copy()
    METAL_PARAMS.update({'ISMEAR': 1, 'SIGMA': 0.2})
    SLAB_PARAMS = BULK_PARAMS.copy()
    SLAB_PARAMS['ISIF'] = 2
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal',
    }
    POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
    POTENTIAL_MAPPING_AG = {'Ag': 'Ag'}
    POTENTIAL_MAPPING_O2 = {'O': 'O'}
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
    """Step 10: Test custom workflow configuration."""

    print("\n" + "="*70)
    print("STEP 10: CUSTOM WORKFLOW CONFIGURATION")
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

    print("\n4. Building custom workflow...")
    print("   Custom configuration (no preset):")
    print("     [x] relax_slabs: True")
    print("     [x] compute_thermodynamics: True")
    print("     [x] compute_cleavage: True")
    print("     [ ] compute_relaxation_energy: False")
    print("     [ ] compute_electronic_properties_bulk: False")
    print("     [ ] compute_electronic_properties_slabs: False")
    print("     [ ] run_aimd: False")

    # Oxygen parameters (ISIF=2 for molecule in box)
    oxygen_params = BULK_PARAMS.copy()
    oxygen_params['ISIF'] = 2

    # Build workgraph with custom flags (no preset, set flags individually)
    wg = build_core_workgraph(
        # NO workflow_preset - this triggers custom workflow mode
        # Individual flags set explicitly
        relax_slabs=True,
        compute_thermodynamics=True,
        compute_cleavage=True,
        compute_relaxation_energy=False,  # Skip relaxation energies
        compute_electronic_properties_bulk=False,
        compute_electronic_properties_slabs=False,
        run_aimd=False,

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
        metal_parameters=METAL_PARAMS,
        metal_options=COMMON_OPTIONS,

        # Oxygen (O2)
        oxygen_potential_mapping=POTENTIAL_MAPPING_O2,
        oxygen_parameters=oxygen_params,
        oxygen_options=COMMON_OPTIONS,

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

        # Thermodynamics sampling
        thermodynamics_sampling=100,

        # Concurrency control
        max_concurrent_jobs=4,

        name=f'Step10_CustomWorkflow_Ag2O_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 10 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, metal_energy, oxygen_energy")
    print(f"  - formation_enthalpy")
    print(f"  - slab_structures (all terminations)")
    print(f"  - slab_energies (relaxed)")
    print(f"  - cleavage_energies")
    print(f"  - surface_energies: gamma(mu_O)")
    print(f"\nNote: Relaxation energies NOT computed (custom choice)")
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
