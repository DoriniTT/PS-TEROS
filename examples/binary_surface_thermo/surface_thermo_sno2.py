#!/usr/bin/env python
"""
Surface Thermodynamics for SnO2 (110)

This script calculates surface thermodynamics for tin dioxide using the
'surface_thermodynamics' preset, which includes:
- Bulk and reference relaxations
- Formation enthalpy calculation
- Slab generation and relaxation
- Surface energies as function of chemical potential

Material: SnO2 (rutile tin dioxide)
Surface: (110) - the most stable rutile surface
References: Sn (beta-tin metal), O2 (molecule)

Note: Cleavage and relaxation energies are OPTIONAL. To include them:
    compute_cleavage=True
    compute_relaxation_energy=True

Usage:
    # Activate your AiiDA environment first
    python surface_thermo_sno2.py
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
        POTENTIAL_MAPPING_SNO2, POTENTIAL_MAPPING_SN, POTENTIAL_MAPPING_O2,
        SLAB_DEFAULTS
    )
except ImportError:
    AIIDA_PROFILE = 'presto'
    VASP_CODE = 'VASP-VTST-6.4.3@bohr'
    POTENTIAL_FAMILY = 'PBE'

    BULK_PARAMS = {
        'PREC': 'Accurate', 'ENCUT': 500, 'EDIFF': 1e-5, 'ISMEAR': 0,
        'SIGMA': 0.05, 'IBRION': 2, 'ISIF': 3, 'NSW': 500, 'EDIFFG': -0.01,
        'ALGO': 'Normal', 'LREAL': 'Auto', 'LWAVE': False, 'LCHARG': False,
    }
    METAL_PARAMS = BULK_PARAMS.copy()
    METAL_PARAMS.update({'ISMEAR': 1, 'SIGMA': 0.2})  # Metallic Sn
    SLAB_PARAMS = BULK_PARAMS.copy()
    SLAB_PARAMS['ISIF'] = 2  # Fix cell, relax ions

    # Bohr cluster: 40 cores per node, PBS Pro scheduler
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 40},
        'queue_name': 'par40',
    }

    # Use Sn_d for better d-electron treatment (14 valence electrons)
    POTENTIAL_MAPPING_SNO2 = {'Sn': 'Sn_d', 'O': 'O'}
    POTENTIAL_MAPPING_SN = {'Sn': 'Sn_d'}
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
    """Run surface thermodynamics workflow for SnO2 (110)."""

    print("\n" + "="*70)
    print("SURFACE THERMODYNAMICS: SnO2 (110)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    # Miller indices for slab generation - (110) is the most stable rutile surface
    miller_indices = [1, 1, 0]

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/sno2.vasp")
    print(f"   Metal:  {structures_dir}/Sn.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")
    print(f"   Surface: ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]})")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")
    print(f"   Sn potential: Sn_d (14 valence electrons)")

    print("\n4. Building workgraph...")
    print("   Preset: 'surface_thermodynamics'")
    print("   Calculates:")
    print("     - Formation enthalpy")
    print("     - Surface energies: gamma(mu_O)")
    print("   Optional (not included by default):")
    print("     - Cleavage energies")
    print("     - Relaxation energies")

    # Oxygen parameters (ISIF=2 for molecule in box)
    oxygen_params = BULK_PARAMS.copy()
    oxygen_params['ISIF'] = 2

    # Build workgraph using surface_thermodynamics preset
    wg = build_core_workgraph(
        workflow_preset='surface_thermodynamics',

        # Structures
        structures_dir=structures_dir,
        bulk_name='sno2.vasp',
        metal_name='Sn.cif',
        oxygen_name='O2.cif',

        # VASP code
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        kpoints_spacing=0.04,
        clean_workdir=False,

        # Bulk (SnO2)
        bulk_potential_mapping=POTENTIAL_MAPPING_SNO2,
        bulk_parameters=BULK_PARAMS,
        bulk_options=COMMON_OPTIONS,

        # Metal (Sn)
        metal_potential_mapping=POTENTIAL_MAPPING_SN,
        metal_parameters=METAL_PARAMS,
        metal_options=COMMON_OPTIONS,

        # Oxygen (O2)
        oxygen_potential_mapping=POTENTIAL_MAPPING_O2,
        oxygen_parameters=oxygen_params,
        oxygen_options=COMMON_OPTIONS,

        # Slab generation - (110) surface
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

        # Thermodynamics sampling (number of chemical potential points)
        thermodynamics_sampling=100,

        # Concurrency control
        max_concurrent_jobs=4,

        name=f'SurfaceThermodynamics_SnO2_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, metal_energy, oxygen_energy")
    print(f"  - formation_enthalpy (eV/formula unit)")
    print(f"  - slab_structures (all (110) terminations)")
    print(f"  - unrelaxed_slab_energies, slab_energies")
    print(f"  - surface_energies: gamma(mu_O) in J/m^2")
    print(f"\nSurface energies show stability vs oxygen chemical potential")
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
