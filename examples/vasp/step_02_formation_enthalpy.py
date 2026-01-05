#!/usr/bin/env python
"""
STEP 2: Formation Enthalpy Calculation

This script calculates the formation enthalpy of Ag2O using:
- Bulk Ag2O relaxation
- Reference calculations (Ag metal, O2 molecule)
- Formation enthalpy: delta_Hf = E(Ag2O) - 2*E(Ag) - 0.5*E(O2)

What this tests:
- Multiple structure calculations in parallel
- Formation enthalpy calculation
- Reference energy handling

Material: Ag2O
References: Ag (fcc metal), O2 (molecule in box)

Usage:
    # Activate your AiiDA environment first
    python step_02_formation_enthalpy.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
try:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, METAL_PARAMS, COMMON_OPTIONS,
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
    METAL_PARAMS = BULK_PARAMS.copy()
    METAL_PARAMS.update({'ISMEAR': 1, 'SIGMA': 0.2})
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal',
    }
    POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
    POTENTIAL_MAPPING_AG = {'Ag': 'Ag'}
    POTENTIAL_MAPPING_O2 = {'O': 'O'}
# ==============================================================================

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 2: Test formation enthalpy calculation."""

    print("\n" + "="*70)
    print("STEP 2: FORMATION ENTHALPY CALCULATION")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    print("\n4. Building workgraph...")
    print("   Preset: 'formation_enthalpy_only'")
    print("   Formula: delta_Hf = E(Ag2O) - 2*E(Ag) - 0.5*E(O2)")

    # Oxygen parameters (ISIF=2 for molecule in box)
    oxygen_params = BULK_PARAMS.copy()
    oxygen_params['ISIF'] = 2

    # Build workgraph
    wg = build_core_workgraph(
        workflow_preset='formation_enthalpy_only',

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

        # Metal (Ag) - uses Methfessel-Paxton smearing
        metal_potential_mapping=POTENTIAL_MAPPING_AG,
        metal_parameters=METAL_PARAMS,
        metal_options=COMMON_OPTIONS,

        # Oxygen (O2 molecule)
        oxygen_potential_mapping=POTENTIAL_MAPPING_O2,
        oxygen_parameters=oxygen_params,
        oxygen_options=COMMON_OPTIONS,

        # Concurrency control
        max_concurrent_jobs=4,

        name='Step02_FormationEnthalpy_Ag2O',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 2 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: E(Ag2O)")
    print(f"  - metal_energy: E(Ag)")
    print(f"  - oxygen_energy: E(O2)")
    print(f"  - formation_enthalpy: delta_Hf in eV/formula unit")
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
