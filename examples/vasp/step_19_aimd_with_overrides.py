#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 19: AIMD with Per-Structure Overrides

Demonstrates the override system in standalone AIMD module:
- Structure-level overrides
- Stage-level overrides
- Matrix-level overrides
- Priority order verification

Material: Ag (silver bulk)
AIMD: 2-stage sequence
Structures: 2 identical structures with different INCAR settings

Usage:
    source ~/envs/aiida/bin/activate
    python step_19_aimd_with_overrides.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.aimd import build_aimd_workgraph
from ase.io import read


def main():
    """Step 19: Test AIMD override system."""

    print("\n" + "="*70)
    print("STEP 19: AIMD WITH PER-STRUCTURE OVERRIDES")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    ag_cif = os.path.join(structures_dir, 'Ag.cif')

    print(f"\n2. Loading structures:")
    print(f"   Bulk: {ag_cif}")

    # Load structures using ASE
    ag_ase = read(ag_cif)
    ag_structure1 = orm.StructureData(ase=ag_ase)
    ag_structure2 = orm.StructureData(ase=ag_ase)

    print(f"   ✓ Loaded 2 Ag structures")
    print(f"     - Structure 1: {len(ag_structure1.sites)} atoms")
    print(f"     - Structure 2: {len(ag_structure2.sites)} atoms")

    # Code configuration
    code_label = 'VASP-6.5.1@cluster02'
    potential_family = 'PBE'

    print(f"\n3. VASP configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # AIMD stages using VASP-native parameter names
    print("\n4. AIMD configuration:")
    aimd_stages = [
        {'TEBEG': 300, 'NSW': 50},   # Equilibration
        {'TEBEG': 300, 'NSW': 100},  # Production
    ]

    for i, stage in enumerate(aimd_stages):
        print(f"   Stage {i}: {stage['TEBEG']} K, {stage['NSW']} steps")

    # Base builder inputs
    builder_inputs = {
        'parameters': {
            'incar': {
                # Basic settings (BASE)
                'PREC': 'Normal',
                'ENCUT': 400,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': 'Auto',

                # AIMD settings
                'IBRION': 0,
                'MDALGO': 2,
                'POTIM': 2.0,
                'SMASS': 0.0,

                # Output control
                'LWAVE': False,
                'LCHARG': False,
            }
        },
        'kpoints_spacing': 0.5,
        'potential_family': potential_family,
        'potential_mapping': {'Ag': 'Ag'},
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'clean_workdir': False,
    }

    print("\n5. Override configuration:")
    print("   Base: ENCUT=400, PREC=Normal, ALGO=Normal")

    # Structure override: ag2 uses higher cutoff
    structure_overrides = {
        'ag2': {'parameters': {'incar': {'ENCUT': 500}}}
    }
    print("   Structure override: ag2 uses ENCUT=500")

    # Stage override: stage 1 (production) uses Accurate
    stage_overrides = {
        1: {'parameters': {'incar': {'PREC': 'Accurate', 'EDIFF': 1e-6}}}
    }
    print("   Stage override: stage 1 uses PREC=Accurate, EDIFF=1e-6")

    # Matrix override: (ag1, stage 1) uses Fast algorithm
    matrix_overrides = {
        ('ag1', 1): {'parameters': {'incar': {'ALGO': 'Fast'}}}
    }
    print("   Matrix override: (ag1, stage 1) uses ALGO=Fast")

    print("\n6. Expected INCAR values:")
    print("   ag1, stage 0: ENCUT=400, PREC=Normal, ALGO=Normal")
    print("   ag1, stage 1: ENCUT=400, PREC=Accurate, ALGO=Fast, EDIFF=1e-6")
    print("   ag2, stage 0: ENCUT=500, PREC=Normal, ALGO=Normal")
    print("   ag2, stage 1: ENCUT=500, PREC=Accurate, ALGO=Normal, EDIFF=1e-6")

    print("\n7. Building workgraph with overrides...")

    # Build AIMD workgraph
    wg = build_aimd_workgraph(
        # Input structures
        structures={
            'ag1': ag_structure1,
            'ag2': ag_structure2,
        },

        # AIMD stages
        aimd_stages=aimd_stages,

        # Code configuration
        code_label=code_label,

        # Base builder inputs
        builder_inputs=builder_inputs,

        # Optional: Create supercells before AIMD
        # Example: supercell_specs={'ag1': [2, 2, 1], 'ag2': [3, 3, 1]}
        supercell_specs={'ag1': [2, 2, 1], 'ag2': [3, 3, 1]},

        # Override system (now functional!)
        structure_overrides=structure_overrides,
        stage_overrides=stage_overrides,
        matrix_overrides=matrix_overrides,

        # Concurrency control
        max_concurrent_jobs=1,

        # Workgraph name
        name='Step19_AIMD_Overrides_Ag',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n8. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 19 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process show {wg.pk}")
    print(f"\nVerify INCAR overrides:")
    print(f"  # Find VASP calculation PKs")
    print(f"  verdi process show {wg.pk}  # Look for VaspWorkChain PKs")
    print(f"  # Check INCAR for each calculation")
    print(f"  verdi calcjob inputcat <VASP_PK> INCAR | grep -E 'ENCUT|PREC|ALGO|EDIFF'")
    print(f"\nExpected to see different INCAR values per structure/stage!")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
