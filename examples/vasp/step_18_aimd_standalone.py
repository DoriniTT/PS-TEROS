#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 18: Standalone AIMD Module

This script demonstrates the standalone AIMD module (teros.core.aimd) which provides
full control over multi-stage AIMD calculations on pre-existing structures.

Features demonstrated:
- Direct structure input (no bulk workflow coupling)
- Sequential multi-stage AIMD with automatic restart chaining
- Optional supercell transformation
- Concurrency control with max_concurrent_jobs
- Simple API focused on AIMD only

Material: Ag (silver bulk)
AIMD: 2-stage sequence (equilibration + production)
Structures: 2 identical Ag structures to demonstrate parallel execution

NOTE: Override parameters (structure_overrides, stage_overrides, matrix_overrides)
      are NOT functional in this version. See module documentation for details.

Usage:
    source ~/envs/aiida/bin/activate
    python step_18_aimd_standalone.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.aimd import build_aimd_workgraph
from ase.io import read


def main():
    """Step 18: Test standalone AIMD module."""

    print("\n" + "="*70)
    print("STEP 18: STANDALONE AIMD MODULE")
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

    # Code configuration (use the code with 24 processors from cluster02)
    code_label = 'VASP6.5.0@cluster02'  # No dash between VASP and version
    potential_family = 'PBE'

    print(f"\n3. VASP configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # AIMD stages
    print("\n4. AIMD configuration:")
    aimd_stages = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    for i, stage in enumerate(aimd_stages):
        print(f"   Stage {i}: {stage['temperature']} K, {stage['steps']} steps")

    # Builder inputs (common for all structures and stages)
    builder_inputs = {
        'parameters': {
            'incar': {
                # Basic settings
                'PREC': 'Normal',
                'ENCUT': 400,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': 'Auto',

                # AIMD settings (temperature and steps will be set by module)
                'IBRION': 0,      # Molecular dynamics
                'MDALGO': 2,      # Nosé-Hoover thermostat
                'POTIM': 2.0,     # Time step in fs
                'SMASS': 0.0,     # Nosé mass

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

    print("\n5. Building standalone AIMD workgraph...")
    print("   Features:")
    print("   - 2 structures (ag1, ag2) running in parallel")
    print("   - 2 sequential stages with restart chaining")
    print("   - Supercell transformation on ag1 only (2x2x1)")
    print("   - Max 2 concurrent VASP jobs")

    # Build AIMD workgraph
    wg = build_aimd_workgraph(
        # Input structures
        structures={
            'ag1': ag_structure1,
            'ag2': ag_structure2,
        },

        # AIMD stages (sequential)
        aimd_stages=aimd_stages,

        # Code configuration
        code_label=code_label,

        # Builder inputs (same for all)
        builder_inputs=builder_inputs,

        # Optional: supercell for ag1 only
        supercell_specs={
            'ag1': [2, 2, 1],  # 2x2x1 supercell
        },

        # Concurrency control
        max_concurrent_jobs=2,  # Limit to 2 concurrent VASP calculations

        # Workgraph name
        name='Step18_AIMD_Standalone_Ag',
    )

    print("   ✓ WorkGraph built successfully")

    # Display workgraph structure
    print("\n6. WorkGraph structure:")
    print("   Expected task flow:")
    print("   1. create_supercell_ag1  (supercell for ag1)")
    print("   2. stage_0_aimd          (equilibration on both structures)")
    print("   3. stage_1_aimd          (production, restarts from stage 0)")

    # Submit
    print("\n7. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 18 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - results: Nested dict with AIMD data per structure per stage")
    print(f"    - results[stage_idx]['structures']: Final structures from stage")
    print(f"    - results[stage_idx]['energies']: Total energies from stage")
    print(f"    - results[stage_idx]['remote_folders']: Remote folders for restart")
    print(f"  - supercells: Dict with supercell structures")
    print(f"    - supercells['ag1']: 2x2x1 supercell of original ag1")
    print(f"\nNote about overrides:")
    print(f"  The override parameters (structure_overrides, stage_overrides,")
    print(f"  matrix_overrides) are NOT functional in this version.")
    print(f"  See teros/core/aimd/README.md for implementation status.")
    print(f"\nAIMD trajectories can be analyzed for:")
    print(f"  - Temperature evolution")
    print(f"  - Atomic diffusion")
    print(f"  - Structural stability")
    print(f"  - Phase transitions")
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
