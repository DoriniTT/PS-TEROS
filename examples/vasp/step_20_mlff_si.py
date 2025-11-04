#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 20: MLFF with Silicon Bulk

Tests basic MLFF workflow: Training → Production
Uses simple Si bulk structure for fast testing.

Workflow:
- Stage 0: Train MLFF (20 steps for quick test)
- Stage 1: Production run with ML forces (50 steps)
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.mlff import build_mlff_workgraph
from ase.io import read


def main():
    print("\n" + "="*70)
    print("STEP 20: MLFF with Silicon Bulk")
    print("="*70)

    # Load profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Load structure
    script_dir = os.path.dirname(os.path.abspath(__file__))
    si_cif = os.path.join(script_dir, 'structures', 'Si.cif')

    print(f"\n2. Loading Si structure: {si_cif}")
    si_ase = read(si_cif)
    si_structure = orm.StructureData(ase=si_ase)
    print(f"   ✓ Loaded: {len(si_structure.sites)} atoms")

    # Code configuration
    code_label = 'VASP6.5.0@cluster02'
    potential_family = 'PBE'

    print(f"\n3. VASP configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential: {potential_family}")

    # Builder inputs
    builder_inputs = {
        'parameters': {
            'incar': {
                # DFT settings
                'PREC': 'Normal',
                'ENCUT': 400,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': 'Auto',

                # MD settings
                'IBRION': 0,      # MD
                'MDALGO': 2,      # Nosé-Hoover
                'POTIM': 1.0,     # 1 fs timestep
                'SMASS': 0.0,

                # Output
                'LWAVE': True,    # Required for restart
                'LCHARG': True,
            }
        },
        'kpoints_spacing': 0.5,
        'potential_family': potential_family,
        'potential_mapping': {'Si': 'Si'},
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'clean_workdir': False,  # Keep ML files!
    }

    print("\n4. MLFF configuration:")
    print("   Training: 20 steps (DFT+ML)")
    print("   Production: 50 steps (pure ML)")
    print("   Temperature: 300 K")

    # Build MLFF workgraph
    print("\n5. Building MLFF workgraph...")
    wg = build_mlff_workgraph(
        structures={'si_bulk': si_structure},
        training_steps=20,      # Short for testing
        production_steps=50,
        temperature=300,
        code_label=code_label,
        builder_inputs=builder_inputs,
        name='MLFF_Si_Test',
    )
    print("   ✓ WorkGraph built")

    print("\n6. Expected workflow:")
    print("   Stage 0: Training (20 steps, ~10 min)")
    print("     → Generates ML_AB, ML_FFN")
    print("   Stage 1: Production (50 steps, ~2 min)")
    print("     → Uses ML forces (should be faster)")

    # Submit
    print("\n7. Submitting...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nCritical checks:")
    print(f"  1. Stage 0 completes (check for ML_ABN, ML_FFN)")
    print(f"  2. Stage 1 finds ML files (check OUTCAR for 'reading ML_FFN')")
    print(f"  3. Stage 1 faster than Stage 0")
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
