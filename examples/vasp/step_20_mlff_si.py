#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 20: MLFF Training with Silicon

Simple MLFF training workflow based on VASP Liquid Si example.
Trains ML force field on-the-fly during MD simulation.

Workflow:
  Stage 0: Training (ML_ISTART=0)
    - VASP runs DFT+MD
    - ML model learns on-the-fly
    - Generates ML_AB (training data) and ML_FFN (model)

Optional Stage 1: Production (ML_ISTART=2)
    - Uses trained ML model
    - Much faster than DFT
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.mlff import build_mlff_workgraph
from ase.io import read


def main():
    print("\n" + "="*70)
    print("STEP 20: MLFF Training with Silicon")
    print("="*70)

    # Load profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Load structure
    script_dir = os.path.dirname(os.path.abspath(__file__))
    si_cif = os.path.join(script_dir, 'structures', 'Si.cif')

    print(f"\n2. Loading structure: {si_cif}")
    si_ase = read(si_cif)
    si_structure = orm.StructureData(ase=si_ase)
    print(f"   ✓ {len(si_structure.sites)} atoms")

    print(f"\n3. MLFF configuration:")
    print(f"   Training: 100 MD steps with ML learning")
    print(f"   Temperature: 300 K")
    print(f"   Code: VASP6.5.0@cluster02")

    # Build MLFF workgraph (training only)
    print("\n4. Building workgraph...")
    wg = build_mlff_workgraph(
        structures={'si_bulk': si_structure},
        training_steps=100,         # Train for 100 steps
        temperature=300,            # 300 K
        code_label='VASP6.5.0@cluster02',
        potential_family='PBE',
        potential_mapping={'Si': 'Si'},
        # Optional: Add production stage
        # production_steps=500,     # Uncomment to add ML production run
    )
    print("   ✓ WorkGraph created")

    print("\n5. Workflow stages:")
    print("   Stage 0: Training (ML_ISTART=0)")
    print("     - DFT+MD with on-the-fly ML learning")
    print("     - Generates ML_AB and ML_FFN files")
    print("     - ~50 minutes estimated")

    # Submit
    print("\n6. Submitting...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("SUBMITTED")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nAfter completion, check for ML files:")
    print(f"  - ML_ABN: Training data")
    print(f"  - ML_FFN: Trained neural network")
    print(f"  - ML_LOGFILE: Training log with errors")
    print(f"\nTo add production stage:")
    print(f"  Uncomment 'production_steps=500' in the script")
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
