#!/home/thiagotd/envs/aiida/bin/python
"""
Test custom calculation module - Multiple structures with same settings.

This example demonstrates running custom VASP calculations on multiple structures
using the same builder inputs for all.
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.custom_calculation import build_custom_calculation_workgraph, get_custom_results

def main():
    """Run custom VASP calculations on multiple structures."""

    print("\n" + "="*70)
    print("CUSTOM VASP CALCULATION - Multiple Structures (Same Settings)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Check daemon
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        return 1
    print("   ✓ Daemon is running")

    # Load multiple structures
    print("\n2. Loading structures...")
    structure_pks = [12345, 12346, 12347]  # MODIFY THESE

    structures = []
    for pk in structure_pks:
        try:
            struct = orm.load_node(pk)
            structures.append(struct)
            print(f"   ✓ Loaded PK {pk}: {struct.get_composition()}")
        except:
            print(f"   ERROR: Could not load structure PK {pk}")
            return 1

    print(f"   Total structures: {len(structures)}")

    # Define builder inputs (same for all structures)
    print("\n3. Defining VASP builder inputs (same for all)...")
    builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
                'NCORE': 6,
                'ISIF': 2,
                'NSW': 100,
                'IBRION': 2,
                'EDIFFG': -0.02,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'kpoints_spacing': 0.3,
        'potential_family': 'PBE.54',
        'potential_mapping': {'Ag': 'Ag', 'O': 'O', 'P': 'P'},
        'clean_workdir': True,
    }

    # Build WorkGraph
    print("\n4. Building WorkGraph...")
    wg = build_custom_calculation_workgraph(
        structure=structures,
        code_label='VASP-6.4.1@cluster02',
        builder_inputs=builder_inputs,  # Single dict for all
        name='test_multiple_custom_calc'
    )

    print(f"   ✓ WorkGraph created: {wg.name}")
    print(f"   Tasks: {len(wg.tasks)} ({list(wg.tasks.keys())})")

    # Submit
    print("\n5. Submitting WorkGraph...")
    wg.submit(wait=False)
    print(f"   ✓ Submitted! PK: {wg.pk}")
    print(f"\n   Monitor with: verdi process show {wg.pk}")

    print("\n" + "="*70)
    print("Test complete!")
    print("="*70 + "\n")

    return 0

if __name__ == '__main__':
    sys.exit(main())
