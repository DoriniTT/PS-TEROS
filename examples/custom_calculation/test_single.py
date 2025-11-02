#!/home/thiagotd/envs/aiida/bin/python
"""
Test custom calculation module - Single structure.

This example demonstrates running a custom VASP calculation on a single structure
with full control over builder inputs.
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.custom_calculation import build_custom_calculation_workgraph, get_custom_results

def main():
    """Run custom VASP calculation on single structure."""

    print("\n" + "="*70)
    print("CUSTOM VASP CALCULATION - Single Structure Test")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   ✓ Daemon is running")

    # Load structure - you can modify this to use your own structure
    print("\n2. Loading structure...")
    # Option 1: Load from file
    structure_file = Path("test_structure.vasp")
    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        structure = orm.StructureData(ase=atoms)
        print(f"   ✓ Loaded from file: {structure_file}")
    else:
        # Option 2: Load from database (modify PK)
        structure_pk = 12345  # MODIFY THIS
        try:
            structure = orm.load_node(structure_pk)
            print(f"   ✓ Loaded from PK: {structure_pk}")
        except:
            print(f"   ERROR: Could not load structure from PK {structure_pk}")
            print("   Please provide a valid structure file or PK")
            return 1

    print(f"   Composition: {structure.get_composition()}")

    # Define builder inputs (full control)
    print("\n3. Defining VASP builder inputs...")
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
        'potential_mapping': {'Ag': 'Ag', 'O': 'O', 'P': 'P'},  # Modify for your structure
        'clean_workdir': True,
    }

    print("   INCAR settings:")
    for key, val in builder_inputs['parameters']['incar'].items():
        print(f"     {key}: {val}")

    # Build WorkGraph
    print("\n4. Building WorkGraph...")
    code_label = 'VASP-6.4.1@cluster02'

    wg = build_custom_calculation_workgraph(
        structure=structure,
        code_label=code_label,
        builder_inputs=builder_inputs,
        name='test_single_custom_calc'
    )

    print(f"   ✓ WorkGraph created: {wg.name}")
    print(f"   Tasks: {list(wg.tasks.keys())}")

    # Submit WorkGraph
    print("\n5. Submitting WorkGraph...")
    wg.submit(wait=False)
    print(f"   ✓ Submitted! PK: {wg.pk}")
    print(f"\n   Monitor with:")
    print(f"     verdi process show {wg.pk}")
    print(f"     verdi process report {wg.pk}")

    # Note: To get results, run with wait=True or check later
    print(f"\n   To get results after completion:")
    print(f"     from aiida import orm")
    print(f"     from teros.core.custom_calculation import get_custom_results")
    print(f"     wg = orm.load_node({wg.pk})")
    print(f"     results = get_custom_results(wg)")
    print(f"     print(results['energies'])")

    print("\n" + "="*70)
    print("Test complete!")
    print("="*70 + "\n")

    return 0

if __name__ == '__main__':
    sys.exit(main())
