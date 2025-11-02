#!/home/thiagotd/envs/aiida/bin/python
"""
Test custom calculation module - Multiple structures with different settings.

This example demonstrates running custom VASP calculations on multiple structures
with different builder inputs for each structure.
"""

import sys
from aiida import orm, load_profile
from teros.core.custom_calculation import build_custom_calculation_workgraph, get_custom_results

def main():
    """Run custom VASP calculations with different settings per structure."""

    print("\n" + "="*70)
    print("CUSTOM VASP CALCULATION - Multiple Structures (Different Settings)")
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

    # Load structures
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

    # Define different builder inputs for each structure
    print("\n3. Defining different VASP builder inputs for each structure...")

    # Structure 1: Relaxation with ALGO=Fast
    builder_1 = {
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

    # Structure 2: Static calculation (NSW=0)
    builder_2 = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-6,  # Tighter convergence
                'ISMEAR': -5,    # Tetrahedron method
                'ALGO': 'Normal',
                'LREAL': False,
                'LWAVE': True,
                'LCHARG': True,
                'NCORE': 6,
                'NSW': 0,        # Static calculation
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'kpoints_spacing': 0.25,  # Denser k-points
        'potential_family': 'PBE.54',
        'potential_mapping': {'Ag': 'Ag', 'O': 'O', 'P': 'P'},
        'clean_workdir': False,  # Keep files for post-processing
    }

    # Structure 3: High-precision relaxation
    builder_3 = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 600,    # Higher cutoff
                'EDIFF': 1e-6,
                'ISMEAR': 0,
                'SIGMA': 0.02,   # Smaller smearing
                'ALGO': 'Normal',
                'LREAL': False,  # No real-space projection
                'LWAVE': False,
                'LCHARG': False,
                'NCORE': 6,
                'ISIF': 2,
                'NSW': 200,      # More steps
                'IBRION': 2,
                'EDIFFG': -0.01, # Tighter forces
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'kpoints_spacing': 0.25,
        'potential_family': 'PBE.54',
        'potential_mapping': {'Ag': 'Ag', 'O': 'O', 'P': 'P'},
        'clean_workdir': True,
    }

    builder_inputs_list = [builder_1, builder_2, builder_3]

    print("   Structure 0: Relaxation with ALGO=Fast")
    print("   Structure 1: Static calculation (NSW=0)")
    print("   Structure 2: High-precision relaxation")

    # Build WorkGraph
    print("\n4. Building WorkGraph...")
    wg = build_custom_calculation_workgraph(
        structure=structures,
        code_label='VASP-6.4.1@cluster02',
        builder_inputs=builder_inputs_list,  # List of dicts
        name='test_different_settings'
    )

    print(f"   ✓ WorkGraph created: {wg.name}")
    print(f"   Tasks: {len(wg.tasks)}")

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
