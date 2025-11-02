#!/home/thiagotd/envs/aiida/bin/python
"""
Custom VASP Calculation - Two Structures with Different Settings

This script demonstrates how to run custom VASP calculations on 2 structures
with different builder inputs (INCAR settings) for each structure using the
PS-TEROS custom calculation module.

This example shows:
- Loading structures from files or by PK
- Defining separate builder_inputs for each structure
- Building and submitting a custom calculation WorkGraph
- Monitoring workflow progress
- Retrieving results after completion

Material: Ag2O and Ag3PO4
Purpose: Compare relaxation strategies (fast vs. high-precision)

Usage:
    source ~/envs/aiida/bin/activate
    python test_custom_calculation_two_structures.py
"""

import sys
import os
from aiida import orm, load_profile
from teros.core.custom_calculation import build_custom_calculation_workgraph, get_custom_results


def load_structures_from_files(structures_dir):
    """
    Load structures from CIF files.

    Args:
        structures_dir: Path to directory containing CIF files

    Returns:
        List of StructureData objects
    """
    from ase.io import read

    # Define structure files
    structure_files = [
        os.path.join(structures_dir, 'ag2o.cif'),
        os.path.join(structures_dir, 'ag3po4.cif'),
    ]

    structures = []
    for filepath in structure_files:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Structure file not found: {filepath}")

        # Read with ASE and convert to StructureData
        ase_atoms = read(filepath)
        structure = orm.StructureData(ase=ase_atoms)
        structures.append(structure)

        print(f"   ✓ Loaded: {os.path.basename(filepath)} ({structure.get_composition()})")

    return structures


def load_structures_from_pks(structure_pks):
    """
    Load structures from AiiDA database by PK.

    Args:
        structure_pks: List of PKs (integers)

    Returns:
        List of StructureData objects
    """
    structures = []
    for pk in structure_pks:
        try:
            struct = orm.load_node(pk)
            if not isinstance(struct, orm.StructureData):
                raise TypeError(f"PK {pk} is not a StructureData node")
            structures.append(struct)
            print(f"   ✓ Loaded PK {pk}: {struct.get_composition()}")
        except Exception as e:
            raise ValueError(f"Could not load structure PK {pk}: {e}")

    return structures


def main():
    """Run custom VASP calculations on 2 structures with different settings."""

    print("\n" + "="*70)
    print("CUSTOM VASP CALCULATION - Two Structures (Different Settings)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded: psteros")

    # Check daemon
    print("\n2. Checking AiiDA daemon...")
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("   ✗ ERROR: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   ✓ Daemon is running")

    # Load structures
    print("\n3. Loading structures...")
    print("   Choose one method:")
    print("   A) Load from CIF files")
    print("   B) Load from AiiDA database by PK")

    # Method A: Load from files (default)
    use_files = True  # SET TO False TO USE PKs INSTEAD

    if use_files:
        print("\n   Using method A: Loading from CIF files...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        structures_dir = os.path.join(script_dir, '../structures')
        structures = load_structures_from_files(structures_dir)
    else:
        print("\n   Using method B: Loading from PKs...")
        # MODIFY THESE PKs TO YOUR STRUCTURES
        structure_pks = [12345, 12346]
        structures = load_structures_from_pks(structure_pks)

    print(f"   Total structures loaded: {len(structures)}")

    # Define builder inputs for each structure
    print("\n4. Defining VASP builder inputs for each structure...")

    # Code configuration (same for both)
    code_label = 'VASP-6.4.1@cluster02'
    potential_family = 'PBE.54'

    # Structure 1: Fast relaxation (Ag2O)
    print("   Structure 0 (Ag2O): Fast relaxation with ALGO=Fast")
    builder_1 = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,           # Standard cutoff
                'EDIFF': 1e-5,          # Standard convergence
                'ISMEAR': 0,            # Gaussian smearing
                'SIGMA': 0.05,
                'ALGO': 'Fast',         # Fast algorithm
                'LREAL': 'Auto',        # Real-space projection
                'LWAVE': False,
                'LCHARG': False,
                'NCORE': 6,
                'ISIF': 2,              # Relax ions only
                'NSW': 100,             # Max ionic steps
                'IBRION': 2,            # CG relaxation
                'EDIFFG': -0.02,        # Force convergence
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'kpoints_spacing': 0.3,
        'potential_family': potential_family,
        'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
        'clean_workdir': True,
    }

    # Structure 2: High-precision relaxation (Ag3PO4)
    print("   Structure 1 (Ag3PO4): High-precision with tighter convergence")
    builder_2 = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 600,           # Higher cutoff
                'EDIFF': 1e-6,          # Tighter convergence
                'ISMEAR': 0,
                'SIGMA': 0.02,          # Smaller smearing
                'ALGO': 'Normal',       # More stable algorithm
                'LREAL': False,         # Reciprocal space (more accurate)
                'LWAVE': False,
                'LCHARG': False,
                'NCORE': 6,
                'ISIF': 2,
                'NSW': 200,             # More ionic steps
                'IBRION': 2,
                'EDIFFG': -0.01,        # Tighter forces
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'kpoints_spacing': 0.25,        # Denser k-points
        'potential_family': potential_family,
        'potential_mapping': {'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        'clean_workdir': True,
    }

    # Create list of builder inputs (one per structure)
    builder_inputs_list = [builder_1, builder_2]

    # Build WorkGraph
    print("\n5. Building custom calculation WorkGraph...")
    wg = build_custom_calculation_workgraph(
        structure=structures,                  # List of 2 StructureData
        code_label=code_label,
        builder_inputs=builder_inputs_list,    # List of 2 builder dicts
        name='custom_calc_two_structures'
    )

    print(f"   ✓ WorkGraph created: {wg.name}")
    print(f"   Tasks: {len(wg.tasks)}")
    print(f"   - vasp_calc_0: Ag2O fast relaxation")
    print(f"   - vasp_calc_1: Ag3PO4 high-precision relaxation")

    # Submit
    print("\n6. Submitting WorkGraph...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")

    print(f"\n--- MONITORING COMMANDS ---")
    print(f"Check status:")
    print(f"  verdi process status {wg.pk}")
    print(f"\nView details:")
    print(f"  verdi process show {wg.pk}")
    print(f"\nView report:")
    print(f"  verdi process report {wg.pk}")
    print(f"\nWatch live (updates every 10s):")
    print(f"  watch -n 10 'verdi process status {wg.pk}'")

    print(f"\n--- EXPECTED OUTPUTS ---")
    print(f"After completion, the workflow will have:")
    print(f"  - energies: [E_Ag2O, E_Ag3PO4] (in eV)")
    print(f"  - structures: [relaxed_Ag2O, relaxed_Ag3PO4]")
    print(f"  - misc: [VASP outputs for each calculation]")

    print(f"\n--- RETRIEVE RESULTS (after completion) ---")
    print(f"In Python:")
    print(f"  from aiida import orm, load_profile")
    print(f"  from teros.core.custom_calculation import get_custom_results")
    print(f"  load_profile('psteros')")
    print(f"  wg = orm.load_node({wg.pk})")
    print(f"  results = get_custom_results(wg)")
    print(f"  print('Energies:', results['energies'])")
    print(f"  print('Structures:', results['structures'])")

    print(f"\n{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
