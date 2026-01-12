#!/home/thiagotd/envs/aiida/bin/python
"""
Test DOS calculation module - Single or multiple structures.

This example demonstrates running DOS calculations using the
build_dos_calculation_workgraph function from the custom_calculation module.

The workflow:
1. Performs SCF calculation to get charge density (CHGCAR)
2. Performs non-SCF DOS calculation using the CHGCAR

Usage:
    source ~/envs/aiida/bin/activate
    python test_dos.py
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.custom_calculation import build_dos_calculation_workgraph, get_dos_results


def main():
    """Run DOS calculation example."""

    print("\n" + "="*70)
    print("DOS CALCULATION - Example Test")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   Daemon is running")

    # Load structure - modify this for your use case
    print("\n2. Loading structure...")
    # Option 1: Load from file
    structure_file = Path("test_structure.vasp")
    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        structure = orm.StructureData(ase=atoms)
        print(f"   Loaded from file: {structure_file}")
    else:
        # Option 2: Load from database (modify PK)
        structure_pk = 12345  # MODIFY THIS
        try:
            structure = orm.load_node(structure_pk)
            print(f"   Loaded from PK: {structure_pk}")
        except Exception:
            print(f"   ERROR: Could not load structure from PK {structure_pk}")
            print("   Please provide a valid structure file or PK")
            return 1

    print(f"   Composition: {structure.get_composition()}")

    # Define SCF inputs
    print("\n3. Defining SCF calculation inputs...")
    scf_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-5,
                'ALGO': 'Fast',
                'NELM': 120,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'LWAVE': True,
                'LCHARG': True,
                'LORBIT': 11,  # For projected DOS
                'NCORE': 6,
                'LREAL': 'Auto',
                'ISPIN': 1,  # Non-spin-polarized (modify for magnetic systems)
            }
        },
        'kpoints_spacing': 0.3,  # SCF k-points spacing
        'settings': {
            'ADDITIONAL_RETRIEVE_LIST': ['vasprun.xml', 'CHGCAR', 'OUTCAR'],
        },
    }

    # Define DOS inputs
    print("\n4. Defining DOS calculation inputs...")
    dos_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-5,
                'ALGO': 'Normal',
                'NELM': 60,
                'ISTART': 1,   # Read WAVECAR
                'ICHARG': 11,  # Non-SCF from CHGCAR
                'ISMEAR': -5,  # Tetrahedron method for DOS
                'LREAL': False,  # Reciprocal space for accurate DOS
                'LWAVE': False,
                'LCHARG': False,
                'LORBIT': 11,  # Projected DOS
                'NEDOS': 2001,  # Number of DOS points
                'NCORE': 6,
                'ISPIN': 1,
            }
        },
        'settings': {
            'ADDITIONAL_RETRIEVE_LIST': ['vasprun.xml', 'DOSCAR', 'OUTCAR'],
        },
    }

    # VASP code and potentials
    code_label = 'VASP-6.4.1@cluster02'  # MODIFY THIS
    potential_family = 'PBE.54'
    potential_mapping = {'Ag': 'Ag', 'O': 'O', 'P': 'P'}  # MODIFY for your structure

    # Calculation options
    options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    }

    # DOS k-points spacing (denser for better resolution)
    dos_kpoints_distance = 0.2

    print(f"\n   VASP code: {code_label}")
    print(f"   SCF k-points spacing: {scf_inputs['kpoints_spacing']}")
    print(f"   DOS k-points spacing: {dos_kpoints_distance}")

    # Build WorkGraph
    print("\n5. Building DOS WorkGraph...")
    wg = build_dos_calculation_workgraph(
        structure=structure,
        code_label=code_label,
        scf_inputs=scf_inputs,
        dos_inputs=dos_inputs,
        dos_kpoints_distance=dos_kpoints_distance,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=options,
        name='test_dos_calculation',
    )

    print(f"   WorkGraph created: {wg.name}")
    print(f"   Tasks: {list(wg.tasks.keys())}")

    # Submit WorkGraph
    print("\n6. Submitting WorkGraph...")
    wg.submit(wait=False)
    print(f"   Submitted! PK: {wg.pk}")
    print(f"\n   Monitor with:")
    print(f"     verdi process show {wg.pk}")
    print(f"     verdi process report {wg.pk}")

    # Note: To get results, run with wait=True or check later
    print(f"\n   To get results after completion:")
    print(f"     from aiida import orm")
    print(f"     from teros.core.custom_calculation import get_dos_results")
    print(f"     wg = orm.load_node({wg.pk})")
    print(f"     results = get_dos_results(wg)")
    print(f"     print(results['dos'])")
    print(f"     print(results['projectors'])  # For projected DOS")

    print("\n" + "="*70)
    print("DOS calculation submitted!")
    print("="*70 + "\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
