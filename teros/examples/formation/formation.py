#!/usr/bin/env python
"""
Example script for calculating formation enthalpy of any ternary oxide using PS-TEROS.

This script demonstrates how to:
1. Load the AiiDA profile
2. Set up separate parameters for bulk, metal, nonmetal, and oxygen calculations
3. Create and run a formation enthalpy WorkGraph
4. The workflow will relax the bulk and all references in parallel

Example: Ag3PO4
    - Bulk: Ag3PO4
    - Metal reference: Ag
    - Nonmetal reference: P
    - Oxygen reference: O2

Usage:
    source ~/envs/aiida/bin/activate && python formation.py
"""

from aiida import load_profile, orm
from teros.workgraph import build_formation_workgraph

def main():
    """Main function to run the formation enthalpy workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile()

    # Define structures directory
    structures_dir = '/home/thiagotd/git/PS-TEROS/teros/structures'

    # Define calculation parameters
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # VASP parameters for bulk relaxation (Ag3PO4)
    bulk_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,  # Conjugate gradient
        'ISIF': 3,    # Relax cell shape, volume, and atoms
        'NSW': 100,   # Max ionic steps
        'EDIFFG': -0.1,  # Force convergence
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options for bulk
    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # VASP parameters for metal reference (Ag)
    metal_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 1,  # Methfessel-Paxton for metals
        'SIGMA': 0.2,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.1,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options for metal reference
    metal_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # VASP parameters for nonmetal reference (P)
    nonmetal_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,  # Gaussian smearing
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.1,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options for nonmetal reference
    nonmetal_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # VASP parameters for oxygen reference (O2)
    # Molecular oxygen requires special treatment
    oxygen_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,  # Gaussian smearing for molecules
        'SIGMA': 0.01,  # Small smearing for molecules
        'IBRION': 2,
        'ISIF': 2,    # Relax ions only, not cell (molecule in box)
        'NSW': 100,
        'EDIFFG': -0.1,
        'ALGO': 'Normal',
        'LREAL': False,  # Use reciprocal space for small systems
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options for oxygen reference
    oxygen_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    print(f"\nSetting up formation enthalpy calculation for Ag3PO4")
    print(f"Structures directory: {structures_dir}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")
    print(f"\nWill relax:")
    print(f"  - Bulk: ag3po4.cif")
    print(f"  - Metal reference: Ag.cif (with metal-specific parameters)")
    print(f"  - Nonmetal reference: P.cif (with nonmetal-specific parameters)")
    print(f"  - Oxygen reference: O2.cif (with molecule-specific parameters)")

    # Create the WorkGraph
    print("\nCreating WorkGraph...")
    wg = build_formation_workgraph(
        structures_dir=structures_dir,
        bulk_name='ag3po4.cif',
        metal_name='Ag.cif',
        nonmetal_name='P.cif',
        oxygen_name='O2.cif',
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        metal_potential_mapping={'Ag': 'Ag'},
        nonmetal_potential_mapping={'P': 'P'},
        oxygen_potential_mapping={'O': 'O'},
        kpoints_spacing=0.3,  # A^-1 * 2pi
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=True,
        name='Ag3PO4_Formation',
    )

    # Optional: Export WorkGraph to HTML for visualization
    try:
        html_file = 'ag3po4_formation_workgraph.html'
        wg.to_html(html_file)
        print(f"WorkGraph visualization saved to: {html_file}")
    except Exception as e:
        print(f"Could not generate HTML visualization: {e}")

    # Submit the WorkGraph
    print("\nSubmitting WorkGraph...")
    wg.submit(wait=False)

    print(f"\nWorkGraph submitted successfully!")
    print(f"WorkGraph PK: {wg.pk}")
    print(f"\nTo check status, run:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")
    print(f"\nTo monitor the workflow:")
    print(f"  verdi process list")
    print(f"\nOutputs available after completion:")
    print(f"  - bulk_energy: Energy of relaxed Ag3PO4")
    print(f"  - metal_energy: Energy of relaxed Ag")
    print(f"  - nonmetal_energy: Energy of relaxed P")
    print(f"  - oxygen_energy: Energy of relaxed O2")
    print(f"  - bulk_structure, metal_structure, nonmetal_structure, oxygen_structure")

    return wg


if __name__ == '__main__':
    """
    Run the formation enthalpy workflow.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    5. Run this script:
       source ~/envs/aiida/bin/activate && python formation.py
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
