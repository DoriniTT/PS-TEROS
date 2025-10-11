#!/home/thiagotd/envs/aiida/bin/python
"""
Example script for running bulk relaxation of Ag3PO4 using PS-TEROS WorkGraph.

This script demonstrates how to:
1. Load the AiiDA profile
2. Set up the structure and calculation parameters
3. Create and run a bulk relaxation WorkGraph
4. Retrieve and analyze results

Usage:
    source ~/envs/aiida/bin/activate && python relaxation.py
"""

import sys
from aiida import load_profile, orm
from teros.core.workgraph import build_relaxation_workgraph


def main():
    """Main function to run the bulk relaxation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define structure path
    structure_path = '/home/thiagotd/git/PS-TEROS/teros/structures/ag3po4.cif'

    # Define calculation parameters
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # VASP parameters for relaxation
    parameters = {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,  # Conjugate gradient
        'ISIF': 3,    # Relax cell shape, volume, and atoms
        'NSW': 100,   # Max ionic steps
        'EDIFFG': -0.01,  # Force convergence
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        }
    }
    # Scheduler options
    options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Potential mapping (element -> potential)
    # For Ag3PO4: Ag, P, O
    potential_mapping = {
        'Ag': 'Ag',
        'P': 'P',
        'O': 'O',
    }

    print(f"\nSetting up bulk relaxation for: {structure_path}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")

    # Create the WorkGraph
    print("\nCreating WorkGraph...")
    wg = build_relaxation_workgraph(
        structure_path=structure_path,
        code_label=code_label,
        kpoints_spacing=0.3,  # A^-1 * 2pi
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        clean_workdir=True,
        name='Ag3PO4_BulkRelaxation',
    )

    # Optional: Export WorkGraph to HTML for visualization
    try:
        html_file = 'ag3po4_relaxation_workgraph.html'
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
    print(f"\nAfter completion, analyze results with:")
    print(f"  verdi process show {wg.pk}")

    return wg


if __name__ == '__main__':
    """
    Run the bulk relaxation workflow.

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
       source ~/envs/aiida/bin/activate && python relaxation.py
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
