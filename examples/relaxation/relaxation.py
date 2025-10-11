#!/home/thiagotd/envs/aiida/bin/python
"""
Example script for running bulk relaxation of Ag3PO4 using PS-TEROS WorkGraph.

This script demonstrates how to:
1. Load the AiiDA profile
2. Set up the structure and calculation parameters
3. Create and run a bulk relaxation WorkGraph
4. Retrieve and analyze results

This example performs bulk relaxation only (no reference structures).
The compute_relaxation_energy, compute_cleavage, and compute_thermodynamics
flags control what gets calculated. Since no metal_name/oxygen_name are provided,
thermodynamics will be skipped even though the flag defaults to True.

Usage:
    source ~/envs/aiida/bin/activate && python relaxation.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def main():
    """Main function to run the bulk relaxation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define structure path (relative to this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    bulk_filename = 'ag3po4.cif'

    # Define calculation parameters
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # VASP parameters for relaxation (INCAR tags)
    bulk_parameters = {
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

    # Scheduler options
    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Potential mapping (element -> potential)
    # For Ag3PO4: Ag, P, O
    bulk_potential_mapping = {
        'Ag': 'Ag',
        'P': 'P',
        'O': 'O',
    }

    print(f"\nSetting up bulk relaxation for: {os.path.join(structures_dir, bulk_filename)}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")
    print(f"Mode: Bulk relaxation only (no formation enthalpy)")

    # Create the WorkGraph
    # By default, all calculation flags are True, but thermodynamics
    # will be skipped since metal_name and oxygen_name are not provided
    print("\nCreating WorkGraph...")
    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name=bulk_filename,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        name='Ag3PO4_BulkRelaxation',
        # Optional: Override default flags if needed
        # compute_relaxation_energy=False,  # Default: True
        # compute_cleavage=False,           # Default: True
        # compute_thermodynamics=False,     # Default: True (but requires metal/oxygen refs)
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
