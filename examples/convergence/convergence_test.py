#!/usr/bin/env python
"""
Example: VASP convergence testing with PS-TEROS.

This script demonstrates how to use the convergence module to determine
optimal ENCUT and k-points spacing for a given structure.

Usage:
    1. Modify the structure PK or load a structure from file
    2. Adjust builder_inputs for your cluster
    3. Adjust conv_settings for your desired scan ranges
    4. Run the script to submit the workflow
    5. After completion, run get_results() to extract recommendations
"""

from aiida import orm, load_profile
from teros.core.convergence import build_convergence_workgraph, get_convergence_results


def submit_convergence_test():
    """Submit a convergence test workflow."""
    # Load AiiDA profile
    load_profile()

    # Load structure - replace with your structure PK or load from file
    # Option 1: Load from database by PK
    # structure = orm.load_node(12345)

    # Option 2: Load from file using pymatgen
    # from pymatgen.core import Structure
    # from aiida.orm import StructureData
    # pmg_structure = Structure.from_file('my_structure.cif')
    # structure = StructureData(pymatgen=pmg_structure)

    # For this example, we'll create a simple Ag structure
    from ase.build import bulk
    from aiida.orm import StructureData
    ase_structure = bulk('Ag', 'fcc', a=4.09)
    structure = StructureData(ase=ase_structure)

    # Define base VASP parameters
    builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ISMEAR': 0,       # Gaussian smearing (good for convergence tests)
                'SIGMA': 0.05,
                'EDIFF': 1e-6,
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
                'NELM': 100,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
            'queue_name': 'normal',
            'walltime_seconds': 3600,  # 1 hour per calculation
        },
        'kpoints_spacing': 0.05,  # Starting value (will be varied in k-points scan)
        'potential_family': 'PBE.54',
        'potential_mapping': {'Ag': 'Ag'},
        'clean_workdir': True,    # Clean up scratch files
    }

    # Define convergence scan ranges
    conv_settings = {
        # Cutoff energy scan parameters
        'cutoff_start': 300,      # Start ENCUT (eV)
        'cutoff_stop': 600,       # End ENCUT (eV)
        'cutoff_step': 50,        # ENCUT step (eV)

        # K-points spacing scan parameters
        'kspacing_start': 0.06,   # Start k-spacing (A^-1) - coarsest
        'kspacing_stop': 0.02,    # End k-spacing (A^-1) - finest
        'kspacing_step': -0.01,   # K-spacing step (negative = decreasing)

        # Fixed values during the other scan
        'cutoff_kconv': 520,      # ENCUT used during k-points scan
        'kspacing_cutconv': 0.03, # K-spacing used during ENCUT scan
    }

    # Build WorkGraph
    wg = build_convergence_workgraph(
        structure=structure,
        code_label='VASP-6.5.1@cluster02',  # Replace with your VASP code label
        builder_inputs=builder_inputs,
        conv_settings=conv_settings,
        convergence_threshold=0.001,  # 1 meV/atom
        name='Ag_convergence_test',
    )

    # Submit (non-blocking)
    wg.submit(wait=False)
    print(f"Submitted convergence test: PK {wg.pk}")
    print(f"Monitor with: verdi process show {wg.pk}")

    return wg


def get_results(pk: int):
    """
    Get results from a completed convergence test.

    Args:
        pk: Process key of the completed WorkGraph
    """
    from aiida import load_profile
    from aiida_workgraph import WorkGraph
    load_profile()

    # Load the completed WorkGraph
    wg = WorkGraph.from_pk(pk)

    # Extract results
    results = get_convergence_results(wg)

    # Print summary
    print("\n" + "=" * 60)
    print("CONVERGENCE TEST RESULTS")
    print("=" * 60)

    if results['recommended_cutoff']:
        print(f"\nRecommended ENCUT: {results['recommended_cutoff']} eV")
    else:
        print("\nENCUT: NOT CONVERGED - increase cutoff_stop")

    if results['recommended_kspacing']:
        print(f"Recommended k-spacing: {results['recommended_kspacing']} A^-1")
    else:
        print("K-spacing: NOT CONVERGED - decrease kspacing_stop")

    if results['convergence_summary']:
        summary = results['convergence_summary']
        print(f"\nThreshold used: {summary['threshold_used'] * 1000:.1f} meV/atom")
        print(f"Cutoff converged at: {summary['cutoff_converged_at']} eV")
        print(f"K-spacing converged at: {summary['kspacing_converged_at']} A^-1")

    # Print detailed analysis if available
    if results['cutoff_analysis']:
        print("\n--- Cutoff Analysis ---")
        analysis = results['cutoff_analysis']
        print(f"Tested values: {analysis['cutoff_values']}")
        print(f"Energy per atom (eV): {[f'{e:.6f}' for e in analysis['energy_per_atom']]}")
        print(f"Diff from max (meV/atom): {[f'{d*1000:.2f}' for d in analysis['energy_diff_from_max']]}")

    if results['kpoints_analysis']:
        print("\n--- K-points Analysis ---")
        analysis = results['kpoints_analysis']
        print(f"Tested values: {analysis['kspacing_values']}")
        print(f"Energy per atom (eV): {[f'{e:.6f}' for e in analysis['energy_per_atom']]}")
        print(f"Diff from densest (meV/atom): {[f'{d*1000:.2f}' for d in analysis['energy_diff_from_densest']]}")

    print("\n" + "=" * 60)

    return results


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        # If PK provided, get results
        pk = int(sys.argv[1])
        get_results(pk)
    else:
        # Otherwise, submit new workflow
        submit_convergence_test()
