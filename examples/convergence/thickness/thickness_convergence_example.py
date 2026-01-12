#!/usr/bin/env python
"""
Example: Slab thickness convergence testing with PS-TEROS.

This script demonstrates how to use the thickness convergence module to determine
the minimum number of atomic layers needed for converged surface energy calculations.

Usage:
    1. Modify the bulk structure path or load a structure from AiiDA
    2. Adjust builder_inputs for your cluster
    3. Set miller_indices and layer_counts for your surface
    4. Run the script to submit the workflow
    5. After completion, run get_results(pk) to extract recommendations
"""

import os
from aiida import orm, load_profile
from teros.core.convergence import (
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
)


def submit_thickness_convergence():
    """Submit a thickness convergence test workflow."""
    # Load AiiDA profile
    load_profile()

    # Example: Create a simple FCC Au structure for testing
    # In practice, you would load from file or database
    from ase.build import bulk
    au_bulk = bulk('Au', 'fcc', a=4.08)
    bulk_structure = orm.StructureData(ase=au_bulk)

    # Alternative: Load from file
    # bulk_structure_path = '/absolute/path/to/bulk.cif'

    # VASP parameters for bulk relaxation
    bulk_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 500,
        'EDIFF': 1e-6,
        'ISMEAR': 1,      # Methfessel-Paxton for metals
        'SIGMA': 0.2,
        'IBRION': 2,
        'ISIF': 3,        # Full relaxation for bulk
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Slab parameters - only relax ionic positions
    slab_parameters = bulk_parameters.copy()
    slab_parameters['ISIF'] = 2  # Fix cell, relax ions

    # Scheduler options
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
        'queue_name': 'normal',
        'walltime_seconds': 86400,  # 24 hours
    }

    # Build WorkGraph
    wg = build_thickness_convergence_workgraph(
        # Structure - provide either path or StructureData
        bulk_structure=bulk_structure,
        # bulk_structure_path='/absolute/path/to/bulk.cif',  # Alternative

        # VASP configuration
        code_label='VASP-6.5.1@cluster02',  # Replace with your code label
        potential_family='PBE.54',
        potential_mapping={'Au': 'Au'},
        kpoints_spacing=0.03,
        clean_workdir=False,

        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,

        # Surface configuration
        miller_indices=[1, 1, 1],           # Test Au(111)
        layer_counts=[3, 5, 7, 9, 11],      # Test these thicknesses
        min_vacuum_thickness=20.0,          # 20 A vacuum
        lll_reduce=True,
        center_slab=True,
        primitive=True,
        termination_index=0,                # First termination

        # Slab relaxation (optional - defaults to bulk parameters)
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.03,

        # Convergence settings
        convergence_threshold=0.01,  # J/m^2

        # Concurrency
        max_concurrent_jobs=4,

        # Workflow name
        name='Au_111_thickness_convergence',
    )

    # Submit (non-blocking)
    wg.submit(wait=False)
    print(f"Submitted thickness convergence test: PK {wg.pk}")
    print(f"Monitor with: verdi process show {wg.pk}")
    print(f"\nAfter completion, get results with:")
    print(f"  python {os.path.basename(__file__)} {wg.pk}")

    return wg


def get_results(pk: int):
    """
    Get results from a completed thickness convergence test.

    Args:
        pk: Process key of the completed WorkGraph
    """
    from aiida_workgraph import WorkGraph
    load_profile()

    # Load the completed WorkGraph
    wg = WorkGraph.load(pk)

    # Extract results
    results = get_thickness_convergence_results(wg)

    # Print summary
    print("\n" + "=" * 70)
    print("THICKNESS CONVERGENCE RESULTS")
    print("=" * 70)

    print(f"\nConverged: {results['converged']}")
    if results['converged']:
        print(f"Recommended layers: {results['recommended_layers']}")
    else:
        print("NOT CONVERGED - consider testing more layers or increasing threshold")

    if results['bulk_energy']:
        print(f"\nBulk energy: {results['bulk_energy']:.6f} eV")

    if results['surface_energies']:
        print("\n--- Surface Energy vs. Thickness ---")
        print(f"{'Layers':<10} {'gamma (J/m^2)':<15}")
        print("-" * 25)
        for layers in sorted(results['surface_energies'].keys()):
            gamma = results['surface_energies'][layers]
            print(f"{layers:<10} {gamma:.4f}")

    # Print convergence analysis
    if results['convergence_results']:
        conv = results['convergence_results']
        summary = conv.get('summary', {})
        print(f"\nConvergence threshold: {summary.get('convergence_threshold', 'N/A')} J/m^2")
        print(f"Max tested layers: {summary.get('max_tested_layers', 'N/A')}")

        miller = conv.get('miller_indices', [])
        if miller:
            print(f"Miller indices: ({miller[0]}{miller[1]}{miller[2]})")

    print("\n" + "=" * 70)

    return results


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        # If PK provided, get results
        pk = int(sys.argv[1])
        get_results(pk)
    else:
        # Otherwise, submit new workflow
        submit_thickness_convergence()
