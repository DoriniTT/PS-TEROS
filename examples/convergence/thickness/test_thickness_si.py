#!/usr/bin/env python
"""
Test: Si slab thickness convergence with minimal parameters.

This is a lightweight test for the thickness convergence workflow
using Si diamond structure with minimal computational settings.
"""

from aiida import orm, load_profile
from teros.core.convergence import (
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
)


def submit_test():
    """Submit a lightweight Si thickness convergence test."""
    load_profile()

    # Simple Si diamond structure (2 atoms in primitive cell)
    from ase.build import bulk
    si_bulk = bulk('Si', 'diamond', a=5.43)
    bulk_structure = orm.StructureData(ase=si_bulk)

    # Minimal VASP parameters for quick testing
    bulk_parameters = {
        'PREC': 'Normal',
        'ENCUT': 300,       # Low cutoff for speed
        'EDIFF': 1e-4,      # Loose electronic convergence
        'ISMEAR': 0,        # Gaussian smearing for semiconductor
        'SIGMA': 0.1,
        'IBRION': 2,
        'ISIF': 3,          # Full relaxation for bulk
        'NSW': 20,          # Few ionic steps
        'EDIFFG': -0.05,    # Loose force convergence
        'ALGO': 'Fast',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Slab parameters - only relax ionic positions
    slab_parameters = bulk_parameters.copy()
    slab_parameters['ISIF'] = 2   # Fix cell, relax ions
    slab_parameters['NSW'] = 15   # Fewer ionic steps for slabs

    # Options for localwork (8 processors)
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
        },
        'max_wallclock_seconds': 3600,  # 1 hour max
    }

    # Build WorkGraph
    wg = build_thickness_convergence_workgraph(
        # Structure
        bulk_structure=bulk_structure,

        # VASP configuration
        code_label='VASP-6.5.1@localwork',
        potential_family='PBE',
        potential_mapping={'Si': 'Si'},
        kpoints_spacing=0.08,        # Coarse k-points for speed
        clean_workdir=False,

        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,

        # Surface configuration - Si(100)
        miller_indices=[1, 0, 0],
        layer_counts=[3, 5, 7],       # Test 3 thicknesses only
        min_vacuum_thickness=12.0,    # 12 A vacuum
        lll_reduce=True,
        center_slab=True,
        primitive=True,
        termination_index=0,

        # Slab relaxation
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.08,

        # Convergence settings
        convergence_threshold=0.05,   # Relaxed threshold for testing

        # Concurrency - localwork runs 1 job at a time
        max_concurrent_jobs=1,

        # Workflow name
        name='Si_100_thickness_test',
    )

    # Submit (non-blocking)
    wg.submit(wait=False)
    print(f"Submitted Si thickness convergence test: PK {wg.pk}")
    print(f"Monitor with: verdi process show {wg.pk}")

    return wg


def get_results(pk: int):
    """Get results from completed workflow."""
    from aiida_workgraph import WorkGraph
    load_profile()

    wg = WorkGraph.load(pk)
    results = get_thickness_convergence_results(wg)

    print("\n" + "=" * 60)
    print("THICKNESS CONVERGENCE RESULTS")
    print("=" * 60)

    print(f"\nConverged: {results['converged']}")
    if results['converged']:
        print(f"Recommended layers: {results['recommended_layers']}")

    if results['bulk_energy']:
        print(f"\nBulk energy: {results['bulk_energy']:.6f} eV")

    if results['surface_energies']:
        print("\n--- Surface Energy vs. Thickness ---")
        for layers in sorted(results['surface_energies'].keys()):
            gamma = results['surface_energies'][layers]
            print(f"  {layers} layers: {gamma:.4f} J/m^2")

    print("=" * 60)
    return results


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        pk = int(sys.argv[1])
        get_results(pk)
    else:
        submit_test()
