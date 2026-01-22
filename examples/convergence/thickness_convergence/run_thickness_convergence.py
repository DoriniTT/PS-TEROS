#!/usr/bin/env python
"""
Slab Thickness Convergence Test for FCC Gold (Au) on obelix cluster.

This script demonstrates the thickness_convergence module for determining
the minimum slab thickness needed for converged surface energy calculations.

The workflow:
1. Relaxes the bulk Au structure with full cell optimization (ISIF=3)
2. Generates Au(111) slabs at multiple thicknesses (3, 5, 7, 9, 11 layers)
3. Relaxes all slabs in parallel (ions only, ISIF=2)
4. Calculates surface energy for each thickness
5. Reports convergence status and recommended thickness

Material: FCC Au
Surface: (111)
Thicknesses: 3, 5, 7, 9, 11 layers
Expected convergence: ~5-7 layers
Expected surface energy: ~0.79 J/m^2

Setup:
    1. Create Au.cif structure file in this directory:
       ```python
       from ase.build import bulk
       from ase.io import write
       au = bulk('Au', 'fcc', a=4.08)
       write('Au.cif', au)
       ```

    2. Activate AiiDA environment:
       source ~/envs/aiida/bin/activate

    3. Run the script:
       python run_thickness_convergence.py

    4. Monitor the workflow:
       verdi process show <PK>

    5. After completion, get results:
       python run_thickness_convergence.py <PK>

Cluster Configuration (obelix):
    - Code: VASP-6.5.1-idefix@obelix
    - Parallelization: Hybrid MPI+OpenMP (4 MPI processes)
    - Nodes: Skylake (88 cores per node)
    - Scheduler: PBS

Usage:
    # Submit new workflow
    python run_thickness_convergence.py

    # Get results from completed workflow
    python run_thickness_convergence.py <PK>
"""

import sys
import os
from aiida import load_profile
from teros.core.convergence import (
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
)


def submit_thickness_convergence():
    """Submit thickness convergence test for Au(111) on obelix cluster."""

    print("\n" + "="*70)
    print("SLAB THICKNESS CONVERGENCE TEST")
    print("Material: FCC Au (111)")
    print("Cluster: obelix (Skylake, PBS)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   Profile loaded: presto")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bulk_structure_path = os.path.join(script_dir, 'Au.cif')

    print(f"\n2. Structure:")
    print(f"   Bulk: {bulk_structure_path}")

    if not os.path.exists(bulk_structure_path):
        print(f"\n   ERROR: Structure file not found!")
        print(f"   Create Au.cif with:")
        print(f"   ```python")
        print(f"   from ase.build import bulk")
        print(f"   from ase.io import write")
        print(f"   au = bulk('Au', 'fcc', a=4.08)")
        print(f"   write('Au.cif', au)")
        print(f"   ```")
        sys.exit(1)

    # Code configuration for obelix cluster
    code_label = 'VASP-6.5.1-idefix@obelix'
    potential_family = 'PBE'

    print(f"\n3. Code configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")
    print(f"   Parallelization: Hybrid MPI+OpenMP (4 MPI procs)")

    # VASP parameters for metals (bulk relaxation)
    bulk_parameters = {
        'prec': 'Accurate',
        'encut': 500,
        'ediff': 1e-6,
        'ismear': 1,      # Methfessel-Paxton for metals
        'sigma': 0.2,
        'ibrion': 2,      # Conjugate gradient
        'isif': 3,        # Full relaxation (ions + cell)
        'nsw': 100,
        'ediffg': -0.01,
        'algo': 'Normal',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }

    # Slab parameters - only relax ionic positions
    slab_parameters = bulk_parameters.copy()
    slab_parameters['isif'] = 2  # Fix cell, relax ions only

    # Scheduler options for obelix (PBS with Skylake nodes)
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 4,  # Hybrid MPI+OpenMP
        },
        'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N Au_conv''',
    }

    # Convergence test parameters
    miller_indices = [1, 1, 1]
    layer_counts = [3, 5, 7, 9, 11]

    print(f"\n4. Convergence test configuration:")
    print(f"   Miller indices: ({miller_indices[0]} {miller_indices[1]} {miller_indices[2]})")
    print(f"   Layer counts: {layer_counts}")
    print(f"   Termination: 0 (first/lowest-index)")
    print(f"   Vacuum thickness: 20.0 Angstrom")
    print(f"   Convergence threshold: 0.01 J/m^2")

    print("\n5. Building WorkGraph...")

    # Build WorkGraph
    wg = build_thickness_convergence_workgraph(
        # Structure
        bulk_structure_path=bulk_structure_path,

        # Code configuration
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping={'Au': 'Au'},
        kpoints_spacing=0.02,  # Fine k-point mesh for bulk
        clean_workdir=False,

        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,

        # Surface configuration
        miller_indices=miller_indices,
        layer_counts=layer_counts,
        min_vacuum_thickness=20.0,
        lll_reduce=True,
        center_slab=False,
        primitive=True,
        termination_index=0,

        # Slab relaxation
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.03,  # Slightly coarser for slabs

        # Convergence settings
        convergence_threshold=0.01,  # J/m^2

        # Concurrency
        max_concurrent_jobs=4,

        # Workflow name
        name='Au_111_thickness_convergence',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n6. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}  # Detailed task hierarchy")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy: Relaxed bulk total energy")
    print(f"  - bulk_structure: Relaxed bulk structure")
    print(f"  - convergence_results: Dict with all surface energies")
    print(f"    - miller_indices: {miller_indices}")
    print(f"    - results: Surface energy for each thickness")
    print(f"    - summary: Convergence status and recommended thickness")
    print(f"\nGet results after completion:")
    print(f"  python {os.path.basename(__file__)} {wg.pk}")
    print(f"\nTypical convergence for Au(111):")
    print(f"  Convergence: ~5-7 layers")
    print(f"  Surface energy: ~0.79 J/m^2")
    print(f"{'='*70}\n")

    return wg


def print_results(pk: int):
    """
    Print results from a completed thickness convergence test.

    Args:
        pk: Process key of the completed WorkGraph
    """
    from aiida_workgraph import WorkGraph
    load_profile()

    print("\n" + "="*70)
    print(f"Loading WorkGraph PK {pk}...")
    print("="*70)

    # Load the completed WorkGraph
    wg = WorkGraph.load(pk)

    # Extract results using the helper function
    results = get_thickness_convergence_results(wg)

    # Print summary
    print("\n" + "="*70)
    print("THICKNESS CONVERGENCE RESULTS")
    print("="*70)

    # Convergence status
    print(f"\nConverged: {results['converged']}")
    if results['converged']:
        print(f"Recommended layers: {results['recommended_layers']}")
    else:
        print("NOT CONVERGED - consider testing more layers or increasing threshold")

    # Bulk energy
    if results['bulk_energy']:
        print(f"\nBulk energy: {results['bulk_energy']:.6f} eV")

    # Surface energies table
    if results['surface_energies']:
        print("\n--- Surface Energy vs. Thickness ---")
        print(f"{'Layers':<10} {'gamma (J/m^2)':<15} {'gamma (eV/Ang^2)':<15}")
        print("-"*40)
        for layers in sorted(results['surface_energies'].keys()):
            gamma_J = results['surface_energies'][layers]
            gamma_eV = gamma_J / 16.0217662  # J/m^2 to eV/Ang^2
            marker = " <-- RECOMMENDED" if layers == results['recommended_layers'] else ""
            print(f"{layers:<10} {gamma_J:>14.4f} {gamma_eV:>14.6f}{marker}")

    # Convergence analysis details
    if results['convergence_results']:
        conv = results['convergence_results']
        summary = conv.get('summary', {})

        print(f"\n--- Convergence Analysis ---")
        print(f"Convergence threshold: {summary.get('convergence_threshold', 'N/A')} J/m^2")
        print(f"Max tested layers: {summary.get('max_tested_layers', 'N/A')}")

        miller = conv.get('miller_indices', [])
        if miller:
            print(f"Miller indices: ({miller[0]} {miller[1]} {miller[2]})")

    # Reference data
    print("\n--- Reference Data (Au 111) ---")
    print("Experimental surface energy: ~0.79 J/m^2")
    print("Typical convergence: 5-7 layers")

    print("\n" + "="*70)
    print()

    return results


if __name__ == '__main__':
    if len(sys.argv) > 1:
        # If PK provided, print results
        try:
            pk = int(sys.argv[1])
            print_results(pk)
        except ValueError:
            print(f"\nError: Invalid PK '{sys.argv[1]}'. Must be an integer.")
            print(f"Usage: python {os.path.basename(__file__)} <PK>")
            sys.exit(1)
        except Exception as e:
            print(f"\nError extracting results: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        # Otherwise, submit new workflow
        try:
            wg = submit_thickness_convergence()
        except Exception as e:
            print(f"\nError submitting workflow: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
