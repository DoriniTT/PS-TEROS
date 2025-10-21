#!/home/thiagotd/envs/aiida/bin/python
"""
Test script for surface_hydroxylation module - Full Workflow with VASP

WARNING: This script submits actual VASP jobs!
Only run this with proper VASP setup and compute resources.

For quick testing without VASP, use: test_structure_generation.py

Tests:
- Complete SurfaceHydroxylationWorkGraph workflow
- Structure generation
- Parallel VASP relaxations with semaphore control
- Result collection and organization

Material: 2x2 Pt(111) with O adlayer
VASP settings: Lightweight for fast testing
Expected runtime: ~10-20 minutes

Usage:
    source ~/envs/aiida/bin/activate
    python test_full_workflow.py
"""

import sys
from aiida import orm, load_profile
from ase.build import fcc111
from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph


def main():
    """Test complete surface hydroxylation workflow."""

    print("\n" + "="*70)
    print("SURFACE HYDROXYLATION - FULL WORKFLOW TEST (WITH VASP)")
    print("="*70)
    print("\nWARNING: This will submit VASP jobs to the cluster!")
    print("Press Ctrl+C to abort, or Enter to continue...")

    try:
        input()
    except KeyboardInterrupt:
        print("\n\nAborted by user.")
        return 0

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded: psteros")

    # Check daemon is running
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   ERROR: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   ✓ Daemon is running")

    # Create test slab (2x2 Pt(111) with O adlayer)
    print("\n2. Creating test structure...")
    print("   Building 2x2 Pt(111) slab with O adlayer")

    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)

    # Add O atoms on top of Pt surface (manually positioned above Pt atoms)
    import numpy as np
    from ase import Atom
    z_max = max(atom.position[2] for atom in slab)
    top_pt_atoms = [atom for atom in slab if abs(atom.position[2] - z_max) < 0.1]

    # Add O atoms 2.0 Å above each top Pt atom
    for pt_atom in top_pt_atoms:
        o_pos = pt_atom.position.copy()
        o_pos[2] += 2.0  # Place O 2 Å above Pt
        slab.append(Atom('O', position=o_pos))

    slab.center(vacuum=10.0, axis=2)

    structure = orm.StructureData(ase=slab)
    print(f"   ✓ Test slab: {len(slab)} atoms ({sum(1 for a in slab if a.symbol == 'O')} O atoms)")

    # Surface generation parameters
    print("\n3. Setting up surface parameters...")
    surface_params = {
        'mode': 'hydrogen',
        'species': 'O',
        'z_window': 0.5,
        'which_surface': 'top',
        'oh_dist': 0.98,
        'include_empty': False,
        'deduplicate_by_coverage': True,
        'coverage_bins': 3
    }

    print("   Parameters:")
    print(f"   - Mode: {surface_params['mode']}")
    print(f"   - Coverage bins: {surface_params['coverage_bins']}")
    print("   - Expected structures: ~3-5")

    # VASP configuration (LIGHTWEIGHT for testing)
    print("\n4. Setting up VASP configuration...")
    print("   Using LIGHTWEIGHT parameters for fast testing")

    code_label = 'VASP-VTST-6.4.3@bohr'
    try:
        code = orm.load_code(code_label)
        print(f"   ✓ VASP code loaded: {code_label}")
    except Exception as e:
        print(f"\n   ERROR: Could not load VASP code '{code_label}'")
        print(f"   {e}")
        print("\n   Available codes:")
        from aiida.orm import QueryBuilder
        qb = QueryBuilder()
        qb.append(orm.Code, project=['label'])
        for (label,) in qb.iterall():
            print(f"     - {label}")
        return 1

    builder_config = {
        'code': code,
        'parameters': {
            'PREC': 'Normal',      # Lower precision for testing
            'ENCUT': 400,          # Lower cutoff
            'EDIFF': 1e-4,         # Looser convergence
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 2,             # Relax atoms only
            'NSW': 10,             # Very few steps for testing
            'EDIFFG': -0.05,       # Looser force convergence
            'ALGO': 'Fast',
            'LREAL': 'Auto',
            'LWAVE': False,
            'LCHARG': False,
        },
        'potential_family': 'PBE',
        'potential_mapping': {},
        'kpoints_spacing': 0.5,    # Coarse k-points
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 40,
            },
            'queue_name': 'par40',
            'max_wallclock_seconds': 3600,  # 1 hour should be enough
        },
        'clean_workdir': False,  # Keep files for debugging
    }

    print("   VASP parameters:")
    print(f"   - ENCUT: {builder_config['parameters']['ENCUT']} eV")
    print(f"   - NSW: {builder_config['parameters']['NSW']} steps")
    print(f"   - K-points spacing: {builder_config['kpoints_spacing']} Å⁻¹")

    # Parallelization
    max_parallel = 2
    print(f"\n5. Parallelization: max {max_parallel} concurrent jobs")

    # Create and submit workflow
    print("\n6. Creating workflow...")

    from aiida_workgraph import WorkGraph

    # The SurfaceHydroxylationWorkGraph is a @task.graph decorator
    # We need to build it into a WorkGraph and submit
    wg = WorkGraph('test_surface_hydroxylation')

    # Add the main task
    main_task = wg.add_task(
        SurfaceHydroxylationWorkGraph,
        name='hydroxylation',
        structure=structure,
        surface_params=surface_params,
        builder_config=builder_config,
        max_parallel_jobs=max_parallel,
    )

    # Set outputs
    wg.add_output('manifest', main_task.outputs.manifest)
    wg.add_output('successful_relaxations', main_task.outputs.successful_relaxations)
    wg.add_output('failed_relaxations', main_task.outputs.failed_relaxations)
    wg.add_output('statistics', main_task.outputs.statistics)

    print("   ✓ Workflow created")

    print("\n7. Submitting workflow...")
    result = wg.submit()
    pk = result.pk

    print(f"   ✓ Workflow submitted: PK = {pk}")

    print("\n" + "="*70)
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print("="*70)

    print(f"\nWorkflow PK: {pk}")
    print("\nMonitoring commands:")
    print(f"  verdi process show {pk}")
    print(f"  verdi process report {pk}")
    print(f"  verdi process list")
    print("\nWatch progress (updates every 30s):")
    print(f"  watch -n 30 verdi process show {pk}")

    print("\nExpected behavior:")
    print("  1. generate_structures creates ~3 variants")
    print("  2. VASP relaxations run (max 2 parallel)")
    print("  3. collect_results organizes outputs")
    print("  4. Workflow completes with exit status [0]")

    print("\nExpected runtime: ~10-20 minutes")

    print("\nChecking results after completion:")
    print("  python -c \"")
    print("from aiida import orm")
    print(f"node = orm.load_node({pk})")
    print("print('Statistics:', node.outputs.statistics.get_dict())")
    print("print('Successful:', node.outputs.successful_relaxations.get_dict())")
    print("  \"")

    print("\n" + "="*70 + "\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
