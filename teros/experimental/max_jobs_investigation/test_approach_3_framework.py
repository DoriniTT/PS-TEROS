#!/usr/bin/env python3
"""
Test Approach 3: Framework Properties

Investigate if WorkGraph has properties that cascade to nested graphs.
Test if setting max_number_jobs on parent affects children.
"""

import typing as t
from aiida import orm, load_profile
from aiida_workgraph import task, WorkGraph, dynamic

load_profile()

from teros.experimental.max_jobs_investigation import process_structures_approach_3


def test_framework_cascade():
    """Test if max_number_jobs cascades from parent to nested workgraphs."""
    print("\n" + "="*80)
    print("TEST: Approach 3 - Framework Properties (Cascading)")
    print("="*80 + "\n")

    # Create test structures
    from ase import Atoms

    structures = {}
    for i in range(5):
        atoms = Atoms('H', positions=[[0, 0, 0]])
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)

    print(f"Created {len(structures)} test structures")

    # Create main workgraph
    with WorkGraph("test_approach_3") as wg:
        # Call the nested @task.graph function
        # Question: Will it inherit max_number_jobs=2 from parent?
        result = process_structures_approach_3(structures=structures)

    # Set max_number_jobs on MAIN workgraph
    wg.max_number_jobs = 2
    print(f"[MAIN] Set main workgraph max_number_jobs = 2")

    print("\n[MAIN] Submitting workgraph...")
    wg.submit(wait=True)

    print(f"\n[MAIN] WorkGraph finished: PK={wg.pk}")
    print(f"[MAIN] State: {wg.state}")

    # Analyze the results
    print("\n" + "="*80)
    print("ANALYSIS")
    print("="*80)
    print("\nIf max_number_jobs cascaded to nested workgraph:")
    print("  → Only 2 VASP tasks should run concurrently")
    print("  → Check process report for timing")
    print("\nIf it DIDN'T cascade:")
    print("  → All 5 VASP tasks started simultaneously")
    print("\nCheck with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")

    return wg.pk


def test_framework_properties():
    """Inspect WorkGraph class for cascade-related properties."""
    print("\n" + "="*80)
    print("TEST: Inspect WorkGraph Properties")
    print("="*80 + "\n")

    import inspect
    from aiida_workgraph import WorkGraph

    # Check WorkGraph class
    print("[INSPECT] WorkGraph class attributes:")
    for attr in dir(WorkGraph):
        if 'max' in attr.lower() or 'concurrent' in attr.lower() or 'limit' in attr.lower():
            print(f"  - {attr}")

    # Check WorkGraph instance
    wg = WorkGraph("inspect_test")
    wg.max_number_jobs = 2

    print(f"\n[INSPECT] After setting max_number_jobs=2:")
    print(f"  wg.max_number_jobs = {wg.max_number_jobs}")

    # Check for other related attributes
    print(f"\n[INSPECT] Checking for cascade-related attributes:")
    for attr in ['cascade_limits', 'inherit_limits', 'propagate_limits', 'global_max_jobs']:
        if hasattr(wg, attr):
            print(f"  ✓ {attr} = {getattr(wg, attr)}")
        else:
            print(f"  ✗ {attr} not found")


if __name__ == "__main__":
    # First, inspect properties
    test_framework_properties()

    # Then run the cascade test
    print("\n" + "="*80 + "\n")
    pk = test_framework_cascade()
    print(f"\nWorkGraph PK: {pk}")
