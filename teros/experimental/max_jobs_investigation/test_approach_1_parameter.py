#!/usr/bin/env python3
"""
Test Approach 1: Parameter Passing

Add max_number_jobs as a parameter to @task.graph functions.
Question: How do we access the WorkGraph instance inside @task.graph to set it?
"""

import typing as t
from aiida import orm, load_profile
from aiida_workgraph import task, WorkGraph, dynamic

load_profile()

from teros.experimental.max_jobs_investigation import process_structures_approach_1


def test_parameter_passing():
    """Test if we can pass max_number_jobs as parameter."""
    print("\n" + "="*80)
    print("TEST: Approach 1 - Parameter Passing")
    print("="*80 + "\n")

    # Create test structures
    from ase import Atoms

    structures = {}
    for i in range(5):
        atoms = Atoms('H', positions=[[0, 0, 0]])
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)

    print(f"Created {len(structures)} test structures")

    # Create main workgraph
    with WorkGraph("test_approach_1") as wg:
        # Call the @task.graph function with max_number_jobs parameter
        result = process_structures_approach_1(
            structures=structures,
            max_number_jobs=2,  # Try to limit to 2 concurrent jobs
        )

    print("\n[MAIN] Submitting workgraph...")
    wg.submit(wait=True)

    print(f"\n[MAIN] WorkGraph finished: PK={wg.pk}")
    print(f"[MAIN] State: {wg.state}")

    # Check if max_number_jobs was respected
    # (We'll need to check the process report to see timing)
    print(f"\n[MAIN] Check timing with: verdi process report {wg.pk}")

    return wg.pk


if __name__ == "__main__":
    pk = test_parameter_passing()
    print(f"\nWorkGraph PK: {pk}")
