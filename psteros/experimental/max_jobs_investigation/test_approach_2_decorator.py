#!/usr/bin/env python3
"""
Test Approach 2: Decorator Configuration

Try passing max_number_jobs directly to the @task.graph decorator.
This is experimental - the documentation doesn't show this syntax.
"""

import typing as t
from aiida import orm, load_profile
from aiida_workgraph import task, WorkGraph, dynamic

load_profile()


# Test 1: Try @task.graph(max_number_jobs=2)
print("\n[APPROACH 2] Testing @task.graph(max_number_jobs=2) syntax...")

try:
    @task.graph(max_number_jobs=2)
    def process_with_decorator_config(
        structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    ) -> t.Annotated[dict, dynamic(orm.Float)]:
        """Try passing max_number_jobs to decorator."""
        from mock_tasks import mock_vasp_calculation

        print("[APPROACH 2] Inside decorator-configured function")

        energies = {}
        for key, structure in structures.items():
            energy = mock_vasp_calculation(structure=structure, label=orm.Str(key))
            energies[key] = energy.result

        return {'energies': energies}

    print("[APPROACH 2] ✓ Decorator syntax accepted!")

except TypeError as e:
    print(f"[APPROACH 2] ✗ Decorator doesn't accept max_number_jobs: {e}")
    process_with_decorator_config = None


# Test 2: Try other possible decorator parameters
print("\n[APPROACH 2] Testing other decorator parameters...")

try:
    # Check what parameters @task.graph accepts
    import inspect
    from aiida_workgraph import task

    # Get the graph decorator
    graph_decorator = task.graph

    print(f"[APPROACH 2] @task.graph type: {type(graph_decorator)}")
    print(f"[APPROACH 2] @task.graph signature: {inspect.signature(graph_decorator)}")

except Exception as e:
    print(f"[APPROACH 2] Error inspecting decorator: {e}")


def test_decorator_configuration():
    """Test if decorator accepts max_number_jobs."""
    print("\n" + "="*80)
    print("TEST: Approach 2 - Decorator Configuration")
    print("="*80 + "\n")

    if process_with_decorator_config is None:
        print("[APPROACH 2] Decorator syntax not supported - test skipped")
        return None

    # Create test structures
    from ase import Atoms

    structures = {}
    for i in range(5):
        atoms = Atoms('H', positions=[[0, 0, 0]])
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)

    print(f"Created {len(structures)} test structures")

    # Create main workgraph
    wg = WorkGraph("test_approach_2")

    # Call the function
    result = process_with_decorator_config(structures=structures)

    print("\n[MAIN] Submitting workgraph...")
    wg.submit(wait=True)

    print(f"\n[MAIN] WorkGraph finished: PK={wg.pk}")
    print(f"[MAIN] Exit status: {wg.exit_status}")

    return wg.pk


if __name__ == "__main__":
    pk = test_decorator_configuration()
    if pk:
        print(f"\nWorkGraph PK: {pk}")
    else:
        print("\nTest could not be run - decorator syntax not supported")
