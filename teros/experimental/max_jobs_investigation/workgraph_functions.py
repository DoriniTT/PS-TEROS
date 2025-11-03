"""WorkGraph functions for max_number_jobs investigation."""

import typing as t
from aiida import orm
from aiida_workgraph import task, dynamic


@task.graph
def process_structures_approach_1(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    max_number_jobs: int = None,
) -> t.Annotated[dict, dynamic(orm.Float)]:
    """
    APPROACH 1: Parameter Passing

    Accept max_number_jobs as a parameter and SET it on the WorkGraph context.
    """
    from aiida_workgraph import get_current_graph  # ← CORRECT API!
    from teros.experimental.max_jobs_investigation.mock_tasks import mock_vasp_calculation

    print(f"[APPROACH 1] Function called with max_number_jobs={max_number_jobs}")

    # Get the current WorkGraph and set max_number_jobs
    if max_number_jobs is not None:
        wg = get_current_graph()
        # Convert to int if it's an AiiDA Int node
        if hasattr(max_number_jobs, 'value'):
            max_jobs_value = max_number_jobs.value
        else:
            max_jobs_value = int(max_number_jobs)

        wg.max_number_jobs = max_jobs_value
        print(f"[APPROACH 1] Successfully set wg.max_number_jobs = {max_jobs_value}")

    # Process all structures
    energies = {}
    for key, structure in structures.items():
        print(f"[APPROACH 1] Creating task for {key}")
        energy = mock_vasp_calculation(structure=structure, label=orm.Str(key))
        energies[key] = energy.result

    print(f"[APPROACH 1] Created {len(energies)} calculation tasks")
    return {'energies': energies}


@task.graph
def process_structures_approach_3(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
) -> t.Annotated[dict, dynamic(orm.Float)]:
    """
    APPROACH 3: Framework Cascade

    Check if this nested workgraph inherits max_number_jobs from parent.
    """
    from teros.experimental.max_jobs_investigation.mock_tasks import mock_vasp_calculation

    print("[APPROACH 3] Nested @task.graph function called")

    energies = {}
    for key, structure in structures.items():
        print(f"[APPROACH 3] Creating task for {key}")
        energy = mock_vasp_calculation(structure=structure, label=orm.Str(key))
        energies[key] = energy.result

    print(f"[APPROACH 3] Created {len(energies)} calculation tasks")
    return {'energies': energies}
