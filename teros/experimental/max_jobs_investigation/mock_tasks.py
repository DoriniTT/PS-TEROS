"""Mock tasks for testing max_number_jobs behavior."""

import time
from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def mock_vasp_calculation(structure: orm.StructureData, label: orm.Str) -> orm.Float:
    """
    Mock VASP calculation that simulates computational work.

    This helps us observe when calculations start/finish to verify
    max_number_jobs behavior.
    """
    print(f"[VASP] Starting calculation: {label.value}")
    time.sleep(2)  # Simulate calculation time
    print(f"[VASP] Finished calculation: {label.value}")

    # Return mock energy
    return orm.Float(42.0)


@task.calcfunction
def create_structures(count: orm.Int):
    """Create multiple mock structures for testing."""
    from ase import Atoms

    structures = {}
    for i in range(count.value):
        # Create simple structure
        atoms = Atoms('H', positions=[[0, 0, 0]])
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)

    return structures
