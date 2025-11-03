"""WorkGraph tasks for AIMD module."""
from aiida import orm
from aiida_workgraph import task


@task
def create_supercell(
    structure: orm.StructureData,
    spec: list[int]
) -> orm.StructureData:
    """
    Create supercell using pymatgen.

    Args:
        structure: Input structure
        spec: [nx, ny, nz] supercell dimensions

    Returns:
        Supercell StructureData
    """
    # Convert to pymatgen
    pmg_struct = structure.get_pymatgen()

    # Create supercell
    pmg_supercell = pmg_struct * spec

    # Convert back to AiiDA
    supercell = orm.StructureData(pymatgen=pmg_supercell)

    return supercell
