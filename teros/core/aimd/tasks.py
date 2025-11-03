"""WorkGraph tasks for AIMD module."""
from aiida import orm
from aiida_workgraph import task


@task.pythonjob()
def create_supercell(
    structure: orm.StructureData,
    spec: list[int]
) -> orm.StructureData:
    """
    Create supercell using pymatgen.

    When running as AiiDA task, structure is automatically unpacked to ASE Atoms.
    We convert ASE -> pymatgen -> supercell -> ASE -> StructureData.

    Args:
        structure: Input structure (unpacked to ASE Atoms in execution)
        spec: [nx, ny, nz] supercell dimensions

    Returns:
        Supercell StructureData
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from aiida import orm

    # When running in AiiDA, structure is unpacked to ASE Atoms
    # Convert ASE -> pymatgen
    adaptor = AseAtomsAdaptor()
    pmg_struct = adaptor.get_structure(structure)

    # Create supercell
    pmg_supercell = pmg_struct * spec

    # Convert back: pymatgen -> ASE -> StructureData
    ase_supercell = adaptor.get_atoms(pmg_supercell)
    supercell_data = orm.StructureData(ase=ase_supercell)

    return supercell_data
