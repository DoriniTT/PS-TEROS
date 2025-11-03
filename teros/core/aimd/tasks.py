"""WorkGraph tasks for AIMD module."""
from aiida import orm
from aiida_workgraph import task


@task.pythonjob()
def create_supercell(
    structure: orm.StructureData,
    spec: list[int]
):
    """
    Create supercell using pymatgen.

    When running as AiiDA PythonJob, structure is unpacked to ASE Atoms.
    We return ASE Atoms, which pythonjob automatically converts to StructureData.

    Args:
        structure: Input structure (unpacked to ASE Atoms in execution)
        spec: [nx, ny, nz] supercell dimensions

    Returns:
        ASE Atoms object (auto-converted to StructureData by pythonjob)
    """
    from pymatgen.io.ase import AseAtomsAdaptor

    # When running in AiiDA, structure is unpacked to ASE Atoms
    # Convert ASE -> pymatgen
    adaptor = AseAtomsAdaptor()
    pmg_struct = adaptor.get_structure(structure)

    # Create supercell
    pmg_supercell = pmg_struct * spec

    # Convert back: pymatgen -> ASE
    # PythonJob will automatically convert ASE Atoms to StructureData
    ase_supercell = adaptor.get_atoms(pmg_supercell)

    return ase_supercell
