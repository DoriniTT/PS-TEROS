"""WorkGraph tasks for AIMD module."""
from aiida import orm
from aiida_workgraph import task
from aiida.engine import calcfunction


@calcfunction
def create_supercell_calcfunc(
    structure: orm.StructureData,
    spec: orm.List
) -> orm.StructureData:
    """
    Create supercell using pymatgen (as calcfunction).

    Args:
        structure: Input structure
        spec: List [nx, ny, nz] supercell dimensions

    Returns:
        StructureData: Supercell structure
    """
    from pymatgen.io.ase import AseAtomsAdaptor

    # Convert StructureData -> ASE -> pymatgen
    ase_atoms = structure.get_ase()
    adaptor = AseAtomsAdaptor()
    pmg_struct = adaptor.get_structure(ase_atoms)

    # Create supercell
    spec_list = spec.get_list()
    pmg_supercell = pmg_struct * spec_list

    # Convert back: pymatgen -> ASE -> StructureData
    ase_supercell = adaptor.get_atoms(pmg_supercell)
    supercell_data = orm.StructureData(ase=ase_supercell)

    return supercell_data


# Wrap calcfunction as WorkGraph task
create_supercell = task(create_supercell_calcfunc)
