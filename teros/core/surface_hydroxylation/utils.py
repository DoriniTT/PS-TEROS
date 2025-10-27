"""Utility functions for surface_hydroxylation module."""

from aiida import orm
from aiida_workgraph import task
from ase import Atoms


def aiida_to_ase(structure: orm.StructureData) -> Atoms:
    """
    Convert AiiDA StructureData to ASE Atoms.

    Args:
        structure: AiiDA StructureData object

    Returns:
        ASE Atoms object
    """
    return structure.get_ase()


def ase_to_aiida(atoms: Atoms) -> orm.StructureData:
    """
    Convert ASE Atoms to AiiDA StructureData.

    Args:
        atoms: ASE Atoms object

    Returns:
        AiiDA StructureData object
    """
    return orm.StructureData(ase=atoms)


@task()
def extract_total_energy_from_misc(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP WorkChain misc output.

    Args:
        misc: Dict node containing various VASP outputs including total_energies

    Returns:
        Float node with final total energy in eV
    """
    misc_dict = misc.get_dict()

    # Get total_energies array from misc
    if 'total_energies' not in misc_dict:
        raise ValueError("total_energies not found in misc output")

    total_energies = misc_dict['total_energies']

    if not total_energies:
        raise ValueError("total_energies array is empty")

    # Return final energy (last ionic step)
    final_energy = float(total_energies[-1])
    return orm.Float(final_energy)
