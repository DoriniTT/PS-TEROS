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


@task.calcfunction
def extract_total_energy_from_misc(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP WorkChain misc output.

    Args:
        misc: Dict node containing various VASP outputs including total_energies

    Returns:
        Float node with final total energy in eV
    """
    misc_dict = misc.get_dict()

    if 'total_energies' not in misc_dict:
        raise ValueError("total_energies not found in misc output")

    energy_dict = misc_dict['total_energies']

    # Try standard keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(float(energy_dict[key]))

    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy. Available keys: {available}')
