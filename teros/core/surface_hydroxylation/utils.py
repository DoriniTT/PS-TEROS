"""Utility functions for surface_hydroxylation module."""

from aiida.orm import StructureData
from ase import Atoms


def aiida_to_ase(structure: StructureData) -> Atoms:
    """
    Convert AiiDA StructureData to ASE Atoms.

    Args:
        structure: AiiDA StructureData object

    Returns:
        ASE Atoms object
    """
    return structure.get_ase()


def ase_to_aiida(atoms: Atoms) -> StructureData:
    """
    Convert ASE Atoms to AiiDA StructureData.

    Args:
        atoms: ASE Atoms object

    Returns:
        AiiDA StructureData object
    """
    return StructureData(ase=atoms)
