# ooh_constructor.py
import numpy as np
from ase import Atom


def construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos,
                         oo_bond_length=1.45, oh_bond_length=0.97):
    """
    Add O and H atoms to form OOH radical on the exposed O atom.

    Args:
        ase_atoms: ASE Atoms object (will be modified)
        exposed_o_idx: int, index of base O atom
        exposed_o_pos: np.array, position of base O atom [x, y, z]
        oo_bond_length: float, O-O bond length in Angstroms (default 1.45)
        oh_bond_length: float, O-H bond length in Angstroms (default 0.97)

    Returns:
        ASE Atoms object with OOH radical added
    """
    # Create a copy to avoid modifying original
    modified_atoms = ase_atoms.copy()

    # Calculate positions for new O and H atoms
    # Place them along positive z direction (away from surface)

    # New O atom position: base_O + (0, 0, oo_bond_length)
    new_o_position = exposed_o_pos + np.array([0.0, 0.0, oo_bond_length])

    # H atom position: new_O + (0, 0, oh_bond_length)
    h_position = new_o_position + np.array([0.0, 0.0, oh_bond_length])

    # Add new O atom
    modified_atoms.append(Atom('O', position=new_o_position))

    # Add H atom
    modified_atoms.append(Atom('H', position=h_position))

    return modified_atoms
