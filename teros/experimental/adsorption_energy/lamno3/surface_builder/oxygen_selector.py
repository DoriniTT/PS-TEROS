# oxygen_selector.py
import numpy as np


def select_most_exposed_oxygen(structure, o_indices):
    """
    Select the most exposed O atom (highest z-coordinate) from a list of O indices.

    Args:
        structure: pymatgen.core.Structure object
        o_indices: list of int, indices of O atoms to consider

    Returns:
        tuple: (exposed_o_index, exposed_o_position)
            exposed_o_index: int, index of most exposed O
            exposed_o_position: np.array, Cartesian coordinates [x, y, z]

    Raises:
        ValueError: If o_indices is empty
    """
    if not o_indices:
        raise ValueError("No oxygen atoms provided")

    # Find O with highest z-coordinate
    max_z = -np.inf
    exposed_o_idx = None

    for idx in o_indices:
        z_coord = structure[idx].coords[2]
        if z_coord > max_z:
            max_z = z_coord
            exposed_o_idx = idx

    exposed_o_position = structure[exposed_o_idx].coords

    return exposed_o_idx, exposed_o_position
