# coordination_finder.py
import numpy as np


def find_coordinated_oxygens(structure, mn_index, cutoff=2.5):
    """
    Find O atoms coordinated to a specific Mn atom.

    Args:
        structure: pymatgen.core.Structure object
        mn_index: int, index of Mn atom in structure
        cutoff: float, maximum Mn-O distance in Angstroms (default 2.5)

    Returns:
        list: indices of O atoms coordinated to Mn, sorted by distance
    """
    mn_position = structure[mn_index].coords

    # Find all O atoms within cutoff distance
    coordinated_o = []

    for i, site in enumerate(structure):
        if str(site.specie) == "O":
            o_position = site.coords
            distance = np.linalg.norm(mn_position - o_position)

            if distance <= cutoff:
                coordinated_o.append((i, distance))

    # Sort by distance
    coordinated_o.sort(key=lambda x: x[1])

    # Return only indices
    return [idx for idx, dist in coordinated_o]
