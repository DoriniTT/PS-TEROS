# surface_detector.py
import numpy as np


def find_surface_mn(structure):
    """
    Find the topmost Mn atom in the structure (highest z-coordinate).

    Args:
        structure: pymatgen.core.Structure object

    Returns:
        tuple: (mn_index, mn_position)
            mn_index: int, index of surface Mn in structure
            mn_position: np.array, Cartesian coordinates [x, y, z]
    """
    # Find all Mn atoms
    mn_sites = [(i, site) for i, site in enumerate(structure) if str(site.specie) == "Mn"]

    if not mn_sites:
        raise ValueError("No Mn atoms found in structure")

    # Find Mn with highest z-coordinate (Cartesian)
    surface_mn_idx = None
    max_z = -np.inf

    for idx, site in mn_sites:
        z_coord = site.coords[2]  # Cartesian z
        if z_coord > max_z:
            max_z = z_coord
            surface_mn_idx = idx

    surface_mn_position = structure[surface_mn_idx].coords

    return surface_mn_idx, surface_mn_position
