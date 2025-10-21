# structure_loader.py
from pathlib import Path
from pymatgen.core import Structure


def load_structure(cif_path):
    """
    Load a structure from a CIF file.

    Args:
        cif_path: Path to CIF file (str or Path)

    Returns:
        pymatgen.core.Structure object

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    cif_path = Path(cif_path)

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    structure = Structure.from_file(str(cif_path))
    return structure


def validate_structure(structure):
    """
    Validate that structure has required elements.

    Args:
        structure: pymatgen.core.Structure object

    Returns:
        dict with keys: valid (bool), elements (list), num_mn (int), num_o (int)
    """
    elements = [str(site.specie) for site in structure]
    element_set = set(elements)

    num_mn = elements.count("Mn")
    num_o = elements.count("O")

    valid = "Mn" in element_set and "O" in element_set

    return {
        "valid": valid,
        "elements": list(element_set),
        "num_mn": num_mn,
        "num_o": num_o
    }
