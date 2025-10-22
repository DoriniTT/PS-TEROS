# structure_exporter.py
from pathlib import Path
from ase.io import write


def export_structure_to_cif(ase_atoms, output_path):
    """
    Export ASE Atoms object to CIF file.

    Args:
        ase_atoms: ASE Atoms object
        output_path: Path to output CIF file (str or Path)
    """
    output_path = Path(output_path)

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write CIF file
    write(str(output_path), ase_atoms, format='cif')
