"""
Utilities to persist generated slab structures to files.

Provides an AiiDA calcfunction that writes a mapping of slab names to StructureData
nodes into a FolderData as CIF files using ASE's writer.
"""
from aiida.engine import calcfunction
from aiida.orm import FolderData, StructureData, Dict
from typing import Dict as TypingDict


@calcfunction
def write_slabs_to_files(slab_structures: Dict) -> FolderData:
    """
    Write slab structures to CIF files and return them bundled in a FolderData.

    Parameters
    ----------
    slab_structures : Dict
        AiiDA Dict-like mapping where the values are StructureData nodes.
        Expected shape matches what `get_slabs` returns under `structures`, e.g.:
        {'s_0': StructureData, 's_1': StructureData, ...}

    Returns
    -------
    FolderData
        A FolderData node containing <name>.cif for each item in the input.

    Notes
    -----
    - Files are named as <key>.cif using the dictionary key, which typically
      follows the pattern 's_0', 's_1', ...
    - ASE is used under the hood via StructureData.get_ase().write(<path>, format="cif").
    """
    folder = FolderData()

    # The input is an AiiDA Dict-like mapping from the WorkGraph port.
    # We iterate through items and write CIF files into the FolderData.
    for name, structure in slab_structures.items():
        if not isinstance(structure, StructureData):
            raise TypeError(f"Value for key '{name}' is not a StructureData: {type(structure)}")
        filename = f"{name}.cif"
        # Write to a temporary local file path managed by FolderData
        with folder.open(filename, 'wb') as handle:
            # ASE write requires a filesystem path; use the underlying name from the open handle
            # We need to close this handle before writing, so instead use a two-step approach:
            pass

    # Re-implement the loop writing via a temporary path and then importing into FolderData
    import tempfile
    import os

    tmpdir = tempfile.mkdtemp(prefix="teros_slabs_")
    try:
        for name, structure in slab_structures.items():
            filename = f"{name}.cif"
            tmp_path = os.path.join(tmpdir, filename)
            structure.get_ase().write(tmp_path, format="cif")
            folder.put_object_from_file(tmp_path, filename)
    finally:
        # Best-effort cleanup of temp dir
        try:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass

    return folder
