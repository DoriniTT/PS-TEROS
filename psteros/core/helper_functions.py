"""
Helper functions for PS-TEROS WorkGraph workflows.

This module contains generic helper functions that are used across
different workflows. These functions do NOT contain specific calculation
parameters - those should be defined in the launch scripts.
"""

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def get_structure_from_file(filepath) -> orm.StructureData:
    """
    Load structure from a file (CIF, POSCAR, XYZ, etc.).

    This is a generic function that uses ASE to read structure files
    in any format supported by ASE.

    Args:
        filepath: Path to the structure file (can be str or orm.Str)

    Returns:
        StructureData node containing the loaded structure
    """
    from ase.io import read

    # Handle both string and orm.Str types
    if isinstance(filepath, orm.Str):
        path_str = filepath.value
    else:
        path_str = str(filepath)

    # Read structure file using ASE (supports CIF, POSCAR, XYZ, etc.)
    atoms = read(path_str)

    # Convert to AiiDA StructureData
    structure = orm.StructureData(ase=atoms)

    return structure


def prepare_vasp_inputs(
    code_label: str,
    parameters: dict,
    options: dict,
    potential_mapping: dict,
    structure: orm.StructureData = None,
    structure_path: str = None,
    kpoints_spacing: float = 0.3,
    potential_family: str = 'PBE.54',
    clean_workdir: bool = True,
):
    """
    Prepare AiiDA nodes from Python dictionaries for VASP inputs.

    This is a generic helper that converts Python types to AiiDA data nodes.
    All specific parameters should be provided by the caller (launch script).

    Args:
        code_label: Label of the VASP code in AiiDA
        parameters: VASP input parameters dictionary
        options: Scheduler options dictionary
        potential_mapping: Mapping of elements to potentials
        structure: Pre-loaded structure (optional)
        structure_path: Path to structure file (optional)
        kpoints_spacing: K-points spacing in A^-1 * 2pi
        potential_family: Name of the potential family
        clean_workdir: Whether to clean work directory after completion

    Returns:
        dict: Dictionary containing all prepared AiiDA nodes

    Raises:
        ValueError: If required inputs are missing
    """
    # Validate inputs
    if structure is None and structure_path is None:
        raise ValueError("Either structure or structure_path must be provided")

    if parameters is None:
        raise ValueError("parameters dictionary must be provided")

    if options is None:
        raise ValueError("options dictionary must be provided")

    if potential_mapping is None:
        raise ValueError("potential_mapping dictionary must be provided")

    # Load code
    code = orm.load_code(code_label)

    # Convert Python dictionaries to AiiDA nodes
    params_dict = orm.Dict(dict=parameters)
    options_dict = orm.Dict(dict=options)
    pot_mapping_dict = orm.Dict(dict=potential_mapping)

    return {
        'code': code,
        'params_dict': params_dict,
        'options_dict': options_dict,
        'pot_mapping_dict': pot_mapping_dict,
        'kpoints_spacing': orm.Float(kpoints_spacing),
        'potential_family': orm.Str(potential_family),
        'clean_workdir': orm.Bool(clean_workdir),
        'structure': structure,
        'structure_path': structure_path,
    }
