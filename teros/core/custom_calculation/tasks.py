"""Helper tasks for custom VASP calculations."""

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def extract_total_energy(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP misc output.

    Args:
        misc: VASP misc output Dict containing energy data

    Returns:
        Total energy as Float (eV)
    """
    misc_dict = misc.get_dict()

    # Navigate to total_energies if present
    energy_dict = misc_dict
    if 'total_energies' in misc_dict:
        energy_dict = misc_dict['total_energies']

    # Try multiple keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(float(energy_dict[key]))

    # If no recognized key found, raise error with available keys
    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in misc output. Available keys: {available}')


@task.calcfunction
def extract_relaxed_structure(misc: orm.Dict) -> orm.StructureData:
    """
    Extract relaxed structure from VASP misc output.

    Args:
        misc: VASP misc output Dict containing structure data

    Returns:
        Relaxed structure as StructureData
    """
    misc_dict = misc.get_dict()

    if 'structure' not in misc_dict:
        raise ValueError(f"No 'structure' key found in misc output. Available keys: {list(misc_dict.keys())}")

    return misc_dict['structure']
