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
    energy_dict = misc_dict.get('total_energies', misc_dict)

    # Try multiple keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    raise ValueError(f"No recognized energy key found in misc output. Available keys: {list(energy_dict.keys())}")
