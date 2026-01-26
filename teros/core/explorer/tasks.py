"""Calcfunction tasks for the explorer module."""

import re
from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def extract_energy(misc: orm.Dict, retrieved: orm.FolderData = None) -> orm.Float:
    """
    Extract total energy from VASP misc output or OUTCAR.

    Args:
        misc: VASP misc output Dict containing energy data
        retrieved: VASP retrieved FolderData (optional) to parse OUTCAR if misc fails

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

    # If no recognized key found in misc, try to parse from retrieved OUTCAR
    if retrieved is not None:
        try:
            content = retrieved.get_object_content('OUTCAR')
            # Look for "free  energy   TOTEN  =       -832.63657516 eV"
            matches = re.findall(r'free\s+energy\s+TOTEN\s+=\s+([-\d.]+)', content)
            if matches:
                return orm.Float(float(matches[-1]))
        except Exception:
            pass

    # If no recognized key found, raise error with available keys
    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in misc output or OUTCAR. Available keys: {available}')
