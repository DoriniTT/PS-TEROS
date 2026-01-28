"""Calcfunction tasks for the explorer module."""

import re
from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def compute_dynamics(
    structure: orm.StructureData,
    fix_type: orm.Str,
    fix_thickness: orm.Float,
    fix_elements: orm.List = None,
) -> orm.Dict:
    """
    Compute selective dynamics (positions_dof) for a structure.

    This calcfunction is used when the structure is not known at build time
    (e.g., from a previous stage's output), so we need to compute the fixed
    atoms at runtime.

    Args:
        structure: Input structure to analyze
        fix_type: Where to fix atoms ('bottom', 'center', 'top')
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to restrict fixing to

    Returns:
        Dict with 'positions_dof' array for VaspWorkChain dynamics input
    """
    from teros.core.fixed_atoms import get_fixed_atoms_list

    fix_type_val = fix_type.value
    fix_thickness_val = fix_thickness.value
    fix_elements_val = fix_elements.get_list() if fix_elements is not None else None

    # Get list of atom indices to fix (1-based)
    fixed_atoms_list = get_fixed_atoms_list(
        structure=structure,
        fix_type=fix_type_val,
        fix_thickness=fix_thickness_val,
        fix_elements=fix_elements_val,
    )

    # Create positions_dof array: True = relax, False = fix
    num_atoms = len(structure.sites)
    positions_dof = []

    for i in range(1, num_atoms + 1):  # 1-based indexing
        if i in fixed_atoms_list:
            positions_dof.append([False, False, False])  # Fix atom
        else:
            positions_dof.append([True, True, True])  # Relax atom

    return orm.Dict(dict={'positions_dof': positions_dof})


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
