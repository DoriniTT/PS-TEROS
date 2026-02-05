"""Calcfunction tasks for the lego module."""

import re
import typing as t

from aiida import orm
from aiida_workgraph import task, dynamic


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


@task.calcfunction
def prepare_poscar_with_velocities(
    structure: orm.StructureData,
    velocities_dict: orm.Dict,
) -> orm.Dict:
    """
    Prepare POSCAR content with velocity vectors from previous AIMD stage.

    This calcfunction creates a VASP POSCAR string that includes atomic
    velocities from a previous CONTCAR, enabling seamless MD continuation
    without thermal re-initialization.

    Args:
        structure: StructureData for POSCAR generation
        velocities_dict: Dict from extract_velocities_from_contcar with velocity data

    Returns:
        orm.Dict with keys:
            - 'poscar_content': POSCAR string with velocities block (if available)
            - 'has_velocities': bool indicating velocity block presence
    """
    from .utils import build_poscar_with_velocities

    vel_dict = velocities_dict.get_dict()
    has_vels = vel_dict.get('has_velocities', False)
    vels_data = vel_dict.get('velocities', [])

    import numpy as np
    velocities = None
    if has_vels and len(vels_data) == len(structure.sites):
        velocities = np.array(vels_data, dtype=float)

    # Build POSCAR with or without velocities
    poscar_str = build_poscar_with_velocities(structure, velocities=velocities)

    return orm.Dict(dict={
        'poscar_content': poscar_str,
        'has_velocities': has_vels,
        'n_atoms': len(structure.sites),
    })


@task.calcfunction
def create_poscar_file_with_velocities(
    structure: orm.StructureData,
    velocities_dict: orm.Dict,
) -> orm.SinglefileData:
    """
    Create a SinglefileData POSCAR with embedded velocities for AimdVaspCalculation.

    This calcfunction creates a POSCAR file (as SinglefileData) that includes
    atomic velocities from a previous AIMD stage. The resulting file can be
    passed to AimdVaspCalculation's poscar_file input for seamless MD continuation.

    Args:
        structure: StructureData for atomic positions.
        velocities_dict: Dict from extract_velocities_from_contcar with velocity data.

    Returns:
        SinglefileData containing POSCAR with velocity block (if velocities available).

    Example:
        >>> poscar_file = create_poscar_file_with_velocities(
        ...     structure=prev_structure,
        ...     velocities_dict=prev_velocities,
        ... )
        >>> calc = AimdVaspCalculation(poscar_file=poscar_file, ...)
    """
    from .utils import build_poscar_with_velocities
    import tempfile
    import os

    vel_dict = velocities_dict.get_dict()
    has_vels = vel_dict.get('has_velocities', False)
    vels_data = vel_dict.get('velocities', [])

    import numpy as np
    velocities = None
    if has_vels and len(vels_data) == len(structure.sites):
        velocities = np.array(vels_data, dtype=float)

    poscar_str = build_poscar_with_velocities(structure, velocities=velocities)

    # Create SinglefileData from POSCAR content
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vasp', delete=False) as f:
        f.write(poscar_str)
        temp_path = f.name

    try:
        result = orm.SinglefileData(file=temp_path, filename='POSCAR')
    finally:
        os.unlink(temp_path)

    return result


@task.calcfunction
def extract_qe_energy(output_parameters: orm.Dict) -> orm.Float:
    """
    Extract total energy from QE output_parameters.

    Args:
        output_parameters: QE output_parameters Dict containing parsed output

    Returns:
        Total energy as Float (eV)

    Raises:
        ValueError: If 'energy' key not found in output_parameters
    """
    params = output_parameters.get_dict()

    if 'energy' in params:
        return orm.Float(float(params['energy']))

    # If no 'energy' key found, raise error with available keys
    available = ', '.join(sorted(params.keys()))
    raise ValueError(f'No energy in QE output. Available keys: {available}')


# Conversion constant: 1 Hartree = 27.211386245988 eV
HARTREE_TO_EV = 27.211386245988


@task.calcfunction
def extract_cp2k_energy(output_parameters: orm.Dict) -> orm.Float:
    """
    Extract total energy from CP2K output_parameters and convert to eV.

    CP2K outputs energy in Hartree. This function converts to eV for
    consistency with VASP and QE energy extractors.

    Args:
        output_parameters: CP2K output_parameters Dict containing parsed output

    Returns:
        Total energy as Float (eV)

    Raises:
        ValueError: If 'energy' key not found in output_parameters
    """
    params = output_parameters.get_dict()

    if 'energy' in params:
        energy_hartree = float(params['energy'])
        energy_ev = energy_hartree * HARTREE_TO_EV
        return orm.Float(energy_ev)

    # If no 'energy' key found, raise error with available keys
    available = ', '.join(sorted(params.keys()))
    raise ValueError(f'No energy in CP2K output. Available keys: {available}')


@task.calcfunction
def compute_cp2k_dynamics(
    structure: orm.StructureData,
    fix_type: orm.Str,
    fix_thickness: orm.Float,
    fix_elements: orm.List = None,
    fix_components: orm.Str = None,
) -> orm.Dict:
    """
    Compute fixed atoms specification for CP2K MOTION.CONSTRAINT.FIXED_ATOMS.

    This calcfunction is used when the structure is not known at build time
    (e.g., from a previous stage's output), so we need to compute the fixed
    atoms at runtime.

    Args:
        structure: Input structure to analyze
        fix_type: Where to fix atoms ('bottom', 'center', 'top')
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to restrict fixing to
        fix_components: Components to fix ('XYZ', 'XY', 'Z', etc.)

    Returns:
        Dict with CP2K FIXED_ATOMS section ready for merging into parameters
    """
    from teros.core.fixed_atoms import get_fixed_atoms_list

    fix_type_val = fix_type.value
    fix_thickness_val = fix_thickness.value
    fix_elements_val = fix_elements.get_list() if fix_elements is not None else None
    fix_components_val = fix_components.value if fix_components is not None else 'XYZ'

    # Get list of atom indices to fix (1-based, which CP2K uses)
    fixed_atoms_list = get_fixed_atoms_list(
        structure=structure,
        fix_type=fix_type_val,
        fix_thickness=fix_thickness_val,
        fix_elements=fix_elements_val,
    )

    if not fixed_atoms_list:
        # No atoms to fix
        return orm.Dict(dict={'fixed_atoms': None})

    # Create FIXED_ATOMS section for CP2K
    # CP2K expects LIST as space-separated string of 1-based indices
    fixed_atoms_section = {
        'LIST': ' '.join(map(str, fixed_atoms_list)),
        'COMPONENTS_TO_FIX': fix_components_val,
    }

    return orm.Dict(dict={'fixed_atoms': fixed_atoms_section})


@task.calcfunction
def concatenate_trajectories(
    trajectories: t.Annotated[dict, dynamic(orm.TrajectoryData)],
) -> orm.TrajectoryData:
    """
    Concatenate multiple AIMD trajectory files from sequential stages.

    This calcfunction merges trajectory data from multiple AIMD stages
    into a single TrajectoryData object, preserving frame order based
    on sorted stage namespace keys.

    Args:
        trajectories: Dict of {stage_namespace: TrajectoryData}.
                     Keys should be sortable (e.g., 's01_...', 's02_...').

    Returns:
        TrajectoryData with all frames from all stages concatenated.

    Note:
        Arrays are concatenated along axis=0 (frame dimension).
        Symbols array is taken from the first trajectory (assumes same atoms).
    """
    import numpy as np

    all_positions = []
    all_cells = []
    all_energies = []
    all_forces = []
    symbols = None

    # Process in sorted order for correct frame sequence
    for stage_key in sorted(trajectories.keys()):
        traj = trajectories[stage_key]

        try:
            all_positions.append(traj.get_array('positions'))
        except KeyError:
            pass

        try:
            all_cells.append(traj.get_array('cells'))
        except KeyError:
            pass

        try:
            all_energies.append(traj.get_array('energies'))
        except KeyError:
            pass

        try:
            all_forces.append(traj.get_array('forces'))
        except KeyError:
            pass

        if symbols is None:
            try:
                symbols = traj.get_array('symbols')
            except KeyError:
                pass

    result = orm.TrajectoryData()

    if all_positions:
        result.set_array('positions', np.concatenate(all_positions, axis=0))
    if all_cells:
        result.set_array('cells', np.concatenate(all_cells, axis=0))
    if all_energies:
        result.set_array('energies', np.concatenate(all_energies, axis=0))
    if all_forces:
        result.set_array('forces', np.concatenate(all_forces, axis=0))
    if symbols is not None:
        result.set_array('symbols', symbols)

    return result
