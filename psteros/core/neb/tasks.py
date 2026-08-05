"""AiiDA calcfunction tasks for NEB calculations.

These tasks handle image interpolation, energy extraction, and barrier calculation
for Nudged Elastic Band (NEB) calculations.
"""

import typing as t

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def interpolate_structures(
    initial: orm.StructureData,
    final: orm.StructureData,
    n_images: orm.Int,
    method: orm.Str,
) -> orm.Dict:
    """
    Generate NEB intermediate images using ASE interpolation.

    Creates intermediate structures between initial and final configurations
    using either linear interpolation or IDPP (Image Dependent Pair Potential).

    Args:
        initial: Initial structure (reactant)
        final: Final structure (product)
        n_images: Number of intermediate images to generate
        method: Interpolation method - 'idpp' (recommended) or 'linear'

    Returns:
        Dict containing intermediate images as StructureData with keys
        'image_01', 'image_02', etc. The output dict has structure:
        {'image_01': StructureData, 'image_02': StructureData, ...}

    Example:
        >>> result = interpolate_structures(
        ...     initial=initial_structure,
        ...     final=final_structure,
        ...     n_images=orm.Int(5),
        ...     method=orm.Str('idpp'),
        ... )
        >>> # Access via: result.get_dict()['image_01'], etc.
    """
    from ase.mep import NEB
    import numpy as np

    # Convert to ASE Atoms objects
    atoms_initial = initial.get_ase()
    atoms_final = final.get_ase()

    n_img = n_images.value
    method_str = method.value

    # Create list of images (including endpoints)
    images = [atoms_initial.copy()]
    for _ in range(n_img):
        images.append(atoms_initial.copy())
    images.append(atoms_final.copy())

    # Set the final image positions
    images[-1].set_positions(atoms_final.get_positions())

    # Create NEB object and interpolate
    neb = NEB(images)

    if method_str == 'idpp':
        # IDPP interpolation - better initial path for transition states
        neb.interpolate('idpp', mic=True)
    else:
        # Linear interpolation with minimum image convention
        neb.interpolate(mic=True)

    # Extract intermediate images (exclude endpoints) and convert to StructureData
    result = {}
    for i, img in enumerate(images[1:-1], start=1):
        # Create StructureData from ASE Atoms
        structure = orm.StructureData(ase=img)
        result[f'image_{i:02d}'] = structure.store()

    return orm.Dict(dict={k: v.pk for k, v in result.items()})


@task.calcfunction
def create_single_neb_image(
    initial: orm.StructureData,
    final: orm.StructureData,
    n_images: orm.Int,
    image_index: orm.Int,
    method: orm.Str,
) -> orm.StructureData:
    """
    Generate a single NEB intermediate image.

    This calcfunction generates one specific intermediate image by performing
    the full interpolation and extracting the requested image. While this
    repeats the interpolation for each image, it allows proper AiiDA provenance.

    Args:
        initial: Initial structure (reactant)
        final: Final structure (product)
        n_images: Total number of intermediate images
        image_index: Index of the image to return (1-based)
        method: Interpolation method - 'idpp' (recommended) or 'linear'

    Returns:
        StructureData for the specified intermediate image
    """
    from ase.mep import NEB

    atoms_initial = initial.get_ase()
    atoms_final = final.get_ase()

    n_img = n_images.value
    idx = image_index.value
    method_str = method.value

    # Create list of images (including endpoints)
    images = [atoms_initial.copy()]
    for _ in range(n_img):
        images.append(atoms_initial.copy())
    images.append(atoms_final.copy())

    # Set the final image positions
    images[-1].set_positions(atoms_final.get_positions())

    # Create NEB object and interpolate
    neb = NEB(images)

    if method_str == 'idpp':
        neb.interpolate('idpp', mic=True)
    else:
        neb.interpolate(mic=True)

    # Return the specific intermediate image (idx is 1-based, images[0] is initial)
    return orm.StructureData(ase=images[idx])


@task.calcfunction
def extract_neb_energies(misc: orm.Dict) -> orm.Dict:
    """
    Extract energies from NEB calculation misc output.

    Parses the misc dictionary from a completed NEB calculation to extract
    the energy of each image along the reaction path.

    Args:
        misc: Dict output from VaspNEBWorkChain containing NEB results

    Returns:
        Dict with image indices as keys and energies (eV) as values:
        {'image_00': energy_0, 'image_01': energy_1, ...}
        Includes endpoints (image_00 = initial, image_N+1 = final)

    Raises:
        ValueError: If energy data cannot be extracted from misc
    """
    misc_dict = misc.get_dict()

    # Try to extract energies from different possible locations
    energies = {}

    # Check for per-image energies in NEB output
    if 'neb' in misc_dict:
        neb_data = misc_dict['neb']
        if 'energies' in neb_data:
            for i, energy in enumerate(neb_data['energies']):
                energies[f'image_{i:02d}'] = float(energy)
            return orm.Dict(dict=energies)

    # Check for total_energies structure (standard VASP output)
    if 'total_energies' in misc_dict:
        energy_dict = misc_dict['total_energies']
        # For NEB, this might be a list or dict per image
        if isinstance(energy_dict, list):
            for i, e_data in enumerate(energy_dict):
                if isinstance(e_data, dict):
                    energy = e_data.get(
                        'energy_extrapolated',
                        e_data.get('energy_no_entropy', e_data.get('energy'))
                    )
                else:
                    energy = float(e_data)
                energies[f'image_{i:02d}'] = float(energy)
            return orm.Dict(dict=energies)
        elif isinstance(energy_dict, dict):
            # Single energy - probably not NEB output
            energy = energy_dict.get(
                'energy_extrapolated',
                energy_dict.get('energy_no_entropy', energy_dict.get('energy'))
            )
            if energy is not None:
                energies['image_00'] = float(energy)
                return orm.Dict(dict=energies)

    # Try to find energies directly in misc_dict
    if 'energies' in misc_dict:
        for i, energy in enumerate(misc_dict['energies']):
            energies[f'image_{i:02d}'] = float(energy)
        return orm.Dict(dict=energies)

    # If we get here, we couldn't find energies
    available_keys = list(misc_dict.keys())
    raise ValueError(
        f"Could not extract NEB energies from misc output. "
        f"Available keys: {available_keys}"
    )


@task.calcfunction
def calculate_barrier(energies: orm.Dict) -> orm.Dict:
    """
    Calculate activation barriers from NEB energies.

    Computes forward and reverse activation barriers, reaction energy,
    and identifies the transition state image.

    Args:
        energies: Dict with image indices as keys and energies (eV) as values
                 Format: {'image_00': E_0, 'image_01': E_1, ...}

    Returns:
        Dict containing:
            - forward_barrier: E_TS - E_initial (eV)
            - reverse_barrier: E_TS - E_final (eV)
            - reaction_energy: E_final - E_initial (eV)
            - saddle_point_index: Index of the transition state image
            - saddle_point_energy: Energy at the saddle point (eV)
            - energies_list: List of all energies in order
            - is_exothermic: True if reaction_energy < 0

    Example:
        >>> barrier_results = calculate_barrier(neb_energies)
        >>> results = barrier_results.get_dict()
        >>> print(f"Forward barrier: {results['forward_barrier']:.3f} eV")
    """
    e_dict = energies.get_dict()

    # Sort by image index to ensure correct order
    sorted_keys = sorted(e_dict.keys())
    e_values = [e_dict[k] for k in sorted_keys]

    e_initial = e_values[0]
    e_final = e_values[-1]
    e_max = max(e_values)
    ts_index = e_values.index(e_max)

    forward_barrier = e_max - e_initial
    reverse_barrier = e_max - e_final
    reaction_energy = e_final - e_initial

    return orm.Dict(dict={
        'forward_barrier': forward_barrier,
        'reverse_barrier': reverse_barrier,
        'reaction_energy': reaction_energy,
        'saddle_point_index': ts_index,
        'saddle_point_energy': e_max,
        'energies_list': e_values,
        'n_images': len(e_values),
        'is_exothermic': reaction_energy < 0,
    })


@task.calcfunction
def extract_neb_trajectory(
    trajectory: orm.TrajectoryData,
) -> orm.ArrayData:
    """
    Extract trajectory data from NEB calculation for visualization.

    Parses the TrajectoryData from a completed NEB calculation to extract
    positions and other data useful for path visualization.

    Args:
        trajectory: TrajectoryData from VaspNEBWorkChain

    Returns:
        ArrayData containing:
            - positions: Atomic positions at each step (n_steps, n_atoms, 3)
            - cells: Cell vectors at each step (n_steps, 3, 3)
            - symbols: Chemical symbols (n_atoms,)
            - step_indices: Array of step indices

    Note:
        The returned ArrayData can be used for visualizing the NEB path
        or for further analysis of atomic movements during the transition.
    """
    import numpy as np

    # Get trajectory arrays
    try:
        positions = trajectory.get_array('positions')
    except KeyError:
        positions = None

    try:
        cells = trajectory.get_array('cells')
    except KeyError:
        cells = None

    try:
        symbols = trajectory.get_array('symbols')
    except KeyError:
        symbols = None

    # Create output ArrayData
    output = orm.ArrayData()

    if positions is not None:
        output.set_array('positions', positions)
        n_steps = positions.shape[0]
        output.set_array('step_indices', np.arange(n_steps))

    if cells is not None:
        output.set_array('cells', cells)

    if symbols is not None:
        output.set_array('symbols', symbols)

    return output


@task.calcfunction
def extract_total_energy(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP misc output.

    This is a general-purpose energy extraction function used for
    endpoint relaxation calculations.

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
def create_neb_summary(
    initial_energy: orm.Float,
    final_energy: orm.Float,
    barrier_results: orm.Dict,
    n_images: orm.Int,
    method: orm.Str,
    climb: orm.Bool,
) -> orm.Dict:
    """
    Create comprehensive summary of NEB calculation results.

    Combines all NEB results into a single summary Dict for easy access.

    Args:
        initial_energy: Energy of initial (relaxed) structure (eV)
        final_energy: Energy of final (relaxed) structure (eV)
        barrier_results: Dict from calculate_barrier with barrier values
        n_images: Number of intermediate images used
        method: Interpolation method used ('idpp' or 'linear')
        climb: Whether climbing image NEB was used

    Returns:
        Dict with complete NEB summary including:
            - All barrier information
            - Endpoint energies
            - Calculation parameters
            - Status information
    """
    barrier_dict = barrier_results.get_dict()

    summary = {
        # Barrier information
        'forward_barrier_eV': barrier_dict['forward_barrier'],
        'reverse_barrier_eV': barrier_dict['reverse_barrier'],
        'reaction_energy_eV': barrier_dict['reaction_energy'],
        'is_exothermic': barrier_dict['is_exothermic'],
        'saddle_point_index': barrier_dict['saddle_point_index'],
        'saddle_point_energy_eV': barrier_dict['saddle_point_energy'],

        # Endpoint energies
        'initial_energy_eV': initial_energy.value,
        'final_energy_eV': final_energy.value,

        # Energy profile
        'energies_list_eV': barrier_dict['energies_list'],
        'n_images_total': barrier_dict['n_images'],

        # Calculation parameters
        'n_intermediate_images': n_images.value,
        'interpolation_method': method.value,
        'climbing_image': climb.value,

        # Unit conversions (useful for analysis)
        'forward_barrier_kJ_mol': barrier_dict['forward_barrier'] * 96.485,  # eV to kJ/mol
        'reverse_barrier_kJ_mol': barrier_dict['reverse_barrier'] * 96.485,
        'reaction_energy_kJ_mol': barrier_dict['reaction_energy'] * 96.485,
    }

    return orm.Dict(dict=summary)
