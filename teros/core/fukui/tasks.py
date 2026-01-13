"""AiiDA calcfunction tasks for Fukui calculations.

These tasks extract energies, collect CHGCAR files, and generate summary
information from the parallel VASP charge-state calculations.
"""

import tempfile
import typing as t
from pathlib import Path

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
def collect_chgcar_files(
    delta_n_list: orm.List,
    labels_list: orm.List,
    **retrieved_kwargs,
) -> orm.FolderData:
    """
    Collect CHGCAR files from multiple VASP calculations into a single FolderData.

    Reads CHGCAR from each retrieved FolderData, renames them to CHGCAR_X.XX format
    (where X.XX is the delta_n value), and writes to a new consolidated FolderData.

    Args:
        delta_n_list: List of delta_n values (e.g., [0.0, 0.05, 0.10, 0.15])
        labels_list: List of labels corresponding to kwargs keys
        **retrieved_kwargs: Dynamic kwargs with FolderData from each calculation
                           Keys are labels like 'delta_0_00', 'delta_0_05', etc.

    Returns:
        FolderData containing CHGCAR files named CHGCAR_0.00, CHGCAR_0.05, etc.

    Raises:
        ValueError: If CHGCAR not found in any retrieved folder
    """
    delta_n_list_py = delta_n_list.get_list()
    labels_list_py = labels_list.get_list()

    # Verify all CHGCARs exist before copying
    missing = []
    for label in labels_list_py:
        if label not in retrieved_kwargs:
            missing.append(f"{label} (not in kwargs)")
            continue
        retrieved = retrieved_kwargs[label]
        try:
            retrieved.get_object_content('CHGCAR')
        except (FileNotFoundError, OSError):
            missing.append(label)

    if missing:
        raise ValueError(
            f"CHGCAR not found in retrieved folders for: {missing}. "
            f"Ensure 'CHGCAR' is in ADDITIONAL_RETRIEVE_LIST settings."
        )

    # Create FolderData with renamed CHGCARs
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        for i, label in enumerate(labels_list_py):
            retrieved = retrieved_kwargs[label]
            delta_n = delta_n_list_py[i]

            # Read CHGCAR from retrieved folder
            chgcar_content = retrieved.get_object_content('CHGCAR')

            # Write with renamed name: CHGCAR_X.XX
            output_name = f"CHGCAR_{delta_n:.2f}"
            output_path = tmppath / output_name

            # Handle both string and bytes content
            if isinstance(chgcar_content, str):
                output_path.write_text(chgcar_content)
            else:
                output_path.write_bytes(chgcar_content)

        # Create FolderData from directory
        folder = orm.FolderData(tree=str(tmppath))

    return folder


@task.calcfunction
def generate_fukui_summary(
    nelect_neutral: orm.Int,
    delta_n_list: orm.List,
    fukui_type: orm.Str,
    **misc_kwargs,
) -> orm.Dict:
    """
    Generate summary dict with calculation metadata.

    Args:
        nelect_neutral: Number of electrons in neutral system
        delta_n_list: List of delta_n values used
        fukui_type: Type of Fukui calculation ('plus' or 'minus')
        **misc_kwargs: Dynamic kwargs with Dict from each VASP calculation
                      Keys are labels like 'delta_0_00', 'delta_0_05', etc.

    Returns:
        Dict with summary information including:
        - nelect_neutral: Base electron count
        - delta_n_values: List of dN values used
        - fukui_type: 'plus' or 'minus'
        - calculations: {label: {nelect, energy, converged}}
        - formula: Description of Fukui function formula used
    """
    nelect_neutral_val = nelect_neutral.value
    delta_n_list_py = delta_n_list.get_list()
    fukui_type_val = fukui_type.value

    # Calculate NELECT for each delta_n
    calculations = {}
    for i, delta_n in enumerate(delta_n_list_py):
        from .utils import make_delta_label
        label = make_delta_label(delta_n)

        if fukui_type_val == 'plus':
            # Fukui+: remove electrons (nucleophilic attack)
            nelect = nelect_neutral_val - delta_n
        else:
            # Fukui-: add electrons (electrophilic attack)
            nelect = nelect_neutral_val + delta_n

        # Extract energy from misc if available
        energy = None
        converged = None
        if label in misc_kwargs:
            misc = misc_kwargs[label]
            if isinstance(misc, orm.Dict):
                misc_dict = misc.get_dict()
            else:
                misc_dict = dict(misc)

            # Get energy
            if 'total_energies' in misc_dict:
                energy_dict = misc_dict['total_energies']
                for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
                    if key in energy_dict:
                        energy = energy_dict[key]
                        break

            # Check convergence
            converged = misc_dict.get('run_status', {}).get('electronic_converged', None)

        calculations[label] = {
            'delta_n': delta_n,
            'nelect': nelect,
            'energy_eV': energy,
            'converged': converged,
        }

    # Build formula description
    if fukui_type_val == 'plus':
        formula = "f+(r) = [rho(N) - rho(N-dN)] / dN (nucleophilic attack - where electrons are added)"
    else:
        formula = "f-(r) = [rho(N+dN) - rho(N)] / dN (electrophilic attack - where electrons are removed)"

    summary = {
        'nelect_neutral': nelect_neutral_val,
        'delta_n_values': delta_n_list_py,
        'fukui_type': fukui_type_val,
        'formula': formula,
        'n_calculations': len(delta_n_list_py),
        'calculations': calculations,
    }

    return orm.Dict(dict=summary)


@task.calcfunction
def collect_chgcar_files_internal(
    delta_n_list: orm.List,
    labels_list: orm.List,
    retrieved_0: orm.FolderData = None,
    retrieved_1: orm.FolderData = None,
    retrieved_2: orm.FolderData = None,
    retrieved_3: orm.FolderData = None,
) -> orm.FolderData:
    """
    Collect CHGCAR files from VASP calculations into a single FolderData.

    This is an internal version using fixed positional arguments instead of
    dynamic **kwargs, to be used within @task.graph functions.

    Args:
        delta_n_list: List of delta_n values (e.g., [0.0, 0.05, 0.10, 0.15])
        labels_list: List of labels corresponding to delta_n values
        retrieved_0: FolderData from first calculation (optional)
        retrieved_1: FolderData from second calculation (optional)
        retrieved_2: FolderData from third calculation (optional)
        retrieved_3: FolderData from fourth calculation (optional)

    Returns:
        FolderData containing CHGCAR files named CHGCAR_0.00, CHGCAR_0.05, etc.
    """
    delta_n_list_py = delta_n_list.get_list()
    labels_list_py = labels_list.get_list()

    # Collect non-None retrieved folders in order
    retrieved_folders = [retrieved_0, retrieved_1, retrieved_2, retrieved_3]

    # Verify we have enough folders for the delta_n values
    available_folders = [f for f in retrieved_folders[:len(delta_n_list_py)] if f is not None]
    if len(available_folders) < len(delta_n_list_py):
        raise ValueError(
            f"Not enough retrieved folders. Expected {len(delta_n_list_py)}, "
            f"got {len(available_folders)}"
        )

    # Verify all CHGCARs exist before copying
    missing = []
    for i, retrieved in enumerate(available_folders):
        try:
            retrieved.get_object_content('CHGCAR')
        except (FileNotFoundError, OSError):
            missing.append(labels_list_py[i])

    if missing:
        raise ValueError(
            f"CHGCAR not found in retrieved folders for: {missing}. "
            f"Ensure 'CHGCAR' is in ADDITIONAL_RETRIEVE_LIST settings."
        )

    # Create FolderData with renamed CHGCARs
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        for i, retrieved in enumerate(available_folders):
            delta_n = delta_n_list_py[i]

            # Read CHGCAR from retrieved folder
            chgcar_content = retrieved.get_object_content('CHGCAR')

            # Write with renamed name: CHGCAR_X.XX
            output_name = f"CHGCAR_{delta_n:.2f}"
            output_path = tmppath / output_name

            # Handle both string and bytes content
            if isinstance(chgcar_content, str):
                output_path.write_text(chgcar_content)
            else:
                output_path.write_bytes(chgcar_content)

        # Create FolderData from directory
        folder = orm.FolderData(tree=str(tmppath))

    return folder


@task.calcfunction
def generate_fukui_summary_internal(
    nelect_neutral: orm.Int,
    delta_n_list: orm.List,
    fukui_type: orm.Str,
    misc_0: orm.Dict = None,
    misc_1: orm.Dict = None,
    misc_2: orm.Dict = None,
    misc_3: orm.Dict = None,
) -> orm.Dict:
    """
    Generate summary dict with calculation metadata.

    This is an internal version using fixed positional arguments instead of
    dynamic **kwargs, to be used within @task.graph functions.

    Args:
        nelect_neutral: Number of electrons in neutral system
        delta_n_list: List of delta_n values used
        fukui_type: Type of Fukui calculation ('plus' or 'minus')
        misc_0: Dict from first VASP calculation (optional)
        misc_1: Dict from second VASP calculation (optional)
        misc_2: Dict from third VASP calculation (optional)
        misc_3: Dict from fourth VASP calculation (optional)

    Returns:
        Dict with summary information
    """
    from .utils import make_delta_label

    nelect_neutral_val = nelect_neutral.value
    delta_n_list_py = delta_n_list.get_list()
    fukui_type_val = fukui_type.value

    # Collect misc dicts in order
    misc_list = [misc_0, misc_1, misc_2, misc_3]

    # Calculate NELECT for each delta_n
    calculations = {}
    for i, delta_n in enumerate(delta_n_list_py):
        label = make_delta_label(delta_n)

        if fukui_type_val == 'plus':
            # Fukui+: remove electrons (nucleophilic attack)
            nelect = nelect_neutral_val - delta_n
        else:
            # Fukui-: add electrons (electrophilic attack)
            nelect = nelect_neutral_val + delta_n

        # Extract energy from misc if available
        energy = None
        converged = None
        if i < len(misc_list) and misc_list[i] is not None:
            misc = misc_list[i]
            if isinstance(misc, orm.Dict):
                misc_dict = misc.get_dict()
            else:
                misc_dict = dict(misc)

            # Get energy
            if 'total_energies' in misc_dict:
                energy_dict = misc_dict['total_energies']
                for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
                    if key in energy_dict:
                        energy = energy_dict[key]
                        break

            # Check convergence
            converged = misc_dict.get('run_status', {}).get('electronic_converged', None)

        calculations[label] = {
            'delta_n': delta_n,
            'nelect': nelect,
            'energy_eV': energy,
            'converged': converged,
        }

    # Build formula description
    if fukui_type_val == 'plus':
        formula = "f+(r) = [rho(N) - rho(N-dN)] / dN (nucleophilic attack - where electrons are added)"
    else:
        formula = "f-(r) = [rho(N+dN) - rho(N)] / dN (electrophilic attack - where electrons are removed)"

    summary = {
        'nelect_neutral': nelect_neutral_val,
        'delta_n_values': delta_n_list_py,
        'fukui_type': fukui_type_val,
        'formula': formula,
        'n_calculations': len(delta_n_list_py),
        'calculations': calculations,
    }

    return orm.Dict(dict=summary)
