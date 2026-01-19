"""AiiDA calcfunction tasks for Fukui calculations.

These tasks extract energies, collect CHGCAR files, and generate summary
information from the parallel VASP charge-state calculations.
"""

import tempfile
import typing as t
from pathlib import Path

from aiida import orm
from aiida_workgraph import task

# Paths for aiida-shell integration (relative to this module)
_MODULE_DIR = Path(__file__).parent
WRAPPER_SCRIPT = str(_MODULE_DIR / 'scripts' / 'fukui_interpolation_wrapper.py')
ELECTRODES_WRAPPER_SCRIPT = str(_MODULE_DIR / 'scripts' / 'fukui_electrodes_wrapper.py')
PERTURBATIVE_WRAPPER_SCRIPT = str(_MODULE_DIR / 'scripts' / 'perturbative_expansion_wrapper.py')
FUKUI_GRID_PATH = str(_MODULE_DIR.parent.parent / 'external' / 'FukuiGrid')


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


# -----------------------------------------------------------------------------
# Internal versions for @task.graph (using positional args instead of **kwargs)
# -----------------------------------------------------------------------------

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

    **IMPORTANT**: This function supports a maximum of 4 delta_n values due to
    the fixed positional argument pattern required by @task.graph. For higher-order
    polynomial fitting requiring more points, use the external `collect_chgcar_files()`
    function instead, or run multiple workflows.

    Args:
        delta_n_list: List of delta_n values (e.g., [0.0, 0.05, 0.10, 0.15]).
                     **Maximum 4 values supported.**
        labels_list: List of labels corresponding to delta_n values
        retrieved_0: FolderData from first calculation (optional)
        retrieved_1: FolderData from second calculation (optional)
        retrieved_2: FolderData from third calculation (optional)
        retrieved_3: FolderData from fourth calculation (optional)

    Returns:
        FolderData containing CHGCAR files named CHGCAR_0.00, CHGCAR_0.05, etc.

    Raises:
        ValueError: If more than 4 delta_n values are provided
    """
    delta_n_list_py = delta_n_list.get_list()
    labels_list_py = labels_list.get_list()

    # Validate maximum number of delta_n values
    if len(delta_n_list_py) > 4:
        raise ValueError(
            f"collect_chgcar_files_internal supports a maximum of 4 delta_n values, "
            f"got {len(delta_n_list_py)}. For more points, use collect_chgcar_files() "
            f"with dynamic kwargs or run multiple workflows."
        )

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

    **IMPORTANT**: This function supports a maximum of 4 delta_n values due to
    the fixed positional argument pattern required by @task.graph.

    Args:
        nelect_neutral: Number of electrons in neutral system
        delta_n_list: List of delta_n values used. **Maximum 4 values supported.**
        fukui_type: Type of Fukui calculation ('plus' or 'minus')
        misc_0: Dict from first VASP calculation (optional)
        misc_1: Dict from second VASP calculation (optional)
        misc_2: Dict from third VASP calculation (optional)
        misc_3: Dict from fourth VASP calculation (optional)

    Returns:
        Dict with summary information

    Raises:
        ValueError: If more than 4 delta_n values are provided
    """
    from .utils import make_delta_label

    nelect_neutral_val = nelect_neutral.value
    delta_n_list_py = delta_n_list.get_list()
    fukui_type_val = fukui_type.value

    # Validate maximum number of delta_n values
    if len(delta_n_list_py) > 4:
        raise ValueError(
            f"generate_fukui_summary_internal supports a maximum of 4 delta_n values, "
            f"got {len(delta_n_list_py)}. For more points, use generate_fukui_summary() "
            f"with dynamic kwargs."
        )

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


# -----------------------------------------------------------------------------
# FukuiGrid interpolation integration
# -----------------------------------------------------------------------------

@task.calcfunction
def run_fukui_interpolation_calcfunc(
    chgcar_files: orm.FolderData,
    delta_n_values: orm.List,
    fukui_type: orm.Str,
) -> orm.SinglefileData:
    """
    Run FukuiGrid interpolation to compute the Fukui function.

    This calcfunction runs the FukuiGrid wrapper script via subprocess to
    compute the Fukui function from CHGCAR files at different charge states.

    Args:
        chgcar_files: FolderData with CHGCAR_0.00, CHGCAR_0.05, etc.
        delta_n_values: List of delta_n values [0.0, 0.05, 0.10, 0.15]
        fukui_type: Str with 'plus' or 'minus'

    Returns:
        SinglefileData containing CHGCAR_FUKUI.vasp
    """
    import subprocess
    import sys
    import os

    # Get Python values from AiiDA types
    dn_list = delta_n_values.get_list()
    ftype = fukui_type.value

    # Prepare file names (sorted by delta_n)
    sorted_dn = sorted(dn_list)
    file_names = [f'CHGCAR_{dn:.2f}' for dn in sorted_dn]
    delta_n_csv = ','.join(str(dn) for dn in sorted_dn)

    # Create temporary directory for processing
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Copy CHGCAR files to temp directory
        for filename in file_names:
            content = chgcar_files.get_object_content(filename, mode='rb')
            filepath = tmppath / filename
            filepath.write_bytes(content)

        # Build command arguments
        args = [
            sys.executable,
            WRAPPER_SCRIPT,
            ftype,
            delta_n_csv,
        ] + file_names

        # Run the wrapper script
        result = subprocess.run(
            args,
            cwd=str(tmppath),
            capture_output=True,
            text=True,
        )

        # Check for errors
        if result.returncode != 0:
            raise RuntimeError(
                f"FukuiGrid interpolation failed:\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Read output file
        output_path = tmppath / 'CHGCAR_FUKUI.vasp'
        if not output_path.exists():
            raise FileNotFoundError(
                f"FukuiGrid did not produce CHGCAR_FUKUI.vasp.\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Create SinglefileData from the output
        output_file = orm.SinglefileData(str(output_path), filename='CHGCAR_FUKUI.vasp')

    return output_file


def run_fukui_interpolation(
    chgcar_files: orm.FolderData,
    delta_n_values: t.List[float],
    fukui_type: str = 'plus',
) -> orm.SinglefileData:
    """
    Run FukuiGrid interpolation using aiida-shell (standalone function).

    This function launches a shell job that runs the FukuiGrid interpolation
    to compute the Fukui function from CHGCAR files at different charge states.

    Note: For use within WorkGraphs, use `run_fukui_interpolation_calcfunc` instead.

    Args:
        chgcar_files: FolderData with CHGCAR_0.00, CHGCAR_0.05, etc.
        delta_n_values: List of delta_n values [0.0, 0.05, 0.10, 0.15]
        fukui_type: 'plus' (nucleophilic f+) or 'minus' (electrophilic f-)

    Returns:
        SinglefileData containing CHGCAR_FUKUI.vasp

    Example:
        >>> from aiida import orm
        >>> chgcar_folder = orm.load_node(31972)
        >>> result = run_fukui_interpolation(
        ...     chgcar_folder,
        ...     delta_n_values=[0.0, 0.05, 0.10, 0.15],
        ...     fukui_type='plus',
        ... )
        >>> print(result.filename)
        'CHGCAR_FUKUI.vasp'
    """
    import sys
    from aiida_shell import launch_shell_job

    # Use the same Python interpreter that's running AiiDA
    python_path = sys.executable

    # Prepare file names (sorted by delta_n)
    sorted_dn = sorted(delta_n_values)
    file_names = [f'CHGCAR_{dn:.2f}' for dn in sorted_dn]
    delta_n_csv = ','.join(str(dn) for dn in sorted_dn)

    # Build arguments: wrapper.py fukui_type delta_n_csv file1 file2 file3 file4
    arguments = f'{WRAPPER_SCRIPT} {fukui_type} {delta_n_csv} ' + ' '.join(file_names)

    # Launch shell job
    # Note: aiida-shell converts dots in output names to underscores
    # So 'CHGCAR_FUKUI.vasp' becomes 'CHGCAR_FUKUI_vasp' in results
    results, node = launch_shell_job(
        python_path,
        arguments=arguments,
        nodes={'folder': chgcar_files},
        filenames={'folder': '.'},  # Extract folder contents to working directory root
        outputs=['CHGCAR_FUKUI.vasp'],
    )

    # Access the output using the underscore version (aiida-shell converts dots to underscores)
    return results['CHGCAR_FUKUI_vasp']


# -----------------------------------------------------------------------------
# Phase 2: Dielectric constant extraction and Fukui potential via electrodes
# -----------------------------------------------------------------------------

@task.calcfunction
def extract_dielectric_constant(
    retrieved: orm.FolderData,
    method: orm.Str = None,
) -> orm.Float:
    """
    Extract dielectric constant from VASP LEPSILON calculation.

    Parses the OUTCAR to get the macroscopic dielectric tensor and computes
    the scalar value using the specified method (default: arithmetic mean).

    Args:
        retrieved: FolderData from VASP LEPSILON calculation (must include OUTCAR)
        method: 'arithmetic' (default) or 'geometric'
                - arithmetic: ε = (εxx + εyy + εzz) / 3
                - geometric: ε = (εxx · εyy · εzz)^(1/3)

    Returns:
        Float with scalar dielectric constant

    Raises:
        FileNotFoundError: If OUTCAR not found in retrieved folder
        ValueError: If dielectric tensor not found in OUTCAR
    """
    from pymatgen.io.vasp import Outcar

    # Default to arithmetic mean
    method_val = method.value if method is not None else 'arithmetic'

    # Read OUTCAR from retrieved folder
    try:
        outcar_content = retrieved.get_object_content('OUTCAR')
    except (FileNotFoundError, OSError) as e:
        raise FileNotFoundError(
            f"OUTCAR not found in retrieved folder. "
            f"Ensure 'OUTCAR' is in ADDITIONAL_RETRIEVE_LIST for DFPT calculation. "
            f"Error: {e}"
        )

    # Write to temp file for pymatgen parsing
    with tempfile.NamedTemporaryFile(mode='w', suffix='_OUTCAR', delete=False) as f:
        if isinstance(outcar_content, bytes):
            f.write(outcar_content.decode('utf-8'))
        else:
            f.write(outcar_content)
        outcar_path = f.name

    try:
        outcar = Outcar(outcar_path)
        outcar.read_lepsilon()

        # Check if dielectric tensor was parsed
        if not hasattr(outcar, 'dielectric_tensor') or outcar.dielectric_tensor is None:
            raise ValueError(
                "Dielectric tensor not found in OUTCAR. "
                "Ensure the VASP calculation used LEPSILON=.TRUE."
            )

        tensor = outcar.dielectric_tensor  # 3x3 matrix

        eps_xx = tensor[0][0]
        eps_yy = tensor[1][1]
        eps_zz = tensor[2][2]

        if method_val == 'geometric':
            epsilon = (eps_xx * eps_yy * eps_zz) ** (1 / 3)
        else:  # arithmetic (default)
            epsilon = (eps_xx + eps_yy + eps_zz) / 3

        return orm.Float(epsilon)
    finally:
        Path(outcar_path).unlink()


@task.calcfunction
def run_fukui_electrodes_calcfunc(
    chgcar_files: orm.FolderData,
    fukui_chgcar: orm.SinglefileData,
    epsilon: orm.Float,
) -> orm.SinglefileData:
    """
    Run FukuiGrid electrodes method to compute Fukui potential.

    The electrodes method applies corrections to the electrostatic potential
    to calculate the Fukui potential for electrode/surface systems. This
    accounts for periodic boundary conditions and dielectric screening.

    Args:
        chgcar_files: FolderData containing CHGCAR_0.00 (neutral charge density)
        fukui_chgcar: SinglefileData with CHGCAR_FUKUI.vasp from Phase 1 interpolation
        epsilon: Float with dielectric constant (from DFPT calculation)

    Returns:
        SinglefileData containing LOCPOT_FUKUI.vasp (Fukui potential)

    Raises:
        FileNotFoundError: If required input files are missing
        RuntimeError: If FukuiGrid electrodes calculation fails
    """
    import subprocess
    import sys

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Copy CHGCAR_0.00 (neutral charge density)
        try:
            neutral_content = chgcar_files.get_object_content('CHGCAR_0.00', mode='rb')
        except (FileNotFoundError, OSError):
            raise FileNotFoundError(
                "CHGCAR_0.00 not found in chgcar_files. "
                "This file should be the neutral charge density from delta_n=0.0 calculation."
            )
        neutral_path = tmppath / 'CHGCAR_0.00'
        neutral_path.write_bytes(neutral_content)

        # Copy CHGCAR_FUKUI.vasp (Fukui function from Phase 1)
        try:
            fukui_content = fukui_chgcar.get_content(mode='rb')
        except (FileNotFoundError, OSError, IOError) as e:
            raise FileNotFoundError(
                f"Could not read Fukui function file: {e}. "
                f"Ensure Phase 1 interpolation completed successfully."
            )
        fukui_path = tmppath / 'CHGCAR_FUKUI.vasp'
        fukui_path.write_bytes(fukui_content)

        # Run the electrodes wrapper script
        result = subprocess.run(
            [
                sys.executable,
                ELECTRODES_WRAPPER_SCRIPT,
                'CHGCAR_0.00',
                'CHGCAR_FUKUI.vasp',
                str(epsilon.value),
            ],
            cwd=str(tmppath),
            capture_output=True,
            text=True,
        )

        # Check for errors
        if result.returncode != 0:
            raise RuntimeError(
                f"FukuiGrid electrodes calculation failed:\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Read output file (FukuiGrid writes LOCPOT_FUKUI.vasp)
        output_path = tmppath / 'LOCPOT_FUKUI.vasp'
        if not output_path.exists():
            raise FileNotFoundError(
                f"FukuiGrid did not produce LOCPOT_FUKUI.vasp.\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Create SinglefileData from the output
        output_file = orm.SinglefileData(str(output_path), filename='LOCPOT_FUKUI.vasp')

    return output_file


# -----------------------------------------------------------------------------
# Phase 4: Perturbative expansion model for interaction energy prediction
# -----------------------------------------------------------------------------

@task.calcfunction
def extract_locpot_from_retrieved(
    retrieved: orm.FolderData,
) -> orm.SinglefileData:
    """
    Extract LOCPOT file from VASP retrieved folder.

    Used in Phase 4 to extract the electrostatic potential from the neutral
    slab calculation (delta_n=0.0), which is needed for the perturbative
    expansion model.

    Args:
        retrieved: FolderData from VASP calculation containing LOCPOT

    Returns:
        SinglefileData containing LOCPOT

    Raises:
        FileNotFoundError: If LOCPOT not found in retrieved folder
    """
    try:
        locpot_content = retrieved.get_object_content('LOCPOT', mode='rb')
    except (FileNotFoundError, OSError) as e:
        raise FileNotFoundError(
            f"LOCPOT not found in retrieved folder. "
            f"Ensure 'LOCPOT' is in ADDITIONAL_RETRIEVE_LIST. Error: {e}"
        )

    with tempfile.NamedTemporaryFile(delete=False, suffix='_LOCPOT') as f:
        f.write(locpot_content)
        temp_path = f.name

    try:
        return orm.SinglefileData(temp_path, filename='LOCPOT')
    finally:
        Path(temp_path).unlink(missing_ok=True)


@task.calcfunction
def run_perturbative_expansion_calcfunc(
    locpot_neutral: orm.SinglefileData,
    fukui_potential: orm.SinglefileData,
    probe_charge: orm.Float,
    electron_transfer: orm.Float,
) -> orm.SinglefileData:
    """
    Run perturbative expansion to compute interaction energy map.

    Uses the c-DFT perturbative expansion formula:
        ΔU(r) = q·Φ(r) - q·ΔN·vf±(r)

    Where:
        Φ(r) = electrostatic potential (LOCPOT from neutral slab)
        vf±(r) = Fukui potential (LOCPOT_FUKUI.vasp from Phase 2)
        q = charge of the probe (point charge model)
        ΔN = electron transfer

    The output MODELPOT_LOCPOT.vasp contains a 3D potential field showing
    favorable (negative) and unfavorable (positive) adsorption sites.

    Args:
        locpot_neutral: LOCPOT from neutral slab (electrostatic potential Φ)
        fukui_potential: LOCPOT_FUKUI.vasp from Phase 2 (Fukui potential vf±)
        probe_charge: Charge q of the probe in |e| (point charge model)
        electron_transfer: Electron transfer ΔN (positive = electron gain,
                          negative = electron donation)

    Returns:
        SinglefileData containing MODELPOT_LOCPOT.vasp

    Raises:
        RuntimeError: If perturbative expansion calculation fails
        FileNotFoundError: If output file not generated
    """
    import subprocess
    import sys

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Copy LOCPOT (neutral electrostatic potential)
        locpot_content = locpot_neutral.get_content(mode='rb')
        locpot_path = tmppath / 'LOCPOT'
        locpot_path.write_bytes(locpot_content)

        # Copy LOCPOT_FUKUI.vasp (Fukui potential from Phase 2)
        fukui_content = fukui_potential.get_content(mode='rb')
        fukui_path = tmppath / 'LOCPOT_FUKUI.vasp'
        fukui_path.write_bytes(fukui_content)

        # Run wrapper script
        result = subprocess.run(
            [
                sys.executable,
                PERTURBATIVE_WRAPPER_SCRIPT,
                'LOCPOT',
                'LOCPOT_FUKUI.vasp',
                str(probe_charge.value),
                str(electron_transfer.value),
            ],
            cwd=str(tmppath),
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"Perturbative expansion calculation failed:\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Read output file
        output_path = tmppath / 'MODELPOT_LOCPOT.vasp'
        if not output_path.exists():
            raise FileNotFoundError(
                f"MODELPOT_LOCPOT.vasp not generated.\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        return orm.SinglefileData(str(output_path), filename='MODELPOT_LOCPOT.vasp')
