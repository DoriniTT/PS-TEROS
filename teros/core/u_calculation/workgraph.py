"""WorkGraph builder for Hubbard U calculation using VASP linear response method.

This module implements the linear response approach to calculate the Hubbard U
parameter for LSDA+U calculations. The workflow follows the VASP wiki method:

1. Ground State: DFT calculation without +U, saves CHGCAR/WAVECAR
2. Non-SCF Response: With LDAUTYPE=3 and ICHARG=11 (fixed charge)
3. SCF Response: With LDAUTYPE=3, no ICHARG (charge evolves)

Steps 2-3 are repeated for multiple potential values, and U is calculated
from linear regression of the occupation changes.

Reference: https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U
"""

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task

from ..utils import get_vasp_parser_settings
from .tasks import (
    extract_d_electron_occupation,
    calculate_occupation_response,
    calculate_hubbard_u_linear_regression,
    gather_responses,
)
from .utils import (
    prepare_ground_state_incar,
    prepare_response_incar,
    validate_target_species,
    get_species_order_from_structure,
    DEFAULT_POTENTIAL_VALUES,
)


def build_u_calculation_workgraph(
    structure: orm.StructureData,
    code_label: str,
    potential_family: str,
    potential_mapping: dict,
    target_species: str,
    potential_values: t.Optional[t.List[float]] = None,
    ldaul: int = 2,
    ldauj: float = 0.0,
    ground_state_parameters: t.Optional[dict] = None,
    response_parameters: t.Optional[dict] = None,
    options: t.Optional[dict] = None,
    kpoints_spacing: float = 0.03,
    clean_workdir: bool = False,
    name: str = 'HubbardUCalculation',
) -> WorkGraph:
    """
    Build a WorkGraph to calculate the Hubbard U parameter using linear response.

    The workflow performs:
    1. Ground state calculation (no +U) to establish baseline d-occupancy
    2. For each potential value V in potential_values:
       - Non-SCF response (ICHARG=11): Fixed charge, apply potential
       - SCF response: Allow charge to evolve
    3. Linear regression to calculate U from multiple responses

    Args:
        structure: AiiDA StructureData for the material
        code_label: VASP code label (e.g., 'VASP-6.5.1@cluster')
        potential_family: POTCAR family name (e.g., 'PBE.54')
        potential_mapping: Element to potential mapping (e.g., {'Fe': 'Fe', 'O': 'O'})
        target_species: Element symbol for Hubbard U (e.g., 'Fe', 'Ni', 'Mn')
        potential_values: List of potentials to apply (eV). Default: [-0.2, -0.1, 0.0, 0.1, 0.2]
        ldaul: Angular momentum (2=d electrons, 3=f electrons). Default: 2
        ldauj: Exchange J parameter. Default: 0.0
        ground_state_parameters: Base INCAR for ground state (ENCUT, EDIFF, etc.)
        response_parameters: Base INCAR for response calculations (can override)
        options: Scheduler options (resources, queue_name, etc.)
        kpoints_spacing: K-point spacing in 1/Angstrom. Default: 0.03
        clean_workdir: Whether to clean remote work directories. Default: False
        name: WorkGraph name. Default: 'HubbardUCalculation'

    Returns:
        WorkGraph ready to submit

    Example:
        >>> from aiida import orm
        >>> structure = orm.load_node(123)  # Your NiO structure
        >>> wg = build_u_calculation_workgraph(
        ...     structure=structure,
        ...     code_label='VASP-6.5.1@cluster',
        ...     potential_family='PBE.54',
        ...     potential_mapping={'Ni': 'Ni', 'O': 'O'},
        ...     target_species='Ni',
        ...     ground_state_parameters={'ENCUT': 520, 'EDIFF': 1e-6, 'ISMEAR': 0, 'SIGMA': 0.05},
        ...     options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        ... )
        >>> wg.submit()
    """
    # Validate inputs
    validate_target_species(structure, target_species)

    # Set defaults
    if potential_values is None:
        potential_values = DEFAULT_POTENTIAL_VALUES

    if ground_state_parameters is None:
        ground_state_parameters = {
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'PREC': 'Accurate',
            'ALGO': 'Normal',
            'NELM': 100,
        }

    if response_parameters is None:
        response_parameters = ground_state_parameters.copy()

    if options is None:
        options = {
            'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
            'max_wallclock_seconds': 86400,
        }

    # Get species order for LDAU arrays
    all_species = get_species_order_from_structure(structure)

    # Load VASP code and wrap as task
    code = orm.load_code(code_label)
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Get parser settings that request orbital data
    settings = get_vasp_parser_settings(
        add_energy=True,
        add_trajectory=True,
        add_structure=True,
        add_kpoints=True,
    )

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # =========================================================================
    # STEP 1: Ground State Calculation
    # =========================================================================
    gs_incar = prepare_ground_state_incar(
        base_params=ground_state_parameters,
        lmaxmix=4 if ldaul == 2 else 6,  # 4 for d, 6 for f electrons
    )

    ground_state = wg.add_task(
        VaspTask,
        name='ground_state',
        structure=structure,
        code=code,
        parameters=orm.Dict(dict={'incar': gs_incar}),
        options=orm.Dict(dict=options),
        potential_family=potential_family,
        potential_mapping=orm.Dict(dict=potential_mapping),
        kpoints_spacing=kpoints_spacing,
        clean_workdir=False,  # MUST keep CHGCAR/WAVECAR
        settings=orm.Dict(dict=settings),
    )

    # Extract ground state d-electron occupation
    gs_occupation = wg.add_task(
        extract_d_electron_occupation,
        name='gs_occupation',
        misc=ground_state.outputs.misc,
        target_species=orm.Str(target_species),
        structure=structure,
    )

    # =========================================================================
    # STEP 2-3: Response Calculations for Each Potential Value
    # =========================================================================
    response_tasks = {}

    for i, V in enumerate(potential_values):
        label = f'V_{i}'
        V_str = f'{V:+.2f}'.replace('.', 'p').replace('-', 'm').replace('+', 'p')

        # ----- Non-SCF Response (ICHARG=11) -----
        nscf_incar = prepare_response_incar(
            base_params=response_parameters,
            potential_value=V,
            target_species=target_species,
            all_species=all_species,
            ldaul=ldaul,
            ldauj=ldauj,
            is_scf=False,  # ICHARG=11
            lmaxmix=4 if ldaul == 2 else 6,
        )

        nscf_task = wg.add_task(
            VaspTask,
            name=f'nscf_{V_str}',
            structure=structure,
            code=code,
            parameters=orm.Dict(dict={'incar': nscf_incar}),
            options=orm.Dict(dict=options),
            potential_family=potential_family,
            potential_mapping=orm.Dict(dict=potential_mapping),
            kpoints_spacing=kpoints_spacing,
            restart_folder=ground_state.outputs.remote_folder,
            clean_workdir=clean_workdir,
            settings=orm.Dict(dict=settings),
        )

        # Extract NSCF d-occupation
        nscf_occ = wg.add_task(
            extract_d_electron_occupation,
            name=f'nscf_occ_{V_str}',
            misc=nscf_task.outputs.misc,
            target_species=orm.Str(target_species),
            structure=structure,
        )

        # ----- SCF Response (no ICHARG) -----
        scf_incar = prepare_response_incar(
            base_params=response_parameters,
            potential_value=V,
            target_species=target_species,
            all_species=all_species,
            ldaul=ldaul,
            ldauj=ldauj,
            is_scf=True,  # No ICHARG
            lmaxmix=4 if ldaul == 2 else 6,
        )

        scf_task = wg.add_task(
            VaspTask,
            name=f'scf_{V_str}',
            structure=structure,
            code=code,
            parameters=orm.Dict(dict={'incar': scf_incar}),
            options=orm.Dict(dict=options),
            potential_family=potential_family,
            potential_mapping=orm.Dict(dict=potential_mapping),
            kpoints_spacing=kpoints_spacing,
            restart_folder=ground_state.outputs.remote_folder,
            clean_workdir=clean_workdir,
            settings=orm.Dict(dict=settings),
        )

        # Extract SCF d-occupation
        scf_occ = wg.add_task(
            extract_d_electron_occupation,
            name=f'scf_occ_{V_str}',
            misc=scf_task.outputs.misc,
            target_species=orm.Str(target_species),
            structure=structure,
        )

        # Calculate response for this potential
        response = wg.add_task(
            calculate_occupation_response,
            name=f'response_{V_str}',
            ground_state_occupation=gs_occupation.outputs.result,
            nscf_occupation=nscf_occ.outputs.result,
            scf_occupation=scf_occ.outputs.result,
            potential_value=orm.Float(V),
        )

        response_tasks[label] = response

    # =========================================================================
    # STEP 4: Gather Responses and Calculate U
    # =========================================================================
    # Gather all responses into a list
    gather_kwargs = {label: task.outputs.result for label, task in response_tasks.items()}
    gathered = wg.add_task(
        gather_responses,
        name='gather_responses',
        **gather_kwargs,
    )

    # Calculate U via linear regression
    calc_u = wg.add_task(
        calculate_hubbard_u_linear_regression,
        name='calculate_u',
        responses=gathered.outputs.result,
    )

    # =========================================================================
    # Set WorkGraph Outputs
    # =========================================================================
    wg.add_output('any', 'hubbard_u_result', from_socket=calc_u.outputs.result)
    wg.add_output('any', 'ground_state_misc', from_socket=ground_state.outputs.misc)
    wg.add_output('any', 'ground_state_remote', from_socket=ground_state.outputs.remote_folder)

    return wg


def get_u_calculation_results(workgraph) -> dict:
    """
    Extract results from a completed U calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_u_calculation_workgraph

    Returns:
        dict with:
            - U: Calculated Hubbard U value (eV)
            - chi_slope: SCF response slope
            - chi_0_slope: NSCF response slope
            - chi_r2: SCF fit R²
            - chi_0_r2: NSCF fit R²
            - potential_values: Potentials used
            - delta_n_scf_values: SCF occupation changes
            - delta_n_nscf_values: NSCF occupation changes

    Example:
        >>> results = get_u_calculation_results(completed_wg)
        >>> print(f"Hubbard U = {results['U']:.2f} eV")
        >>> print(f"R² = {results['chi_r2']:.4f}")
    """
    if hasattr(workgraph.outputs, 'hubbard_u_result'):
        return workgraph.outputs.hubbard_u_result.get_dict()
    else:
        # Try to get from task directly
        calc_u_task = workgraph.tasks.get('calculate_u')
        if calc_u_task and hasattr(calc_u_task.outputs, 'result'):
            return calc_u_task.outputs.result.value.get_dict()
        raise ValueError(
            "Could not find hubbard_u_result in workgraph outputs. "
            "Ensure the workgraph has completed successfully."
        )
