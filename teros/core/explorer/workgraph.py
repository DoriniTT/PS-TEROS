"""WorkGraph builders for the explorer module.

This module provides lightweight, incremental VASP calculation wrappers
for exploratory work. The goal is simplicity: submit a calculation,
check results, decide next step, optionally restart.
"""

import typing as t
import time

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task

from .tasks import extract_energy
from .utils import prepare_restart_settings, get_status
from ..utils import deep_merge_dicts


def quick_vasp(
    structure: t.Union[orm.StructureData, int] = None,
    code_label: str = None,
    incar: dict = None,
    kpoints_spacing: float = 0.03,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    retrieve: t.List[str] = None,
    restart_from: int = None,
    copy_wavecar: bool = True,
    copy_chgcar: bool = False,
    name: str = 'quick_vasp',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> int:
    """
    Submit a single VASP calculation with minimal boilerplate.

    This is the primary entry point for exploratory VASP calculations.
    Specify INCAR parameters directly - no presets, maximum flexibility.

    Args:
        structure: StructureData or PK. If restart_from is provided, this is optional.
        code_label: VASP code label (e.g., 'VASP-6.5.1@localwork')
        incar: INCAR parameters dict (e.g., {'NSW': 100, 'IBRION': 2})
        kpoints_spacing: K-points spacing in A^-1 (default: 0.03)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping (e.g., {'Sn': 'Sn_d', 'O': 'O'})
        options: Scheduler options dict with 'resources' key
        retrieve: List of files to retrieve (e.g., ['CONTCAR', 'CHGCAR'])
        restart_from: PK of previous calculation to restart from
        copy_wavecar: Copy WAVECAR from restart_from (sets ISTART=1)
        copy_chgcar: Copy CHGCAR from restart_from (sets ICHARG=1)
        name: WorkGraph name for identification
        wait: If True, block until calculation finishes (default: False)
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean the work directory after completion

    Returns:
        PK of the submitted WorkGraph

    Example:
        >>> pk = quick_vasp(
        ...     structure=my_structure,
        ...     code_label='VASP-6.5.1@localwork',
        ...     incar={'NSW': 100, 'IBRION': 2, 'ISIF': 3},
        ...     kpoints_spacing=0.03,
        ...     potential_family='PBE',
        ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
        ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ...     retrieve=['CONTCAR', 'CHGCAR'],
        ...     name='sno2_relax',
        ... )

        >>> # Restart from previous calculation
        >>> pk2 = quick_vasp(
        ...     restart_from=pk,
        ...     code_label='VASP-6.5.1@localwork',
        ...     incar={'NSW': 0, 'NEDOS': 2000},
        ...     retrieve=['DOSCAR'],
        ...     name='sno2_dos',
        ... )
    """
    # Handle restart
    restart_folder = None
    incar_additions = {}

    if restart_from is not None:
        restart_structure, restart_settings = prepare_restart_settings(
            restart_from,
            copy_wavecar=copy_wavecar,
            copy_chgcar=copy_chgcar,
        )
        if structure is None:
            structure = restart_structure
        restart_folder = restart_settings['folder']
        incar_additions = restart_settings['incar_additions']

    # Validate required inputs
    if structure is None:
        raise ValueError("structure is required (or provide restart_from)")
    if code_label is None:
        raise ValueError("code_label is required")
    if incar is None:
        raise ValueError("incar is required - always specify INCAR parameters explicitly")
    if options is None:
        raise ValueError("options is required - specify scheduler resources")

    # Load structure if PK
    if isinstance(structure, int):
        structure = orm.load_node(structure)

    # Merge restart INCAR additions
    if incar_additions:
        incar = deep_merge_dicts(incar, incar_additions)

    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Build WorkGraph
    wg = WorkGraph(name=name)

    # Prepare builder inputs
    builder_inputs = _prepare_builder_inputs(
        incar=incar,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=potential_mapping or {},
        options=options,
        retrieve=retrieve,
        restart_folder=restart_folder,
        clean_workdir=clean_workdir,
    )

    # Add VASP task
    vasp_task = wg.add_task(
        VaspTask,
        name='vasp_calc',
        structure=structure,
        code=code,
        **builder_inputs
    )

    # Add energy extraction task
    energy_task = wg.add_task(
        extract_energy,
        name='extract_energy',
        misc=vasp_task.outputs.misc,
        retrieved=vasp_task.outputs.retrieved,
    )

    # Set WorkGraph outputs
    wg.outputs.energy = energy_task.outputs.result
    wg.outputs.structure = vasp_task.outputs.structure
    wg.outputs.misc = vasp_task.outputs.misc
    wg.outputs.retrieved = vasp_task.outputs.retrieved

    # Submit
    wg.submit()

    # Wait if requested
    if wait:
        _wait_for_completion(wg.pk, poll_interval)

    return wg.pk


def quick_vasp_batch(
    structures: t.Dict[str, t.Union[orm.StructureData, int]],
    code_label: str,
    incar: dict,
    kpoints_spacing: float = 0.03,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    retrieve: t.List[str] = None,
    incar_overrides: t.Dict[str, dict] = None,
    max_concurrent_jobs: int = None,
    name: str = 'quick_vasp_batch',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> t.Dict[str, int]:
    """
    Submit multiple VASP calculations with the same base settings.

    Useful for Fukui-style calculations or comparing different structures
    with consistent computational parameters.

    Args:
        structures: Dict mapping keys to StructureData or PKs
                   (e.g., {'clean': s1, 'defect': s2})
        code_label: VASP code label
        incar: Base INCAR parameters dict (applied to all calculations)
        kpoints_spacing: K-points spacing in A^-1 (default: 0.03)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options dict
        retrieve: List of files to retrieve
        incar_overrides: Per-structure INCAR overrides
                        (e.g., {'delta_0.05': {'NELECT': 191.95}})
        max_concurrent_jobs: Maximum parallel VASP jobs (default: unlimited)
        name: WorkGraph name
        wait: If True, block until all calculations finish
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean work directories after completion

    Returns:
        Dict mapping structure keys to individual WorkGraph PKs

    Example:
        >>> # Same settings for all structures
        >>> pks = quick_vasp_batch(
        ...     structures={'clean': s1, 'defect': s2},
        ...     code_label='VASP-6.5.1@localwork',
        ...     incar={'NSW': 100},
        ...     retrieve=['CONTCAR'],
        ...     max_concurrent_jobs=2,
        ... )

        >>> # Fukui-style with per-structure INCAR overrides
        >>> pks = quick_vasp_batch(
        ...     structures={'delta_0.00': s, 'delta_0.05': s, 'delta_0.10': s},
        ...     incar={'NSW': 0, 'ALGO': 'All'},
        ...     incar_overrides={
        ...         'delta_0.05': {'NELECT': 191.95},
        ...         'delta_0.10': {'NELECT': 191.90},
        ...     },
        ...     retrieve=['CHGCAR'],
        ...     max_concurrent_jobs=4,
        ... )
    """
    # Validate inputs
    if not structures:
        raise ValueError("structures dict cannot be empty")
    if code_label is None:
        raise ValueError("code_label is required")
    if incar is None:
        raise ValueError("incar is required - always specify INCAR parameters explicitly")
    if options is None:
        raise ValueError("options is required - specify scheduler resources")

    if incar_overrides is None:
        incar_overrides = {}

    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Build WorkGraph
    wg = WorkGraph(name=name)

    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    # Track task names for each key
    task_map = {}

    # Process each structure
    for key, struct_input in structures.items():
        # Load structure if PK
        if isinstance(struct_input, int):
            struct = orm.load_node(struct_input)
        else:
            struct = struct_input

        # Merge base INCAR with per-structure overrides
        if key in incar_overrides:
            merged_incar = deep_merge_dicts(incar, incar_overrides[key])
        else:
            merged_incar = incar

        # Prepare builder inputs
        builder_inputs = _prepare_builder_inputs(
            incar=merged_incar,
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping or {},
            options=options,
            retrieve=retrieve,
            restart_folder=None,
            clean_workdir=clean_workdir,
        )

        # Add VASP task for this structure
        task_name = f'vasp_{key}'
        vasp_task = wg.add_task(
            VaspTask,
            name=task_name,
            structure=struct,
            code=code,
            **builder_inputs
        )

        # Add energy extraction task
        energy_task_name = f'energy_{key}'
        wg.add_task(
            extract_energy,
            name=energy_task_name,
            misc=vasp_task.outputs.misc,
            retrieved=vasp_task.outputs.retrieved,
        )

        task_map[key] = {
            'vasp_task': task_name,
            'energy_task': energy_task_name,
        }

    # Submit
    wg.submit()

    # Wait if requested
    if wait:
        _wait_for_completion(wg.pk, poll_interval)

    # Return the WorkGraph PK with keys for reference
    # Users can extract individual results using get_batch_results_from_workgraph
    return {
        '__workgraph_pk__': wg.pk,
        '__task_map__': task_map,
        **{key: wg.pk for key in structures.keys()},
    }


def _prepare_builder_inputs(
    incar: dict,
    kpoints_spacing: float,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    retrieve: t.List[str] = None,
    restart_folder=None,
    clean_workdir: bool = False,
) -> dict:
    """
    Prepare builder inputs for VaspWorkChain.

    Args:
        incar: INCAR parameters dict
        kpoints_spacing: K-points spacing
        potential_family: POTCAR family
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options
        retrieve: Additional files to retrieve
        restart_folder: RemoteData for restart
        clean_workdir: Whether to clean work directory

    Returns:
        Dict of prepared inputs for VaspWorkChain
    """
    prepared = {}

    # Parameters (INCAR)
    prepared['parameters'] = orm.Dict(dict={'incar': incar})

    # K-points spacing
    prepared['kpoints_spacing'] = float(kpoints_spacing)

    # Potentials
    prepared['potential_family'] = potential_family
    prepared['potential_mapping'] = orm.Dict(dict=potential_mapping)

    # Options
    prepared['options'] = orm.Dict(dict=options)

    # Clean workdir
    prepared['clean_workdir'] = clean_workdir

    # Settings (for file retrieval)
    settings = {}
    if retrieve:
        settings['additional_retrieve_list'] = retrieve
    if settings:
        prepared['settings'] = orm.Dict(dict=settings)

    # Restart folder
    if restart_folder is not None:
        prepared['restart'] = {'folder': restart_folder}

    return prepared


def _wait_for_completion(pk: int, poll_interval: float) -> None:
    """
    Block until a WorkGraph completes.

    Args:
        pk: WorkGraph PK
        poll_interval: Seconds between status checks
    """
    print(f"Waiting for WorkGraph PK {pk} to complete...")

    while True:
        status = get_status(pk)

        if status in ('finished', 'failed', 'excepted', 'killed'):
            print(f"WorkGraph PK {pk} completed with status: {status}")
            break

        time.sleep(poll_interval)


def get_batch_results_from_workgraph(batch_result: dict) -> t.Dict[str, dict]:
    """
    Extract results from a quick_vasp_batch result.

    Args:
        batch_result: Return value from quick_vasp_batch

    Returns:
        Dict mapping structure keys to result dicts
    """
    from .results import get_results

    wg_pk = batch_result['__workgraph_pk__']
    task_map = batch_result['__task_map__']

    wg_node = orm.load_node(wg_pk)

    results = {}
    for key, task_info in task_map.items():
        # Extract results for this key
        result = {
            'energy': None,
            'structure': None,
            'misc': None,
            'files': None,
            'pk': wg_pk,
            'key': key,
        }

        # Access task outputs
        vasp_task_name = task_info['vasp_task']
        energy_task_name = task_info['energy_task']

        if hasattr(wg_node, 'tasks'):
            # Live WorkGraph
            tasks = wg_node.tasks
            if vasp_task_name in tasks:
                vasp_task = tasks[vasp_task_name]
                if hasattr(vasp_task.outputs, 'misc'):
                    misc_val = vasp_task.outputs.misc.value
                    if hasattr(misc_val, 'get_dict'):
                        result['misc'] = misc_val.get_dict()
                if hasattr(vasp_task.outputs, 'structure'):
                    struct_val = vasp_task.outputs.structure.value
                    if struct_val is not None:
                        result['structure'] = struct_val
                if hasattr(vasp_task.outputs, 'retrieved'):
                    result['files'] = vasp_task.outputs.retrieved.value

            if energy_task_name in tasks:
                energy_task = tasks[energy_task_name]
                if hasattr(energy_task.outputs, 'result'):
                    energy_val = energy_task.outputs.result.value
                    if energy_val is not None:
                        result['energy'] = energy_val.value if hasattr(energy_val, 'value') else float(energy_val)

        # Extract energy from misc if not found
        if result['energy'] is None and result['misc'] is not None:
            from .results import _extract_energy_from_misc
            result['energy'] = _extract_energy_from_misc(result['misc'])

        results[key] = result

    return results


def quick_dos(
    structure: t.Union[orm.StructureData, int] = None,
    code_label: str = None,
    scf_incar: dict = None,
    dos_incar: dict = None,
    kpoints_spacing: float = 0.03,
    dos_kpoints_spacing: float = None,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    retrieve: t.List[str] = None,
    name: str = 'quick_dos',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> int:
    """
    Submit a DOS calculation using BandsWorkChain with only_dos=True.

    This function uses AiiDA-VASP's BandsWorkChain which handles the
    SCF -> DOS workflow internally with proper CHGCAR/WAVECAR passing.

    Args:
        structure: StructureData or PK of structure to calculate DOS for.
        code_label: VASP code label (e.g., 'VASP-6.5.1@localwork')
        scf_incar: INCAR parameters for SCF stage (e.g., {'encut': 400, 'ediff': 1e-5}).
                   lwave and lcharg are forced to True internally by BandsWorkChain.
        dos_incar: INCAR parameters for DOS stage (e.g., {'nedos': 2000, 'lorbit': 11}).
                   Passed to the 'dos' namespace of BandsWorkChain.
        kpoints_spacing: K-points spacing in A^-1 for SCF (default: 0.03)
        dos_kpoints_spacing: K-points spacing for DOS (default: 80% of kpoints_spacing)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping (e.g., {'Sn': 'Sn_d', 'O': 'O'})
        options: Scheduler options dict with 'resources' key
        retrieve: List of additional files to retrieve (e.g., ['DOSCAR'])
        name: Calculation label for identification
        wait: If True, block until calculation finishes (default: False)
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean the work directory after completion

    Returns:
        PK of the submitted BandsWorkChain

    Example:
        >>> pk = quick_dos(
        ...     structure=my_structure,
        ...     code_label='VASP-6.5.1@localwork',
        ...     scf_incar={'encut': 400, 'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05},
        ...     dos_incar={'nedos': 2000, 'lorbit': 11, 'ismear': -5},
        ...     kpoints_spacing=0.03,
        ...     dos_kpoints_spacing=0.02,
        ...     potential_family='PBE',
        ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
        ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ...     retrieve=['DOSCAR'],
        ...     name='sno2_dos',
        ... )

    Note:
        AiiDA-VASP requires lowercase INCAR keys (e.g., 'encut' not 'ENCUT').
    """
    from aiida.engine import submit

    # Validate required inputs
    if structure is None:
        raise ValueError("structure is required")
    if code_label is None:
        raise ValueError("code_label is required")
    if scf_incar is None:
        raise ValueError("scf_incar is required - always specify SCF INCAR parameters explicitly")
    if dos_incar is None:
        raise ValueError("dos_incar is required - always specify DOS INCAR parameters explicitly")
    if options is None:
        raise ValueError("options is required - specify scheduler resources")

    # Load structure if PK
    if isinstance(structure, int):
        structure = orm.load_node(structure)

    # Default DOS k-points spacing to 80% of SCF spacing (denser)
    if dos_kpoints_spacing is None:
        dos_kpoints_spacing = kpoints_spacing * 0.8

    # Load code and BandsWorkChain
    code = orm.load_code(code_label)
    BandsWorkChain = WorkflowFactory('vasp.v2.bands')

    # Prepare INCAR parameters
    # Force lwave and lcharg for non-SCF DOS calculation
    scf_incar_final = dict(scf_incar)
    scf_incar_final['lwave'] = True
    scf_incar_final['lcharg'] = True
    # Ensure static calculation
    if 'nsw' not in scf_incar_final:
        scf_incar_final['nsw'] = 0
    if 'ibrion' not in scf_incar_final:
        scf_incar_final['ibrion'] = -1

    # DOS INCAR
    dos_incar_final = dict(dos_incar)
    if 'nsw' not in dos_incar_final:
        dos_incar_final['nsw'] = 0
    if 'ibrion' not in dos_incar_final:
        dos_incar_final['ibrion'] = -1

    # Band settings - only_dos mode
    # dos_kpoints_distance uses same units as kpoints_spacing (Å⁻¹)
    band_settings = {
        'only_dos': True,
        'run_dos': True,
        'dos_kpoints_distance': float(dos_kpoints_spacing),
    }

    # Build overrides for protocol-based builder
    # Note: settings must be {} not None to avoid aiida-vasp bug
    overrides = {
        'scf': {
            'parameters': {'incar': scf_incar_final},
            'kpoints_spacing': float(kpoints_spacing),
            'potential_family': potential_family,
            'potential_mapping': potential_mapping or {},
            'clean_workdir': False,  # Keep for DOS restart
            'settings': {},  # Must be dict, not None
        },
        'dos': {
            'parameters': {'incar': dos_incar_final},
            'settings': {},  # Must be dict, not None
        },
        'band_settings': band_settings,
    }

    # Add settings for file retrieval
    if retrieve:
        overrides['dos']['settings'] = {
            'additional_retrieve_list': retrieve,
            'parser_settings': {
                'include_node': ['dos'],
            },
        }

    # Clean workdir option
    if clean_workdir:
        overrides['clean_children_workdir'] = 'all'

    # Use protocol-based builder with run_relax=False to exclude relax namespace
    # Note: band_settings must be in overrides, not as separate parameter
    # (to avoid recursive_merge error with None)
    builder = BandsWorkChain.get_builder_from_protocol(
        code=code,
        structure=structure,
        run_relax=False,  # Critical: skip relaxation
        overrides=overrides,
        options=options,
    )

    # Set metadata label
    builder.metadata.label = name

    # Note: The `retrieve` parameter is passed through protocol overrides.
    # If DOSCAR doesn't appear in retrieved files, use export_files() to
    # manually download it from the remote folder after calculation completes.

    # Submit the builder directly (not unpacked)
    node = submit(builder)
    pk = node.pk

    # Wait if requested
    if wait:
        _wait_for_completion(pk, poll_interval)

    return pk
