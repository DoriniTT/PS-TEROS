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

from .tasks import extract_energy, compute_dynamics
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


def _builder_to_dict(builder) -> dict:
    """
    Recursively convert a ProcessBuilder to a plain dict.

    ProcessBuilderNamespace objects need to be converted to regular dicts
    for use with WorkGraph add_task().

    Args:
        builder: ProcessBuilder or ProcessBuilderNamespace

    Returns:
        Plain dict with all nested namespaces converted
    """
    from aiida.engine.processes.builder import ProcessBuilderNamespace

    result = {}
    for key, value in builder.items():
        if isinstance(value, ProcessBuilderNamespace):
            # Recursively convert nested namespaces
            nested = _builder_to_dict(value)
            if nested:  # Only include non-empty namespaces
                result[key] = nested
        elif value is not None:
            result[key] = value
    return result


def _prepare_builder_inputs(
    incar: dict,
    kpoints_spacing: float,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    retrieve: t.List[str] = None,
    restart_folder=None,
    clean_workdir: bool = False,
    kpoints_mesh: t.List[int] = None,
    structure: orm.StructureData = None,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
) -> dict:
    """
    Prepare builder inputs for VaspWorkChain.

    Args:
        incar: INCAR parameters dict
        kpoints_spacing: K-points spacing (used if kpoints_mesh not provided)
        potential_family: POTCAR family
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options
        retrieve: Additional files to retrieve
        restart_folder: RemoteData for restart
        clean_workdir: Whether to clean work directory
        kpoints_mesh: Explicit k-points mesh [nx, ny, nz] (overrides kpoints_spacing)
        structure: StructureData (required for selective dynamics)
        fix_type: Where to fix atoms ('bottom', 'center', 'top', or None)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix

    Returns:
        Dict of prepared inputs for VaspWorkChain
    """
    from ..fixed_atoms import get_fixed_atoms_list

    prepared = {}

    # Parameters (INCAR)
    prepared['parameters'] = orm.Dict(dict={'incar': incar})

    # K-points: explicit mesh or spacing
    if kpoints_mesh is not None:
        kpoints = orm.KpointsData()
        kpoints.set_kpoints_mesh(kpoints_mesh)
        prepared['kpoints'] = kpoints
    else:
        prepared['kpoints_spacing'] = float(kpoints_spacing)

    # Potentials
    prepared['potential_family'] = potential_family
    prepared['potential_mapping'] = orm.Dict(dict=potential_mapping)

    # Options
    prepared['options'] = orm.Dict(dict=options)

    # Clean workdir
    prepared['clean_workdir'] = clean_workdir

    # Settings (for file retrieval)
    # Note: aiida-vasp expects UPPERCASE keys for settings
    settings = {}
    if retrieve:
        settings['ADDITIONAL_RETRIEVE_LIST'] = retrieve
    if settings:
        prepared['settings'] = orm.Dict(dict=settings)

    # Restart folder
    if restart_folder is not None:
        prepared['restart'] = {'folder': restart_folder}

    # Selective dynamics (fix atoms)
    if structure is not None and fix_type is not None and fix_thickness > 0.0:
        fixed_atoms_list = get_fixed_atoms_list(
            structure=structure,
            fix_type=fix_type,
            fix_thickness=fix_thickness,
            fix_elements=fix_elements,
        )

        if fixed_atoms_list:
            # Create positions_dof array: True = relax, False = fix
            num_atoms = len(structure.sites)
            positions_dof = []

            for i in range(1, num_atoms + 1):  # 1-based indexing
                if i in fixed_atoms_list:
                    positions_dof.append([False, False, False])  # Fix atom
                else:
                    positions_dof.append([True, True, True])  # Relax atom

            prepared['dynamics'] = orm.Dict(dict={'positions_dof': positions_dof})

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


def _validate_stages(stages: t.List[dict]) -> None:
    """
    Validate sequential stage configuration.

    Args:
        stages: List of stage configuration dicts

    Raises:
        ValueError: If validation fails
    """
    if not stages:
        raise ValueError("stages list cannot be empty")

    stage_names = set()
    for i, stage in enumerate(stages):
        # Require name
        if 'name' not in stage:
            raise ValueError(f"Stage {i} missing required 'name' field")

        name = stage['name']
        if name in stage_names:
            raise ValueError(f"Duplicate stage name: '{name}'")
        stage_names.add(name)

        # Get stage type (default to 'vasp')
        stage_type = stage.get('type', 'vasp')
        valid_types = ('vasp', 'dos', 'batch')
        if stage_type not in valid_types:
            raise ValueError(
                f"Stage '{name}' type='{stage_type}' must be one of {valid_types}"
            )

        # Type-specific validation
        if stage_type == 'dos':
            # DOS stages require scf_incar and dos_incar
            if 'scf_incar' not in stage:
                raise ValueError(f"Stage '{name}': DOS stages require 'scf_incar'")
            if 'dos_incar' not in stage:
                raise ValueError(f"Stage '{name}': DOS stages require 'dos_incar'")
            if 'structure_from' not in stage:
                raise ValueError(f"Stage '{name}': DOS stages require 'structure_from'")

            # Validate structure_from for DOS stages
            structure_from = stage['structure_from']
            if structure_from not in stage_names:
                raise ValueError(
                    f"Stage '{name}' structure_from='{structure_from}' must reference "
                    f"a previous stage name"
                )

        elif stage_type == 'batch':
            # Batch stages require structure_from, base_incar, calculations
            if 'structure_from' not in stage:
                raise ValueError(f"Stage '{name}': batch stages require 'structure_from'")
            if 'base_incar' not in stage:
                raise ValueError(f"Stage '{name}': batch stages require 'base_incar'")
            if 'calculations' not in stage or not stage['calculations']:
                raise ValueError(
                    f"Stage '{name}': batch stages require non-empty 'calculations' dict"
                )

            # Validate structure_from references a previous stage
            structure_from = stage['structure_from']
            if structure_from not in stage_names:
                raise ValueError(
                    f"Stage '{name}' structure_from='{structure_from}' must reference "
                    f"a previous stage name"
                )

        else:  # stage_type == 'vasp'
            # VASP stages require incar
            if 'incar' not in stage:
                raise ValueError(f"Stage '{name}' missing required 'incar' field")

            # Require restart (must be None or a previous stage name)
            if 'restart' not in stage:
                raise ValueError(f"Stage '{name}' missing required 'restart' field (use None or a stage name)")

            restart = stage['restart']
            if restart is not None:
                if restart not in stage_names:
                    raise ValueError(
                        f"Stage '{name}' restart='{restart}' references unknown or "
                        f"later stage (must be defined before this stage)"
                    )

            # Validate structure_from for VASP stages
            structure_from = stage.get('structure_from', 'previous')
            if structure_from not in ('previous', 'input') and structure_from not in stage_names:
                raise ValueError(
                    f"Stage '{name}' structure_from='{structure_from}' must be 'previous', "
                    f"'input', or a previous stage name"
                )

            # Validate supercell spec
            if 'supercell' in stage:
                spec = stage['supercell']
                if not isinstance(spec, (list, tuple)) or len(spec) != 3:
                    raise ValueError(
                        f"Stage '{name}' supercell must be [nx, ny, nz], got: {spec}"
                    )
                for val in spec:
                    if not isinstance(val, int) or val < 1:
                        raise ValueError(
                            f"Stage '{name}' supercell values must be positive integers, "
                            f"got: {spec}"
                        )

            # Validate fix_type
            fix_type = stage.get('fix_type', None)
            if fix_type is not None:
                valid_fix_types = ('bottom', 'center', 'top')
                if fix_type not in valid_fix_types:
                    raise ValueError(
                        f"Stage '{name}' fix_type='{fix_type}' must be one of {valid_fix_types}"
                    )

                # If fix_type is set, fix_thickness must be positive
                fix_thickness = stage.get('fix_thickness', 0.0)
                if fix_thickness <= 0.0:
                    raise ValueError(
                        f"Stage '{name}' has fix_type='{fix_type}' but fix_thickness={fix_thickness}. "
                        f"fix_thickness must be > 0 when fix_type is set."
                    )


def _create_dos_stage_tasks(
    wg: WorkGraph,
    stage: dict,
    stage_name: str,
    input_structure,
    code,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    base_kpoints_spacing: float,
    clean_workdir: bool,
) -> dict:
    """
    Create DOS task using BandsWorkChain (vasp.v2.bands).

    BandsWorkChain handles the SCF → non-SCF DOS workflow internally with
    proper ICHARG settings and CHGCAR passing.

    Args:
        wg: WorkGraph to add tasks to
        stage: Stage configuration dict with keys:
            - scf_incar: INCAR dict for SCF step
            - dos_incar: INCAR dict for DOS step
            - kpoints_spacing: K-points spacing for SCF (default: base_kpoints_spacing)
            - kpoints: Explicit k-points mesh [nx, ny, nz] for SCF (overrides kpoints_spacing)
            - dos_kpoints_spacing: K-points spacing for DOS (default: scf_spacing * 0.8)
            - dos_kpoints: Explicit k-points mesh [nx, ny, nz] for DOS (overrides dos_kpoints_spacing)
            - retrieve: List of files to retrieve from DOS calc (default: ['DOSCAR'])
        stage_name: Unique stage identifier
        input_structure: Structure socket from previous stage
        code: VASP code node
        potential_family: POTCAR family
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options
        base_kpoints_spacing: Default k-points spacing
        clean_workdir: Whether to clean work directories

    Returns:
        Dict with task references for later stages
    """
    # Get BandsWorkChain and wrap as task
    BandsWorkChain = WorkflowFactory('vasp.v2.bands')
    BandsTask = task(BandsWorkChain)

    # Handle SCF k-points: explicit mesh or spacing
    scf_kpoints_mesh = stage.get('kpoints', None)
    scf_kpoints_spacing = stage.get('kpoints_spacing', base_kpoints_spacing)

    # Handle DOS k-points: explicit mesh or spacing
    dos_kpoints_mesh = stage.get('dos_kpoints', None)
    dos_kpoints_spacing = stage.get('dos_kpoints_spacing', scf_kpoints_spacing * 0.8)

    # Prepare SCF INCAR (BandsWorkChain handles lwave/lcharg internally)
    scf_incar = dict(stage['scf_incar'])
    scf_incar.update({
        'nsw': 0,        # Static calculation
        'ibrion': -1,    # No ionic relaxation
    })

    # Prepare DOS INCAR (BandsWorkChain handles ICHARG internally)
    dos_incar = dict(stage['dos_incar'])
    dos_incar.setdefault('ismear', -5)   # Tetrahedron method for DOS
    dos_incar.setdefault('lorbit', 11)   # Projected DOS
    dos_incar.setdefault('nedos', 2000)  # Number of DOS points
    dos_incar.update({
        'nsw': 0,        # Static calculation
        'ibrion': -1,    # No ionic relaxation
    })

    # Files to retrieve from DOS calc
    dos_retrieve = stage.get('retrieve', ['DOSCAR'])

    # DOS settings for file retrieval
    dos_settings = {}
    if dos_retrieve:
        dos_settings = {
            'ADDITIONAL_RETRIEVE_LIST': dos_retrieve,
        }

    # Prepare SCF input dict
    scf_input = {
        'code': code,
        'parameters': orm.Dict({'incar': scf_incar}),
        'potential_family': potential_family,
        'potential_mapping': orm.Dict(potential_mapping),
        'options': orm.Dict(options),
        'settings': orm.Dict({}),
        'clean_workdir': False,  # Keep for DOS restart
    }

    # SCF k-points: explicit mesh or spacing
    if scf_kpoints_mesh is not None:
        scf_kpoints = orm.KpointsData()
        scf_kpoints.set_kpoints_mesh(scf_kpoints_mesh)
        scf_input['kpoints'] = scf_kpoints
    else:
        scf_input['kpoints_spacing'] = float(scf_kpoints_spacing)

    # Prepare DOS input dict
    dos_input = {
        'code': code,
        'parameters': orm.Dict({'incar': dos_incar}),
        'potential_family': potential_family,
        'potential_mapping': orm.Dict(potential_mapping),
        'options': orm.Dict(options),
        'settings': orm.Dict(dos_settings),
    }

    # DOS k-points: explicit mesh or spacing
    if dos_kpoints_mesh is not None:
        dos_kpoints = orm.KpointsData()
        dos_kpoints.set_kpoints_mesh(dos_kpoints_mesh)
        dos_input['kpoints'] = dos_kpoints
        # Set dummy spacing in band_settings (required but ignored when kpoints is set)
        band_settings = orm.Dict({
            'only_dos': True,
            'run_dos': True,
            'dos_kpoints_distance': 0.03,  # Will be overridden by explicit kpoints
        })
    else:
        # Use spacing-based k-points for DOS
        band_settings = orm.Dict({
            'only_dos': True,
            'run_dos': True,
            'dos_kpoints_distance': float(dos_kpoints_spacing),
        })

    # Add BandsWorkChain task
    bands_task = wg.add_task(
        BandsTask,
        name=f'bands_{stage_name}',
        structure=input_structure,
        scf=scf_input,
        dos=dos_input,
        band_settings=band_settings,
        clean_children_workdir=orm.Str('all') if clean_workdir else None,
    )

    return {
        'bands_task': bands_task,
        'structure': input_structure,  # DOS doesn't modify structure
    }


def _create_batch_stage_tasks(
    wg: WorkGraph,
    stage: dict,
    stage_name: str,
    input_structure,
    code,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    base_kpoints_spacing: float,
    clean_workdir: bool,
) -> dict:
    """
    Create batch stage tasks: multiple parallel VASP calculations with varying parameters.

    For each entry in the stage's 'calculations' dict, a separate VaspTask and
    extract_energy task are created. Per-calculation overrides for incar, kpoints,
    and retrieve are merged with stage-level defaults.

    Args:
        wg: WorkGraph to add tasks to
        stage: Stage configuration dict with keys:
            - base_incar: Base INCAR dict applied to all calculations
            - calculations: Dict of {label: {incar: {overrides}, kpoints: [...], ...}}
            - kpoints_spacing: Default k-points spacing (fallback to base_kpoints_spacing)
            - kpoints: Default explicit k-points mesh [nx, ny, nz]
            - retrieve: Default files to retrieve
        stage_name: Unique stage identifier
        input_structure: Structure socket from previous stage
        code: VASP code node
        potential_family: POTCAR family
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options
        base_kpoints_spacing: Default k-points spacing
        clean_workdir: Whether to clean work directories

    Returns:
        Dict with:
            - calc_tasks: {label: vasp_task} for each calculation
            - energy_tasks: {label: energy_task} for each calculation
            - structure: input structure (unchanged, batch is static)
    """
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    base_incar = stage['base_incar']
    calculations = stage['calculations']

    # Stage-level defaults
    stage_kpoints_spacing = stage.get('kpoints_spacing', base_kpoints_spacing)
    stage_kpoints_mesh = stage.get('kpoints', None)
    stage_retrieve = stage.get('retrieve', None)

    calc_tasks = {}
    energy_tasks = {}

    for calc_label, calc_config in calculations.items():
        # Deep-merge base_incar with per-calculation incar overrides
        calc_incar_overrides = calc_config.get('incar', {})
        if calc_incar_overrides:
            merged_incar = deep_merge_dicts(base_incar, calc_incar_overrides)
        else:
            merged_incar = dict(base_incar)

        # Per-calc kpoints or fall back to stage-level
        calc_kpoints_mesh = calc_config.get('kpoints', stage_kpoints_mesh)
        calc_kpoints_spacing = calc_config.get('kpoints_spacing', stage_kpoints_spacing)

        # Per-calc retrieve or fall back to stage-level
        calc_retrieve = calc_config.get('retrieve', stage_retrieve)

        # Prepare builder inputs
        builder_inputs = _prepare_builder_inputs(
            incar=merged_incar,
            kpoints_spacing=calc_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            retrieve=calc_retrieve,
            restart_folder=None,
            clean_workdir=clean_workdir,
            kpoints_mesh=calc_kpoints_mesh,
        )

        # Add VASP task
        vasp_task_name = f'vasp_{stage_name}_{calc_label}'
        vasp_task = wg.add_task(
            VaspTask,
            name=vasp_task_name,
            structure=input_structure,
            code=code,
            **builder_inputs
        )

        # Add energy extraction task
        energy_task_name = f'energy_{stage_name}_{calc_label}'
        energy_task = wg.add_task(
            extract_energy,
            name=energy_task_name,
            misc=vasp_task.outputs.misc,
            retrieved=vasp_task.outputs.retrieved,
        )

        calc_tasks[calc_label] = vasp_task
        energy_tasks[calc_label] = energy_task

    return {
        'calc_tasks': calc_tasks,
        'energy_tasks': energy_tasks,
        'structure': input_structure,
    }


def quick_vasp_sequential(
    structure: t.Union[orm.StructureData, int] = None,
    stages: t.List[dict] = None,
    code_label: str = None,
    kpoints_spacing: float = 0.03,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    name: str = 'quick_vasp_sequential',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> dict:
    """
    Submit a multi-stage sequential VASP calculation with automatic restart chaining.

    Each stage runs after the previous one completes, automatically using the
    previous stage's output structure and remote_folder for restart. Supercell
    transformations can be inserted between stages.

    Args:
        structure: Initial StructureData or PK
        stages: List of stage configuration dicts (see Stage Configuration below)
        code_label: VASP code label (e.g., 'VASP-6.5.1@localwork')
        kpoints_spacing: Default k-points spacing in A^-1 (default: 0.03)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping (e.g., {'Sn': 'Sn_d', 'O': 'O'})
        options: Scheduler options dict with 'resources' key
        name: WorkGraph name for identification
        wait: If True, block until calculation finishes (default: False)
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean work directories after completion

    Returns:
        Dict with:
            - __workgraph_pk__: WorkGraph PK
            - __stage_names__: List of stage names in order
            - __stage_types__: Dict mapping stage names to types ('vasp' or 'dos')
            - <stage_name>: WorkGraph PK (for each stage)

    Stage Configuration (VASP stages, type='vasp' or omitted):
        - name (required): Unique stage identifier
        - type: 'vasp' (default) - standard VASP calculation
        - incar (required): INCAR parameters for this stage
        - restart (required): None or stage name to restart from
        - structure_from: 'previous' (default), 'input', or specific stage name
        - supercell: [nx, ny, nz] to create supercell
        - kpoints_spacing: Override k-points spacing for this stage
        - kpoints: Explicit k-points mesh [nx, ny, nz] (overrides kpoints_spacing)
        - retrieve: ADDITIONAL_RETRIEVE_LIST for this stage (e.g., ['CONTCAR', 'WAVECAR'])
        - fix_type: Where to fix atoms ('bottom', 'center', 'top', or None)
        - fix_thickness: Thickness in Angstroms for fixing region (required if fix_type set)
        - fix_elements: Optional list of element symbols to restrict fixing to

    Stage Configuration (DOS stages, type='dos'):
        - name (required): Unique stage identifier
        - type: 'dos' - DOS calculation (SCF + DOS)
        - structure_from (required): Stage name to get structure from
        - scf_incar (required): INCAR for SCF step (lwave/lcharg forced to True)
        - dos_incar (required): INCAR for DOS step (ismear defaults to -5, lorbit to 11)
        - kpoints_spacing: K-points for SCF (default: base value)
        - dos_kpoints_spacing: K-points for DOS (default: kpoints_spacing * 0.8)
        - retrieve: Files to retrieve from DOS step (default: ['DOSCAR'])

    Example:
        >>> stages = [
        ...     {
        ...         'name': 'relax_1x1_rough',
        ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'EDIFF': 1e-4, 'ENCUT': 400},
        ...         'restart': None,
        ...         'kpoints_spacing': 0.06,
        ...         'retrieve': ['CONTCAR', 'OUTCAR'],
        ...     },
        ...     {
        ...         'name': 'relax_1x1_fine',
        ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'EDIFF': 1e-6, 'ENCUT': 520},
        ...         'restart': 'relax_1x1_rough',
        ...         'kpoints_spacing': 0.03,
        ...         'retrieve': ['CONTCAR', 'OUTCAR'],
        ...     },
        ...     {
        ...         'name': 'relax_2x2_rough',
        ...         'supercell': [2, 2, 1],
        ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'EDIFF': 1e-4, 'ENCUT': 400},
        ...         'restart': None,
        ...         'kpoints': [2, 2, 1],
        ...         'retrieve': ['CONTCAR', 'OUTCAR'],
        ...     },
        ...     {
        ...         'name': 'relax_2x2_fine',
        ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'EDIFF': 1e-6, 'ENCUT': 520},
        ...         'restart': 'relax_2x2_rough',
        ...         'kpoints': [6, 6, 1],
        ...         'retrieve': ['CONTCAR', 'OUTCAR', 'CHGCAR'],
        ...     },
        ... ]
        >>>
        >>> result = quick_vasp_sequential(
        ...     structure=sno2_structure,
        ...     stages=stages,
        ...     code_label='VASP-6.5.1@localwork',
        ...     potential_family='PBE',
        ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
        ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ... )
        >>>
        >>> # Monitor: verdi process show result['__workgraph_pk__']
        >>> # Results: print_sequential_results(result)

    Note:
        - When supercell is specified, restart is automatically set to None
        - VaspWorkChain handles WAVECAR/CHGCAR copying from restart.folder automatically
    """
    from teros.core.aimd.tasks import create_supercell

    # Validate required inputs
    if structure is None:
        raise ValueError("structure is required")
    if stages is None:
        raise ValueError("stages is required - provide list of stage configurations")
    if code_label is None:
        raise ValueError("code_label is required")
    if options is None:
        raise ValueError("options is required - specify scheduler resources")

    # Validate stages
    _validate_stages(stages)

    # Load structure if PK
    if isinstance(structure, int):
        structure = orm.load_node(structure)

    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Build WorkGraph
    wg = WorkGraph(name=name)

    # Track structures and remote folders across stages
    stage_tasks = {}  # name -> {'vasp': task, 'energy': task, 'supercell': task, ...}
    stage_names = []  # Ordered list
    stage_types = {}  # name -> 'vasp' or 'dos'

    # Get the initial structure (may need to be wrapped for first supercell)
    current_structure = structure

    for i, stage in enumerate(stages):
        stage_name = stage['name']
        stage_type = stage.get('type', 'vasp')
        stage_names.append(stage_name)
        stage_types[stage_name] = stage_type

        if stage_type == 'dos':
            # DOS stage: get structure from specified stage
            structure_from = stage['structure_from']
            # Get structure from specified previous stage
            prev_stage_type = stage_types.get(structure_from, 'vasp')
            if prev_stage_type in ('dos', 'batch'):
                # DOS/batch stages use the input structure directly
                stage_structure = stage_tasks[structure_from]['structure']
            else:
                stage_structure = stage_tasks[structure_from]['vasp'].outputs.structure

            # Create DOS tasks
            dos_tasks = _create_dos_stage_tasks(
                wg=wg,
                stage=stage,
                stage_name=stage_name,
                input_structure=stage_structure,
                code=code,
                potential_family=potential_family,
                potential_mapping=potential_mapping or {},
                options=options,
                base_kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
            )

            # Store tasks for later reference
            stage_tasks[stage_name] = dos_tasks

            # Expose outputs for DOS stage (BandsWorkChain)
            bands_task = dos_tasks['bands_task']

            # Expose BandsWorkChain outputs that are available
            # These are optional - BandsWorkChain may or may not produce them
            # depending on parser settings and calculation success
            try:
                setattr(wg.outputs, f'{stage_name}_dos', bands_task.outputs.dos)
            except AttributeError:
                pass  # Will be extracted via link traversal
            try:
                setattr(wg.outputs, f'{stage_name}_projectors', bands_task.outputs.projectors)
            except AttributeError:
                pass  # Will be extracted via link traversal

            # Expose internal SCF workchain outputs (always available)
            try:
                setattr(wg.outputs, f'{stage_name}_scf_misc', bands_task.outputs.scf_misc)
            except AttributeError:
                pass
            try:
                setattr(wg.outputs, f'{stage_name}_scf_remote', bands_task.outputs.scf_remote_folder)
            except AttributeError:
                pass
            try:
                setattr(wg.outputs, f'{stage_name}_scf_retrieved', bands_task.outputs.scf_retrieved)
            except AttributeError:
                pass

            # Expose internal DOS workchain outputs (always available)
            try:
                setattr(wg.outputs, f'{stage_name}_dos_misc', bands_task.outputs.dos_misc)
            except AttributeError:
                pass
            try:
                setattr(wg.outputs, f'{stage_name}_dos_remote', bands_task.outputs.dos_remote_folder)
            except AttributeError:
                pass
            try:
                setattr(wg.outputs, f'{stage_name}_dos_retrieved', bands_task.outputs.dos_retrieved)
            except AttributeError:
                pass

        elif stage_type == 'batch':
            # Batch stage: multiple parallel VASP calculations on same structure
            structure_from = stage['structure_from']
            # Resolve structure from referenced stage
            ref_stage_type = stage_types.get(structure_from, 'vasp')
            if ref_stage_type == 'dos':
                stage_structure = stage_tasks[structure_from]['structure']
            elif ref_stage_type == 'batch':
                stage_structure = stage_tasks[structure_from]['structure']
            else:
                stage_structure = stage_tasks[structure_from]['vasp'].outputs.structure

            # Create batch tasks
            batch_tasks = _create_batch_stage_tasks(
                wg=wg,
                stage=stage,
                stage_name=stage_name,
                input_structure=stage_structure,
                code=code,
                potential_family=potential_family,
                potential_mapping=potential_mapping or {},
                options=options,
                base_kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
            )

            # Store tasks for later reference
            stage_tasks[stage_name] = batch_tasks

            # Expose outputs for each calculation in the batch
            for calc_label, vasp_task in batch_tasks['calc_tasks'].items():
                energy_task = batch_tasks['energy_tasks'][calc_label]
                setattr(wg.outputs, f'{stage_name}_{calc_label}_energy',
                        energy_task.outputs.result)
                setattr(wg.outputs, f'{stage_name}_{calc_label}_misc',
                        vasp_task.outputs.misc)
                setattr(wg.outputs, f'{stage_name}_{calc_label}_remote',
                        vasp_task.outputs.remote_folder)
                setattr(wg.outputs, f'{stage_name}_{calc_label}_retrieved',
                        vasp_task.outputs.retrieved)

        else:  # stage_type == 'vasp'
            # VASP stage: existing logic
            # Determine structure source
            structure_from = stage.get('structure_from', 'previous')
            if i == 0:
                # First stage always uses input structure
                stage_structure = current_structure
            elif structure_from == 'input':
                stage_structure = structure
            elif structure_from == 'previous':
                # Use output structure from previous stage
                prev_name = stage_names[i - 1]
                prev_stage_type = stage_types[prev_name]
                if prev_stage_type in ('dos', 'batch'):
                    # DOS/batch stages don't modify structure, use their input structure
                    stage_structure = stage_tasks[prev_name]['structure']
                else:
                    stage_structure = stage_tasks[prev_name]['vasp'].outputs.structure
            else:
                # Use output structure from specific stage
                ref_stage_type = stage_types.get(structure_from, 'vasp')
                if ref_stage_type in ('dos', 'batch'):
                    stage_structure = stage_tasks[structure_from]['structure']
                else:
                    stage_structure = stage_tasks[structure_from]['vasp'].outputs.structure

            # Handle supercell transformation
            supercell_task = None
            if 'supercell' in stage:
                supercell_spec = stage['supercell']
                supercell_task = wg.add_task(
                    create_supercell,
                    name=f'supercell_{stage_name}',
                    structure=stage_structure,
                    spec=orm.List(list=supercell_spec),
                )
                stage_structure = supercell_task.outputs.result

            # Determine restart source (None or stage name)
            restart = stage['restart']
            restart_folder = None
            if restart is not None:
                restart_folder = stage_tasks[restart]['vasp'].outputs.remote_folder

            # Prepare builder inputs for this stage
            stage_incar = stage['incar']
            stage_kpoints_spacing = stage.get('kpoints_spacing', kpoints_spacing)
            stage_kpoints_mesh = stage.get('kpoints', None)
            stage_retrieve = stage.get('retrieve', None)

            # Get fix parameters for this stage
            stage_fix_type = stage.get('fix_type', None)
            stage_fix_thickness = stage.get('fix_thickness', 0.0)
            stage_fix_elements = stage.get('fix_elements', None)

            # Determine if we can compute dynamics at build time
            # (only possible if structure is an actual StructureData, not a socket)
            is_structure_socket = not isinstance(stage_structure, orm.StructureData)

            # Prepare builder inputs
            # If fix_type is set and we have an actual structure, compute dynamics at build time
            if stage_fix_type is not None and not is_structure_socket:
                builder_inputs = _prepare_builder_inputs(
                    incar=stage_incar,
                    kpoints_spacing=stage_kpoints_spacing,
                    potential_family=potential_family,
                    potential_mapping=potential_mapping or {},
                    options=options,
                    retrieve=stage_retrieve,
                    restart_folder=None,  # We'll pass restart separately
                    clean_workdir=clean_workdir,
                    kpoints_mesh=stage_kpoints_mesh,
                    structure=stage_structure,
                    fix_type=stage_fix_type,
                    fix_thickness=stage_fix_thickness,
                    fix_elements=stage_fix_elements,
                )
            else:
                builder_inputs = _prepare_builder_inputs(
                    incar=stage_incar,
                    kpoints_spacing=stage_kpoints_spacing,
                    potential_family=potential_family,
                    potential_mapping=potential_mapping or {},
                    options=options,
                    retrieve=stage_retrieve,
                    restart_folder=None,  # We'll pass restart separately
                    clean_workdir=clean_workdir,
                    kpoints_mesh=stage_kpoints_mesh,
                )

            # If fix_type is set and structure is a socket, compute dynamics at runtime
            dynamics_task = None
            if stage_fix_type is not None and is_structure_socket:
                dynamics_task = wg.add_task(
                    compute_dynamics,
                    name=f'dynamics_{stage_name}',
                    structure=stage_structure,
                    fix_type=orm.Str(stage_fix_type),
                    fix_thickness=orm.Float(stage_fix_thickness),
                    fix_elements=orm.List(list=stage_fix_elements) if stage_fix_elements else None,
                )

            # Add VASP task
            vasp_task_kwargs = {
                'name': f'vasp_{stage_name}',
                'structure': stage_structure,
                'code': code,
                **builder_inputs
            }

            # Add restart if available
            if restart_folder is not None:
                vasp_task_kwargs['restart'] = {'folder': restart_folder}

            # Add dynamics from calcfunction if computed at runtime
            if dynamics_task is not None:
                vasp_task_kwargs['dynamics'] = dynamics_task.outputs.result

            vasp_task = wg.add_task(VaspTask, **vasp_task_kwargs)

            # Add energy extraction task
            energy_task = wg.add_task(
                extract_energy,
                name=f'energy_{stage_name}',
                misc=vasp_task.outputs.misc,
                retrieved=vasp_task.outputs.retrieved,
            )

            # Store tasks for later reference
            stage_tasks[stage_name] = {
                'vasp': vasp_task,
                'energy': energy_task,
                'supercell': supercell_task,
            }

            # Expose outputs for this stage
            setattr(wg.outputs, f'{stage_name}_energy', energy_task.outputs.result)
            setattr(wg.outputs, f'{stage_name}_structure', vasp_task.outputs.structure)
            setattr(wg.outputs, f'{stage_name}_misc', vasp_task.outputs.misc)
            setattr(wg.outputs, f'{stage_name}_remote', vasp_task.outputs.remote_folder)
            setattr(wg.outputs, f'{stage_name}_retrieved', vasp_task.outputs.retrieved)

    # Submit
    wg.submit()

    # Wait if requested
    if wait:
        _wait_for_completion(wg.pk, poll_interval)

    # Return result dict
    return {
        '__workgraph_pk__': wg.pk,
        '__stage_names__': stage_names,
        '__stage_types__': stage_types,
        **{name: wg.pk for name in stage_names},
    }


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


def quick_dos_batch(
    structures: t.Dict[str, t.Union[orm.StructureData, int]],
    code_label: str,
    scf_incar: dict,
    dos_incar: dict,
    kpoints_spacing: float = 0.03,
    dos_kpoints_spacing: float = None,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    retrieve: t.List[str] = None,
    scf_incar_overrides: t.Dict[str, dict] = None,
    dos_incar_overrides: t.Dict[str, dict] = None,
    max_concurrent_jobs: int = None,
    name: str = 'quick_dos_batch',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> t.Dict[str, int]:
    """
    Submit multiple DOS calculations in parallel using BandsWorkChain.

    Each structure runs through SCF -> DOS workflow with optional per-structure
    INCAR overrides. Useful for comparing DOS across different structures
    (e.g., pristine vs defects, different terminations).

    Args:
        structures: Dict mapping keys to StructureData or PKs
                   (e.g., {'pristine': s1, 'vacancy': s2})
        code_label: VASP code label
        scf_incar: Base INCAR for SCF stage (lowercase keys, e.g., {'encut': 400})
        dos_incar: Base INCAR for DOS stage (lowercase keys, e.g., {'nedos': 2000})
        kpoints_spacing: K-points spacing in A^-1 for SCF (default: 0.03)
        dos_kpoints_spacing: K-points spacing for DOS (default: 80% of kpoints_spacing)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping
        options: Scheduler options dict
        retrieve: List of files to retrieve (e.g., ['DOSCAR'])
        scf_incar_overrides: Per-structure SCF INCAR overrides
                            (e.g., {'vacancy': {'ismear': 0, 'sigma': 0.02}})
        dos_incar_overrides: Per-structure DOS INCAR overrides
        max_concurrent_jobs: Maximum parallel DOS jobs (default: unlimited)
        name: WorkGraph name
        wait: If True, block until all calculations finish
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean work directories after completion

    Returns:
        Dict with:
            - __workgraph_pk__: WorkGraph PK
            - __task_map__: Dict mapping keys to task names
            - <key>: WorkGraph PK (for each structure key)

    Example:
        >>> # Compare DOS for different structures
        >>> result = quick_dos_batch(
        ...     structures={'pristine': s1, 'vacancy': s2},
        ...     code_label='VASP-6.5.1@localwork',
        ...     scf_incar={'encut': 400, 'ediff': 1e-6, 'ismear': 0},
        ...     dos_incar={'nedos': 2000, 'lorbit': 11, 'ismear': -5},
        ...     kpoints_spacing=0.03,
        ...     max_concurrent_jobs=2,
        ... )
        >>> print(f"WorkGraph PK: {result['__workgraph_pk__']}")

        >>> # With per-structure INCAR overrides
        >>> result = quick_dos_batch(
        ...     structures={'metal': s1, 'insulator': s2},
        ...     scf_incar={'encut': 400, 'ediff': 1e-6},
        ...     dos_incar={'nedos': 2000, 'lorbit': 11},
        ...     scf_incar_overrides={
        ...         'metal': {'ismear': 1, 'sigma': 0.2},
        ...         'insulator': {'ismear': 0, 'sigma': 0.05},
        ...     },
        ...     dos_incar_overrides={
        ...         'metal': {'ismear': 1},
        ...         'insulator': {'ismear': -5},
        ...     },
        ... )

    Note:
        AiiDA-VASP requires lowercase INCAR keys (e.g., 'encut' not 'ENCUT').

    Exposed Outputs:
        For each structure key, the following outputs are exposed on the WorkGraph:
        - {key}_scf_misc: Dict with SCF calculation results
        - {key}_scf_remote: RemoteData for SCF calculation
        - {key}_scf_retrieved: FolderData with SCF retrieved files
        - {key}_dos_misc: Dict with DOS calculation results
        - {key}_dos_remote: RemoteData for DOS calculation
        - {key}_dos_retrieved: FolderData with DOS retrieved files (includes DOSCAR)
    """
    # Validate inputs
    if not structures:
        raise ValueError("structures dict cannot be empty")
    if code_label is None:
        raise ValueError("code_label is required")
    if scf_incar is None:
        raise ValueError("scf_incar is required - always specify SCF INCAR explicitly")
    if dos_incar is None:
        raise ValueError("dos_incar is required - always specify DOS INCAR explicitly")
    if options is None:
        raise ValueError("options is required - specify scheduler resources")

    if scf_incar_overrides is None:
        scf_incar_overrides = {}
    if dos_incar_overrides is None:
        dos_incar_overrides = {}

    # Default DOS k-points spacing to 80% of SCF spacing (denser)
    if dos_kpoints_spacing is None:
        dos_kpoints_spacing = kpoints_spacing * 0.8

    # Load code and wrap VaspWorkChain as task
    code = orm.load_code(code_label)
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
        if key in scf_incar_overrides:
            merged_scf_incar = deep_merge_dicts(scf_incar, scf_incar_overrides[key])
        else:
            merged_scf_incar = dict(scf_incar)

        if key in dos_incar_overrides:
            merged_dos_incar = deep_merge_dicts(dos_incar, dos_incar_overrides[key])
        else:
            merged_dos_incar = dict(dos_incar)

        # Prepare SCF INCAR - force lwave and lcharg for DOS restart
        scf_incar_final = dict(merged_scf_incar)
        scf_incar_final['lwave'] = True
        scf_incar_final['lcharg'] = True
        if 'nsw' not in scf_incar_final:
            scf_incar_final['nsw'] = 0
        if 'ibrion' not in scf_incar_final:
            scf_incar_final['ibrion'] = -1

        # Prepare DOS INCAR
        # Note: Don't set ISTART/ICHARG - the restart.folder mechanism in
        # VaspWorkChain handles WAVECAR/CHGCAR copying automatically
        dos_incar_final = dict(merged_dos_incar)
        if 'nsw' not in dos_incar_final:
            dos_incar_final['nsw'] = 0
        if 'ibrion' not in dos_incar_final:
            dos_incar_final['ibrion'] = -1

        # Prepare SCF builder inputs
        scf_builder_inputs = _prepare_builder_inputs(
            incar=scf_incar_final,
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping or {},
            options=options,
            retrieve=None,  # No special retrieval for SCF
            restart_folder=None,
            clean_workdir=False,  # Keep for DOS restart
        )

        # Add SCF task
        scf_task_name = f'scf_{key}'
        scf_task = wg.add_task(
            VaspTask,
            name=scf_task_name,
            structure=struct,
            code=code,
            **scf_builder_inputs
        )

        # Prepare DOS builder inputs
        # Note: We prepare the base inputs here, then add restart in add_task
        dos_builder_inputs = _prepare_builder_inputs(
            incar=dos_incar_final,
            kpoints_spacing=dos_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping or {},
            options=options,
            retrieve=retrieve,
            restart_folder=None,  # Will be passed directly below
            clean_workdir=clean_workdir,
        )

        # Add DOS task with restart from SCF
        # Pass restart directly in add_task to avoid potential issues with set_inputs
        dos_task_name = f'dos_{key}'
        dos_task = wg.add_task(
            VaspTask,
            name=dos_task_name,
            structure=struct,
            code=code,
            restart={'folder': scf_task.outputs.remote_folder},  # Wire restart here
            **dos_builder_inputs
        )

        # Expose outputs for this structure
        # SCF outputs
        setattr(wg.outputs, f'{key}_scf_misc', scf_task.outputs.misc)
        setattr(wg.outputs, f'{key}_scf_remote', scf_task.outputs.remote_folder)
        setattr(wg.outputs, f'{key}_scf_retrieved', scf_task.outputs.retrieved)
        # DOS outputs
        setattr(wg.outputs, f'{key}_dos_misc', dos_task.outputs.misc)
        setattr(wg.outputs, f'{key}_dos_remote', dos_task.outputs.remote_folder)
        setattr(wg.outputs, f'{key}_dos_retrieved', dos_task.outputs.retrieved)

        task_map[key] = {
            'scf_task': scf_task_name,
            'dos_task': dos_task_name,
        }

    # Submit
    wg.submit()

    # Wait if requested
    if wait:
        _wait_for_completion(wg.pk, poll_interval)

    # Return the WorkGraph PK with keys for reference
    return {
        '__workgraph_pk__': wg.pk,
        '__task_map__': task_map,
        **{key: wg.pk for key in structures.keys()},
    }


def quick_dos(
    structure: t.Union[orm.StructureData, int] = None,
    code_label: str = None,
    scf_incar: dict = None,
    dos_incar: dict = None,
    kpoints_spacing: float = 0.03,
    kpoints: t.List[int] = None,
    dos_kpoints_spacing: float = None,
    dos_kpoints: t.List[int] = None,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    retrieve: t.List[str] = None,
    name: str = 'quick_dos',
    wait: bool = False,
    poll_interval: float = 10.0,
    clean_workdir: bool = False,
) -> dict:
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
        kpoints: Explicit k-points mesh for SCF [nx, ny, nz] (overrides kpoints_spacing)
        dos_kpoints_spacing: K-points spacing for DOS (default: 80% of kpoints_spacing)
        dos_kpoints: Explicit k-points mesh for DOS [nx, ny, nz] (overrides dos_kpoints_spacing)
        potential_family: POTCAR family (default: 'PBE')
        potential_mapping: Element to POTCAR mapping (e.g., {'Sn': 'Sn_d', 'O': 'O'})
        options: Scheduler options dict with 'resources' key
        retrieve: List of additional files to retrieve (e.g., ['DOSCAR'])
        name: Calculation label for identification
        wait: If True, block until calculation finishes (default: False)
        poll_interval: Seconds between status checks when wait=True
        clean_workdir: Whether to clean the work directory after completion

    Returns:
        Dict with '__workgraph_pk__' key containing the BandsWorkChain PK

    Example (using spacing):
        >>> result = quick_dos(
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

    Example (using explicit k-points - recommended for slabs):
        >>> result = quick_dos(
        ...     structure=my_structure,
        ...     code_label='VASP-6.5.1@localwork',
        ...     scf_incar={'encut': 400, 'ediff': 1e-6},
        ...     dos_incar={'nedos': 3000, 'lorbit': 11, 'ismear': -5},
        ...     kpoints=[6, 6, 1],          # SCF k-points mesh
        ...     dos_kpoints=[8, 8, 1],      # DOS k-points mesh (denser)
        ...     potential_family='PBE',
        ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
        ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ...     retrieve=['DOSCAR', 'PROCAR', 'vasprun.xml'],
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
    # Only used if explicit dos_kpoints mesh is not provided
    if dos_kpoints_spacing is None and dos_kpoints is None:
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
    # dos_kpoints_distance is required but will be ignored if explicit kpoints is set
    band_settings = {
        'only_dos': True,
        'run_dos': True,
        'dos_kpoints_distance': float(dos_kpoints_spacing) if dos_kpoints_spacing else 0.03,
    }

    # Build overrides for protocol-based builder
    # Note: settings must be {} not None to avoid aiida-vasp bug
    scf_overrides = {
        'parameters': {'incar': scf_incar_final},
        'potential_family': potential_family,
        'potential_mapping': potential_mapping or {},
        'clean_workdir': False,  # Keep for DOS restart
        'settings': {},  # Must be dict, not None
    }

    # SCF k-points: explicit mesh or spacing
    if kpoints is not None:
        scf_kpoints_data = orm.KpointsData()
        scf_kpoints_data.set_kpoints_mesh(kpoints)
        scf_overrides['kpoints'] = scf_kpoints_data
    else:
        scf_overrides['kpoints_spacing'] = float(kpoints_spacing)

    # DOS overrides
    dos_overrides = {
        'parameters': {'incar': dos_incar_final},
        'settings': {},  # Must be dict, not None
    }

    # DOS k-points: explicit mesh or spacing
    if dos_kpoints is not None:
        dos_kpoints_data = orm.KpointsData()
        dos_kpoints_data.set_kpoints_mesh(dos_kpoints)
        dos_overrides['kpoints'] = dos_kpoints_data

    overrides = {
        'scf': scf_overrides,
        'dos': dos_overrides,
        'band_settings': band_settings,
    }

    # Add settings for file retrieval
    # Note: aiida-vasp expects UPPERCASE keys for settings
    if retrieve:
        overrides['dos']['settings'] = {
            'ADDITIONAL_RETRIEVE_LIST': retrieve,
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

    # Return dict for consistency with other quick_* functions
    return {
        '__workgraph_pk__': pk,
    }
