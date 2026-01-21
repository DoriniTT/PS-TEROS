"""WorkGraph builder for NEB (Nudged Elastic Band) calculations.

This module provides the main workflow builder for NEB calculations using
the aiida-vasp NEB plugin. It supports:
- Optional endpoint relaxation
- IDPP or linear interpolation
- Two-stage NEB (regular NEB → CI-NEB) for accurate barrier calculation
- Single-stage NEB for quick estimates

The workflow follows the PS-TEROS architectural patterns for consistency
with other modules.
"""

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task

from .tasks import (
    interpolate_structures_output,
    extract_neb_energies,
    calculate_barrier,
    extract_total_energy,
    create_neb_summary,
)
from .utils import (
    get_neb_incar_parameters,
    validate_neb_structures,
    format_neb_results,
    get_endpoint_relax_incar,
    get_stage2_parameters,
)
from ..utils import deep_merge_dicts


def _prepare_relax_inputs(
    builder_inputs: dict,
    structure: orm.StructureData,
    code: orm.InstalledCode,
) -> dict:
    """
    Prepare inputs for endpoint relaxation calculation.

    Args:
        builder_inputs: Base builder inputs dict
        structure: Structure to relax
        code: VASP code

    Returns:
        dict with AiiDA-compatible types ready for VaspTask
    """
    # Get default relaxation INCAR
    relax_incar = get_endpoint_relax_incar()

    # Create override dict
    relax_overrides = {
        'parameters': {
            'incar': relax_incar
        }
    }

    # Merge with user inputs (overrides take precedence for INCAR keys)
    merged = deep_merge_dicts(builder_inputs, relax_overrides)

    # Convert to AiiDA types
    prepared = {
        'structure': structure,
        'code': code,
    }

    # Convert dict-type parameters to orm.Dict
    for key in ('parameters', 'options', 'potential_mapping', 'settings'):
        if key in merged:
            if isinstance(merged[key], dict):
                prepared[key] = orm.Dict(dict=merged[key])
            else:
                prepared[key] = merged[key]

    # Handle kpoints_spacing
    if 'kpoints_spacing' in merged:
        kps = merged['kpoints_spacing']
        if isinstance(kps, (int, float)):
            prepared['kpoints_spacing'] = float(kps)
        else:
            prepared['kpoints_spacing'] = kps

    # Copy string/bool values directly
    for key in ('potential_family', 'clean_workdir'):
        if key in merged:
            prepared[key] = merged[key]

    return prepared


def _prepare_neb_inputs(
    builder_inputs: dict,
    initial_structure: t.Union[orm.StructureData, t.Any],
    final_structure: t.Union[orm.StructureData, t.Any],
    neb_images_dict: dict,
    code: orm.InstalledCode,
    n_images: int,
    climb: bool,
    spring_constant: float,
    neb_optimizer: int,
    force_convergence: float,
    max_steps: int,
    restart_folder: orm.RemoteData = None,
) -> dict:
    """
    Prepare inputs for NEB calculation.

    Args:
        builder_inputs: Base builder inputs dict
        initial_structure: Initial structure (or socket from relaxation)
        final_structure: Final structure (or socket from relaxation)
        neb_images_dict: Dictionary of intermediate images
        code: VASP code
        n_images: Number of intermediate images
        climb: Enable climbing image NEB
        spring_constant: NEB spring constant
        neb_optimizer: VTST optimizer (1=LBFGS, 3=FIRE)
        force_convergence: EDIFFG value
        max_steps: Maximum ionic steps (NSW)
        restart_folder: Optional RemoteData for restart

    Returns:
        dict with inputs ready for VaspNEBWorkChain
    """
    # Get NEB-specific INCAR parameters
    neb_incar = get_neb_incar_parameters(
        n_images=n_images,
        climb=climb,
        spring_constant=spring_constant,
        neb_optimizer=neb_optimizer,
        force_convergence=force_convergence,
        max_steps=max_steps,
    )

    # Create override dict
    neb_overrides = {
        'parameters': {
            'incar': neb_incar
        }
    }

    # Merge with user inputs
    merged = deep_merge_dicts(builder_inputs, neb_overrides)

    # Convert to AiiDA types
    prepared = {
        'initial_structure': initial_structure,
        'final_structure': final_structure,
        'code': code,
    }

    # Add NEB images as namespace inputs
    # VaspNEBWorkChain expects neb_images as a namespace with keys '01', '02', etc.
    for key, struct in neb_images_dict.items():
        # Convert 'image_01' to '01'
        image_num = key.replace('image_', '')
        if f'neb_images__{image_num}' not in prepared:
            prepared[f'neb_images__{image_num}'] = struct

    # Convert dict-type parameters to orm.Dict
    for key in ('parameters', 'options', 'potential_mapping', 'settings'):
        if key in merged:
            if isinstance(merged[key], dict):
                prepared[key] = orm.Dict(dict=merged[key])
            else:
                prepared[key] = merged[key]

    # Handle kpoints_spacing
    if 'kpoints_spacing' in merged:
        kps = merged['kpoints_spacing']
        if isinstance(kps, (int, float)):
            prepared['kpoints_spacing'] = orm.Float(kps)
        else:
            prepared['kpoints_spacing'] = kps

    # Copy string/bool values
    if 'potential_family' in merged:
        prepared['potential_family'] = orm.Str(merged['potential_family'])

    if 'clean_workdir' in merged:
        prepared['clean_workdir'] = orm.Bool(merged['clean_workdir'])

    # Add restart folder if provided
    if restart_folder is not None:
        prepared['restart_folder'] = restart_folder

    return prepared


def build_neb_workgraph(
    initial_structure: t.Union[orm.StructureData, int],
    final_structure: t.Union[orm.StructureData, int],
    n_images: int,
    code_label: str,
    builder_inputs: dict,
    relax_endpoints: bool = True,
    interpolation_method: str = 'idpp',
    climb: bool = True,
    spring_constant: float = -5.0,
    neb_optimizer: int = 1,
    force_convergence: float = -0.05,
    max_steps: int = 500,
    restart_folder: t.Union[orm.RemoteData, int] = None,
    name: str = 'NEBWorkGraph',
    max_concurrent_jobs: int = None,
) -> WorkGraph:
    """
    Build a WorkGraph for NEB (Nudged Elastic Band) calculations.

    Creates a workflow that calculates the minimum energy path (MEP) and
    activation barrier between two structures (reactant and product).

    Workflow stages:
    1. (Optional) Relax endpoint structures using vasp.v2.vasp
    2. Generate intermediate images via IDPP/linear interpolation
    3. Stage 1: Run NEB (LCLIMB=False) to get initial MEP
    4. (Optional) Stage 2: Run CI-NEB (LCLIMB=True) for accurate saddle point
    5. Calculate forward/reverse barriers from final energies

    Args:
        initial_structure: Initial structure (StructureData or PK)
        final_structure: Final structure (StructureData or PK)
        n_images: Number of intermediate NEB images (typically 3-7)
        code_label: VASP code label (must be VTST-compiled for full NEB support)
        builder_inputs: VASP builder configuration dict containing:
            - parameters: {'incar': {...}}
            - options: {'resources': {...}, ...}
            - kpoints_spacing: float (Å⁻¹)
            - potential_family: str
            - potential_mapping: dict
        relax_endpoints: If True, relax initial/final structures first.
                        Uses ISIF=2, NSW=100, IBRION=2. Default: True
        interpolation_method: 'idpp' (recommended) or 'linear'. Default: 'idpp'
        climb: If True, run two-stage NEB→CI-NEB workflow.
              If False, run single NEB stage only. Default: True
        spring_constant: NEB spring constant. Negative = variable nudging.
                        Default: -5.0 (VASP VTST default)
        neb_optimizer: VTST optimizer (IOPT). 1=LBFGS (default), 3=FIRE
        force_convergence: EDIFFG value. Negative = force-based (eV/Å).
                          Default: -0.05 eV/Å
        max_steps: Maximum ionic steps per NEB stage (NSW). Default: 500
        restart_folder: RemoteData (or PK) from previous NEB for restart
        name: WorkGraph name. Default: 'NEBWorkGraph'
        max_concurrent_jobs: Limit concurrent calculations (optional)

    Returns:
        WorkGraph ready for submission

    Outputs (available after completion):
        - barrier: Dict with forward_barrier, reverse_barrier, reaction_energy,
                  saddle_point_index, saddle_point_energy, energies_list
        - summary: Dict with comprehensive NEB results
        - misc: Dict with raw NEB output from final stage
        - remote_folder: RemoteData for restarts/analysis
        - stage1_remote_folder: RemoteData from Stage 1 (if climb=True)
        - relaxed_initial: StructureData (if relax_endpoints=True)
        - relaxed_final: StructureData (if relax_endpoints=True)

    Example:
        >>> from teros.core.neb import build_neb_workgraph, print_neb_summary
        >>> wg = build_neb_workgraph(
        ...     initial_structure=initial,
        ...     final_structure=final,
        ...     n_images=5,
        ...     code_label='VASP-6.4.3@bohr',
        ...     builder_inputs={
        ...         'parameters': {'incar': {'encut': 520, 'ismear': 0}},
        ...         'options': {
        ...             'resources': {'num_machines': 3, 'num_cores_per_machine': 40},
        ...             'queue_name': 'par120',
        ...         },
        ...         'kpoints_spacing': 0.03,
        ...         'potential_family': 'PBE',
        ...         'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
        ...     },
        ...     relax_endpoints=True,
        ...     climb=True,
        ... )
        >>> wg.submit(wait=False)
        >>> # After completion:
        >>> print_neb_summary(wg.pk)

    Notes:
        - VASP must be compiled with VTST extensions for IOPT optimizer
        - For accurate barriers, use climb=True (two-stage workflow)
        - IDPP interpolation provides better initial paths than linear
        - Endpoint relaxation ensures consistent reference energies
        - Spring constant < 0 uses variable nudged elastic band
    """
    # Load structures if PKs provided
    if isinstance(initial_structure, int):
        initial_structure = orm.load_node(initial_structure)
    if isinstance(final_structure, int):
        final_structure = orm.load_node(final_structure)
    if isinstance(restart_folder, int):
        restart_folder = orm.load_node(restart_folder)

    # Validate structures are compatible
    validate_neb_structures(initial_structure, final_structure)

    # Load code
    code = orm.load_code(code_label)

    # Set defaults for builder_inputs
    if builder_inputs is None:
        builder_inputs = {}

    # Get VaspWorkChain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspNEBWorkChain = WorkflowFactory('vasp.neb')
    VaspTask = task(VaspWorkChain)
    NEBTask = task(VaspNEBWorkChain)

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # Track structures for NEB (may be from relaxation or direct input)
    initial_for_neb = initial_structure
    final_for_neb = final_structure

    # Step 1: Optional endpoint relaxation
    if relax_endpoints:
        # Prepare relaxation inputs
        relax_initial_inputs = _prepare_relax_inputs(
            builder_inputs, initial_structure, code
        )
        relax_final_inputs = _prepare_relax_inputs(
            builder_inputs, final_structure, code
        )

        # Add relaxation tasks (run in parallel)
        relax_initial_task = wg.add_task(
            VaspTask,
            name='relax_initial',
            **relax_initial_inputs
        )
        relax_final_task = wg.add_task(
            VaspTask,
            name='relax_final',
            **relax_final_inputs
        )

        # Use relaxed structures for NEB
        initial_for_neb = relax_initial_task.outputs.structure
        final_for_neb = relax_final_task.outputs.structure

        # Expose relaxed structures as outputs
        wg.outputs.relaxed_initial = relax_initial_task.outputs.structure
        wg.outputs.relaxed_final = relax_final_task.outputs.structure

        # Extract energies from relaxation for summary
        initial_energy_task = wg.add_task(
            extract_total_energy,
            name='extract_initial_energy',
            misc=relax_initial_task.outputs.misc,
        )
        final_energy_task = wg.add_task(
            extract_total_energy,
            name='extract_final_energy',
            misc=relax_final_task.outputs.misc,
        )
    else:
        # No relaxation - need to do static calculations for endpoint energies
        # or use the NEB endpoints
        initial_energy_task = None
        final_energy_task = None

    # Step 2: Generate intermediate images via interpolation
    interp_task = wg.add_task(
        interpolate_structures_output,
        name='interpolate_structures',
        initial=initial_for_neb,
        final=final_for_neb,
        n_images=orm.Int(n_images),
        method=orm.Str(interpolation_method),
    )

    # Step 3: Stage 1 - Regular NEB (LCLIMB=False)
    # Get NEB-specific INCAR parameters for Stage 1
    neb_incar_stage1 = get_neb_incar_parameters(
        n_images=n_images,
        climb=False,  # Stage 1: no climbing image
        spring_constant=spring_constant,
        neb_optimizer=neb_optimizer,
        force_convergence=force_convergence,
        max_steps=max_steps,
    )

    # Merge with user parameters
    neb_params_stage1 = deep_merge_dicts(
        builder_inputs.get('parameters', {}),
        {'incar': neb_incar_stage1}
    )

    # Prepare NEB inputs
    neb_inputs_stage1 = {
        'initial_structure': initial_for_neb,
        'final_structure': final_for_neb,
        'code': code,
        'parameters': orm.Dict(dict=neb_params_stage1),
    }

    # Add common builder inputs
    if 'options' in builder_inputs:
        neb_inputs_stage1['options'] = orm.Dict(dict=builder_inputs['options'])
    if 'potential_family' in builder_inputs:
        neb_inputs_stage1['potential_family'] = orm.Str(builder_inputs['potential_family'])
    if 'potential_mapping' in builder_inputs:
        neb_inputs_stage1['potential_mapping'] = orm.Dict(dict=builder_inputs['potential_mapping'])
    if 'kpoints_spacing' in builder_inputs:
        neb_inputs_stage1['kpoints_spacing'] = orm.Float(builder_inputs['kpoints_spacing'])
    if 'settings' in builder_inputs:
        neb_inputs_stage1['settings'] = orm.Dict(dict=builder_inputs['settings'])

    # Add restart folder if provided
    if restart_folder is not None:
        neb_inputs_stage1['restart_folder'] = restart_folder

    # Create Stage 1 NEB task
    # Note: neb_images are passed dynamically from interpolation output
    neb_stage1_task = wg.add_task(
        NEBTask,
        name='neb_stage1',
        **neb_inputs_stage1,
    )

    # Link interpolated images to NEB task
    # The interpolation task returns a dict with 'image_01', 'image_02', etc.
    # We need to connect these to the neb_images namespace
    for i in range(1, n_images + 1):
        image_key = f'image_{i:02d}'
        neb_key = f'neb_images__{i:02d}'
        # Connect interpolation output to NEB input
        wg.add_link(
            interp_task.outputs[image_key],
            neb_stage1_task.inputs[neb_key],
        )

    # Expose Stage 1 remote folder for restarts
    wg.outputs.stage1_remote_folder = neb_stage1_task.outputs.remote_folder

    # Step 4: Stage 2 - CI-NEB (optional, LCLIMB=True)
    if climb:
        # Get Stage 2 parameters (tighter convergence)
        neb_incar_stage2 = get_stage2_parameters(
            neb_incar_stage1,
            tighter_convergence=True,
        )

        neb_params_stage2 = deep_merge_dicts(
            builder_inputs.get('parameters', {}),
            {'incar': neb_incar_stage2}
        )

        # Stage 2 uses restart from Stage 1
        neb_inputs_stage2 = {
            'initial_structure': initial_for_neb,
            'final_structure': final_for_neb,
            'code': code,
            'parameters': orm.Dict(dict=neb_params_stage2),
            'restart_folder': neb_stage1_task.outputs.remote_folder,
        }

        # Add common builder inputs
        if 'options' in builder_inputs:
            neb_inputs_stage2['options'] = orm.Dict(dict=builder_inputs['options'])
        if 'potential_family' in builder_inputs:
            neb_inputs_stage2['potential_family'] = orm.Str(builder_inputs['potential_family'])
        if 'potential_mapping' in builder_inputs:
            neb_inputs_stage2['potential_mapping'] = orm.Dict(dict=builder_inputs['potential_mapping'])
        if 'kpoints_spacing' in builder_inputs:
            neb_inputs_stage2['kpoints_spacing'] = orm.Float(builder_inputs['kpoints_spacing'])
        if 'settings' in builder_inputs:
            neb_inputs_stage2['settings'] = orm.Dict(dict=builder_inputs['settings'])

        # Create Stage 2 CI-NEB task
        neb_stage2_task = wg.add_task(
            NEBTask,
            name='neb_cineb',
            **neb_inputs_stage2,
        )

        # Link interpolated images to Stage 2 as well
        for i in range(1, n_images + 1):
            image_key = f'image_{i:02d}'
            neb_key = f'neb_images__{i:02d}'
            wg.add_link(
                interp_task.outputs[image_key],
                neb_stage2_task.inputs[neb_key],
            )

        final_neb_task = neb_stage2_task
    else:
        final_neb_task = neb_stage1_task

    # Step 5: Extract energies and calculate barrier
    energies_task = wg.add_task(
        extract_neb_energies,
        name='extract_neb_energies',
        misc=final_neb_task.outputs.misc,
    )

    barrier_task = wg.add_task(
        calculate_barrier,
        name='calculate_barrier',
        energies=energies_task.outputs.result,
    )

    # Step 6: Create summary
    # Only create summary if we have endpoint energies
    if initial_energy_task is not None and final_energy_task is not None:
        summary_task = wg.add_task(
            create_neb_summary,
            name='create_summary',
            initial_energy=initial_energy_task.outputs.result,
            final_energy=final_energy_task.outputs.result,
            barrier_results=barrier_task.outputs.result,
            n_images=orm.Int(n_images),
            method=orm.Str(interpolation_method),
            climb=orm.Bool(climb),
        )
        wg.outputs.summary = summary_task.outputs.result

    # Expose outputs
    wg.outputs.barrier = barrier_task.outputs.result
    wg.outputs.energies = energies_task.outputs.result
    wg.outputs.misc = final_neb_task.outputs.misc
    wg.outputs.remote_folder = final_neb_task.outputs.remote_folder

    # Expose trajectory if available
    # Note: TrajectoryData may not be present in all NEB outputs
    # wg.outputs.trajectory = final_neb_task.outputs.trajectory

    # Set max_concurrent_jobs
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg


def get_neb_results(workgraph) -> dict:
    """
    Extract results from completed NEB WorkGraph.

    Args:
        workgraph: Completed WorkGraph node (PK, UUID, or WorkGraphNode)

    Returns:
        dict with:
            - barrier: Dict with barrier information (forward, reverse, reaction energy)
            - summary: Dict with comprehensive NEB results (if available)
            - energies: Dict with per-image energies
            - relaxed_initial: StructureData (if relax_endpoints was True)
            - relaxed_final: StructureData (if relax_endpoints was True)
            - remote_folder: RemoteData for further analysis
            - stage1_remote_folder: RemoteData from Stage 1 (if climb was True)

    Example:
        >>> results = get_neb_results(wg.pk)
        >>> print(f"Forward barrier: {results['barrier']['forward_barrier']:.3f} eV")
    """
    # Load workgraph if PK or UUID provided
    if isinstance(workgraph, (int, str)):
        workgraph = orm.load_node(workgraph)

    results = {}

    # Get barrier results
    if hasattr(workgraph.outputs, 'barrier'):
        barrier_node = workgraph.outputs.barrier
        if isinstance(barrier_node, orm.Dict):
            results['barrier'] = barrier_node.get_dict()
        else:
            results['barrier'] = dict(barrier_node)
    else:
        results['barrier'] = None

    # Get summary results
    if hasattr(workgraph.outputs, 'summary'):
        summary_node = workgraph.outputs.summary
        if isinstance(summary_node, orm.Dict):
            results['summary'] = summary_node.get_dict()
        else:
            results['summary'] = dict(summary_node)
    else:
        results['summary'] = None

    # Get energies
    if hasattr(workgraph.outputs, 'energies'):
        energies_node = workgraph.outputs.energies
        if isinstance(energies_node, orm.Dict):
            results['energies'] = energies_node.get_dict()
        else:
            results['energies'] = dict(energies_node)
    else:
        results['energies'] = None

    # Get relaxed structures
    if hasattr(workgraph.outputs, 'relaxed_initial'):
        results['relaxed_initial'] = workgraph.outputs.relaxed_initial
    else:
        results['relaxed_initial'] = None

    if hasattr(workgraph.outputs, 'relaxed_final'):
        results['relaxed_final'] = workgraph.outputs.relaxed_final
    else:
        results['relaxed_final'] = None

    # Get remote folders
    if hasattr(workgraph.outputs, 'remote_folder'):
        results['remote_folder'] = workgraph.outputs.remote_folder
    else:
        results['remote_folder'] = None

    if hasattr(workgraph.outputs, 'stage1_remote_folder'):
        results['stage1_remote_folder'] = workgraph.outputs.stage1_remote_folder
    else:
        results['stage1_remote_folder'] = None

    # Get misc
    if hasattr(workgraph.outputs, 'misc'):
        misc_node = workgraph.outputs.misc
        if isinstance(misc_node, orm.Dict):
            results['misc'] = misc_node.get_dict()
        else:
            results['misc'] = dict(misc_node)
    else:
        results['misc'] = None

    return results


def print_neb_summary(workgraph) -> None:
    """
    Print formatted summary of NEB calculation results.

    Args:
        workgraph: Completed WorkGraph node (PK, UUID, or WorkGraphNode)

    Example:
        >>> print_neb_summary(wg.pk)
        ============================================================
        NEB CALCULATION RESULTS
        ============================================================

        Activation Barriers:
        ----------------------------------------
          Forward barrier:  0.4521 eV (43.62 kJ/mol)
          Reverse barrier:  0.3245 eV (31.31 kJ/mol)
          Reaction energy:  -0.1276 eV (-12.31 kJ/mol)
          Reaction type:    Exothermic
        ...
    """
    results = get_neb_results(workgraph)
    formatted = format_neb_results(results)
    print(formatted)
