"""
AIMD Module for PS-TEROS - CP2K Implementation

Ab initio molecular dynamics calculations on slab structures using CP2K.
Sequential AIMD stages with automatic restart chaining.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace

# Get CP2K workchain for type annotation
Cp2kBaseWorkChain = WorkflowFactory('cp2k.base')


@task.graph
def aimd_single_stage_scatter_cp2k(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    temperature: float,
    steps: int,
    code: orm.Code,
    aimd_parameters: dict,
    basis_file: orm.SinglefileData,
    pseudo_file: orm.SinglefileData,
    options: dict,  # Scheduler options - plain dict, wrapped with dict() when used
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    fixed_atoms_lists: dict = None,
    fix_components: str = "XYZ",
    # Dynamic fixed atoms parameters (for auto-generated slabs)
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list = None,
) -> t.Annotated[dict, namespace(
    structures=dynamic(orm.StructureData),
    remote_folders=dynamic(orm.RemoteData),
    parameters=dynamic(orm.Dict),
    trajectories=dynamic(orm.TrajectoryData),
    retrieved=dynamic(orm.FolderData)
)]:
    """
    Run single CP2K AIMD stage on all slabs in parallel using scatter-gather pattern.

    This function handles ONE temperature/timestep stage for all slabs.
    Call it multiple times sequentially to build multi-stage AIMD workflows.

    Args:
        slabs: Dictionary of slab structures to run AIMD on
        temperature: Target temperature in K
        steps: Number of MD steps for this stage
        code: CP2K code
        aimd_parameters: Base CP2K AIMD parameters (from get_aimd_defaults_cp2k)
        basis_file: SinglefileData containing BASIS_MOLOPT
        pseudo_file: SinglefileData containing GTH_POTENTIALS
        options: Scheduler options (metadata)
        clean_workdir: Whether to clean work directory
        restart_folders: Optional dict of RemoteData for restart (from previous stage)
        fixed_atoms_lists: Optional dict mapping slab_label -> list of fixed atom indices
        fix_components: Components to fix ("XYZ", "XY", "Z")
        fix_type: Optional - for dynamic calculation: 'bottom', 'top', 'center'
        fix_thickness: Optional - for dynamic calculation: thickness in Angstroms
        fix_elements: Optional - for dynamic calculation: list of element symbols

    Returns:
        Dictionary with outputs per slab:
            - structures: Output structures from this stage
            - remote_folders: RemoteData nodes for next stage restart
            - parameters: Output parameters from CP2K
            - trajectories: Trajectory data from MD
            - retrieved: Retrieved folder data
    """
    from teros.core.builders.aimd_builder_cp2k import prepare_aimd_parameters_cp2k
    from teros.core.fixed_atoms import add_fixed_atoms_to_cp2k_parameters, get_fixed_atoms_list

    # Get CP2K workchain
    Cp2kWorkChain = WorkflowFactory('cp2k.base')
    Cp2kTask = task(Cp2kWorkChain)

    structures_out = {}
    remote_folders_out = {}
    parameters_out = {}
    trajectories_out = {}
    retrieved_out = {}

    # If fixed_atoms_lists not provided but fix_type is, calculate dynamically
    computed_fixed_atoms = {}
    if fixed_atoms_lists is None and fix_type is not None and fix_thickness > 0:
        for slab_label, slab_structure in slabs.items():
            fixed_list = get_fixed_atoms_list(
                slab_structure,
                fix_type=str(fix_type) if fix_type else None,
                fix_thickness=float(fix_thickness) if fix_thickness else 0.0,
                fix_elements=list(fix_elements) if fix_elements else None,
            )
            computed_fixed_atoms[slab_label] = fixed_list
    elif fixed_atoms_lists:
        # Use provided fixed_atoms_lists, unwrapping any TaggedValues
        computed_fixed_atoms = {k: list(v) if v else [] for k, v in fixed_atoms_lists.items()}

    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Prepare parameters for this stage (unwrap TaggedValues with dict/float/int)
        stage_params = prepare_aimd_parameters_cp2k(dict(aimd_parameters), float(temperature), int(steps))

        # Add fixed atoms if available for this slab
        if slab_label in computed_fixed_atoms and computed_fixed_atoms[slab_label]:
            stage_params = add_fixed_atoms_to_cp2k_parameters(
                stage_params,
                computed_fixed_atoms[slab_label],
                str(fix_components) if fix_components else "XYZ",
            )

        # Build CP2K inputs
        cp2k_inputs = {
            'structure': slab_structure,
            'parameters': orm.Dict(dict=stage_params),
            'code': code,
            'metadata': {'options': dict(options)},  # Wrap scheduler options (dict() unwraps TaggedValue)
            'file': {
                'basis': basis_file,
                'pseudo': pseudo_file,
            },
            'settings': orm.Dict(dict={
                'additional_retrieve_list': [
                    "aiida-1.ener",
                    "aiida-1.restart",
                    "aiida-pos-1.xyz"
                ]
            }),
        }

        # Add restart folder if provided for this slab
        if restart_folders and slab_label in restart_folders:
            cp2k_inputs['parent_calc_folder'] = restart_folders[slab_label]

        # Create CP2K task
        aimd_task = Cp2kTask(
            cp2k=cp2k_inputs,
            max_iterations=orm.Int(3),
            clean_workdir=orm.Bool(bool(clean_workdir)),  # Unwrap TaggedValue
        )

        # Store CP2K outputs (note: different from VASP!)
        structures_out[slab_label] = aimd_task.output_structure
        remote_folders_out[slab_label] = aimd_task.remote_folder
        parameters_out[slab_label] = aimd_task.output_parameters
        trajectories_out[slab_label] = aimd_task.output_trajectory
        retrieved_out[slab_label] = aimd_task.retrieved

    # Gather: return collected results
    return {
        'structures': structures_out,
        'remote_folders': remote_folders_out,
        'parameters': parameters_out,
        'trajectories': trajectories_out,
        'retrieved': retrieved_out,
    }
