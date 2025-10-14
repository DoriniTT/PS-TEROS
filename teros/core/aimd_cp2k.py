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
    options: t.Annotated[dict, Cp2kBaseWorkChain.spec().inputs['metadata']],
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    fixed_atoms_lists: dict = None,
    fix_components: str = "XYZ",
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

    Returns:
        Dictionary with outputs per slab:
            - structures: Output structures from this stage
            - remote_folders: RemoteData nodes for next stage restart
            - parameters: Output parameters from CP2K
            - trajectories: Trajectory data from MD
            - retrieved: Retrieved folder data
    """
    from teros.core.builders.aimd_builder_cp2k import prepare_aimd_parameters_cp2k
    from teros.core.fixed_atoms import add_fixed_atoms_to_cp2k_parameters

    # Get CP2K workchain
    Cp2kWorkChain = WorkflowFactory('cp2k.base')
    Cp2kTask = task(Cp2kWorkChain)

    structures_out = {}
    remote_folders_out = {}
    parameters_out = {}
    trajectories_out = {}
    retrieved_out = {}

    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters_cp2k(aimd_parameters, temperature, steps)

        # Add fixed atoms if provided for this slab
        if fixed_atoms_lists and slab_label in fixed_atoms_lists:
            stage_params = add_fixed_atoms_to_cp2k_parameters(
                stage_params,
                fixed_atoms_lists[slab_label],
                fix_components,
            )

        # Build CP2K inputs
        cp2k_inputs = {
            'structure': slab_structure,
            'parameters': orm.Dict(dict=stage_params),
            'code': code,
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

        # Create CP2K task - pass metadata and cp2k separately
        aimd_task = Cp2kTask(
            cp2k=cp2k_inputs,
            metadata=dict(options),  # top-level metadata for the workchain
            max_iterations=orm.Int(3),
            clean_workdir=orm.Bool(clean_workdir),
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
