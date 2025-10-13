"""
AIMD Module for PS-TEROS

Ab initio molecular dynamics calculations on slab structures.
Sequential AIMD stages with automatic restart chaining.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace
from teros.core.slabs import extract_total_energy


def get_settings():
    """
    Parser settings for aiida-vasp (copied from workgraph.py for consistency).
    """
    return {
        'parser_settings': {
            'add_energy': True,
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }


def prepare_aimd_parameters(
    base_parameters: dict,
    temperature: float,
    steps: int
) -> dict:
    """
    Prepare INCAR parameters for a single AIMD stage.

    Takes base AIMD parameters and injects stage-specific values:
    - TEBEG and TEEND (both set to temperature for isothermal)
    - NSW (number of MD steps for this stage)

    Args:
        base_parameters: Base AIMD INCAR dict (IBRION=0, MDALGO, POTIM, etc.)
        temperature: Target temperature in K
        steps: Number of MD steps

    Returns:
        Complete INCAR dict for this stage
    """
    params = base_parameters.copy()
    params['TEBEG'] = temperature
    params['TEEND'] = temperature  # Isothermal
    params['NSW'] = steps
    return params


def _aimd_single_stage_scatter_impl(
    slabs: dict,
    temperature: float,
    steps: int,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
    restart_folders: dict = None,
):
    """
    Internal implementation of AIMD single stage scatter.
    This is the actual logic without @task.graph decorator.
    """
    from aiida.plugins import WorkflowFactory
    
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    structures_out = {}
    remote_folders_out = {}
    energies_out = {}

    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters(aimd_parameters, temperature, steps)

        # Build VASP inputs
        vasp_inputs = {
            'structure': slab_structure,
            'code': code,
            'parameters': {'incar': stage_params},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }

        # Add restart folder if provided for this slab
        if restart_folders and slab_label in restart_folders:
            vasp_inputs['restart_folder'] = restart_folders[slab_label]

        # Create VASP task
        aimd_task = VaspTask(**vasp_inputs)

        # Store outputs for this slab
        structures_out[slab_label] = aimd_task.structure
        remote_folders_out[slab_label] = aimd_task.remote_folder
        energies_out[slab_label] = extract_total_energy(energies=aimd_task.misc).result

    # Gather: return collected results
    return {
        'structures': structures_out,
        'remote_folders': remote_folders_out,
        'energies': energies_out,
    }


@task.graph
def aimd_single_stage_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    temperature: float,
    steps: int,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
) -> t.Annotated[dict, namespace(structures=dynamic(orm.StructureData), remote_folders=dynamic(orm.RemoteData), energies=dynamic(orm.Float))]:
    """
    Run single AIMD stage on all slabs in parallel using scatter-gather pattern.

    This function handles ONE temperature/timestep stage for all slabs.
    Call it multiple times sequentially to build multi-stage AIMD workflows.

    Args:
        slabs: Dictionary of slab structures to run AIMD on
        temperature: Target temperature in K
        steps: Number of MD steps for this stage
        code: VASP code
        aimd_parameters: Base AIMD INCAR parameters (IBRION=0, MDALGO, etc.)
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory
        restart_folders: Optional dict of RemoteData for restart (from previous stage)

    Returns:
        Dictionary with outputs per slab:
            - structures: Output structures from this stage
            - remote_folders: RemoteData nodes for potential next stage restart
            - energies: Total energies from this stage
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    structures_out = {}
    remote_folders_out = {}
    energies_out = {}

    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters(aimd_parameters, temperature, steps)

        # Build VASP inputs
        vasp_inputs = {
            'structure': slab_structure,
            'code': code,
            'parameters': {'incar': stage_params},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }

        # Add restart folder if provided for this slab
        if restart_folders and slab_label in restart_folders:
            vasp_inputs['restart_folder'] = restart_folders[slab_label]

        # Create VASP task
        aimd_task = VaspTask(**vasp_inputs)

        # Store outputs for this slab
        structures_out[slab_label] = aimd_task.structure
        remote_folders_out[slab_label] = aimd_task.remote_folder
        energies_out[slab_label] = extract_total_energy(energies=aimd_task.misc).result

    # Gather: return collected results
    return {
        'structures': structures_out,
        'remote_folders': remote_folders_out,
        'energies': energies_out,
    }
