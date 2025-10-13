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


@task.graph
def aimd_sequential_slab(
    slab_label: str,
    structure: orm.StructureData,
    aimd_sequence: list,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
) -> t.Annotated[dict, dynamic(namespace())]:
    """
    Run sequential AIMD stages for a single slab termination.

    Creates a chain of VASP calculations where each stage:
    1. Uses structure from previous stage (or initial for first)
    2. Uses restart_folder from previous stage for WAVECAR/CHGCAR
    3. Runs MD at specified temperature for specified steps

    Args:
        slab_label: Label for this slab (e.g., 'term_0')
        structure: Initial slab structure
        aimd_sequence: List of dicts with 'temperature' and 'steps' keys
        code: VASP code
        aimd_parameters: Base AIMD INCAR parameters
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory

    Returns:
        Dictionary with flattened outputs for all stages:
            - stage_00_300K_structure, stage_00_300K_trajectory, etc.
            - final_structure, final_remote, final_trajectory
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    prev_structure = structure
    prev_remote = None
    all_outputs = {}

    # Loop through AIMD sequence
    for stage_idx, stage_config in enumerate(aimd_sequence):
        temp = stage_config['temperature']
        steps = stage_config['steps']

        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters(
            aimd_parameters, temp, steps
        )

        # Build VASP inputs
        vasp_inputs = {
            'structure': prev_structure,
            'code': code,
            'parameters': {'incar': stage_params},
            'options': options,
            'potential_family': potential_family,
            'potential_mapping': potential_mapping,
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }

        # Add restart folder if available
        if prev_remote is not None:
            vasp_inputs['restart_folder'] = prev_remote

        # Create VASP task
        aimd_task = VaspTask(**vasp_inputs)

        # Store outputs with descriptive keys
        stage_prefix = f"stage_{stage_idx:02d}_{temp:03.0f}K"
        all_outputs[f"{stage_prefix}_structure"] = aimd_task.structure
        all_outputs[f"{stage_prefix}_trajectory"] = aimd_task.trajectory
        all_outputs[f"{stage_prefix}_energy"] = extract_total_energy(
            energies=aimd_task.misc
        ).result
        all_outputs[f"{stage_prefix}_remote"] = aimd_task.remote_folder
        all_outputs[f"{stage_prefix}_retrieved"] = aimd_task.retrieved

        # Update for next stage
        prev_structure = aimd_task.structure
        prev_remote = aimd_task.remote_folder

    # Add final outputs
    all_outputs['final_structure'] = prev_structure
    all_outputs['final_remote'] = prev_remote
    all_outputs['final_trajectory'] = aimd_task.trajectory

    return all_outputs


@task.graph
def aimd_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    aimd_sequence: list,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
) -> t.Annotated[dict, namespace(**{})]:
    """
    Run AIMD on all slabs in parallel using scatter-gather pattern.

    Each slab gets the same AIMD sequence but runs independently in parallel.

    Args:
        slabs: Dictionary of slab structures (from slab generation)
        aimd_sequence: List of {'temperature': K, 'steps': N} dicts
        code: VASP code
        aimd_parameters: Base AIMD INCAR parameters
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory

    Returns:
        Dictionary with AIMD results for each slab:
            - term_0: {stage_00_300K_structure, ..., final_structure}
            - term_1: {stage_00_300K_structure, ..., final_structure}
            - ...
    """
    aimd_outputs = {}

    for slab_label, slab_structure in slabs.items():
        # Run sequential AIMD for this slab
        slab_aimd_result = aimd_sequential_slab(
            slab_label=slab_label,
            structure=slab_structure,
            aimd_sequence=aimd_sequence,
            code=code,
            aimd_parameters=aimd_parameters,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            kpoints_spacing=kpoints_spacing,
            clean_workdir=clean_workdir,
        )

        # Store outputs for this slab
        aimd_outputs[slab_label] = slab_aimd_result

    return aimd_outputs
