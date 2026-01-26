"""
Stage 9: AIMD Stage.

This module provides functions for adding Ab Initio Molecular Dynamics (AIMD)
calculations to the workgraph. It supports both VASP and CP2K calculators with
multi-stage sequential AIMD simulations.

The main functions are:
- ``resolve_aimd_parameters``: Resolve AIMD parameters with fallback chain.
- ``prepare_fixed_atoms``: Prepare fixed atoms lists for each slab structure.
- ``determine_initial_slabs_source``: Determine the source of initial structures.
- ``add_aimd_stage``: Add the complete AIMD stage to the workgraph.

These functions implement the AIMD stage (Stage 9) of the
``build_core_workgraph()`` function.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Union

from aiida import orm

if TYPE_CHECKING:
    from aiida_workgraph import WorkGraph

__all__ = [
    "resolve_aimd_parameters",
    "prepare_fixed_atoms",
    "determine_initial_slabs_source",
    "add_aimd_stage",
]

logger = logging.getLogger(__name__)


def resolve_aimd_parameters(
    aimd_parameters: Optional[Dict[str, Any]],
    aimd_options: Optional[Dict[str, Any]],
    aimd_potential_mapping: Optional[Dict[str, str]],
    aimd_kpoints_spacing: Optional[float],
    slab_parameters: Optional[Dict[str, Any]],
    slab_options: Optional[Dict[str, Any]],
    slab_potential_mapping: Optional[Dict[str, str]],
    slab_kpoints_spacing: Optional[float],
    bulk_parameters: Dict[str, Any],
    bulk_options: Dict[str, Any],
    bulk_potential_mapping: Dict[str, str],
    kpoints_spacing: float,
) -> tuple:
    """Resolve AIMD parameters with fallback chain.

    Applies the fallback hierarchy: aimd -> slab -> bulk for all parameter types.

    Args:
        aimd_parameters: AIMD-specific INCAR parameters or None.
        aimd_options: AIMD-specific scheduler options or None.
        aimd_potential_mapping: AIMD-specific potential mapping or None.
        aimd_kpoints_spacing: AIMD-specific k-points spacing or None.
        slab_parameters: Slab-specific INCAR parameters or None.
        slab_options: Slab-specific scheduler options or None.
        slab_potential_mapping: Slab-specific potential mapping or None.
        slab_kpoints_spacing: Slab-specific k-points spacing or None.
        bulk_parameters: Bulk INCAR parameters (required).
        bulk_options: Bulk scheduler options (required).
        bulk_potential_mapping: Bulk potential mapping (required).
        kpoints_spacing: Default k-points spacing (required).

    Returns:
        Tuple of (params, opts, pot_map, kpts) with resolved values:
            - params: Resolved INCAR parameters dict
            - opts: Resolved scheduler options dict
            - pot_map: Resolved potential mapping dict
            - kpts: Resolved k-points spacing float

    Example:
        >>> params, opts, pot_map, kpts = resolve_aimd_parameters(
        ...     aimd_parameters=None,
        ...     aimd_options={'walltime': 48*3600},
        ...     aimd_potential_mapping=None,
        ...     aimd_kpoints_spacing=0.1,
        ...     slab_parameters={'incar': {'ENCUT': 400}},
        ...     slab_options=None,
        ...     slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     slab_kpoints_spacing=0.04,
        ...     bulk_parameters={'incar': {'ENCUT': 520}},
        ...     bulk_options={'resources': {'num_machines': 1}},
        ...     bulk_potential_mapping={'Ag': 'Ag'},
        ...     kpoints_spacing=0.03,
        ... )
        >>> kpts
        0.1
    """
    # First resolve slab parameters from bulk if not provided
    slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
    slab_opts = slab_options if slab_options is not None else bulk_options
    slab_pot_map = (
        slab_potential_mapping
        if slab_potential_mapping is not None
        else bulk_potential_mapping
    )
    slab_kpts = (
        slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing
    )

    # Then resolve AIMD parameters from slab
    params = aimd_parameters if aimd_parameters is not None else slab_params
    opts = aimd_options if aimd_options is not None else slab_opts
    pot_map = (
        aimd_potential_mapping if aimd_potential_mapping is not None else slab_pot_map
    )
    kpts = aimd_kpoints_spacing if aimd_kpoints_spacing is not None else slab_kpts

    return params, opts, pot_map, kpts


def prepare_fixed_atoms(
    input_slabs: Optional[Dict[str, Any]],
    fix_atoms: bool,
    fix_type: Optional[str],
    fix_thickness: float,
    fix_elements: Optional[List[str]],
    log: logging.Logger = logger,
) -> Dict[str, List[int]]:
    """Prepare fixed atoms lists for each slab structure.

    Calculates the list of atoms to fix for each input slab based on the
    specified fixing criteria. Only applicable when input_slabs is provided.

    Args:
        input_slabs: Dictionary of input slab structures (label -> StructureData).
            If None, returns empty dict (fixed atoms will be calculated dynamically).
        fix_atoms: Whether atom fixing is enabled.
        fix_type: Type of fixing ('bottom', 'top', 'center', or None).
        fix_thickness: Thickness in Angstroms for the fixed region.
        fix_elements: Optional list of element symbols to fix. If None, all
            elements in the region are fixed.
        log: Logger instance for messages.

    Returns:
        Dictionary mapping slab labels to lists of 1-based atom indices to fix.
        Empty dict if fix_atoms is False, fix_type is None, or input_slabs is None.

    Example:
        >>> fixed_atoms = prepare_fixed_atoms(
        ...     input_slabs={'term_0': slab_0, 'term_1': slab_1},
        ...     fix_atoms=True,
        ...     fix_type='bottom',
        ...     fix_thickness=7.0,
        ...     fix_elements=['Ag'],
        ... )
        >>> fixed_atoms['term_0']
        [1, 2, 3, 4, 5]
    """
    fixed_atoms_lists: Dict[str, List[int]] = {}

    if not fix_atoms or fix_type is None:
        return fixed_atoms_lists

    log.info("  -> Preparing fixed atoms constraints")
    log.info("     Type: %s", fix_type)
    log.info("     Thickness: %s A", fix_thickness)
    log.info("     Elements: %s", fix_elements if fix_elements else "all")

    # For input_slabs, calculate fixed atoms now (pre-computed)
    if input_slabs is not None:
        from teros.core.fixed_atoms import get_fixed_atoms_list

        for label, slab_struct in input_slabs.items():
            fixed_list = get_fixed_atoms_list(
                slab_struct,
                fix_type=fix_type,
                fix_thickness=fix_thickness,
                fix_elements=fix_elements,
            )
            fixed_atoms_lists[label] = fixed_list
            log.info("     %s: %d atoms fixed", label, len(fixed_list))

    return fixed_atoms_lists


def determine_initial_slabs_source(
    wg: WorkGraph,
    relax_slabs: bool,
    input_slabs: Optional[Dict[str, Any]],
    log: logging.Logger = logger,
) -> Any:
    """Determine the source of initial structures for AIMD.

    Finds the appropriate source of slab structures based on workflow configuration.
    Priority order:
    1. Relaxed slabs from relax_slabs_scatter task (if relax_slabs=True)
    2. Relaxed slabs from collect_slab_outputs_restart task (restart scenario)
    3. Unrelaxed input slabs (if input_slabs provided)
    4. Generated slabs from generate_slab_structures task (auto-generation)

    Args:
        wg: WorkGraph instance being built.
        relax_slabs: Whether slab relaxation is enabled in the workflow.
        input_slabs: Dictionary of user-provided input slab structures or None.
        log: Logger instance for messages.

    Returns:
        Either a dict of StructureData (for input_slabs) or a task output socket
        that will be connected to the AIMD task inputs.

    Note:
        The returned value can be:
        - A plain dict (for input_slabs case)
        - A task output socket (relax_task.outputs.relaxed_structures, etc.)

    Example:
        >>> source = determine_initial_slabs_source(
        ...     wg=wg,
        ...     relax_slabs=True,
        ...     input_slabs=None,
        ... )
        >>> # source is now relax_slabs_scatter.outputs.relaxed_structures
    """
    if relax_slabs and "relax_slabs_scatter" in wg.tasks:
        relax_task = wg.tasks["relax_slabs_scatter"]
        source = relax_task.outputs.relaxed_structures
        log.info("     Using relaxed slabs from relax_slabs_scatter task")
    elif relax_slabs and "collect_slab_outputs_restart" in wg.tasks:
        collector_task = wg.tasks["collect_slab_outputs_restart"]
        source = collector_task.outputs.structures
        log.info("     Using relaxed slabs from restart collector")
    elif input_slabs is not None:
        source = input_slabs
        log.info("     Using unrelaxed input slabs as initial structures")
    else:
        gen_task = wg.tasks["generate_slab_structures"]
        source = gen_task.outputs.slabs
        log.info("     Using unrelaxed generated slabs as initial structures")

    return source


def add_aimd_stage(
    wg: WorkGraph,
    calculator: str,
    aimd_sequence: List[Dict[str, Any]],
    code_label: str,
    aimd_code_label: Optional[str],
    aimd_parameters: Dict[str, Any],
    aimd_options: Dict[str, Any],
    aimd_potential_mapping: Dict[str, str],
    aimd_kpoints_spacing: float,
    potential_family: str,
    clean_workdir: bool,
    max_concurrent_jobs: Optional[int],
    fix_atoms: bool,
    fix_type: Optional[str],
    fix_thickness: float,
    fix_elements: Optional[List[str]],
    fix_components: str,
    aimd_supercell: Optional[List[int]],
    relax_slabs: bool,
    input_slabs: Optional[Dict[str, Any]],
    basis_file: Optional[Any] = None,
    pseudo_file: Optional[Any] = None,
    log: logging.Logger = logger,
) -> List[Any]:
    """Add the complete AIMD stage to the workgraph.

    This is the main function for Stage 9. It handles:
    1. Determining the calculator type (VASP or CP2K)
    2. Loading the appropriate code
    3. Handling fixed atoms if requested
    4. Determining the initial slab source
    5. Creating supercells if requested
    6. Creating sequential AIMD stage tasks
    7. Wiring outputs between stages

    Args:
        wg: WorkGraph instance being built.
        calculator: Calculator type ('vasp' or 'cp2k').
        aimd_sequence: List of stage configurations. Each dict contains:
            - For VASP: 'TEBEG' (temperature), 'NSW' (steps), optional 'TEEND', 'POTIM'
            - For CP2K: 'temperature', 'steps'
        code_label: Main code label (fallback if aimd_code_label not provided).
        aimd_code_label: AIMD-specific code label or None.
        aimd_parameters: Resolved INCAR/CP2K parameters dict.
        aimd_options: Resolved scheduler options dict.
        aimd_potential_mapping: Resolved potential mapping dict.
        aimd_kpoints_spacing: Resolved k-points spacing.
        potential_family: POTCAR family name (e.g., 'PBE').
        clean_workdir: Whether to clean work directories.
        max_concurrent_jobs: Maximum parallel AIMD calculations or None.
        fix_atoms: Whether atom fixing is enabled.
        fix_type: Type of fixing ('bottom', 'top', 'center', or None).
        fix_thickness: Thickness in Angstroms for the fixed region.
        fix_elements: Optional list of element symbols to fix.
        fix_components: Components to fix ('XYZ', 'XY', 'Z').
        aimd_supercell: Supercell dimensions [nx, ny, nz] or None.
        relax_slabs: Whether slab relaxation is enabled.
        input_slabs: Dictionary of input slab structures or None.
        basis_file: CP2K basis file (SinglefileData) or None. Required for CP2K.
        pseudo_file: CP2K pseudopotential file (SinglefileData) or None. Required for CP2K.
        log: Logger instance for messages.

    Returns:
        List of stage task objects for reference.

    Raises:
        ValueError: If calculator is not 'vasp' or 'cp2k'.

    Note:
        - For VASP, uses aimd_single_stage_scatter from teros.core.aimd_functions
        - For CP2K, uses aimd_single_stage_scatter_cp2k from teros.core.aimd_cp2k
        - Stages are wired sequentially: outputs of stage N feed into stage N+1
        - Access outputs via: wg.tasks['aimd_stage_XX_XXXK'].outputs

    Example:
        >>> stage_tasks = add_aimd_stage(
        ...     wg=wg,
        ...     calculator='vasp',
        ...     aimd_sequence=[
        ...         {'TEBEG': 300, 'NSW': 1000},
        ...         {'TEBEG': 500, 'NSW': 2000},
        ...     ],
        ...     code_label='VASP-6.5.1@cluster',
        ...     aimd_code_label=None,
        ...     aimd_parameters={'incar': {'IBRION': 0, 'MDALGO': 2}},
        ...     ...
        ... )
        >>> len(stage_tasks)
        2
    """
    from aiida.orm import load_code

    log.info("  -> Adding AIMD stages to workflow")
    log.info("     Calculator: %s", calculator)
    log.info("     Number of stages: %d", len(aimd_sequence))

    # Import appropriate scatter function based on calculator
    if calculator == "vasp":
        from teros.core.aimd_functions import aimd_single_stage_scatter

        aimd_scatter_func = aimd_single_stage_scatter
    elif calculator == "cp2k":
        from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k

        aimd_scatter_func = aimd_single_stage_scatter_cp2k
    else:
        raise ValueError(f"Unknown calculator: {calculator}")

    # Load code - use aimd_code_label if provided, otherwise use code_label
    aimd_code_to_use = aimd_code_label if aimd_code_label is not None else code_label
    code = load_code(aimd_code_to_use)
    log.info("     Using code: %s", aimd_code_to_use)

    # Handle fixed atoms if requested
    fixed_atoms_lists = prepare_fixed_atoms(
        input_slabs=input_slabs,
        fix_atoms=fix_atoms,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        fix_elements=fix_elements,
        log=log,
    )

    # Determine which slabs to use as initial structures
    initial_slabs_source = determine_initial_slabs_source(
        wg=wg,
        relax_slabs=relax_slabs,
        input_slabs=input_slabs,
        log=log,
    )

    # Handle supercell creation
    if aimd_supercell is not None:
        log.info("     Supercell: %s", aimd_supercell)
        from teros.core.aimd.tasks import create_supercells_scatter

        # Add supercell task
        sc_task = wg.add_task(
            create_supercells_scatter,
            name="create_aimd_supercells",
            slabs=initial_slabs_source,
            spec=orm.List(list=aimd_supercell),
        )
        initial_slabs_source = sc_task.outputs.result
        log.info("     Created supercell generation task")

    # Sequential AIMD stages
    current_structures = initial_slabs_source
    current_remotes: Dict[str, Any] = {}
    stage_tasks: List[Any] = []

    for stage_idx, stage_config in enumerate(aimd_sequence):
        # Get temperature for stage naming (support both old and new format)
        stage_temp = stage_config.get("TEBEG", stage_config.get("temperature", 0))
        stage_name = f"aimd_stage_{stage_idx:02d}_{stage_temp}K"

        # Print stage info (support both old and new format)
        stage_steps = stage_config.get("NSW", stage_config.get("steps", 0))
        log.info("     Stage %d: %sK x %s steps", stage_idx, stage_temp, stage_steps)

        # Build stage inputs based on calculator
        if calculator == "vasp":
            stage_inputs = {
                "slabs": current_structures,
                "stage_config": stage_config,
                "code": code,
                "base_aimd_parameters": aimd_parameters,
                "structure_aimd_overrides": None,
                "potential_family": potential_family,
                "potential_mapping": aimd_potential_mapping,
                "options": aimd_options,
                "kpoints_spacing": aimd_kpoints_spacing,
                "clean_workdir": clean_workdir,
                "restart_folders": current_remotes,
                "max_number_jobs": (
                    orm.Int(max_concurrent_jobs)
                    if max_concurrent_jobs is not None
                    else None
                ),
            }
        elif calculator == "cp2k":
            stage_inputs = {
                "slabs": current_structures,
                "temperature": stage_config["temperature"],
                "steps": stage_config["steps"],
                "code": code,
                "base_aimd_parameters": aimd_parameters,
                "structure_aimd_overrides": None,
                "basis_file": basis_file,
                "pseudo_file": pseudo_file,
                "options": aimd_options,  # Pass scheduler options directly
                "clean_workdir": clean_workdir,
                "restart_folders": current_remotes,
                "max_number_jobs": (
                    orm.Int(max_concurrent_jobs)
                    if max_concurrent_jobs is not None
                    else None
                ),
            }

            # Add fixed atoms for CP2K
            if fixed_atoms_lists:
                # Pre-computed fixed atoms (for input_slabs)
                stage_inputs["fixed_atoms_lists"] = fixed_atoms_lists
                stage_inputs["fix_components"] = fix_components
            elif fix_atoms and fix_type is not None:
                # Dynamic fixed atoms calculation (for auto-generated slabs)
                stage_inputs["fix_type"] = fix_type
                stage_inputs["fix_thickness"] = fix_thickness if fix_thickness else 0.0
                stage_inputs["fix_elements"] = fix_elements
                stage_inputs["fix_components"] = fix_components

        # Add task
        stage_task = wg.add_task(
            aimd_scatter_func,
            name=stage_name,
            **stage_inputs,
        )

        stage_tasks.append(stage_task)

        # Wire outputs to next stage inputs
        current_structures = stage_task.outputs.structures
        current_remotes = stage_task.outputs.remote_folders

    log.info("  AIMD calculation enabled (%d sequential stages)", len(aimd_sequence))
    log.info("     Access AIMD outputs via: wg.tasks['aimd_stage_XX_XXXK'].outputs")

    return stage_tasks
