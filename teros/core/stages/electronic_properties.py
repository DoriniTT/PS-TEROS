"""Stage 7-8: Electronic properties calculations for bulk and slabs.

This module provides functions for adding electronic properties (DOS and bands)
calculations to the workgraph. These stages are optional and are enabled by
the `compute_electronic_properties_bulk` and `compute_electronic_properties_slabs`
flags respectively.

Stage 7 handles slab electronic properties via scatter-gather pattern, supporting
per-slab parameter overrides. Stage 8 handles bulk electronic properties using
the BandsWorkChain.

Functions:
    add_slab_electronic_properties_stage: Add slab DOS/bands calculations.
    add_bulk_electronic_properties_stage: Add bulk DOS/bands calculations.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from aiida import orm
    from aiida_workgraph import WorkGraph

__all__ = [
    "add_slab_electronic_properties_stage",
    "add_bulk_electronic_properties_stage",
]


def add_slab_electronic_properties_stage(
    wg: WorkGraph,
    code_label: str,
    potential_family: str,
    bulk_options: Dict[str, Any],
    slab_potential_mapping: Optional[Dict[str, str]],
    slab_bands_parameters: Optional[Dict[str, Any]],
    slab_bands_options: Optional[Dict[str, Any]],
    slab_band_settings: Optional[Dict[str, Any]],
    slab_electronic_properties: Dict[str, Dict[str, Any]],
    max_concurrent_jobs: Optional[int],
    use_restart_mode: bool,
    logger: logging.Logger,
    scatter_task: Optional[Any] = None,
    collector: Optional[Any] = None,
) -> None:
    """Add slab electronic properties scatter task to the workgraph.

    This function adds a scatter-gather task that calculates DOS and bands
    for selected slab terminations. It supports per-slab parameter overrides
    through the `slab_electronic_properties` dictionary.

    The source of relaxed slab structures depends on the workflow mode:
    - Restart mode: Uses collector.outputs.structures
    - Non-restart input_slabs: Uses scatter_task.outputs.relaxed_structures

    Args:
        wg: WorkGraph instance to add the task to.
        code_label: AiiDA code label for VASP calculations.
        potential_family: POTCAR family name (e.g., 'PBE', 'PBE.54').
        bulk_options: Default scheduler options (fallback when slab_bands_options
            is None).
        slab_potential_mapping: Element to potential mapping for slabs. If None,
            must be provided separately.
        slab_bands_parameters: VASP parameters for slab DOS/bands calculations.
            Can contain 'scf', 'bands', 'dos', and 'scf_kpoints_distance' keys.
        slab_bands_options: Scheduler options for slab bands calculations. Falls
            back to bulk_options if None.
        slab_band_settings: Band workflow settings (e.g., seekpath parameters).
        slab_electronic_properties: Dictionary mapping slab labels to parameter
            configs. Each config can contain:
            - 'bands_parameters': Dict with 'scf', 'bands', 'dos' keys
            - 'bands_options': Scheduler options
            - 'band_settings': Band workflow settings
            Only slabs present in this dictionary will have electronic properties
            calculated.
        max_concurrent_jobs: Maximum number of concurrent VASP calculations.
            None means unlimited.
        use_restart_mode: Whether restart mode is being used. Determines the
            source of relaxed slab structures.
        logger: Logger instance for status messages.
        scatter_task: The scatter task from normal input_slabs mode. Used when
            use_restart_mode is False.
        collector: The collector task from restart mode. Used when
            use_restart_mode is True.

    Returns:
        None. Modifies wg in place by adding:
        - Task: 'calculate_electronic_properties_slabs'
        - Outputs: slab_bands, slab_dos, slab_primitive_structures,
          slab_seekpath_parameters

    Raises:
        ValueError: If neither scatter_task nor collector is provided based on
            the restart mode.

    Example:
        >>> add_slab_electronic_properties_stage(
        ...     wg=workgraph,
        ...     code_label='VASP-6.5.1@cluster',
        ...     potential_family='PBE',
        ...     bulk_options={'resources': {'num_machines': 1}},
        ...     slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     slab_bands_parameters={'scf': {'ENCUT': 520}},
        ...     slab_bands_options=None,
        ...     slab_band_settings=None,
        ...     slab_electronic_properties={'term_0': {}, 'term_2': {}},
        ...     max_concurrent_jobs=4,
        ...     use_restart_mode=False,
        ...     logger=logger,
        ...     scatter_task=scatter_task,
        ...     collector=None,
        ... )
    """
    from aiida.orm import load_code
    from teros.core.slabs import calculate_electronic_properties_slabs_scatter

    # Get default parameters (per-slab overrides handled in scatter function)
    default_params = slab_bands_parameters if slab_bands_parameters else {}
    default_opts = slab_bands_options if slab_bands_options else bulk_options
    default_settings = slab_band_settings if slab_band_settings else {}

    # Get code
    code = load_code(code_label)

    # Determine which slabs output to use based on restart mode
    if use_restart_mode:
        # Restart mode: use collector outputs
        if collector is None:
            raise ValueError("collector must be provided when use_restart_mode is True")
        relaxed_slabs_source = collector.outputs.structures
    else:
        # Non-restart input_slabs: use scatter task outputs
        if scatter_task is None:
            raise ValueError(
                "scatter_task must be provided when use_restart_mode is False"
            )
        relaxed_slabs_source = scatter_task.outputs.relaxed_structures

    # Add electronic properties task
    slab_elec_task = wg.add_task(
        calculate_electronic_properties_slabs_scatter,
        name="calculate_electronic_properties_slabs",
        slabs=relaxed_slabs_source,
        slab_electronic_properties=slab_electronic_properties,
        code=code,
        potential_family=potential_family,
        potential_mapping=slab_potential_mapping,
        clean_workdir=False,  # Default to preserve remote folders
        default_bands_parameters=default_params,
        default_bands_options=default_opts,
        default_band_settings=default_settings,
        max_number_jobs=max_concurrent_jobs,
    )

    # Connect outputs
    wg.outputs.slab_bands = slab_elec_task.outputs.slab_bands
    wg.outputs.slab_dos = slab_elec_task.outputs.slab_dos
    wg.outputs.slab_primitive_structures = (
        slab_elec_task.outputs.slab_primitive_structures
    )
    wg.outputs.slab_seekpath_parameters = (
        slab_elec_task.outputs.slab_seekpath_parameters
    )

    logger.info(
        "  ✓ Slab electronic properties calculation enabled for %d slabs",
        len(slab_electronic_properties),
    )


def add_bulk_electronic_properties_stage(
    wg: WorkGraph,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: Dict[str, str],
    bulk_options: Dict[str, Any],
    bands_parameters: Optional[Dict[str, Any]],
    bands_options: Optional[Dict[str, Any]],
    band_settings: Optional[Dict[str, Any]],
    clean_workdir: bool,
    logger: logging.Logger,
) -> None:
    """Add bulk electronic properties (DOS/bands) calculation to the workgraph.

    This function adds a BandsWorkChain task that calculates the electronic
    density of states (DOS) and band structure for the relaxed bulk structure.
    The bulk structure is taken from the output of the bulk VaspWorkChain task.

    The BandsWorkChain supports three calculation namespaces:
    - scf: Self-consistent field calculation (required)
    - bands: Band structure calculation (optional)
    - dos: Density of states calculation (optional)

    Args:
        wg: WorkGraph instance to add the task to.
        code_label: AiiDA code label for VASP calculations.
        potential_family: POTCAR family name (e.g., 'PBE', 'PBE.54').
        bulk_potential_mapping: Element to potential mapping for bulk calculations.
        bulk_options: Default scheduler options (fallback when bands_options
            is None).
        bands_parameters: VASP parameters for DOS/bands calculations. Can contain:
            - 'scf': Dict with INCAR parameters for SCF calculation
            - 'bands': Dict with INCAR parameters for bands calculation
            - 'dos': Dict with INCAR parameters for DOS calculation
            - 'scf_kpoints_distance': K-points spacing for SCF (optional)
        bands_options: Scheduler options for bands calculations. Falls back
            to bulk_options if None.
        band_settings: Band workflow settings (e.g., seekpath parameters for
            high-symmetry path generation).
        clean_workdir: Whether to clean remote directories after calculations.
        logger: Logger instance for status messages.

    Returns:
        None. Modifies wg in place by adding:
        - Task: 'BandsWorkChain_bulk'
        - Outputs: bulk_bands, bulk_dos, bulk_primitive_structure,
          bulk_seekpath_parameters

    Note:
        This function assumes a task named 'VaspWorkChain' exists in the
        workgraph and provides the relaxed bulk structure via its 'structure'
        output.

    Example:
        >>> add_bulk_electronic_properties_stage(
        ...     wg=workgraph,
        ...     code_label='VASP-6.5.1@cluster',
        ...     potential_family='PBE',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     bulk_options={'resources': {'num_machines': 1}},
        ...     bands_parameters={
        ...         'scf': {'ENCUT': 520, 'ISMEAR': -5},
        ...         'bands': {'ENCUT': 520, 'ISMEAR': 0, 'LORBIT': 11},
        ...         'dos': {'ENCUT': 520, 'ISMEAR': -5, 'NEDOS': 3001},
        ...     },
        ...     bands_options=None,
        ...     band_settings=None,
        ...     clean_workdir=False,
        ...     logger=logger,
        ... )
    """
    from aiida.plugins import WorkflowFactory
    from aiida.orm import load_code
    from aiida_workgraph import task as wg_task

    # Get BandsWorkChain and wrap it as a task (same pattern as VaspWorkChain)
    BandsWorkChain = WorkflowFactory("vasp.v2.bands")
    BandsTask = wg_task(BandsWorkChain)

    # Get bulk task from workgraph
    bulk_vasp_task = wg.tasks["VaspWorkChain"]

    # Get code
    code = load_code(code_label)

    # Use bands-specific options or fall back to bulk options
    bands_opts = bands_options if bands_options is not None else bulk_options

    # Build inputs dictionary
    # BandsWorkChain expects inputs in namespaces: scf, bands, dos
    bands_inputs: Dict[str, Any] = {
        "structure": bulk_vasp_task.outputs.structure,  # Socket from bulk task
        "metadata": {
            "label": "Bulk_Electronic_Properties",
            "description": "DOS and bands calculation for relaxed bulk structure",
        },
    }

    # Add band_settings if provided
    if band_settings:
        bands_inputs["band_settings"] = band_settings

    # Build SCF namespace inputs (required)
    scf_inputs: Dict[str, Any] = {
        "code": code,
        "potential_family": potential_family,
        "potential_mapping": bulk_potential_mapping,
        "options": bands_opts,
        "clean_workdir": clean_workdir,
    }
    if bands_parameters and "scf" in bands_parameters:
        scf_inputs["parameters"] = {"incar": bands_parameters["scf"]}
    if bands_parameters and "scf_kpoints_distance" in bands_parameters:
        scf_inputs["kpoints_spacing"] = bands_parameters["scf_kpoints_distance"]

    bands_inputs["scf"] = scf_inputs

    # Build Bands namespace inputs (optional)
    if bands_parameters and "bands" in bands_parameters:
        bands_inputs["bands"] = {
            "code": code,
            "potential_family": potential_family,
            "potential_mapping": bulk_potential_mapping,
            "options": bands_opts,
            "clean_workdir": clean_workdir,
            "parameters": {"incar": bands_parameters["bands"]},
        }

    # Build DOS namespace inputs (optional)
    if bands_parameters and "dos" in bands_parameters:
        bands_inputs["dos"] = {
            "code": code,
            "potential_family": potential_family,
            "potential_mapping": bulk_potential_mapping,
            "options": bands_opts,
            "clean_workdir": clean_workdir,
            "parameters": {"incar": bands_parameters["dos"]},
        }

    # Add the bands task to the workgraph
    bands_task = wg.add_task(BandsTask, name="BandsWorkChain_bulk", **bands_inputs)

    # Connect all outputs from BandsWorkChain
    wg.outputs.bulk_bands = bands_task.outputs.band_structure
    wg.outputs.bulk_dos = bands_task.outputs.dos
    wg.outputs.bulk_primitive_structure = bands_task.outputs.primitive_structure
    wg.outputs.bulk_seekpath_parameters = bands_task.outputs.seekpath_parameters

    logger.info("  ✓ Electronic properties calculation enabled (DOS and bands)")
