"""
Stage 10: Adsorption Energy Stage.

This module provides functions for adding adsorption energy calculations
to the PS-TEROS core workgraph.

The main functions are:
- ``resolve_adsorption_parameters``: Resolve adsorption parameters with fallback chain.
- ``extract_builder_params``: Extract INCAR parameters from builder_inputs.
- ``add_adsorption_energy_task``: Add the adsorption energy scatter task.
- ``add_adsorption_energy_stage``: Main orchestrator for the adsorption stage.

These functions implement the adsorption energy stage (Stage 10) of the
``build_core_workgraph()`` function.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from aiida import orm
    from aiida_workgraph import WorkGraph

__all__ = [
    "resolve_adsorption_parameters",
    "extract_builder_params",
    "add_adsorption_energy_task",
    "add_adsorption_energy_stage",
]


def resolve_adsorption_parameters(
    adsorption_parameters: Optional[Dict[str, Any]],
    adsorption_options: Optional[Dict[str, Any]],
    adsorption_potential_mapping: Optional[Dict[str, str]],
    adsorption_kpoints_spacing: Optional[float],
    slab_parameters: Optional[Dict[str, Any]],
    slab_options: Optional[Dict[str, Any]],
    slab_potential_mapping: Optional[Dict[str, str]],
    slab_kpoints_spacing: Optional[float],
    bulk_parameters: Dict[str, Any],
    bulk_options: Dict[str, Any],
    bulk_potential_mapping: Dict[str, str],
    kpoints_spacing: float,
) -> tuple[Dict[str, Any], Dict[str, Any], Dict[str, str], float]:
    """Resolve adsorption parameters with fallback chain.

    Implements the fallback chain: adsorption -> slab -> bulk for each
    parameter type.

    Args:
        adsorption_parameters: Adsorption-specific VASP INCAR parameters.
        adsorption_options: Adsorption-specific scheduler options.
        adsorption_potential_mapping: Adsorption-specific potential mapping.
        adsorption_kpoints_spacing: Adsorption-specific k-points spacing.
        slab_parameters: Slab VASP INCAR parameters (intermediate fallback).
        slab_options: Slab scheduler options (intermediate fallback).
        slab_potential_mapping: Slab potential mapping (intermediate fallback).
        slab_kpoints_spacing: Slab k-points spacing (intermediate fallback).
        bulk_parameters: Bulk VASP INCAR parameters (final fallback).
        bulk_options: Bulk scheduler options (final fallback).
        bulk_potential_mapping: Bulk potential mapping (final fallback).
        kpoints_spacing: Main k-points spacing (final fallback).

    Returns:
        Tuple of (resolved_params, resolved_opts, resolved_pot_map, resolved_kpts):
            - resolved_params: VASP INCAR parameters for adsorption
            - resolved_opts: Scheduler options for adsorption
            - resolved_pot_map: Potential mapping for adsorption
            - resolved_kpts: K-points spacing for adsorption

    Example:
        >>> params, opts, pot_map, kpts = resolve_adsorption_parameters(
        ...     adsorption_parameters=None,
        ...     adsorption_options={'queue_name': 'fast'},
        ...     adsorption_potential_mapping=None,
        ...     adsorption_kpoints_spacing=None,
        ...     slab_parameters={'ENCUT': 400},
        ...     slab_options=None,
        ...     slab_potential_mapping={'Ag': 'Ag'},
        ...     slab_kpoints_spacing=0.04,
        ...     bulk_parameters={'ENCUT': 520},
        ...     bulk_options={'queue_name': 'normal'},
        ...     bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     kpoints_spacing=0.03,
        ... )
        >>> params  # From slab (intermediate fallback)
        {'ENCUT': 400}
        >>> opts  # From adsorption (specified)
        {'queue_name': 'fast'}
    """
    # First level fallback: slab -> bulk
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

    # Second level fallback: adsorption -> slab
    ads_params = (
        adsorption_parameters if adsorption_parameters is not None else slab_params
    )
    ads_opts = adsorption_options if adsorption_options is not None else slab_opts
    ads_pot_map = (
        adsorption_potential_mapping
        if adsorption_potential_mapping is not None
        else slab_pot_map
    )
    ads_kpts = (
        adsorption_kpoints_spacing
        if adsorption_kpoints_spacing is not None
        else slab_kpts
    )

    return ads_params, ads_opts, ads_pot_map, ads_kpts


def extract_builder_params(
    relax_builder_inputs: Optional[Dict[str, Any]],
    scf_builder_inputs: Optional[Dict[str, Any]],
    fallback_params: Optional[Dict[str, Any]],
) -> tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """Extract INCAR parameters from builder_inputs for backward compatibility.

    This function extracts INCAR parameters from the new-style builder_inputs
    format while maintaining backward compatibility with the old-style
    direct parameter passing.

    Args:
        relax_builder_inputs: New-style builder inputs for relaxation.
            Expected format: {'parameters': {'incar': {...}}, ...}
        scf_builder_inputs: New-style builder inputs for SCF.
            Expected format: {'parameters': {'incar': {...}}, ...}
        fallback_params: Fallback parameters to use if builder_inputs
            not provided (old-style direct parameters).

    Returns:
        Tuple of (relax_params, scf_params):
            - relax_params: Extracted relaxation parameters or fallback
            - scf_params: Extracted SCF parameters or fallback

    Example:
        >>> # New-style with builder_inputs
        >>> relax_inputs = {'parameters': {'incar': {'NSW': 100, 'IBRION': 2}}}
        >>> scf_inputs = {'parameters': {'incar': {'NSW': 0}}}
        >>> relax_p, scf_p = extract_builder_params(relax_inputs, scf_inputs, None)
        >>> relax_p
        {'NSW': 100, 'IBRION': 2}

        >>> # Old-style with direct parameters
        >>> fallback = {'ENCUT': 520, 'EDIFF': 1e-6}
        >>> relax_p, scf_p = extract_builder_params(None, None, fallback)
        >>> relax_p
        {'ENCUT': 520, 'EDIFF': 1e-6}
    """
    relax_params = None
    scf_params = None

    if relax_builder_inputs is not None:
        # New style: extract INCAR from builder_inputs (for backward compat)
        relax_params = relax_builder_inputs.get("parameters", {})
        if "incar" in relax_params:
            relax_params = relax_params["incar"]
    elif fallback_params is not None:
        # Old style: use ads_params directly
        relax_params = dict(fallback_params)

    if scf_builder_inputs is not None:
        # New style: extract INCAR from builder_inputs (for backward compat)
        scf_params = scf_builder_inputs.get("parameters", {})
        if "incar" in scf_params:
            scf_params = scf_params["incar"]
    elif fallback_params is not None:
        # Old style: use ads_params directly
        scf_params = dict(fallback_params)

    return relax_params, scf_params


def add_adsorption_energy_task(
    wg: WorkGraph,
    code: orm.Code,
    adsorption_structures: Dict[str, Any],
    adsorption_formulas: Dict[str, str],
    potential_family: str,
    potential_mapping: Dict[str, str],
    relax_before_adsorption: bool,
    relax_params: Optional[Dict[str, Any]],
    scf_params: Optional[Dict[str, Any]],
    options: Dict[str, Any],
    kpoints_spacing: float,
    clean_workdir: bool,
    relax_builder_inputs: Optional[Dict[str, Any]],
    scf_builder_inputs: Optional[Dict[str, Any]],
    structure_specific_relax_builder_inputs: Optional[Dict[str, Any]],
    structure_specific_scf_builder_inputs: Optional[Dict[str, Any]],
    structure_component_specific_scf_builder_inputs: Optional[Dict[str, Any]],
    fix_atoms: bool,
    fix_type: Optional[str],
    fix_thickness: float,
    fix_elements: Optional[List[str]],
    max_concurrent_jobs: Optional[int],
) -> None:
    """Add adsorption energy scatter task and connect outputs.

    This function adds the compute_adsorption_energies_scatter task to the
    workgraph and connects all outputs to the workgraph outputs.

    Args:
        wg: WorkGraph instance to add the task to.
        code: Loaded AiiDA Code for VASP calculations.
        adsorption_structures: Dictionary of adsorption structures.
        adsorption_formulas: Dictionary of adsorbate formulas.
        potential_family: POTCAR family name.
        potential_mapping: Element to potential mapping.
        relax_before_adsorption: Whether to relax before SCF.
        relax_params: Relaxation INCAR parameters (old-style).
        scf_params: SCF INCAR parameters (old-style).
        options: Scheduler options.
        kpoints_spacing: K-points spacing in A^-1.
        clean_workdir: Whether to clean work directory.
        relax_builder_inputs: New-style builder inputs for relaxation.
        scf_builder_inputs: New-style builder inputs for SCF.
        structure_specific_relax_builder_inputs: Per-structure relax overrides.
        structure_specific_scf_builder_inputs: Per-structure SCF overrides.
        structure_component_specific_scf_builder_inputs: Per-structure,
            per-component SCF overrides.
        fix_atoms: Whether to apply atom fixing constraints.
        fix_type: Type of fixing ('bottom', 'top', 'center').
        fix_thickness: Thickness for fixing region in Angstrom.
        fix_elements: Elements to fix (None = all).
        max_concurrent_jobs: Maximum concurrent jobs (None = unlimited).

    Returns:
        None. The task is added to wg and outputs are connected.

    Note:
        This function modifies the workgraph in-place by adding the
        adsorption task and connecting its outputs to wg.outputs.
    """
    from aiida import orm

    from teros.core.adsorption_energy import compute_adsorption_energies_scatter

    # Add adsorption energy scatter task (with NEW builder inputs support)
    adsorption_task = wg.add_task(
        compute_adsorption_energies_scatter,
        name="compute_adsorption_energies_scatter",
        structures=adsorption_structures,
        adsorbate_formulas=adsorption_formulas,
        code=code,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        # Relaxation parameters (old-style for backward compat)
        relax_before_adsorption=relax_before_adsorption,
        relax_parameters=relax_params,
        scf_parameters=scf_params,
        # Common settings (old-style for backward compat)
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
        # NEW: Builder inputs (full control)
        relax_builder_inputs=relax_builder_inputs,
        scf_builder_inputs=scf_builder_inputs,
        structure_specific_relax_builder_inputs=structure_specific_relax_builder_inputs,
        structure_specific_scf_builder_inputs=structure_specific_scf_builder_inputs,
        structure_component_specific_scf_builder_inputs=structure_component_specific_scf_builder_inputs,
        # Atom fixing parameters
        fix_atoms=fix_atoms,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        fix_elements=fix_elements,
        # Concurrency control
        max_number_jobs=(
            orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None
        ),
    )

    # Connect outputs
    wg.outputs.relaxed_complete_structures = (
        adsorption_task.outputs.relaxed_complete_structures
    )
    wg.outputs.separated_structures = adsorption_task.outputs.separated_structures
    wg.outputs.substrate_energies = adsorption_task.outputs.substrate_energies
    wg.outputs.molecule_energies = adsorption_task.outputs.molecule_energies
    wg.outputs.complete_energies = adsorption_task.outputs.complete_energies
    wg.outputs.adsorption_energies = adsorption_task.outputs.adsorption_energies


def add_adsorption_energy_stage(
    wg: WorkGraph,
    code_label: str,
    adsorption_structures: Dict[str, Any],
    adsorption_formulas: Dict[str, str],
    potential_family: str,
    # Adsorption-specific parameters (may be None for fallback)
    adsorption_potential_mapping: Optional[Dict[str, str]] = None,
    adsorption_parameters: Optional[Dict[str, Any]] = None,
    adsorption_options: Optional[Dict[str, Any]] = None,
    adsorption_kpoints_spacing: Optional[float] = None,
    # Relaxation control
    relax_before_adsorption: bool = False,
    # Builder inputs (new-style API)
    adsorption_relax_builder_inputs: Optional[Dict[str, Any]] = None,
    adsorption_scf_builder_inputs: Optional[Dict[str, Any]] = None,
    adsorption_structure_specific_relax_builder_inputs: Optional[Dict[str, Any]] = None,
    adsorption_structure_specific_scf_builder_inputs: Optional[Dict[str, Any]] = None,
    adsorption_structure_component_specific_scf_builder_inputs: Optional[
        Dict[str, Any]
    ] = None,
    # Atom fixing
    adsorption_fix_atoms: bool = False,
    adsorption_fix_type: Optional[str] = None,
    adsorption_fix_thickness: float = 0.0,
    adsorption_fix_elements: Optional[List[str]] = None,
    # Common settings
    clean_workdir: bool = False,
    max_concurrent_jobs: Optional[int] = None,
    # Fallback parameters (slab -> bulk chain)
    slab_parameters: Optional[Dict[str, Any]] = None,
    bulk_parameters: Optional[Dict[str, Any]] = None,
    slab_options: Optional[Dict[str, Any]] = None,
    bulk_options: Optional[Dict[str, Any]] = None,
    slab_potential_mapping: Optional[Dict[str, str]] = None,
    bulk_potential_mapping: Optional[Dict[str, str]] = None,
    slab_kpoints_spacing: Optional[float] = None,
    kpoints_spacing: float = 0.03,
    # Logging
    logger: Optional[logging.Logger] = None,
) -> None:
    """Add adsorption energy calculation stage to the workgraph.

    This is the main orchestrator function for Stage 10 (Adsorption Energy)
    of build_core_workgraph(). It:

    1. Loads the VASP code
    2. Resolves parameters with fallback chain (adsorption -> slab -> bulk)
    3. Extracts relax and SCF parameters from builder_inputs
    4. Adds the compute_adsorption_energies_scatter task
    5. Connects all outputs to wg.outputs

    Args:
        wg: WorkGraph instance to add the adsorption stage to.
        code_label: Label of the AiiDA Code for VASP calculations.
        adsorption_structures: Dictionary mapping labels to adsorption
            StructureData nodes. Keys should match adsorption_formulas.
        adsorption_formulas: Dictionary mapping labels to adsorbate
            chemical formulas (e.g., {'term_0': 'OH'}).
        potential_family: POTCAR family name (e.g., 'PBE', 'PBE.54').

        adsorption_potential_mapping: Element to potential mapping for
            adsorption calculations. Falls back to slab -> bulk if None.
        adsorption_parameters: VASP INCAR parameters for adsorption.
            Falls back to slab -> bulk if None.
        adsorption_options: Scheduler options for adsorption calculations.
            Falls back to slab -> bulk if None.
        adsorption_kpoints_spacing: K-points spacing for adsorption.
            Falls back to slab -> bulk if None.

        relax_before_adsorption: If True, relax complete structures
            before separation and SCF calculations.

        adsorption_relax_builder_inputs: New-style builder inputs for
            relaxation. Format: {'parameters': {'incar': {...}}, ...}.
        adsorption_scf_builder_inputs: New-style builder inputs for SCF.
            Format: {'parameters': {'incar': {...}}, ...}.
        adsorption_structure_specific_relax_builder_inputs: Per-structure
            overrides for relaxation. Format: {0: {...}, 1: {...}}.
        adsorption_structure_specific_scf_builder_inputs: Per-structure
            overrides for SCF calculations.
        adsorption_structure_component_specific_scf_builder_inputs:
            Per-structure, per-component SCF overrides.
            Format: {0: {'substrate': {...}, 'molecule': {...}}}.

        adsorption_fix_atoms: If True, apply atom fixing constraints
            during relaxation.
        adsorption_fix_type: Type of fixing ('bottom', 'top', 'center').
        adsorption_fix_thickness: Thickness in Angstrom for fixing region.
        adsorption_fix_elements: Elements to fix (None = all elements).

        clean_workdir: Whether to clean work directory after calculations.
        max_concurrent_jobs: Maximum concurrent VASP jobs (None = unlimited).

        slab_parameters: Slab VASP parameters (intermediate fallback).
        bulk_parameters: Bulk VASP parameters (final fallback).
        slab_options: Slab scheduler options (intermediate fallback).
        bulk_options: Bulk scheduler options (final fallback).
        slab_potential_mapping: Slab potential mapping (intermediate fallback).
        bulk_potential_mapping: Bulk potential mapping (final fallback).
        slab_kpoints_spacing: Slab k-points spacing (intermediate fallback).
        kpoints_spacing: Main k-points spacing (final fallback).

        logger: Logger instance for status messages. If None, uses module logger.

    Returns:
        None. The workgraph is modified in-place.

    Example:
        >>> from aiida_workgraph import WorkGraph
        >>> wg = WorkGraph('my_workflow')
        >>> add_adsorption_energy_stage(
        ...     wg=wg,
        ...     code_label='VASP-6.5.1@cluster',
        ...     adsorption_structures={'sys1': struct1, 'sys2': struct2},
        ...     adsorption_formulas={'sys1': 'OH', 'sys2': 'OOH'},
        ...     potential_family='PBE',
        ...     bulk_potential_mapping={'La': 'La', 'Mn': 'Mn', 'O': 'O', 'H': 'H'},
        ...     bulk_parameters={'ENCUT': 520, 'EDIFF': 1e-6},
        ...     bulk_options={'resources': {'num_machines': 1}},
        ...     kpoints_spacing=0.03,
        ...     relax_before_adsorption=True,
        ... )
    """
    from aiida.orm import load_code

    # Use module logger if none provided
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("  -> Adding adsorption energy calculation")
    logger.info("     Number of structures: %d", len(adsorption_structures))
    logger.info("     Structure keys: %s", list(adsorption_structures.keys()))
    logger.info("     Adsorbate formulas: %s", adsorption_formulas)
    if relax_before_adsorption:
        logger.info("     Relaxation: ENABLED (relax -> separate -> SCF)")
    else:
        logger.info("     Relaxation: DISABLED (separate -> SCF)")

    # Load code
    code = load_code(code_label)

    # Resolve parameters with fallback chain (adsorption -> slab -> bulk)
    # Provide empty dicts as final fallbacks if bulk params are None
    _bulk_params = bulk_parameters if bulk_parameters is not None else {}
    _bulk_opts = bulk_options if bulk_options is not None else {}
    _bulk_pot_map = bulk_potential_mapping if bulk_potential_mapping is not None else {}

    ads_params, ads_opts, ads_pot_map, ads_kpts = resolve_adsorption_parameters(
        adsorption_parameters=adsorption_parameters,
        adsorption_options=adsorption_options,
        adsorption_potential_mapping=adsorption_potential_mapping,
        adsorption_kpoints_spacing=adsorption_kpoints_spacing,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        bulk_parameters=_bulk_params,
        bulk_options=_bulk_opts,
        bulk_potential_mapping=_bulk_pot_map,
        kpoints_spacing=kpoints_spacing,
    )

    # Extract relax and SCF parameters from builder_inputs (for backward compat)
    relax_params, scf_params = extract_builder_params(
        relax_builder_inputs=adsorption_relax_builder_inputs,
        scf_builder_inputs=adsorption_scf_builder_inputs,
        fallback_params=ads_params,
    )

    # Add adsorption energy task with all parameters
    add_adsorption_energy_task(
        wg=wg,
        code=code,
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        potential_family=potential_family,
        potential_mapping=ads_pot_map,
        relax_before_adsorption=relax_before_adsorption,
        relax_params=relax_params,
        scf_params=scf_params,
        options=ads_opts,
        kpoints_spacing=ads_kpts,
        clean_workdir=clean_workdir,
        relax_builder_inputs=adsorption_relax_builder_inputs,
        scf_builder_inputs=adsorption_scf_builder_inputs,
        structure_specific_relax_builder_inputs=adsorption_structure_specific_relax_builder_inputs,
        structure_specific_scf_builder_inputs=adsorption_structure_specific_scf_builder_inputs,
        structure_component_specific_scf_builder_inputs=adsorption_structure_component_specific_scf_builder_inputs,
        fix_atoms=adsorption_fix_atoms,
        fix_type=adsorption_fix_type,
        fix_thickness=adsorption_fix_thickness,
        fix_elements=adsorption_fix_elements,
        max_concurrent_jobs=max_concurrent_jobs,
    )

    logger.info("  + Adsorption energy calculation enabled")
    if relax_before_adsorption:
        logger.info(
            "     Access relaxed structures via: wg.outputs.relaxed_complete_structures"
        )
    logger.info("     Access adsorption energies via: wg.outputs.adsorption_energies")
