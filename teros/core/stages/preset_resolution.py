"""
Stage 1: Workflow Preset Resolution.

This module handles the resolution of workflow presets and validation of inputs
for the PS-TEROS build_core_workgraph() function.

The preset resolution stage performs:
1. Check for deprecated old-style API usage (flags without preset)
2. Resolve preset and apply user overrides
3. Validate preset requirements (required parameters)
4. Check flag dependencies and emit warnings/errors

Example:
    >>> from teros.core.stages.preset_resolution import resolve_workflow_preset
    >>> preset_name, flags = resolve_workflow_preset(
    ...     workflow_preset='surface_thermodynamics',
    ...     relax_slabs=True,
    ...     compute_thermodynamics=True,
    ...     compute_cleavage=None,  # Use preset default
    ...     compute_relaxation_energy=None,
    ...     compute_electronic_properties_bulk=None,
    ...     compute_electronic_properties_slabs=None,
    ...     run_aimd=None,
    ...     run_adsorption_energy=None,
    ...     metal_name='Ag.cif',
    ...     oxygen_name='O2.cif',
    ... )
    >>> print(preset_name)
    'surface_thermodynamics'
    >>> print(flags['relax_slabs'])
    True
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

from teros.core.workflow_presets import (
    check_old_style_api,
    resolve_preset,
    validate_flag_dependencies,
    validate_preset_inputs,
)

__all__ = ["resolve_workflow_preset"]


def resolve_workflow_preset(
    workflow_preset: Optional[str],
    relax_slabs: Optional[bool],
    compute_thermodynamics: Optional[bool],
    compute_cleavage: Optional[bool],
    compute_relaxation_energy: Optional[bool],
    compute_electronic_properties_bulk: Optional[bool],
    compute_electronic_properties_slabs: Optional[bool],
    run_aimd: Optional[bool],
    run_adsorption_energy: Optional[bool],
    # Validation inputs
    metal_name: Optional[str] = None,
    oxygen_name: Optional[str] = None,
    nonmetal_name: Optional[str] = None,
    miller_indices: Optional[list] = None,
    input_slabs: Optional[dict] = None,
    bands_parameters: Optional[dict] = None,
    bands_options: Optional[dict] = None,
    band_settings: Optional[dict] = None,
    aimd_sequence: Optional[list] = None,
    aimd_parameters: Optional[dict] = None,
    aimd_options: Optional[dict] = None,
    aimd_potential_mapping: Optional[dict] = None,
    aimd_kpoints_spacing: Optional[float] = None,
    slab_bands_parameters: Optional[dict] = None,
    slab_bands_options: Optional[dict] = None,
    slab_band_settings: Optional[dict] = None,
    slab_electronic_properties: Optional[dict] = None,
    adsorption_structures: Optional[dict] = None,
    adsorption_formulas: Optional[dict] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[str, dict]:
    """Resolve workflow preset and apply user overrides.

    This function is Stage 1 of the build_core_workgraph() decomposition.
    It performs the complete preset resolution process:

    1. **Check old-style API**: Emit deprecation warning if user specifies flags
       without an explicit workflow_preset parameter.

    2. **Resolve preset**: Load the preset configuration and apply any user
       overrides for individual flags.

    3. **Validate preset inputs**: Check that all required parameters for the
       selected preset are provided.

    4. **Validate flag dependencies**: Check that flag dependencies are satisfied
       (e.g., compute_cleavage requires relax_slabs).

    Args:
        workflow_preset: Name of workflow preset to use. If None, uses the default
            preset ('surface_thermodynamics'). Available presets:
            - 'surface_thermodynamics': Complete surface thermodynamics workflow
            - 'surface_thermodynamics_unrelaxed': Surface thermodynamics without relaxation
            - 'bulk_only': Bulk structure optimization only
            - 'formation_enthalpy_only': Formation enthalpy calculation
            - 'cleavage_only': Cleavage energy calculations
            - 'relaxation_energy_only': Relaxation energy calculations
            - 'electronic_structure_bulk_only': DOS/bands for bulk
            - 'electronic_structure_slabs_only': DOS/bands for slabs
            - 'electronic_structure_bulk_and_slabs': DOS/bands for both
            - 'aimd_only': AIMD simulation
            - 'adsorption_energy': Adsorption energy calculations
            - 'comprehensive': Everything enabled

        relax_slabs: Override preset default for slab relaxation. If None,
            uses the preset default.

        compute_thermodynamics: Override preset default for thermodynamics
            calculation. If None, uses the preset default.

        compute_cleavage: Override preset default for cleavage energy
            calculation. If None, uses the preset default.

        compute_relaxation_energy: Override preset default for relaxation
            energy calculation. If None, uses the preset default.

        compute_electronic_properties_bulk: Override preset default for bulk
            electronic properties (DOS/bands). If None, uses the preset default.

        compute_electronic_properties_slabs: Override preset default for slab
            electronic properties (DOS/bands). If None, uses the preset default.

        run_aimd: Override preset default for AIMD simulation. If None,
            uses the preset default.

        run_adsorption_energy: Override preset default for adsorption energy
            calculation. If None, uses the preset default.

        metal_name: Filename of metal reference structure. Required for
            thermodynamics calculations.

        oxygen_name: Filename of oxygen reference structure (e.g., 'O2.cif').
            Required for thermodynamics calculations.

        nonmetal_name: Filename of nonmetal reference structure (optional).
            Required for ternary/quaternary oxides.

        miller_indices: List of Miller indices for slab generation.
            E.g., [[1, 0, 0], [1, 1, 0], [1, 1, 1]].

        input_slabs: Dictionary of pre-generated slab structures. If provided,
            slabs are not generated from miller_indices.

        bands_parameters: VASP INCAR parameters for bulk electronic structure
            calculation. Required for compute_electronic_properties_bulk.

        bands_options: Scheduler options for bulk electronic structure
            calculation.

        band_settings: Band structure calculation settings (e.g., k-path).

        aimd_sequence: List of AIMD stage configurations. Required for run_aimd.
            Each element is a dict with 'temperature', 'timesteps', etc.

        aimd_parameters: VASP INCAR parameters for AIMD. Required for run_aimd.

        aimd_options: Scheduler options for AIMD calculations.

        aimd_potential_mapping: Element to potential mapping for AIMD.

        aimd_kpoints_spacing: K-points spacing for AIMD calculations.

        slab_bands_parameters: VASP INCAR parameters for slab electronic
            structure calculation. Required for compute_electronic_properties_slabs.

        slab_bands_options: Scheduler options for slab electronic structure.

        slab_band_settings: Slab band structure calculation settings.

        slab_electronic_properties: Dictionary specifying which slabs to
            calculate electronic properties for, with optional overrides.

        adsorption_structures: Dictionary mapping labels to adsorption structure
            file paths. Required for run_adsorption_energy.

        adsorption_formulas: Dictionary mapping labels to adsorbate chemical
            formulas (e.g., {'term_0': 'OH'}). Required for run_adsorption_energy.

        logger: Logger instance for emitting warnings. If None, uses module logger.

    Returns:
        Tuple of (resolved_preset_name, resolved_flags_dict):
            - resolved_preset_name: Name of the resolved preset (string)
            - resolved_flags_dict: Dictionary with all resolved flags:
                - 'relax_slabs': bool
                - 'compute_thermodynamics': bool
                - 'compute_cleavage': bool
                - 'compute_relaxation_energy': bool
                - 'compute_electronic_properties_bulk': bool
                - 'compute_electronic_properties_slabs': bool
                - 'run_aimd': bool
                - 'run_adsorption_energy': bool

    Raises:
        ValueError: If:
            - workflow_preset is not a valid preset name
            - Required parameters for the preset are missing
            - Flag dependencies are not satisfied (errors, not warnings)

    Example:
        Basic usage with default preset:

        >>> preset_name, flags = resolve_workflow_preset(
        ...     workflow_preset=None,  # Uses default
        ...     relax_slabs=None,
        ...     compute_thermodynamics=None,
        ...     compute_cleavage=None,
        ...     compute_relaxation_energy=None,
        ...     compute_electronic_properties_bulk=None,
        ...     compute_electronic_properties_slabs=None,
        ...     run_aimd=None,
        ...     run_adsorption_energy=None,
        ...     metal_name='Ag.cif',
        ...     oxygen_name='O2.cif',
        ... )
        >>> preset_name
        'surface_thermodynamics'

        With explicit preset and override:

        >>> preset_name, flags = resolve_workflow_preset(
        ...     workflow_preset='surface_thermodynamics',
        ...     relax_slabs=True,
        ...     compute_thermodynamics=True,
        ...     compute_cleavage=True,  # Override default (False)
        ...     compute_relaxation_energy=None,
        ...     compute_electronic_properties_bulk=None,
        ...     compute_electronic_properties_slabs=None,
        ...     run_aimd=None,
        ...     run_adsorption_energy=None,
        ...     metal_name='Ag.cif',
        ...     oxygen_name='O2.cif',
        ... )
        >>> flags['compute_cleavage']
        True

    Note:
        This function does not modify any global state. It returns the resolved
        configuration which can then be used by subsequent stages.
    """
    # Use module logger if none provided
    if logger is None:
        logger = logging.getLogger(__name__)

    # ========================================================================
    # Step 1: Check for deprecated old-style API usage
    # ========================================================================
    check_old_style_api(
        workflow_preset,
        relax_slabs,
        compute_thermodynamics,
        compute_cleavage,
        compute_relaxation_energy,
        compute_electronic_properties_bulk,
        compute_electronic_properties_slabs,
        run_aimd,
        run_adsorption_energy,
    )

    # ========================================================================
    # Step 2: Resolve preset and apply user overrides
    # ========================================================================
    resolved_preset_name, resolved_flags = resolve_preset(
        workflow_preset,
        relax_slabs,
        compute_thermodynamics,
        compute_cleavage,
        compute_relaxation_energy,
        compute_electronic_properties_bulk,
        compute_electronic_properties_slabs,
        run_aimd,
        run_adsorption_energy,
    )

    # ========================================================================
    # Step 3: Validate preset requirements
    # ========================================================================
    validation_errors = validate_preset_inputs(
        resolved_preset_name,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        nonmetal_name=nonmetal_name,
        miller_indices=miller_indices,
        bands_parameters=bands_parameters,
        bands_options=bands_options,
        band_settings=band_settings,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        aimd_options=aimd_options,
        aimd_potential_mapping=aimd_potential_mapping,
        aimd_kpoints_spacing=aimd_kpoints_spacing,
        slab_bands_parameters=slab_bands_parameters,
        slab_bands_options=slab_bands_options,
        slab_band_settings=slab_band_settings,
        slab_electronic_properties=slab_electronic_properties,
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
    )

    if validation_errors:
        error_msg = "\n".join(validation_errors)
        raise ValueError(f"Preset validation failed:\n{error_msg}")

    # ========================================================================
    # Step 4: Check flag dependencies and emit warnings
    # ========================================================================
    dependency_warnings = validate_flag_dependencies(
        resolved_flags,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        miller_indices=miller_indices,
        input_slabs=input_slabs,
        bands_parameters=bands_parameters,
        slab_bands_parameters=slab_bands_parameters,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
    )

    for warning in dependency_warnings:
        if warning.startswith("ERROR:"):
            raise ValueError(warning)
        logger.warning(warning)

    return resolved_preset_name, resolved_flags
