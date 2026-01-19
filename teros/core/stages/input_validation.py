"""Stage 2: Input validation and configuration logging.

This module provides functions for validating workflow inputs based on the
resolved workflow configuration, and for logging the configuration in a
formatted manner.

The input validation stage determines whether bulk structure files are needed
based on the enabled workflow features, and raises clear error messages when
required inputs are missing.

Functions:
    validate_workflow_inputs: Validate required inputs and determine bulk needs.
    log_workflow_configuration: Log the workflow configuration.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from aiida import orm

__all__ = [
    "validate_workflow_inputs",
    "log_workflow_configuration",
]

logger = logging.getLogger(__name__)


def validate_workflow_inputs(
    resolved_preset_name: str,
    resolved_flags: Dict[str, bool],
    structures_dir: Optional[str],
    bulk_name: Optional[str],
    metal_name: Optional[str],
    oxygen_name: Optional[str],
    miller_indices: Optional[List[tuple]],
    input_slabs: Optional[Dict[str, Any]],
) -> bool:
    """Validate that required inputs are provided based on workflow configuration.

    This function performs the following:

    1. Determines if bulk structure is needed based on workflow flags:
       - Formation enthalpy calculation requires bulk (metal_name and oxygen_name provided)
       - Surface thermodynamics requires bulk
       - Bulk electronic properties calculation requires bulk
       - Slab generation via miller_indices requires bulk (when input_slabs not provided)

    2. Validates that structures_dir and bulk_name are provided when needed.

    3. Logs the workflow configuration via ``log_workflow_configuration()``.

    4. Returns whether bulk structure is needed for downstream processing.

    Args:
        resolved_preset_name: Name of the resolved workflow preset (e.g.,
            'surface_thermodynamics', 'aimd_only').
        resolved_flags: Dictionary of resolved workflow flags. Expected keys:
            - 'compute_thermodynamics': Whether surface thermodynamics is enabled.
            - 'compute_electronic_properties_bulk': Whether bulk electronic
              properties (DOS/bands) calculation is enabled.
        structures_dir: Path to directory containing structure files. Required
            when bulk structure is needed.
        bulk_name: Filename of the bulk structure within structures_dir.
            Required when bulk structure is needed.
        metal_name: Filename of the metal reference structure. If provided along
            with oxygen_name, formation enthalpy will be calculated, requiring
            bulk structure.
        oxygen_name: Filename of the oxygen reference structure. If provided
            along with metal_name, formation enthalpy will be calculated,
            requiring bulk structure.
        miller_indices: List of Miller indices for slab generation. If provided
            without input_slabs, bulk structure is needed to generate slabs.
        input_slabs: Dictionary of pre-generated slab structures. If provided,
            bulk structure is not needed for slab generation.

    Returns:
        True if bulk structure is needed for the workflow, False otherwise.
        This value can be used by downstream stages to conditionally skip
        bulk-related processing.

    Raises:
        ValueError: If bulk structure is needed but structures_dir or bulk_name
            are not provided. The error message includes details about which
            features require the bulk structure.

    Example:
        >>> # Workflow that needs bulk (surface thermodynamics)
        >>> needs_bulk = validate_workflow_inputs(
        ...     resolved_preset_name='surface_thermodynamics',
        ...     resolved_flags={'compute_thermodynamics': True,
        ...                     'compute_electronic_properties_bulk': False},
        ...     structures_dir='/path/to/structures',
        ...     bulk_name='Ag2O_bulk.vasp',
        ...     metal_name=None,
        ...     oxygen_name=None,
        ...     miller_indices=[(1, 1, 0)],
        ...     input_slabs=None,
        ... )
        >>> needs_bulk
        True

        >>> # Workflow that doesn't need bulk (aimd_only with input slabs)
        >>> needs_bulk = validate_workflow_inputs(
        ...     resolved_preset_name='aimd_only',
        ...     resolved_flags={'compute_thermodynamics': False,
        ...                     'compute_electronic_properties_bulk': False},
        ...     structures_dir=None,
        ...     bulk_name=None,
        ...     metal_name=None,
        ...     oxygen_name=None,
        ...     miller_indices=None,
        ...     input_slabs={'slab_0': some_structure},
        ... )
        >>> needs_bulk
        False
    """
    # Determine if formation enthalpy calculation is enabled
    compute_formation_enthalpy = metal_name is not None and oxygen_name is not None

    # Extract flags (with safe defaults)
    compute_thermodynamics = resolved_flags.get("compute_thermodynamics", False)
    compute_electronic_properties_bulk = resolved_flags.get(
        "compute_electronic_properties_bulk", False
    )

    # Determine if bulk structure is needed
    needs_bulk = (
        compute_formation_enthalpy
        or compute_thermodynamics
        or compute_electronic_properties_bulk
        or (miller_indices is not None and input_slabs is None)
    )

    # Only require structures_dir and bulk_name if we actually need the bulk structure
    if needs_bulk:
        if structures_dir is None or bulk_name is None:
            # Build detailed error message showing which features require bulk
            reasons = []
            if compute_formation_enthalpy:
                reasons.append("  - Formation enthalpy calculation")
            if compute_thermodynamics:
                reasons.append("  - Surface thermodynamics")
            if compute_electronic_properties_bulk:
                reasons.append("  - Electronic properties (bulk)")
            if miller_indices is not None and input_slabs is None:
                reasons.append("  - Slab generation (miller_indices provided)")

            raise ValueError(
                "structures_dir and bulk_name are required for workflows that "
                "need bulk structure.\n"
                f"Current preset '{resolved_preset_name}' requires bulk structure for:\n"
                + "\n".join(reasons)
                + "\n\nEither provide structures_dir and bulk_name, or use "
                "input_slabs with aimd_only preset."
            )

    # Log workflow configuration
    log_workflow_configuration(resolved_preset_name, resolved_flags)

    return needs_bulk


def log_workflow_configuration(
    preset_name: str,
    resolved_flags: Dict[str, bool],
) -> None:
    """Log the workflow configuration in a formatted manner.

    Outputs a summary of the workflow preset and which features are enabled
    or disabled. Uses checkmarks and crosses for visual clarity.

    Args:
        preset_name: Name of the resolved workflow preset.
        resolved_flags: Dictionary mapping flag names to their boolean values.
            Common keys include:
            - 'relax_slabs'
            - 'compute_thermodynamics'
            - 'compute_cleavage'
            - 'compute_relaxation_energy'
            - 'compute_electronic_properties_bulk'
            - 'compute_electronic_properties_slabs'
            - 'run_aimd'
            - 'run_adsorption_energy'

    Example:
        >>> log_workflow_configuration(
        ...     'surface_thermodynamics',
        ...     {'relax_slabs': True, 'compute_thermodynamics': True,
        ...      'compute_cleavage': False, 'run_aimd': False}
        ... )
        # Output (via logger.info):
        # ======================================================================
        # WORKFLOW CONFIGURATION
        # ======================================================================
        # Preset: surface_thermodynamics
        # Active Components:
        #   [check] relax_slabs: True
        #   [check] compute_thermodynamics: True
        #   [cross] compute_cleavage: False
        #   [cross] run_aimd: False
        # ======================================================================
    """
    separator = "=" * 70

    logger.info(separator)
    logger.info("WORKFLOW CONFIGURATION")
    logger.info(separator)
    logger.info("Preset: %s", preset_name)
    logger.info("Active Components:")

    for flag_name, flag_value in resolved_flags.items():
        # Use Unicode checkmark/cross for visual indication
        status = "\u2713" if flag_value else "\u2717"
        logger.info("  %s %s: %s", status, flag_name, flag_value)

    logger.info(separator)
