"""
PS-TEROS Workflow Stages Module.

This module provides composable stage functions for building the core workgraph.
Each stage is independently testable and maintains clear responsibilities.

Stage Builder Pattern
---------------------

The ``build_core_workgraph()`` function orchestrates 10 distinct stages of workflow
construction. Each stage is implemented as a separate module with focused responsibilities:

1. **Preset Resolution** (``preset_resolution.py``): Resolve workflow presets and apply
   user overrides to determine which workflow features are enabled.

2. **Input Validation** (``input_validation.py``): Validate required inputs based on
   the resolved workflow configuration and log the configuration.

3. **Restart Handling** (``restart_handling.py``): Extract data from a previous
   workgraph for restart scenarios (folders and slab structures).

4. **CP2K Setup** (``cp2k_setup.py``): Create CP2K basis and pseudopotential files
   when using the CP2K calculator.

5. **Core Build** (``core_build.py``): Build the core workgraph by calling the
   @task.graph decorated function.

6. **Input Slabs** (``input_slabs.py``): Handle manual/restart slab relaxation,
   thermodynamics, and cleavage calculations.

7. **Slab Electronic** (``slab_electronic.py``): Add electronic properties (DOS/bands)
   calculations for selected slabs.

8. **Bulk Electronic** (``bulk_electronic.py``): Add electronic properties (DOS/bands)
   calculations for the relaxed bulk structure.

9. **AIMD Stage** (``aimd_stage.py``): Add multi-stage AIMD simulation tasks.

10. **Adsorption Stage** (``adsorption_stage.py``): Add adsorption energy calculation
    tasks.

StageContext
------------

The ``StageContext`` dataclass holds shared state between stages, providing a clean
interface for passing configuration and intermediate results through the workflow
construction pipeline. This avoids passing dozens of individual parameters between
stage functions.

Usage Example
-------------

.. code-block:: python

    from teros.core.stages import StageContext, resolve_workflow_preset

    # Create context with resolved configuration
    ctx = StageContext(
        resolved_flags={'relax_slabs': True, ...},
        resolved_preset_name='surface_thermodynamics',
        code_label='VASP-6.5.1@cluster',
        potential_family='PBE',
        ...
    )

    # Each stage function receives and may modify the context
    preset_name, flags = resolve_workflow_preset(workflow_preset, ...)
    ctx.resolved_preset_name = preset_name
    ctx.resolved_flags = flags

Public API
----------

Stage Functions:
    - ``resolve_workflow_preset``: Stage 1 - Preset resolution
    - ``validate_workflow_inputs``: Stage 2 - Input validation and bulk needs detection
    - ``log_workflow_configuration``: Stage 2 - Configuration logging
    - ``extract_restart_data``: Stage 3 - Restart handling
    - ``prepare_input_slabs``: Stage 3 - Input slabs preparation
    - ``setup_cp2k_files``: Stage 4 - CP2K setup
    - ``check_cleavage_compatibility``: Stage 4 - Cleavage compatibility check
    - ``build_core_workgraph_base``: Stage 5 - Core workgraph construction
    - ``resolve_slab_parameters``: Stage 6 - Slab parameter resolution
    - ``add_restart_slab_tasks``: Stage 6 - Restart slab tasks
    - ``add_normal_slab_tasks``: Stage 6 - Normal slab tasks
    - ``add_thermodynamics_task``: Stage 6 - Thermodynamics task
    - ``add_cleavage_task``: Stage 6 - Cleavage task
    - ``add_slab_electronic_properties``: Stage 7 - Slab electronic properties
    - ``add_bulk_electronic_properties``: Stage 8 - Bulk electronic properties
    - ``resolve_aimd_parameters``: Stage 9 - AIMD parameter resolution
    - ``prepare_fixed_atoms``: Stage 9 - Fixed atoms preparation
    - ``determine_initial_slabs_source``: Stage 9 - Initial slabs source
    - ``add_aimd_stages``: Stage 9 - AIMD stages
    - ``resolve_adsorption_parameters``: Stage 10 - Adsorption parameter resolution
    - ``extract_builder_params``: Stage 10 - Builder parameter extraction
    - ``add_adsorption_energy_task``: Stage 10 - Adsorption energy task

Type Definitions:
    - ``StageContext``: Shared context between stages
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from aiida import orm
    from aiida_workgraph import WorkGraph

__all__ = [
    # Types
    "StageContext",
    "RestartData",
    # Stage 1: Preset resolution
    "resolve_workflow_preset",
    # Stage 2: Input validation
    "validate_workflow_inputs",
    "log_workflow_configuration",
    # Stage 3: Restart handling
    "extract_restart_data",
    "handle_restart_from_node",
    "prepare_input_slabs",
    # Stage 4: CP2K setup
    "setup_cp2k_files",
    "check_cleavage_compatibility",
    # Stage 5: Core build
    "build_core_workgraph_base",
    # Stage 6: Input slabs handling
    "resolve_slab_parameters",
    "add_restart_slab_tasks",
    "add_normal_slab_tasks",
    "add_thermodynamics_task",
    "add_cleavage_task",
    # Stage 7: Slab electronic properties
    "add_slab_electronic_properties",
    # Stage 8: Bulk electronic properties
    "add_bulk_electronic_properties",
    # Stage 9: AIMD
    "resolve_aimd_parameters",
    "prepare_fixed_atoms",
    "determine_initial_slabs_source",
    "add_aimd_stages",
    # Stage 10: Adsorption energy
    "resolve_adsorption_parameters",
    "extract_builder_params",
    "add_adsorption_energy_task",
]


@dataclass
class StageContext:
    """Shared context passed between workflow stages.

    This dataclass holds the shared state that is built up during workflow
    construction. Each stage may read from and write to this context, avoiding
    the need to pass dozens of individual parameters between functions.

    The context is organized into logical groups:
    - Workflow configuration (resolved flags, preset name)
    - Code configuration (code labels for different calculation types)
    - Potential configuration (family and element mappings)
    - Parameter configuration (INCAR parameters for different calculations)
    - Options configuration (scheduler options for different calculations)
    - K-points configuration
    - Restart/input data
    - Runtime state (logger, workgraph instance)

    Attributes:
        wg: WorkGraph instance being built. May be None before Stage 5.
        resolved_flags: Dictionary of resolved workflow flags (relax_slabs,
            compute_thermodynamics, etc.).
        resolved_preset_name: Name of the resolved workflow preset.

        code_label: Main VASP/CP2K code label (used as fallback).
        bulk_code_label: Code label for bulk calculations (or None to use code_label).
        slab_code_label: Code label for slab calculations (or None to use code_label).
        metal_code_label: Code label for metal reference (or None to use code_label).
        oxygen_code_label: Code label for oxygen reference (or None to use code_label).
        nonmetal_code_label: Code label for nonmetal reference (or None to use code_label).
        aimd_code_label: Code label for AIMD calculations (or None to use code_label).

        potential_family: POTCAR family name (e.g., 'PBE', 'PBE.54').
        bulk_potential_mapping: Element to potential mapping for bulk calculations.
        slab_potential_mapping: Element to potential mapping for slabs (or None to use bulk).
        metal_potential_mapping: Element to potential mapping for metal reference.
        oxygen_potential_mapping: Element to potential mapping for oxygen reference.
        nonmetal_potential_mapping: Element to potential mapping for nonmetal reference.
        aimd_potential_mapping: Element to potential mapping for AIMD (or None for fallback).
        adsorption_potential_mapping: Element to potential mapping for adsorption.

        bulk_parameters: VASP INCAR parameters for bulk calculations.
        slab_parameters: VASP INCAR parameters for slabs (or None to use bulk).
        metal_parameters: VASP INCAR parameters for metal reference.
        oxygen_parameters: VASP INCAR parameters for oxygen reference.
        nonmetal_parameters: VASP INCAR parameters for nonmetal reference.
        aimd_parameters: VASP INCAR parameters for AIMD (or None for fallback).
        adsorption_parameters: VASP INCAR parameters for adsorption.
        bands_parameters: VASP INCAR parameters for bulk DOS/bands.
        slab_bands_parameters: VASP INCAR parameters for slab DOS/bands.

        bulk_options: Scheduler options for bulk calculations.
        slab_options: Scheduler options for slabs (or None to use bulk).
        metal_options: Scheduler options for metal reference.
        oxygen_options: Scheduler options for oxygen reference.
        nonmetal_options: Scheduler options for nonmetal reference.
        aimd_options: Scheduler options for AIMD (or None for fallback).
        adsorption_options: Scheduler options for adsorption.
        bands_options: Scheduler options for bulk DOS/bands.
        slab_bands_options: Scheduler options for slab DOS/bands.

        kpoints_spacing: K-points spacing in A^-1 * 2pi (main setting).
        slab_kpoints_spacing: K-points spacing for slabs (or None to use main).
        aimd_kpoints_spacing: K-points spacing for AIMD (or None for fallback).
        adsorption_kpoints_spacing: K-points spacing for adsorption.

        clean_workdir: Whether to clean work directory after calculations.

        logger: Logger instance for workflow construction messages.

        restart_folders: Dict of RemoteData for restart (or None).
        restart_slabs: Dict of StructureData from restart (or None).
        input_slabs: Dict of user-provided input slab structures (or None).
        use_input_slabs: Whether to use input_slabs instead of generating slabs.

    Example:
        >>> ctx = StageContext(
        ...     resolved_flags={'relax_slabs': True, 'compute_thermodynamics': True},
        ...     resolved_preset_name='surface_thermodynamics',
        ...     code_label='VASP-6.5.1@cluster',
        ...     potential_family='PBE',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     bulk_parameters={'incar': {'ENCUT': 520}},
        ...     bulk_options={'resources': {'num_machines': 1}},
        ...     kpoints_spacing=0.03,
        ...     clean_workdir=False,
        ... )
        >>> ctx.resolved_flags['relax_slabs']
        True
    """

    # WorkGraph instance (may be None before Stage 5)
    wg: Optional[WorkGraph] = None

    # Workflow configuration
    resolved_flags: Dict[str, bool] = field(default_factory=dict)
    resolved_preset_name: str = "surface_thermodynamics"

    # Code configuration
    code_label: str = "VASP-6.5.1@cluster"
    bulk_code_label: Optional[str] = None
    slab_code_label: Optional[str] = None
    metal_code_label: Optional[str] = None
    oxygen_code_label: Optional[str] = None
    nonmetal_code_label: Optional[str] = None
    aimd_code_label: Optional[str] = None

    # Potential configuration
    potential_family: str = "PBE"
    bulk_potential_mapping: Dict[str, str] = field(default_factory=dict)
    slab_potential_mapping: Optional[Dict[str, str]] = None
    metal_potential_mapping: Optional[Dict[str, str]] = None
    oxygen_potential_mapping: Optional[Dict[str, str]] = None
    nonmetal_potential_mapping: Optional[Dict[str, str]] = None
    aimd_potential_mapping: Optional[Dict[str, str]] = None
    adsorption_potential_mapping: Optional[Dict[str, str]] = None

    # Parameter configuration (INCAR parameters)
    bulk_parameters: Dict[str, Any] = field(default_factory=dict)
    slab_parameters: Optional[Dict[str, Any]] = None
    metal_parameters: Optional[Dict[str, Any]] = None
    oxygen_parameters: Optional[Dict[str, Any]] = None
    nonmetal_parameters: Optional[Dict[str, Any]] = None
    aimd_parameters: Optional[Dict[str, Any]] = None
    adsorption_parameters: Optional[Dict[str, Any]] = None
    bands_parameters: Optional[Dict[str, Any]] = None
    slab_bands_parameters: Optional[Dict[str, Any]] = None

    # Options configuration (scheduler options)
    bulk_options: Dict[str, Any] = field(default_factory=dict)
    slab_options: Optional[Dict[str, Any]] = None
    metal_options: Optional[Dict[str, Any]] = None
    oxygen_options: Optional[Dict[str, Any]] = None
    nonmetal_options: Optional[Dict[str, Any]] = None
    aimd_options: Optional[Dict[str, Any]] = None
    adsorption_options: Optional[Dict[str, Any]] = None
    bands_options: Optional[Dict[str, Any]] = None
    slab_bands_options: Optional[Dict[str, Any]] = None

    # K-points configuration
    kpoints_spacing: float = 0.03
    slab_kpoints_spacing: Optional[float] = None
    aimd_kpoints_spacing: Optional[float] = None
    adsorption_kpoints_spacing: Optional[float] = None

    # Workdir configuration
    clean_workdir: bool = False

    # Logger
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger(__name__))

    # Restart/input data
    restart_folders: Optional[Dict[str, Any]] = None  # Dict[str, orm.RemoteData]
    restart_slabs: Optional[Dict[str, Any]] = None  # Dict[str, orm.StructureData]
    input_slabs: Optional[Dict[str, Any]] = None  # Dict[str, orm.StructureData]
    use_input_slabs: bool = False

    def get_effective_code_label(self, calculation_type: str) -> str:
        """Get the effective code label for a calculation type.

        Args:
            calculation_type: Type of calculation ('bulk', 'slab', 'metal',
                'oxygen', 'nonmetal', 'aimd').

        Returns:
            The code label to use, falling back to code_label if the
            specific code label is not set.
        """
        code_map = {
            "bulk": self.bulk_code_label,
            "slab": self.slab_code_label,
            "metal": self.metal_code_label,
            "oxygen": self.oxygen_code_label,
            "nonmetal": self.nonmetal_code_label,
            "aimd": self.aimd_code_label,
        }
        specific_code = code_map.get(calculation_type)
        return specific_code if specific_code is not None else self.code_label

    def get_effective_kpoints_spacing(self, calculation_type: str) -> float:
        """Get the effective k-points spacing for a calculation type.

        Args:
            calculation_type: Type of calculation ('slab', 'aimd', 'adsorption').

        Returns:
            The k-points spacing to use, falling back to kpoints_spacing if the
            specific spacing is not set.
        """
        spacing_map = {
            "slab": self.slab_kpoints_spacing,
            "aimd": self.aimd_kpoints_spacing,
            "adsorption": self.adsorption_kpoints_spacing,
        }
        specific_spacing = spacing_map.get(calculation_type)
        return (
            specific_spacing if specific_spacing is not None else self.kpoints_spacing
        )

    def get_effective_potential_mapping(
        self, calculation_type: str
    ) -> Optional[Dict[str, str]]:
        """Get the effective potential mapping for a calculation type.

        Args:
            calculation_type: Type of calculation ('bulk', 'slab', 'metal',
                'oxygen', 'nonmetal', 'aimd', 'adsorption').

        Returns:
            The potential mapping to use, falling back through the chain:
            specific -> slab -> bulk.
        """
        mapping_map = {
            "bulk": self.bulk_potential_mapping,
            "slab": self.slab_potential_mapping,
            "metal": self.metal_potential_mapping,
            "oxygen": self.oxygen_potential_mapping,
            "nonmetal": self.nonmetal_potential_mapping,
            "aimd": self.aimd_potential_mapping,
            "adsorption": self.adsorption_potential_mapping,
        }
        specific_mapping = mapping_map.get(calculation_type)
        if specific_mapping is not None:
            return specific_mapping
        # Fallback chain: specific -> slab -> bulk
        if calculation_type in ("aimd", "adsorption"):
            if self.slab_potential_mapping is not None:
                return self.slab_potential_mapping
        return self.bulk_potential_mapping

    def get_effective_parameters(
        self, calculation_type: str
    ) -> Optional[Dict[str, Any]]:
        """Get the effective INCAR parameters for a calculation type.

        Args:
            calculation_type: Type of calculation ('bulk', 'slab', 'metal',
                'oxygen', 'nonmetal', 'aimd', 'adsorption', 'bands', 'slab_bands').

        Returns:
            The parameters to use, falling back through the chain:
            specific -> slab -> bulk.
        """
        params_map = {
            "bulk": self.bulk_parameters,
            "slab": self.slab_parameters,
            "metal": self.metal_parameters,
            "oxygen": self.oxygen_parameters,
            "nonmetal": self.nonmetal_parameters,
            "aimd": self.aimd_parameters,
            "adsorption": self.adsorption_parameters,
            "bands": self.bands_parameters,
            "slab_bands": self.slab_bands_parameters,
        }
        specific_params = params_map.get(calculation_type)
        if specific_params is not None:
            return specific_params
        # Fallback chain: specific -> slab -> bulk
        if calculation_type in ("aimd", "adsorption"):
            if self.slab_parameters is not None:
                return self.slab_parameters
        return self.bulk_parameters

    def get_effective_options(self, calculation_type: str) -> Optional[Dict[str, Any]]:
        """Get the effective scheduler options for a calculation type.

        Args:
            calculation_type: Type of calculation ('bulk', 'slab', 'metal',
                'oxygen', 'nonmetal', 'aimd', 'adsorption', 'bands', 'slab_bands').

        Returns:
            The options to use, falling back through the chain:
            specific -> slab -> bulk.
        """
        options_map = {
            "bulk": self.bulk_options,
            "slab": self.slab_options,
            "metal": self.metal_options,
            "oxygen": self.oxygen_options,
            "nonmetal": self.nonmetal_options,
            "aimd": self.aimd_options,
            "adsorption": self.adsorption_options,
            "bands": self.bands_options,
            "slab_bands": self.slab_bands_options,
        }
        specific_options = options_map.get(calculation_type)
        if specific_options is not None:
            return specific_options
        # Fallback chain: specific -> slab -> bulk
        if calculation_type in ("aimd", "adsorption"):
            if self.slab_options is not None:
                return self.slab_options
        return self.bulk_options


# =============================================================================
# Stage function imports
# =============================================================================
# Import from individual stage modules as they are implemented.

# Stage 1: Preset resolution - implemented in preset_resolution.py
from teros.core.stages.preset_resolution import resolve_workflow_preset


# Stage 2: Input validation - implemented in input_validation.py
from teros.core.stages.input_validation import (
    validate_workflow_inputs,
    log_workflow_configuration,
)


# Stage 3: Restart handling - implemented in restart_handling.py
from .restart_handling import (
    RestartData,
    handle_restart_from_node,
    prepare_input_slabs,
)

# Alias for backward compatibility with analysis.md naming
extract_restart_data = handle_restart_from_node


# Stage 7: Slab electronic properties - implemented in electronic_properties.py
from .electronic_properties import (
    add_slab_electronic_properties_stage,
)

# Alias for __all__ compatibility
add_slab_electronic_properties = add_slab_electronic_properties_stage


# Stage 8: Bulk electronic properties - implemented in electronic_properties.py
from .electronic_properties import (
    add_bulk_electronic_properties_stage,
)

# Alias for __all__ compatibility
add_bulk_electronic_properties = add_bulk_electronic_properties_stage


# Stage 9: AIMD - implemented in aimd_stage.py
from .aimd_stage import (
    resolve_aimd_parameters,
    prepare_fixed_atoms,
    determine_initial_slabs_source,
    add_aimd_stage,
)

# Alias for __all__ compatibility
add_aimd_stages = add_aimd_stage


# Stage 10: Adsorption energy - implemented in adsorption_stage.py
from .adsorption_stage import (
    resolve_adsorption_parameters,
    extract_builder_params,
    add_adsorption_energy_task,
    add_adsorption_energy_stage,
)


# =============================================================================
# Placeholder functions for stages not yet implemented
# =============================================================================


def setup_cp2k_files(*args, **kwargs):
    """Stage 4: Create CP2K basis and pseudopotential files.

    Placeholder - will be implemented in cp2k_setup.py.
    """
    raise NotImplementedError(
        "setup_cp2k_files not yet implemented. "
        "Import from teros.core.stages.cp2k_setup when available."
    )


def check_cleavage_compatibility(*args, **kwargs):
    """Stage 4: Check and potentially disable cleavage for manual slabs.

    Placeholder - will be implemented in cp2k_setup.py.
    """
    raise NotImplementedError(
        "check_cleavage_compatibility not yet implemented. "
        "Import from teros.core.stages.cp2k_setup when available."
    )


def build_core_workgraph_base(*args, **kwargs):
    """Stage 5: Build the core workgraph.

    Placeholder - will be implemented in core_build.py.
    """
    raise NotImplementedError(
        "build_core_workgraph_base not yet implemented. "
        "Import from teros.core.stages.core_build when available."
    )


def resolve_slab_parameters(*args, **kwargs):
    """Stage 6: Resolve slab parameters with fallback to bulk.

    Placeholder - will be implemented in input_slabs.py.
    """
    raise NotImplementedError(
        "resolve_slab_parameters not yet implemented. "
        "Import from teros.core.stages.input_slabs when available."
    )


def add_restart_slab_tasks(*args, **kwargs):
    """Stage 6: Add VASP tasks with restart folders for each slab.

    Placeholder - will be implemented in input_slabs.py.
    """
    raise NotImplementedError(
        "add_restart_slab_tasks not yet implemented. "
        "Import from teros.core.stages.input_slabs when available."
    )


def add_normal_slab_tasks(*args, **kwargs):
    """Stage 6: Add SCF + relax + relaxation energy tasks for input slabs.

    Placeholder - will be implemented in input_slabs.py.
    """
    raise NotImplementedError(
        "add_normal_slab_tasks not yet implemented. "
        "Import from teros.core.stages.input_slabs when available."
    )


def add_thermodynamics_task(*args, **kwargs):
    """Stage 6: Add surface thermodynamics calculation task.

    Placeholder - will be implemented in input_slabs.py.
    """
    raise NotImplementedError(
        "add_thermodynamics_task not yet implemented. "
        "Import from teros.core.stages.input_slabs when available."
    )


def add_cleavage_task(*args, **kwargs):
    """Stage 6: Add cleavage energy calculation task.

    Placeholder - will be implemented in input_slabs.py.
    """
    raise NotImplementedError(
        "add_cleavage_task not yet implemented. "
        "Import from teros.core.stages.input_slabs when available."
    )
