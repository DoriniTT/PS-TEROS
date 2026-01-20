"""PS-TEROS modules for AiiDA-WorkGraph workflows."""

from .helper_functions import (
    get_structure_from_file,
    prepare_vasp_inputs,
)
from .hf import calculate_formation_enthalpy
from .slabs import (
    generate_slab_structures,
    relax_slabs_scatter,
    extract_total_energy,
    deep_merge_dicts,  # Re-export from slabs for backward compatibility
)
from .thermodynamics import (
    identify_oxide_type,
    calculate_surface_energy_ternary,
    calculate_surface_energy_binary,
    compute_surface_energies_scatter,
)
from .cleavage import (
    calculate_cleavage_energy,
    compute_cleavage_energies_scatter,
)
from .workflow_presets import (
    list_workflow_presets,
    get_preset_config,
    get_preset_summary,
    WORKFLOW_PRESETS,
    DEFAULT_PRESET,
)
from .adsorption_energy import (
    separate_adsorbate_structure,
    calculate_adsorption_energy,
    compute_adsorption_energies_scatter,
)
from .surface_energy import (
    calculate_metal_surface_energy,
    build_metal_surface_energy_workgraph,
    identify_compound_type,
)
from .convergence import (
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
    generate_thickness_series,
    extract_recommended_layers,
)
from .u_calculation import (
    build_u_calculation_workgraph,
    get_u_calculation_results,
)
from .fukui import (
    build_fukui_workgraph,
    get_fukui_results,
    print_fukui_summary,
)
from .vasp_parallelization import (
    build_parallelization_benchmark_workgraph,
    get_benchmark_results,
    print_benchmark_summary,
    generate_benchmark_combinations,
)
from .utils import (
    deep_merge_dicts,
    get_vasp_parser_settings,
    extract_max_jobs_value,
    ensure_python_dict,
    ensure_python_float,
    ensure_python_int,
    ensure_python_bool,
    ensure_python_str,
    ensure_python_list,
    TaskOutputPlaceholder,
    EnergyOutputPlaceholder,
    FormationEnthalpyPlaceholder,
    # Structure analysis utilities
    calculate_surface_area,
    get_atom_counts,
    get_formula_units,
    get_reduced_stoichiometry,
    get_metal_elements,
)
from .constants import (
    EV_PER_ANGSTROM2_TO_J_PER_M2,
    EV_TO_KJ_PER_MOL,
    EV_TO_JOULE,
    STOICHIOMETRY_RTOL,
    DEFAULT_KPOINTS_SPACING,
)
from .exceptions import (
    TerosError,
    ValidationError,
    StructureError,
    ConvergenceError,
    ConfigurationError,
    WorkflowError,
)

__all__ = [
    # Helper functions
    "get_structure_from_file",
    "prepare_vasp_inputs",
    # Formation enthalpy
    "calculate_formation_enthalpy",
    # Slab generation and relaxation
    "generate_slab_structures",
    "relax_slabs_scatter",
    "extract_total_energy",
    # Thermodynamics
    "identify_oxide_type",
    "calculate_surface_energy_ternary",
    "calculate_surface_energy_binary",
    "compute_surface_energies_scatter",
    # Cleavage
    "calculate_cleavage_energy",
    "compute_cleavage_energies_scatter",
    # Workflow presets
    "list_workflow_presets",
    "get_preset_config",
    "get_preset_summary",
    "WORKFLOW_PRESETS",
    "DEFAULT_PRESET",
    # Adsorption
    "separate_adsorbate_structure",
    "calculate_adsorption_energy",
    "compute_adsorption_energies_scatter",
    # Metal surface energy
    "calculate_metal_surface_energy",
    "build_metal_surface_energy_workgraph",
    "identify_compound_type",
    # Thickness convergence
    "build_thickness_convergence_workgraph",
    "get_thickness_convergence_results",
    "generate_thickness_series",
    "extract_recommended_layers",
    # Hubbard U calculation
    "build_u_calculation_workgraph",
    "get_u_calculation_results",
    # Fukui function calculation
    "build_fukui_workgraph",
    "get_fukui_results",
    "print_fukui_summary",
    # VASP parallelization benchmarking
    "build_parallelization_benchmark_workgraph",
    "get_benchmark_results",
    "print_benchmark_summary",
    "generate_benchmark_combinations",
    # Utilities
    "deep_merge_dicts",
    "get_vasp_parser_settings",
    "extract_max_jobs_value",
    "ensure_python_dict",
    "ensure_python_float",
    "ensure_python_int",
    "ensure_python_bool",
    "ensure_python_str",
    "ensure_python_list",
    # Placeholder classes
    "TaskOutputPlaceholder",
    "EnergyOutputPlaceholder",
    "FormationEnthalpyPlaceholder",
    # Structure analysis utilities
    "calculate_surface_area",
    "get_atom_counts",
    "get_formula_units",
    "get_reduced_stoichiometry",
    "get_metal_elements",
    # Constants
    "EV_PER_ANGSTROM2_TO_J_PER_M2",
    "EV_TO_KJ_PER_MOL",
    "EV_TO_JOULE",
    "STOICHIOMETRY_RTOL",
    "DEFAULT_KPOINTS_SPACING",
    # Exceptions
    "TerosError",
    "ValidationError",
    "StructureError",
    "ConvergenceError",
    "ConfigurationError",
    "WorkflowError",
]
