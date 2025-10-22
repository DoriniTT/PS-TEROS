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

__all__ = [
    'get_structure_from_file',
    'prepare_vasp_inputs',
    'calculate_formation_enthalpy',
    'generate_slab_structures',
    'relax_slabs_scatter',
    'extract_total_energy',
    'identify_oxide_type',
    'calculate_surface_energy_ternary',
    'calculate_surface_energy_binary',
    'compute_surface_energies_scatter',
    'calculate_cleavage_energy',
    'compute_cleavage_energies_scatter',
    'list_workflow_presets',
    'get_preset_config',
    'get_preset_summary',
    'WORKFLOW_PRESETS',
    'DEFAULT_PRESET',
    'separate_adsorbate_structure',
    'calculate_adsorption_energy',
    'compute_adsorption_energies_scatter',
]
