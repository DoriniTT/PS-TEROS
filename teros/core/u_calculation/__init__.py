"""Hubbard U calculation module for PS-TEROS.

This module provides functionality to calculate the Hubbard U parameter
for LSDA+U calculations using VASP's linear response method.

The workflow follows the method described in the VASP wiki:
https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U

Main entry point:
    build_u_calculation_workgraph() - Builds a WorkGraph for U calculation

Example:
    >>> from teros.core.u_calculation import build_u_calculation_workgraph
    >>> wg = build_u_calculation_workgraph(
    ...     structure=nio_structure,
    ...     code_label='VASP-6.5.1@cluster',
    ...     potential_family='PBE.54',
    ...     potential_mapping={'Ni': 'Ni', 'O': 'O'},
    ...     target_species='Ni',
    ... )
    >>> wg.submit()
"""

from .workgraph import (
    build_u_calculation_workgraph,
    get_u_calculation_results,
)
from .tasks import (
    extract_d_electron_occupation,
    calculate_occupation_response,
    calculate_hubbard_u_linear_regression,
    calculate_hubbard_u_single_point,
    gather_responses,
)
from .utils import (
    linear_regression,
    prepare_ground_state_incar,
    prepare_response_incar,
    build_ldau_arrays,
    validate_target_species,
    get_species_order_from_structure,
    DEFAULT_POTENTIAL_VALUES,
)

__all__ = [
    # Main entry points
    'build_u_calculation_workgraph',
    'get_u_calculation_results',
    # Task functions
    'extract_d_electron_occupation',
    'calculate_occupation_response',
    'calculate_hubbard_u_linear_regression',
    'calculate_hubbard_u_single_point',
    'gather_responses',
    # Utility functions
    'linear_regression',
    'prepare_ground_state_incar',
    'prepare_response_incar',
    'build_ldau_arrays',
    'validate_target_species',
    'get_species_order_from_structure',
    'DEFAULT_POTENTIAL_VALUES',
]
