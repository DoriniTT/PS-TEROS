"""Custom VASP calculation module for PS-TEROS."""

from .workgraph import (
    build_custom_calculation_workgraph,
    get_custom_results,
    build_dos_calculation_workgraph,
    get_dos_results,
)

__all__ = [
    # Standard VASP calculations
    'build_custom_calculation_workgraph',
    'get_custom_results',
    # DOS calculations
    'build_dos_calculation_workgraph',
    'get_dos_results',
]
