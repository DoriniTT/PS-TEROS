"""
Compatibility module for teros.workgraph imports.

This module re-exports functions from teros.CORE.workgraph for backward compatibility.
"""

from teros.CORE.workgraph import (
    build_core_workgraph,
    build_core_workgraph_with_map,
    core_workgraph,
    load_structure,
    extract_total_energy,
)

__all__ = [
    'build_core_workgraph',
    'build_core_workgraph_with_map',
    'core_workgraph',
    'load_structure',
    'extract_total_energy',
]
