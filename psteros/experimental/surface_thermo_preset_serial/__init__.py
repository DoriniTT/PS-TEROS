"""
Serial Surface Thermodynamics Preset

Experimental flat-graph implementation where all VASP nodes exist at the same
graph level, allowing max_concurrent_jobs to control concurrent execution.
"""

from .workgraph import (
    surface_thermodynamics_serial_workgraph,
    generate_slabs_from_relaxed_bulk_workgraph,
)

__all__ = [
    'surface_thermodynamics_serial_workgraph',
    'generate_slabs_from_relaxed_bulk_workgraph',
]
