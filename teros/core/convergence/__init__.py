"""Convergence testing module for PS-TEROS.

This module provides:
1. Automated ENCUT and k-points convergence testing using aiida-vasp
2. Slab thickness convergence testing for surface energy calculations

Usage:
    # ENCUT/k-points convergence
    from teros.core.convergence import build_convergence_workgraph
    wg = build_convergence_workgraph(structure=..., code_label=..., ...)

    # Thickness convergence
    from teros.core.convergence import build_thickness_convergence_workgraph
    wg = build_thickness_convergence_workgraph(
        bulk_structure_path='/path/to/bulk.cif',
        miller_indices=[1, 1, 1],
        layer_counts=[3, 5, 7, 9, 11],
        ...
    )
"""

from .workgraph import (
    # ENCUT/k-points convergence
    build_convergence_workgraph,
    get_convergence_results,
    # Thickness convergence
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
    extract_total_energy,
    calculate_surface_energy,
    analyze_thickness_convergence,
    relax_thickness_series,
    compute_surface_energies,
    gather_surface_energies,
)
from .slabs import (
    generate_thickness_series,
    extract_recommended_layers,
)

__all__ = [
    # ENCUT/k-points convergence
    'build_convergence_workgraph',
    'get_convergence_results',
    # Thickness convergence
    'build_thickness_convergence_workgraph',
    'get_thickness_convergence_results',
    'generate_thickness_series',
    'extract_recommended_layers',
    'extract_total_energy',
    'calculate_surface_energy',
    'analyze_thickness_convergence',
    'relax_thickness_series',
    'compute_surface_energies',
    'gather_surface_energies',
]
