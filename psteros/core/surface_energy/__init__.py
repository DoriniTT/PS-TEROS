"""
Metal and Intermetallic Surface Energy Module

This module computes surface energies for elemental metals and stoichiometric
intermetallics using the simple formula:
γ = (E_slab - N·E_bulk/atom) / (2A)

where:
- E_slab: Total energy of the slab
- N: Number of atoms in the slab
- E_bulk/atom: Bulk energy per atom
- A: Surface area
- Factor of 2 accounts for two surfaces

Supported materials:
- Elemental metals: Au, Ag, Cu, Pt, Pd, Ni, Fe, etc.
- Stoichiometric intermetallics: PdIn, AuCu, NiAl, Cu3Au, etc.

For stoichiometric and symmetric surfaces, no chemical potential dependencies
are needed. For non-stoichiometric surfaces, use the oxide thermodynamics
approach instead.
"""

from .surface_energy import (
    calculate_metal_surface_energy,
    compute_metal_surface_energies_scatter,
    identify_compound_type,
)

from .workgraph import build_metal_surface_energy_workgraph

from .wulff import (
    build_wulff_shape,
    get_symmetrically_equivalent_miller_indices,
    expand_surface_energies_with_symmetry,
    visualize_wulff_shape,
    get_wulff_shape_summary,
)

from .stoichiometric_finder import (
    # Data classes
    SlabSearchResult,
    MillerFeasibilityReport,
    NoStoichiometricSymmetricSurfaceError,
    # Core functions
    find_stoichiometric_symmetric_slabs,
    analyze_miller_feasibility,
    generate_stoichiometric_symmetric_slabs,
    check_slab_properties,
    detect_bonds,
    get_bulk_composition,
    get_feasibility_summary,
    # Constants
    DEFAULT_STRATEGIES,
)

__all__ = [
    # Surface energy calculations
    'calculate_metal_surface_energy',
    'compute_metal_surface_energies_scatter',
    'identify_compound_type',
    'build_metal_surface_energy_workgraph',
    # Wulff shape
    'build_wulff_shape',
    'get_symmetrically_equivalent_miller_indices',
    'expand_surface_energies_with_symmetry',
    'visualize_wulff_shape',
    'get_wulff_shape_summary',
    # Stoichiometric+Symmetric Surface Finder (EXPERIMENTAL)
    'SlabSearchResult',
    'MillerFeasibilityReport',
    'NoStoichiometricSymmetricSurfaceError',
    'find_stoichiometric_symmetric_slabs',
    'analyze_miller_feasibility',
    'generate_stoichiometric_symmetric_slabs',
    'check_slab_properties',
    'detect_bonds',
    'get_bulk_composition',
    'get_feasibility_summary',
    'DEFAULT_STRATEGIES',
]
