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

Wulff Shape Generation:
The module also supports Wulff shape visualization, which shows the equilibrium
crystal shape based on the computed surface energies. Enable with
`generate_wulff_shape=True` in build_metal_surface_energy_workgraph().
"""

from .surface_energy import (
    calculate_metal_surface_energy,
    compute_metal_surface_energies_scatter,
    identify_compound_type,
)

from .workgraph import build_metal_surface_energy_workgraph

from .wulff import (
    generate_wulff_shape_data,
    select_lowest_energy_terminations,
    create_wulff_shape,
)

__all__ = [
    # Surface energy calculations
    'calculate_metal_surface_energy',
    'compute_metal_surface_energies_scatter',
    'identify_compound_type',
    'build_metal_surface_energy_workgraph',
    # Wulff shape generation
    'generate_wulff_shape_data',
    'select_lowest_energy_terminations',
    'create_wulff_shape',
]
