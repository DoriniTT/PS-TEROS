"""
Metal Surface Energy Module

This module computes surface energies for elemental metals using the simple formula:
γ = (E_slab - N·E_bulk/atom) / (2A)

where:
- E_slab: Total energy of the slab
- N: Number of atoms in the slab
- E_bulk/atom: Bulk energy per atom
- A: Surface area
- Factor of 2 accounts for two surfaces

This is much simpler than oxide thermodynamics since there are no chemical potential
dependencies for pure metals.
"""

from .surface_energy import (
    calculate_metal_surface_energy,
    compute_metal_surface_energies_scatter,
)

from .workgraph import build_metal_surface_energy_workgraph

__all__ = [
    'calculate_metal_surface_energy',
    'compute_metal_surface_energies_scatter',
    'build_metal_surface_energy_workgraph',
]
