"""
Thermodynamics module for ab initio atomistic thermodynamics calculations.

Contains functions for calculating formation enthalpies and surface energies
for binary and ternary oxide systems.
"""

from teros.functions.thermodynamics.formation import calculate_formation_enthalpy
from teros.functions.thermodynamics.binary import calculate_surface_energy_binary
from teros.functions.thermodynamics.ternary import calculate_surface_energy_ternary

__all__ = [
    'calculate_formation_enthalpy',
    'calculate_surface_energy_binary',
    'calculate_surface_energy_ternary',
]