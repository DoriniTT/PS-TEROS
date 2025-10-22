#!/usr/bin/env python3
"""
Visualize LaMnO3 surface structure with OOH using ASE.
"""

from ase.io import read
from ase.visualize import view

# Read the CIF file
structure = read('LaMnO3_100_A4_surface_OOH.cif')

# Print basic information
print(f"Chemical formula: {structure.get_chemical_formula()}")
print(f"Number of atoms: {len(structure)}")
print(f"Cell parameters:")
print(f"  a = {structure.cell.lengths()[0]:.3f} Å")
print(f"  b = {structure.cell.lengths()[1]:.3f} Å")
print(f"  c = {structure.cell.lengths()[2]:.3f} Å")
print(f"  α = {structure.cell.angles()[0]:.2f}°")
print(f"  β = {structure.cell.angles()[1]:.2f}°")
print(f"  γ = {structure.cell.angles()[2]:.2f}°")

# Visualize the structure
view(structure)
