#!/usr/bin/env python3
"""Verify surface detector implementation without pytest"""

from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
import numpy as np

# Test 1: Finding topmost Mn
print("Test 1: Finding topmost Mn atom...")
cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
structure = load_structure(cif_path)

mn_index, mn_position = find_surface_mn(structure)

print(f"  Found Mn at index: {mn_index}")
print(f"  Position (Cartesian): {mn_position}")
print(f"  Z-coordinate: {mn_position[2]:.6f}")

# Verify it's actually a Mn atom
assert str(structure[mn_index].specie) == "Mn", "Not a Mn atom!"
print("  ✓ Confirmed it's a Mn atom")

# Verify it's the highest z among all Mn atoms
mn_sites = [(i, site) for i, site in enumerate(structure) if str(site.specie) == "Mn"]
print(f"  Found {len(mn_sites)} Mn atoms total")

max_z_found = mn_position[2]
for idx, site in mn_sites:
    z_coord = site.coords[2]
    if z_coord > max_z_found + 1e-6:
        raise AssertionError(f"Found Mn at index {idx} with higher z={z_coord:.6f}")

print("  ✓ Confirmed it has the highest z-coordinate among all Mn atoms")

# Test 2: Verify Cartesian coordinates
print("\nTest 2: Verifying Cartesian coordinates...")
expected_coords = structure[mn_index].coords
if not np.allclose(mn_position, expected_coords, atol=1e-6):
    raise AssertionError("Position doesn't match structure's Cartesian coords")

print("  ✓ Position matches structure's Cartesian coordinates")

# Additional validation
print("\nAdditional details:")
print(f"  Surface Mn species: {structure[mn_index].specie}")
print(f"  Surface Mn fractional coords: {structure[mn_index].frac_coords}")

# Show z-distribution of all Mn atoms
z_coords = [site.coords[2] for _, site in mn_sites]
print(f"\n  Z-coordinate distribution of Mn atoms:")
print(f"    Min z: {min(z_coords):.6f}")
print(f"    Max z: {max(z_coords):.6f}")
print(f"    Range: {max(z_coords) - min(z_coords):.6f}")

print("\n✅ All tests passed successfully!")