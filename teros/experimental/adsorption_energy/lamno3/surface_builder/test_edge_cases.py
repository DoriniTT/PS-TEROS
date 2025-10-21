#!/usr/bin/env python3
"""Test edge cases for surface detector"""

from pymatgen.core import Structure, Lattice
from surface_detector import find_surface_mn
import numpy as np

# Test with no Mn atoms
print("Testing with structure containing no Mn atoms...")
try:
    lattice = Lattice.cubic(10.0)
    no_mn_structure = Structure(
        lattice,
        ["O", "La"],
        [[0, 0, 0], [0.5, 0.5, 0.5]]
    )
    find_surface_mn(no_mn_structure)
    print("  ERROR: Should have raised ValueError!")
except ValueError as e:
    print(f"  ✓ Correctly raised ValueError: {e}")
except Exception as e:
    print(f"  ERROR: Unexpected exception: {e}")

# Test with single Mn atom
print("\nTesting with single Mn atom...")
single_mn_structure = Structure(
    Lattice.cubic(10.0),
    ["Mn"],
    [[0.25, 0.25, 0.75]]
)
idx, pos = find_surface_mn(single_mn_structure)
print(f"  Found Mn at index: {idx}")
print(f"  Position: {pos}")
assert idx == 0, "Should return index 0 for single Mn"
assert np.allclose(pos, [2.5, 2.5, 7.5]), "Position incorrect"
print("  ✓ Correctly handles single Mn atom")

# Test with multiple Mn at same z
print("\nTesting with multiple Mn atoms at same z-coordinate...")
lattice = Lattice.cubic(10.0)
same_z_structure = Structure(
    lattice,
    ["Mn", "Mn", "Mn"],
    [[0, 0, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]]  # All at z=0.5
)
idx, pos = find_surface_mn(same_z_structure)
print(f"  Found Mn at index: {idx}")
print(f"  Position: {pos}")
# Should return one of them (implementation dependent on iteration order)
assert idx in [0, 1, 2], "Should return valid index"
assert abs(pos[2] - 5.0) < 1e-6, "Should have z=5.0"
print("  ✓ Correctly handles multiple Mn at same z-level")

print("\n✅ All edge case tests passed!")