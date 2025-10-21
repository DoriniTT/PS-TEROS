#!/usr/bin/env python
"""Verification script to test the structure loader implementation."""

from pathlib import Path
from structure_loader import load_structure, validate_structure


def test_load_structure():
    """Test loading CIF file."""
    print("\n=== Testing load_structure() ===")
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"

    try:
        structure = load_structure(cif_path)
        print(f"✓ Successfully loaded CIF file")
        print(f"  - Number of atoms: {len(structure)}")
        print(f"  - Lattice parameters: a={structure.lattice.a:.3f}, b={structure.lattice.b:.3f}, c={structure.lattice.c:.3f}")
        return structure
    except Exception as e:
        print(f"✗ Failed to load CIF file: {e}")
        return None


def test_validate_structure(structure):
    """Test structure validation."""
    print("\n=== Testing validate_structure() ===")

    if structure is None:
        print("✗ Cannot test validation without structure")
        return

    result = validate_structure(structure)

    print(f"Validation results:")
    print(f"  - Valid: {result['valid']}")
    print(f"  - Elements present: {sorted(result['elements'])}")
    print(f"  - Number of Mn atoms: {result['num_mn']}")
    print(f"  - Number of O atoms: {result['num_o']}")

    # Check requirements
    if result['valid']:
        print("✓ Structure validation passed - contains Mn and O")
    else:
        print("✗ Structure validation failed - missing required elements")

    if result['num_mn'] == 7 and result['num_o'] == 22:
        print("✓ Element counts match expected values (7 Mn, 22 O)")
    else:
        print(f"⚠ Element counts differ from expected (expected 7 Mn, 22 O)")


def test_error_handling():
    """Test error handling for invalid file."""
    print("\n=== Testing error handling ===")

    try:
        load_structure("nonexistent.cif")
        print("✗ Failed to raise FileNotFoundError")
    except FileNotFoundError as e:
        print(f"✓ Correctly raised FileNotFoundError: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")


if __name__ == "__main__":
    print("Verifying structure_loader implementation...")

    # Test 1: Load structure
    structure = test_load_structure()

    # Test 2: Validate structure
    test_validate_structure(structure)

    # Test 3: Error handling
    test_error_handling()

    print("\n=== Verification complete ===")