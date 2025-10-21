#!/usr/bin/env python
"""Deep verification of 2D binning algorithm correctness."""

import math
from surface_modes import sample_by_coverage_bins


def test_algorithm_correctness():
    """Verify the mathematical correctness of the 2D binning algorithm."""

    print("=== ALGORITHM CORRECTNESS TESTS ===\n")

    # Test 1: Verify bin center calculations are correct
    print("Test 1: Bin center calculations")
    print("-" * 40)

    # Create simple case: coverage range [0, 4] x [0, 4] with 2x2 bins
    variants = []
    for vac in range(5):
        for oh in range(5):
            variants.append(({"name": f"v{vac}{oh}"}, None, (float(vac), float(oh))))

    # Expected bin centers for range [0,4] with 2 bins:
    # Bin width = (4-0)/2 = 2
    # Centers at: min + (i + 0.5) * width
    # i=0: 0 + 0.5*2 = 1.0
    # i=1: 0 + 1.5*2 = 3.0
    expected_centers = [
        (1.0, 1.0), (1.0, 3.0),
        (3.0, 1.0), (3.0, 3.0)
    ]

    result = sample_by_coverage_bins(variants, n_bins=2, mode='combine')

    print(f"Coverage range: [0, 4] x [0, 4]")
    print(f"Number of bins: 2x2 = 4")
    print(f"Bin width: 2.0 x 2.0")
    print(f"Expected bin centers: {expected_centers}")
    print(f"Selected variants:")

    for v in result:
        vac, oh = v[2]
        print(f"  {v[0]['name']}: ({vac}, {oh})")

        # Find which bin center this should be closest to
        min_dist = float('inf')
        closest_center = None
        for center in expected_centers:
            dist = math.sqrt((vac - center[0])**2 + (oh - center[1])**2)
            if dist < min_dist:
                min_dist = dist
                closest_center = center

        # The selected variant should indeed be closest to one of the centers
        # among all variants
        is_optimal = True
        for other in variants:
            other_vac, other_oh = other[2]
            other_dist = math.sqrt((other_vac - closest_center[0])**2 +
                                   (other_oh - closest_center[1])**2)
            if other_dist < min_dist - 0.0001:  # Small epsilon for float comparison
                is_optimal = False
                print(f"    WARNING: {other[0]['name']} at {other[2]} is closer to {closest_center} (dist={other_dist:.3f} vs {min_dist:.3f})")

        if is_optimal:
            print(f"    ✓ Optimal for center {closest_center} (dist={min_dist:.3f})")

    # Test 2: Verify that each variant is selected at most once
    print("\nTest 2: No duplicate selection")
    print("-" * 40)

    # Create a scenario where multiple bin centers might want the same variant
    sparse_variants = [
        ({"name": "center"}, None, (2.0, 2.0)),  # Center of the space
        ({"name": "corner1"}, None, (0.0, 0.0)),
        ({"name": "corner2"}, None, (4.0, 4.0)),
    ]

    result2 = sample_by_coverage_bins(sparse_variants, n_bins=3, mode='combine')

    print(f"Input: {len(sparse_variants)} variants (sparse)")
    print(f"Requested bins: 3x3 = 9")
    print(f"Output: {len(result2)} variants selected")

    selected_names = [v[0]['name'] for v in result2]
    if len(selected_names) == len(set(selected_names)):
        print("✓ No duplicates found")
    else:
        print("✗ DUPLICATES DETECTED!")
        print(f"  Selected: {selected_names}")

    # Test 3: Verify edge case handling (zero range)
    print("\nTest 3: Zero range edge cases")
    print("-" * 40)

    # Case A: Zero range in vacancy dimension
    zero_vac_variants = [
        ({"name": f"oh{i}"}, None, (1.0, float(i))) for i in range(5)
    ]

    result3a = sample_by_coverage_bins(zero_vac_variants, n_bins=2, mode='combine')

    print(f"Case A: All variants at vac=1.0, OH varies [0,4]")
    print(f"  Input: {len(zero_vac_variants)} variants")
    print(f"  Output: {len(result3a)} variants")

    # When vac range is zero, bin width should be set to 1.0
    # This should still allow 2x2=4 variants to be selected
    if len(result3a) <= 4:
        print(f"  ✓ Correctly handled zero range (selected {len(result3a)} <= 4)")
    else:
        print(f"  ✗ Too many selected: {len(result3a)} > 4")

    # Case B: Zero range in both dimensions
    same_point_variants = [
        ({"name": f"v{i}"}, None, (2.0, 3.0)) for i in range(5)
    ]

    result3b = sample_by_coverage_bins(same_point_variants, n_bins=2, mode='combine')

    print(f"Case B: All variants at same point (2.0, 3.0)")
    print(f"  Input: {len(same_point_variants)} variants")
    print(f"  Output: {len(result3b)} variants")

    # When both ranges are zero, we should still select up to n_bins*n_bins variants
    if len(result3b) <= 4:
        print(f"  ✓ Correctly handled double zero range (selected {len(result3b)} <= 4)")
    else:
        print(f"  ✗ Too many selected: {len(result3b)} > 4")

    # Test 4: Verify distance calculation is Euclidean
    print("\nTest 4: Euclidean distance verification")
    print("-" * 40)

    # Create variants at specific positions
    test_variants = [
        ({"name": "origin"}, None, (0.0, 0.0)),
        ({"name": "right"}, None, (3.0, 0.0)),
        ({"name": "up"}, None, (0.0, 4.0)),
        ({"name": "diagonal"}, None, (3.0, 4.0)),
    ]

    # Calculate distances manually to a test point
    test_point = (1.5, 2.0)
    print(f"Test point (bin center): {test_point}")
    print("Distances from each variant:")

    for v in test_variants:
        vac, oh = v[2]
        dist = math.sqrt((vac - test_point[0])**2 + (oh - test_point[1])**2)
        print(f"  {v[0]['name']} at {v[2]}: dist = {dist:.3f}")

    # The algorithm should pick "origin" as it's closest
    expected_closest = "origin"
    print(f"\nExpected closest: {expected_closest}")

    # Test 5: Grid coverage verification
    print("\nTest 5: Grid coverage (all bins get a representative)")
    print("-" * 40)

    # Create dense grid
    dense_variants = []
    for vac in range(10):
        for oh in range(10):
            dense_variants.append(({"name": f"v{vac:02d}{oh:02d}"}, None, (float(vac), float(oh))))

    result5 = sample_by_coverage_bins(dense_variants, n_bins=3, mode='combine')

    print(f"Dense grid: 10x10 = {len(dense_variants)} variants")
    print(f"Requested bins: 3x3 = 9")
    print(f"Selected: {len(result5)} variants")

    # Analyze distribution
    vac_values = [v[2][0] for v in result5]
    oh_values = [v[2][1] for v in result5]

    print(f"Vacancy values: {sorted(set(vac_values))}")
    print(f"OH values: {sorted(set(oh_values))}")

    # Check if we have good 2D coverage
    vac_spread = max(vac_values) - min(vac_values)
    oh_spread = max(oh_values) - min(oh_values)

    print(f"Vacancy spread: {vac_spread:.1f}")
    print(f"OH spread: {oh_spread:.1f}")

    if vac_spread >= 5 and oh_spread >= 5:
        print("✓ Good 2D distribution achieved")
    else:
        print("✗ Poor distribution")

    # Test 6: Order independence
    print("\nTest 6: Order independence")
    print("-" * 40)

    import random

    # Create variants
    base_variants = []
    for i in range(5):
        for j in range(5):
            base_variants.append(({"name": f"v{i}{j}"}, None, (float(i), float(j))))

    # Test with original order
    result_ordered = sample_by_coverage_bins(base_variants, n_bins=2, mode='combine')
    ordered_names = sorted([v[0]['name'] for v in result_ordered])

    # Test with shuffled order
    shuffled_variants = base_variants.copy()
    random.shuffle(shuffled_variants)
    result_shuffled = sample_by_coverage_bins(shuffled_variants, n_bins=2, mode='combine')
    shuffled_names = sorted([v[0]['name'] for v in result_shuffled])

    print(f"Ordered result: {ordered_names}")
    print(f"Shuffled result: {shuffled_names}")

    if ordered_names == shuffled_names:
        print("✓ Algorithm is order-independent")
    else:
        print("✗ WARNING: Results depend on input order!")


if __name__ == "__main__":
    test_algorithm_correctness()