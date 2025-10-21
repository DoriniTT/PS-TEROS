#!/usr/bin/env python
"""Detailed verification of 2D binning algorithm for code review."""

import numpy as np
from surface_modes import sample_by_coverage_bins


def test_2d_binning_edge_cases():
    """Test edge cases in 2D binning algorithm."""

    print("\n=== Test 1: Regular 5x5 grid to 2x2 bins ===")
    variants = []
    for vac in range(5):
        for oh in range(5):
            variants.append(({"name": f"v{vac}_{oh}"}, None, (float(vac), float(oh))))

    result = sample_by_coverage_bins(variants, n_bins=2, mode='combine')

    print(f"Input: {len(variants)} variants")
    print(f"Output: {len(result)} sampled variants")
    print("Selected variants:")
    for v in result:
        print(f"  {v[0]['name']}: coverage={v[2]}")

    # Verify the bin centers are approximately correct
    # For range [0,4] with 2 bins: centers should be at 1.0 and 3.0
    vac_coverages = sorted([v[2][0] for v in result])
    oh_coverages = sorted([v[2][1] for v in result])

    print(f"Vacancy coverages: {vac_coverages}")
    print(f"OH coverages: {oh_coverages}")

    print("\n=== Test 2: Single dimension with zero range ===")
    # All variants have same vacancy coverage
    variants2 = []
    for oh in range(5):
        variants2.append(({"name": f"v2_{oh}"}, None, (2.0, float(oh))))

    result2 = sample_by_coverage_bins(variants2, n_bins=2, mode='combine')

    print(f"Input: {len(variants2)} variants (all vac=2.0)")
    print(f"Output: {len(result2)} sampled variants")
    print("Selected variants:")
    for v in result2:
        print(f"  {v[0]['name']}: coverage={v[2]}")

    print("\n=== Test 3: Both dimensions with zero range ===")
    # All variants have same coverage
    variants3 = []
    for i in range(5):
        variants3.append(({"name": f"v{i}"}, None, (1.5, 2.5)))

    result3 = sample_by_coverage_bins(variants3, n_bins=2, mode='combine')

    print(f"Input: {len(variants3)} variants (all at 1.5, 2.5)")
    print(f"Output: {len(result3)} sampled variants")
    print("Selected variants:")
    for v in result3:
        print(f"  {v[0]['name']}: coverage={v[2]}")

    print("\n=== Test 4: Verify no duplicates are selected ===")
    # Create a sparse grid where some bins might try to select same variant
    variants4 = [
        ({"name": "v00"}, None, (0.0, 0.0)),
        ({"name": "v01"}, None, (0.0, 1.0)),
        ({"name": "v10"}, None, (1.0, 0.0)),
        ({"name": "v11"}, None, (1.0, 1.0)),
        ({"name": "v_center"}, None, (0.5, 0.5)),
    ]

    result4 = sample_by_coverage_bins(variants4, n_bins=3, mode='combine')

    print(f"Input: {len(variants4)} sparse variants")
    print(f"Output: {len(result4)} sampled variants")
    print("Selected variants:")
    for v in result4:
        print(f"  {v[0]['name']}: coverage={v[2]}")

    # Check no duplicates
    selected_names = [v[0]['name'] for v in result4]
    assert len(selected_names) == len(set(selected_names)), "Duplicates found!"
    print("✓ No duplicates in selection")

    print("\n=== Test 5: Verify Euclidean distance calculation ===")
    # Create variants and manually verify which should be selected
    variants5 = []
    # Create a 3x3 grid
    for vac in [0, 1, 2]:
        for oh in [0, 1, 2]:
            variants5.append(({"name": f"v{vac}{oh}"}, None, (float(vac), float(oh))))

    result5 = sample_by_coverage_bins(variants5, n_bins=2, mode='combine')

    print(f"Input: 3x3 grid of variants")
    print(f"Output: {len(result5)} sampled variants (2x2 bins)")

    # For range [0,2] with 2 bins:
    # Bin centers should be at (0.5, 0.5), (0.5, 1.5), (1.5, 0.5), (1.5, 1.5)
    expected_centers = [(0.5, 0.5), (0.5, 1.5), (1.5, 0.5), (1.5, 1.5)]

    print("\nExpected bin centers and closest variants:")
    for center in expected_centers:
        # Find closest variant to this center
        min_dist = float('inf')
        closest = None
        for v in variants5:
            vac, oh = v[2]
            dist = ((vac - center[0])**2 + (oh - center[1])**2)**0.5
            if dist < min_dist:
                min_dist = dist
                closest = v
        print(f"  Center {center} -> {closest[0]['name']} at {closest[2]} (dist={min_dist:.3f})")

    print("\nActually selected:")
    for v in result5:
        print(f"  {v[0]['name']}: coverage={v[2]}")

    print("\n=== Test 6: Large grid performance ===")
    # Create a large grid to test performance
    variants6 = []
    for vac in range(20):
        for oh in range(20):
            variants6.append(({"name": f"v{vac:02d}{oh:02d}"}, None, (float(vac), float(oh))))

    import time
    start = time.time()
    result6 = sample_by_coverage_bins(variants6, n_bins=5, mode='combine')
    elapsed = time.time() - start

    print(f"Input: {len(variants6)} variants (20x20 grid)")
    print(f"Output: {len(result6)} sampled variants (5x5 bins)")
    print(f"Time: {elapsed:.3f} seconds")

    # Verify good distribution
    vac_coverages = [v[2][0] for v in result6]
    oh_coverages = [v[2][1] for v in result6]

    print(f"Vacancy coverage range: [{min(vac_coverages):.1f}, {max(vac_coverages):.1f}]")
    print(f"OH coverage range: [{min(oh_coverages):.1f}, {max(oh_coverages):.1f}]")
    print(f"Unique vac values: {len(set(vac_coverages))}")
    print(f"Unique OH values: {len(set(oh_coverages))}")


if __name__ == "__main__":
    test_2d_binning_edge_cases()