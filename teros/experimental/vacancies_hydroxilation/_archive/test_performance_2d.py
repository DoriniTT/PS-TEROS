#!/usr/bin/env python
"""Performance testing for 2D binning implementation."""

import time
import random
from surface_modes import sample_by_coverage_bins


def test_performance():
    """Test performance with various input sizes."""

    print("=== PERFORMANCE TESTING FOR 2D BINNING ===\n")

    test_cases = [
        (10, 10, 3),     # 100 variants, 9 bins
        (20, 20, 5),     # 400 variants, 25 bins
        (50, 50, 7),     # 2500 variants, 49 bins
        (100, 100, 10),  # 10000 variants, 100 bins
    ]

    results = []

    for grid_size_x, grid_size_y, n_bins in test_cases:
        # Create variants
        variants = []
        for vac in range(grid_size_x):
            for oh in range(grid_size_y):
                coverage = (float(vac) + random.random() * 0.1,
                            float(oh) + random.random() * 0.1)
                variants.append(({"name": f"v{vac}_{oh}"}, None, coverage))

        # Time the binning
        start_time = time.perf_counter()
        result = sample_by_coverage_bins(variants, n_bins=n_bins, mode='combine')
        elapsed = time.perf_counter() - start_time

        n_variants = len(variants)
        n_selected = len(result)
        max_possible = min(n_bins * n_bins, n_variants)

        print(f"Test: {n_variants} variants → {n_bins}×{n_bins} bins")
        print(f"  Input size: {n_variants}")
        print(f"  Bins requested: {n_bins}×{n_bins} = {n_bins*n_bins}")
        print(f"  Selected: {n_selected}/{max_possible} possible")
        print(f"  Time: {elapsed:.4f} seconds")
        print(f"  Time per variant: {elapsed/n_variants*1000:.3f} ms")
        print(f"  Time per bin: {elapsed/(n_bins*n_bins)*1000:.3f} ms")
        print()

        results.append({
            'variants': n_variants,
            'bins': n_bins * n_bins,
            'time': elapsed,
            'selected': n_selected
        })

    # Performance analysis
    print("=" * 50)
    print("PERFORMANCE ANALYSIS")
    print("=" * 50)

    # Check if time complexity is reasonable
    # Should be O(n_bins² × n_variants)
    print("\nTime Complexity Check:")
    for i in range(1, len(results)):
        prev = results[i-1]
        curr = results[i]

        # Calculate expected time based on complexity
        expected_factor = (curr['bins'] * curr['variants']) / (prev['bins'] * prev['variants'])
        actual_factor = curr['time'] / prev['time']

        print(f"  {prev['variants']} → {curr['variants']} variants:")
        print(f"    Expected time factor: {expected_factor:.2f}x")
        print(f"    Actual time factor: {actual_factor:.2f}x")

        if actual_factor < expected_factor * 1.5:  # Allow 50% overhead
            print(f"    ✓ Performance scaling is good")
        else:
            print(f"    ⚠ Performance degradation detected")

    # Check selection completeness
    print("\nSelection Completeness:")
    for r in results:
        max_possible = min(r['bins'], r['variants'])
        completeness = r['selected'] / max_possible * 100
        print(f"  {r['variants']} variants, {r['bins']} bins: "
              f"{r['selected']}/{max_possible} = {completeness:.1f}%")

    # Memory estimation
    print("\nMemory Usage Estimation:")
    for r in results:
        # Each variant tuple: ~200 bytes (dict + atoms object placeholder + coverage tuple)
        # sampled_indices set: 8 bytes per index
        # sampled list: 8 bytes per pointer
        memory_input = r['variants'] * 200  # bytes
        memory_indices = r['selected'] * 8  # bytes
        memory_output = r['selected'] * 8  # bytes
        total_memory = memory_input + memory_indices + memory_output

        print(f"  {r['variants']} variants: ~{total_memory/1024:.1f} KB")

    return all(r['time'] < 1.0 for r in results)  # All should complete in < 1 second


def test_worst_case():
    """Test worst-case scenario: all variants at same point."""

    print("\n" + "=" * 50)
    print("WORST CASE TESTING")
    print("=" * 50)

    # Create many variants at the same coverage point
    variants = []
    for i in range(1000):
        variants.append(({"name": f"v{i}"}, None, (5.0, 5.0)))

    print(f"\nTest: {len(variants)} variants all at same point (5.0, 5.0)")

    start_time = time.perf_counter()
    result = sample_by_coverage_bins(variants, n_bins=10, mode='combine')
    elapsed = time.perf_counter() - start_time

    print(f"  Bins requested: 10×10 = 100")
    print(f"  Selected: {len(result)}")
    print(f"  Time: {elapsed:.4f} seconds")

    if elapsed < 0.5:
        print("  ✓ Handles worst case efficiently")
    else:
        print("  ⚠ Performance issue in worst case")


if __name__ == "__main__":
    print("Running performance tests...\n")

    success = test_performance()
    test_worst_case()

    print("\n" + "=" * 50)
    if success:
        print("RESULT: Performance is ACCEPTABLE ✓")
    else:
        print("RESULT: Performance needs optimization")
    print("=" * 50)