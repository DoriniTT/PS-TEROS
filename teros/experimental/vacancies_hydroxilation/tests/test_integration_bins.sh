#!/bin/bash
set -e

echo "=== Integration Test: Coverage Bins Sampling with st2.vasp ==="

# Test with vacancies mode and 5 bins
echo "Test 1: Vacancies mode with 5 bins"
python surface_modes.py st2.vasp --mode vacancies \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 5 \
    --outdir output/test_bins_vac

# Check that output exists
if [ ! -f output/test_bins_vac/manifest.json ]; then
    echo "FAIL: manifest.json not created"
    exit 1
fi

# Check that binning stats exist
if ! grep -q "binning_stats" output/test_bins_vac/manifest.json; then
    echo "FAIL: binning_stats not found in manifest"
    exit 1
fi

# Count variants - should be <= 5
variant_count=$(grep -o '"name"' output/test_bins_vac/manifest.json | wc -l)
if [ "$variant_count" -gt 5 ]; then
    echo "FAIL: Too many variants ($variant_count > 5)"
    exit 1
fi

echo "PASS: Vacancies mode produced $variant_count variants (≤5)"

# Test with hydrogen mode and 5 bins
echo ""
echo "Test 2: Hydrogen mode with 5 bins"
python surface_modes.py st2.vasp --mode hydrogen \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 5 \
    --outdir output/test_bins_H

variant_count=$(grep -o '"name"' output/test_bins_H/manifest.json | wc -l)
if [ "$variant_count" -gt 5 ]; then
    echo "FAIL: Too many variants ($variant_count > 5)"
    exit 1
fi

echo "PASS: Hydrogen mode produced $variant_count variants (≤5)"

# Test with combined mode and 3 bins (3x3=9 max)
echo ""
echo "Test 3: Combined mode with 3 bins (3x3 grid)"
python surface_modes.py st2.vasp --mode combine \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 3 \
    --outdir output/test_bins_combine

variant_count=$(grep -o '"name"' output/test_bins_combine/manifest.json | wc -l)
if [ "$variant_count" -gt 9 ]; then
    echo "FAIL: Too many variants ($variant_count > 9)"
    exit 1
fi

echo "PASS: Combined mode produced $variant_count variants (≤9)"

# Verify total reduction
echo ""
echo "=== Summary ==="
echo "All integration tests passed!"
echo "Coverage bins sampling working correctly"

# Cleanup
rm -rf output/test_bins_*
