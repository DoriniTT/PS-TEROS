#!/bin/bash
# Integration test for st2.vasp with 2x2x1 supercell and coverage deduplication

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check st2.vasp exists
if [ ! -f "../experimental/st2.vasp" ]; then
    echo "ERROR: st2.vasp not found at ../experimental/st2.vasp"
    exit 1
fi

echo "=== Integration Test: st2.vasp with 2x2x1 supercell and coverage deduplication ==="

# Test 1: Complete mode with supercell and deduplication
echo ""
echo "Test 1: Complete mode with 2x2x1 supercell and deduplication"
rm -rf test_st2_complete
python surface_modes.py ../experimental/st2.vasp \
    --mode complete \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --which-surface top \
    --outdir test_st2_complete

if [ ! -f "test_st2_complete/manifest.json" ]; then
    echo "FAIL: manifest.json not created"
    exit 1
fi

if ! grep -q "surface_area_nm2" test_st2_complete/manifest_vacancies.json; then
    echo "FAIL: Coverage metadata missing from manifest"
    exit 1
fi

echo "PASS: Complete mode with deduplication"

# Test 2: Vacancies mode only
echo ""
echo "Test 2: Vacancies mode with 2x2x1 supercell"
rm -rf test_st2_vac
python surface_modes.py ../experimental/st2.vasp \
    --mode vacancies \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --which-surface top \
    --outdir test_st2_vac

if [ ! -f "test_st2_vac/manifest.json" ]; then
    echo "FAIL: manifest.json not created"
    exit 1
fi

echo "PASS: Vacancies mode"

# Test 3: Combined mode
echo ""
echo "Test 3: Combined mode with 2x2x1 supercell"
rm -rf test_st2_combo
python surface_modes.py ../experimental/st2.vasp \
    --mode combine \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --which-surface top \
    --outdir test_st2_combo

if [ ! -f "test_st2_combo/manifest.json" ]; then
    echo "FAIL: manifest.json not created"
    exit 1
fi

if ! grep -q "vacancy_coverage" test_st2_combo/manifest.json; then
    echo "FAIL: vacancy_coverage missing in combined mode"
    exit 1
fi

if ! grep -q "OH_coverage" test_st2_combo/manifest.json; then
    echo "FAIL: OH_coverage missing in combined mode"
    exit 1
fi

echo "PASS: Combined mode"

# Cleanup
echo ""
echo "Cleaning up test outputs..."
rm -rf test_st2_complete test_st2_vac test_st2_combo

echo ""
echo "=== All integration tests PASSED ==="
