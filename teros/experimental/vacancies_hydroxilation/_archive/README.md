# Archive

This directory contains old test outputs, analysis files, and experimental code from previous development iterations.

## Contents

### `old_test_outputs/`
Old test output directories from coverage calculations and binning verification:
- `test_vac_cov/`, `test_h_cov/`, `test_combo_cov/` - Coverage calculation test outputs
- `test_vac_dedup/`, `test_h_dedup/`, `test_combo_dedup/` - Deduplication test outputs
- `all_possibilities/` - Full enumeration without optimization
- `output/`, `outputs/`, `test_output/` - Various test run outputs
- `test_temp/`, `surf_variants/` - Temporary files from analysis

### Analysis Files
- `TASK3_CODE_REVIEW.md` - Code review notes from binning implementation
- `test_2d_algorithm_correctness.py` - 2D binning verification tests
- `test_2d_binning_verification.py` - Binning correctness validation
- `test_performance_2d.py` - Performance analysis
- `static_analysis_2d_binning.py` - Static analysis of 2D binning algorithm

## When to Use

Refer to archived content if:
- Debugging specific test scenarios from earlier runs
- Reviewing the development history of binning algorithms
- Understanding the full enumeration results (before optimization)

Otherwise, use the active production code:
- **Tests**: See `tests/test_surface_modes.py` (35 unit tests, all passing)
- **Integration tests**: See `tests/test_integration_*.sh`
- **Main script**: See `surface_modes.py` in root directory
