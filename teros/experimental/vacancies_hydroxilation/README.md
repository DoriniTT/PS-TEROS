# Surface Modification with Coverage-Driven Optimization

A Python tool for generating and optimizing surface modification variants for computational chemistry studies. Reduces combinatorial explosion from thousands to dozens of structures while preserving all unique surface chemistries.

## 🎯 Overview

`surface_modes.py` generates surface-modified structures (vacancies, hydroxylation, combined modifications) with **two-stage optimization**:

1. **Coverage Deduplication**: Reduces identical surface coverages (species/nm²)
2. **Bin Sampling**: Strategically samples coverage space (0-100%)

**Result**: 99.5% reduction (6560 → 35 structures) while maintaining full coverage diversity.

## ✨ Features

### Core Modes
- **Vacancies**: Remove surface oxygen atoms
- **Hydrogen (OH groups)**: Add OH termination to surface oxygen
- **Combined**: Vacancies + OH groups simultaneously
- **Complete**: Runs all three modes

### Optimization Features
- **Supercell Generation**: Automatic N×M×L supercell creation
- **Coverage Deduplication**: Keep one structure per unique coverage value
- **Bin Sampling**: Sample N representative structures across coverage range
  - 1D binning for single modifications (vacancies, OH)
  - 2D grid binning for combined modifications
- **Configurable Precision**: Control coverage rounding (default: 4 decimal places)

### Output Features
- Structures in VASP, CIF, or XYZ format
- JSON manifests with complete metadata
- Coverage values (species/nm²) for all variants
- Bin metadata (bin_id, bin_center)
- Deduplication and binning statistics

## 📋 Requirements

```bash
python >= 3.10
ase >= 3.22
numpy
```

## 🚀 Quick Start

### Basic Usage

```bash
# Generate all modifications without optimization
python surface_modes.py input.vasp complete --which-surface top

# With coverage deduplication (96% reduction)
python surface_modes.py input.vasp complete --supercell 2 2 1 \
    --deduplicate-by-coverage

# With deduplication + binning (99.5% reduction)
python surface_modes.py input.vasp complete --supercell 2 2 1 \
    --deduplicate-by-coverage --coverage-bins 5
```

### Single Mode Examples

```bash
# Vacancies only, 5 bins
python surface_modes.py structure.vasp vacancies \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 5 \
    --outdir output/vac

# Hydroxylation, 10 bins
python surface_modes.py structure.vasp hydrogen \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 10

# Combined mode, 3×3 grid (9 bins max)
python surface_modes.py structure.vasp combine \
    --supercell 2 2 1 \
    --deduplicate-by-coverage \
    --coverage-bins 3
```

## 🔬 How It Works

### Algorithm Pipeline

```
Input Structure
    ↓
[1] Supercell Generation (optional)
    ↓
[2] Generate All Combinations
    ↓ (e.g., 6560 structures)
[3] Coverage Deduplication
    ↓ (→ 44 unique coverages, 99.3% reduction)
[4] Bin Sampling (optional)
    ↓ (→ 35 structures, 99.5% total reduction)
Output Representative Structures
```

### Coverage Calculation

Coverage is calculated as species per unit area:

```
coverage = n_species / surface_area_nm²
```

Where `surface_area_nm²` is the xy-plane projection calculated via cross product:

```python
area = |cell[0] × cell[1]| / 100.0  # Convert Ų to nm²
```

### Deduplication

Structures with identical coverage (rounded to precision) are considered duplicates. The first occurrence is kept.

**Example**:
- Structure A: 2 vacancies on 10.5 nm² = 0.1905 vac/nm²
- Structure B: 2 vacancies on 10.5 nm² = 0.1905 vac/nm² ← **duplicate**
- Structure C: 3 vacancies on 10.5 nm² = 0.2857 vac/nm²

Result: Keep A and C, discard B.

### Bin Sampling

**1D Binning** (vacancies, OH modes):
1. Divide coverage range [min, max] into N equal-width bins
2. Calculate bin centers: `center_i = min + (i + 0.5) × width`
3. For each bin, select structure with coverage closest to center

**2D Binning** (combined mode):
1. Create N×N grid in (vacancy_coverage, OH_coverage) space
2. For each grid cell, select structure closest to cell center
3. Uses Euclidean distance: `dist = √((Δvac)² + (ΔOH)²)`

**Example** (5 bins, coverage range 0-10):
- Bin 1: center = 1.0, selects structure with coverage 1.2
- Bin 2: center = 3.0, selects structure with coverage 2.9
- Bin 3: center = 5.0, selects structure with coverage 5.1
- Bin 4: center = 7.0, selects structure with coverage 7.0
- Bin 5: center = 9.0, selects structure with coverage 8.8

Result: 5 evenly-distributed representative structures.

## 📊 Results

### Test Case: st2.vasp with 2×2×1 Supercell

#### Without Optimization
```
Vacancies:  255 structures
Hydrogen:   255 structures
Combined:  6050 structures
Total:     6560 structures
```

#### With Coverage Deduplication Only
```
Vacancies:  255 → 8 structures (96.9% reduction)
Hydrogen:   255 → 8 structures (96.9% reduction)
Combined:  6050 → 28 structures (99.5% reduction)
Total:     6560 → 44 structures (99.3% reduction)
```

#### With Deduplication + 5-Bin Sampling
```
Vacancies:  255 → 8 → 5 structures
Hydrogen:   255 → 8 → 5 structures
Combined:  6050 → 28 → 25 structures (5×5 grid)
Total:     6560 → 44 → 35 structures (99.5% reduction)
```

**Performance**: ~1.8 seconds for complete mode with binning

### Coverage Distribution

Binning ensures even sampling across the coverage range:

```
Vacancies (5 bins):
  0.95, 1.43, 1.90, 2.38, 2.85 vac/nm²

Hydroxylation (5 bins):
  0.95, 1.43, 1.90, 2.38, 2.85 OH/nm²

Combined (3×3 grid):
  9 structures covering the (vac, OH) space evenly
```

## 📖 Command-Line Reference

### Required Arguments
```
input                 Input structure file (VASP, CIF, XYZ, etc.)
--mode MODE           vacancies|hydrogen|combine|complete
```

### Mode Selection
```
--which-surface       top|bottom|both (default: top)
--species SPECIES     Target species (default: O)
```

### Optimization Arguments
```
--supercell NX NY NZ           Create NX×NY×NZ supercell
--deduplicate-by-coverage      Enable coverage deduplication
--coverage-precision N         Decimal places for rounding (default: 4)
--coverage-bins N              Number of bins for sampling
```

### Output Options
```
--outdir DIR          Output directory (default: surf_variants)
--fmt FORMAT          Output format: vasp|cif|xyz (default: vasp)
--prefix PREFIX       Filename prefix (default: surf)
--tag-original        Include copy of original structure
--include-empty       Include no-modification baseline
--no-manifest         Skip JSON manifest generation
```

### Surface Parameters
```
--z-window ANGSTROM   Surface detection window (default: 0.5 Å)
--oh-dist ANGSTROM    O-H bond length (default: 0.98 Å)
```

## 📁 Manifest Format

Each mode generates a `manifest.json` with complete metadata:

```json
{
  "mode": "vacancies",
  "surface_area_nm2": 10.5,
  "deduplication_stats": {
    "total_generated": 255,
    "unique_coverages": 8,
    "kept": 8
  },
  "binning_stats": {
    "bins_requested": 5,
    "before_binning": 8,
    "after_binning": 5
  },
  "variants": [
    {
      "name": "surf_top_vac_n1_54",
      "file": "output/surf_top_vac_n1_54.vasp",
      "n_removed": 1,
      "removed_indices": [54],
      "which_surface": "top",
      "vacancy_coverage": 0.9524,
      "bin_id": 0,
      "bin_center": 0.9524
    },
    ...
  ]
}
```

### Combined Mode Manifest

For combined mode, coverage is a tuple `(vac_coverage, OH_coverage)`:

```json
{
  "mode": "combine",
  "binning_stats": {
    "bins_requested": 3,
    "binning_type": "2D_grid"
  },
  "variants": [
    {
      "name": "surf_top_combo_v1_54__H1_58",
      "vacancy_coverage": 0.9524,
      "OH_coverage": 0.9524,
      "bin_id_vac": 0,
      "bin_id_oh": 0,
      "bin_center_vac": 0.9524,
      "bin_center_oh": 0.9524
    },
    ...
  ]
}
```

## 🧪 Testing

### Unit Tests

```bash
# Run all unit tests (29 tests)
pytest tests/test_surface_modes.py -v

# Run specific test
pytest tests/test_surface_modes.py::test_sample_by_coverage_bins_1d -v
```

**Test Coverage**:
- Supercell generation
- Surface area calculation (cubic, rectangular, non-orthogonal)
- Coverage calculation and precision
- Deduplication (single and dual coverage)
- 1D and 2D binning algorithms
- Integration with run_vacancies(), run_hydrogen(), run_combine()

### Integration Tests

```bash
# Test deduplication with st2.vasp
./tests/test_integration_st2.sh

# Test binning with st2.vasp
./tests/test_integration_bins.sh
```

### Test Structure

```
tests/
├── test_surface_modes.py        # Unit tests (29 tests)
├── test_integration_st2.sh      # Deduplication integration
├── test_integration_bins.sh     # Binning integration
├── outputs/                     # Test output storage
└── [test output directories]    # Generated during testing
```

## 🏗️ Architecture

### Core Functions

```python
# Supercell and coverage utilities
make_supercell(atoms, nx, ny, nz) -> Atoms
calculate_surface_area(atoms) -> float  # nm²
calculate_coverage(n_species, area_nm2, precision) -> float

# Deduplication and binning
deduplicate_by_coverage(variants_list, mode) -> (list, stats)
sample_by_coverage_bins(variants_list, n_bins, mode) -> list
add_bin_metadata(variants_list, n_bins, mode) -> list
```

### SurfaceModifier Class

```python
class SurfaceModifier:
    def __init__(self, atoms, species='O', z_window=0.5,
                 which_surface='top', oh_dist=0.98,
                 supercell=None,
                 deduplicate_by_coverage=False,
                 coverage_precision=4,
                 coverage_bins=None, ...):
        ...

    def run_vacancies() -> manifest
    def run_hydrogen() -> manifest
    def run_combine() -> manifest
    def run_complete() -> manifest
```

### Processing Pipeline

1. **Surface Selection**: Identify surface atoms by z-coordinate
2. **Combination Generation**: Generate all non-empty subsets
3. **Modification Application**: Apply vacancies/OH groups
4. **Coverage Calculation**: Calculate species/nm² for each variant
5. **Deduplication**: Keep one per unique coverage
6. **Binning**: Sample N representatives
7. **Metadata Addition**: Add bin_id and bin_center
8. **File Writing**: Write structures and manifest

## 🎓 Advanced Usage

### Custom Coverage Precision

```bash
# Higher precision (5 decimals) - less aggressive deduplication
python surface_modes.py input.vasp complete --supercell 2 2 1 \
    --deduplicate-by-coverage --coverage-precision 5

# Lower precision (2 decimals) - more aggressive deduplication
python surface_modes.py input.vasp complete --supercell 2 2 1 \
    --deduplicate-by-coverage --coverage-precision 2
```

### Both Surfaces

```bash
# Operate on top and bottom surfaces independently
python surface_modes.py input.vasp complete --which-surface both \
    --supercell 2 2 1 --deduplicate-by-coverage --coverage-bins 5
```

### Custom Output Format

```bash
# Export as CIF files
python surface_modes.py input.vasp vacancies --fmt cif \
    --supercell 2 2 1 --deduplicate-by-coverage
```

## 📚 Additional Documentation

- `SUMMARY.txt` - Coverage deduplication feature summary
- `BINNING_RESULTS.txt` - Bin sampling results and performance
- `docs/plans/` - Design documents and implementation plans
- `tests/` - Complete test suite with examples

## 🤝 Contributing

The implementation follows strict TDD (Test-Driven Development):
1. Write failing test
2. Implement minimal code to pass test
3. Refactor while keeping tests green
4. Commit

All 29 unit tests must pass before merging changes.

## 📄 License

Part of the fosfato computational chemistry toolkit.

## 📮 Support

For issues or questions, please refer to the test files in `tests/` for usage examples, or consult the planning documents in `docs/plans/`.

---

**Generated with** [Claude Code](https://claude.com/claude-code)
