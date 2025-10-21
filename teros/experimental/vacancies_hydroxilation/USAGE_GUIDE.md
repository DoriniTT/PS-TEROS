# How to Run surface_modes.py

A step-by-step guide to using the surface modification tool with practical examples.

## Prerequisites

### Environment Setup

```bash
# Ensure you have the required packages
pip install ase numpy

# Activate your environment (if using conda/venv)
conda activate aiida  # or your environment name
```

### Input Structure File

You need a structure file in ASE-readable format:
- **VASP**: POSCAR, CONTCAR
- **CIF**: crystallographic information file
- **XYZ**: extended xyz format
- Other formats: ASE supports many more (see `ase.io` documentation)

Example: `st2.vasp` (included in repository)

---

## Quick Start

### 1. Basic Command - Generate All Variants

```bash
python surface_modes.py st2.vasp --mode complete
```

**What it does:**
- Generates vacancies, hydrogen (OH), and combined modifications
- Creates structures without optimization
- Outputs to `all_possibilities/` directory

**Output:**
```
Mode: complete. Wrote 6560 variants total to 'all_possibilities/'.
  - Vacancies: 255 variants
  - Hydrogen: 255 variants
  - Combine: 6050 variants
Manifest: all_possibilities/manifest.json
```

### 2. With Coverage Deduplication (96% reduction)

```bash
python surface_modes.py st2.vasp --mode complete \
  --supercell 2 2 1 \
  --deduplicate-by-coverage
```

**What it does:**
- Creates 2×2×1 supercell
- Keeps only one structure per unique surface coverage
- Reduces 6560 → 44 structures

**Output:**
```
Mode: complete. Wrote 44 variants total to 'all_possibilities/'.
  - Vacancies: 255 -> 8 variants
  - Hydrogen: 255 -> 8 variants
  - Combine: 6050 -> 28 variants
  - Vacancies deduplication: 255 -> 8 variants
  - Hydrogen deduplication: 255 -> 8 variants
  - Combine deduplication: 6050 -> 28 variants
```

### 3. With Deduplication + Binning (99.5% reduction)

```bash
python surface_modes.py st2.vasp --mode complete \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-bins 5
```

**What it does:**
- Deduplicates by coverage (44 structures)
- Samples into 5 bins (1D) or 5×5 grid (2D for combined)
- Reduces 6560 → 35 structures total

**Output:**
```
Mode: complete. Wrote 35 variants total to 'all_possibilities/'.
  - Vacancies: 255 -> 8 -> 5 variants
  - Hydrogen: 255 -> 8 -> 5 variants
  - Combine: 6050 -> 28 -> 25 variants
  - Vacancies deduplication: 255 -> 8 variants
  - Vacancies binning: 8 -> 5 (5 bins)
  - Hydrogen deduplication: 255 -> 8 variants
  - Hydrogen binning: 8 -> 5 (5 bins)
  - Combine deduplication: 6050 -> 28 variants
  - Combine binning: 28 -> 25 (5 bins)
```

---

## Mode-Specific Examples

### Vacancies Only

```bash
# Generate all vacancy combinations on top surface
python surface_modes.py st2.vasp --mode vacancies --which-surface top

# With deduplication and custom bins
python surface_modes.py st2.vasp --mode vacancies \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-bins 10 \
  --outdir output/vacancies
```

### Hydrogen (OH) Only

```bash
# Generate all hydroxylation combinations on bottom surface
python surface_modes.py st2.vasp --mode hydrogen \
  --which-surface bottom \
  --oh-dist 1.00

# With deduplication
python surface_modes.py st2.vasp --mode hydrogen \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --outdir output/hydroxyl
```

### Combined (Vacancies + OH)

```bash
# Generate all combined modifications on both surfaces
python surface_modes.py st2.vasp --mode combine \
  --which-surface both \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-bins 3 \
  --outdir output/combined
```

In combined mode with `--coverage-bins 3`, you get a 3×3 grid (up to 9 structures).

---

## Command-Line Options Reference

### Essential Arguments

```bash
python surface_modes.py INPUT_FILE --mode MODE [options]
```

| Argument | Values | Default | Description |
|----------|--------|---------|-------------|
| `INPUT_FILE` | path/to/structure | — | Input structure (POSCAR, CIF, XYZ, etc.) |
| `--mode` | vacancies, hydrogen, combine, complete | — | **REQUIRED** - which modifications to generate |

### Surface Selection

```bash
--which-surface {top, bottom, both}  # Default: top
  # top    - Modify only top surface
  # bottom - Modify only bottom surface
  # both   - Modify both surfaces independently
```

### Structure Modification

```bash
--species ATOM                      # Default: O
  # Target atom species (typically oxygen)

--z-window ANGSTROM                 # Default: 0.5
  # Distance window in Å from surface extreme z-coordinate
  # Larger values include more atoms as "surface"

--oh-dist ANGSTROM                  # Default: 0.98
  # O-H bond length in Ångströms for hydrogen mode
  # Typical values: 0.90-1.10 Å
```

### Supercell and Optimization

```bash
--supercell NX NY NZ                # Default: none (1×1×1)
  # Create supercell before modifications
  # Example: --supercell 2 2 1

--deduplicate-by-coverage           # Default: OFF
  # Keep only one structure per unique coverage (species/nm²)
  # Reduces by 90-99%

--coverage-precision N              # Default: 4
  # Decimal places for coverage rounding
  # Higher = stricter deduplication (fewer structures)
  # Lower = more relaxed (more structures kept)

--coverage-bins N                   # Default: none (no binning)
  # Sample N representative structures across coverage range
  # Requires --deduplicate-by-coverage
  # For combined mode: creates N×N grid (up to N² structures)
```

### Output Options

```bash
--outdir PATH                       # Default: surf_variants or all_possibilities
  # Output directory for structures and manifest

--fmt {vasp, cif, xyz}              # Default: vasp
  # Output structure format
  # VASP: most common for DFT codes
  # CIF: crystallographic format
  # XYZ: simple atomic format

--prefix NAME                       # Default: surf
  # Prefix for output filenames
  # Example: --prefix my_struct → my_struct_top_vac_n1_54.vasp

--tag-original                      # Default: OFF
  # Write original structure as 000_original.*
  # Useful for reference

--include-empty                     # Default: OFF
  # Include baseline (no modification) variant
  # Useful to have unmodified structure in results

--no-manifest                       # Default: OFF
  # Skip JSON manifest generation
  # Usually keep this OFF for tracking metadata
```

---

## Practical Examples

### Example 1: Quick Test

```bash
python surface_modes.py st2.vasp --mode vacancies --which-surface top
```

**Result:** 255 VASP files in `surf_variants/`

**Use case:** Quickly test if script runs with your structure

---

### Example 2: Optimized for Computation

```bash
python surface_modes.py st2.vasp --mode complete \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-bins 5 \
  --fmt vasp \
  --tag-original \
  --outdir my_surfaces
```

**Result:** 35 optimized structures in `my_surfaces/`

**Use case:** Reduce computational cost while maintaining coverage diversity

---

### Example 3: High-Precision Deduplication

```bash
python surface_modes.py st2.vasp --mode combine \
  --supercell 3 3 1 \
  --deduplicate-by-coverage \
  --coverage-precision 3 \
  --outdir high_precision
```

**Result:** Many unique structures (3 decimal places = strict)

**Use case:** Explore fine details of coverage landscape

---

### Example 4: Low-Precision Aggressive Reduction

```bash
python surface_modes.py st2.vasp --mode combine \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-precision 1 \
  --coverage-bins 3 \
  --outdir aggressive
```

**Result:** Very few structures (1 decimal place = relaxed, 3×3 grid)

**Use case:** Minimal set for quick screening

---

### Example 5: Both Surfaces with Binning

```bash
python surface_modes.py st2.vasp --mode combine \
  --which-surface both \
  --supercell 2 2 1 \
  --deduplicate-by-coverage \
  --coverage-bins 4 \
  --outdir both_surfaces
```

**Result:** Up to 64 structures (4×4 grid for each surface pair)

**Use case:** Study surface-surface interactions

---

## Understanding the Output

### Manifest File Format

Each run generates `manifest.json` with complete metadata:

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
      "deleted_indices": [54],
      "surface": "top",
      "count_deleted": 1,
      "surface_area_nm2": 10.5,
      "vacancy_coverage": 0.0952,
      "bin_id": 0,
      "bin_center": 0.0952
    }
    // ... more variants ...
  ]
}
```

**Key fields:**
- `name`: Structure identifier
- `file`: Path to VASP/CIF/XYZ file
- `*_coverage`: Surface coverage in species/nm²
- `bin_id`: Which bin (0 to n_bins-1)
- `bin_center`: Target coverage for that bin

---

## Tips and Tricks

### 1. Estimate Output Size Before Running

```bash
# Dry run with just first few structures
python surface_modes.py st2.vasp --mode combine --no-manifest
```

### 2. Test with Small Supercell First

```bash
# Test with 1×1×1 (original size)
python surface_modes.py st2.vasp --mode vacancies --outdir test_1x1x1

# Then scale up
python surface_modes.py st2.vasp --mode vacancies --supercell 2 2 1
```

### 3. Use Reasonable Bin Counts

- **1D modes (vacancies, hydrogen):** 5-10 bins typical
- **2D mode (combined):** 3-5 bins typical (creates 9-25 structures)
- **More bins** = more structures = more computation

```bash
# Conservative (few structures)
--coverage-bins 3

# Balanced (medium structures)
--coverage-bins 5

# Comprehensive (many structures)
--coverage-bins 10
```

### 4. Check Structure Validity After Generation

```bash
# Validate one generated structure
python -c "from ase.io import read; atoms = read('output/structure.vasp'); print(f'Atoms: {len(atoms)}, Formula: {atoms.get_chemical_formula()}')"
```

### 5. Monitor Coverage Distribution

```bash
# Parse manifest to check coverage values
python << 'EOF'
import json
with open('output/manifest.json') as f:
    data = json.load(f)
for v in data['variants']:
    print(f"{v['name']}: coverage={v.get('vacancy_coverage', v.get('OH_coverage'))}")
EOF
```

---

## Common Issues and Solutions

### Issue: "No atoms of species 'O' found"

**Problem:** Script can't find oxygen atoms

**Solution:**
```bash
# Check what atoms are in structure
python -c "from ase.io import read; atoms = read('st2.vasp'); print([a.symbol for a in atoms])"

# Use correct species
python surface_modes.py st2.vasp --mode vacancies --species Ti
```

### Issue: "No 'O' atoms found within z_window"

**Problem:** Surface atoms are outside z_window

**Solution:**
```bash
# Increase z_window
python surface_modes.py st2.vasp --mode vacancies --z-window 1.0

# Or check structure orientation
python -c "from ase.io import read; atoms = read('st2.vasp'); print('Z coords:', atoms.positions[:, 2])"
```

### Issue: Too Many/Too Few Structures

**Too many:**
```bash
# More aggressive deduplication
--coverage-precision 2  # Fewer decimal places
--coverage-bins 3       # Fewer bins
```

**Too few:**
```bash
# Less aggressive deduplication
--coverage-precision 5  # More decimal places
--coverage-bins 10      # More bins
```

### Issue: Out of Memory

**Problem:** Too many structures in memory at once

**Solution:** Can't currently fix (single mode only), but output directory shows progress

---

## Performance Guidelines

| Supercell | Mode | Structures | Time* | Memory |
|-----------|------|-----------|-------|--------|
| 1×1×1 | vacancies | ~15 | <0.1s | ~10MB |
| 1×1×1 | combined | ~100 | <1s | ~50MB |
| 2×2×1 | combined | 6560 | ~2s | ~500MB |
| 2×2×1 | combined + dedup | 44 | ~2s | ~500MB |
| 2×2×1 | combined + dedup + bins | 35 | ~2s | ~500MB |

*Approximate on typical laptop

---

## Next Steps After Generation

1. **Review manifest.json** to understand generated structures

2. **Pick representative structures** for DFT calculations based on coverage bins

3. **Convert format if needed:**
   ```bash
   python -c "from ase.io import read, write; write('output.cif', read('input.vasp'))"
   ```

4. **Run DFT calculations** on selected structures

5. **Analyze results** using coverage metadata from manifest

---

## Questions?

Refer to:
- `README.md` - Comprehensive feature documentation
- `SUMMARY.txt` - Deduplication results overview
- `BINNING_RESULTS.txt` - Binning effectiveness data
- `tests/test_integration_*.sh` - Example workflows

