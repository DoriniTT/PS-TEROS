# Thickness Convergence Example: Au(111)

This example demonstrates how to use the thickness convergence module to determine the minimum slab thickness for converged surface energy calculations.

## Overview

The workflow tests Au(111) slabs at different thicknesses (3, 5, 7, 9, 11 layers) to find the minimum thickness where the surface energy converges.

**Expected Results:**
- Convergence: ~5-7 layers
- Surface energy: ~0.79 J/m²

## Prerequisites

1. **AiiDA environment:**
   ```bash
   source ~/envs/aiida/bin/activate
   verdi profile set-default presto
   ```

2. **Cluster configuration:** This example is configured for the **obelix** cluster with:
   - Code: `VASP-6.5.1-idefix@obelix`
   - Potential family: `PBE`
   - PBS scheduler with Skylake nodes
   - Hybrid MPI+OpenMP parallelization (4 MPI processes)

## Setup

### Step 1: Create Au.cif Structure File

Create the bulk Au structure file in this directory:

```python
from ase.build import bulk
from ase.io import write

# Create FCC Au with experimental lattice parameter
au = bulk('Au', 'fcc', a=4.08)

# Write to CIF file
write('Au.cif', au)
```

Or run this one-liner:
```bash
python -c "from ase.build import bulk; from ase.io import write; write('Au.cif', bulk('Au', 'fcc', a=4.08))"
```

### Step 2: Verify Structure

Check that `Au.cif` exists in this directory:
```bash
ls -lh Au.cif
```

## Usage

### Submit Workflow

```bash
python run_thickness_convergence.py
```

This will:
1. Load the Au.cif structure
2. Build a thickness convergence WorkGraph
3. Submit to the AiiDA daemon
4. Print the WorkGraph PK for monitoring

### Monitor Progress

```bash
# Check overall status
verdi process show <PK>

# View detailed task hierarchy
verdi process report <PK>

# List all sub-processes
verdi process list -a -p 1
```

### Get Results

After the workflow completes (status shows `Finished [0]`):

```bash
python run_thickness_convergence.py <PK>
```

This will print:
- Convergence status
- Recommended number of layers
- Surface energy for each thickness
- Comparison with experimental reference data

## Expected Output

```
========================================================================
THICKNESS CONVERGENCE RESULTS
========================================================================

Converged: True
Recommended layers: 7

Bulk energy: -3.756234 eV

--- Surface Energy vs. Thickness ---
Layers     gamma (J/m^2)   gamma (eV/Ang^2)
----------------------------------------
3               0.8245         0.051468
5               0.7968         0.049740
7               0.7896         0.049291  <-- RECOMMENDED
9               0.7903         0.049335
11              0.7899         0.049310

--- Convergence Analysis ---
Convergence threshold: 0.01 J/m^2
Max tested layers: 11
Miller indices: (1 1 1)

--- Reference Data (Au 111) ---
Experimental surface energy: ~0.79 J/m^2
Typical convergence: 5-7 layers
========================================================================
```

## Workflow Structure

The thickness convergence workflow performs these steps:

1. **Bulk Relaxation**
   - Full cell optimization (ISIF=3)
   - High-precision parameters (ENCUT=500, ISMEAR=1)

2. **Slab Generation**
   - Creates slabs at specified thicknesses
   - Same termination for all slabs
   - 20 Å vacuum thickness

3. **Slab Relaxation** (parallel)
   - Ionic relaxation only (ISIF=2, fixed cell)
   - Same parameters as bulk
   - 4 jobs run concurrently

4. **Surface Energy Calculation**
   - γ = (E_slab - N_atoms * E_bulk_per_atom) / (2 * A)
   - Converted to J/m²

5. **Convergence Analysis**
   - Checks if consecutive thicknesses differ by < threshold
   - Reports recommended thickness

## Customization

### Test Different Surfaces

Edit `run_thickness_convergence.py`:

```python
# For Au(100)
miller_indices = [1, 0, 0]

# For Au(110)
miller_indices = [1, 1, 0]
```

### Adjust Layer Counts

```python
# Finer sampling around expected convergence
layer_counts = [3, 5, 7, 9, 11, 13, 15]

# Quick test with fewer layers
layer_counts = [3, 5, 7]
```

### Change Convergence Threshold

```python
# Stricter convergence (more layers needed)
convergence_threshold = 0.005  # J/m^2

# Looser convergence (fewer layers)
convergence_threshold = 0.02  # J/m^2
```

## Troubleshooting

### Structure file not found
```
ERROR: Structure file not found!
```
**Solution:** Create `Au.cif` as described in Step 1 above.

### Code not found
```
NotExistent: Code 'VASP-6.5.1-idefix@obelix' does not exist
```
**Solution:** Check available codes with `verdi code list` and update `code_label`.

### Potential family not found
```
ValueError: Potential family 'PBE' not found
```
**Solution:** Check available families with `verdi data core.potcar listfamilies` and update `potential_family`.

### Workflow stuck in Waiting
- Check daemon status: `verdi daemon status`
- Check PBS queue: `qstat -u $USER` (if you have SSH access)
- View detailed report: `verdi process report <PK>`

### VASP calculation failed
- Check VASP output: `verdi calcjob outputcat <CALC_PK>`
- Common issues:
  - POTCAR mismatch: Verify `potential_mapping={'Au': 'Au'}`
  - Resource limits: Check PBS queue limits on obelix

## Files

- `run_thickness_convergence.py` - Main script for submission and results
- `Au.cif` - Bulk Au structure (you create this)
- `README.md` - This file

## References

- Au(111) surface energy: Tyson, W. R., & Miller, W. A. (1977). Surface free energies of solid metals. *Surface Science*, 62(1), 267-276.
- PS-TEROS documentation: `/home/thiagotd/git/PS-TEROS/CLAUDE.md`
