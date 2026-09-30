# Relaxation Energy Example

## Overview

This example demonstrates how to calculate the relaxation energy for Ag2O (100) surface slabs.

## What is Relaxation Energy?

The relaxation energy measures the energetic stabilization when atoms relax from their bulk-truncated positions to their equilibrium surface geometry:

```
E_relax = E_relaxed - E_unrelaxed
```

A negative value (typical) indicates that relaxation stabilizes the surface.

## Example Script

**File**: `ag2o_100_relaxation_energy.py`

This complete working example:
1. Relaxes bulk Ag2O and reference structures (Ag, O2)
2. Generates slab structures for the (100) orientation
3. **Automatically** performs SCF on unrelaxed slabs (NSW=0, IBRION=-1)
4. Relaxes the slabs (user-specified parameters)
5. **Automatically** calculates relaxation energies
6. Outputs all results

## Running the Example

```bash
# Activate the psteros environment
source ~/envs/psteros/bin/activate

# Navigate to the examples directory
cd /home/thiagotd/git/worktree/PS-TEROS/feature-relax-energy/examples/slabs

# Run the script
python ag2o_100_relaxation_energy.py
```

## What Happens

### Step 1: Workflow Setup
The script configures all parameters and builds a WorkGraph with `relax_slabs=True`.

### Step 2: Submission
The WorkGraph is submitted to AiiDA. You'll see output like:
```
WORKGRAPH PK: 12345
```

### Step 3: Execution
The workflow runs these tasks in parallel:
- Bulk relaxation (Ag2O)
- Reference relaxations (Ag, O2)
- Slab generation from relaxed bulk
- **SCF calculations on unrelaxed slabs** (automatic)
- **Slab relaxations** (your parameters)
- **Relaxation energy calculation** (automatic)

### Step 4: Monitoring
While the workflow runs:
```bash
# Check status
verdi process show 12345

# Watch in real-time
watch -n 5 'verdi process show 12345'
```

### Step 5: Results
After completion (typically 1-4 hours depending on system size):

```python
from aiida import orm

# Load the workgraph
node = orm.load_node(12345)  # Replace with your PK

# Access relaxation energies
print("\nRelaxation Energies:")
print("-" * 40)
for label, energy in node.outputs.relaxation_energies.items():
    print(f"  {label}: {energy.value:+.4f} eV")
```

Example output:
```
Relaxation Energies:
----------------------------------------
  term_0: -2.3456 eV
  term_1: -1.8932 eV
  term_2: -3.1245 eV
```

## Available Outputs

After the workflow completes, you can access:

| Output | Description |
|--------|-------------|
| `slab_structures` | Unrelaxed slab structures (generated from bulk) |
| `unrelaxed_slab_energies` | Total energies from SCF (NSW=0, IBRION=-1) |
| `unrelaxed_slab_remote` | RemoteData for SCF calculations |
| `relaxed_slabs` | Relaxed slab structures |
| `slab_energies` | Total energies from relaxation |
| `slab_remote` | RemoteData for relaxation calculations |
| **`relaxation_energies`** | **E_relaxed - E_unrelaxed for each slab** |

## Accessing All Data

```python
from aiida import orm

node = orm.load_node(PK)

# Slab structures (unrelaxed)
for label, structure in node.outputs.slab_structures.items():
    print(f"{label}: {structure.pk} - {len(structure.sites)} atoms")

# Unrelaxed energies (from SCF)
for label, energy in node.outputs.unrelaxed_slab_energies.items():
    print(f"{label}: E_unrelaxed = {energy.value:.4f} eV")

# Relaxed energies
for label, energy in node.outputs.slab_energies.items():
    print(f"{label}: E_relaxed = {energy.value:.4f} eV")

# Relaxation energies (the difference)
for label, energy in node.outputs.relaxation_energies.items():
    print(f"{label}: E_relax = {energy.value:+.4f} eV")
```

## Customizing the Example

### Change Surface Orientation
```python
miller_indices = [1, 1, 0]  # (110) surface
```

### Adjust Slab Size
```python
min_slab_thickness = 15.0   # Thicker slab
min_vacuum_thickness = 20.0  # More vacuum
```

### Modify Relaxation Parameters
```python
slab_parameters = {
    # ... other parameters ...
    "NSW": 100,        # More relaxation steps
    "EDIFFG": -0.01,   # Tighter convergence
}
```

### Use Different Material
```python
bulk_name="cuo.cif",        # CuO instead of Ag2O
metal_name="Cu.cif",
# ... update potential mappings ...
```

## Key Points

1. **Automatic**: Just set `relax_slabs=True` - the SCF and relaxation energy calculation are automatic!

2. **Parallel**: All slabs are processed in parallel for maximum efficiency

3. **Complete Provenance**: All calculations are tracked in AiiDA

4. **Flexible**: Works with any oxide system (binary or ternary)

## Computational Cost

For Ag2O (100) with typical settings:
- Bulk + references: ~30 minutes
- Slab generation: <1 minute
- SCF per slab: ~5-10 minutes
- Relaxation per slab: ~30-60 minutes

**Total**: ~1-2 hours for 2-3 terminations (running in parallel)

## Troubleshooting

### Issue: Workflow fails during slab generation
**Check**: 
- Bulk relaxation completed successfully
- Miller indices are valid for the structure
- Slab thickness parameters are reasonable

### Issue: SCF calculation fails
**Check**:
- Unrelaxed slab structure is physically reasonable
- No overlapping atoms
- Sufficient vacuum

### Issue: Very large relaxation energy (< -5 eV)
**Possible causes**:
- Poor initial structure
- Insufficient slab thickness
- Convergence issues

**Solution**: Review structure, increase thickness, tighten convergence

## Next Steps

1. Run the example and monitor progress
2. After completion, analyze the relaxation energies
3. Compare different terminations
4. Try different surface orientations
5. Apply to your own material systems

## Additional Resources

- User Guide: `../../docs/RELAXATION_ENERGY.md`
- Implementation Details: `../../RELAXATION_ENERGY_IMPLEMENTATION.md`
- Quick Start: `../../QUICKSTART_RELAXATION_ENERGY.md`

## Questions?

The relaxation energy calculation is automatically integrated into PS-TEROS. Just use `relax_slabs=True` and you'll get:
- SCF energies (automatic)
- Relaxed energies (your parameters)
- Relaxation energies (automatic)

No additional configuration needed!
