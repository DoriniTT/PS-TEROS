# Quick Start: Relaxation Energy Calculation

## What is the Relaxation Energy?

The relaxation energy measures how much energy is gained when atoms in a surface slab relax from their bulk-truncated positions to their equilibrium surface geometry:

```
E_relax = E_relaxed - E_unrelaxed
```

Negative values (typical) indicate that relaxation stabilizes the surface.

## How to Enable

Simply set `relax_slabs=True` when building your workgraph. The relaxation energy calculation is **automatically included** – no additional parameters needed!

```python
from teros.core.workgraph import build_core_workgraph
from aiida_workgraph import submit
from aiida import load_profile, orm

load_profile()

wg = build_core_workgraph(
    structures_dir="/path/to/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    bulk_parameters={...},
    bulk_options={...},
    slab_parameters={...},
    slab_options={...},
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,  # ← This enables relaxation energy calculation!
)

result, node, process = submit(wg)
print(f"WorkGraph PK: {node.pk}")
```

## What Happens Behind the Scenes

When `relax_slabs=True`, the workflow automatically:

1. **Performs SCF calculation** on each unrelaxed slab (NSW=0, IBRION=-1)
2. **Performs relaxation** on each slab (user-defined NSW and IBRION)
3. **Calculates relaxation energy** = E_relaxed - E_unrelaxed

All three steps run in parallel for all slab terminations using the scatter-gather pattern.

## Accessing Results

After the workgraph completes:

```python
from aiida import orm

# Load the completed workgraph
node = orm.load_node(PK)  # Replace PK with your workgraph PK

# Access relaxation energies
relaxation_energies = node.outputs.relaxation_energies

# Print results
print("\nRelaxation Energies:")
print("-" * 40)
for label, energy in relaxation_energies.items():
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

## All Available Outputs

With `relax_slabs=True`, you get access to:

| Output | Description |
|--------|-------------|
| `slab_structures` | Unrelaxed slab structures |
| `unrelaxed_slab_energies` | Energies from SCF (NSW=0) |
| `unrelaxed_slab_remote` | RemoteData for SCF calculations |
| `relaxed_slabs` | Relaxed slab structures |
| `slab_energies` | Energies from relaxation |
| `slab_remote` | RemoteData for relaxation |
| **`relaxation_energies`** | **E_relaxed - E_unrelaxed** ⭐ |

## Example: Complete Workflow

```python
#!/usr/bin/env python
"""Example: Ag2O (100) surface with relaxation energy calculation."""

from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph
from aiida_workgraph import submit

load_profile()

# Parameters
structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"

bulk_params = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

slab_params = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,  # Relax atoms, fix cell
    "NSW": 50,
    "EDIFFG": -0.05,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

options = {
    "resources": {
        "num_machines": 1,
        "num_cores_per_machine": 40,
    },
    "queue_name": "par40",
}

# Build workgraph
wg = build_core_workgraph(
    structures_dir=structures_dir,
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    bulk_parameters=bulk_params,
    bulk_options=options,
    metal_parameters=bulk_params,
    metal_options=options,
    oxygen_parameters=bulk_params,
    oxygen_options=options,
    slab_parameters=slab_params,
    slab_options=options,
    slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    kpoints_spacing=0.3,
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,  # Enable relaxation energy
    name='Ag2O_100_RelaxEnergy',
)

# Submit
result, node, process = submit(wg)
print(f"Submitted! WorkGraph PK: {node.pk}")
print(f"Check status: verdi process show {node.pk}")

# After completion (run separately):
# node = orm.load_node(PK)
# for label, energy in node.outputs.relaxation_energies.items():
#     print(f"{label}: {energy.value:+.4f} eV")
```

## Monitoring Progress

```bash
# Check status
verdi process show <PK>

# Watch progress
watch -n 5 'verdi process show <PK>'

# Check daemon
verdi daemon status
```

## Interpreting Results

### Typical Values
- Small slabs (< 10 Å): -0.5 to -2.0 eV
- Medium slabs (10-15 Å): -1.0 to -3.0 eV
- Large slabs (> 15 Å): -2.0 to -5.0 eV

### Warning Signs
- **Positive values**: Check convergence (NSW, EDIFFG)
- **Very large negative values (< -5 eV)**: May indicate:
  - Poor initial structure
  - Insufficient slab thickness
  - Convergence issues

### Comparison
Different terminations can be compared:
- Lower (more negative) = more atomic rearrangement
- Higher (less negative) = structure closer to bulk-like

## Works With All Features

The relaxation energy calculation is compatible with:

✅ **Slab Generation**: Automatically generates and calculates
✅ **Input Slabs**: Works with pre-provided slabs
✅ **Restart Mode**: Can restart from previous calculations
✅ **Surface Energies**: Use with `compute_thermodynamics=True`
✅ **Cleavage Energies**: Use with `compute_cleavage=True`
✅ **Binary Oxides**: Ag2O, CuO, etc.
✅ **Ternary Oxides**: Ag3PO4, Fe2WO6, etc.

## Computational Cost

- **SCF calculation**: ~5-10% of total slab calculation time
- **Relaxation**: Main computational cost (unchanged)
- **Energy calculation**: Negligible (<1 second)

Total overhead: **~5-10% increase** in wall time, but all calculations run in parallel.

## Need Help?

- Full documentation: `docs/RELAXATION_ENERGY.md`
- Implementation details: `RELAXATION_ENERGY_IMPLEMENTATION.md`
- Example script: `examples/slabs/test_relaxation_energy.py`

## Summary

To calculate relaxation energies:
1. Set `relax_slabs=True` (that's it!)
2. Submit your workgraph
3. Access `node.outputs.relaxation_energies` after completion

The feature is automatic, parallel, and fully integrated with existing PS-TEROS functionality.
