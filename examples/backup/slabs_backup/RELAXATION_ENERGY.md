# Relaxation Energy Calculation in PS-TEROS

## Overview

The relaxation energy module calculates the energy difference between relaxed and unrelaxed slab structures. This is an important quantity for understanding the energetic stabilization due to atomic relaxation at surfaces.

## Theory

The relaxation energy is defined as:

```
E_relax = E_relaxed - E_unrelaxed
```

Where:
- `E_unrelaxed`: Total energy of the slab structure in its initial (unrelaxed) configuration, obtained from a single-point SCF calculation (NSW=0, IBRION=-1)
- `E_relaxed`: Total energy of the slab structure after atomic relaxation
- `E_relax`: Relaxation energy (negative values indicate stabilization)

A negative relaxation energy indicates that the relaxation process stabilizes the system, which is typically the case for surfaces.

## Implementation

### Workflow Steps

When `relax_slabs=True` is enabled in the workflow, the following steps are executed for each slab termination:

1. **SCF Calculation (Unrelaxed)**: Performs a single-point energy calculation on the unrelaxed slab structure
   - VASP parameters: `NSW=0`, `IBRION=-1`
   - Output: `unrelaxed_slab_energies` (total energy)

2. **Relaxation Calculation**: Performs atomic relaxation on the slab structure
   - VASP parameters: User-defined (typically `IBRION=2`, `NSW>0`)
   - Output: `relaxed_slabs` (structure) and `slab_energies` (total energy)

3. **Relaxation Energy Calculation**: Computes the energy difference
   - Calculation: `E_relax = E_relaxed - E_unrelaxed`
   - Output: `relaxation_energies` (energy difference in eV)

### Key Functions

#### `scf_slabs_scatter`
Located in `teros/core/slabs.py`

Performs SCF calculations on unrelaxed slab structures in parallel using the scatter-gather pattern.

**Arguments:**
- `slabs`: Dictionary of unrelaxed slab structures
- `code`: AiiDA code for VASP
- `potential_family`, `potential_mapping`: Pseudopotential settings
- `parameters`: VASP parameters (NSW and IBRION are overridden)
- `options`: Computation resources
- `kpoints_spacing`: K-points spacing
- `clean_workdir`: Whether to clean remote directories

**Returns:**
- `energies`: Dictionary of total energies for each slab
- `remote_folders`: Dictionary of RemoteData nodes

#### `calculate_relaxation_energies_scatter`
Located in `teros/core/slabs.py`

Calculates relaxation energies for all slab terminations using the scatter-gather pattern.

**Arguments:**
- `unrelaxed_energies`: Dictionary of energies from SCF calculations
- `relaxed_energies`: Dictionary of energies from relaxation calculations

**Returns:**
- `relaxation_energies`: Dictionary of relaxation energies (E_relaxed - E_unrelaxed)

### Integration in WorkGraph

The relaxation energy calculation is automatically integrated into the `core_workgraph` when `relax_slabs=True`:

```python
from teros.core.workgraph import build_core_workgraph

wg = build_core_workgraph(
    structures_dir="/path/to/structures",
    bulk_name="oxide.cif",
    metal_name="metal.cif",
    oxygen_name="O2.cif",
    # ... other parameters ...
    relax_slabs=True,  # Enable relaxation and relaxation energy calculation
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
)
```

## Accessing Results

After the workflow completes, you can access the relaxation energies:

```python
from aiida import orm

# Load the workgraph node
node = orm.load_node(PK)

# Access relaxation energies
relaxation_energies = node.outputs.relaxation_energies

# Print results
for label, energy in relaxation_energies.items():
    print(f"{label}: {energy.value:.4f} eV")
```

Example output:
```
term_0: -2.3456 eV
term_1: -1.8932 eV
term_2: -3.1245 eV
```

## Available Outputs

When `relax_slabs=True`, the following outputs are available:

| Output | Type | Description |
|--------|------|-------------|
| `slab_structures` | StructureData | Unrelaxed slab structures |
| `unrelaxed_slab_energies` | Float | Total energies from SCF calculations |
| `unrelaxed_slab_remote` | RemoteData | Remote folders for SCF calculations |
| `relaxed_slabs` | StructureData | Relaxed slab structures |
| `slab_energies` | Float | Total energies from relaxation |
| `slab_remote` | RemoteData | Remote folders for relaxation |
| `relaxation_energies` | Float | Relaxation energies (E_relaxed - E_unrelaxed) |

All outputs use dynamic namespaces, so each slab termination (e.g., `term_0`, `term_1`) has its own entry.

## Example Usage

### Full Example with Slab Generation

```python
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph
from aiida_workgraph import submit

load_profile()

# Define parameters
structures_dir = "/path/to/structures"
bulk_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    # ... other INCAR tags ...
}

slab_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "IBRION": 2,
    "ISIF": 2,  # Relax atoms, fix cell
    "NSW": 50,
    # ... other INCAR tags ...
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
    bulk_parameters=bulk_parameters,
    bulk_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}},
    metal_parameters=bulk_parameters,
    metal_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}},
    oxygen_parameters=bulk_parameters,
    oxygen_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}},
    slab_parameters=slab_parameters,
    slab_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}},
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,  # Enable relaxation energy calculation
    name='RelaxationEnergyCalculation',
)

# Submit
result, node, process = submit(wg)
print(f"Submitted WorkGraph PK: {node.pk}")

# After completion, access results
# node = orm.load_node(PK)
# for label, energy in node.outputs.relaxation_energies.items():
#     print(f"{label}: {energy.value:.4f} eV")
```

### Using Pre-generated Slabs

```python
from aiida import orm

# Load pre-generated slabs
input_slabs = {
    'term_0': orm.load_node(12345),
    'term_1': orm.load_node(12346),
}

wg = build_core_workgraph(
    # ... other parameters ...
    input_slabs=input_slabs,  # Use pre-generated slabs
    relax_slabs=True,
)
```

## VASP Parameters

### SCF Calculation (Unrelaxed)
The SCF calculation automatically sets:
- `NSW = 0`: No ionic steps
- `IBRION = -1`: No relaxation algorithm

All other INCAR parameters are taken from `slab_parameters` (or `bulk_parameters` if not specified).

### Relaxation Calculation
Uses the parameters specified in `slab_parameters`. Typical settings:
- `IBRION = 2`: Conjugate gradient algorithm
- `ISIF = 2`: Relax atoms, keep cell fixed
- `NSW = 50-200`: Number of ionic steps
- `EDIFFG = -0.01 to -0.05`: Force convergence criterion

## Performance Notes

- Both SCF and relaxation calculations run in **parallel** for all slab terminations using the scatter-gather pattern
- The SCF calculation is typically much faster than relaxation (single point vs. many ionic steps)
- Total computational cost ≈ Cost(SCF) + Cost(Relaxation) for all slabs in parallel

## Physical Interpretation

The relaxation energy provides insights into:

1. **Surface Stability**: Larger negative values indicate more substantial atomic rearrangement and stabilization
2. **Termination Comparison**: Different terminations can be compared by their relaxation energies
3. **Convergence Verification**: Large relaxation energies may indicate need for:
   - More relaxation steps (higher NSW)
   - Tighter convergence criteria (lower EDIFFG)
   - Thicker slabs

## Troubleshooting

### Issue: SCF calculation fails
**Solution**: Check that the unrelaxed slab structure is physically reasonable (no overlapping atoms, appropriate vacuum, etc.)

### Issue: Relaxation energy is positive
**Cause**: This can happen if:
- Relaxation didn't converge (check NSW, EDIFFG)
- Initial structure was already optimized
- Numerical precision issues

**Solution**: Check convergence, increase NSW, or tighten EDIFFG

### Issue: Very large negative relaxation energy (< -5 eV)
**Cause**: May indicate:
- Poor initial structure
- Insufficient slab thickness
- Convergence issues

**Solution**: Review initial structure, increase slab thickness, check convergence

## References

- Relaxation energy is commonly reported in surface science literature
- Related to surface energy but focuses specifically on atomic rearrangement contribution
- Important for comparing DFT results with experimental surface preparation methods
