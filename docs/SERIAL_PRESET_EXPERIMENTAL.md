# Serial Surface Thermodynamics Preset (EXPERIMENTAL)

**Status:** ⚠️ EXPERIMENTAL - In Active Development

**Module:** `teros.experimental.surface_thermo_preset_serial`

**Available since:** Development version (not yet released)

**Example:** `examples/vasp/step_16_surface_thermodynamics_serial.py`

---

## Overview

The serial preset creates a flat-graph architecture for surface thermodynamics workflows. All VASP calculation nodes exist at the same graph level, allowing `max_number_jobs` to control concurrent execution across all calculations.

**This solves a limitation in the standard preset where `max_number_jobs` does not propagate to nested sub-workgraphs.**

---

## Why This Exists

### Problem: Nested Sub-Workgraphs

The standard surface thermodynamics preset uses `@task.graph` decorators to create nested sub-workgraphs:

```python
@task.graph
def relax_slabs(slabs, ...):
    """Creates isolated graph context."""
    for slab in slabs:
        vasp_node = VaspWorkChain(...)  # In nested context
    return results
```

**Issue:** Each `@task.graph` creates an isolated graph context. Setting `max_number_jobs` on the main graph does not propagate to nested graphs.

**Result:** All VASP jobs submit simultaneously regardless of the limit.

### Solution: Flat Graph Architecture

The serial preset replaces nested sub-workgraphs with direct node addition:

```python
def build_relax_slabs_nodes(wg, slabs, ...):
    """Adds nodes to main graph."""
    for slab in slabs:
        wg.add_task(VaspWorkChain, ...)  # Same level as main graph
    return nodes
```

**Result:** All VASP nodes exist in the main graph. Setting `wg.max_number_jobs=2` limits concurrent execution to 2 VASP jobs.

---

## Current Status

✅ **Implemented:**
- Core flat-graph architecture
- All five workflow phases (bulk, references, slabs, calculations, thermodynamics)
- Concurrency control with `max_number_jobs`
- Basic testing and verification

⚠️ **Limitations:**
- Pre-generated slabs required (dynamic generation not supported)
- Limited production testing
- API may change
- Experimental module location

❌ **Not Yet Implemented:**
- Dynamic slab generation from bulk structure
- Two-stage execution with bulk-then-slabs
- Comprehensive error handling

---

## Installation and Import

The serial preset is in the experimental module:

```python
from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph
)
```

**Note:** Do not use `from teros.experimental import *`. Users must explicitly import experimental modules.

---

## Basic Usage

### 1. Prepare Pre-Generated Slabs

The serial preset requires slabs as input:

```python
from aiida import orm
from ase.io import read

slabs_dict = {
    'slab_100_term_0': orm.StructureData(ase=read('slab_100_term_0.cif')),
    'slab_100_term_1': orm.StructureData(ase=read('slab_100_term_1.cif')),
    'slab_110': orm.StructureData(ase=read('slab_110.cif')),
}
```

### 2. Build Workgraph

```python
from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph
)

wg = surface_thermodynamics_serial_workgraph(
    # Structures
    structures_dir='structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',

    # Code and parameters
    code_label='VASP-6.5.0@bohr-new',
    potential_family='PBE',
    kpoints_spacing=0.4,

    # Bulk calculations
    bulk_parameters={'PREC': 'Accurate', 'ENCUT': 420, ...},
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    bulk_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 40}},

    # Reference calculations
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 420, ...},
    metal_potential_mapping={'Ag': 'Ag'},
    metal_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 40}},

    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 420, ...},
    oxygen_potential_mapping={'O': 'O'},
    oxygen_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 40}},

    # Slab calculations (REQUIRED)
    input_slabs=slabs_dict,
    slab_parameters={'PREC': 'Accurate', 'ISIF': 2, ...},
    slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    slab_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 40}},

    # Control flags
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_relaxation_energy=True,
)
```

### 3. Set Concurrency Limit

```python
# Limit to 2 concurrent VASP jobs
wg.max_number_jobs = 2
```

### 4. Submit

```python
result = wg.submit()
pk = result.pk
```

---

## Parameter Differences from Standard Preset

| Parameter | Standard Preset | Serial Preset |
|-----------|----------------|---------------|
| Slab generation | `miller_indices=(1,0,0)` | Not supported |
| Slab input | Optional | **Required** (`input_slabs`) |
| `@task.graph` | Yes (nested graphs) | No (flat graph) |
| `max_number_jobs` | Limited effect | Controls all VASP |

---

## VASP Parameter Format

**Important:** VASP parameters must be wrapped with `incar` key:

```python
# Correct
bulk_parameters = {
    'PREC': 'Accurate',
    'ENCUT': 420,
    'EDIFF': 1e-4,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 3,
    'NSW': 100,
}

# The preset internally wraps this as: {'incar': bulk_parameters}
```

This is handled automatically by the `prepare_vasp_parameters()` utility function.

---

## Verification

### Monitor Concurrent Jobs

```bash
# Watch process list (should show max 2 VASP jobs)
watch -n 2 'verdi process list'

# Check workgraph status
verdi process show <PK>
verdi process status <PK>
```

### Expected Behavior

With `max_number_jobs=2`:

**Phase 1 (Bulk):**
- 1 VASP job: bulk relaxation

**Phase 2 (References):**
- 2 VASP jobs: metal + oxygen (both run concurrently)

**Phase 3 (Formation Enthalpy):**
- Calculation nodes (no VASP)

**Phase 4 (Slab Calculations):**
- Maximum 2 concurrent VASP jobs
- If 3 slabs: slab_1 and slab_2 run, slab_3 queued
- When slab_1 finishes, slab_3 launches

**Phase 5 (Thermodynamics):**
- Calculation nodes (no VASP)

---

## Workflow Phases

The serial preset implements five phases:

### Phase 1: Bulk Calculations
- Load bulk structure from file
- Relax bulk with VASP
- Extract total energy

### Phase 2: Reference Calculations
- Relax metal reference
- Relax oxygen reference
- Optional: Relax nonmetal reference (for ternary oxides)
- Extract energies for all references

### Phase 3: Thermodynamics Preparation
- Calculate formation enthalpy
- Build reference energies dictionary
- Identify oxide type (binary or ternary)

### Phase 4: Slab Calculations
- SCF calculations for all slabs (unrelaxed energies)
- Relaxation calculations for all slabs (optional)
- Extract energies
- Calculate relaxation energies (optional)

### Phase 5: Thermodynamics Calculations
- Calculate surface energies for each slab
- Select appropriate result (binary vs ternary)

---

## Output Structure

```python
outputs = {
    # Bulk results
    'bulk_energy': Float,
    'bulk_structure': StructureData,

    # Reference results
    'metal_energy': Float,
    'metal_structure': StructureData,
    'oxygen_energy': Float,
    'oxygen_structure': StructureData,
    'nonmetal_energy': Float,  # If ternary oxide
    'nonmetal_structure': StructureData,  # If ternary oxide

    # Thermodynamics results
    'formation_enthalpy': Dict,
    'reference_energies': Dict,
    'oxide_type': Str,  # 'binary' or 'ternary'

    # Slab results
    'unrelaxed_slab_energies': {slab_id: Float},
    'relaxed_slabs': {slab_id: StructureData},  # If relax_slabs=True
    'slab_energies': {slab_id: Float},  # If relax_slabs=True
    'relaxation_energies': {slab_id: Float},  # If compute_relaxation_energy=True

    # Final results
    'surface_energies': {slab_id: Dict},
}
```

---

## Comparison with Standard Preset

| Feature | Standard Preset | Serial Preset |
|---------|----------------|---------------|
| **Graph Structure** | Nested sub-workgraphs | Flat single-level graph |
| **Concurrency Control** | `max_number_jobs` limited effect | `max_number_jobs` controls all VASP |
| **Slab Generation** | Dynamic from bulk | Pre-generated (required) |
| **Node Builders** | `@task.graph` decorators | Direct `wg.add_task()` |
| **Provenance** | Multiple graph levels | Single graph level |
| **Module Location** | `teros.core.workgraph` | `teros.experimental.surface_thermo_preset_serial` |
| **Stability** | Stable production code | Experimental development |

---

## Known Limitations

### 1. No Dynamic Slab Generation

**Issue:** Cannot generate slabs dynamically from bulk structure.

**Reason:** The workflow cannot wait for bulk relaxation to complete before generating slabs in the same graph.

**Workaround:** Pre-generate slabs before building the workflow.

### 2. Pre-Generated Slabs Required

The `input_slabs` parameter is **mandatory**:

```python
wg = surface_thermodynamics_serial_workgraph(
    input_slabs=slabs_dict,  # REQUIRED
    # ...
)
```

**Future:** Two-stage execution may enable dynamic generation.

### 3. Experimental Status

- API may change
- Limited production testing
- Not feature-complete compared to standard preset

---

## Module Structure

```
teros/experimental/surface_thermo_preset_serial/
├── __init__.py              # Exports surface_thermodynamics_serial_workgraph
├── workgraph.py             # Main workgraph (no @task.graph decorator)
├── utils.py                 # Parameter preparation (prepare_vasp_parameters)
├── slab_operations.py       # Slab node builders (SCF, relax, energy extraction)
└── thermodynamics_operations.py  # Thermodynamics node builders
```

### Design Pattern: Node Builders

Node builders are regular Python functions (not decorated with `@task.graph`):

```python
def build_scf_slabs_nodes(wg, slabs, code, parameters, ...):
    """
    Add SCF nodes for each slab directly to main graph.

    Args:
        wg: WorkGraph instance (main graph)
        slabs: Dict of {slab_id: StructureData}
        ...

    Returns:
        Dict of {slab_id: vasp_node}
    """
    scf_nodes = {}
    for slab_id, structure in slabs.items():
        scf_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"scf_slab_{slab_id}",
            structure=structure,
            code=code,
            **parameters,
        )
    return scf_nodes
```

**Key:** Takes `wg` as first parameter and calls `wg.add_task()` directly.

---

## Testing

### Run Example Script

```bash
source ~/envs/aiida/bin/activate
python examples/vasp/step_16_surface_thermodynamics_serial.py
```

The script demonstrates:
1. Loading AiiDA profile
2. Creating test slab structures
3. Building serial workgraph with `max_number_jobs=2`
4. Submitting and monitoring execution

### Expected Timeline

- Bulk and reference calculations: ~15-30 minutes
- Slab calculations (6 jobs with limit of 2): ~45-60 minutes
- Thermodynamics calculations: ~1-2 minutes
- **Total: ~60-90 minutes**

---

## Troubleshooting

### Q: All jobs still submit simultaneously

**A:** Check that you set `wg.max_number_jobs` (not `wg.max_concurrent_jobs`):

```python
# Correct
wg.max_number_jobs = 2

# Wrong (does nothing)
wg.max_concurrent_jobs = 2
```

### Q: Error about missing incar tags

**A:** The serial preset automatically wraps parameters in `{'incar': ...}`. If you see this error, you may have passed parameters already wrapped:

```python
# Correct
bulk_parameters = {'PREC': 'Accurate', ...}

# Wrong (double wrapping)
bulk_parameters = {'incar': {'PREC': 'Accurate', ...}}
```

### Q: Cannot dynamically generate slabs

**A:** This is a known limitation. Pre-generate slabs before building the workflow:

```python
from aiida import orm
from ase.io import read

slabs_dict = {
    'slab_100': orm.StructureData(ase=read('slab_100.cif')),
    'slab_110': orm.StructureData(ase=read('slab_110.cif')),
}
```

---

## Future Development

Planned improvements:

1. **Two-stage execution:** Support dynamic slab generation with bulk-then-slabs approach
2. **Error handling:** Add comprehensive validation and retry logic
3. **API stabilization:** Finalize parameter names and structure
4. **Testing:** Expand production testing across different materials
5. **Migration path:** Define upgrade path from experimental to production

---

## When to Use This Preset

**Use serial preset when:**
- You need strict control over concurrent VASP jobs
- You work on resource-constrained clusters
- You have pre-generated slabs
- You need better provenance tracking

**Use standard preset when:**
- You need dynamic slab generation from bulk
- You prefer stable production code
- You don't need strict concurrency control
- You use unlimited parallel execution

---

## Contributing

This is experimental code actively seeking feedback:

1. Test with your materials systems
2. Report unexpected behavior
3. Propose architectural improvements
4. Submit fixes via pull request

File issues at: (repository issue tracker)

---

## References

- **Full documentation:** `examples/vasp/step_16_DOCUMENTATION.md`
- **Example script:** `examples/vasp/step_16_surface_thermodynamics_serial.py`
- **Module README:** `teros/experimental/surface_thermo_preset_serial/README.md`
- **Standard preset:** `teros/core/workgraph.build_core_workgraph()`
- **Concurrency control:** `docs/CONCURRENCY_CONTROL.md`

---

## API Stability Warning

⚠️ **This is experimental code. The API will change.**

Expect changes to:
- Function names and signatures
- Parameter names
- Output structure
- Module location

Do not use in production workflows without careful testing and version pinning.
