# Serial Surface Thermodynamics Preset (EXPERIMENTAL)

**Status:** Experimental - In Active Development
**Location:** `teros.experimental.surface_thermo_preset_serial`
**Example:** `examples/vasp/step_16_surface_thermodynamics_serial.py`

## Overview

The serial preset creates a flat-graph architecture where all VASP nodes exist at the same level. This allows `max_number_jobs` to control concurrent execution across all calculations.

## Problem with Standard Preset

The standard surface thermodynamics preset uses nested sub-workgraphs created by `@task.graph` decorators. These create isolated graph contexts where:

- Each sub-workgraph spawns its own VASP calculations
- `max_number_jobs` set on the main graph does not propagate to nested graphs
- All VASP jobs submit simultaneously regardless of the limit

Result: Setting `max_number_jobs=2` on the main graph does not limit concurrent VASP jobs.

## Solution: Flat-Graph Serial Preset

The serial preset replaces nested sub-workgraphs with direct node addition:

```python
# Standard preset (nested graphs)
@task.graph
def relax_slabs(slabs, ...):
    for slab in slabs:
        vasp_node = VaspWorkChain(...)  # Nested context
    return results

# Serial preset (flat graph)
def build_relax_slabs_nodes(wg, slabs, ...):
    for slab in slabs:
        wg.add_task(VaspWorkChain, ...)  # Same level
    return nodes
```

All VASP nodes now exist in the main graph. Setting `wg.max_number_jobs=2` limits concurrent execution to 2 VASP jobs.

## Current Implementation

The preset implements five phases:

1. **Bulk calculations**: Load and relax bulk structure, extract energy
2. **Reference calculations**: Relax metal, oxygen, and optional nonmetal references
3. **Slab handling**: Requires pre-generated slabs (dynamic generation not yet supported)
4. **Slab calculations**: SCF and relaxation for each slab
5. **Thermodynamics**: Calculate formation enthalpy and surface energies

## Usage

```python
from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph
)

# Build workgraph
wg = surface_thermodynamics_serial_workgraph(
    structures_dir='structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    code_label='VASP-6.5.0@bohr-new',
    potential_family='PBE',
    kpoints_spacing=0.4,
    bulk_parameters=bulk_params,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    bulk_options=options,
    metal_parameters=metal_params,
    metal_potential_mapping={'Ag': 'Ag'},
    metal_options=options,
    oxygen_parameters=oxygen_params,
    oxygen_potential_mapping={'O': 'O'},
    oxygen_options=options,
    input_slabs=slabs_dict,  # REQUIRED
    slab_parameters=slab_params,
    slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    slab_options=options,
    relax_slabs=True,
    compute_thermodynamics=True,
)

# Set concurrency limit
wg.max_number_jobs = 2

# Submit
wg.submit()
```

## Required Input: Pre-Generated Slabs

The serial preset requires pre-generated slabs as `input_slabs` dictionary:

```python
from aiida import orm
from ase.io import read

slabs_dict = {
    'slab_100_term_0': orm.StructureData(ase=read('slab_100_term_0.cif')),
    'slab_100_term_1': orm.StructureData(ase=read('slab_100_term_1.cif')),
    'slab_110': orm.StructureData(ase=read('slab_110.cif')),
}
```

Dynamic slab generation is not yet supported. The workflow cannot wait for bulk relaxation to complete before generating slabs in the same graph.

## Verification

Monitor concurrent jobs to verify `max_number_jobs` works:

```bash
# Watch process list (should show max 2 VASP jobs running)
watch -n 2 'verdi process list'

# Check workgraph status
verdi process show <PK>
verdi process status <PK>
```

Expected behavior:
- Initial phase: 1-2 VASP jobs (bulk + metal, or bulk + metal + oxygen)
- Slab SCF phase: Maximum 2 concurrent VASP jobs
- Slab relaxation phase: Maximum 2 concurrent VASP jobs

## Known Limitations

1. **No dynamic slab generation**: Must provide pre-generated slabs
2. **Experimental status**: Implementation may change
3. **Limited testing**: Requires more production validation
4. **Parameter format**: Must wrap VASP parameters as `{'incar': {...}}`

## Differences from Standard Preset

| Feature | Standard Preset | Serial Preset |
|---------|----------------|---------------|
| Graph structure | Nested sub-workgraphs | Flat single-level |
| `max_number_jobs` | Does not propagate | Controls all VASP jobs |
| Slab generation | Dynamic (from bulk) | Pre-generated (required) |
| Node builders | `@task.graph` decorators | Direct `wg.add_task()` calls |
| Provenance | Multiple graph levels | Single graph level |

## Future Development

Planned improvements:

- Support dynamic slab generation with two-stage execution
- Add error handling and validation
- Implement retry logic for failed calculations
- Add more comprehensive testing
- Document migration path from standard preset

## Example Output

The workflow returns:

```python
outputs = {
    'bulk_energy': Float,
    'bulk_structure': StructureData,
    'metal_energy': Float,
    'oxygen_energy': Float,
    'metal_structure': StructureData,
    'oxygen_structure': StructureData,
    'formation_enthalpy': Dict,
    'reference_energies': Dict,
    'oxide_type': Str,
    'surface_energies': {slab_id: Dict},
    'relaxed_slabs': {slab_id: StructureData},
    'slab_energies': {slab_id: Float},
    'unrelaxed_slab_energies': {slab_id: Float},
    'relaxation_energies': {slab_id: Float},
}
```

## Testing

Run the example script:

```bash
source ~/envs/aiida/bin/activate
python examples/vasp/step_16_surface_thermodynamics_serial.py
```

The script:
1. Loads the presto profile
2. Creates 3 test slab structures
3. Builds the serial workgraph with `max_number_jobs=2`
4. Submits and monitors execution

Expected timeline:
- Bulk and reference calculations: ~15-30 minutes
- Slab calculations (6 jobs with limit of 2): ~45-60 minutes
- Thermodynamics calculations: ~1-2 minutes
- Total: ~60-90 minutes

## Technical Implementation

### Module Structure

```
teros/experimental/surface_thermo_preset_serial/
├── __init__.py              # Exports main workgraph function
├── workgraph.py             # Main workgraph logic
├── utils.py                 # Parameter preparation helpers
├── slab_operations.py       # Slab calculation node builders
└── thermodynamics_operations.py  # Thermodynamics node builders
```

### Node Builder Pattern

Node builders are regular Python functions that take the workgraph as first parameter:

```python
def build_scf_slabs_nodes(wg, slabs, code, parameters, ...):
    """Add SCF nodes directly to main graph."""
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

This pattern keeps all nodes at the same graph level.

## Contributing

This is experimental code. Report issues or suggest improvements:

1. Test with your materials system
2. Document unexpected behavior
3. Propose architectural improvements
4. Submit fixes via pull request

## References

- Standard preset: `teros.core.workgraph.build_core_workgraph()`
- Core thermodynamics: `teros.core.thermodynamics`
- Core formation enthalpy: `teros.core.hf`
- Example: `examples/vasp/step_16_surface_thermodynamics_serial.py`
