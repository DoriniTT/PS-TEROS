# Serial Surface Thermodynamics Preset (EXPERIMENTAL)

⚠️ **WARNING: EXPERIMENTAL CODE IN ACTIVE DEVELOPMENT** ⚠️

This module is experimental and under active development. The API may change without notice.

## Purpose

This preset creates a flat-graph architecture for surface thermodynamics workflows where all VASP calculation nodes exist at the same graph level. This allows `max_number_jobs` to control concurrent execution across all calculations.

## Why This Exists

The standard surface thermodynamics preset uses nested sub-workgraphs created by `@task.graph` decorators. Setting `max_number_jobs` on the main graph does not limit concurrent VASP jobs because each nested sub-workgraph spawns its own calculations independently.

This serial preset solves that problem by replacing nested sub-workgraphs with direct node addition to the main graph.

## Status

- ✅ Core functionality implemented
- ✅ Basic testing completed (step_16 example)
- ✅ `max_number_jobs` concurrency control verified
- ⚠️ Dynamic slab generation not yet supported
- ⚠️ Limited production testing
- ⚠️ API may change

## Documentation

See comprehensive documentation:
- **Full guide**: `examples/vasp/step_16_DOCUMENTATION.md`
- **Example script**: `examples/vasp/step_16_surface_thermodynamics_serial.py`

## Quick Start

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
    input_slabs=slabs_dict,  # Pre-generated slabs (required)
    # ... other parameters
)

# Set concurrency limit
wg.max_number_jobs = 2

# Submit
wg.submit()
```

## Known Limitations

1. **Pre-generated slabs required**: Cannot dynamically generate slabs from bulk structure
2. **Experimental status**: API and implementation may change
3. **Limited validation**: Needs more production testing
4. **Parameter format**: VASP parameters must be wrapped as `{'incar': {...}}`

## Module Structure

```
surface_thermo_preset_serial/
├── __init__.py              # Module exports
├── workgraph.py             # Main workgraph logic (no @task.graph)
├── utils.py                 # Parameter preparation
├── slab_operations.py       # Slab calculation node builders
└── thermodynamics_operations.py  # Thermodynamics node builders
```

## Design Pattern

All node builders are regular Python functions that take the workgraph as first parameter and add nodes directly using `wg.add_task()`:

```python
def build_scf_slabs_nodes(wg, slabs, code, parameters, ...):
    """Add SCF nodes directly to main graph."""
    nodes = {}
    for slab_id, structure in slabs.items():
        nodes[slab_id] = wg.add_task(VaspWorkChain, ...)
    return nodes
```

This keeps all nodes at the same graph level, allowing `max_number_jobs` to control all VASP calculations.

## Contributing

Report issues, test with your systems, and propose improvements. This is experimental code actively seeking feedback.

## Migration Path

When this module stabilizes, it may:
- Replace the standard preset
- Become an option in the standard preset (`use_flat_graph=True`)
- Remain as an alternative preset

The decision depends on production testing results and community feedback.
