# Surface Thermodynamics Serial Preset Design

**Date:** 2025-11-02
**Author:** Claude Code
**Status:** Design Phase

## Purpose

The `max_concurrent_jobs` parameter fails to limit concurrent VASP calculations in the current `surface_thermodynamics` preset because nested sub-workgraphs created by `@task.graph` decorators isolate calculations from the main graph's concurrency control.

This design proposes a serial variant that flattens the graph structure. All VASP calculation nodes will exist at the same graph level, allowing `max_concurrent_jobs` to control their execution.

## Design Constraints

- **API compatibility required**: Users must call the serial preset the same way they call the current preset from `build_core_workgraph()`
- **Scope**: Implement core thermodynamics workflow (bulk + slabs + surface energies), design for future extension
- **Success criteria**: Setting `max_concurrent_jobs=2` must limit concurrent VASP jobs to 2

## Architecture

### Core Principle: Direct Node Addition

The current implementation uses `@task.graph` decorators that create nested sub-workgraphs:

```python
# Current approach (nested)
@task.graph
def relax_slabs_scatter(slabs, code, ...):
    for slab_id, slab in slabs.items():
        vasp_task = VaspWorkChain(structure=slab, ...)
    return namespace
```

The serial approach removes the `@task.graph` decorator and converts these to regular Python functions that add nodes directly to the main graph:

```python
# Serial approach (flat)
def build_relax_slabs_nodes(wg, slabs, code, ...):
    """Add VASP relax nodes directly to wg"""
    vasp_nodes = {}
    for slab_id, slab_data in slabs.items():
        vasp_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"relax_{slab_id}",
            structure=slab_data,
            ...
        )
    return vasp_nodes
```

All VASP nodes exist at the same graph level. The WorkGraph engine applies `max_concurrent_jobs` to all of them.

### Module Structure

Create the new module in `teros/experimental/surface_thermo_preset_serial/`:

```
teros/experimental/surface_thermo_preset_serial/
├── __init__.py                      # Export surface_thermodynamics_serial_workgraph
├── workgraph.py                     # Main @task.graph entry point
├── slab_operations.py               # Node builders for slab calculations
├── thermodynamics_operations.py    # Node builders for thermodynamics
└── utils.py                         # Shared utilities
```

**File responsibilities:**

- `workgraph.py`: Contains `surface_thermodynamics_serial_workgraph()`, the main orchestrator
- `slab_operations.py`: Functions like `build_relax_slabs_nodes()`, `build_scf_slabs_nodes()`, `build_energy_extraction_nodes()`
- `thermodynamics_operations.py`: Functions like `build_surface_energy_nodes()`
- `utils.py`: Parameter preparation helpers

### Workflow Sequence

**Phase 1: Bulk calculations**
1. Load bulk structure
2. Add VASP bulk relax node → `bulk_structure`, `bulk_energy`
3. Add reference calculations (metal, oxygen, optional nonmetal) → energies
4. Add formation_enthalpy calcfunction node

**Phase 2: Slab generation and relaxation**
5. Generate or wrap input slabs → `slab_structures` dict
6. For each slab:
   - Add VASP SCF node (unrelaxed) → `unrelaxed_energy`
   - Add VASP relax node → `relaxed_structure`, `relaxed_energy`
   - Add relaxation_energy calcfunction node

**Phase 3: Thermodynamics calculations**
7. Add oxide_type identification node
8. For each relaxed slab, add surface_energy calcfunction node
9. Collect results

### Output Collection

The serial approach builds output dictionaries manually:

```python
def surface_thermodynamics_serial_workgraph(...):
    # Build all nodes
    bulk_node = build_bulk_node(wg, ...)
    relax_nodes = build_relax_slabs_nodes(wg, slabs, ...)
    energy_nodes = build_energy_extraction_nodes(wg, relax_nodes, ...)
    surface_nodes = build_surface_energy_nodes(wg, bulk_node, energy_nodes, ...)

    # Collect outputs
    relaxed_slabs = {}
    slab_energies = {}
    surface_energies = {}

    for slab_id in slabs.keys():
        relaxed_slabs[slab_id] = relax_nodes[slab_id].outputs.structure
        slab_energies[slab_id] = energy_nodes[slab_id].outputs.result
        surface_energies[slab_id] = surface_nodes[slab_id].outputs.result

    # Return outputs matching current preset API
    return {
        'bulk_energy': bulk_node.outputs.energy,
        'bulk_structure': bulk_node.outputs.structure,
        'slab_structures': slab_structures,
        'relaxed_slabs': relaxed_slabs,
        'slab_energies': slab_energies,
        'surface_energies': surface_energies,
        # ... remaining outputs
    }
```

### Code Reuse

The serial implementation will reuse existing code:

- **Calcfunctions**: Reuse all calcfunctions from `teros.core.thermodynamics` (`calculate_surface_energy_ternary`, `calculate_surface_energy_binary`, etc.)
- **VASP setup**: Reuse parameter preparation logic from current implementation
- **Slab generation**: Reuse `generate_slab_structures` from `teros.core.slabs`

Only the graph construction pattern changes. The underlying calculations remain identical.

### Extensibility Design

The design supports adding features (cleavage, electronic properties, AIMD) later:

1. **Modular node builders**: Each feature gets its own file:
   - Current: `slab_operations.py`, `thermodynamics_operations.py`
   - Future: `cleavage_operations.py`, `electronic_operations.py`, `aimd_operations.py`

2. **Conditional inclusion** in main workgraph:
   ```python
   def surface_thermodynamics_serial_workgraph(
       ...,
       compute_cleavage=False,
       compute_electronic_properties=False,
       ...
   ):
       # Core workflow always runs
       bulk_nodes = build_bulk_nodes(wg, ...)
       slab_nodes = build_relax_slabs_nodes(wg, ...)
       thermo_nodes = build_surface_energy_nodes(wg, ...)

       # Optional features
       if compute_cleavage:
           cleavage_nodes = build_cleavage_nodes(wg, slab_nodes, ...)

       if compute_electronic_properties:
           bands_nodes = build_bands_nodes(wg, bulk_nodes, slab_nodes, ...)
   ```

3. **Consistent outputs**: Always return the same output dictionary keys. Use `None` or empty dicts for disabled features.

## Testing Strategy

Test the serial preset with `max_concurrent_jobs=2`:

1. Create test script in `examples/vasp/test_surface_thermo_serial/`
2. Configure workflow with multiple slabs (at least 4 slabs to see concurrent behavior)
3. Set `max_concurrent_jobs=2` on the WorkGraph
4. Run workflow: `verdi daemon restart && python test_script.py`
5. Monitor execution: `verdi process list` should show maximum 2 VASP jobs running simultaneously
6. Verify completion: Main node returns `[0]`, surface energies calculated correctly

## Implementation Notes

- Keep the current `surface_thermodynamics` preset unchanged
- Place serial variant in `teros/experimental/` for testing
- If successful, consider migrating other presets to flat structure
- Document differences between nested and flat approaches

## References

- AiiDA-WorkGraph graph_task concept: https://aiida-workgraph.readthedocs.io/en/latest/concept/autogen/graph_task_concept.html
- Current max_concurrent_jobs implementation: `docs/plans/2025-11-02-max-concurrent-jobs-implementation.md`

---

## Implementation Status

**Date:** 2025-11-02
**Status:** Implemented

### Completed Components

- [x] Module structure created
- [x] Utilities module (parameter preparation)
- [x] Slab operations module (node builders)
- [x] Thermodynamics operations module (node builders)
- [x] Main workgraph implementation
- [x] Test script created

### Known Limitations

1. **Dynamic slab generation not supported**: Current implementation requires `input_slabs` to be provided as a dictionary. Dynamic slab generation from `miller_indices` requires additional work to handle the dynamic namespace properly.

2. **Simplified oxide type handling**: Currently creates both binary and ternary surface energy nodes and selects the appropriate one. Could be optimized to create only the needed type.

### Testing

Test script location: `examples/vasp/test_surface_thermo_serial/test_serial_preset.py`

Verification steps:
1. Run test with `max_concurrent_jobs=2`
2. Monitor with `verdi process list` to confirm max 2 concurrent VASP jobs
3. Verify successful completion with exit code [0]
4. Check surface energies calculated correctly

### Next Steps

1. Test with real data
2. Add support for dynamic slab generation
3. Optimize oxide type handling
4. Add support for optional features (cleavage, electronic properties, AIMD)
5. If successful, consider migrating other presets to flat structure
