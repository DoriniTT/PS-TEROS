# AIMD Implementation Errors and Investigation

## Date
2025-10-13

## Overview
Attempting to implement AIMD (Ab Initio Molecular Dynamics) functionality in PS-TEROS following the scatter-gather pattern used by other modules (slabs, thermodynamics). The implementation encounters a persistent WorkGraph serialization error.

## Primary Error

### Error Message
```
yaml.constructor.ConstructorError: found unconstructable recursive node
in "<unicode string>", line 1454, column 20:
    _widget: &id004 !!python/object:node_grap ...
             ^
```

### Error Context
- **When**: During WorkGraph construction/serialization, before any calculations run
- **Where**: In `build_core_workgraph()` when creating the workflow
- **State**: WorkGraph state = "Excepted"
- **Observation**: The error occurs even with simplified implementations (single AIMD stage, no nested loops)

## Implementation Attempts

### Attempt 1: Nested @task.graph functions
**Approach**: Created `aimd_sequential_slab` as @task.graph, called from `aimd_slabs_scatter`

**Error**: 
```
Invalid assignment into namespace socket: graph_outputs.outputs
Field 'term_0' is not defined and this namespace is not dynamic.
```

**Reason**: Nested @task.graph functions don't properly expose their outputs to parent graph

---

### Attempt 2: Single @task.graph with nested loops
**Approach**: Inlined all logic into `aimd_slabs_scatter` with for loops (slabs loop + stages loop)

**Error**: 
```
yaml.constructor.ConstructorError: found unconstructable recursive node
```

**Reason**: Unknown - possibly circular dependencies from chaining `prev_structure` and `prev_remote`

---

### Attempt 3: Simplified to single stage
**Approach**: Removed inner loop, only run first AIMD stage per slab

**Error**: 
```
yaml.constructor.ConstructorError: found unconstructable recursive node
```

**Reason**: Same error persists even without nested loops

---

### Attempt 4: Changed parameter types to AiiDA types
**Approach**: Changed function signature to use `orm.List`, `orm.Dict` instead of Python `list`, `dict`

**Status**: In progress

## Code Structure

### Current Implementation Pattern

```python
@task.graph
def aimd_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    aimd_sequence: list,  # or orm.List
    code: orm.Code,
    aimd_parameters: dict,  # or orm.Dict
    ...
) -> t.Annotated[dict, namespace(final_structures=dynamic(...), ...)]:
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)
    
    final_structures = {}
    
    for slab_label, slab_structure in slabs.items():
        # Create VASP task for this slab
        aimd_task = VaspTask(**inputs)
        final_structures[slab_label] = aimd_task.structure
    
    return {'final_structures': final_structures}
```

### Working Reference Pattern (from relax_slabs_scatter)

```python
@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    ...
) -> t.Annotated[dict, namespace(relaxed_structures=dynamic(orm.StructureData), ...)]:
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_wc)
    
    relaxed: dict[str, orm.StructureData] = {}
    
    for label, structure in slabs.items():
        relaxation = relax_task_cls(**inputs)
        relaxed[label] = relaxation.structure
    
    return {'relaxed_structures': relaxed}
```

## Key Differences from Working Code

1. **Parameter Types**: `relax_slabs_scatter` uses `t.Mapping` for dict parameters, not raw `dict`
2. **No Complex Objects**: `relax_slabs_scatter` doesn't take list-of-dicts as input
3. **No Sequential Dependencies**: Each relaxation is independent; AIMD stages chain together
4. **Called Directly**: `relax_slabs_scatter` is called directly in `core_workgraph`, not nested

## Hypotheses

### Hypothesis 1: Python complex objects can't be serialized
**Evidence**: The error mentions YAML serialization of Python objects
**Test**: Change all parameters to AiiDA types (orm.List, orm.Dict)

### Hypothesis 2: Circular reference in task graph
**Evidence**: Error says "recursive node"
**Test**: Remove all chaining (prev_structure, prev_remote)

### Hypothesis 3: WorkGraph can't handle function with many parameters
**Evidence**: `aimd_slabs_scatter` has 9 parameters vs 7-8 for working functions
**Test**: Reduce parameters, use defaults

### Hypothesis 4: Issue with how @task.graph decorator processes the function
**Evidence**: Error occurs during graph construction
**Test**: Try without @task.graph decorator, call inline

## Questions for Deep Investigation

1. What exactly triggers the YAML serialization error?
2. How does WorkGraph handle @task.graph functions with complex parameter types?
3. What creates the "recursive node" - is it the loop, the task chaining, or something else?
4. Why does the same pattern work for relax_slabs_scatter but not aimd_slabs_scatter?
5. Is there a limit to nesting depth or complexity in WorkGraph task graphs?

## Next Steps

1. Compare line-by-line differences between `aimd_slabs_scatter` and `relax_slabs_scatter`
2. Create minimal reproduction case outside PS-TEROS
3. Check WorkGraph documentation for serialization limitations
4. Test with absolutely minimal implementation (1 slab, hardcoded parameters)
5. Consider alternative architectures (separate AIMD workflow, not integrated)
