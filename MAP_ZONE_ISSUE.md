# Map/Zone Issue with Dynamic Task Outputs in AiiDA-WorkGraph

## Problem Summary

We are unable to use the `Map` or `Zone` context managers with **dynamic task outputs** as the source data in the latest version of AiiDA-WorkGraph. The legacy `active_map_zone` from an older version worked perfectly for this use case, but we cannot find the equivalent implementation in the current version.

## Our Use Case

We need to:
1. Generate multiple slab structures from a relaxed bulk structure (dynamic number of slabs)
2. Relax **all slab structures in parallel** using VASP
3. Extract energies from each relaxed slab

The number of slabs is not known until runtime (it depends on the symmetry and terminations of the bulk structure).

## What Worked: Legacy `active_map_zone`

In our legacy code (AiiDA-WorkGraph older version), we successfully used `active_map_zone` with task outputs:

```python
# From: /home/thiagotd/git/PS-TEROS/legacy/v2/teros/core/workgraph.py

from aiida_workgraph import WorkGraph, active_map_zone

with WorkGraph(workgraph_name) as wg:
    # Generate slabs from relaxed bulk
    generate_slabs_task = wg.add_task(
        get_slabs,
        name="generate_slabs",
        relaxed_structure=bulk_task.outputs.structure,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        # ... other parameters
    )

    # Use task output as source for active_map_zone
    slab_structures = generate_slabs_task.outputs.structures

    # Map over dynamically generated slabs
    with active_map_zone(slab_structures) as map_zone:
        slab_task = map_zone.add_task(
            VaspWorkChain,
            name="slab_relaxation",
            structure=map_zone.item,  # Each slab structure
        )
        slab_task.set_from_builder(builder_slab)

    # Use outputs from map zone
    slab_structures_out = slab_task.outputs.structure
    slab_parameters_out = slab_task.outputs.misc
```

**This worked perfectly** - it:
- ✅ Accepted task outputs as the source (`generate_slabs_task.outputs.structures`)
- ✅ Created separate VASP tasks for each generated slab
- ✅ Ran all slab relaxations in parallel
- ✅ Collected results in a dynamic namespace

## What We've Tried with Current Version

### Attempt 1: Using `Map` with task outputs

Following the documentation pattern, we tried:

```python
from aiida_workgraph import WorkGraph, Map, task

VaspTask = task(WorkflowFactory('vasp.v2.vasp'))

with WorkGraph(name='SlabRelax') as wg:
    # Generate slabs
    slabs_task = wg.add_task(
        get_slabs,
        name="generate_slabs",
        relaxed_structure=bulk_vasp.outputs.structure,
        # ... parameters
    )

    # Try to use Map with task output
    with Map(slabs_task.outputs.slabs) as map_zone:
        vasp_result = VaspTask(
            structure=map_zone.item,  # or map_zone.item.value
            code=code,
            # ... parameters
        )
```

**Result:** `KeyError: 'source'`

```
File "/home/thiagotd/envs/aiida/lib/python3.10/site-packages/aiida_workgraph/engine/task_manager.py", line 325, in execute_map_task
    source = kwargs["source"]
KeyError: 'source'
```

### Attempt 2: Using `Zone` instead of `Map`

```python
from aiida_workgraph import WorkGraph, Zone, task

with WorkGraph(name='SlabRelax') as wg:
    slabs_task = wg.add_task(get_slabs, ...)

    with Zone() as slab_zone:
        # Try to iterate over task outputs
        for slab_id, slab_structure in slabs_task.outputs.slabs.items():
            vasp_task = VaspTask(
                structure=slab_structure,
                code=code,
                # ... parameters
            )
```

**Result:** `TypeError: 'illegal operation on a future value (Socket): iteration'`

Cannot iterate over task outputs in `@task.graph` or context manager.

### Attempt 3: Using `@task.graph` with iteration

```python
@task.graph(outputs=['relaxed_slabs', 'slab_energies'])
def relax_slabs(slabs, code, ...):
    VaspTask = task(WorkflowFactory('vasp.v2.vasp'))

    relaxed_slabs = {}
    vasp_tasks = {}

    # Extract and iterate over slabs
    if hasattr(slabs.slabs, "_get_keys"):
        slab_items = [(key, getattr(slabs.slabs, key)) for key in slabs.slabs._get_keys()]
    else:
        slab_items = list(slabs.slabs.items())

    for slab_id, slab_structure in slab_items:
        vasp_tasks[slab_id] = VaspTask(
            structure=slab_structure,
            code=code,
            # ... parameters
        )

    return {'relaxed_slabs': relaxed_slabs, 'slab_energies': energies}
```

**Result:** Tasks are created in Python but **never submitted to AiiDA**. The workflow finishes without running any slab relaxations.

```
2025-10-04 13:15:22 [2157 | REPORT]: [6339|WorkGraphEngine|update_task_state]: Task: get_slabs, type: PYFUNCTION, finished.
2025-10-04 13:15:22 [2158 | REPORT]: [6339|WorkGraphEngine|continue_workgraph]: tasks ready to run:
2025-10-04 13:15:22 [2159 | REPORT]: [6339|WorkGraphEngine|finalize]: Finalize workgraph.
```

No slab relaxation tasks appear in the "Called" section of `verdi process show`.

### Attempt 4: Pre-loading data before WorkGraph

The **only pattern that works** with current Map is pre-defined data:

```python
from aiida_workgraph import Map, task

# Pre-defined data (NOT from task outputs)
len_list = 4
data = {f'data_{i}': {'x': i, 'y': i} for i in range(len_list)}

with WorkGraph('AddMap') as wg:
    with Map(data) as map_zone:  # Works with pre-defined data
        result = add(
            x=get_value(map_zone.item.value, 'x').result,
            y=get_value(map_zone.item.value, 'y').result,
        ).result
        map_zone.gather({'result': result})
```

This works, but **requires knowing the data before building the WorkGraph**, which defeats the purpose for our use case.

## Key Differences

| Feature | Legacy `active_map_zone` | Current `Map`/`Zone` |
|---------|-------------------------|----------------------|
| Task outputs as source | ✅ Works | ❌ `KeyError: 'source'` |
| Pre-defined data as source | ✅ Works | ✅ Works |
| Integration with WorkGraph | ✅ Seamless | ❌ Broken for task outputs |
| Parallel execution | ✅ Automatic | ❌ Never launches tasks |

## Questions

1. **How do we use `Map` or `Zone` with task outputs as the source?**
   - The documentation only shows examples with pre-defined data
   - We need to map over data generated by an upstream task

2. **Why was `active_map_zone` removed?**
   - It worked perfectly for our use case
   - Is there a migration guide from `active_map_zone` to the new Map/Zone?

3. **Is there a different pattern we should use?**
   - Should we use `While` loops instead?
   - Is there a `graph_builder` pattern that supports this?
   - Should we create two separate workflows (one to generate slabs, one to relax them)?

## Simple Reproducible Example

```python
from aiida import orm, load_profile
from aiida_workgraph import WorkGraph, Map, task
from aiida.plugins import WorkflowFactory

load_profile()

# Task that generates dynamic data (simulates slab generation)
@task
def generate_structures():
    """Generate multiple structures (simulates get_slabs)."""
    from ase.build import bulk
    structures = {}
    for i in range(3):
        atoms = bulk('Cu', 'fcc', a=3.6 + i*0.1)
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)
    return orm.Dict(dict={'structures': structures})

# Main workflow
VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
VaspTask = task(VaspWorkChain)
code = orm.load_code('VASP-VTST-6.4.3@bohr')

with WorkGraph(name='TestMapWithTaskOutputs') as wg:
    # Generate structures dynamically
    gen_task = wg.add_task(
        generate_structures,
        name='generate_structures'
    )

    # Try to map over task outputs
    with Map(gen_task.outputs.result) as map_zone:  # This fails
        vasp = VaspTask(
            structure=map_zone.item,
            code=code,
            # ... parameters
        )

# This gives: KeyError: 'source'
wg.submit()
```

## Expected Behavior

We expect the `Map` context manager to:
1. Accept task outputs as source (like `active_map_zone` did)
2. Wait for the upstream task to complete
3. Extract the dynamic namespace/dictionary from the task output
4. Create and launch one task per item in the source
5. Run all tasks in parallel
6. Collect results in a dynamic namespace

## Current Behavior

- ❌ `Map` with task outputs: `KeyError: 'source'`
- ❌ `Zone` with iteration: Cannot iterate over task outputs
- ❌ `@task.graph` with iteration: Tasks created but never submitted
- ✅ `Map` with pre-defined data: Works, but not useful for our case

## Environment

- **AiiDA-WorkGraph version:** (latest, installed via pip)
- **Python version:** 3.10
- **AiiDA-core version:** Latest
- **Legacy code location:** `/home/thiagotd/git/PS-TEROS/legacy/v2/teros/core/workgraph.py`
- **Current code location:** `/home/thiagotd/git/PS-TEROS/teros/workgraph.py` and `workgraph_map.py`

## Request for Help

We need guidance on the correct pattern for mapping over dynamic task outputs in the current version of AiiDA-WorkGraph. The legacy `active_map_zone` worked perfectly, but we cannot find the equivalent implementation in the documentation or codebase.

Any help would be greatly appreciated!

---

## UPDATE: Forum Response and Further Testing

After posting this on the AiiDA forum, we received the response:

> "In the latest version, map_zone.item has key and value outputs. For your example code, you only need to change map_zone.item to map_zone.item.value in both cases (static and dynamic)."

### What We Tried

We updated our code in `/home/thiagotd/git/PS-TEROS/teros/test_modules/zone_approach/workgraph.py` to use `map_zone.item.value`:

```python
with Map(slab_namespace) as map_zone:
    relaxation = relax_single_slab(
        structure=map_zone.item.value,  # Changed from map_zone.item
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
    )

    map_zone.gather({
        'relaxed_slabs': relaxation.relaxed_structure,
        'slab_energies': relaxation.energy,
    })
```

**Result: Still produces the same `KeyError: 'source'` error!**

### Additional Attempts

We also tried:
- `with Map(source_socket=slab_namespace)` - same error
- Both positional and keyword arguments - same error  
- Different import statements - same error

### Error Traceback

```
File "/home/thiagotd/envs/aiida/lib/python3.10/site-packages/aiida_workgraph/engine/task_manager.py", line 298, in execute_map_task
    source = kwargs['source']
KeyError: 'source'
```

This suggests that the Map context manager is not properly registering task outputs as the source parameter internally in version 1.0.0b3.

### Confirmed Working vs Not Working

**Confirmed Working:**
- ✅ Map with pre-defined data (dict created before WorkGraph)
- ✅ `map_zone.item.value` syntax (when source is pre-defined data)

**Confirmed NOT Working (all produce `KeyError: 'source'`):**
- ❌ Map with task outputs in version 1.0.0b3
- ❌ Both `Map(task_output)` and `Map(source_socket=task_output)`
- ❌ Both `map_zone.item` and `map_zone.item.value` (when source is task output)
- ❌ Using `map_zone.gather()` with task outputs

### Test File

We created a test file at `/home/thiagotd/git/PS-TEROS/teros/test_modules/zone_approach/slabs_relax.py` that demonstrates the issue. Running it produces:

```
[6583|WorkGraphEngine|update_task_state]: Task: generate_slab_structures, type: PYFUNCTION, finished.
[6583|WorkGraphEngine|continue_workgraph]: tasks ready to run: map_zone
[6583|WorkGraphEngine|on_except]: KeyError: 'source'
```

### Conclusion

**This appears to be a limitation or bug in AiiDA-WorkGraph 1.0.0b3.**

The Map context manager cannot accept task outputs as the source, despite:
1. The forum response suggesting it should work with `map_zone.item.value`
2. The code being structured correctly according to documentation examples
3. Static data working fine with the same pattern

### Next Steps

**Options to consider:**
1. **Downgrade to version 0.5.2** where `active_map_zone` worked
2. **Wait for a fix** in a future aiida-workgraph version  
3. **Open a GitHub issue** with this reproducible example
4. **Use a different pattern** (two-stage workflow, While loops, etc.)

The legacy `active_map_zone` from version 0.5.2 worked perfectly for this exact use case, so downgrading might be the most practical solution until the Map zone is fixed.
