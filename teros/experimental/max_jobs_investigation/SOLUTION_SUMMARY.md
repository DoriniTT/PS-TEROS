# Solution Summary: Making max_number_jobs Work in Nested WorkGraphs

## The Problem You Had

Setting `max_number_jobs` on the main WorkGraph didn't limit concurrent VASP calculations because nested sub-workgraphs (created by `@task.graph`) spawn their own calculations independently.

## The Solution We Found

✅ **Use `get_current_graph()` to access and configure the WorkGraph from within `@task.graph` functions!**

## How to Fix Your Code

###  Step 1: Add max_number_jobs Parameter

Add `max_number_jobs` as a parameter to your `@task.graph` functions:

```python
@task.graph
def relax_slabs_scatter(
    slabs,
    code,
    potential_family,
    potential_mapping,
    parameters,
    options,
    kpoints_spacing=None,
    clean_workdir=True,
    max_number_jobs: int = None,  # ← Add this parameter
):
```

### Step 2: Set max_number_jobs Using get_current_graph()

At the START of the function, get the current WorkGraph and set max_number_jobs:

```python
@task.graph
def relax_slabs_scatter(..., max_number_jobs: int = None):
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph
    if max_number_jobs is not None:
        wg = get_current_graph()
        # Convert AiiDA Int node to Python int
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # Rest of your function unchanged...
    for slab_id, structure in slabs.items():
        vasp_calc = VaspWorkChain(...)
    return results
```

### Step 3: Pass max_number_jobs from Main WorkGraph

In your main workgraph, pass `max_number_jobs` to the `@task.graph` function:

```python
# In surface_thermodynamics_workgraph function
relaxation_outputs = relax_slabs_scatter(
    slabs=slab_namespace,
    code=code,
    potential_family=potential_family,
    potential_mapping=slab_pot_map,
    parameters=slab_params,
    options=slab_opts,
    kpoints_spacing=slab_kpts,
    clean_workdir=clean_workdir,
    max_number_jobs=orm.Int(max_concurrent_jobs),  # ← Pass it here!
)
```

## Why This Works

When `@task.graph` executes your function, it creates a WorkGraph and sets it as the "current graph" using a context manager. The `get_current_graph()` function gives you access to that WorkGraph instance, allowing you to configure it (like setting `max_number_jobs`).

## Which Functions Need This Fix?

Apply this pattern to ALL your `@task.graph` functions that create VASP calculations:

1. `relax_slabs_scatter` in `teros/core/slabs.py`
2. `scf_slabs_scatter` in `teros/core/slabs.py`
3. `compute_surface_energies_scatter` in `teros/core/surface_energy.py`
4. `scf_relax_and_calculate_relaxation_energy` in `teros/core/slabs.py`
5. Any other `@task.graph` functions that create child workgraphs

## Testing

After applying the fix:

1. Set `max_number_jobs` on your main workgraph:
   ```python
   wg = surface_thermodynamics_workgraph(...)
   wg.max_number_jobs = 4  # Limit to 4 concurrent jobs
   wg.submit()
   ```

2. Monitor the running calculations:
   ```bash
   watch -n 2 'verdi process list -a -p 1 | head -30'
   ```

3. You should see no more than 4 VASP calculations running simultaneously!

## Complete Example

```python
# In teros/core/slabs.py

@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    restart_folders: dict = None,
    max_number_jobs: int = None,  # ← NEW PARAMETER
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
    remote_folders=dynamic(orm.RemoteData),
)]:
    """Relax multiple slab structures with concurrency control."""
    from aiida_workgraph import get_current_graph
    from aiida.plugins import WorkflowFactory

    # ← NEW: Set max_number_jobs
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # Rest is unchanged...
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    vasp_task_cls = task(vasp_wc)

    structures_dict: dict[str, orm.StructureData] = {}
    energies_dict: dict[str, orm.Float] = {}
    remote_dict: dict[str, orm.RemoteData] = {}

    for label, structure in slabs.items():
        inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            inputs['kpoints_spacing'] = kpoints_spacing
        if restart_folders is not None and label in restart_folders:
            inputs['restart_folder'] = restart_folders[label]

        calc = vasp_task_cls(**inputs)
        structures_dict[label] = calc.structure
        energies_dict[label] = extract_total_energy(energies=calc.misc).result
        remote_dict[label] = calc.remote_folder

    return {
        'relaxed_structures': structures_dict,
        'energies': energies_dict,
        'remote_folders': remote_dict,
    }
```

## API Reference

- **Function**: `get_current_graph()`
- **Module**: `aiida_workgraph` or `aiida_workgraph.manager`
- **Returns**: The currently active WorkGraph instance
- **Use**: Inside `@task.graph` functions to access the WorkGraph being built

## Benefits

✅ Simple - just 3 lines of code per function
✅ No architecture changes needed
✅ Works with existing `@task.graph` pattern
✅ Can propagate through multiple nesting levels
✅ Uses official AiiDA WorkGraph API

---

**Investigation**: `/home/thiagotd/git/PS-TEROS/teros/experimental/max_jobs_investigation/`
**Test Results**: See `SUCCESS.md` for full details
**Date**: 2025-11-02
