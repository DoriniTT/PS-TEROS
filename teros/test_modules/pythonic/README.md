# Pythonic Scatter-Gather Slab Relaxation Workflow

This directory demonstrates the **scatter-gather pattern** for parallelizing slab relaxation calculations using AiiDA-WorkGraph's pythonic approach.

## Overview

The scatter-gather pattern is implemented using `@task.graph` decorators, which allow you to:
1. **Scatter**: Distribute work across multiple independent tasks that run in parallel
2. **Gather**: Collect results from all parallel tasks

This approach is more flexible than the Map-zone pattern and uses plain Python control flow to build workflows.

## Key Differences from Map-Zone Approach

| Feature | Pythonic Scatter-Gather | Map-Zone |
|---------|------------------------|----------|
| Pattern | `@task.graph` with Python loops | `with Map(...)` context manager |
| Code Style | Plain Python | DSL-style context manager |
| Flexibility | High - full Python power | Limited to Map iteration |
| Provenance | Tasks nested in graph | Flat workflow structure |
| Learning Curve | Familiar Python patterns | New DSL to learn |

## Files

- **`slabs_relax.py`**: CLI driver that launches the workflow (supports `--mock` mode for testing)
- **`workgraph.py`**: Workflow implementation using scatter-gather pattern
- **`README.md`**: This file

## Implementation Details

### Core Components

1. **`generate_slab_structures`** (`@task.calcfunction`)
   - Generates slab terminations from bulk structure
   - Returns dynamic namespace with all generated slabs
   - Uses pymatgen's `SlabGenerator`

2. **`relax_slabs_scatter`** (`@task.graph`)
   - **Scatter phase**: Creates independent VASP WorkChain tasks for each slab
   - All relaxations run in parallel automatically (no dependencies)
   - **Gather phase**: Collects relaxed structures and energies
   - Uses `@task.graph` to handle dynamic slab dictionary

3. **`extract_total_energy`** (`@task.calcfunction`)
   - Extracts total energy from VASP outputs
   - Tries multiple keys for compatibility

4. **`slab_relaxation_scatter_gather`** (`@task.graph`)
   - Top-level workflow orchestrating the entire process
   - Loads code, generates slabs, then scatters relaxations

### Key Patterns

#### Using `@task.calcfunction` for AiiDA Data Types

When functions return AiiDA Data types (e.g., `orm.Float`, `orm.StructureData`), use `@task.calcfunction` instead of `@task`:

```python
@task.calcfunction
def generate_slab_structures(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    ...
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    # Generate slabs
    return {'slabs': slab_dict}
```

#### Accessing Task Outputs

When a `@task` or `@task.calcfunction` returns a namespace annotation, the outputs are accessible directly (no `.result` wrapper):

```python
# This task returns namespace(values=dynamic(orm.Float))
data = generate_mock_scalars(count=orm.Int(3)).values

# NOT: data = generate_mock_scalars(...).result.values
```

#### Scatter-Gather with `@task.graph`

The key pattern is wrapping the loop in a `@task.graph` because you're iterating over future values:

```python
@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    ...
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float)
)]:
    relaxed = {}
    energies_ns = {}
    
    # Scatter: create independent tasks (runs in parallel)
    for label, structure in slabs.items():
        relaxation = relax_task_cls(structure=structure, ...)
        relaxed[label] = relaxation.structure
        energies_ns[label] = extract_total_energy(relaxation.misc).result
    
    # Gather: return collected results
    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
    }
```

## Running the Workflow

### Prerequisites

```bash
source ~/envs/psteros/bin/activate
verdi profile set-default psteros
verdi daemon restart
```

### Mock Mode (Testing)

Test the scatter-gather pattern without VASP:

```bash
python slabs_relax.py --mock --mock-count 3 --mock-delta 0.5
```

Output:
```
Mock scatter-gather workflow finished:
  Source values:
  value_00: ... value: 0.0
  value_01: ... value: 1.0
  value_02: ... value: 2.0
  Shifted values:
  value_00: ... value: 0.5
  value_01: ... value: 1.5
  value_02: ... value: 2.5
```

### Full VASP Workflow

Run the complete slab relaxation workflow:

```bash
python slabs_relax.py
```

This will:
1. Load the Ag₃PO₄ bulk structure
2. Generate all symmetrically distinct (100) slab terminations
3. Launch VASP relaxations for all slabs **in parallel**
4. Extract and store total energies
5. Return all results with full provenance

## Example Output

```
Workflow finished:
  PK: 11739
  Generated slabs: 4 terminations
    Labels: ['slab_00', 'slab_01', 'slab_02', 'slab_03']
```

Check results:
```bash
verdi process show 11739
```

Outputs:
- `generated_slabs`: Original slab structures (4 terminations)
- `relaxed_slabs`: Relaxed structures after VASP optimization
- `slab_energies`: Total energies for each slab

## Advantages of Scatter-Gather Pattern

1. **Natural Python**: Uses familiar Python patterns (loops, dictionaries, function calls)
2. **Flexible Control Flow**: Can use conditionals, nested loops, etc.
3. **Reusable Components**: Graph tasks can be composed and reused
4. **Type Safety**: Full Python type annotations with IDE support
5. **Debugging**: Standard Python debugging tools work

## Comparison with Zone Approach

The zone approach (`test_modules/zone_approach`) uses the `Map` context manager:

```python
with Map(slab_namespace) as map_zone:
    relaxation = relax_task_cls(structure=map_zone.item.value, ...)
    energy = extract_total_energy(energies=relaxation.misc).result
    map_zone.gather({'relaxed_slabs': relaxation.structure, 'slab_energies': energy})
```

The scatter-gather approach uses plain Python:

```python
@task.graph
def relax_slabs_scatter(slabs, ...):
    relaxed = {}
    energies = {}
    for label, structure in slabs.items():
        relaxation = relax_task_cls(structure=structure, ...)
        relaxed[label] = relaxation.structure
        energies[label] = extract_total_energy(relaxation.misc).result
    return {'relaxed_structures': relaxed, 'energies': energies}
```

Both approaches achieve the same parallelization, but the scatter-gather pattern is more "pythonic" and flexible.

## Tips and Tricks

1. **Load external resources early**: Load codes, structures, etc. outside nested `@task.graph` functions
2. **Use orm types in calcfunctions**: All inputs/outputs of `@task.calcfunction` should be AiiDA Data types
3. **Dynamic namespaces**: Use `dynamic(Type)` annotation for collections with unknown size
4. **Parallel execution is automatic**: Tasks without dependencies run in parallel automatically
5. **Accessing outputs**: Remember that namespace annotations expose outputs directly, not via `.result`

## Documentation References

- [AiiDA-WorkGraph Scatter-Gather Guide](https://aiida-workgraph.readthedocs.io/en/latest/howto/run_tasks_in_parallel.html)
- [AiiDA-WorkGraph Task Graphs](https://aiida-workgraph.readthedocs.io/en/latest/concept/task.html#graph-task)
- [Dynamic Namespaces](https://aiida-workgraph.readthedocs.io/en/latest/concept/socket.html#dynamic-namespaces)
