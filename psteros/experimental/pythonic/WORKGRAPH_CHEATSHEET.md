# AiiDA-WorkGraph Cheat Sheet

Quick reference for building workflows with AiiDA-WorkGraph's pythonic approach, based on the scatter-gather implementation in PS-TEROS.

---

## Task Types

### `@task` - Simple Python Function
```python
@task
def add(x: int, y: int) -> int:
    return x + y
```
- No provenance for the operation itself
- Fast, lightweight
- Use for temporary calculations

### `@task.calcfunction` - AiiDA Calculation Function
```python
@task.calcfunction
def extract_energy(data: orm.Dict) -> orm.Float:
    return orm.Float(data['energy'])
```
- **All inputs/outputs MUST be AiiDA types** (`orm.Float`, `orm.StructureData`, etc.)
- Full provenance tracking
- Creates CalcFunctionNode in database
- Use when you need to track this specific operation

### `@task.graph` - Workflow Builder
```python
@task.graph
def parallel_workflow(
    inputs: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)]
) -> t.Annotated[dict, namespace(results=dynamic(orm.Float))]:
    results = {}
    for key, value in inputs.items():
        results[key] = process(value).result
    return {'results': results}
```
- **Builds subworkflows dynamically**
- Use for scatter-gather patterns
- Handles "future" values (outputs from other tasks)
- Python control flow evaluated at graph build time

### `task(WorkChain)` - Wrap Existing WorkChain
```python
from aiida.plugins import WorkflowFactory

VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
vasp_task = task(VaspWorkChain)

# Use it
@task.graph
def run_vasp(structure):
    return vasp_task(structure=structure, code=code, ...).structure
```
- No need for builders
- Pass parameters directly as kwargs

---

## Type Annotations

### Basic Types
```python
from aiida import orm

x: orm.Int          # Single integer
y: orm.Float        # Single float
s: orm.StructureData  # Single structure
```

### Dynamic Namespaces (Unknown Number of Items)
```python
import typing as t
from aiida_workgraph import dynamic, namespace

# Input: dict with unknown keys
inputs: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)]

# Output: namespace with unknown keys  
-> t.Annotated[dict, namespace(results=dynamic(orm.Float))]
```

### Named Namespace Outputs
```python
-> t.Annotated[dict, namespace(
    energy=orm.Float,
    structure=orm.StructureData,
    forces=orm.ArrayData
)]
```

---

## Scatter-Gather Pattern

**Problem**: Process N items in parallel, where N is unknown at workflow definition time.

**Solution**: Use `@task.graph` with Python loop.

```python
@task.calcfunction
def process_item(item: orm.Float) -> orm.Float:
    return orm.Float(item.value * 2)

@task.graph
def parallel_processing(
    items: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)]
) -> t.Annotated[dict, namespace(results=dynamic(orm.Float))]:
    """
    Scatter-gather: process all items in parallel.
    Each iteration creates an independent task.
    """
    results = {}
    
    # Scatter: create tasks (no dependencies between iterations)
    for key, item in items.items():
        results[key] = process_item(item).result
    
    # Gather: return all results
    return {'results': results}
```

**Why it works**: Each `process_item` call is independent → WorkGraph executes them in parallel automatically.

---

## Building and Running Workflows

### Method 1: Build then Run
```python
# Build the graph
workflow = my_workflow.build(param1=value1, param2=value2)

# Inspect before running
workflow.to_html()  # Visualize

# Run
workflow.run()

# Access results
results = workflow.process.outputs
```

### Method 2: Submit (Non-blocking)
```python
workflow = my_workflow.build(...)
workflow.submit()  # Returns immediately

# Check later
from aiida import orm
node = orm.load_node(workflow.pk)
if node.is_finished_ok:
    results = node.outputs
```

---

## Accessing Outputs

### During Graph Building (Before .run())
```python
@task.graph
def workflow():
    result1 = task1(x=10).result
    result2 = task2(y=result1).output  # Use .result or specific output name
    return result2
```

### After Running (In Script)
```python
workflow.run()
outputs = workflow.process.outputs

# Two ways to access:
value1 = outputs.key_name        # Attribute access (cleaner)
value2 = outputs['key_name']     # Dict access (also works)

# For dynamic namespaces, iterate:
for key in outputs._sockets:
    if not key.startswith('_'):
        item = getattr(outputs, key)
```

### Check if Outputs Exist
```python
if hasattr(outputs, '_sockets'):
    # It's a TaskSocketNamespace (during build)
    keys = [k for k in outputs._sockets if not k.startswith('_')]
else:
    # It's an AttributeDict (after run)
    keys = [k for k in outputs.keys() if not k.startswith('_')]
```

---

## Common Patterns

### Pattern 1: Generate Dynamic Collection
```python
@task.calcfunction
def generate_items(n: orm.Int) -> t.Annotated[
    dict, namespace(items=dynamic(orm.Int))
]:
    """Generate n items."""
    return {
        'items': {
            f'item_{i}': orm.Int(i)
            for i in range(n.value)
        }
    }
```

### Pattern 2: Process Collection in Parallel
```python
@task.graph
def process_all(
    items: t.Annotated[dict[str, orm.Int], dynamic(orm.Int)]
) -> t.Annotated[dict, namespace(processed=dynamic(orm.Int))]:
    """Process each item independently."""
    processed = {}
    for key, item in items.items():
        processed[key] = process_one(item).result
    return {'processed': processed}
```

### Pattern 3: Conditional Workflow Steps
```python
@task.graph
def conditional_workflow(data, do_extra_step=False):
    result = step1(data).output
    
    if do_extra_step:
        result = step2(result).output  # Only added if True
    
    return result
```

**Note**: Condition evaluated at graph build time, not execution time.

### Pattern 4: Nested Graphs
```python
@task.graph
def inner_workflow(item):
    return process(item).result

@task.graph
def outer_workflow(items):
    results = {}
    for key, item in items.items():
        results[key] = inner_workflow(item)  # Nested graph
    return {'results': results}
```

### Pattern 5: Combining Results
```python
@task.calcfunction
def sum_values(
    values: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)]
) -> orm.Float:
    """Gather: combine all parallel results."""
    total = sum(v.value for v in values.values())
    return orm.Float(total)

@task.graph
def scatter_gather_sum(items):
    # Scatter
    squared = {}
    for key, value in items.items():
        squared[key] = square(value).result
    
    # Gather
    total = sum_values(squared).result
    return total
```

---

## Best Practices

### ✅ DO

1. **Use dynamic namespaces for collections**
   ```python
   items: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)]
   ```

2. **Keep calcfunctions pure (no side effects)**
   ```python
   @task.calcfunction
   def calc(x: orm.Float) -> orm.Float:
       return orm.Float(x.value * 2)  # Pure function
   ```

3. **Use descriptive names**
   ```python
   relax_slabs_scatter  # Good
   process_data         # Too generic
   ```

4. **Check workflow status before accessing results**
   ```python
   workflow.run()
   if workflow.process.is_finished_ok:
       results = workflow.process.outputs
   ```

5. **Use type annotations everywhere**
   ```python
   def func(x: orm.Float) -> orm.Float:  # Clear types
   ```

### ❌ DON'T

1. **Don't use Python types in calcfunctions**
   ```python
   @task.calcfunction
   def bad(x: int) -> int:  # ❌ Won't work!
       return x * 2
   ```

2. **Don't try to iterate at execution time in @task.graph**
   ```python
   @task.graph
   def bad(items):
       # This won't work - items is a future, not actual data
       for i in range(len(items)):  # ❌ Can't get len of future
           ...
   ```

3. **Don't modify inputs in calcfunctions**
   ```python
   @task.calcfunction
   def bad(data: orm.Dict) -> orm.Dict:
       data['new_key'] = 123  # ❌ Don't mutate inputs!
       return data
   ```

4. **Don't forget to return dictionaries for namespace outputs**
   ```python
   @task.graph
   def bad(x):
       result = process(x).result
       return result  # ❌ Should return {'result': result}
   ```

5. **Don't use .value on future outputs**
   ```python
   @task.graph
   def bad(x):
       result = process(x).result
       value = result.value  # ❌ Can't access .value of future!
       return value
   ```

---

## Debugging Tips

### Visualize the Graph
```python
workflow = my_workflow.build(...)
workflow.to_html()  # Opens in browser
```

### Check Process Status
```bash
verdi process show <PK>
verdi process report <PK>
```

### Get Failed Job Logs
```bash
verdi process report <PK> | grep -i error
```

### Test with Mock Mode
```python
# Add a mock mode to your workflow
def build_workflow(use_mock=False):
    if use_mock:
        return build_mock_workflow()
    else:
        return build_real_workflow()
```

### Print During Graph Building
```python
@task.graph
def debug_workflow(items):
    print(f"Building graph with {len(items)} items")  # This prints during .build()
    results = {}
    for key, item in items.items():
        print(f"Adding task for {key}")  # This also prints during .build()
        results[key] = process(item).result
    return {'results': results}
```

---

## Provenance Inspection

### View Workflow Tree
```bash
verdi node graph generate <PK>  # Generate graph image
```

### Query Outputs
```python
from aiida import orm

workflow = orm.load_node(<PK>)

# Get outputs
outputs = workflow.outputs

# Get specific output
energy = workflow.outputs.energy

# Check called processes
for child in workflow.called:
    print(f"{child.process_label}: {child.pk}")
```

### Find All Workflows of Type
```bash
verdi process list -S FINISHED -p WorkGraph
```

---

## Real-World Example: Slab Relaxation

From the PS-TEROS implementation:

```python
@task.graph
def slab_relaxation_scatter_gather(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    code: orm.Code,
    compute_thermodynamics: bool = False,
    bulk_energy: orm.Float | None = None,
    reference_energies: orm.Dict | None = None,
    ...
):
    """
    Complete workflow:
    1. Generate slab terminations (unknown number)
    2. Relax all slabs in parallel with VASP
    3. Extract energies
    4. Optionally compute thermodynamics in parallel
    """
    
    # Step 1: Generate slabs
    slabs = generate_slab_structures(
        structure=bulk_structure,
        miller_indices=miller_indices,
        ...
    ).slabs
    
    # Step 2: Relax all slabs in parallel
    relaxation = relax_slabs_scatter(
        slabs=slabs,
        code=code,
        ...
    )
    
    output = {
        'generated_slabs': slabs,
        'relaxed_slabs': relaxation.relaxed_structures,
        'slab_energies': relaxation.energies,
    }
    
    # Step 3: Optional thermodynamics
    if compute_thermodynamics:
        thermo = compute_surface_energies_scatter(
            slabs=relaxation.relaxed_structures,
            energies=relaxation.energies,
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            reference_energies=reference_energies,
            ...
        )
        output['surface_energies'] = thermo.surface_energies
    
    return output
```

**Key points:**
- Dynamic number of slabs handled automatically
- Two levels of parallelization (relaxations + thermodynamics)
- Conditional step using Python `if`
- Full type annotations with dynamic namespaces
- Clean, readable code using Python patterns

---

## Quick Command Reference

### Workflow Management
```bash
verdi daemon restart              # Restart after code changes
verdi process list                # List running processes
verdi process show <PK>           # Show process details
verdi process report <PK>         # Show process report
verdi process kill <PK>           # Kill running process
```

### Node Inspection
```bash
verdi node show <PK>              # Show node details
verdi node graph generate <PK>    # Generate provenance graph
verdi data structure show <PK>    # Show structure
```

### Database
```bash
verdi profile list                # List profiles
verdi profile set-default <NAME>  # Set default profile
verdi status                      # Check AiiDA status
```

---

## Additional Resources

- **Official Docs**: https://aiida-workgraph.readthedocs.io/
- **Scatter-Gather Guide**: https://aiida-workgraph.readthedocs.io/en/latest/howto/run_tasks_in_parallel.html
- **Task Types**: https://aiida-workgraph.readthedocs.io/en/latest/concept/task.html
- **Dynamic Namespaces**: https://aiida-workgraph.readthedocs.io/en/latest/concept/socket.html#dynamic-namespaces
- **PS-TEROS Implementation**: See `README.md`, `THERMODYNAMICS_COMPLETE.md`, and `AIAT_IMPLEMENTATION.md` in this directory

---

**Last Updated**: 2025-10-06
**Version**: 1.0 (Based on aiida-workgraph 1.0+)
