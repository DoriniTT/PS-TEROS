# max_number_jobs Investigation - Findings

## Summary

We investigated three approaches to make nested workgraphs (created by `@task.graph`) respect `max_number_jobs` concurrency limits.

### Results

| Approach | Can Pass Parameter? | Does it Work? | Status |
|----------|---------------------|---------------|--------|
| 1. Parameter Passing | ✅ YES | ❌ NOT YET | Partial Success |
| 2. Decorator Configuration | ❌ NO | ❌ NO | Not Supported |
| 3. Framework Cascade | N/A | ❌ NO | Not Supported |

## Detailed Findings

### Approach 1: Parameter Passing ✅ (Partial Success)

**Hypothesis**: Add `max_number_jobs` as a parameter to `@task.graph` functions.

**Test**: Added `max_number_jobs: int = None` parameter to decorated function.

**Result**:
- ✅ **Parameter CAN be passed** to @task.graph functions
- ✅ **Parameter APPEARS as input** to the nested workgraph (verified with `verdi process show`)
- ❌ **Parameter is NOT automatically used** by the nested workgraph
- ❌ **Cannot access WorkGraph instance** inside @task.graph to set it manually

**Evidence**:
```
# From verdi process show 122259
Inputs                      PK      Type
--------------------------  ------  -------------
graph_inputs
    max_number_jobs         122258  Int  # ← VALUE WAS PASSED!
```

```python
# Checking the value
>>> node = orm.load_node(122258)
>>> node.value
2  # ← Correct value was passed
```

**Timing Analysis**:
```
# From verdi process report 122259
18:45:36 - All 5 tasks ready
18:45:39 - 4 tasks ready (1 completed) # +3 sec
18:45:42 - 3 tasks ready (2 completed) # +3 sec
18:45:44 - 2 tasks ready (3 completed) # +2 sec
18:45:47 - 1 task ready (4 completed)  # +3 sec
18:45:49 - 0 tasks ready (5 completed) # +2 sec
```

Tasks ran **one at a time**, not 2 at a time as `max_number_jobs=2` would require.

**Attempted Fix**:
Tried to access WorkGraph context to set `max_number_jobs` programmatically:
```python
from aiida_workgraph import get_workgraph_context  # ← Doesn't exist!
wg = get_workgraph_context()
wg.max_number_jobs = max_number_jobs
```

**Error**: `ImportError: cannot import name 'get_workgraph_context'`

**Conclusion**:
- Parameter passing WORKS for getting the value into the nested workgraph
- But we need to find the right API to actually SET `max_number_jobs` on the WorkGraph instance created by `@task.graph`

---

### Approach 2: Decorator Configuration ❌ (Not Supported)

**Hypothesis**: Pass `max_number_jobs` directly to the `@task.graph` decorator.

**Test**:
```python
@task.graph(max_number_jobs=2)  # ← Try this syntax
def process_structures(...):
    ...
```

**Result**:
```
TypeError: TaskDecoratorCollection.decorator_graph() got an unexpected keyword argument 'max_number_jobs'
```

**Decorator Inspection**:
```python
>>> import inspect
>>> from aiida_workgraph import task
>>> inspect.signature(task.graph)
(*args, **kwargs)  # ← No documented parameters
```

**Conclusion**: The `@task.graph` decorator does NOT accept configuration parameters.

---

### Approach 3: Framework Cascade ❌ (Not Supported)

**Hypothesis**: Setting `max_number_jobs` on the parent workgraph automatically propagates to nested workgraphs.

**Test**:
```python
wg = WorkGraph("main")
wg.max_number_jobs = 2  # ← Set on parent

# Call @task.graph function (creates nested workgraph)
result = process_structures(...)
```

**Result**: Nested workgraph did NOT inherit `max_number_jobs` from parent.

**Framework Inspection**:
```python
>>> wg = WorkGraph("test")
>>> wg.max_number_jobs = 2

# Check for cascade-related attributes
>>> hasattr(wg, 'cascade_limits')
False
>>> hasattr(wg, 'inherit_limits')
False
>>> hasattr(wg, 'propagate_limits')
False
>>> hasattr(wg, 'global_max_jobs')
False
```

**Documentation Check**:
> "The `max_number_jobs` attribute only governs child processes created by **this specific** WorkGraph instance."

The word "specific" indicates no cascading behavior.

**Conclusion**: WorkGraph does NOT have cascade/inheritance of `max_number_jobs`.

---

## Key Discovery

**WE CAN PASS `max_number_jobs` AS A PARAMETER**, but we need to find how to **USE it** inside the `@task.graph` function.

The nested workgraph IS created with the parameter as an input, but it doesn't automatically apply it.

---

## Next Steps to Investigate

### 1. Find the WorkGraph Instance API

Within a `@task.graph` function, we need to find how to:
- Access the WorkGraph instance being created
- Set `max_number_jobs` on that instance

**Possible approaches**:
- Check if there's a `self` or `this` reference
- Look for WorkGraph builder patterns in source code
- Check if WorkGraph() can be instantiated explicitly inside @task.graph
- Look for aiida_workgraph context managers

### 2. Alternative: Return WorkGraph Directly

Instead of using `@task.graph`, maybe we can:
```python
def create_slab_calculations(..., max_number_jobs=None):
    wg = WorkGraph("slab_calcs")
    wg.max_number_jobs = max_number_jobs  # ← Set it explicitly

    # Add tasks
    for slab in slabs:
        wg.add_task(VaspWorkChain, ...)

    return wg
```

Then call it as a regular task (not @task.graph).

### 3. Check AiiDA WorkGraph Source Code

Look at:
- `aiida_workgraph/tasks/graph_task.py` - How @task.graph creates WorkGraphs
- `aiida_workgraph/engine/workgraph.py` - WorkGraph execution engine
- Find where `max_number_jobs` is actually used

### 4. Ask AiiDA WorkGraph Developers

This might be a feature request:
- "Allow @task.graph to accept max_number_jobs parameter"
- "Provide API to set WorkGraph properties from within @task.graph"

---

## Test Files Created

All test files are in: `/home/thiagotd/git/PS-TEROS/teros/experimental/max_jobs_investigation/`

- `workgraph_functions.py` - @task.graph functions for testing
- `mock_tasks.py` - Mock VASP calculations
- `test_approach_1_parameter.py` - Test parameter passing
- `test_approach_2_decorator.py` - Test decorator configuration
- `test_approach_3_framework.py` - Test framework cascade
- `run_investigation.py` - Run all tests

### Example Test Results

**Approach 1** - WorkGraph PK: 122252
- Nested WorkGraph: 122259
- max_number_jobs input: 122258 (value=2)
- 5 VASP tasks created
- Tasks ran sequentially (not respecting max_number_jobs=2)

**Approach 3** - WorkGraph PK: 122275
- Similar results

---

## Recommendation

**Short term**: We confirmed that parameters CAN be passed to @task.graph functions. We just need to find the right API to USE them.

**Next action**: Investigate aiida_workgraph source code or contact developers to find how to access/set the WorkGraph instance properties from within @task.graph.

**Alternative**: Consider the two-stage workflow approach from `surface_thermo_preset_serial` as a working solution while we investigate the proper API.

---

## Questions for AiiDA WorkGraph Developers

1. How do we access the WorkGraph instance created by `@task.graph` from within the decorated function?
2. Is there an API to set `max_number_jobs` programmatically inside `@task.graph`?
3. Should `max_number_jobs` be supported as a standard parameter for `@task.graph` functions?
4. Is there a recommended pattern for controlling concurrency in nested workgraphs?

---

**Date**: 2025-11-02
**Status**: Investigation Complete - Partial Success
**Follow-up**: Need to find WorkGraph instance API or request feature from developers
