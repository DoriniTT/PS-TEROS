# max_number_jobs Investigation

## Problem

The `max_number_jobs` parameter on the main WorkGraph doesn't limit concurrent VASP calculations because nested sub-workgraphs created by `@task.graph` decorators spawn their own independent calculations.

**Current behavior:**
```
Main WorkGraph (max_number_jobs=2)
  └─> @task.graph function creates sub-workgraph
       └─> Spawns 10+ VASP calculations (ignores parent's max_number_jobs)
```

## Hypothesis

There may be a way to pass `max_number_jobs` to the sub-workgraphs created by `@task.graph` decorators, so they respect the concurrency limit.

## Investigation Approaches

### Approach 1: Parameter Passing
Add `max_number_jobs` as a parameter to `@task.graph` decorated functions, then set it on the WorkGraph instance created inside.

**Question:** How do we access the WorkGraph instance inside `@task.graph`?

### Approach 2: Decorator Configuration
Try passing `max_number_jobs` directly to the `@task.graph` decorator.

**Example:** `@task.graph(max_number_jobs=2)`

### Approach 3: Framework Properties
Check if WorkGraph has global/cascading concurrency control that automatically propagates to nested graphs.

## Test Files

- `test_approach_1_parameter.py` - Test parameter passing
- `test_approach_2_decorator.py` - Test decorator configuration
- `test_approach_3_framework.py` - Test framework properties
- `run_all_tests.py` - Run all tests and compare results

## Expected Outcome

Find a simple way to make nested workgraphs respect `max_number_jobs` without restructuring the entire codebase.
