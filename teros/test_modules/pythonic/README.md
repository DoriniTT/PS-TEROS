# Pythonic Scatter-Gather Pattern: Complete Implementation Guide

This directory demonstrates the **scatter-gather pattern** for parallelizing computational workflows using AiiDA-WorkGraph's pythonic approach. It is applied to slab surface generation and parallel VASP relaxation calculations.

## Table of Contents

1. [Overview](#overview)
2. [Core Concepts](#core-concepts)
3. [Architecture](#architecture)
4. [Implementation Details](#implementation-details)
5. [Slab Relaxation Use Case](#slab-relaxation-use-case)
6. [Task Connections and Data Flow](#task-connections-and-data-flow)
7. [Running the Workflow](#running-the-workflow)
8. [Best Practices](#best-practices)
9. [Troubleshooting](#troubleshooting)

---

## Overview

The **scatter-gather pattern** is a workflow design pattern for parallel execution:

1. **Scatter Phase**: Distribute a collection of inputs to multiple independent tasks that execute in parallel
2. **Gather Phase**: Collect and aggregate results from all parallel tasks

This implementation uses AiiDA-WorkGraph's `@task.graph` decorator, which provides a **pythonic** alternative to context managers like `Map`, allowing you to use plain Python control flow (loops, conditionals, etc.) to build complex workflows.

### Key Benefits

- **Natural Python Syntax**: Use familiar Python patterns instead of learning DSL
- **Full Flexibility**: Conditionals, nested loops, early returns, etc.
- **Type Safety**: Complete Python type annotations with IDE support
- **Automatic Parallelization**: Tasks without dependencies run concurrently
- **Composable**: Graph tasks can be nested and reused
- **Full Provenance**: Complete AiiDA provenance tracking

### Comparison: Scatter-Gather vs Map-Zone

| Feature | Pythonic Scatter-Gather | Map-Zone |
|---------|------------------------|----------|
| **Pattern** | `@task.graph` with Python loops | `with Map(...)` context manager |
| **Code Style** | Plain Python | DSL-style context manager |
| **Flexibility** | High - full Python power | Limited to Map iteration |
| **Conditionals** | Native Python `if`/`else` | Requires `If` zone |
| **Nested Loops** | Native Python | Requires nested `Map` |
| **Provenance** | Tasks nested in graph | Flat workflow structure |
| **Learning Curve** | Familiar Python patterns | New DSL to learn |
| **When to Use** | Complex logic, reusable components | Simple iterations |

## Files

- **`slabs_relax.py`**: CLI driver that launches the workflow (supports `--mock` mode for testing)
- **`workgraph.py`**: Complete workflow implementation using scatter-gather pattern
- **`README.md`**: This comprehensive guide
---

## Core Concepts

### 1. Tasks: The Building Blocks

A **task** is the fundamental unit of work in AiiDA-WorkGraph. There are several types:

#### a) `@task` - Python Function Task

Basic Python function wrapped as a task:

```python
@task
def add(x: int, y: int) -> int:
    """Simple addition - no AiiDA provenance for intermediates."""
    return x + y
```

**When to use**: Lightweight operations, temporary calculations, or when you don't need provenance for the specific operation.

#### b) `@task.calcfunction` - AiiDA Calculation Function

Python function that stores full provenance in AiiDA database:

```python
@task.calcfunction
def extract_energy(results: orm.Dict) -> orm.Float:
    """Stores this operation in provenance graph."""
    return orm.Float(results['energy'])
```

**When to use**: 
- ✅ When inputs/outputs are AiiDA Data types (`orm.Float`, `orm.StructureData`, etc.)
- ✅ When you need provenance tracking for this specific operation
- ✅ When working with dynamic namespaces that return AiiDA types

**Critical Rule**: All inputs and outputs MUST be AiiDA Data types (`orm.*`). Regular Python types will cause errors.

#### c) `@task.graph` - Graph Task (Workflow Builder)

A task that builds and returns a WorkGraph dynamically:

```python
@task.graph
def parallel_calculations(
    inputs: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
) -> t.Annotated[dict, namespace(results=dynamic(orm.Float))]:
    """Builds a graph at runtime - enables scatter-gather pattern."""
    results = {}
    for key, value in inputs.items():
        results[key] = process_value(value).result
    return {'results': results}
```

**When to use**:
- ✅ **Scatter-gather patterns** (iterate over dynamic collections)
- ✅ When you need to build subgraphs based on runtime data
- ✅ When inputs are "future values" (outputs from other tasks)
- ✅ Complex conditional logic that would be awkward with zones

**Key Property**: The loop/logic runs at **graph construction time**, not at data processing time. This is why it can handle future values.

#### d) `task(WorkChain)` - WorkChain Wrapper

Wraps an existing AiiDA WorkChain as a task:

```python
from aiida.plugins import WorkflowFactory

vasp_wc = WorkflowFactory('vasp.v2.vasp')
relax_task = task(vasp_wc)

# Use it like any task
result = relax_task(structure=slab, code=code, ...)
```

**When to use**: Integrating existing AiiDA WorkChains into WorkGraph workflows.

### 2. WorkGraph: The Container

A **WorkGraph** is a container that holds tasks and defines their connections:

```python
from aiida_workgraph import WorkGraph

# Method 1: Using context manager
with WorkGraph(name='my_workflow') as wg:
    task1 = add(x=1, y=2)
    task2 = multiply(value=task1.result, factor=3)

# Method 2: Using @task.graph (recommended for scatter-gather)
@task.graph
def my_workflow(n: int) -> int:
    task1 = add(x=1, y=n)
    task2 = multiply(value=task1.result, factor=3)
    return task2.result
```

**WorkGraph responsibilities**:
- Manages task execution order based on dependencies
- Handles parallel execution (tasks without dependencies run concurrently)
- Stores workflow provenance in AiiDA database
- Provides input/output management

### 3. Sockets and Data Flow

**Sockets** are the connection points between tasks:

```python
@task
def add(x: int, y: int) -> int:
    return x + y

# When you call the task in a graph:
result_task = add(x=5, y=3)

# Access outputs via sockets:
value = result_task.result  # This is a "future value" socket
```

#### Socket Types

1. **Regular Result Socket**: Single output
   ```python
   @task
   def func() -> int:
       return 42
   
   # Access: task.result
   ```

2. **Named Namespace Socket**: Multiple named outputs
   ```python
   @task.calcfunction
   def func() -> t.Annotated[dict, namespace(a=int, b=int)]:
       return {'a': 1, 'b': 2}
   
   # Access: task.a, task.b (NOT task.result.a)
   ```

3. **Dynamic Namespace Socket**: Variable number of outputs
   ```python
   @task.calcfunction
   def generate_slabs(...) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
       return {'slabs': {'slab_00': struct1, 'slab_01': struct2, ...}}
   
   # Access: task.slabs (which is a namespace with slab_00, slab_01, etc.)
   ```

### 4. Future Values

A **future value** is a socket that represents data that will be available after a task executes:

```python
@task.graph
def workflow():
    # At this point, 'result' is a "future value" - data doesn't exist yet
    result = expensive_calculation(x=100).result
    
    # You CANNOT do this (will raise error):
    # if result > 50:  # ❌ Can't evaluate future value
    
    # But you CAN do this (pass future to another task):
    next_result = process_value(value=result).result  # ✅ OK
    
    # You CAN iterate over future namespaces in @task.graph:
    for key, value in result.items():  # ✅ OK in @task.graph
        process(value)
```

**Key Insight**: `@task.graph` functions can work with future values because they build the graph structure, they don't process the actual data.

---

## Architecture

### Workflow Structure Hierarchy

```
WorkGraph (Top Level)
├── @task.graph: slab_relaxation_scatter_gather
│   ├── Loads external resources (code)
│   ├── @task.calcfunction: generate_slab_structures
│   │   └── Returns: dynamic namespace of slabs
│   └── @task.graph: relax_slabs_scatter (nested graph)
│       ├── Loop over slabs (scatter phase)
│       ├── task(WorkChain): VaspWorkChain (parallel instances)
│       │   └── Runs VASP calculations
│       ├── @task.calcfunction: extract_total_energy (parallel instances)
│       └── Returns: collected results (gather phase)
└── Outputs: generated_slabs, relaxed_slabs, slab_energies
```

### Data Flow Diagram

```
Input: bulk_structure
         ↓
┌─────────────────────────────────────────────────────────┐
│ generate_slab_structures (@task.calcfunction)          │
│ - Generates N slab terminations                        │
└─────────────────────────────────────────────────────────┘
         ↓
    Dynamic Namespace: {slab_00: struct, slab_01: struct, ...}
         ↓
┌─────────────────────────────────────────────────────────┐
│ relax_slabs_scatter (@task.graph)                      │
│                                                          │
│  ╔═══════════════════════════════════════════════════╗ │
│  ║ SCATTER PHASE (loop creates N tasks in parallel) ║ │
│  ╚═══════════════════════════════════════════════════╝ │
│                                                          │
│  for each slab:                                         │
│    ┌──────────────────────────┐                        │
│    │ VaspWorkChain            │ ← Parallel Instance 1  │
│    └──────────────────────────┘                        │
│            ↓                                            │
│    ┌──────────────────────────┐                        │
│    │ extract_total_energy     │                        │
│    └──────────────────────────┘                        │
│                                                          │
│    ┌──────────────────────────┐                        │
│    │ VaspWorkChain            │ ← Parallel Instance 2  │
│    └──────────────────────────┘                        │
│            ↓                                            │
│    ┌──────────────────────────┐                        │
│    │ extract_total_energy     │                        │
│    └──────────────────────────┘                        │
│                                                          │
│    ... (N instances total)                              │
│                                                          │
│  ╔═══════════════════════════════════════════════════╗ │
│  ║ GATHER PHASE (collect all results)               ║ │
│  ╚═══════════════════════════════════════════════════╝ │
│                                                          │
│  Returns: {relaxed_structures: {...},                  │
│            energies: {...}}                             │
└─────────────────────────────────────────────────────────┘
         ↓
    Outputs: relaxed_slabs, slab_energies
```

---

## Implementation Details


---

## Implementation Details

### Component 1: Slab Generation

**Function**: `generate_slab_structures`

```python
@task.calcfunction
def generate_slab_structures(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    min_slab_thickness: orm.Float,
    min_vacuum_thickness: orm.Float,
    lll_reduce: orm.Bool,
    center_slab: orm.Bool,
    symmetrize: orm.Bool,
    primitive: orm.Bool,
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """Generate slab terminations from bulk structure."""
    # Implementation uses pymatgen SlabGenerator
    # ...
    return {'slabs': slab_nodes}
```

**Key Points**:

1. **Decorator**: Uses `@task.calcfunction` because:
   - Returns AiiDA Data types (`orm.StructureData`)
   - We want provenance for slab generation
   - Dynamic namespace requires proper AiiDA node handling

2. **Inputs**: ALL inputs are AiiDA Data types (e.g., `orm.Float` not `float`)
   - This is a requirement for `@task.calcfunction`
   - Allows AiiDA to track provenance

3. **Return Annotation**: `t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]`
   - `namespace(...)`: Creates a named namespace output
   - `slabs=dynamic(...)`: The namespace contains a dynamic collection
   - `dynamic(orm.StructureData)`: Unknown number of `StructureData` nodes
   - Results in output socket: `task.slabs` (NOT `task.result.slabs`)

4. **Return Value**: Dictionary with key matching annotation (`'slabs'`)

### Component 2: Energy Extraction

**Function**: `extract_total_energy`

```python
@task.calcfunction
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """Extract total energy from VASP output dictionary."""
    energy_dict = energies.get_dict()
    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']
    
    # Try multiple keys for compatibility
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])
    
    raise ValueError(f'Unable to find total energy')
```

**Key Points**:

1. **Simple extraction logic** with fallbacks for different VASP output formats
2. **Returns `orm.Float`** for provenance tracking
3. **Will be called N times** (once per slab) in parallel

### Component 3: Scatter-Gather Orchestrator

**Function**: `relax_slabs_scatter`

```python
@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float)
)]:
    """Scatter-gather: relax all slabs in parallel."""
    
    # Load WorkChain class
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_wc)
    
    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}

    # ═══════════════════════════════════════
    # SCATTER PHASE: Create parallel tasks
    # ═══════════════════════════════════════
    for label, structure in slabs.items():
        # Create VASP relaxation task
        relaxation = relax_task_cls(
            structure=structure,
            code=code,
            parameters={'incar': dict(parameters)},
            options=dict(options),
            potential_family=potential_family,
            potential_mapping=dict(potential_mapping),
            clean_workdir=clean_workdir,
            kpoints_spacing=kpoints_spacing,
        )
        
        # Extract energy from VASP output
        energy = extract_total_energy(energies=relaxation.misc).result
        
        # Store future values
        relaxed[label] = relaxation.structure
        energies_ns[label] = energy

    # ═══════════════════════════════════════
    # GATHER PHASE: Return collected results
    # ═══════════════════════════════════════
    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
    }
```

**Key Points**:

1. **`@task.graph` Decorator**: Required because we're:
   - Iterating over a future value (`slabs` doesn't exist yet at graph-build time)
   - Creating multiple tasks dynamically based on runtime data
   - Building a subgraph that will be executed later

2. **Scatter Phase** (the loop):
   - Each iteration creates **two tasks**: `relax_task_cls` and `extract_total_energy`
   - These tasks have **no dependencies on each other** (between slabs)
   - AiiDA automatically executes them **in parallel**
   - We're collecting **future values** (sockets), not actual data

3. **Gather Phase** (the return):
   - We return dictionaries of future values
   - The annotation tells AiiDA these are dynamic namespaces
   - Results will be available after all tasks complete

4. **Why This Works**:
   ```python
   # At graph-build time:
   for label, structure in slabs.items():  # ← This iterates over socket names
       relaxation = relax_task_cls(...)    # ← Creates task nodes
       relaxed[label] = relaxation.structure  # ← Stores socket references
   
   # At runtime (later):
   # AiiDA executes all VaspWorkChain tasks in parallel
   # Then executes all extract_total_energy tasks in parallel
   # Finally populates the output sockets with actual data
   ```

5. **Code Loading**: `code = orm.load_code(code_label)` happens in **parent graph**
   - Can't load within this function (would try to load from future value)
   - Must be loaded where `code_label` is a concrete string

### Component 4: Top-Level Workflow

**Function**: `slab_relaxation_scatter_gather`

```python
@task.graph
def slab_relaxation_scatter_gather(
    bulk_structure: orm.StructureData,
    *,
    miller_indices: tuple[int, int, int] | list[int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    code_label: str,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    generated_slabs=dynamic(orm.StructureData),
    relaxed_slabs=dynamic(orm.StructureData),
    slab_energies=dynamic(orm.Float),
)]:
    """Top-level workflow orchestrating slab generation and relaxation."""
    
    # ═══════════════════════════════════════
    # SETUP: Load resources that need concrete values
    # ═══════════════════════════════════════
    code = orm.load_code(code_label)  # ← Must happen here!
    
    if isinstance(miller_indices, tuple):
        miller_indices = list(miller_indices)

    # ═══════════════════════════════════════
    # STEP 1: Generate slabs
    # ═══════════════════════════════════════
    slab_namespace = generate_slab_structures(
        bulk_structure=bulk_structure,
        miller_indices=orm.List(list=miller_indices),
        min_slab_thickness=orm.Float(min_slab_thickness),
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        symmetrize=orm.Bool(symmetrize),
        primitive=orm.Bool(primitive),
    ).slabs  # ← Access namespace directly (no .result)

    # ═══════════════════════════════════════
    # STEP 2: Scatter-gather relaxations
    # ═══════════════════════════════════════
    relaxation_outputs = relax_slabs_scatter(
        slabs=slab_namespace,  # ← Pass future namespace
        code=code,  # ← Pass concrete code node
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
    )

    # ═══════════════════════════════════════
    # STEP 3: Return all outputs
    # ═══════════════════════════════════════
    return {
        'generated_slabs': slab_namespace,
        'relaxed_slabs': relaxation_outputs.relaxed_structures,
        'slab_energies': relaxation_outputs.energies,
    }
```

**Key Points**:

1. **Setup Section**:
   - Loads code with `orm.load_code()` using the string label
   - This MUST happen here because `code_label` is a concrete value
   - Inside nested `@task.graph`, it would be a future value

2. **Step 1 - Slab Generation**:
   - Convert Python types to AiiDA types (`orm.Float`, `orm.List`, etc.)
   - Required because `generate_slab_structures` is a `@task.calcfunction`
   - Access output via `.slabs` (not `.result.slabs`)

3. **Step 2 - Scatter-Gather**:
   - Pass the future `slab_namespace` to nested graph
   - Pass the concrete `code` object
   - The nested graph will build tasks for each slab

4. **Step 3 - Return**:
   - Expose all three output namespaces
   - These become the workflow outputs accessible via `workflow.outputs.*`

### Building and Running

**Building the workflow**:

```python
def build_pythonic_workgraph(**kwargs) -> WorkGraph:
    """Convenience wrapper returning a built workgraph."""
    if 'miller_indices' in kwargs and isinstance(kwargs['miller_indices'], tuple):
        kwargs = dict(kwargs)
        kwargs['miller_indices'] = list(kwargs['miller_indices'])
    return slab_relaxation_scatter_gather.build(**kwargs)
```

**Execution**:

```python
workflow = build_pythonic_workgraph(
    bulk_structure=structure,
    miller_indices=(1, 0, 0),
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    code_label='VASP-VTST-6.4.3@bohr',
    # ... other parameters
)

# Execute and wait for completion
workflow.run()

# Access results
for label, energy in workflow.outputs.slab_energies._sockets.items():
    if not label.startswith('_'):
        print(f'{label}: {energy.value} eV')
```

---

## Slab Relaxation Use Case

### Problem Statement

**Goal**: Relax all symmetrically distinct surface terminations of a material in parallel.

**Challenges**:
1. Number of terminations is unknown until runtime
2. Each relaxation is independent and can run in parallel
3. Need to track provenance for all calculations
4. Must handle variable numbers of inputs/outputs

**Why Scatter-Gather?**

Traditional approaches would require:
- Writing separate tasks for each slab (not scalable)
- Sequential processing (slow)
- Complex WorkChain logic (error-prone)

Scatter-gather provides:
- ✅ Automatic parallelization
- ✅ Dynamic handling of N slabs
- ✅ Clean, readable code
- ✅ Full provenance tracking

### Workflow Execution Timeline

```
Time →

t0: Submit workflow
    └─→ Start generate_slab_structures
    
t1: Slabs generated (4 terminations found)
    └─→ Enter relax_slabs_scatter
        └─→ Build graph with 4 VaspWorkChain tasks
            └─→ Submit all 4 tasks to daemon
    
t2-t5: Parallel execution
    ├─→ VaspWorkChain (slab_00) running
    ├─→ VaspWorkChain (slab_01) running  } All running
    ├─→ VaspWorkChain (slab_02) running  } simultaneously
    └─→ VaspWorkChain (slab_03) running
    
t6: First slab finishes
    └─→ extract_total_energy (slab_00) runs
    
t7-t8: Other slabs finish
    └─→ extract_total_energy tasks run in parallel
    
t9: All tasks complete
    └─→ Gather results
    └─→ Workflow finishes
    └─→ Outputs available
```

### Real Example Output

Running on Ag₃PO₄ (100) surface:

```bash
$ python slabs_relax.py

# Output:
10/06/2025 10:27:04 PM: Task: generate_slab_structures, finished.
10/06/2025 10:27:05 PM: tasks ready to run: relax_slabs_scatter
10/06/2025 10:27:06 PM: tasks ready to run: VaspWorkChain,VaspWorkChain1,VaspWorkChain2,VaspWorkChain3

# (4 VASP calculations run in parallel)

10/06/2025 10:27:45 PM: Workflow finished

Workflow finished:
  PK: 11739
  Generated slabs: 4 terminations
    Labels: ['slab_00', 'slab_01', 'slab_02', 'slab_03']
```

**Provenance Graph**:
```bash
$ verdi process show 11739

Outputs:
  generated_slabs
    slab_00      11741  StructureData
    slab_01      11742  StructureData
    slab_02      11743  StructureData
    slab_03      11744  StructureData
  relaxed_slabs
    slab_00      11780  StructureData
    slab_01      11787  StructureData
    slab_02      11794  StructureData
    slab_03      11807  StructureData
  slab_energies
    slab_00      11809  Float  (-143.05 eV)
    slab_01      11811  Float  (-162.43 eV)
    slab_02      11813  Float  (-126.11 eV)
    slab_03      11815  Float  (-132.46 eV)
```

---

## Task Connections and Data Flow

### Connection Mechanisms

#### 1. Direct Value Passing (Concrete Data)

Simple case - passing known values:

```python
@task.graph
def workflow():
    # 'code' is a concrete AiiDA node
    code = orm.load_code('vasp@cluster')
    
    # Pass it directly to tasks
    result = vasp_calculation(code=code, ...)  # ✅ Works
```

#### 2. Socket Linking (Future Values)

Connecting task outputs to inputs:

```python
@task.graph
def workflow():
    # 'structure' is a future value (socket)
    structure = generate_structure(...).result
    
    # Pass future value to next task
    relaxed = relax_structure(structure=structure)  # ✅ Works
    
    # AiiDA creates a link: generate_structure.result → relax_structure.structure
```

#### 3. Namespace Linking (Future Namespaces)

Passing entire dynamic namespaces:

```python
@task.graph
def workflow():
    # 'slabs' is a future namespace (TaskSocketNamespace)
    slabs = generate_slabs(...).slabs
    
    # Pass entire namespace to nested graph
    results = process_all_slabs(slabs=slabs)  # ✅ Works
    
    # Inside process_all_slabs:
    for label, slab in slabs.items():  # Iterates over socket names
        ...  # Creates tasks for each slab
```

#### 4. Iteration Over Future Namespaces

The key pattern for scatter-gather:

```python
@task.graph  # ← Required to iterate over future values
def process_multiple(items: t.Annotated[dict[str, T], dynamic(T)]):
    results = {}
    
    # This loop runs at GRAPH-BUILD time
    for key, item in items.items():
        # 'key' is a string (known now)
        # 'item' is a future value (socket)
        result = process_one(input=item).result  # ← Creates task
        results[key] = result  # ← Stores socket
    
    return {'results': results}
```

**What's happening**:
1. `items.items()` returns socket names and socket objects
2. The loop creates one `process_one` task per item
3. Tasks are added to the graph structure
4. At runtime, all tasks execute (in parallel if independent)

### Data Type Conversions

#### Python Types → AiiDA Types

For `@task.calcfunction`:

```python
# Input conversions
python_int = 42
aiida_int = orm.Int(python_int)

python_float = 3.14
aiida_float = orm.Float(python_float)

python_list = [1, 2, 3]
aiida_list = orm.List(list=python_list)

python_dict = {'a': 1, 'b': 2}
aiida_dict = orm.Dict(dict=python_dict)

python_bool = True
aiida_bool = orm.Bool(python_bool)

# Pass to calcfunction
result = my_calcfunction(
    count=aiida_int,
    threshold=aiida_float,
    indices=aiida_list,
    params=aiida_dict,
    flag=aiida_bool,
)
```

#### AiiDA Types → Python Types

Extracting values:

```python
int_node = orm.Int(42)
python_value = int_node.value  # → 42

list_node = orm.List(list=[1, 2, 3])
python_list = list_node.get_list()  # → [1, 2, 3]

dict_node = orm.Dict(dict={'a': 1})
python_dict = dict_node.get_dict()  # → {'a': 1}
```

### Common Connection Patterns

#### Pattern 1: Simple Chain

```python
@task.graph
def chain_workflow(x: int):
    result1 = step1(x=x).result
    result2 = step2(y=result1).result
    result3 = step3(z=result2).result
    return result3

# Data flow: x → step1 → step2 → step3 → output
# Execution: Sequential (each waits for previous)
```

#### Pattern 2: Parallel Split

```python
@task.graph
def parallel_workflow(x: int):
    # All three run in parallel (no dependencies)
    result1 = process_a(x=x).result
    result2 = process_b(x=x).result
    result3 = process_c(x=x).result
    
    # Combine results
    final = combine(a=result1, b=result2, c=result3).result
    return final

# Data flow:      ┌→ process_a ┐
#            x →  │→ process_b ├→ combine → output
#                 └→ process_c ┘
# Execution: process_a/b/c in parallel, then combine
```

#### Pattern 3: Scatter-Gather (Dynamic)

```python
@task.graph
def scatter_gather(items: t.Annotated[dict[str, T], dynamic(T)]):
    results = {}
    
    # Scatter: create N parallel tasks
    for key, item in items.items():
        results[key] = process(item=item).result
    
    # Gather: return collected results
    return {'results': results}

# Data flow:      ┌→ process(item1) ┐
#         items → │→ process(item2) ├→ collected results
#                 └→ process(item3) ┘
# Execution: All process tasks run in parallel
```

#### Pattern 4: Nested Scatter-Gather

```python
@task.graph
def nested_workflow(groups: t.Annotated[dict, dynamic(dynamic(T))]):
    all_results = {}
    
    # Outer scatter
    for group_key, group_items in groups.items():
        # Inner scatter (via nested graph)
        group_results = scatter_gather(items=group_items).results
        all_results[group_key] = group_results
    
    return {'all_results': all_results}

# Handles hierarchical parallelization
```

### Understanding Task Dependencies

AiiDA-WorkGraph automatically determines dependencies:

```python
@task.graph
def auto_dependencies():
    # Step 1: No dependencies → runs immediately
    a = task_a(x=1).result
    
    # Step 2: Depends on task_a → waits for a
    b = task_b(y=a).result
    
    # Step 3: Depends on task_a → waits for a (parallel with task_b)
    c = task_c(z=a).result
    
    # Step 4: Depends on both b and c → waits for both
    d = task_d(b=b, c=c).result
    
    return d

# Execution timeline:
# t0: task_a starts
# t1: task_a completes → task_b and task_c start in parallel
# t2: Both complete → task_d starts
# t3: task_d completes → workflow finishes
```

**Rule**: A task runs when ALL its input sockets have data.

### Socket Access Patterns

#### Common Socket Access Mistakes and Solutions

**❌ Mistake 1: Trying to access `.result` on namespace outputs**

```python
@task.calcfunction
def func() -> t.Annotated[dict, namespace(data=dynamic(int))]:
    return {'data': {...}}

# Wrong:
result = func().result.data  # ❌ Error: no 'result' socket

# Correct:
result = func().data  # ✅ Direct namespace access
```

**❌ Mistake 2: Trying to evaluate future values**

```python
@task.graph
def workflow():
    value = expensive_calc().result
    
    # Wrong:
    if value > 10:  # ❌ Error: can't evaluate future
        do_something()
    
    # Correct: pass to another task
    result = conditional_task(threshold=10, value=value)  # ✅
```

**❌ Mistake 3: Iterating without `@task.graph`**

```python
# Wrong:
def workflow():
    items = generate_items().items
    for key, item in items.items():  # ❌ Error: future namespace
        process(item)

# Correct:
@task.graph  # ← Add this decorator
def workflow():
    items = generate_items().items
    for key, item in items.items():  # ✅ Works in @task.graph
        process(item)
```

**❌ Mistake 4: Using Python types in `@task.calcfunction`**

```python
# Wrong:
@task.calcfunction
def func(x: int) -> float:  # ❌ Should use orm types
    return x * 2.5

# Correct:
@task.calcfunction
def func(x: orm.Int) -> orm.Float:  # ✅ AiiDA types
    return orm.Float(x.value * 2.5)
```

---

## Running the Workflow

### Prerequisites

```bash
# Activate AiiDA environment
source ~/envs/psteros/bin/activate

# Set the correct profile
verdi profile set-default psteros

# Check AiiDA status
verdi status

# Restart daemon to load any code changes
verdi daemon restart
```

### Mock Mode: Testing Without VASP

Test the scatter-gather pattern with lightweight calculations:

```bash
python slabs_relax.py --mock --mock-count 3 --mock-delta 0.5
```

**What it does**:
1. Generates 3 mock scalar values (0.0, 1.0, 2.0)
2. Applies a shift operation to each in parallel (+0.5)
3. Returns shifted values (0.5, 1.5, 2.5)

**Output**:
```
Mock scatter-gather workflow finished:
  Source values:
  value_00: uuid: ... value: 0.0
  value_01: uuid: ... value: 1.0
  value_02: uuid: ... value: 2.0
  Shifted values:
  value_00: uuid: ... value: 0.5
  value_01: uuid: ... value: 1.5
  value_02: uuid: ... value: 2.5
```

**Use cases**:
- ✅ Testing workflow logic
- ✅ Debugging connection issues
- ✅ Verifying parallelization
- ✅ CI/CD testing
- ✅ Training/demos

### Full VASP Workflow

Run the complete slab relaxation workflow:

```bash
python slabs_relax.py
```

**What it does**:
1. Loads Ag₃PO₄ bulk structure from `structures/ag3po4.cif`
2. Generates (100) surface slab terminations (typically 2-6 terminations)
3. Submits VASP relaxation WorkChains for all slabs in parallel
4. Extracts total energy from each completed calculation
5. Returns all structures and energies with full provenance

**Configuration** (in `slabs_relax.py`):

```python
# Structure
bulk_path = structures_dir / 'ag3po4.cif'

# Compute resources
code_label = 'VASP-VTST-6.4.3@bohr'
potential_family = 'PBE'
potential_mapping = {'Ag': 'Ag', 'P': 'P', 'O': 'O'}

# VASP parameters
slab_parameters = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 2,
    'NSW': 100,
    'EDIFFG': -0.02,
    # ...
}

# HPC options
slab_options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 40,
    },
    'queue_name': 'par40',
}

# Surface generation
miller_indices = (1, 0, 0)
min_slab_thickness = 10.0  # Å
min_vacuum_thickness = 15.0  # Å
kpoints_spacing = 0.3  # Å⁻¹
```

### Analyzing Results

#### Check Workflow Status

```bash
# Show workflow details
verdi process show <PK>

# Example output:
Property     Value
-----------  -----------------------------------------
type         WorkGraph<slab_relaxation_scatter_gather>
state        Finished [0]
pk           11739

Outputs:
  generated_slabs
    slab_00      11741  StructureData
    slab_01      11742  StructureData
    slab_02      11743  StructureData
    slab_03      11744  StructureData
  relaxed_slabs
    slab_00      11780  StructureData
    slab_01      11787  StructureData
    slab_02      11794  StructureData
    slab_03      11807  StructureData
  slab_energies
    slab_00      11809  Float
    slab_01      11811  Float
    slab_02      11813  Float
    slab_03      11815  Float
```

#### Extract Energies

```python
from aiida import orm, load_profile
load_profile()

workflow_pk = 11739
workflow = orm.load_node(workflow_pk)

# Access via outputs
for key in workflow.outputs.slab_energies._sockets:
    if not key.startswith('_'):
        energy_node = getattr(workflow.outputs.slab_energies, key)
        print(f'{key}: {energy_node.value:.6f} eV')

# Output:
# slab_00: -143.052465 eV
# slab_01: -162.426982 eV
# slab_02: -126.111982 eV
# slab_03: -132.455038 eV
```

#### Find Most Stable Termination

```python
from aiida import orm, load_profile
load_profile()

workflow_pk = 11739
workflow = orm.load_node(workflow_pk)

# Collect energies
energies = {}
for key in workflow.outputs.slab_energies._sockets:
    if not key.startswith('_'):
        energy = getattr(workflow.outputs.slab_energies, key)
        energies[key] = energy.value

# Find minimum
most_stable = min(energies, key=energies.get)
min_energy = energies[most_stable]

print(f'Most stable termination: {most_stable}')
print(f'Energy: {min_energy:.6f} eV')

# Get the corresponding structure
stable_structure = getattr(workflow.outputs.relaxed_slabs, most_stable)
print(f'Structure PK: {stable_structure.pk}')
```

#### Visualize Structure

```bash
# Export to file
verdi data core.structure export --format xsf <structure_pk> > structure.xsf

# View with your favorite visualization tool (VESTA, XCrySDen, etc.)
```

#### Check Calculation Details

```bash
# Show VASP calculation details
verdi process show <vasp_calc_pk>

# View output files
verdi calcjob outputls <vasp_calc_pk>

# Get specific output
verdi calcjob outputcat <vasp_calc_pk> OUTCAR
```

---

## Best Practices

### 1. Task Design

**✅ DO: Keep tasks focused and composable**

```python
@task.calcfunction
def extract_energy(results: orm.Dict) -> orm.Float:
    """Single responsibility: extract energy."""
    return orm.Float(results['energy'])

@task.calcfunction
def extract_forces(results: orm.Dict) -> orm.ArrayData:
    """Single responsibility: extract forces."""
    return orm.ArrayData(results['forces'])
```

**❌ DON'T: Create monolithic tasks**

```python
@task.calcfunction
def extract_everything(results: orm.Dict) -> orm.Dict:
    """Too many responsibilities."""
    return orm.Dict(dict={
        'energy': results['energy'],
        'forces': results['forces'],
        'stress': results['stress'],
        # ... many more
    })
```

### 2. Type Annotations

**✅ DO: Use complete type annotations**

```python
import typing as t
from aiida_workgraph import namespace, dynamic

@task.graph
def workflow(
    data: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)]
) -> t.Annotated[dict, namespace(results=dynamic(orm.Float))]:
    ...
```

Benefits:
- IDE autocomplete and type checking
- Clear documentation
- Catches errors early

**❌ DON'T: Skip annotations**

```python
@task.graph
def workflow(data):  # What type is data?
    ...
```

### 3. Error Handling

**✅ DO: Validate inputs and provide clear errors**

```python
@task.calcfunction
def extract_energy(results: orm.Dict) -> orm.Float:
    """Extract energy with clear error messages."""
    data = results.get_dict()
    
    # Try multiple keys
    for key in ('energy', 'total_energy', 'final_energy'):
        if key in data:
            return orm.Float(data[key])
    
    # Provide helpful error
    available_keys = ', '.join(sorted(data.keys()))
    raise ValueError(
        f'Could not find energy in results. '
        f'Available keys: {available_keys}'
    )
```

**❌ DON'T: Let obscure errors propagate**

```python
@task.calcfunction
def extract_energy(results: orm.Dict) -> orm.Float:
    # Will raise KeyError with no context
    return orm.Float(results['energy'])
```

### 4. Resource Management

**✅ DO: Load heavy resources in outer scope**

```python
@task.graph
def workflow(code_label: str, ...):
    # Load once
    code = orm.load_code(code_label)
    
    # Reuse in all tasks
    for item in items:
        calculate(code=code, item=item)
```

**❌ DON'T: Reload resources repeatedly**

```python
@task.graph
def workflow(code_label: str, ...):
    for item in items:
        code = orm.load_code(code_label)  # ❌ Inefficient
        calculate(code=code, item=item)
```

### 5. Naming Conventions

**✅ DO: Use descriptive, consistent names**

```python
# Clear purpose
@task.calcfunction
def extract_total_energy(...) -> orm.Float:
    ...

@task.graph
def relax_slabs_scatter(...):
    ...

# Consistent naming within workflow
slab_00, slab_01, slab_02  # ✅ Sequential
relaxed_slabs, slab_energies  # ✅ Descriptive namespaces
```

**❌ DON'T: Use vague names**

```python
@task.calcfunction
def process(...):  # ❌ What does it process?
    ...

@task.graph
def do_stuff(...):  # ❌ What stuff?
    ...
```

### 6. Testing Strategy

**✅ DO: Create mock versions for testing**

```python
# Production version
@task.calcfunction
def expensive_calculation(structure: orm.StructureData) -> orm.Float:
    """Runs actual DFT calculation."""
    ...

# Mock version for testing
@task.calcfunction
def mock_calculation(structure: orm.StructureData) -> orm.Float:
    """Returns fake result for testing."""
    return orm.Float(42.0)

# Workflow accepts either
@task.graph
def workflow(calc_function=expensive_calculation):
    result = calc_function(structure=...).result
```

**Benefits**:
- Fast testing without HPC
- CI/CD integration
- Development iteration

### 7. Documentation

**✅ DO: Document task behavior clearly**

```python
@task.calcfunction
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP output dictionary.
    
    Tries the following keys in order:
    1. 'total_energies.energy_extrapolated'
    2. 'total_energies.energy_no_entropy'
    3. 'total_energies.energy'
    
    Args:
        energies: VASP misc output dictionary
        
    Returns:
        Total energy in eV
        
    Raises:
        ValueError: If no compatible energy key is found
    """
    ...
```

### 8. Provenance Considerations

**✅ DO: Use `@task.calcfunction` for important operations**

```python
@task.calcfunction  # ← Stores provenance
def calculate_surface_energy(
    slab_energy: orm.Float,
    bulk_energy: orm.Float,
    area: orm.Float,
) -> orm.Float:
    """Calculate surface energy - important result!"""
    n_atoms_slab = ...
    n_atoms_bulk = ...
    return orm.Float(
        (slab_energy.value - n_atoms_slab/n_atoms_bulk * bulk_energy.value) 
        / area.value
    )
```

**❌ DON'T: Use `@task` for results you want to query later**

```python
@task  # ← No provenance stored
def calculate_surface_energy(...) -> float:
    """Can't query this result later in database!"""
    ...
```

### 9. Performance Optimization

**✅ DO: Maximize parallelization**

```python
@task.graph
def optimized_workflow(items):
    # All tasks created in one loop → run in parallel
    results = {}
    for key, item in items.items():
        results[key] = process(item).result
    return {'results': results}
```

**❌ DON'T: Create artificial dependencies**

```python
@task.graph
def slow_workflow(items):
    results = {}
    prev_result = None
    for key, item in items.items():
        # Each task waits for previous → sequential!
        result = process(item, prev=prev_result).result
        prev_result = result  # ❌ Unnecessary dependency
        results[key] = result
    return {'results': results}
```

---

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: "AttributeError: no sub-socket 'result'"

**Error**:
```python
AttributeError: TaskSocketNamespace: 'task.outputs' has no sub-socket 'result'.
Available: data, _outputs, _wait
```

**Cause**: Trying to access `.result` on a task that returns a namespace.

**Solution**: Access the namespace field directly:
```python
# Wrong:
data = generate_data().result.data  # ❌

# Correct:
data = generate_data().data  # ✅
```

#### Issue 2: "TypeError: missing required positional argument"

**Error**:
```python
TypeError: scatter_function() missing 1 required positional argument: 'data'
```

**Cause**: Data isn't being passed from previous task (link not created).

**Solutions**:

1. Check that you're accessing the correct output socket:
```python
# Wrong:
data = generate().values  # If annotation doesn't have 'values'

# Correct:
data = generate().result  # Or whatever the annotation specifies
```

2. Verify the annotation matches the return:
```python
@task.calcfunction
def generate() -> t.Annotated[dict, namespace(values=dynamic(int))]:
    return {'values': {...}}  # ✅ Key matches annotation
```

#### Issue 3: "NotExistent: no Code found with LABEL"

**Error**:
```python
aiida.common.exceptions.NotExistent: no Code found with LABEL<...>
```

**Cause**: Trying to load code inside `@task.graph` from a future value.

**Solution**: Load code in outer scope where label is concrete:
```python
# Wrong:
@task.graph
def inner_workflow(code_label: str):
    code = orm.load_code(code_label)  # ❌ code_label is future

# Correct:
@task.graph
def outer_workflow(code_label: str):
    code = orm.load_code(code_label)  # ✅ code_label is concrete
    inner_workflow(code=code)  # Pass code object
```

#### Issue 4: "GraphDeferredIllegalOperationError"

**Error**:
```python
GraphDeferredIllegalOperationError: Illegal operation on a future value
```

**Cause**: Trying to evaluate or subscript a future value outside `@task.graph`.

**Solution**: Wrap in `@task.graph` or pass to another task:
```python
# Wrong:
def workflow():
    value = calc().result
    if value > 10:  # ❌ Can't evaluate future
        ...

# Correct Option 1: Use @task.graph
@task.graph
def workflow():
    value = calc().result
    return conditional_task(threshold=10, value=value).result  # ✅

# Correct Option 2: Handle in separate task
@task
def conditional_processing(value: float, threshold: float):
    if value > threshold:  # ✅ Now value is concrete
        ...
```

#### Issue 5: AiiDA Type Errors in calcfunction

**Error**:
```python
ValueError: Cannot store object of type <class 'float'>
```

**Cause**: Passing Python types to `@task.calcfunction` instead of AiiDA types.

**Solution**: Convert all inputs to `orm.*` types:
```python
# Wrong:
@task.calcfunction
def func(x: float) -> orm.Float:  # ❌ float input
    return orm.Float(x * 2)

# Correct:
@task.calcfunction
def func(x: orm.Float) -> orm.Float:  # ✅ orm.Float input
    return orm.Float(x.value * 2)

# When calling:
result = func(x=orm.Float(3.14))  # ✅ Pass AiiDA type
```

#### Issue 6: Tasks Not Running in Parallel

**Symptom**: Tasks that should run in parallel are running sequentially.

**Diagnosis**:
```bash
# Check task dependencies
verdi process show <workflow_pk>

# Look for unnecessary links between tasks that should be independent
```

**Common causes**:

1. Accidental data dependency:
```python
# Wrong: Creates dependency chain
@task.graph
def workflow(items):
    prev = None
    for item in items:
        result = process(item, prev=prev).result  # ❌ prev creates dependency
        prev = result
```

2. Shared mutable state:
```python
# Wrong: Tasks modify shared object
shared_list = []
for item in items:
    process(item, output_list=shared_list)  # ❌ Shared state
```

**Solution**: Ensure tasks are truly independent:
```python
# Correct: No dependencies between iterations
@task.graph
def workflow(items):
    results = {}
    for key, item in items.items():
        results[key] = process(item).result  # ✅ Independent
    return {'results': results}
```

### Debugging Strategies

#### 1. Check WorkGraph Structure

```bash
# Show workflow structure
verdi process show <pk>

# Check task states
verdi process report <pk>

# View provenance graph
verdi node graph generate <pk>
```

#### 2. Enable Debug Logging

```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Run workflow
workflow.run()
```

#### 3. Test with Mock Data

Use `--mock` flag to test workflow logic without expensive calculations:

```bash
python slabs_relax.py --mock
```

#### 4. Inspect Intermediate Results

```python
from aiida import orm

# Load workflow
workflow = orm.load_node(workflow_pk)

# Check each called process
for called in workflow.called:
    print(f'{called.process_label}: {called.process_state}')
    if called.is_failed:
        print(f'  Failed: {called.exit_message}')
```

#### 5. Python Debugger

For graph-building issues, use standard Python debugger:

```python
@task.graph
def workflow(...):
    import pdb; pdb.set_trace()  # Set breakpoint
    
    # Step through graph building logic
    ...
```

### Getting Help

If you encounter issues not covered here:

1. **Check AiiDA-WorkGraph docs**: https://aiida-workgraph.readthedocs.io/
2. **AiiDA Discourse**: https://aiida.discourse.group/
3. **GitHub Issues**: https://github.com/aiidateam/aiida-workgraph/issues
4. **Local experts**: Ask team members familiar with the codebase

---

## Summary and Key Takeaways

### The Scatter-Gather Pattern

**Core Idea**: Use `@task.graph` to dynamically create multiple parallel tasks:

```python
@task.graph
def scatter_gather(items: dynamic_namespace):
    results = {}
    for key, item in items.items():  # ← Scatter
        results[key] = process(item).result
    return {'results': results}  # ← Gather
```

### When to Use Scatter-Gather

✅ **Use scatter-gather when**:
- Number of items is unknown until runtime
- Items can be processed independently
- You want automatic parallelization
- You need flexible Python control flow

❌ **Don't use scatter-gather when**:
- Number of items is fixed and small (use explicit tasks)
- Items must be processed sequentially
- Simple iteration suffices (consider `Map` zone)

### Critical Patterns

1. **`@task.calcfunction`** for AiiDA provenance
2. **`@task.graph`** for dynamic task creation
3. **`namespace(...)` annotations** for structured outputs
4. **`dynamic(Type)`** for variable-length collections
5. **Direct namespace access** (not `.result`)

### Implementation Checklist

When implementing a new scatter-gather workflow:

- [ ] Identify what varies at runtime (→ dynamic input)
- [ ] Design reusable atomic tasks
- [ ] Create scatter graph with `@task.graph`
- [ ] Use `@task.calcfunction` for provenance
- [ ] Add type annotations
- [ ] Create mock version for testing
- [ ] Document task behavior
- [ ] Test with small dataset first
- [ ] Verify parallel execution
- [ ] Check provenance graph

### This Implementation

Successfully demonstrates:
- ✅ Dynamic slab generation (2-6+ terminations)
- ✅ Parallel VASP relaxations
- ✅ Energy extraction and collection
- ✅ Full AiiDA provenance tracking
- ✅ Clean, readable code
- ✅ Mock mode for testing

**Production ready** for PS-TEROS slab screening workflows!

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

---

## Appendix: Complete Code Reference

### Minimal Working Example

Complete minimal scatter-gather workflow:

```python
import typing as t
from aiida import orm, load_profile
from aiida_workgraph import task, dynamic, namespace

load_profile()

# Step 1: Define atomic tasks
@task.calcfunction
def square(x: orm.Int) -> orm.Int:
    """Square a number."""
    return orm.Int(x.value ** 2)

@task.calcfunction
def generate_numbers(n: orm.Int) -> t.Annotated[
    dict, namespace(numbers=dynamic(orm.Int))
]:
    """Generate numbers 1 to n."""
    return {
        'numbers': {
            f'num_{i}': orm.Int(i)
            for i in range(1, n.value + 1)
        }
    }

# Step 2: Create scatter-gather graph
@task.graph
def parallel_square(
    numbers: t.Annotated[dict[str, orm.Int], dynamic(orm.Int)]
) -> t.Annotated[dict, namespace(squares=dynamic(orm.Int))]:
    """Square all numbers in parallel."""
    squares = {}
    
    # Scatter: create tasks
    for key, num in numbers.items():
        squares[key] = square(x=num).result
    
    # Gather: return results
    return {'squares': squares}

# Step 3: Build top-level workflow
@task.graph
def main_workflow(n: int):
    """Generate and square numbers."""
    nums = generate_numbers(n=orm.Int(n)).numbers
    results = parallel_square(numbers=nums).squares
    return results

# Step 4: Run
wg = main_workflow.build(n=5)
wg.run()

# Step 5: Access results
for key in wg.outputs._sockets:
    if not key.startswith('_'):
        value = getattr(wg.outputs, key)
        print(f'{key}: {value.value}')
```

---

## References and Resources

### Official Documentation

- **AiiDA-WorkGraph**: https://aiida-workgraph.readthedocs.io/
  - [Scatter-Gather Guide](https://aiida-workgraph.readthedocs.io/en/latest/howto/run_tasks_in_parallel.html)
  - [Task Types](https://aiida-workgraph.readthedocs.io/en/latest/concept/task.html)
  - [Dynamic Namespaces](https://aiida-workgraph.readthedocs.io/en/latest/concept/socket.html#dynamic-namespaces)
  - [Graph Tasks](https://aiida-workgraph.readthedocs.io/en/latest/concept/task.html#graph-task)

- **AiiDA Core**: https://aiida.readthedocs.io/
  - [CalcFunctions](https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/concepts.html#calculation-functions)
  - [Data Types](https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/data_types.html)

### Related Files in PS-TEROS

- **Zone Approach**: `teros/test_modules/zone_approach/`
  - Alternative implementation using Map context manager
  - Compare with this scatter-gather approach

---

**Last Updated**: 2025-10-06  
**Author**: PS-TEROS Development Team  
**Version**: 1.0
