# Working with PS-TEROS and AiiDA-WorkGraph

This is a comprehensive guide for developing and working with the **PS-TEROS** code using **AiiDA** and **AiiDA-WorkGraph**.

---

## Quick Start Commands

* **Before you begin**

Make sure that you are in the 'psteros' profile in AiiDA by doing 'verdi profile set-default psteros', then do 'verdi status' to check if everything is working fine.
If the daemon is not running, you can start it by doing 'verdi daemon start'.

* **After every modification**

  ```bash
  verdi daemon restart
  ```

* **Python binary**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python
  ```

* **Launch workflow**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python /home/thiagotd/git/PS-TEROS/examples/vasp/update_psteros/psteros_vasp.py
  ```

* **Wait for results**: Use `sleep 15` (or more) before analyzing nodes

* **Analyze results**

  ```bash
  verdi process show <PK>
  verdi process report <PK>
  ```

* **Clear Python cache**

  ```bash
  find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
  ```

---

## Essential Documentation

* [AiiDA Core](https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/index.html)
* [AiiDA-WorkGraph](https://aiida-workgraph.readthedocs.io/en/latest/)
* [Materials Science Tutorial](https://aiida-workgraph.readthedocs.io/en/latest/tutorial/autogen/materials_science_ase.html)

---

## Key Concepts — AiiDA-WorkGraph (Latest Version)

### Task

A **Task** is the basic building block of a WorkGraph. It has **inputs**, **outputs**, and an **executor**. The executor can be:

* A Python function
* AiiDA components: `calcfunction`, `workfunction`, `calcjob`, `WorkChain`, or `ProcessBuilder`

You can create a Task in three ways:

#### 1. Decorator

```python
from aiida_workgraph import task
from aiida import orm

@task
def add(x, y):
    return x + y

@task.calcfunction
def multiply(x, y):
    return orm.Float(x + y)

# Export to HTML visualization
add()._node.to_html()
```

You can inspect the input and output sockets of a task:

```python
add1 = add()._node
print('Inputs:', add1.get_input_names())
print('Outputs:', add1.get_output_names())
```

#### 2. Build from Callable

You can also build a task from an existing Python function:

```python
from scipy.linalg import norm
from aiida_workgraph import task, WorkGraph

NormTask = task()(norm)
wg = WorkGraph()
norm_task = wg.add_task(NormTask, name='norm1')
norm_task.to_html()
```

For multiple outputs, define them explicitly:

```python
def calculate_stats(data):
    import numpy as np
    return np.mean(data), np.std(data)

StatsTask = task(outputs=['mean', 'std'])(calculate_stats)
```

#### 3. Define a Task Class

You can create a custom Task class:

```python
from aiida_workgraph.task import Task

class MyAdd(Task):
    identifier = 'MyAdd'
    name = 'MyAdd'
    node_type = 'calcfunction'
    catalog = 'Test'

    _executor = {
        'module_path': 'aiida_workgraph.executors.test',
        'callable_name': 'add',
    }

    def update_sockets(self):
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.Any', 'x')
        self.add_input('workgraph.Any', 'y')
        self.add_output('workgraph.Any', 'sum')
```

You can then use `MyAdd` in a WorkGraph by its identifier:

```python
wg = WorkGraph()
add1_task = wg.add_task(MyAdd, name='add1')
```

---

### WorkGraph

A **WorkGraph** is a collection of tasks and links. You can define tasks, connect their inputs and outputs, and execute them.

#### Creating and Running a WorkGraph

```python
from aiida_workgraph import WorkGraph, task

wg = WorkGraph(name='my_first_workgraph')

@task()
def add(x, y):
    return x + y

add1 = wg.add_task(add, name='add1')
add2 = wg.add_task(add, name='add2')

wg.add_link(add1.outputs.result, add2.inputs.x)

wg.to_html()

from aiida import load_profile
load_profile()
wg.run(inputs={'add1': {'x': 1, 'y': 2}, 'add2': {'y': 3}})
```

You can define **graph-level inputs/outputs** for cleaner workflows, reuse inputs, and expose only key results:

```python
wg = WorkGraph('graph_inputs_outputs')

wg.inputs.x = 2
wg.add_task(add, 'add1', x=wg.inputs.x, y=3)
wg.add_task(add, 'add2', x=wg.inputs.x, y=wg.tasks.add1.outputs.result)

wg.outputs.sum1 = wg.tasks.add1.outputs.result
wg.outputs.sum2 = wg.tasks.add2.outputs.result
wg.run()
```

---

### Sockets

**Sockets** are the input and output connection points between tasks. WorkGraph automatically creates sockets from function arguments and return values.

Example:

```python
@task()
def multiply(x, y):
    return x * y

wg = WorkGraph()
task1 = wg.add_task(multiply, x=3, y=4)
print('Input sockets:', task1.get_input_names())
print('Output sockets:', task1.get_output_names())
```

You can define **custom output sockets**:

```python
from aiida_workgraph import spec
from typing import Any

@task
def add_and_subtract(x, y) -> spec.namespace(sum=Any, difference=Any):
    return {'sum': x + y, 'difference': x - y}
```

Sockets support Python type hints (e.g., `int → orm.Int`, `float → orm.Float`). Default values can be set via function arguments. You can also organize sockets with **namespaces** and even create **dynamic namespaces** for outputs generated at runtime.

---

### Graph Task

A **graph task** is a Python function decorated with `@task.graph` that builds and returns a WorkGraph. It allows you to use Python logic (conditionals, loops, etc.) to construct a workflow dynamically.

```python
from aiida_workgraph import task, WorkGraph

@task
def add(x, y):
    return x + y

@task.graph
def my_workflow(x, y):
    outputs = add(x=x, y=y)
    return outputs.result

wg = my_workflow.build(x=1, y=2)
wg.run()
print('Workflow outputs:', wg.outputs.result)
```

**Graph tasks vs Context Managers**:

* **Graph tasks** use plain Python to build subgraphs before execution.
* **Context managers** (`If`, `While`, `Map`) insert logic nodes into the provenance graph.

This allows you to choose between flexibility and full provenance depending on your needs.

---

## Summary Table: Graph Task vs Context Managers

| Feature         | Graph Task                         | Context Managers (`If`, `While`, `Map`) |
| --------------- | ---------------------------------- | --------------------------------------- |
| Flexibility     | High (full Python power)           | DSL-limited                             |
| Provenance      | Less detailed (setup logic hidden) | Fully detailed                          |
| Structure       | Nested workflows                   | Flat workflows                          |
| Iterative Logic | Not direct (no while loop)         | Native with `While`/`Map`               |
| Best For        | Reusable, complex setup logic      | Strict provenance, iterative workflows  |

---

This markdown file preserves the structure of the original document but uses proper Markdown syntax, code blocks, and tables to make it more readable and portable.
