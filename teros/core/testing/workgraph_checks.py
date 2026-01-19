"""WorkGraph structure validation.

This module provides utilities to inspect and validate WorkGraph
structures before submission. Helps catch wiring issues early.

Example:
    >>> from teros.core.testing import check_workgraph_wiring, print_workgraph_structure
    >>> from teros.core.fukui import build_fukui_workgraph
    >>> wg = build_fukui_workgraph(...)  # Don't submit
    >>> print_workgraph_structure(wg)
    >>> warnings = check_workgraph_wiring(wg)
"""

from typing import List, Set, Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class WorkGraphAnalysis:
    """Analysis result for a WorkGraph.

    Attributes:
        name: WorkGraph name
        num_tasks: Number of tasks
        task_names: List of task names
        warnings: List of potential issues
        connections: List of (from_task, from_socket, to_task, to_socket) tuples
        orphan_tasks: Tasks with no connections
        root_tasks: Tasks with no input connections (entry points)
        leaf_tasks: Tasks with no output connections (exit points)
    """

    name: str
    num_tasks: int
    task_names: List[str]
    warnings: List[str] = field(default_factory=list)
    connections: List[tuple] = field(default_factory=list)
    orphan_tasks: List[str] = field(default_factory=list)
    root_tasks: List[str] = field(default_factory=list)
    leaf_tasks: List[str] = field(default_factory=list)
    info: Dict[str, Any] = field(default_factory=dict)


def _get_tasks_dict(workgraph) -> Dict[str, Any]:
    """Get tasks as a dictionary, handling different WorkGraph API versions.

    The aiida-workgraph API changed - tasks is now a NodeCollection
    that doesn't have .keys() or .items() methods.
    """
    tasks_dict = {}

    # Try new API first (NodeCollection)
    try:
        for task in workgraph.tasks:
            tasks_dict[task.name] = task
        return tasks_dict
    except TypeError:
        pass

    # Try old API (.keys()/.items())
    try:
        return dict(workgraph.tasks.items())
    except (AttributeError, TypeError):
        pass

    # Fallback: try direct dict conversion
    try:
        return dict(workgraph.tasks)
    except (TypeError, ValueError):
        pass

    return tasks_dict


def analyze_workgraph(workgraph) -> WorkGraphAnalysis:
    """Perform comprehensive analysis of WorkGraph structure.

    Args:
        workgraph: Built (but not submitted) WorkGraph

    Returns:
        WorkGraphAnalysis with structure information and warnings

    Example:
        >>> analysis = analyze_workgraph(wg)
        >>> print(f"Tasks: {analysis.num_tasks}")
        >>> print(f"Warnings: {analysis.warnings}")
    """
    # Get tasks as dictionary
    tasks_dict = _get_tasks_dict(workgraph)

    analysis = WorkGraphAnalysis(
        name=workgraph.name,
        num_tasks=len(tasks_dict),
        task_names=list(tasks_dict.keys()),
    )

    # Track which tasks have connections
    has_input_connection: Set[str] = set()
    has_output_connection: Set[str] = set()

    # Internal/system tasks to ignore
    system_tasks = {'graph_inputs', 'graph_outputs', 'graph_ctx',
                    'WORKGRAPH_INPUT', 'WORKGRAPH_OUTPUT'}

    # Analyze each task
    for name, task in tasks_dict.items():
        if name in system_tasks:
            continue

        # Get inputs - handle both dict-like and list-like access
        inputs = {}
        if hasattr(task, 'inputs'):
            try:
                # Try dict-like access
                inputs = dict(task.inputs.items()) if hasattr(task.inputs, 'items') else {}
            except (TypeError, AttributeError):
                # Try iteration
                try:
                    for inp in task.inputs:
                        inp_name = getattr(inp, 'name', str(inp))
                        inputs[inp_name] = inp
                except TypeError:
                    pass

        # Get outputs - handle both dict-like and list-like access
        outputs = {}
        if hasattr(task, 'outputs'):
            try:
                outputs = dict(task.outputs.items()) if hasattr(task.outputs, 'items') else {}
            except (TypeError, AttributeError):
                try:
                    for out in task.outputs:
                        out_name = getattr(out, 'name', str(out))
                        outputs[out_name] = out
                except TypeError:
                    pass

        # Check inputs for connections
        for input_name, input_socket in inputs.items():
            links = getattr(input_socket, 'links', None)
            if links:
                has_input_connection.add(name)
                for link in links:
                    if hasattr(link, 'node') and hasattr(link, 'socket'):
                        source_name = getattr(link.node, 'name', 'unknown')
                        source_socket = getattr(link.socket, 'name', 'unknown')
                        analysis.connections.append(
                            (source_name, source_socket, name, input_name)
                        )

        # Check outputs for connections
        for output_name, output_socket in outputs.items():
            links = getattr(output_socket, 'links', None)
            if links:
                has_output_connection.add(name)

        # Check for required but unconnected inputs
        for input_name, input_socket in inputs.items():
            is_required = getattr(input_socket, 'required', False)
            has_links = bool(getattr(input_socket, 'links', None))
            has_value = getattr(input_socket, 'value', None) is not None

            if is_required and not has_links and not has_value:
                analysis.warnings.append(
                    f"Task '{name}' has unconnected required input '{input_name}'"
                )

    # Identify orphans, roots, and leaves
    user_tasks = set(tasks_dict.keys()) - system_tasks

    for task_name in user_tasks:
        has_in = task_name in has_input_connection
        has_out = task_name in has_output_connection

        if not has_in and not has_out:
            analysis.orphan_tasks.append(task_name)
            analysis.warnings.append(
                f"Task '{task_name}' appears orphaned (no connections)"
            )
        elif not has_in:
            analysis.root_tasks.append(task_name)
        elif not has_out:
            analysis.leaf_tasks.append(task_name)

    return analysis


def check_workgraph_wiring(workgraph) -> List[str]:
    """Check WorkGraph for common wiring issues.

    A simpler interface that just returns warnings.

    Args:
        workgraph: Built (but not submitted) WorkGraph

    Returns:
        List of warning messages

    Example:
        >>> warnings = check_workgraph_wiring(wg)
        >>> for w in warnings:
        ...     print(f"Warning: {w}")
    """
    analysis = analyze_workgraph(workgraph)
    return analysis.warnings


def print_workgraph_structure(
    workgraph,
    show_connections: bool = True,
    show_inputs: bool = False,
    show_outputs: bool = False,
) -> None:
    """Print WorkGraph structure for inspection.

    Args:
        workgraph: Built (but not submitted) WorkGraph
        show_connections: Show task connections
        show_inputs: Show WorkGraph-level inputs
        show_outputs: Show WorkGraph-level outputs

    Example:
        >>> print_workgraph_structure(wg)
        WorkGraph: fukui_workgraph
        ============================================================
        Tasks (5):
          - scf_neutral: VaspWorkChainTask
          - scf_charged_0_05: VaspWorkChainTask
          ...
    """
    analysis = analyze_workgraph(workgraph)
    tasks_dict = _get_tasks_dict(workgraph)

    # System tasks to skip
    system_tasks = {'graph_inputs', 'graph_outputs', 'graph_ctx',
                    'WORKGRAPH_INPUT', 'WORKGRAPH_OUTPUT'}

    print(f"\nWorkGraph: {analysis.name}")
    print("=" * 60)

    # Tasks summary (excluding system tasks)
    user_tasks = [n for n in analysis.task_names if n not in system_tasks]
    print(f"\nTasks ({len(user_tasks)}):")
    for name in user_tasks:
        task = tasks_dict.get(name)
        if task:
            task_type = task.__class__.__name__
            print(f"  - {name}: {task_type}")

    # Root and leaf tasks
    if analysis.root_tasks:
        print(f"\nEntry points (no input connections):")
        for name in analysis.root_tasks:
            print(f"  - {name}")

    if analysis.leaf_tasks:
        print(f"\nExit points (no output connections):")
        for name in analysis.leaf_tasks:
            print(f"  - {name}")

    # Connections
    if show_connections and analysis.connections:
        print("\nConnections:")
        for src_task, src_socket, dst_task, dst_socket in analysis.connections:
            print(f"  {src_task}.{src_socket} -> {dst_task}.{dst_socket}")

    # WorkGraph-level inputs
    if show_inputs and hasattr(workgraph, 'inputs'):
        print("\nWorkGraph inputs:")
        try:
            for inp in workgraph.inputs:
                inp_name = getattr(inp, 'name', str(inp))
                value = getattr(inp, 'value', None)
                if value is not None:
                    print(f"  - {inp_name}: {type(value).__name__}")
                else:
                    print(f"  - {inp_name}: (not set)")
        except TypeError:
            pass

    # WorkGraph-level outputs
    if show_outputs and hasattr(workgraph, 'outputs'):
        print("\nWorkGraph outputs:")
        try:
            for out in workgraph.outputs:
                out_name = getattr(out, 'name', str(out))
                print(f"  - {out_name}")
        except TypeError:
            pass

    # Warnings
    if analysis.warnings:
        print("\nWarnings:")
        for w in analysis.warnings:
            print(f"  ! {w}")

    print("=" * 60)


def get_workgraph_execution_order(workgraph) -> List[List[str]]:
    """Determine task execution order (topological sort with levels).

    Returns tasks grouped by execution level - tasks in the same
    level can run in parallel.

    Args:
        workgraph: Built WorkGraph

    Returns:
        List of lists, where each inner list contains task names
        that can execute in parallel

    Example:
        >>> order = get_workgraph_execution_order(wg)
        >>> for i, level in enumerate(order):
        ...     print(f"Level {i}: {level}")
        Level 0: ['scf_neutral', 'scf_charged_0_05']
        Level 1: ['interpolation']
    """
    analysis = analyze_workgraph(workgraph)

    # System tasks to skip
    system_tasks = {'graph_inputs', 'graph_outputs', 'graph_ctx',
                    'WORKGRAPH_INPUT', 'WORKGRAPH_OUTPUT'}

    # Build dependency graph
    dependencies: Dict[str, Set[str]] = {}
    for task_name in analysis.task_names:
        if task_name in system_tasks:
            continue
        dependencies[task_name] = set()

    for src_task, _, dst_task, _ in analysis.connections:
        if dst_task in dependencies and src_task in dependencies:
            dependencies[dst_task].add(src_task)

    # Kahn's algorithm for topological sort with levels
    in_degree = {task: len(deps) for task, deps in dependencies.items()}
    levels: List[List[str]] = []

    remaining = set(dependencies.keys())

    while remaining:
        # Find all tasks with no remaining dependencies
        ready = [t for t in remaining if in_degree[t] == 0]

        if not ready:
            # Cycle detected
            break

        levels.append(sorted(ready))

        for task in ready:
            remaining.remove(task)
            # Update in-degrees
            for other in remaining:
                if task in dependencies[other]:
                    in_degree[other] -= 1

    return levels


def estimate_workgraph_complexity(workgraph) -> Dict[str, Any]:
    """Estimate computational complexity of WorkGraph.

    Provides rough estimates based on structure.

    Args:
        workgraph: Built WorkGraph

    Returns:
        Dictionary with complexity metrics

    Example:
        >>> info = estimate_workgraph_complexity(wg)
        >>> print(f"VASP jobs: {info['vasp_tasks']}")
        >>> print(f"Max parallel: {info['max_parallel']}")
    """
    analysis = analyze_workgraph(workgraph)
    order = get_workgraph_execution_order(workgraph)
    tasks_dict = _get_tasks_dict(workgraph)

    # System tasks to skip
    system_tasks = {'graph_inputs', 'graph_outputs', 'graph_ctx',
                    'WORKGRAPH_INPUT', 'WORKGRAPH_OUTPUT'}

    # Count VASP-related tasks
    vasp_tasks = 0
    calcfunctions = 0
    graph_tasks = 0
    other_tasks = 0

    for name in analysis.task_names:
        if name in system_tasks:
            continue
        task = tasks_dict.get(name)
        if task:
            task_type = task.__class__.__name__.lower()

            if 'vasp' in task_type:
                vasp_tasks += 1
            elif 'calcfunction' in task_type or 'calcfunc' in task_type or 'function' in task_type:
                calcfunctions += 1
            elif 'graph' in task_type:
                graph_tasks += 1
            else:
                other_tasks += 1

    # Count user tasks (excluding system)
    user_task_count = len([n for n in analysis.task_names if n not in system_tasks])

    return {
        'total_tasks': user_task_count,
        'vasp_tasks': vasp_tasks,
        'calcfunctions': calcfunctions,
        'graph_tasks': graph_tasks,
        'other_tasks': other_tasks,
        'execution_levels': len(order),
        'max_parallel': max(len(level) for level in order) if order else 0,
        'root_tasks': len(analysis.root_tasks),
        'leaf_tasks': len(analysis.leaf_tasks),
        'has_warnings': len(analysis.warnings) > 0,
    }
