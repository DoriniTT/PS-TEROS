"""WorkGraph builder for VASP parallelization benchmarking.

This module provides the main workflow for benchmarking VASP parallelization
parameters (NCORE, KPAR) to find optimal settings for a given structure
and computational resources.

Example:
    >>> from teros.core.vasp_parallelization import (
    ...     build_parallelization_benchmark_workgraph,
    ...     print_benchmark_summary,
    ... )
    >>> wg = build_parallelization_benchmark_workgraph(
    ...     structure=my_structure,
    ...     code_label='VASP-6.5.1@localwork',
    ...     potential_family='PBE',
    ...     potential_mapping={'Si': 'Si'},
    ...     num_procs=8,
    ... )
    >>> wg.submit()
    >>> # After completion:
    >>> print_benchmark_summary(wg.pk)
"""

from __future__ import annotations

import typing as t

from aiida import orm
from aiida.common.links import LinkType
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task

from teros.core.utils import get_vasp_parser_settings

from .tasks import (
    compile_benchmark_results,
    extract_benchmark_metrics,
    extract_recommended_parameters,
)
from .utils import (
    estimate_kpoints_count,
    generate_benchmark_combinations,
    prepare_benchmark_incar,
)


def _accept_benchmark_vasp_exit(task, engine):
    """Error handler that accepts VaspWorkChain exit codes for benchmarking.

    For benchmarking, we intentionally run non-converging calculations (limited NELM).
    VaspWorkChain will exit with error code 401 (max iterations exceeded) because
    the underlying VaspCalculation returns non-zero exit codes (e.g., 701 for
    electronic not converged).

    This handler raises AcceptError to signal to WorkGraph that the task should
    be marked as FINISHED despite the error, allowing child tasks (metrics
    extraction) to run.

    Args:
        task: The WorkGraph task that failed.
        engine: The ErrorHandlerManager instance (provides access to WorkGraph internals).

    Raises:
        AcceptError: Always raised to accept the task as finished.
    """
    from aiida_workgraph.engine.error_handler_manager import AcceptError

    raise AcceptError(
        "Benchmark VASP task accepted (non-convergence expected with limited NELM)"
    )


def build_parallelization_benchmark_workgraph(
    structure: orm.StructureData,
    code_label: str,
    potential_family: str,
    potential_mapping: dict[str, str],
    num_procs: int,
    # Benchmark configuration
    num_kpoints: int | None = None,
    ncore_values: list[int] | None = None,
    kpar_values: list[int] | None = None,
    nelm: int = 3,
    # VASP configuration
    base_incar: dict[str, t.Any] | None = None,
    options: dict[str, t.Any] | None = None,
    kpoints_spacing: float = 0.03,
    clean_workdir: bool = True,
    # Workflow control
    max_concurrent_jobs: int | None = None,
    weight_time: float = 0.7,
    weight_memory: float = 0.3,
    # WorkGraph settings
    name: str = "VaspParallelizationBenchmark",
) -> WorkGraph:
    """Build WorkGraph for benchmarking VASP parallelization parameters.

    Creates a workflow that runs multiple short VASP calculations with
    different NCORE/KPAR combinations, extracts timing/memory metrics,
    and ranks the results to find optimal parallelization settings.

    Args:
        structure: Structure to benchmark with.
        code_label: AiiDA code label for VASP (e.g., 'VASP-6.5.1@localwork').
        potential_family: POTCAR family name (e.g., 'PBE', 'PBE.54').
        potential_mapping: Element to POTCAR mapping (e.g., {'Si': 'Si'}).
        num_procs: Total number of MPI processes available.

        num_kpoints: Number of k-points. Auto-estimated from structure if None.
        ncore_values: Explicit NCORE values to test. Auto-generated if None.
        kpar_values: Explicit KPAR values to test. Auto-generated if None.
        nelm: Number of electronic steps for benchmark (1-3 recommended).

        base_incar: Base INCAR parameters (ENCUT, ISMEAR, etc.).
        options: Scheduler options (resources, queue, walltime).
        kpoints_spacing: K-points spacing in Angstrom^-1.
        clean_workdir: Whether to clean VASP working directory after completion.

        max_concurrent_jobs: Maximum number of parallel benchmark jobs.
        weight_time: Weight for elapsed time in ranking (0-1).
        weight_memory: Weight for memory usage in ranking (0-1).

        name: Name for the WorkGraph.

    Returns:
        WorkGraph ready to be submitted.

        WorkGraph outputs (after completion):
        - benchmark_results: Dict with full results including:
            - all_results: List of all benchmark runs with timing/memory
            - best_by_time: Parameters with fastest elapsed time
            - best_by_memory: Parameters with lowest memory usage
            - recommended: Recommended parameters (best combined score)
        - summary: Dict with key metrics for quick access:
            - recommended_ncore: Int - recommended NCORE value
            - recommended_kpar: Int - recommended KPAR value
            - recommended_elapsed_time_s: Float - best elapsed time
            - recommended_memory_mb: Float - memory usage
            - num_successful: Int - successful benchmark count
            - num_failed: Int - failed benchmark count

    Example:
        >>> wg = build_parallelization_benchmark_workgraph(
        ...     structure=structure,
        ...     code_label='VASP-6.5.1@localwork',
        ...     potential_family='PBE',
        ...     potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     num_procs=8,
        ...     nelm=3,
        ...     base_incar={'ENCUT': 400, 'ISMEAR': 0, 'SIGMA': 0.05},
        ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ...     max_concurrent_jobs=2,
        ... )
        >>> wg.submit()
    """
    # Load VASP code
    code = orm.load_code(code_label)

    # Estimate k-points if not provided
    if num_kpoints is None:
        num_kpoints = estimate_kpoints_count(structure, kpoints_spacing)

    # Generate benchmark combinations
    combinations = generate_benchmark_combinations(
        num_procs=num_procs,
        num_kpoints=num_kpoints,
        ncore_values=ncore_values,
        kpar_values=kpar_values,
    )

    if not combinations:
        raise ValueError(
            f"No valid NCORE/KPAR combinations found for num_procs={num_procs}"
        )

    # Default options
    if options is None:
        options = {
            "resources": {
                "num_machines": 1,
                "num_mpiprocs_per_machine": num_procs,
            },
        }

    # Default base INCAR
    if base_incar is None:
        base_incar = {
            "ENCUT": 400,
            "ISMEAR": 0,
            "SIGMA": 0.05,
            "EDIFF": 1e-4,  # Looser for benchmark
        }

    # Wrap VaspWorkChain as task
    VaspWorkChain = WorkflowFactory("vasp.v2.vasp")
    VaspTask = task(VaspWorkChain)

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # Set max concurrent jobs
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    # Create benchmark tasks for each combination
    metrics_tasks = {}

    for combo in combinations:
        label = combo["label"]
        ncore = combo["ncore"]
        kpar = combo["kpar"]

        # Prepare INCAR with parallelization parameters
        benchmark_incar = prepare_benchmark_incar(
            base_incar=base_incar,
            ncore=ncore,
            kpar=kpar,
            nelm=nelm,
            force_short_scf=True,
        )

        # Create VASP task with benchmark-specific settings:
        # - max_iterations=1: Prevent VaspWorkChain from relaunching on non-convergence
        # - handler_overrides: Disable handlers that would restart the calculation
        #   and enable handler_always_attach_outputs to expose outputs despite errors
        #
        # This is critical for benchmarking because we intentionally use few
        # electronic steps (NELM=3) which won't converge. Without these settings,
        # VaspWorkChain would detect non-convergence and restart with ALGO=normal
        # and NELM=150, defeating the purpose of quick benchmarking.
        #
        # The handler_always_attach_outputs ensures that even when VaspCalculation
        # returns a non-zero exit code (e.g., 701 for electronic not converged),
        # the outputs (misc, retrieved, remote_folder) are still attached to the
        # VaspWorkChain, making them available for the metrics extraction task.
        handler_overrides = {
            # Disable handlers that restart calculations on non-convergence
            "handler_electronic_conv": {"enabled": False},
            "handler_electronic_conv_alt": {"enabled": False},
            "handler_unfinished_calc_generic": {"enabled": False},
            "handler_unfinished_calc_generic_alt": {"enabled": False},
            "handler_unfinished_calc_ionic": {"enabled": False},
            "handler_unfinished_calc_ionic_alt": {"enabled": False},
            # Disable convergence checks that would mark calculation as failed
            "check_electronic_converged": {"enabled": False},
            "check_ionic_converged": {"enabled": False},
            # Enable handler that attaches outputs even when calculation "fails"
            # This is critical: without it, VaspWorkChain won't expose outputs
            # and downstream tasks can't access the retrieved OUTCAR
            "handler_always_attach_outputs": {"enabled": True},
        }

        vasp_task = wg.add_task(
            VaspTask,
            name=f"vasp_{label}",
            structure=structure,
            code=code,
            parameters=orm.Dict(dict={"incar": benchmark_incar}),
            options=orm.Dict(dict=options),
            potential_family=potential_family,
            potential_mapping=orm.Dict(dict=potential_mapping),
            kpoints_spacing=kpoints_spacing,
            clean_workdir=clean_workdir,
            settings=orm.Dict(dict=get_vasp_parser_settings()),
            # Benchmark-specific: prevent restarts and handler interventions
            max_iterations=orm.Int(1),
            handler_overrides=orm.Dict(dict=handler_overrides),
        )

        # Add error handler to accept expected VaspWorkChain exit codes for benchmarking.
        # This allows the workflow to continue even though VaspWorkChain "fails"
        # due to the intentionally non-converging benchmark calculation.
        #
        # Expected exit codes:
        # - 401: ERROR_MAXIMUM_ITERATION_EXCEEDED (when handlers are disabled)
        # - 503: ERROR_ELECTRONIC_STRUCTURE_NOT_CONVERGED (when using handler_always_attach_outputs)
        # - 504: ERROR_IONIC_NOT_CONVERGED (if testing with ionic relaxation)
        vasp_task.add_error_handler(
            {
                "accept_benchmark_exit": {
                    "executor": _accept_benchmark_vasp_exit,
                    "exit_codes": [401, 503, 504],
                    "max_retries": 1,  # Must be >= 1 for handler to run
                }
            }
        )

        # Create metrics extraction task
        metrics_task = wg.add_task(
            extract_benchmark_metrics,
            name=f"metrics_{label}",
            retrieved=vasp_task.outputs.retrieved,
            label=orm.Str(label),
            ncore=orm.Int(ncore),
            kpar=orm.Int(kpar),
        )

        metrics_tasks[label] = metrics_task

    # Create compilation task
    compile_kwargs = {
        label: metrics_task.outputs.result
        for label, metrics_task in metrics_tasks.items()
    }
    compile_kwargs["weight_time"] = orm.Float(weight_time)
    compile_kwargs["weight_memory"] = orm.Float(weight_memory)

    compile_task = wg.add_task(
        compile_benchmark_results,
        name="compile_results",
        **compile_kwargs,
    )

    # Extract recommended parameters as a summary
    extract_task = wg.add_task(
        extract_recommended_parameters,
        name="extract_recommended",
        benchmark_results=compile_task.outputs.result,
    )

    # Expose outputs - the summary contains all key results in an easy format
    # Full results: benchmark_results.all_results, .best_by_time, .best_by_memory
    # Quick access: summary has recommended_ncore, recommended_kpar, etc.
    wg.outputs.benchmark_results = compile_task.outputs.result
    wg.outputs.summary = extract_task.outputs.result

    return wg


def get_benchmark_results(workgraph) -> dict[str, t.Any]:
    """Extract benchmark results from a completed WorkGraph.

    Args:
        workgraph: WorkGraph node (can be PK, UUID string, or node object).

    Returns:
        Dict with benchmark results including:
        - 'all_results': List of all benchmark runs with timing/memory
        - 'best_by_time': Fastest parameters
        - 'best_by_memory': Lowest memory parameters
        - 'recommended': Best combined score parameters

    Raises:
        ValueError: If results cannot be found.
    """
    # Load node if PK or UUID provided
    if isinstance(workgraph, (int, str)):
        workgraph = orm.load_node(workgraph)

    # Try direct output access for live WorkGraph
    if hasattr(workgraph, "outputs") and hasattr(
        workgraph.outputs, "benchmark_results"
    ):
        try:
            return workgraph.outputs.benchmark_results.get_dict()
        except Exception:
            pass

    # For stored nodes, traverse CALL_CALC links
    if hasattr(workgraph, "base") and hasattr(workgraph.base, "links"):
        # Look for compile_results calcfunction
        calcfuncs = workgraph.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
        for link in calcfuncs.all():
            if "compile" in link.link_label.lower():
                for out_link in link.node.base.links.get_outgoing(
                    link_type=LinkType.CREATE
                ).all():
                    if out_link.link_label == "result":
                        return out_link.node.get_dict()

        # Also check RETURN links for WorkGraph outputs
        returns = workgraph.base.links.get_outgoing(link_type=LinkType.RETURN)
        for link in returns.all():
            if "benchmark" in link.link_label.lower():
                if isinstance(link.node, orm.Dict):
                    return link.node.get_dict()

    raise ValueError(
        "Could not find benchmark results. "
        "Ensure the WorkGraph has completed successfully."
    )


def get_benchmark_summary(workgraph) -> dict[str, t.Any]:
    """Extract the summary Dict from a completed WorkGraph.

    This returns the simplified summary with just the key metrics:
    recommended_ncore, recommended_kpar, recommended_elapsed_time_s, etc.

    Args:
        workgraph: WorkGraph node (can be PK, UUID string, or node object).

    Returns:
        Dict with summary metrics.

    Raises:
        ValueError: If summary cannot be found.
    """
    # Load node if PK or UUID provided
    if isinstance(workgraph, (int, str)):
        workgraph = orm.load_node(workgraph)

    # For stored nodes, traverse RETURN links
    if hasattr(workgraph, "base") and hasattr(workgraph.base, "links"):
        returns = workgraph.base.links.get_outgoing(link_type=LinkType.RETURN)
        for link in returns.all():
            if "summary" in link.link_label.lower():
                if isinstance(link.node, orm.Dict):
                    return link.node.get_dict()

        # Also try CALL_CALC for extract_recommended
        calcfuncs = workgraph.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
        for link in calcfuncs.all():
            if "extract" in link.link_label.lower():
                for out_link in link.node.base.links.get_outgoing(
                    link_type=LinkType.CREATE
                ).all():
                    if out_link.link_label == "result":
                        return out_link.node.get_dict()

    raise ValueError(
        "Could not find benchmark summary. "
        "Ensure the WorkGraph has completed successfully."
    )


def print_benchmark_summary(workgraph, show_all: bool = True) -> None:
    """Print a formatted summary of benchmark results.

    Args:
        workgraph: WorkGraph node (can be PK, UUID string, or node object).
        show_all: Whether to show all results or just top 5.

    Example output:
        VASP Parallelization Benchmark Results
        ======================================
        Rank | NCORE | KPAR | Time (s) | Memory (MB) | Score
        -----|-------|------|----------|-------------|------
          1  |   4   |   1  |   12.34  |     512.0   | 0.000
          2  |   8   |   1  |   14.56  |     384.0   | 0.156
          3  |   2   |   1  |   15.78  |     640.0   | 0.234

        Recommended: NCORE=4, KPAR=1
          Best balance of time (12.34s) and memory (512 MB)
    """
    try:
        results = get_benchmark_results(workgraph)
    except ValueError as e:
        print(f"Error: {e}")
        return

    print("\nVASP Parallelization Benchmark Results")
    print("=" * 60)

    # Table header
    print(
        f"{'Rank':^5} | {'NCORE':^5} | {'KPAR':^4} | {'Time (s)':^10} | "
        f"{'Memory (MB)':^11} | {'Score':^6} | {'Status':^7}"
    )
    print(
        "-" * 5
        + "-+-"
        + "-" * 5
        + "-+-"
        + "-" * 4
        + "-+-"
        + "-" * 10
        + "-+-"
        + "-" * 11
        + "-+-"
        + "-" * 6
        + "-+-"
        + "-" * 7
    )

    # Table rows
    all_results = results.get("all_results", [])
    if not show_all:
        all_results = all_results[:5]

    for r in all_results:
        rank = r.get("rank", "-")
        ncore = r.get("ncore", "-")
        kpar = r.get("kpar", "-")
        time_s = r.get("elapsed_time")
        time_str = f"{time_s:.2f}" if time_s is not None else "-"
        mem_mb = r.get("memory_mb")
        mem_str = f"{mem_mb:.1f}" if mem_mb is not None else "-"
        score = r.get("score")
        score_str = f"{score:.3f}" if score is not None else "-"
        status = r.get("status", "unknown")[:7]

        print(
            f"{rank:^5} | {ncore:^5} | {kpar:^4} | {time_str:^10} | "
            f"{mem_str:^11} | {score_str:^6} | {status:^7}"
        )

    # Summary
    print()
    recommended = results.get("recommended")
    if recommended:
        print(f"Recommended: NCORE={recommended['ncore']}, KPAR={recommended['kpar']}")
        time_s = recommended.get("elapsed_time")
        mem_mb = recommended.get("memory_mb")
        time_str = f"{time_s:.2f}s" if time_s is not None else "N/A"
        mem_str = f"{mem_mb:.0f} MB" if mem_mb is not None else "N/A"
        print(f"  Best balance of time ({time_str}) and memory ({mem_str})")
    else:
        print("No recommendation available (all benchmarks may have failed).")

    # Stats
    num_ok = results.get("num_successful", 0)
    num_fail = results.get("num_failed", 0)
    print(f"\nBenchmarks: {num_ok} successful, {num_fail} failed")

    weights = results.get("weights", {})
    if weights:
        print(
            f"Ranking weights: time={weights.get('time', 0.7):.1%}, "
            f"memory={weights.get('memory', 0.3):.1%}"
        )
