"""AiiDA calcfunction tasks for VASP parallelization benchmarking.

This module provides @task.calcfunction helpers for:
- Extracting timing/memory metrics from VASP OUTCAR
- Compiling and ranking benchmark results
"""

from __future__ import annotations

import typing as t

from aiida import orm
from aiida_workgraph import task

from .utils import parse_outcar_timing, rank_benchmark_results


@task.calcfunction
def extract_benchmark_metrics(
    retrieved: orm.FolderData,
    label: orm.Str,
    ncore: orm.Int,
    kpar: orm.Int,
) -> orm.Dict:
    """Extract timing and memory metrics from VASP calculation.

    Parses the OUTCAR file from a completed VASP calculation to extract
    performance metrics for benchmarking parallelization parameters.

    Note:
        This function works even when the VASP calculation did not converge
        electronically (which is expected in benchmarking with limited NELM).
        It reads directly from the retrieved OUTCAR file, which is available
        regardless of convergence status, as long as VASP ran successfully.

    Args:
        retrieved: VASP retrieved FolderData containing OUTCAR.
        label: Label identifying this benchmark run (e.g., 'ncore4_kpar1').
        ncore: NCORE value used in this calculation.
        kpar: KPAR value used in this calculation.

    Returns:
        Dict containing:
        - 'label': Benchmark label
        - 'ncore': NCORE value
        - 'kpar': KPAR value
        - 'elapsed_time': Total elapsed time in seconds
        - 'cpu_time': Total CPU time in seconds
        - 'memory_kb': Maximum memory used in KB
        - 'memory_mb': Maximum memory used in MB (convenience)
        - 'loop_times': List of per-electronic-step timing
        - 'status': 'success' or 'failed'
        - 'error': Error message if failed
    """
    result = {
        "label": label.value,
        "ncore": ncore.value,
        "kpar": kpar.value,
        "elapsed_time": None,
        "cpu_time": None,
        "memory_kb": None,
        "memory_mb": None,
        "loop_times": [],
        "status": "failed",
        "error": None,
    }

    # Try to read OUTCAR
    try:
        if "OUTCAR" in retrieved.list_object_names():
            outcar_content = retrieved.get_object_content("OUTCAR")
        elif "OUTCAR.gz" in retrieved.list_object_names():
            import gzip

            outcar_bytes = retrieved.get_object_content("OUTCAR.gz", mode="rb")
            outcar_content = gzip.decompress(outcar_bytes).decode("utf-8")
        else:
            result["error"] = "OUTCAR not found in retrieved files"
            return orm.Dict(dict=result)
    except Exception as e:
        result["error"] = f"Error reading OUTCAR: {str(e)}"
        return orm.Dict(dict=result)

    # Parse timing information
    timing = parse_outcar_timing(outcar_content)

    if timing["parsing_success"]:
        result["elapsed_time"] = timing["elapsed_time"]
        result["cpu_time"] = timing["cpu_time"]
        result["memory_kb"] = timing["memory_kb"]
        if timing["memory_kb"]:
            result["memory_mb"] = timing["memory_kb"] / 1024.0
        result["loop_times"] = timing["loop_times"]
        result["status"] = "success"
    else:
        result["error"] = "Failed to parse timing from OUTCAR"

    return orm.Dict(dict=result)


@task.calcfunction
def compile_benchmark_results(
    weight_time: orm.Float = None,
    weight_memory: orm.Float = None,
    **metrics_dicts: orm.Dict,
) -> orm.Dict:
    """Compile and rank benchmark results from multiple runs.

    Aggregates timing/memory metrics from all benchmark calculations,
    ranks them by a weighted score of time and memory, and identifies
    the best parameters.

    Args:
        weight_time: Weight for elapsed time in ranking (0-1). Default 0.7.
        weight_memory: Weight for memory in ranking (0-1). Default 0.3.
        **metrics_dicts: Keyword arguments of orm.Dict with benchmark metrics.
            Keys are the benchmark labels (e.g., 'ncore4_kpar1').

    Returns:
        Dict containing:
        - 'all_results': List of all benchmark results, sorted by rank
        - 'best_by_time': Parameters with fastest elapsed time
        - 'best_by_memory': Parameters with lowest memory usage
        - 'recommended': Recommended parameters (best combined score)
        - 'num_successful': Number of successful benchmarks
        - 'num_failed': Number of failed benchmarks
        - 'weights': {'time': float, 'memory': float} used for ranking
    """
    # Set default weights
    wt = weight_time.value if weight_time is not None else 0.7
    wm = weight_memory.value if weight_memory is not None else 0.3

    # Collect all results
    all_results: list[dict[str, t.Any]] = []
    for key, metrics_node in metrics_dicts.items():
        if isinstance(metrics_node, orm.Dict):
            metrics = metrics_node.get_dict()
            all_results.append(metrics)

    # Count successes and failures
    successful = [r for r in all_results if r.get("status") == "success"]
    failed = [r for r in all_results if r.get("status") != "success"]

    # Initialize output
    output = {
        "all_results": [],
        "best_by_time": None,
        "best_by_memory": None,
        "recommended": None,
        "num_successful": len(successful),
        "num_failed": len(failed),
        "weights": {"time": wt, "memory": wm},
    }

    if not successful:
        output["all_results"] = all_results
        return orm.Dict(dict=output)

    # Rank results
    ranked = rank_benchmark_results(successful, weight_time=wt, weight_memory=wm)
    output["all_results"] = ranked + failed  # Successful first, then failed

    # Best by time
    by_time = sorted(successful, key=lambda x: x.get("elapsed_time", float("inf")))
    if by_time:
        best_time = by_time[0]
        output["best_by_time"] = {
            "ncore": best_time["ncore"],
            "kpar": best_time["kpar"],
            "elapsed_time": best_time["elapsed_time"],
            "label": best_time["label"],
        }

    # Best by memory
    valid_memory = [r for r in successful if r.get("memory_kb") is not None]
    if valid_memory:
        by_memory = sorted(valid_memory, key=lambda x: x.get("memory_kb", float("inf")))
        best_mem = by_memory[0]
        output["best_by_memory"] = {
            "ncore": best_mem["ncore"],
            "kpar": best_mem["kpar"],
            "memory_kb": best_mem["memory_kb"],
            "memory_mb": best_mem.get("memory_mb"),
            "label": best_mem["label"],
        }

    # Recommended (best combined score)
    if ranked:
        best = ranked[0]
        output["recommended"] = {
            "ncore": best["ncore"],
            "kpar": best["kpar"],
            "elapsed_time": best["elapsed_time"],
            "memory_kb": best.get("memory_kb"),
            "memory_mb": best.get("memory_mb"),
            "score": best.get("score"),
            "label": best["label"],
        }

    return orm.Dict(dict=output)


@task.calcfunction
def extract_recommended_parameters(benchmark_results: orm.Dict) -> orm.Dict:
    """Extract recommended NCORE/KPAR as a summary Dict.

    This calcfunction takes the compiled benchmark results and returns
    the key values in a single summary Dict for easy access.

    Args:
        benchmark_results: Dict from compile_benchmark_results containing
            'recommended', 'best_by_time', 'best_by_memory' keys.

    Returns:
        Dict with key metrics:
        - 'recommended_ncore': recommended NCORE value
        - 'recommended_kpar': recommended KPAR value
        - 'recommended_elapsed_time_s': best elapsed time in seconds
        - 'recommended_memory_mb': memory usage in MB
        - 'num_successful': number of successful benchmarks
        - 'num_failed': number of failed benchmarks
        - 'fastest_*': info about fastest configuration
        - 'lowest_memory_*': info about lowest memory configuration
    """
    results = benchmark_results.get_dict()

    # Get recommended values (default to 1 if not found)
    recommended = results.get("recommended", {})
    ncore = recommended.get("ncore", 1)
    kpar = recommended.get("kpar", 1)
    elapsed_time = recommended.get("elapsed_time")

    # Get counts
    num_successful = results.get("num_successful", 0)
    num_failed = results.get("num_failed", 0)

    # Build summary
    summary = {
        "recommended_ncore": ncore,
        "recommended_kpar": kpar,
        "recommended_elapsed_time_s": elapsed_time,
        "recommended_memory_mb": recommended.get("memory_mb"),
        "num_benchmarks_total": num_successful + num_failed,
        "num_successful": num_successful,
        "num_failed": num_failed,
    }

    # Add best by time info
    best_time = results.get("best_by_time", {})
    if best_time:
        summary["fastest_ncore"] = best_time.get("ncore")
        summary["fastest_kpar"] = best_time.get("kpar")
        summary["fastest_time_s"] = best_time.get("elapsed_time")

    # Add best by memory info
    best_mem = results.get("best_by_memory", {})
    if best_mem:
        summary["lowest_memory_ncore"] = best_mem.get("ncore")
        summary["lowest_memory_kpar"] = best_mem.get("kpar")
        summary["lowest_memory_mb"] = best_mem.get("memory_mb")

    return orm.Dict(dict=summary)
