"""VASP parallelization benchmarking module.

This module provides tools to benchmark and optimize VASP parallelization
parameters (NCORE, KPAR) for a given structure and computational resources.

Main entry point:
    build_parallelization_benchmark_workgraph(): Create benchmark workflow

Helper functions:
    get_benchmark_results(): Extract results from completed workflow
    print_benchmark_summary(): Print formatted results table

Utility functions:
    generate_ncore_values(): Generate NCORE values for testing
    generate_kpar_values(): Generate KPAR values for testing
    generate_benchmark_combinations(): Create NCORE/KPAR test matrix
    estimate_kpoints_count(): Estimate k-points from structure
    parse_outcar_timing(): Parse timing from VASP OUTCAR

Example:
    >>> from teros.core.vasp_parallelization import (
    ...     build_parallelization_benchmark_workgraph,
    ...     print_benchmark_summary,
    ... )
    >>>
    >>> wg = build_parallelization_benchmark_workgraph(
    ...     structure=my_structure,
    ...     code_label='VASP-6.5.1@localwork',
    ...     potential_family='PBE',
    ...     potential_mapping={'Si': 'Si'},
    ...     num_procs=8,
    ...     max_concurrent_jobs=2,
    ... )
    >>> wg.submit()
    >>>
    >>> # After completion:
    >>> print_benchmark_summary(wg.pk)
"""

from .workgraph import (
    build_parallelization_benchmark_workgraph,
    get_benchmark_results,
    get_benchmark_summary,
    print_benchmark_summary,
)
from .tasks import (
    compile_benchmark_results,
    extract_benchmark_metrics,
    extract_recommended_parameters,
)
from .utils import (
    estimate_kpoints_count,
    generate_benchmark_combinations,
    generate_kpar_values,
    generate_ncore_values,
    get_divisors,
    parse_outcar_timing,
    prepare_benchmark_incar,
    rank_benchmark_results,
)

__all__ = [
    # Main entry point
    "build_parallelization_benchmark_workgraph",
    # Helper functions
    "get_benchmark_results",
    "print_benchmark_summary",
    # Task functions
    "extract_benchmark_metrics",
    "compile_benchmark_results",
    "extract_recommended_parameters",
    # Utility functions
    "generate_ncore_values",
    "generate_kpar_values",
    "generate_benchmark_combinations",
    "estimate_kpoints_count",
    "parse_outcar_timing",
    "prepare_benchmark_incar",
    "rank_benchmark_results",
    "get_divisors",
]
