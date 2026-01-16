"""Fukui function calculation module for PS-TEROS.

This module implements Fukui function calculations via the interpolation method.
It generates CHGCAR files at different charge states for later processing with
FukuiGrid.py.

Public API:
    build_fukui_workgraph: Main workflow builder
    get_fukui_results: Extract results from completed workflow
    print_fukui_summary: Print formatted summary of results
    FukuiCalculationScatter: @task.graph for parallel charge calculations
    collect_chgcar_files: @task.calcfunction to collect CHGCAR files
    generate_fukui_summary: @task.calcfunction to generate summary dict

Example:
    >>> from teros.core.fukui import build_fukui_workgraph, get_fukui_results
    >>> wg = build_fukui_workgraph(
    ...     structure=my_structure,
    ...     nelect_neutral=192,
    ...     code_label='VASP-6.5.1@localwork',
    ...     builder_inputs={...},
    ...     fukui_type='plus',
    ... )
    >>> wg.submit(wait=False)
    >>> # After completion:
    >>> results = get_fukui_results(wg.pk)
    >>> print(results['file_names'])
    ['CHGCAR_0.00', 'CHGCAR_0.05', 'CHGCAR_0.10', 'CHGCAR_0.15']
"""

from .workgraph import (
    build_fukui_workgraph,
    get_fukui_results,
    print_fukui_summary,
    FukuiCalculationScatter,
)

from .tasks import (
    collect_chgcar_files,
    generate_fukui_summary,
    extract_total_energy,
    collect_chgcar_files_internal,
    generate_fukui_summary_internal,
    run_fukui_interpolation,
    run_fukui_interpolation_calcfunc,
)

from .utils import (
    make_delta_label,
    validate_fukui_inputs,
    DEFAULT_DELTA_N_VALUES,
    calculate_nelect,
    print_nelect_breakdown,
)

__all__ = [
    # Main workflow
    'build_fukui_workgraph',
    'get_fukui_results',
    'print_fukui_summary',

    # Task graphs
    'FukuiCalculationScatter',

    # Calcfunctions
    'collect_chgcar_files',
    'generate_fukui_summary',
    'extract_total_energy',
    'collect_chgcar_files_internal',
    'generate_fukui_summary_internal',
    'run_fukui_interpolation',
    'run_fukui_interpolation_calcfunc',

    # Utilities
    'make_delta_label',
    'validate_fukui_inputs',
    'DEFAULT_DELTA_N_VALUES',
    'calculate_nelect',
    'print_nelect_breakdown',
]
