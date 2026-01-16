"""Fukui function calculation module for PS-TEROS.

This module implements Fukui function calculations via the interpolation method
(Phase 1), Fukui potential via the electrodes method (Phase 2), and perturbative
expansion for interaction energy prediction (Phase 4).

Phase 1 (Interpolation): Generates CHGCAR files at different charge states
and computes the Fukui function via polynomial interpolation using FukuiGrid.

Phase 2 (Electrodes): Computes the Fukui potential by applying electrostatic
corrections for electrode/surface systems using the dielectric constant
obtained from DFPT calculations.

Phase 4 (Perturbative Expansion): Computes the interaction energy map using the
c-DFT perturbative expansion: ΔU(r) = q·Φ(r) - q·ΔN·vf±(r). This predicts
favorable/unfavorable adsorption sites for charged adsorbates.

Public API:
    build_fukui_workgraph: Main workflow builder (supports all phases)
    get_fukui_results: Extract results from completed workflow
    print_fukui_summary: Print formatted summary of results
    FukuiCalculationScatter: @task.graph for parallel charge calculations
    extract_dielectric_constant: @task.calcfunction to extract epsilon from DFPT
    run_fukui_electrodes_calcfunc: @task.calcfunction for electrodes method
    run_perturbative_expansion_calcfunc: @task.calcfunction for perturbative expansion
    extract_locpot_from_retrieved: @task.calcfunction to extract LOCPOT from VASP

Example (Phase 1 only):
    >>> from teros.core.fukui import build_fukui_workgraph, get_fukui_results
    >>> wg = build_fukui_workgraph(
    ...     structure=my_structure,
    ...     nelect_neutral=192,
    ...     code_label='VASP-6.5.1@localwork',
    ...     builder_inputs={...},
    ...     fukui_type='plus',
    ...     compute_fukui=True,
    ... )
    >>> wg.submit(wait=False)

Example (Phase 1 + Phase 2):
    >>> wg = build_fukui_workgraph(
    ...     structure=slab_structure,
    ...     nelect_neutral=192,
    ...     code_label='VASP-6.5.1@localwork',
    ...     builder_inputs={...},
    ...     fukui_type='plus',
    ...     compute_fukui=True,
    ...     compute_fukui_potential=True,  # Enable Phase 2
    ...     bulk_structure=bulk_structure,  # Required for DFPT
    ... )

Example (Full workflow with Phase 4):
    >>> wg = build_fukui_workgraph(
    ...     structure=slab_structure,
    ...     nelect_neutral=192,
    ...     code_label='VASP-6.5.1@localwork',
    ...     builder_inputs={...},
    ...     fukui_type='plus',
    ...     compute_fukui=True,
    ...     compute_fukui_potential=True,
    ...     bulk_structure=bulk_structure,
    ...     compute_perturbative_expansion=True,  # Enable Phase 4
    ...     probe_charge=0.3,        # Charge of adsorbate in |e|
    ...     electron_transfer=-0.3,  # Electron donation to surface
    ... )
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
    # Phase 2: Electrodes method
    extract_dielectric_constant,
    run_fukui_electrodes_calcfunc,
    # Phase 4: Perturbative expansion
    extract_locpot_from_retrieved,
    run_perturbative_expansion_calcfunc,
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

    # Phase 1 Calcfunctions (Interpolation)
    'collect_chgcar_files',
    'generate_fukui_summary',
    'extract_total_energy',
    'collect_chgcar_files_internal',
    'generate_fukui_summary_internal',
    'run_fukui_interpolation',
    'run_fukui_interpolation_calcfunc',

    # Phase 2 Calcfunctions (Electrodes)
    'extract_dielectric_constant',
    'run_fukui_electrodes_calcfunc',

    # Phase 4 Calcfunctions (Perturbative Expansion)
    'extract_locpot_from_retrieved',
    'run_perturbative_expansion_calcfunc',

    # Utilities
    'make_delta_label',
    'validate_fukui_inputs',
    'DEFAULT_DELTA_N_VALUES',
    'calculate_nelect',
    'print_nelect_breakdown',
]
