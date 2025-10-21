"""Main WorkGraph for surface hydroxylation workflow."""

import typing as t
from aiida import orm
from aiida_workgraph import task, namespace, dynamic

from .tasks import (
    generate_structures,
    collect_results,
    extract_manifest,
    extract_successful_relaxations,
    extract_failed_relaxations,
    extract_statistics,
)
from .relaxations import relax_slabs_with_semaphore


@task.graph(outputs=[
    'manifest',
    'successful_relaxations',
    'failed_relaxations',
    'statistics',
])
def SurfaceHydroxylationWorkGraph(
    structure: orm.StructureData,
    surface_params: dict,
    builder_config: dict,
    max_parallel_jobs: int = 2,
) -> dict:
    """
    Main workflow for surface hydroxylation studies.

    This WorkGraph orchestrates the complete hydroxylation workflow:
    1. Generate surface structure variants (hydroxylation/vacancies) from input slab
    2. Relax all variants in parallel with semaphore-based concurrency control
    3. Collect and organize results

    Args:
        structure: Input relaxed slab structure (StructureData)
        surface_params: Parameters for surface_modes.py (dict with keys):
            - mode: str ('vacancies'/'hydrogen'/'combine')
            - species: str (default 'O')
            - z_window: float (default 0.5)
            - which_surface: str ('top'/'bottom'/'both')
            - oh_dist: float (default 0.98)
            - include_empty: bool (default False)
            - supercell: list[int] or None
            - deduplicate_by_coverage: bool
            - coverage_bins: int or None
        builder_config: Complete VASP relaxation configuration for vasp.v2.relax (dict):
            - code: AiiDA Code for VASP
            - relax: Dict with relaxation settings (positions, shape, volume, force_cutoff, steps, algo)
            - base: Dict with INCAR parameters (PREC, ENCUT, EDIFF, etc.)
            - kpoints_distance: Float (Angstrom^-1, typically 0.3-0.5)
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - potential_mapping: Dict (optional, element->potential)
            - options: Dict with scheduler settings (resources, queue, time)
            - clean_workdir: Bool (default: False)
        max_parallel_jobs: Number of structures to process in this run (default: 2)
                          Uses simple batch approach: processes first N structures only

    Returns:
        Dictionary with outputs:
            - manifest: Dict with surface_modes metadata
            - successful_relaxations: Dict with list of successful results
                Each result: {name, structure_pk, energy, coverage, metadata}
            - failed_relaxations: Dict with list of failed results
                Each result: {name, coverage, exit_status, error_message}
            - statistics: Dict with {total, succeeded, failed} counts

    Example:
        >>> from aiida import orm
        >>> from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph
        >>>
        >>> # Load relaxed slab structure
        >>> structure = orm.load_node(1234)
        >>>
        >>> # Define parameters
        >>> surface_params = {
        ...     'mode': 'hydrogen',
        ...     'species': 'O',
        ...     'coverage_bins': 5,
        ... }
        >>> builder_config = {
        ...     'code': orm.load_code('vasp@cluster'),
        ...     'relax': {'positions': True, 'force_cutoff': 0.02, 'steps': 200},
        ...     'base': {'EDIFF': 1e-6, 'ENCUT': 520, 'PREC': 'Accurate'},
        ...     'kpoints_distance': 0.3,
        ...     'potential_family': 'PBE',
        ...     'options': {'resources': {'num_machines': 1}},
        ... }
        >>>
        >>> # Create and submit workflow
        >>> wg = SurfaceHydroxylationWorkGraph(
        ...     structure=structure,
        ...     surface_params=surface_params,
        ...     builder_config=builder_config,
        ...     max_parallel_jobs=3,
        ... )
        >>> wg.submit()
    """
    # Validate builder_config (for vasp.v2.relax)
    required_keys = ['code', 'options']
    missing_keys = [key for key in required_keys if key not in builder_config]
    if missing_keys:
        raise ValueError(f"builder_config missing required keys: {missing_keys}")

    # Warn if using old format
    if 'parameters' in builder_config:
        print("WARNING: builder_config uses old 'parameters' key.")
        print("For vasp.v2.relax, use 'base' for INCAR parameters and 'relax' for relaxation settings.")
        print("See examples/surface_hydroxylation/run_production.py for correct format.")

    # Convert inputs to AiiDA nodes if needed
    if not isinstance(surface_params, orm.Dict):
        surface_params = orm.Dict(dict=surface_params)

    # Task 1: Generate surface structure variants
    # Wrap @calcfunction with task() for use in WorkGraph
    # Returns: {manifest: Dict, structure_0: StructureData, structure_1: ..., structure_N: ...}
    gen_task = task(generate_structures)(
        structure=structure,
        params=surface_params,
    )

    # Task 1b: Extract manifest from result dict
    manifest_task = task(extract_manifest)(result=gen_task.result)

    # Task 2: Prepare VASP configuration
    # Extract code and convert other components to AiiDA Dict nodes
    code = builder_config['code']
    options = orm.Dict(dict=builder_config.get('options', {}))

    # VASP config without code and options
    vasp_config_dict = {
        'relax': builder_config.get('relax', {}),
        'base': builder_config.get('base', {}),
        'kpoints_distance': builder_config.get('kpoints_distance', 0.5),
        'potential_family': builder_config.get('potential_family', 'PBE'),
        'potential_mapping': builder_config.get('potential_mapping', {}),
        'clean_workdir': builder_config.get('clean_workdir', False),
    }
    vasp_config = orm.Dict(dict=vasp_config_dict)

    # Task 3: Relax all generated structures in parallel
    # Pass the full result dict from gen_task
    # relax_slabs_with_semaphore will extract structure_* keys
    relax_outputs = relax_slabs_with_semaphore(
        structures=gen_task.result,
        code=code,
        vasp_config=vasp_config,
        options=options,
        max_parallel=max_parallel_jobs,
    )

    # Task 4: Collect and organize results
    # Wrap collect_results with task() to use it in the WorkGraph
    collect_task = task(collect_results)(
        manifest=manifest_task.result,
        structures=relax_outputs.structures,
        energies=relax_outputs.energies,
        exit_statuses=relax_outputs.exit_statuses,
        errors=relax_outputs.errors,
    )

    # Task 5: Extract individual outputs from collect_results result
    successful_task = task(extract_successful_relaxations)(result=collect_task.result)
    failed_task = task(extract_failed_relaxations)(result=collect_task.result)
    statistics_task = task(extract_statistics)(result=collect_task.result)

    # Return outputs
    return {
        'manifest': manifest_task.result,
        'successful_relaxations': successful_task.result,
        'failed_relaxations': failed_task.result,
        'statistics': statistics_task.result,
    }
