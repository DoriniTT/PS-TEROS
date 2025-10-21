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
    code: orm.InstalledCode,
    vasp_config: dict,
    options: dict,
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
        code: AiiDA Code for VASP (InstalledCode)
        vasp_config: VASP relaxation configuration (dict) for vasp.v2.relax:
            - relax: Dict with relaxation settings (positions, shape, volume, force_cutoff, steps, algo)
            - base: Dict with INCAR parameters (PREC, ENCUT, EDIFF, etc.)
            - kpoints_distance: Float (Angstrom^-1, typically 0.3-0.5)
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - potential_mapping: Dict (optional, element->potential)
            - clean_workdir: Bool (default: False)
        options: Scheduler options (dict) with resources, queue, time, etc.
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
        >>> code = orm.load_code('vasp@cluster')
        >>>
        >>> # Define parameters
        >>> surface_params = {
        ...     'mode': 'hydrogen',
        ...     'species': 'O',
        ...     'coverage_bins': 5,
        ... }
        >>> vasp_config = {
        ...     'relax': {'positions': True, 'force_cutoff': 0.02, 'steps': 200},
        ...     'base': {'EDIFF': 1e-6, 'ENCUT': 520, 'PREC': 'Accurate'},
        ...     'kpoints_distance': 0.3,
        ...     'potential_family': 'PBE',
        ...     'clean_workdir': False,
        ... }
        >>> options = {'resources': {'num_machines': 1}}
        >>>
        >>> # Create and submit workflow
        >>> wg = SurfaceHydroxylationWorkGraph(
        ...     structure=structure,
        ...     surface_params=surface_params,
        ...     code=code,
        ...     vasp_config=vasp_config,
        ...     options=options,
        ...     max_parallel_jobs=3,
        ... )
        >>> wg.submit()
    """

    # Convert inputs to AiiDA nodes if needed
    if not isinstance(surface_params, orm.Dict):
        surface_params = orm.Dict(dict=surface_params)

    # Task 1: Generate surface structure variants
    # Use task() wrapper to create task connections, not execute immediately
    # Returns: {manifest: Dict, structure_0: StructureData, structure_1: ..., structure_N: ...}
    from aiida_workgraph import task as wgtask
    gen_task = wgtask(generate_structures)(
        structure=structure,
        params=surface_params,
    )

    # Task 1b: Extract manifest from result dict
    manifest_task = wgtask(extract_manifest)(result=gen_task.outputs['result'])

    # Task 2: Relax all generated structures in parallel
    # Pass the full result dict from gen_task
    # relax_slabs_with_semaphore will extract structure_* keys
    # Pass plain values (int/dict) instead of AiiDA nodes to avoid serialization issues
    relax_outputs = relax_slabs_with_semaphore(
        structures=gen_task.outputs['result'],
        code_pk=code.pk,
        vasp_config=vasp_config,
        options=options,
        max_parallel=max_parallel_jobs,
    )

    # Task 4: Collect and organize results
    # Use task() wrapper to create task connections
    collect_task = wgtask(collect_results)(
        manifest=manifest_task.outputs['result'],
        structures=relax_outputs.structures,
        energies=relax_outputs.energies,
        exit_statuses=relax_outputs.exit_statuses,
        errors=relax_outputs.errors,
    )

    # Task 5: Extract individual outputs from collect_results result
    successful_task = wgtask(extract_successful_relaxations)(result=collect_task.outputs['result'])
    failed_task = wgtask(extract_failed_relaxations)(result=collect_task.outputs['result'])
    statistics_task = wgtask(extract_statistics)(result=collect_task.outputs['result'])

    # Return outputs - reference task outputs, not actual values
    return {
        'manifest': manifest_task.outputs['result'],
        'successful_relaxations': successful_task.outputs['result'],
        'failed_relaxations': failed_task.outputs['result'],
        'statistics': statistics_task.outputs['result'],
    }
