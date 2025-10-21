"""Main WorkGraph for surface hydroxylation workflow."""

import typing as t
from aiida import orm
from aiida_workgraph import task, namespace, dynamic

from .tasks import generate_structures, collect_results, package_namespace_to_dict
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
        builder_config: Complete VASP builder configuration (dict with keys):
            - code: AiiDA Code for VASP
            - parameters: Dict with INCAR parameters
            - potential_family: Pseudopotential family name
            - potential_mapping: Dict mapping elements to potentials
            - kpoints_spacing: K-points spacing (optional)
            - options: Scheduler options dict
            - clean_workdir: bool
            - settings: Parser settings dict (optional)
        max_parallel_jobs: Maximum number of concurrent VASP relaxations (default: 2)

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
        ...     'parameters': {'EDIFF': 1e-5, 'ENCUT': 520},
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
    # Convert inputs to AiiDA nodes if needed
    if not isinstance(surface_params, orm.Dict):
        surface_params = orm.Dict(dict=surface_params)

    # Task 1: Generate surface structure variants
    # Returns: {manifest: Dict, structure_0: StructureData, structure_1: ..., structure_N: ...}
    gen_task = generate_structures(
        structure=structure,
        params=surface_params,
    )

    # Task 2: Relax all generated structures in parallel with semaphore
    # Input: Need to extract structure_* outputs from gen_task and pass as dict
    # In WorkGraph's @task.graph, we can directly access the returned dict items
    # generate_structures returns a dict with manifest and structure_N keys
    # We filter to get only the structure keys and pass them

    # Extract structures from gen_task outputs
    # gen_task returns a dict, so we can access items like: gen_task['structure_0']
    # However, in @task.graph context, we need to build the structures dict at graph definition time
    # The actual filtering happens at runtime through the namespace

    # Pass the full gen_task namespace (contains manifest + structure_0, structure_1, ...)
    # relax_slabs_with_semaphore will filter out the structure_* keys
    relax_outputs = relax_slabs_with_semaphore(
        structures=gen_task.namespace,
        builder_config=builder_config,
        max_parallel=max_parallel_jobs,
    )

    # Task 3: Package namespace outputs into Dict nodes
    # relax_slabs_with_semaphore returns namespaces with dynamic outputs
    # We need to convert them to Dict nodes for collect_results

    # Package each namespace into a Dict
    structures_dict = package_namespace_to_dict(relax_outputs.structures)
    energies_dict = package_namespace_to_dict(relax_outputs.energies)
    exit_statuses_dict = package_namespace_to_dict(relax_outputs.exit_statuses)
    errors_dict = package_namespace_to_dict(relax_outputs.errors)

    # Task 4: Collect and organize results
    collect_task = collect_results(
        manifest=gen_task.manifest,
        structures=structures_dict.result,
        energies=energies_dict.result,
        exit_statuses=exit_statuses_dict.result,
        errors=errors_dict.result,
    )

    # Return outputs
    return {
        'manifest': gen_task.manifest,
        'successful_relaxations': collect_task.successful_relaxations,
        'failed_relaxations': collect_task.failed_relaxations,
        'statistics': collect_task.statistics,
    }
