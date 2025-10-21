"""Main WorkGraph for surface hydroxylation workflow."""

import typing as t
from aiida import orm
from aiida_workgraph import task, namespace, dynamic, WorkGraph

from .tasks import (
    generate_structures,
    extract_manifest,
    extract_successful_relaxations,
    extract_failed_relaxations,
    extract_statistics,
)
from .relaxations import relax_slabs_with_semaphore


@task.graph(outputs=[
    'manifest',
    'structures',
    'energies',
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
        vasp_config: VASP configuration (dict) for vasp.v2.vasp:
            - parameters: Dict with INCAR parameters (PREC, ENCUT, EDIFF, ISIF, NSW, IBRION, EDIFFG, etc.)
            - kpoints_spacing: Float (Angstrom^-1, typically 0.3-0.5)
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
        ...     'parameters': {
        ...         'PREC': 'Accurate',
        ...         'ENCUT': 520,
        ...         'EDIFF': 1e-6,
        ...         'ISIF': 2,
        ...         'NSW': 200,
        ...         'IBRION': 2,
        ...         'EDIFFG': -0.02,
        ...     },
        ...     'kpoints_spacing': 0.3,
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
    # Returns namespace with manifest (Dict) and structures (dynamic dict of StructureData)
    gen_outputs = generate_structures(
        structure=structure,
        params=surface_params,
    )

    # Task 2: Relax all generated structures in parallel
    # Pass structures dict from generate_structures
    relax_outputs = relax_slabs_with_semaphore(
        structures=gen_outputs.structures,
        code_pk=code.pk,
        vasp_config=vasp_config,
        options=options,
        max_parallel=max_parallel_jobs,
    )

    # Return raw outputs
    # Note: Result organization (successful/failed/statistics) must be done
    # in post-processing due to WorkGraph namespace serialization limitations.
    # Use the organize_hydroxylation_results() helper function after workflow completion.
    return {
        'manifest': gen_outputs.manifest,
        'structures': relax_outputs.structures,
        'energies': relax_outputs.energies,
    }


def build_surface_hydroxylation_workgraph(
    structure: orm.StructureData = None,
    structure_pk: int = None,
    surface_params: dict = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    vasp_config: dict = None,
    options: dict = None,
    max_parallel_jobs: int = 2,
    name: str = 'SurfaceHydroxylation',
) -> WorkGraph:
    """
    Build a WorkGraph for surface hydroxylation/vacancy calculations.

    This builder function creates a ready-to-submit WorkGraph following the same
    architectural pattern as build_core_workgraph() in teros.core.workgraph.

    The workflow performs:
    1. Generation of surface structure variants (hydroxylation/vacancies)
    2. Parallel VASP relaxation of all variants with batch control
    3. Collection and organization of results with success/failure tracking

    Args:
        structure: Input relaxed slab structure (StructureData).
                  Either this or structure_pk must be provided.
        structure_pk: PK of relaxed slab structure.
                     Either this or structure must be provided.
        surface_params: Parameters for surface modification (dict):
            - mode: str ('vacancies'/'hydrogen'/'combine')
            - species: str (default 'O')
            - z_window: float (default 0.5 Angstrom)
            - which_surface: str ('top'/'bottom'/'both')
            - oh_dist: float (default 0.98 Angstrom, O-H bond distance)
            - include_empty: bool (default False)
            - supercell: list[int] or None
            - deduplicate_by_coverage: bool (default True)
            - coverage_bins: int or None
        code_label: Label of VASP code (e.g., 'VASP-VTST-6.4.3@bohr')
        vasp_config: VASP configuration (dict) for vasp.v2.vasp workflow:
            - parameters: Dict with INCAR parameters
              * PREC: str ('Normal'/'Accurate')
              * ENCUT: float (plane-wave cutoff in eV)
              * EDIFF: float (electronic convergence in eV)
              * ISMEAR: int (smearing method)
              * SIGMA: float (smearing width in eV)
              * ALGO: str (electronic minimization algorithm)
              * LREAL: bool or str (real-space projection)
              * NELM: int (max electronic steps)
              * LWAVE: bool (write WAVECAR)
              * LCHARG: bool (write CHGCAR)
              * ISIF: int (relaxation type: 2=positions only, 3=positions+cell)
              * NSW: int (max ionic steps)
              * IBRION: int (ionic relaxation: 2=CG, 1=RMM-DIIS)
              * EDIFFG: float (force convergence in eV/Å, negative = force criterion)
            - kpoints_spacing: Float (Angstrom^-1, typically 0.3-0.5)
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - potential_mapping: Dict (optional, element->potential mapping)
            - clean_workdir: Bool (default False)
        options: Scheduler options (dict):
            - resources: Dict with 'num_machines' and 'num_mpiprocs_per_machine'
            - queue_name: str (queue name)
            - max_wallclock_seconds: int (walltime in seconds)
            - account: str (optional, cluster account)
            - prepend_text: str (optional, commands before VASP)
            - custom_scheduler_commands: str (optional, scheduler directives)
        max_parallel_jobs: Number of structures to process in this run (default 2).
                          Uses simple batch approach: processes first N structures only.
                          Increase this value in subsequent runs to process more structures.
        name: Name for the workflow (default 'SurfaceHydroxylation')

    Returns:
        WorkGraph: Ready-to-submit workflow instance

    Raises:
        ValueError: If neither structure nor structure_pk is provided
        ValueError: If code_label cannot be loaded

    Example:
        >>> from aiida import orm
        >>> from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph
        >>>
        >>> # Define parameters
        >>> surface_params = {
        ...     'mode': 'hydrogen',
        ...     'species': 'O',
        ...     'coverage_bins': 5,
        ... }
        >>> vasp_config = {
        ...     'parameters': {
        ...         'PREC': 'Accurate',
        ...         'ENCUT': 520,
        ...         'EDIFF': 1e-6,
        ...         'ISIF': 2,
        ...         'NSW': 200,
        ...         'IBRION': 2,
        ...         'EDIFFG': -0.02,
        ...     },
        ...     'kpoints_spacing': 0.3,
        ...     'potential_family': 'PBE',
        ... }
        >>> options = {
        ...     'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
        ...     'queue_name': 'normal',
        ...     'max_wallclock_seconds': 3600 * 10,
        ... }
        >>>
        >>> # Build and submit workflow
        >>> wg = build_surface_hydroxylation_workgraph(
        ...     structure_pk=1234,
        ...     surface_params=surface_params,
        ...     code_label='vasp@cluster',
        ...     vasp_config=vasp_config,
        ...     options=options,
        ...     max_parallel_jobs=3,
        ... )
        >>> wg.submit()
    """
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================

    # Validate structure input
    if structure is None and structure_pk is None:
        raise ValueError(
            "Either 'structure' or 'structure_pk' must be provided.\n"
            "Example: structure_pk=1234 or structure=orm.load_node(1234)"
        )

    # Load structure from PK if needed
    if structure is None:
        try:
            structure = orm.load_node(structure_pk)
            if not isinstance(structure, orm.StructureData):
                raise ValueError(
                    f"Node {structure_pk} is not a StructureData.\n"
                    f"Found: {type(structure).__name__}"
                )
        except Exception as e:
            raise ValueError(
                f"Could not load structure from PK {structure_pk}.\n"
                f"Error: {e}\n"
                f"Use 'verdi data core.structure list' to see available structures."
            )

    # Load code
    try:
        code = orm.load_code(code_label)
    except Exception as e:
        raise ValueError(
            f"Could not load VASP code '{code_label}'.\n"
            f"Error: {e}\n\n"
            f"Available codes:\n"
            f"  Use 'verdi code list' to see configured codes.\n"
            f"  Or setup VASP code with: verdi code create core.code.installed"
        )

    # Set defaults for optional parameters
    if surface_params is None:
        surface_params = {
            'mode': 'hydrogen',
            'species': 'O',
            'z_window': 0.5,
            'which_surface': 'top',
            'oh_dist': 0.98,
            'include_empty': False,
            'deduplicate_by_coverage': True,
            'coverage_bins': 10,
        }

    if vasp_config is None:
        # Default lightweight VASP configuration
        vasp_config = {
            'parameters': {
                'PREC': 'Normal',
                'ENCUT': 400,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': False,
                'NELM': 100,
                'LWAVE': False,
                'LCHARG': False,
                'ISIF': 2,
                'NSW': 100,
                'IBRION': 2,
                'EDIFFG': -0.02,
            },
            'kpoints_spacing': 0.5,
            'potential_family': 'PBE',
            'potential_mapping': {},
            'clean_workdir': False,
        }

    if options is None:
        # Default scheduler options
        options = {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 16,
            },
            'queue_name': 'normal',
            'max_wallclock_seconds': 3600,
        }

    # ========================================================================
    # BUILD WORKGRAPH
    # ========================================================================

    # Build the WorkGraph using the @task.graph function
    wg = SurfaceHydroxylationWorkGraph.build(
        structure=structure,
        surface_params=surface_params,
        code=code,
        vasp_config=vasp_config,
        options=options,
        max_parallel_jobs=max_parallel_jobs,
    )

    # Set the workflow name
    wg.name = name

    return wg


def organize_hydroxylation_results(workflow_node):
    """
    Organize hydroxylation workflow results into successful/failed/statistics.

    This helper function processes the raw outputs from a completed
    SurfaceHydroxylation workflow and organizes them into structured results.

    Args:
        workflow_node: Completed WorkGraph node (orm.WorkChainNode)

    Returns:
        dict with:
            - successful_relaxations: list of dicts with successful results
            - failed_relaxations: list of dicts with failed results
            - statistics: dict with total/succeeded/failed counts

    Example:
        >>> from aiida import orm
        >>> from teros.core.surface_hydroxylation import organize_hydroxylation_results
        >>>
        >>> node = orm.load_node(12345)  # Completed workflow
        >>> results = organize_hydroxylation_results(node)
        >>> print(f"Succeeded: {results['statistics']['succeeded']}")
        >>> print(f"Failed: {results['statistics']['failed']}")
    """
    # Get outputs
    manifest = workflow_node.outputs.manifest.get_dict()
    structures = workflow_node.outputs.structures
    energies = workflow_node.outputs.energies

    variants = manifest['variants']
    successful = []
    failed = []

    # Process each variant
    for idx, variant in enumerate(variants):
        # Construct the same descriptive key format used in outputs
        # Format: "idx_variantname" (e.g., "0_oh_000_3_7572")
        # Replace dots with underscores for AiiDA link label compatibility
        variant_name_safe = variant['name'].replace('.', '_')
        key = f"{idx}_{variant_name_safe}"

        # Check if relaxation succeeded (has outputs)
        try:
            structure_node = structures[key]
            energy_node = energies[key]

            successful.append({
                'name': variant['name'],
                'structure_pk': structure_node.pk,
                'energy': energy_node.value,
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'metadata': variant
            })
        except (KeyError, AttributeError):
            # Missing from outputs - relaxation failed
            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'error_message': 'Relaxation failed - no structure/energy output'
            })

    return {
        'successful_relaxations': successful,
        'failed_relaxations': failed,
        'statistics': {
            'total': len(variants),
            'succeeded': len(successful),
            'failed': len(failed)
        }
    }
