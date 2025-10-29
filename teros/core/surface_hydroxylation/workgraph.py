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
from .surface_energy_workgraph import create_surface_energy_task, calculate_all_surface_energies


@task.graph(outputs=[
    'manifest',
    'structures',
    'energies',
    'bulk_structure',
    'bulk_energy',
    'pristine_structure',
    'pristine_energy',
])
def SurfaceHydroxylationWorkGraph(
    structure: orm.StructureData,
    surface_params: dict,
    code: orm.InstalledCode,
    builder_inputs: dict,
    bulk_structure: orm.StructureData,
    bulk_builder_inputs: dict,
    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    structure_specific_builder_inputs: dict = None,
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
        builder_inputs: Complete builder configuration (dict) for vasp.v2.vasp WorkChain.
            This dict should contain ALL builder parameters as you would set them directly:
            - parameters: Dict with nested 'incar' dict containing VASP tags
            - kpoints_spacing: Float (or provide 'kpoints' KpointsData)
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - potential_mapping: Dict or orm.Dict (element->potential mapping)
            - options: Dict or orm.Dict (scheduler options: resources, queue, walltime)
            - settings: Dict or orm.Dict (parser settings, optional)
            - dynamics: Dict (selective dynamics, optional - will be added automatically if fix_type is set)
            - clean_workdir: Bool (default: False)
            Note: 'code' and 'structure' will be set automatically by the workflow
        bulk_structure: Bulk crystal structure for reference calculations (StructureData).
                       Required for surface energy calculations per Section S2.
        bulk_builder_inputs: VASP parameters for bulk relaxation (dict). Must include
                            ISIF=3 for cell relaxation. Same format as builder_inputs.
        max_parallel_jobs: Number of structures to process in this run (default: 2)
                          Uses simple batch approach: processes first N structures only
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None). Default: None (no fixing)
        fix_thickness: Thickness in Angstroms for fixing region. Default: 0.0
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O']).
                     If None, all elements in the region are fixed. Default: None
        structure_specific_builder_inputs: Optional dict to override builder_inputs for specific
                     structure indices. Keys are integer indices (0, 1, 2, ...) corresponding to
                     generated structures. Values are PARTIAL builder_inputs dicts that will be
                     MERGED into the default builder_inputs. You only need to specify the parameters
                     you want to override.
                     Example: {0: {'parameters': {'incar': {'ALGO': 'Normal'}}},
                              2: {'kpoints_spacing': 0.2, 'options': {'resources': {...}}}}
                     Structure 0 will use ALGO='Normal' with all other default parameters.
                     Structure 2 will use kpoints_spacing=0.2 and custom resources with all other defaults.
                     Structures 1, 3, etc. will use the complete default builder_inputs.
                     Default: None (use default builder_inputs for all)

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

    # =========================================================================
    # Task 0: Bulk Relaxation (NEW)
    # =========================================================================
    from aiida.plugins import WorkflowFactory

    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Prepare settings with parser configuration
    bulk_settings = bulk_builder_inputs.get('settings', {})
    if not bulk_settings:
        bulk_settings = {
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }
    if not isinstance(bulk_settings, orm.Dict):
        bulk_settings = orm.Dict(dict=bulk_settings)

    # Create bulk relaxation task
    bulk_vasp = VaspTask(
        structure=bulk_structure,
        code=code,
        parameters=orm.Dict(dict=bulk_builder_inputs['parameters']),
        kpoints_spacing=orm.Float(bulk_builder_inputs.get('kpoints_spacing', 0.5)),
        potential_family=orm.Str(bulk_builder_inputs['potential_family']),
        potential_mapping=orm.Dict(dict=bulk_builder_inputs.get('potential_mapping', {})),
        options=orm.Dict(dict=bulk_builder_inputs['options']),
        settings=bulk_settings,
        clean_workdir=orm.Bool(bulk_builder_inputs.get('clean_workdir', False)),
    )

    # Extract energy from bulk relaxation
    from .utils import extract_total_energy_from_misc
    bulk_energy = extract_total_energy_from_misc(misc=bulk_vasp.misc)

    # =========================================================================
    # Task 0.5: Pristine Slab Relaxation (NEW)
    # =========================================================================

    # Prepare settings for pristine slab
    pristine_settings = builder_inputs.get('settings', {})
    if not pristine_settings:
        pristine_settings = {
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }
    if not isinstance(pristine_settings, orm.Dict):
        pristine_settings = orm.Dict(dict=pristine_settings)

    # Create pristine slab relaxation task
    pristine_vasp = VaspTask(
        structure=structure,  # Use input slab structure
        code=code,
        parameters=orm.Dict(dict=builder_inputs['parameters']),
        kpoints_spacing=orm.Float(builder_inputs.get('kpoints_spacing', 0.5)),
        potential_family=orm.Str(builder_inputs['potential_family']),
        potential_mapping=orm.Dict(dict=builder_inputs.get('potential_mapping', {})),
        options=orm.Dict(dict=builder_inputs['options']),
        settings=pristine_settings,
        clean_workdir=orm.Bool(builder_inputs.get('clean_workdir', False)),
    )

    pristine_energy = extract_total_energy_from_misc(misc=pristine_vasp.misc)

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
        builder_inputs=builder_inputs,
        max_parallel=max_parallel_jobs,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        fix_elements=fix_elements,
        structure_specific_builder_inputs=structure_specific_builder_inputs,
    )

    # Return raw outputs
    # Note: Result organization (successful/failed/statistics) must be done
    # in post-processing due to WorkGraph namespace serialization limitations.
    # Use the organize_hydroxylation_results() helper function after workflow completion.
    return {
        'manifest': gen_outputs.manifest,
        'structures': relax_outputs.structures,
        'energies': relax_outputs.energies,
        'bulk_structure': bulk_vasp.structure,
        'bulk_energy': bulk_energy,
        'pristine_structure': pristine_vasp.structure,
        'pristine_energy': pristine_energy,
    }


def build_surface_hydroxylation_workgraph(
    structure: orm.StructureData = None,
    structure_pk: int = None,
    surface_params: dict = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    builder_inputs: dict = None,
    bulk_structure: orm.StructureData = None,
    bulk_structure_pk: int = None,
    bulk_cif_path: str = None,
    bulk_builder_inputs: dict = None,
    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    structure_specific_builder_inputs: dict = None,
    name: str = 'SurfaceHydroxylation',
    # NEW: Surface energy calculation control
    calculate_surface_energies: bool = True,
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
        builder_inputs: Complete builder configuration (dict) for vasp.v2.vasp WorkChain.
            Provide ALL builder parameters you want to control:

            Required keys:
            - parameters: Dict with nested 'incar' dict
                {'incar': {'PREC': 'Accurate', 'ENCUT': 500, 'EDIFF': 1e-6, ...}}
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - options: Dict with scheduler settings
                {'resources': {'num_machines': 1, 'num_cores_per_machine': 16},
                 'queue_name': 'normal', 'max_wallclock_seconds': 3600}
        bulk_structure: Bulk crystal structure (StructureData).
                       Provide one of: bulk_structure, bulk_structure_pk, or bulk_cif_path.
        bulk_structure_pk: PK of bulk StructureData node (int).
        bulk_cif_path: Path to CIF file for bulk structure (str).
                      Example: 'ag3po4.cif'
        bulk_builder_inputs: VASP parameters for bulk relaxation (dict).
                            Must include ISIF=3 for cell relaxation.
                            If None, uses sensible defaults.

            Optional keys:
            - kpoints_spacing: Float (Å⁻¹, default: 0.5)
            - potential_mapping: Dict (element->potential, default: {})
            - settings: Dict (parser settings, default: add_trajectory/structure/kpoints)
            - dynamics: Dict (selective dynamics - auto-added if fix_type set)
            - clean_workdir: Bool (default: False)
            - max_iterations: Int (max restart attempts, default: 5)
            - verbose: Bool (detailed logging, default: False)

            Example:
                builder_inputs = {
                    'parameters': {'incar': {'ENCUT': 500, 'EDIFF': 1e-6, ...}},
                    'kpoints_spacing': 0.3,
                    'potential_family': 'PBE',
                    'potential_mapping': {},
                    'options': {'resources': {...}, 'queue_name': 'normal'},
                    'clean_workdir': False,
                }
        max_parallel_jobs: Number of structures to process in this run (default 2).
                          Uses simple batch approach: processes first N structures only.
                          Increase this value in subsequent runs to process more structures.
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None). Default: None (no fixing)
                 Example: 'bottom' to fix bottom 5Å of atoms in slab
        fix_thickness: Thickness in Angstroms for fixing region. Default: 0.0
                      Example: 5.0 to fix 5Å region
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O']).
                     If None, all elements in the region are fixed. Default: None
        structure_specific_builder_inputs: Optional dict to override builder_inputs for specific
                     structure indices. Keys are integer indices (0, 1, 2, ...) matching the
                     generated structure order. Values are PARTIAL builder_inputs dicts that will
                     be MERGED into the default builder. Only specify parameters you want to override.
                     Example: {0: {'parameters': {'incar': {'ALGO': 'Normal'}}, 'kpoints_spacing': 0.2}}
                              overrides only ALGO and kpoints_spacing for structure 0, keeping all
                              other parameters (potential_mapping, options, etc.) from default builder.
                     Structures not listed will use the complete default builder_inputs.
                     Useful for rerunning only failed calculations with different parameters.
                     Default: None (use default builder_inputs for all structures)
        name: Name for the workflow (default 'SurfaceHydroxylation')
        calculate_surface_energies: Enable automatic surface energy calculations at workflow end.
                     If True (default), calculates γ for all 3 formation reactions (H2O/O2, H2/H2O, H2/O2)
                     and exposes outputs: surface_energies_reaction1/2/3.
                     Set to False to disable surface energy calculations (v2 backward compat mode).
                     Default: True

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

    # ========================================================================
    # BULK STRUCTURE VALIDATION
    # ========================================================================

    # Validate exactly one bulk input is provided
    bulk_inputs_provided = sum([
        bulk_structure is not None,
        bulk_structure_pk is not None,
        bulk_cif_path is not None,
    ])

    if bulk_inputs_provided == 0:
        raise ValueError(
            "Bulk structure is required for surface energy calculations.\n"
            "Provide one of: bulk_structure, bulk_structure_pk, or bulk_cif_path.\n\n"
            "Examples:\n"
            "  bulk_cif_path='ag3po4.cif'\n"
            "  bulk_structure_pk=1234\n"
            "  bulk_structure=orm.load_node(1234)"
        )

    if bulk_inputs_provided > 1:
        raise ValueError(
            "Provide exactly ONE bulk structure input.\n"
            f"Found: bulk_structure={bulk_structure is not None}, "
            f"bulk_structure_pk={bulk_structure_pk is not None}, "
            f"bulk_cif_path={bulk_cif_path is not None}"
        )

    # Load/convert bulk structure to StructureData
    if bulk_cif_path is not None:
        try:
            from ase.io import read
            from pathlib import Path
            cif_path = Path(bulk_cif_path).resolve()  # Resolve symlinks and .. paths
            working_dir = Path.cwd().resolve()

            # Validate the file is within allowed directory
            if not str(cif_path).startswith(str(working_dir)):
                raise ValueError(
                    f"CIF file must be in current working directory.\n"
                    f"Requested: {cif_path}\n"
                    f"Working dir: {working_dir}\n"
                    f"Use relative paths only (e.g., 'structures/ag3po4.cif')"
                )

            if not cif_path.exists():
                raise FileNotFoundError(f"CIF file not found: {bulk_cif_path}")
            atoms = read(str(cif_path))
            bulk_structure = orm.StructureData(ase=atoms)
            print(f"   ✓ Loaded bulk structure from CIF: {bulk_cif_path}")
            print(f"   Formula: {atoms.get_chemical_formula()}")
        except Exception as e:
            raise ValueError(
                f"Could not read CIF file: {bulk_cif_path}\n"
                f"Error: {e}\n"
                f"Check file exists and is valid CIF format."
            )

    elif bulk_structure_pk is not None:
        try:
            bulk_structure = orm.load_node(bulk_structure_pk)
            if not isinstance(bulk_structure, orm.StructureData):
                raise ValueError(
                    f"Node {bulk_structure_pk} is not a StructureData.\n"
                    f"Found: {type(bulk_structure).__name__}"
                )
            print(f"   ✓ Loaded bulk structure from PK: {bulk_structure_pk}")
        except Exception as e:
            raise ValueError(
                f"Could not load bulk structure from PK {bulk_structure_pk}.\n"
                f"Error: {e}\n"
                f"Use 'verdi data core.structure list' to see available structures."
            )

    # If bulk_structure provided directly, validate type
    elif bulk_structure is not None:
        if not isinstance(bulk_structure, orm.StructureData):
            raise ValueError(
                f"bulk_structure must be StructureData, got {type(bulk_structure).__name__}"
            )

    # Set default bulk_builder_inputs if not provided
    if bulk_builder_inputs is None:
        print("   ⚠ Using default bulk_builder_inputs (ISIF=3, ENCUT=500)")
        # Default bulk relaxation parameters (based on default_ag3po4_builders.py)
        bulk_builder_inputs = {
            'parameters': {
                'incar': {
                    'PREC': 'Accurate',
                    'ENCUT': 500,
                    'EDIFF': 1e-6,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ALGO': 'Fast',
                    'LREAL': False,
                    'NELM': 200,
                    'LWAVE': False,
                    'LCHARG': False,
                    'ISIF': 3,        # Cell + ionic relaxation for bulk
                    'NSW': 500,
                    'IBRION': 2,
                    'EDIFFG': -0.01,  # Tighter convergence for bulk
                }
            },
            'kpoints_spacing': 0.3,  # Denser k-points for bulk
            'potential_family': 'PBE',
            'potential_mapping': {},
            'options': {
                'resources': {
                    'num_machines': 1,
                    'num_mpiprocs_per_machine': 16,
                },
                'queue_name': 'normal',
                'max_wallclock_seconds': 3600 * 4,
            },
            'clean_workdir': False,
        }
    # Validate ISIF=3 for user-provided bulk_builder_inputs (CRITICAL)
    else:
        # User provided bulk_builder_inputs - validate ISIF=3
        isif_value = bulk_builder_inputs.get('parameters', {}).get('incar', {}).get('ISIF')
        if isif_value != 3:
            raise ValueError(
                f"bulk_builder_inputs MUST include ISIF=3 for cell relaxation.\n"
                f"Found: ISIF={isif_value}\n\n"
                f"Bulk reference calculations require full cell relaxation (ISIF=3).\n"
                f"This is critical for accurate surface energy calculations per Section S2.\n\n"
                f"Fix: bulk_builder_inputs['parameters']['incar']['ISIF'] = 3"
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

    if builder_inputs is None:
        # Default lightweight VASP builder configuration
        builder_inputs = {
            'parameters': {
                'incar': {
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
                }
            },
            'kpoints_spacing': 0.5,
            'potential_family': 'PBE',
            'potential_mapping': {},
            'options': {
                'resources': {
                    'num_machines': 1,
                    'num_mpiprocs_per_machine': 16,
                },
                'queue_name': 'normal',
                'max_wallclock_seconds': 3600,
            },
            'clean_workdir': False,
        }

    # ========================================================================
    # BUILD WORKGRAPH
    # ========================================================================

    # Convert structure_specific_builder_inputs integer keys to strings
    # (AiiDA Dict requires string keys for serialization)
    if structure_specific_builder_inputs is not None:
        structure_specific_builder_inputs_str = {}
        for key, value in structure_specific_builder_inputs.items():
            structure_specific_builder_inputs_str[str(key)] = value
        structure_specific_builder_inputs = structure_specific_builder_inputs_str

    # Build the WorkGraph using the @task.graph function
    wg = SurfaceHydroxylationWorkGraph.build(
        structure=structure,
        surface_params=surface_params,
        code=code,
        builder_inputs=builder_inputs,
        bulk_structure=bulk_structure,
        bulk_builder_inputs=bulk_builder_inputs,
        max_parallel_jobs=max_parallel_jobs,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        fix_elements=fix_elements,
        structure_specific_builder_inputs=structure_specific_builder_inputs,
    )

    # Set the workflow name
    wg.name = name

    # Surface energy calculations (if enabled)
    if calculate_surface_energies:
        # Add surface energy calculation task
        # Connect directly to the WorkGraph outputs
        se_task = wg.add_task(
            calculate_all_surface_energies,
            name='surface_energy_calc',
            structures_dict=wg.outputs['structures'].value,
            energies_dict=wg.outputs['energies'].value,
            bulk_structure=wg.outputs['bulk_structure'].value,
            bulk_energy=wg.outputs['bulk_energy'].value,
            temperature=298.0,
            pressures={'H2O': 0.023, 'O2': 0.21, 'H2': 1.0},
        )

        # Expose outputs
        wg.outputs.new(
            'surface_energies_reaction1',
            value=se_task.outputs.reaction1_results
        )
        wg.outputs.new(
            'surface_energies_reaction2',
            value=se_task.outputs.reaction2_results
        )
        wg.outputs.new(
            'surface_energies_reaction3',
            value=se_task.outputs.reaction3_results
        )

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
            - reference_data: dict with bulk and pristine reference calculations

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

    # Extract reference data from bulk and pristine calculations
    reference_data = {
        'bulk_structure_pk': workflow_node.outputs.bulk_structure.pk,
        'bulk_energy': workflow_node.outputs.bulk_energy.value,
        'pristine_structure_pk': workflow_node.outputs.pristine_structure.pk,
        'pristine_energy': workflow_node.outputs.pristine_energy.value,
    }

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
        },
        'reference_data': reference_data,
    }
