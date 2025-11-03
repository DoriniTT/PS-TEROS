"""
Slab Generation Module

This module provides functions to generate surface slabs from bulk structures
using Pymatgen's SlabGenerator with scatter-gather pattern for parallel relaxation.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def extract_restart_folders_from_node(node_pk: int) -> dict[str, orm.RemoteData]:
    """
    Extract restart folders from a previous PS-TEROS workgraph run.
    
    This function loads a completed (or failed) PS-TEROS workgraph node and extracts
    the RemoteData nodes from its slab_remote outputs. These can then be used to
    restart the slab relaxation calculations from where they left off.
    
    Args:
        node_pk: PK of the previous PS-TEROS workgraph node
    
    Returns:
        Dictionary mapping slab labels (e.g., 'term_0', 'term_1') to RemoteData nodes
    
    Raises:
        ValueError: If the node doesn't have slab_remote outputs or is invalid
    
    Example:
        >>> restart_folders = extract_restart_folders_from_node(19774)
        >>> print(restart_folders.keys())
        dict_keys(['term_0', 'term_1'])
        >>> # Use in a new workgraph
        >>> wg = build_core_workgraph(..., restart_from_node=19774)
    """
    try:
        node = orm.load_node(node_pk)
    except Exception as e:
        raise ValueError(f"Could not load node {node_pk}: {e}")
    
    # Check if the node has slab_remote outputs
    if not hasattr(node.outputs, 'slab_remote'):
        raise ValueError(
            f"Node {node_pk} does not have 'slab_remote' outputs. "
            f"This might not be a PS-TEROS workgraph with slab relaxations, "
            f"or it was run with an older version of the code. "
            f"Available outputs: {list(node.outputs)}"
        )
    
    # Extract the RemoteData nodes
    slab_remote = node.outputs.slab_remote
    restart_folders = {}
    
    for label in slab_remote.keys():
        remote_data = slab_remote[label]
        if not isinstance(remote_data, orm.RemoteData):
            raise ValueError(
                f"Expected RemoteData for {label}, got {type(remote_data).__name__}"
            )
        restart_folders[label] = remote_data
    
    return restart_folders


@task.graph
def wrap_input_slabs(
    **slabs: orm.StructureData
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """
    Wrap user-provided slab structures into the proper namespace format.
    
    This function takes pre-generated slab structures and returns them
    in the same format as generate_slab_structures, allowing user-provided
    slabs to work with the scatter-gather relaxation pattern.
    
    This is a @task.graph (not calcfunction) to avoid creating cycles in the 
    provenance graph when using already-stored StructureData nodes. This matches
    the scatter-gather pattern used in relax_slabs_scatter.
    
    Args:
        **slabs: Keyword arguments where each key is a slab identifier
                 (e.g., "term_0") and value is a StructureData node
    
    Returns:
        Dictionary with key 'slabs' containing a dict of slab structures.
        Each slab is keyed by termination identifier (e.g., "term_0", "term_1")
        and contains AiiDA StructureData nodes.
    """
    # Simply return the slabs in the same format as generate_slab_structures
    # Using @task.graph instead of @task.calcfunction avoids provenance cycles
    return {'slabs': dict(slabs)}


@task.graph
def collect_slab_outputs(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    remote_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)],
) -> t.Annotated[dict, namespace(
    structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
    remote_folders=dynamic(orm.RemoteData)
)]:
    """
    Collect and passthrough slab outputs in the proper format.
    
    This is a passthrough task that takes dictionaries of slab outputs
    and returns them in the format expected by WorkGraph dynamic outputs.
    Used in restart mode to properly expose outputs from individual VASP tasks.
    
    Args:
        structures: Dictionary of relaxed slab structures
        energies: Dictionary of slab energies
        remote_folders: Dictionary of RemoteData nodes
    
    Returns:
        Dictionary with structures, energies, and remote_folders namespaces
    """
    return {
        'structures': structures,
        'energies': energies,
        'remote_folders': remote_folders,
    }


@task.calcfunction
def generate_slab_structures(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    min_slab_thickness: orm.Float,
    min_vacuum_thickness: orm.Float,
    lll_reduce: orm.Bool,
    center_slab: orm.Bool,
    symmetrize: orm.Bool,
    primitive: orm.Bool,
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """
    Generate slab structures from a bulk crystal structure using Pymatgen's SlabGenerator.

    This function wraps Pymatgen's slab generation capabilities to produce various
    surface terminations for a given bulk material and Miller index. All slabs are
    generated with symmetric, reduced, c-axis orthogonal cells.

    Args:
        bulk_structure: AiiDA StructureData of the bulk crystal
        miller_indices: AiiDA List of Miller indices for slab generation (e.g., [1, 0, 0])
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        lll_reduce: Reduce cell using LLL algorithm before slab generation
        center_slab: Center the slab in the c direction of the cell
        symmetrize: Generate symmetrically distinct terminations
        primitive: Find primitive cell before generating slabs

    Returns:
        Dictionary with key 'slabs' containing a dict of slab structures.
        Each slab is keyed by termination identifier (e.g., "term_0", "term_1")
        and contains AiiDA StructureData nodes.
    """
    adaptor = AseAtomsAdaptor()

    ase_structure = bulk_structure.get_ase()
    pymatgen_structure = adaptor.get_structure(ase_structure)

    if primitive.value:
        analyzer = SpacegroupAnalyzer(pymatgen_structure)
        pymatgen_structure = analyzer.get_primitive_standard_structure()

    generator = SlabGenerator(
        pymatgen_structure,
        tuple(miller_indices.get_list()),
        min_slab_thickness.value,
        min_vacuum_thickness.value,
        center_slab=center_slab.value,
        in_unit_planes=False,
        max_normal_search=None,
        lll_reduce=lll_reduce.value,
    )

    slab_nodes: dict[str, orm.StructureData] = {}
    for index, slab in enumerate(generator.get_slabs(symmetrize=symmetrize.value)):
        orthogonal_slab = slab.get_orthogonal_c_slab()
        ase_slab = adaptor.get_atoms(orthogonal_slab)
        slab_nodes[f"term_{index}"] = orm.StructureData(ase=ase_slab)

    if not slab_nodes:
        raise ValueError(
            f"No slabs were generated for the specified surface orientation "
            f"({miller_index.value[0]}, {miller_index.value[1]}, {miller_index.value[2]}). "
            f"This may occur if the Miller indices are invalid for the given structure, "
            f"or if the slab generation parameters (min_slab_size={min_slab_size.value}, "
            f"min_vacuum_size={min_vacuum_size.value}) are too restrictive."
        )

    return {'slabs': slab_nodes}


@task.calcfunction
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP energies output.
    
    Args:
        energies: Dictionary containing energy outputs from VASP (from misc output)
    
    Returns:
        Total energy as Float
    """
    energy_dict = energies.get_dict()
    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in VASP outputs. Available keys: {available}')


def get_settings():
    """
    Parser settings for aiida-vasp.
    """
    return {
        'parser_settings': {
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }


@task.graph
def scf_relax_and_calculate_relaxation_energy(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    restart_folders: dict = None,
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    unrelaxed_energies=dynamic(orm.Float),
    unrelaxed_remote_folders=dynamic(orm.RemoteData),
    relaxed_structures=dynamic(orm.StructureData),
    relaxed_energies=dynamic(orm.Float),
    relaxed_remote_folders=dynamic(orm.RemoteData),
    relaxation_energies=dynamic(orm.Float),
)]:
    """
    Combined workflow: SCF → Relax → Calculate relaxation energy for all slabs.
    
    This function performs three steps in sequence:
    1. SCF calculation on unrelaxed slabs (NSW=0, IBRION=-1)
    2. Relaxation of slabs (user-specified parameters)
    3. Calculation of relaxation energy (E_relaxed - E_unrelaxed)
    
    Args:
        slabs: Dictionary of slab structures
        code: AiiDA code for VASP
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        parameters: VASP parameters (for relaxation; NSW/IBRION overridden for SCF)
        options: Computation resources
        kpoints_spacing: K-points spacing (optional)
        clean_workdir: Whether to clean remote directories
        restart_folders: Dictionary of RemoteData for restarting relaxations (optional)
        max_number_jobs: Maximum number of concurrent VASP calculations (optional)

    Returns:
        Dictionary with all outputs:
        - unrelaxed_energies: Energies from SCF
        - unrelaxed_remote_folders: RemoteData from SCF
        - relaxed_structures: Relaxed structures
        - relaxed_energies: Energies from relaxation
        - relaxed_remote_folders: RemoteData from relaxation
        - relaxation_energies: E_relaxed - E_unrelaxed
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    vasp_task_cls = task(vasp_wc)
    
    # Dictionaries for outputs
    unrelaxed_energies_dict: dict[str, orm.Float] = {}
    unrelaxed_remote_dict: dict[str, orm.RemoteData] = {}
    relaxed_structures_dict: dict[str, orm.StructureData] = {}
    relaxed_energies_dict: dict[str, orm.Float] = {}
    relaxed_remote_dict: dict[str, orm.RemoteData] = {}
    relaxation_energies_dict: dict[str, orm.Float] = {}
    
    # Process each slab
    for label, structure in slabs.items():
        # ==== STEP 1: SCF on unrelaxed slab ====
        scf_params = dict(parameters)
        scf_params['NSW'] = 0
        scf_params['IBRION'] = -1
        
        scf_inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': scf_params},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            scf_inputs['kpoints_spacing'] = kpoints_spacing
        
        scf_calc = vasp_task_cls(**scf_inputs)
        unrelaxed_energies_dict[label] = extract_total_energy(energies=scf_calc.misc).result
        unrelaxed_remote_dict[label] = scf_calc.remote_folder
        
        # ==== STEP 2: Relaxation ====
        relax_inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            relax_inputs['kpoints_spacing'] = kpoints_spacing
        if restart_folders is not None and label in restart_folders:
            relax_inputs['restart_folder'] = restart_folders[label]
        
        relax_calc = vasp_task_cls(**relax_inputs)
        relaxed_structures_dict[label] = relax_calc.structure
        relaxed_energies_dict[label] = extract_total_energy(energies=relax_calc.misc).result
        relaxed_remote_dict[label] = relax_calc.remote_folder
        
        # ==== STEP 3: Calculate relaxation energy ====
        relax_energy = calculate_energy_difference(
            energy_final=relaxed_energies_dict[label],
            energy_initial=unrelaxed_energies_dict[label],
        ).result
        relaxation_energies_dict[label] = relax_energy
    
    # Return all outputs
    return {
        'unrelaxed_energies': unrelaxed_energies_dict,
        'unrelaxed_remote_folders': unrelaxed_remote_dict,
        'relaxed_structures': relaxed_structures_dict,
        'relaxed_energies': relaxed_energies_dict,
        'relaxed_remote_folders': relaxed_remote_dict,
        'relaxation_energies': relaxation_energies_dict,
    }


@task.graph
def scf_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    max_number_jobs: int = None,
    # NEW: Builder inputs (full control over VASP builder)
    scf_builder_inputs: dict | None = None,
) -> t.Annotated[dict, namespace(energies=dynamic(orm.Float), remote_folders=dynamic(orm.RemoteData))]:
    """
    Scatter-gather phase: perform SCF calculations on unrelaxed slab structures in parallel.

    This task graph performs single-point energy calculations (NSW=0, IBRION=-1)
    on unrelaxed slab structures to obtain the initial energy before relaxation.
    This is needed for calculating relaxation energies.

    Args:
        slabs: Dictionary of slab structures (unrelaxed)
        code: AiiDA code for VASP
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        parameters: VASP parameters (old-style, NSW and IBRION will be overridden)
        options: Computation resources (old-style)
        kpoints_spacing: K-points spacing (optional, old-style)
        clean_workdir: Whether to clean remote directories
        max_number_jobs: Maximum number of concurrent VASP calculations (optional)
        scf_builder_inputs: Dict of builder parameters for full control over VASP.
                           Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                           NSW and IBRION will be overridden to 0 and -1 respectively.
                           If provided, overrides parameters, options, kpoints_spacing.
                           Default: None

    Returns:
        Dictionary with energies and remote_folders namespaces
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    scf_task_cls = task(vasp_wc)
    
    energies_ns: dict[str, orm.Float] = {}
    remote_folders: dict[str, orm.RemoteData] = {}

    # NEW: Use builder_inputs if provided (full control), otherwise use old-style params
    if scf_builder_inputs is not None:
        # New-style: extract from builder_inputs
        actual_params = scf_builder_inputs.get('parameters', {'incar': dict(parameters)})
        if 'incar' not in actual_params:
            # Assume the dict IS the INCAR
            actual_params = {'incar': actual_params}
        actual_options = scf_builder_inputs.get('options', dict(options))
        actual_kpoints = scf_builder_inputs.get('kpoints_spacing', kpoints_spacing)
    else:
        # Old-style: use provided params
        actual_params = {'incar': dict(parameters)}
        actual_options = dict(options)
        actual_kpoints = kpoints_spacing

    # Scatter: create independent SCF tasks for each slab (runs in parallel)
    for label, structure in slabs.items():
        # Create SCF parameters by overriding NSW and IBRION
        scf_params = dict(actual_params['incar'])
        scf_params['NSW'] = 0
        scf_params['IBRION'] = -1

        inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': scf_params},
            'options': actual_options,
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if actual_kpoints is not None:
            inputs['kpoints_spacing'] = actual_kpoints
        
        scf_calc = scf_task_cls(**inputs)
        energies_ns[label] = extract_total_energy(energies=scf_calc.misc).result
        remote_folders[label] = scf_calc.remote_folder

    # Gather: return collected results
    return {
        'energies': energies_ns,
        'remote_folders': remote_folders,
    }


@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    restart_folders: dict = None,
    max_number_jobs: int = None,
    # NEW: Builder inputs (full control over VASP builder)
    relax_builder_inputs: dict | None = None,
) -> t.Annotated[dict, namespace(relaxed_structures=dynamic(orm.StructureData), energies=dynamic(orm.Float), remote_folders=dynamic(orm.RemoteData))]:
    """
    Scatter-gather phase: relax each slab structure in parallel.

    This task graph iterates over the input slabs dictionary and creates
    independent relaxation tasks that run in parallel. The loop wrapping
    is necessary because slabs is a future output that isn't available
    at graph construction time.

    Args:
        slabs: Dictionary of slab structures to relax
        code: AiiDA code for VASP
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        parameters: VASP parameters (old-style, for backward compat)
        options: Computation resources (old-style, for backward compat)
        kpoints_spacing: K-points spacing (optional, old-style)
        clean_workdir: Whether to clean remote directories
        restart_folders: Dictionary of RemoteData nodes for restarting calculations (optional).
                        Keys must match slab labels (e.g., 'term_0', 'term_1').
                        If provided, calculations will restart from these remote folders.
        max_number_jobs: Maximum number of concurrent VASP calculations (optional)
        relax_builder_inputs: Dict of builder parameters for full control over VASP.
                             Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                             If provided, overrides parameters, options, kpoints_spacing.
                             Default: None

    Returns:
        Dictionary with relaxed_structures, energies, and remote_folders namespaces
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_wc)
    
    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}
    remote_folders: dict[str, orm.RemoteData] = {}

    # NEW: Use builder_inputs if provided (full control), otherwise use old-style params
    if relax_builder_inputs is not None:
        # New-style: extract from builder_inputs
        actual_params = relax_builder_inputs.get('parameters', {'incar': dict(parameters)})
        if 'incar' not in actual_params:
            # Assume the dict IS the INCAR
            actual_params = {'incar': actual_params}
        actual_options = relax_builder_inputs.get('options', dict(options))
        actual_kpoints = relax_builder_inputs.get('kpoints_spacing', kpoints_spacing)
    else:
        # Old-style: use provided params
        actual_params = {'incar': dict(parameters)}
        actual_options = dict(options)
        actual_kpoints = kpoints_spacing

    # Scatter: create independent tasks for each slab (runs in parallel)
    for label, structure in slabs.items():
        inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': actual_params,
            'options': actual_options,
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if actual_kpoints is not None:
            inputs['kpoints_spacing'] = actual_kpoints
        
        # Add restart_folder if provided for this slab
        if restart_folders is not None and label in restart_folders:
            inputs['restart_folder'] = restart_folders[label]
        
        relaxation = relax_task_cls(**inputs)
        relaxed[label] = relaxation.structure
        energies_ns[label] = extract_total_energy(energies=relaxation.misc).result
        remote_folders[label] = relaxation.remote_folder

    # Gather: return collected results
    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
        'remote_folders': remote_folders,
    }


@task.graph
def calculate_relaxation_energies_scatter(
    unrelaxed_energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    relaxed_energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
) -> t.Annotated[dict, namespace(relaxation_energies=dynamic(orm.Float))]:
    """
    Calculate relaxation energy for each slab termination.
    
    The relaxation energy is defined as:
        E_relax = E_relaxed - E_unrelaxed
    
    A negative relaxation energy indicates that the relaxation stabilizes the system.
    
    Args:
        unrelaxed_energies: Dictionary of energies from SCF calculations on unrelaxed slabs
        relaxed_energies: Dictionary of energies from relaxation calculations
    
    Returns:
        Dictionary with relaxation_energies namespace containing the energy difference
        for each slab termination
    """
    relaxation_energies: dict[str, orm.Float] = {}
    
    for label in unrelaxed_energies.keys():
        if label in relaxed_energies:
            relax_energy = calculate_energy_difference(
                energy_final=relaxed_energies[label],
                energy_initial=unrelaxed_energies[label],
            ).result
            relaxation_energies[label] = relax_energy
    
    return {'relaxation_energies': relaxation_energies}


@task.calcfunction
def calculate_energy_difference(
    energy_final: orm.Float,
    energy_initial: orm.Float,
) -> orm.Float:
    """
    Calculate energy difference: E_final - E_initial
    
    Args:
        energy_final: Final energy (e.g., relaxed)
        energy_initial: Initial energy (e.g., unrelaxed)
    
    Returns:
        Energy difference as Float (in eV)
    """
    return orm.Float(energy_final.value - energy_initial.value)


@task.graph
def calculate_electronic_properties_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    slab_electronic_properties: t.Mapping[str, t.Mapping[str, t.Any]],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    clean_workdir: bool,
    default_bands_parameters: t.Mapping[str, t.Any] = None,
    default_bands_options: t.Mapping[str, t.Any] = None,
    default_band_settings: t.Mapping[str, t.Any] = None,
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    slab_bands=dynamic(orm.BandsData),
    slab_dos=dynamic(orm.ArrayData),
    slab_primitive_structures=dynamic(orm.StructureData),
    slab_seekpath_parameters=dynamic(orm.Dict),
)]:
    """
    Scatter-gather phase: calculate electronic properties (DOS and bands) for selected slabs.

    This task graph creates independent BandsWorkChain tasks for each selected slab,
    allowing per-slab parameter overrides. Only slabs present in slab_electronic_properties
    dictionary will have their electronic properties calculated.

    Args:
        slabs: Dictionary of relaxed slab structures
        slab_electronic_properties: Dictionary mapping slab labels to parameter configs.
                                    Each config can contain:
                                    - 'bands_parameters': Dict with 'scf', 'bands', 'dos' keys
                                    - 'bands_options': Scheduler options
                                    - 'band_settings': Band workflow settings
                                    If keys are missing, defaults are used.
        code: AiiDA code for VASP
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        clean_workdir: Whether to clean remote directories
        default_bands_parameters: Default parameters (fallback)
        default_bands_options: Default scheduler options (fallback)
        default_band_settings: Default band settings (fallback)
        max_number_jobs: Maximum number of concurrent VASP calculations (None = unlimited)

    Returns:
        Dictionary with four namespaces:
        - slab_bands: BandsData for each selected slab
        - slab_dos: DOS ArrayData for each selected slab
        - slab_primitive_structures: Primitive structures used for band paths
        - slab_seekpath_parameters: Seekpath parameters for each slab

    Example:
        >>> scatter_outputs = calculate_electronic_properties_slabs_scatter(
        ...     slabs=relaxed_slabs_dict,
        ...     slab_electronic_properties={
        ...         'term_0': {'bands_parameters': defaults},
        ...         'term_2': {'bands_parameters': custom_params, 'bands_options': high_mem},
        ...     },
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     clean_workdir=True,
        ... )
    """
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # Get BandsWorkChain and wrap as task
    BandsWorkChain = WorkflowFactory('vasp.v2.bands')
    BandsTask = task(BandsWorkChain)

    # Output dictionaries
    bands_dict: dict[str, orm.BandsData] = {}
    dos_dict: dict[str, orm.ArrayData] = {}
    primitive_structures_dict: dict[str, orm.StructureData] = {}
    seekpath_parameters_dict: dict[str, orm.Dict] = {}

    # Scatter: create independent tasks for each selected slab
    for label, slab_config in slab_electronic_properties.items():
        # Skip if slab doesn't exist
        if label not in slabs:
            continue

        structure = slabs[label]

        # Extract per-slab parameters with fallback to defaults
        bands_params = slab_config.get('bands_parameters', default_bands_parameters or {})
        bands_opts = slab_config.get('bands_options', default_bands_options or {})
        band_settings = slab_config.get('band_settings', default_band_settings or {})

        # Build BandsWorkChain inputs (same structure as bulk)
        bands_inputs: dict[str, t.Any] = {
            'structure': structure,
            'metadata': {
                'label': f'Slab_{label}_Electronic_Properties',
                'description': f'DOS and bands calculation for slab termination {label}',
            }
        }

        # Add band_settings if provided
        if band_settings:
            bands_inputs['band_settings'] = band_settings

        # Build SCF namespace inputs (required)
        scf_inputs = {
            'code': code,
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'options': dict(bands_opts),
            'clean_workdir': clean_workdir,
        }
        if bands_params and 'scf' in bands_params:
            scf_inputs['parameters'] = {'incar': bands_params['scf']}
        if bands_params and 'scf_kpoints_distance' in bands_params:
            scf_inputs['kpoints_spacing'] = bands_params['scf_kpoints_distance']

        bands_inputs['scf'] = scf_inputs

        # Build Bands namespace inputs (optional)
        if bands_params and 'bands' in bands_params:
            bands_inputs['bands'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'options': dict(bands_opts),
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_params['bands']},
            }

        # Build DOS namespace inputs (optional)
        if bands_params and 'dos' in bands_params:
            bands_inputs['dos'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'options': dict(bands_opts),
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_params['dos']},
            }

        # Create BandsWorkChain task
        bands_task = BandsTask(**bands_inputs)

        # Collect outputs
        bands_dict[label] = bands_task.band_structure
        dos_dict[label] = bands_task.dos
        primitive_structures_dict[label] = bands_task.primitive_structure
        seekpath_parameters_dict[label] = bands_task.seekpath_parameters

    # Gather: return collected results
    return {
        'slab_bands': bands_dict,
        'slab_dos': dos_dict,
        'slab_primitive_structures': primitive_structures_dict,
        'slab_seekpath_parameters': seekpath_parameters_dict,
    }
