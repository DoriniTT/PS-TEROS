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
        parameters: VASP parameters
        options: Computation resources
        kpoints_spacing: K-points spacing (optional)
        clean_workdir: Whether to clean remote directories
        restart_folders: Dictionary of RemoteData nodes for restarting calculations (optional).
                        Keys must match slab labels (e.g., 'term_0', 'term_1'). 
                        If provided, calculations will restart from these remote folders.
    
    Returns:
        Dictionary with relaxed_structures, energies, and remote_folders namespaces
    """
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_wc)
    
    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}
    remote_folders: dict[str, orm.RemoteData] = {}

    # Scatter: create independent tasks for each slab (runs in parallel)
    for label, structure in slabs.items():
        inputs: dict[str, t.Any] = {
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
            inputs['kpoints_spacing'] = kpoints_spacing
        
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
