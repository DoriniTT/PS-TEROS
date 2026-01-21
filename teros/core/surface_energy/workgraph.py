"""
Metal and Intermetallic Surface Energy WorkGraph Builder.

This module provides the main workflow builder for computing surface energies
of elemental metals and stoichiometric intermetallics. It supports:
- Multiple Miller indices in a single parent workflow
- Single bulk calculation shared across all orientations
- Multiple terminations per orientation
- Bulk relaxation → Slab generation → Slab relaxation → Surface energy calculation

Supported materials:
- Elemental metals: Au, Ag, Cu, Pt, Pd, Ni, Fe, etc.
- Stoichiometric intermetallics: PdIn, AuCu, NiAl, Cu3Au, etc.

For stoichiometric and symmetric intermetallic surfaces, the simple formula
γ = (E_slab - N·E_bulk/atom) / (2A) is used, same as for elemental metals.
"""

from __future__ import annotations

import typing as t
import os

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, dynamic, namespace
from ase.io import read

from teros.core.slabs import (
    generate_slab_structures,
    relax_slabs_scatter,
    extract_total_energy,
)
from teros.core.surface_energy.surface_energy import (
    calculate_metal_surface_energy,
    calculate_bulk_energy_per_atom,
)
from teros.core.surface_energy.wulff import build_wulff_shape
from teros.core.surface_energy.stoichiometric_finder import (
    generate_stoichiometric_symmetric_slabs,
    NoStoichiometricSymmetricSurfaceError,
)


def load_structure_from_file(filepath: str) -> orm.StructureData:
    """
    Load structure from a file using ASE.

    Args:
        filepath: Path to structure file

    Returns:
        StructureData node
    """
    atoms = read(filepath)
    return orm.StructureData(ase=atoms)


def get_settings():
    """Parser settings for aiida-vasp."""
    return {
        'parser_settings': {
            'add_energy': True,
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }


def miller_to_key(miller: list) -> str:
    """
    Convert Miller indices list to a key string.
    
    Args:
        miller: List like [1, 1, 1]
        
    Returns:
        String like 'hkl_111'
    """
    return f'hkl_{miller[0]}{miller[1]}{miller[2]}'


@task.calcfunction
def gather_surface_energies(**kwargs) -> orm.Dict:
    """
    Gather surface energies from multiple orientations into a single dictionary.
    
    Input kwargs are expected to be namespaces like surface_energies_hkl_111.
    Each input is a dynamic namespace of per-termination results.
    
    Returns:
        Dict with structure:
        {
            'hkl_111': {
                'term_0': {...},
                'term_1': {...}
            },
            'hkl_100': {
                'term_0': {...}
            }
        }
    """
    results = {}
    for key, val in kwargs.items():
        # key corresponds to the task output name (e.g., surface_energies_hkl_111)
        # val is the plain Python dict (since it's a dynamic namespace) or Dict node
        
        # Strip the prefix to get the hkl key
        hkl_key = key.replace('surface_energies_', '')
        
        # Process the value
        if isinstance(val, orm.Dict):
             results[hkl_key] = val.get_dict()
        elif isinstance(val, dict):
            # It's a namespace dict of nodes (or values if already loaded)
            results[hkl_key] = {}
            for term_key, term_node in val.items():
                if isinstance(term_node, orm.Dict):
                    results[hkl_key][term_key] = term_node.get_dict()
                else:
                    results[hkl_key][term_key] = term_node
        else:
            results[hkl_key] = val
            
    return orm.Dict(dict=results)


@task.graph
def compute_metal_surface_energies_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
) -> t.Annotated[dict, namespace(surface_energies=dynamic(orm.Dict))]:
    """
    Scatter-gather pattern for computing surface energies for multiple metal slabs.
    """
    surface_results = {}
    
    for key, slab_structure in slabs.items():
        slab_energy = energies[key]
        
        surface_data = calculate_metal_surface_energy(
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energy,
        ).result
        
        surface_results[key] = surface_data
    
    return {'surface_energies': surface_results}


@task.graph(outputs=[
    'slab_structures', 'relaxed_slabs', 'slab_energies', 'surface_energies', 'slab_remote',
])
def metal_surface_energy_single_hkl(
    # Bulk inputs (already relaxed!)
    bulk_structure: orm.StructureData = None,
    bulk_energy: orm.Float = None,

    # VASP configuration
    code: orm.Code = None,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    clean_workdir: bool = False,

    # Slab generation - single Miller index
    miller_indices: list = None,
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,

    # Slab relaxation parameters
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_kpoints_spacing: float = None,

    # Concurrency
    max_concurrent_jobs: int = None,

    # Sequential execution - creates data dependency for chaining tasks
    wait_for: orm.Data = None,

    # Stoichiometric+Symmetric Surface Finder (EXPERIMENTAL)
    require_stoichiometric_symmetric: bool = False,
    stoichiometric_strategies: list = None,
    stoichiometric_max_thickness: float = 30.0,
):
    """
    Slab generation and relaxation for a SINGLE Miller index.

    Receives already-relaxed bulk structure and energy from parent.

    If require_stoichiometric_symmetric=True, uses the experimental stoichiometric
    finder to generate only stoichiometric AND symmetric slabs.
    """
    # ===== SLAB GENERATION =====
    if require_stoichiometric_symmetric:
        # Use stoichiometric+symmetric finder (EXPERIMENTAL)
        slab_generation_outputs = generate_stoichiometric_symmetric_slabs(
            bulk_structure=bulk_structure,
            miller_indices=orm.List(list=[miller_indices]),  # Single Miller index in list
            min_slab_thickness=orm.Float(min_slab_thickness),
            min_vacuum_thickness=orm.Float(min_vacuum_thickness),
            lll_reduce=orm.Bool(lll_reduce),
            center_slab=orm.Bool(center_slab),
            primitive=orm.Bool(primitive),
            strategies=orm.List(list=stoichiometric_strategies) if stoichiometric_strategies else None,
            max_thickness=orm.Float(stoichiometric_max_thickness),
        )
        # Extract slabs for the single Miller index
        hkl_key = f"hkl_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}"
        slab_namespace = slab_generation_outputs['slabs'][hkl_key]
    else:
        # Standard slab generation
        slab_namespace = generate_slab_structures(
            bulk_structure=bulk_structure,
            miller_indices=orm.List(list=miller_indices),
            min_slab_thickness=orm.Float(min_slab_thickness),
            min_vacuum_thickness=orm.Float(min_vacuum_thickness),
            lll_reduce=orm.Bool(lll_reduce),
            center_slab=orm.Bool(center_slab),
            symmetrize=orm.Bool(symmetrize),
            primitive=orm.Bool(primitive),
        ).slabs
    
    # ===== SLAB RELAXATION =====
    relaxation_outputs = relax_slabs_scatter(
        slabs=slab_namespace,
        code=code,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=slab_parameters,
        options=slab_options,
        kpoints_spacing=slab_kpoints_spacing,
        clean_workdir=clean_workdir,
        max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
    )
    
    # ===== SURFACE ENERGY CALCULATION =====
    surface_outputs = compute_metal_surface_energies_scatter(
        slabs=relaxation_outputs.relaxed_structures,
        energies=relaxation_outputs.energies,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
    )
    
    return {
        'slab_structures': slab_namespace,
        'relaxed_slabs': relaxation_outputs.relaxed_structures,
        'slab_energies': relaxation_outputs.energies,
        'surface_energies': surface_outputs.surface_energies,
        'slab_remote': relaxation_outputs.remote_folders,
    }


def build_metal_surface_energy_workgraph(
    # Structure input (file path or StructureData)
    bulk_structure_path: str = None,
    bulk_structure: orm.StructureData = None,
    
    # VASP configuration
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    clean_workdir: bool = False,
    
    # Bulk parameters
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    
    # Slab generation - multiple Miller indices!
    miller_indices: list = None,  # List of Miller indices: [[1,1,1], [1,0,0], [1,1,0]]
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    
    # Slab relaxation
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_kpoints_spacing: float = None,
    
    # Concurrency
    max_concurrent_jobs: int = 4,

    # Stoichiometric+Symmetric Surface Finder (EXPERIMENTAL)
    require_stoichiometric_symmetric: bool = False,
    stoichiometric_strategies: list = None,
    stoichiometric_bonds: dict = None,
    stoichiometric_auto_detect_bonds: bool = False,
    stoichiometric_max_thickness: float = 30.0,

    # Workflow name
    name: str = 'MetalSurfaceEnergy',
) -> WorkGraph:
    """
    Build a single WorkGraph for metal surface energy calculations with multiple Miller indices.
    
    This creates ONE parent WorkGraph that:
    1. Runs bulk relaxation ONCE
    2. Spawns sub-tasks for each Miller index
    3. Consolidates all surface energies into ONE output Dict
    
    Args:
        bulk_structure_path: Path to bulk structure file (e.g., 'au.cif')
        bulk_structure: Bulk structure as StructureData (alternative to path)
        code_label: VASP code label in AiiDA
        potential_family: Potential family name
        potential_mapping: Element to potential mapping (e.g., {'Au': 'Au'})
        kpoints_spacing: K-points spacing for bulk
        clean_workdir: Clean work directory after completion
        bulk_parameters: VASP INCAR parameters for bulk
        bulk_options: Scheduler options for bulk
        miller_indices: List of Miller indices, e.g., [[1,1,1], [1,0,0], [1,1,0]]
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        lll_reduce: Use LLL reduction for slab cell
        center_slab: Center slab in c direction
        symmetrize: Generate symmetric terminations
        primitive: Find primitive cell before slab generation
        slab_parameters: VASP INCAR parameters for slab relaxation
        slab_options: Scheduler options for slab calculations
        slab_kpoints_spacing: K-points spacing for slabs
        max_concurrent_jobs: Maximum concurrent VASP calculations
        require_stoichiometric_symmetric: (EXPERIMENTAL) If True, use stoichiometric+symmetric
            surface finder to pre-filter slabs before DFT. Only slabs that are BOTH
            stoichiometric AND symmetric will be generated. Raises
            NoStoichiometricSymmetricSurfaceError if no valid surfaces found.
        stoichiometric_strategies: List of search strategies for stoichiometric finder.
            Options: 'filter_first', 'symmetrize_check', 'thickness_scan', 'bond_preservation'.
            Default: all strategies in order.
        stoichiometric_bonds: Bond dictionary for 'bond_preservation' strategy.
            Format: {('element1', 'element2'): bond_length}
        stoichiometric_auto_detect_bonds: Auto-detect bonds from covalent radii for
            the 'bond_preservation' strategy.
        stoichiometric_max_thickness: Maximum thickness for 'thickness_scan' strategy.
        name: Name for the WorkGraph

    Returns:
        WorkGraph ready for submission
        
    Outputs (accessible via task names):
        - bulk_relax: VaspWorkChain with bulk relaxation results
        - bulk_energy: extract_total_energy result
        - bulk_energy_per_atom: calculate_bulk_energy_per_atom result
        - surface_energies: Dict containing consolidated results for ALL orientations
        - wulff_shape: Dict containing Wulff shape analysis (shape_factor, anisotropy,
          weighted_surface_energy, facet_fractions, dominant_facet, etc.)
        - relaxed_slabs_hkl_XXX: Per-orientation slab structures
        - slab_energies_hkl_XXX: Per-orientation slab energies
    """
    # Validate inputs
    if bulk_structure_path is None and bulk_structure is None:
        raise ValueError("Either bulk_structure_path or bulk_structure must be provided")
    
    if miller_indices is None or len(miller_indices) == 0:
        raise ValueError("At least one Miller index must be specified")
    
    if potential_mapping is None:
        raise ValueError("potential_mapping must be provided (e.g., {'Au': 'Au'})")
    
    # Load structure from file if path provided
    if bulk_structure is None:
        if not os.path.isabs(bulk_structure_path):
            raise ValueError(f"bulk_structure_path must be absolute: {bulk_structure_path}")
        bulk_structure = load_structure_from_file(bulk_structure_path)
    
    # Load code
    code = orm.load_code(code_label)
    
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)
    
    # Build the main workflow
    wg = WorkGraph(name=name)
    
    # ===== SINGLE BULK RELAXATION =====
    bulk_task = wg.add_task(
        VaspTask,
        name='bulk_relax',
        structure=bulk_structure,
        code=code,
        parameters={'incar': bulk_parameters},
        options=bulk_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        clean_workdir=clean_workdir,
        settings=orm.Dict(dict=get_settings()),
    )
    
    bulk_energy_task = wg.add_task(
        extract_total_energy,
        name='bulk_energy',
        energies=bulk_task.outputs.misc,
    )
    
    bulk_energy_per_atom_task = wg.add_task(
        calculate_bulk_energy_per_atom,
        name='bulk_energy_per_atom',
        bulk_energy=bulk_energy_task.outputs.result,
        bulk_structure=bulk_task.outputs.structure,
    )
    
    # Set parameters for slabs
    slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
    slab_opts = slab_options if slab_options is not None else bulk_options
    slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing
    
    # To collect surface energy outputs for consolidation
    gather_inputs = {}

    # Track previous task for sequential chaining
    previous_hkl_task = None

    # ===== ADD TASK FOR EACH MILLER INDEX =====
    for miller in miller_indices:
        # Normalize to list of ints
        if isinstance(miller, tuple):
            miller = list(miller)
        miller_list = [int(m) for m in miller]
        hkl_key = miller_to_key(miller_list)
        task_name = f'surface_{hkl_key}'

        # Build task kwargs - optionally include wait_for for sequential chaining
        task_kwargs = dict(
            bulk_structure=bulk_task.outputs.structure,
            bulk_energy=bulk_energy_task.outputs.result,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            clean_workdir=clean_workdir,
            miller_indices=miller_list,
            min_slab_thickness=min_slab_thickness,
            min_vacuum_thickness=min_vacuum_thickness,
            lll_reduce=lll_reduce,
            center_slab=center_slab,
            symmetrize=symmetrize,
            primitive=primitive,
            slab_parameters=slab_params,
            slab_options=slab_opts,
            slab_kpoints_spacing=slab_kpts,
            max_concurrent_jobs=max_concurrent_jobs,
            # Stoichiometric+Symmetric Surface Finder (EXPERIMENTAL)
            require_stoichiometric_symmetric=require_stoichiometric_symmetric,
            stoichiometric_strategies=stoichiometric_strategies,
            stoichiometric_max_thickness=stoichiometric_max_thickness,
        )

        # Chain Miller index tasks sequentially to respect max_concurrent_jobs=1
        # Each HKL task waits for the previous one via data dependency
        if previous_hkl_task is not None:
            task_kwargs['wait_for'] = previous_hkl_task.outputs.surface_energies

        hkl_task = wg.add_task(
            metal_surface_energy_single_hkl,
            name=task_name,
            **task_kwargs,
        )

        # Chain: Bulk energy >> HKL task
        wg.add_link(bulk_energy_task.outputs.result, hkl_task.inputs.bulk_energy)

        previous_hkl_task = hkl_task

        # Expose individual outputs
        setattr(wg.outputs, f'relaxed_slabs_{hkl_key}', hkl_task.outputs.relaxed_slabs)
        setattr(wg.outputs, f'slab_energies_{hkl_key}', hkl_task.outputs.slab_energies)
        setattr(wg.outputs, f'slab_structures_{hkl_key}', hkl_task.outputs.slab_structures)

        # Collect for gathered output
        gather_inputs[f'surface_energies_{hkl_key}'] = hkl_task.outputs.surface_energies
    
    # ===== GATHER SURFACE ENERGIES =====
    gather_task = wg.add_task(
        gather_surface_energies,
        name='gather_surface_energies',
        **gather_inputs
    )

    # ===== BUILD WULFF SHAPE =====
    wulff_task = wg.add_task(
        build_wulff_shape,
        name='build_wulff_shape',
        bulk_structure=bulk_task.outputs.structure,
        surface_energies=gather_task.outputs.result,
    )

    # Expose consolidated output
    wg.outputs.surface_energies = gather_task.outputs.result
    wg.outputs.wulff_shape = wulff_task.outputs.result
    
    # Expose bulk outputs
    wg.outputs.bulk_energy = bulk_energy_task.outputs.result
    wg.outputs.bulk_energy_per_atom = bulk_energy_per_atom_task.outputs.result
    wg.outputs.bulk_structure = bulk_task.outputs.structure
    
    # Set concurrency limit
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs
    
    return wg
