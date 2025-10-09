"""
PS-TEROS Core WorkGraph

This module contains the core workflow for PS-TEROS calculations using
the pythonic scatter-gather pattern for parallel execution.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, dynamic, namespace
from ase.io import read
from teros.core.hf import calculate_formation_enthalpy
from teros.core.slabs import (
    generate_slab_structures,
    relax_slabs_scatter,
    extract_total_energy,
)
from teros.core.thermodynamics import (
    identify_oxide_type,
    compute_surface_energies_scatter,
)
from teros.core.cleavage import compute_cleavage_energies_scatter

def load_structure_from_file(filepath: str) -> orm.StructureData:
    """
    Load structure from a file using ASE (helper function for graph construction).

    Args:
        filepath: Path to structure file

    Returns:
        StructureData node
    """
    atoms = read(filepath)
    return orm.StructureData(ase=atoms)


def get_settings():
    """
    Parser settings for aiida-vasp.
    """
    return {
        'parser_settings': {
            'add_energy': True,
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }

# NOTE: Using scatter-gather pattern for slab relaxations
# All VASP tasks created in parallel within @task.graph

@task.graph(outputs=[
    'bulk_energy', 'metal_energy', 'nonmetal_energy', 'oxygen_energy',
    'bulk_structure', 'metal_structure', 'nonmetal_structure', 'oxygen_structure',
    'formation_enthalpy', 'slab_structures', 'relaxed_slabs', 'slab_energies',
    'surface_energies', 'cleavage_energies',
])
def core_workgraph(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict,
    metal_potential_mapping: dict,
    nonmetal_potential_mapping: dict,
    oxygen_potential_mapping: dict,
    kpoints_spacing: float,
    bulk_parameters: dict,
    bulk_options: dict,
    metal_parameters: dict,
    metal_options: dict,
    nonmetal_parameters: dict,
    nonmetal_options: dict,
    oxygen_parameters: dict,
    oxygen_options: dict,
    clean_workdir: bool,
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = False,
    center_slab: bool = True,
    symmetrize: bool = False,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = False,
    compute_thermodynamics: bool = False,
    thermodynamics_sampling: int = 100,
    compute_cleavage: bool = False,
):
    """
    Core WorkGraph for formation enthalpy calculations of ternary oxides with slab generation.

    This workflow relaxes the bulk compound and all reference elements in parallel,
    then calculates the formation enthalpy and generates slab structures from the relaxed bulk.
    Optionally relaxes slabs and computes surface energies with thermodynamics.
    The workflow is general and works for any ternary oxide system (e.g., Ag3PO4, Fe2WO6, etc.).
    Each reference (metal, nonmetal, oxygen) can have its own specific calculation parameters.

    Args:
        structures_dir: Directory containing all structure files
        bulk_name: Filename of bulk structure (e.g., 'ag3po4.cif')
        metal_name: Filename of metal reference structure (e.g., 'Ag.cif')
        nonmetal_name: Filename of nonmetal reference structure (e.g., 'P.cif')
        oxygen_name: Filename of oxygen reference structure (e.g., 'O2.cif')
        code_label: Label of the VASP code in AiiDA
        potential_family: Name of the potential family
        bulk_potential_mapping: Mapping for bulk (e.g., {'Ag': 'Ag', 'P': 'P', 'O': 'O'})
        metal_potential_mapping: Mapping for metal (e.g., {'Ag': 'Ag'})
        nonmetal_potential_mapping: Mapping for nonmetal (e.g., {'P': 'P'})
        oxygen_potential_mapping: Mapping for oxygen (e.g., {'O': 'O'})
        kpoints_spacing: K-points spacing in A^-1 * 2pi
        bulk_parameters: VASP parameters for bulk
        bulk_options: Scheduler options for bulk
        metal_parameters: VASP parameters for metal reference
        metal_options: Scheduler options for metal reference
        nonmetal_parameters: VASP parameters for nonmetal reference
        nonmetal_options: Scheduler options for nonmetal reference
        oxygen_parameters: VASP parameters for oxygen reference
        oxygen_options: Scheduler options for oxygen reference
        clean_workdir: Whether to clean work directory
        miller_indices: Miller indices for slab generation (e.g., [1, 0, 0])
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        slab_parameters: VASP parameters for slab relaxation. Default: None (uses bulk_parameters)
        slab_options: Scheduler options for slab calculations. Default: None (uses bulk_options)
        slab_potential_mapping: Potential mapping for slabs. Default: None (uses bulk_potential_mapping)
        slab_kpoints_spacing: K-points spacing for slabs. Default: None (uses kpoints_spacing)
        lll_reduce: Reduce cell using LLL algorithm. Default: False
        center_slab: Center the slab in c direction. Default: True
        symmetrize: Generate symmetrically distinct terminations. Default: False
        primitive: Find primitive cell before slab generation. Default: True
        in_unit_planes: Restrict Miller indices to unit planes. Default: False
        max_normal_search: Max normal search for Miller indices. Default: None
        relax_slabs: Whether to relax the generated slabs with VASP. Default: False
        compute_thermodynamics: Whether to compute surface energies. Default: False (requires relax_slabs=True)
        thermodynamics_sampling: Number of grid points for chemical potential sampling. Default: 100
        compute_cleavage: Whether to compute cleavage energies for complementary slab pairs. Default: False (requires relax_slabs=True)

    Returns:
        Dictionary with energies, structures, formation enthalpy, slab structures, and optionally surface energies:
            - bulk_energy, metal_energy, nonmetal_energy, oxygen_energy: Total energies (Float)
            - bulk_structure, metal_structure, nonmetal_structure, oxygen_structure: Relaxed structures (StructureData)
            - formation_enthalpy: Formation enthalpy and related data (Dict)
            - slab_structures: Dynamic namespace with all generated slab terminations (StructureData)
            - relaxed_slabs: Dynamic namespace with relaxed slab structures (if relax_slabs=True)
            - slab_energies: Dynamic namespace with slab energies (if relax_slabs=True)
            - surface_energies: Dynamic namespace with surface energy data (if compute_thermodynamics=True)
            - cleavage_energies: Dynamic namespace with cleavage energy data (if compute_cleavage=True)
    """
    # Load the code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap it as a task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # ===== BULK RELAXATION =====
    bulk_struct = load_structure_from_file(f"{structures_dir}/{bulk_name}")

    bulk_vasp = VaspTask(
        structure=bulk_struct,
        code=code,
        parameters={'incar': bulk_parameters},
        options=bulk_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=bulk_potential_mapping,
        clean_workdir=clean_workdir,
        settings=orm.Dict(dict=get_settings()),
    )

    bulk_energy = extract_total_energy(energies=bulk_vasp.misc)

    # ===== METAL RELAXATION =====
    metal_struct = load_structure_from_file(f"{structures_dir}/{metal_name}")

    metal_vasp = VaspTask(
        structure=metal_struct,
        code=code,
        parameters={'incar': metal_parameters},
        options=metal_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=metal_potential_mapping,
        clean_workdir=clean_workdir,
        settings=orm.Dict(dict=get_settings()),
    )

    metal_energy = extract_total_energy(energies=metal_vasp.misc)

    # ===== NONMETAL RELAXATION =====
    nonmetal_struct = load_structure_from_file(f"{structures_dir}/{nonmetal_name}")

    nonmetal_vasp = VaspTask(
        structure=nonmetal_struct,
        code=code,
        parameters={'incar': nonmetal_parameters},
        options=nonmetal_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=nonmetal_potential_mapping,
        clean_workdir=clean_workdir,
        settings=orm.Dict(dict=get_settings()),
    )

    nonmetal_energy = extract_total_energy(energies=nonmetal_vasp.misc)

    # ===== OXYGEN RELAXATION =====
    oxygen_struct = load_structure_from_file(f"{structures_dir}/{oxygen_name}")

    oxygen_vasp = VaspTask(
        structure=oxygen_struct,
        code=code,
        parameters={'incar': oxygen_parameters},
        options=oxygen_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=oxygen_potential_mapping,
        clean_workdir=clean_workdir,
        settings=orm.Dict(dict=get_settings()),
    )

    oxygen_energy = extract_total_energy(energies=oxygen_vasp.misc)

    # ===== FORMATION ENTHALPY CALCULATION =====
    formation_hf = calculate_formation_enthalpy(
        bulk_structure=bulk_vasp.structure,
        bulk_energy=bulk_energy.result,
        metal_structure=metal_vasp.structure,
        metal_energy=metal_energy.result,
        nonmetal_structure=nonmetal_vasp.structure,
        nonmetal_energy=nonmetal_energy.result,
        oxygen_structure=oxygen_vasp.structure,
        oxygen_energy=oxygen_energy.result,
    )

    # ===== SLAB GENERATION =====
    # Generate slabs using the pythonic scatter-gather pattern
    slab_namespace = generate_slab_structures(
        bulk_structure=bulk_vasp.structure,
        miller_indices=orm.List(list=miller_indices),
        min_slab_thickness=orm.Float(min_slab_thickness),
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        symmetrize=orm.Bool(symmetrize),
        primitive=orm.Bool(primitive),
    ).slabs

    # ===== SLAB RELAXATION (OPTIONAL) =====
    relaxed_slabs_output = None
    slab_energies_output = None

    if relax_slabs:
        # Use slab-specific parameters or fall back to bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        # Use scatter-gather pattern for parallel slab relaxation
        relaxation_outputs = relax_slabs_scatter(
            slabs=slab_namespace,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            parameters=slab_params,
            options=slab_opts,
            kpoints_spacing=slab_kpts,
            clean_workdir=clean_workdir,
        )

        relaxed_slabs_output = relaxation_outputs.relaxed_structures
        slab_energies_output = relaxation_outputs.energies

    # ===== THERMODYNAMICS CALCULATION (OPTIONAL) =====
    surface_energies_output = None
    if compute_thermodynamics and relax_slabs:
        # Identify oxide type (binary or ternary)
        oxide_type_result = identify_oxide_type(bulk_structure=bulk_vasp.structure)
        
        # Compute surface energies for all slabs
        # The formation_hf.result Dict contains all needed reference energies
        surface_outputs = compute_surface_energies_scatter(
            slabs=relaxed_slabs_output,
            energies=slab_energies_output,
            bulk_structure=bulk_vasp.structure,
            bulk_energy=bulk_energy.result,
            reference_energies=formation_hf.result,  # Dict with reference energies
            formation_enthalpy=formation_hf.result,  # Dict with formation_enthalpy_ev key
            oxide_type=oxide_type_result.result,  # Pass the Str node directly
            sampling=thermodynamics_sampling,
        )
        
        surface_energies_output = surface_outputs.surface_energies

    # ===== CLEAVAGE ENERGY CALCULATION (OPTIONAL) =====
    cleavage_energies_output = None
    if compute_cleavage and relax_slabs:
        # Compute cleavage energies for all complementary slab pairs
        cleavage_outputs = compute_cleavage_energies_scatter(
            slabs=relaxed_slabs_output,
            energies=slab_energies_output,
            bulk_structure=bulk_vasp.structure,
            bulk_energy=bulk_energy.result,
        )
        
        cleavage_energies_output = cleavage_outputs.cleavage_energies

    # Return all outputs
    return {
        'bulk_energy': bulk_energy.result,
        'bulk_structure': bulk_vasp.structure,
        'metal_energy': metal_energy.result,
        'metal_structure': metal_vasp.structure,
        'nonmetal_energy': nonmetal_energy.result,
        'nonmetal_structure': nonmetal_vasp.structure,
        'oxygen_energy': oxygen_energy.result,
        'oxygen_structure': oxygen_vasp.structure,
        'formation_enthalpy': formation_hf.result,
        'slab_structures': slab_namespace,
        'relaxed_slabs': relaxed_slabs_output,
        'slab_energies': slab_energies_output,
        'surface_energies': surface_energies_output,
        'cleavage_energies': cleavage_energies_output,
    }


def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    metal_potential_mapping: dict = None,
    nonmetal_potential_mapping: dict = None,
    oxygen_potential_mapping: dict = None,
    kpoints_spacing: float = 0.3,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    clean_workdir: bool = True,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = False,
    center_slab: bool = True,
    symmetrize: bool = False,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = False,
    compute_thermodynamics: bool = False,
    thermodynamics_sampling: int = 100,
    compute_cleavage: bool = False,
    name: str = 'FormationEnthalpy',
):
    """
    Build a formation enthalpy WorkGraph for any ternary oxide with slab generation and optional relaxation.

    This is a convenience wrapper that builds and returns a WorkGraph
    ready to calculate formation enthalpy, generate slab structures, and optionally relax them.

    Args:
        relax_slabs: Whether to relax generated slabs with VASP. Default: False
        compute_thermodynamics: Whether to compute surface energies. Default: False
        thermodynamics_sampling: Grid resolution for chemical potential sampling. Default: 100
        compute_cleavage: Whether to compute cleavage energies. Default: False
        slab_parameters: VASP parameters for slab relaxation. Default: None (uses bulk_parameters)
        slab_options: Scheduler options for slab calculations. Default: None (uses bulk_options)
        slab_potential_mapping: Potential mapping for slabs. Default: None (uses bulk_potential_mapping)
        slab_kpoints_spacing: K-points spacing for slabs. Default: None (uses kpoints_spacing)
        ... (other parameters as before)

    Returns:
        WorkGraph instance ready to be submitted
    """
    # Build the workgraph
    # Note: parameters will be wrapped with {'incar': ...} inside the graph
    wg = core_workgraph.build(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        metal_name=metal_name,
        nonmetal_name=nonmetal_name,
        oxygen_name=oxygen_name,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping or {},
        metal_potential_mapping=metal_potential_mapping or {},
        nonmetal_potential_mapping=nonmetal_potential_mapping or {},
        oxygen_potential_mapping=oxygen_potential_mapping or {},
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters or {},
        bulk_options=bulk_options or {},
        metal_parameters=metal_parameters or {},
        metal_options=metal_options or {},
        nonmetal_parameters=nonmetal_parameters or {},
        nonmetal_options=nonmetal_options or {},
        oxygen_parameters=oxygen_parameters or {},
        oxygen_options=oxygen_options or {},
        clean_workdir=clean_workdir,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        relax_slabs=relax_slabs,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        compute_cleavage=compute_cleavage,
    )

    # Set the name
    wg.name = name

    return wg


def build_core_workgraph_with_map(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    metal_potential_mapping: dict = None,
    nonmetal_potential_mapping: dict = None,
    oxygen_potential_mapping: dict = None,
    kpoints_spacing: float = 0.3,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    clean_workdir: bool = True,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = False,
    center_slab: bool = True,
    symmetrize: bool = False,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = False,
    compute_thermodynamics: bool = False,
    thermodynamics_sampling: int = 100,
    compute_cleavage: bool = False,
    name: str = 'FormationEnthalpy_ScatterGather',
) -> WorkGraph:
    """
    DEPRECATED: Now uses scatter-gather pattern instead of Map zone.

    Build a formation enthalpy WorkGraph using scatter-gather pattern for slab relaxations.

    This function creates a WorkGraph that:
    1. Relaxes bulk and reference structures in parallel
    2. Calculates formation enthalpy
    3. Generates slab structures from relaxed bulk
    4. Uses scatter-gather pattern (@task.graph) to relax all slabs in parallel (if relax_slabs=True)

    Note: This function is now an alias for build_core_workgraph with updated
    implementation using the pythonic scatter-gather pattern instead of Map zones.

    Args:
        Same as build_core_workgraph

    Returns:
        WorkGraph instance ready to be submitted
    """
    # Simply forward to build_core_workgraph which now uses scatter-gather
    return build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        metal_name=metal_name,
        nonmetal_name=nonmetal_name,
        oxygen_name=oxygen_name,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping,
        metal_potential_mapping=metal_potential_mapping,
        nonmetal_potential_mapping=nonmetal_potential_mapping,
        oxygen_potential_mapping=oxygen_potential_mapping,
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=clean_workdir,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        relax_slabs=relax_slabs,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        compute_cleavage=compute_cleavage,
        name=name,
    )


if __name__ == '__main__':
    """
    Simple test/example of the WorkGraph.
    For full examples, see examples/formation/
    """
    print("PS-TEROS Core WorkGraph Module")
    print("=" * 50)
    print("\nThis creates a WorkGraph using @task.graph decorator:")
    print("  - load_structure tasks (using @task)")
    print("  - VASP relaxation tasks (run in parallel)")
    print("  - extract_total_energy tasks (using @task)")
    print("\nFor examples, see:")
    print("  - examples/formation/formation.py (formation enthalpy)")
