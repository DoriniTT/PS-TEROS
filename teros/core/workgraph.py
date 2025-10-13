"""
PS-TEROS Core WorkGraph

This module contains the core workflow for PS-TEROS calculations using
the pythonic scatter-gather pattern for parallel execution.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, dynamic, namespace
from aiida.engine import submit
from ase.io import read
from teros.core.hf import calculate_formation_enthalpy
from teros.core.slabs import (
    generate_slab_structures,
    wrap_input_slabs,
    scf_slabs_scatter,
    relax_slabs_scatter,
    calculate_relaxation_energies_scatter,
    extract_total_energy,
)
from teros.core.thermodynamics import (
    identify_oxide_type,
    compute_surface_energies_scatter,
)
from teros.core.cleavage import compute_cleavage_energies_scatter
from teros.core.aimd import aimd_slabs_scatter

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
    'unrelaxed_slab_energies', 'relaxation_energies',
    'surface_energies', 'cleavage_energies',
    'slab_remote', 'unrelaxed_slab_remote',
    # Electronic properties outputs
    'bulk_bands', 'bulk_dos', 'bulk_primitive_structure', 'bulk_seekpath_parameters',
    'aimd_results',
])
def core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
    bulk_parameters: dict,
    bulk_options: dict,
    clean_workdir: bool,
    metal_name: str = None,
    oxygen_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    miller_indices: list = None,
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,  # Calculate relaxation energy by default
    input_slabs: dict = None,
    use_input_slabs: bool = False,  # Signal to skip slab generation
    compute_cleavage: bool = True,
    run_aimd: bool = False,
    aimd_sequence: list = None,
    aimd_parameters: dict = None,
    aimd_options: dict = None,
    aimd_potential_mapping: dict = None,
    aimd_kpoints_spacing: float = None,
    # Note: Electronic properties parameters removed - handled in build_core_workgraph
):
    """
    Core WorkGraph for formation enthalpy calculations of binary and ternary oxides with optional slab calculations.

    This workflow relaxes the bulk compound and all reference elements in parallel,
    then calculates the formation enthalpy. Optionally generates and relaxes slab structures,
    and computes surface energies, cleavage energies, and relaxation energies based on flags.

    The workflow is general and works for any binary oxide (e.g., Ag2O, CuO) or ternary oxide
    system (e.g., Ag3PO4, Fe2WO6, etc.). For binary oxides, nonmetal parameters should be None.
    Each reference (metal, nonmetal, oxygen) can have its own specific calculation parameters.

    Calculations controlled by boolean flags (all default to True):
        - compute_relaxation_energy: Calculate energy difference between unrelaxed and relaxed slabs
        - compute_cleavage: Calculate cleavage energies for complementary slab pairs
        - compute_thermodynamics: Calculate surface energies with chemical potential sampling
          (requires metal_name and oxygen_name to be provided)

    Note: Electronic properties (DOS/bands) are handled in build_core_workgraph, not here

    Args:
        structures_dir: Directory containing all structure files
        bulk_name: Filename of bulk structure (e.g., 'ag2o.cif' or 'ag3po4.cif')
        code_label: Label of the VASP code in AiiDA
        potential_family: Name of the potential family
        bulk_potential_mapping: Mapping for bulk (e.g., {'Ag': 'Ag', 'O': 'O'} for binary or {'Ag': 'Ag', 'P': 'P', 'O': 'O'} for ternary)
        kpoints_spacing: K-points spacing in A^-1 * 2pi
        bulk_parameters: VASP parameters for bulk
        bulk_options: Scheduler options for bulk
        clean_workdir: Whether to clean work directory
        metal_name: Filename of metal reference structure (e.g., 'Ag.cif'). Optional, default: None
        oxygen_name: Filename of oxygen reference structure (e.g., 'O2.cif'). Optional, default: None
        metal_potential_mapping: Mapping for metal (e.g., {'Ag': 'Ag'}). Optional, default: None
        metal_parameters: VASP parameters for metal reference. Optional, default: None
        metal_options: Scheduler options for metal reference. Optional, default: None
        oxygen_potential_mapping: Mapping for oxygen (e.g., {'O': 'O'}). Optional, default: None
        oxygen_parameters: VASP parameters for oxygen reference. Optional, default: None
        oxygen_options: Scheduler options for oxygen reference. Optional, default: None
        nonmetal_name: Filename of nonmetal reference structure (e.g., 'P.cif'). Optional, default: None
        nonmetal_potential_mapping: Mapping for nonmetal (e.g., {'P': 'P'}). Optional, default: None
        nonmetal_parameters: VASP parameters for nonmetal reference. Optional, default: None
        nonmetal_options: Scheduler options for nonmetal reference. Optional, default: None
        miller_indices: Miller indices for slab generation (e.g., [1, 0, 0]). Optional, default: None (skips slab generation if None and input_slabs not provided)
        min_slab_thickness: Minimum slab thickness in Angstroms. Default: 15.0
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms. Default: 15.0
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
        compute_thermodynamics: Whether to compute surface energies. Default: True (requires metal_name and oxygen_name)
        thermodynamics_sampling: Number of grid points for chemical potential sampling. Default: 100
        compute_relaxation_energy: Whether to compute relaxation energies. Default: True (requires relax_slabs=True)
        input_slabs: Dictionary of pre-generated slab structures (e.g., {'term_0': StructureData, ...}).
                     If provided, slab generation is skipped. Default: None
        compute_cleavage: Whether to compute cleavage energies for complementary slab pairs. Default: True (requires relax_slabs=True)

    Returns:
        Dictionary with energies, structures, formation enthalpy, slab structures, and optionally surface energies:
            - bulk_energy, metal_energy, nonmetal_energy, oxygen_energy: Total energies (Float)
            - bulk_structure, metal_structure, nonmetal_structure, oxygen_structure: Relaxed structures (StructureData)
            - formation_enthalpy: Formation enthalpy and related data (Dict)
            - slab_structures: Dynamic namespace with all generated slab terminations (StructureData)
            - relaxed_slabs: Dynamic namespace with relaxed slab structures (if relax_slabs=True)
            - slab_energies: Dynamic namespace with relaxed slab energies (if relax_slabs=True)
            - unrelaxed_slab_energies: Dynamic namespace with unrelaxed slab energies from SCF (if relax_slabs=True)
            - relaxation_energies: Dynamic namespace with relaxation energies (E_relaxed - E_unrelaxed) (if relax_slabs=True)
            - slab_remote: Dynamic namespace with RemoteData nodes for each slab relaxation (if relax_slabs=True)
            - unrelaxed_slab_remote: Dynamic namespace with RemoteData nodes for each SCF calculation (if relax_slabs=True)
            - surface_energies: Dynamic namespace with surface energy data (if compute_thermodynamics=True)
            - cleavage_energies: Dynamic namespace with cleavage energy data (if compute_cleavage=True)

    Note: Electronic properties outputs (bulk_bands, bulk_dos, bulk_electronic_properties_misc)
          are added in build_core_workgraph when compute_electronic_properties_bulk=True
    """
    # Load the code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap it as a task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Determine if we should compute formation enthalpy
    # (requires both metal_name and oxygen_name)
    compute_formation_enthalpy = (metal_name is not None and oxygen_name is not None)

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

    # ===== REFERENCE RELAXATIONS AND FORMATION ENTHALPY (CONDITIONAL) =====
    if compute_formation_enthalpy:
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

        # ===== NONMETAL RELAXATION (ONLY FOR TERNARY OXIDES) =====
        if nonmetal_name is not None:
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
        else:
            # For binary oxides, use dummy values (will be ignored by calculate_formation_enthalpy)
            nonmetal_vasp = type('obj', (object,), {'structure': metal_vasp.structure})()
            nonmetal_energy = type('obj', (object,), {'result': metal_energy.result})()

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
    else:
        # Bulk-only mode: Create dummy outputs with proper types
        # Use bulk structure as placeholder for reference structures
        metal_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
        metal_energy = type('obj', (object,), {'result': bulk_energy.result})()
        nonmetal_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
        nonmetal_energy = type('obj', (object,), {'result': bulk_energy.result})()
        oxygen_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
        oxygen_energy = type('obj', (object,), {'result': bulk_energy.result})()
        formation_hf = type('obj', (object,), {'result': bulk_energy.result})()  # Placeholder

    # ===== ELECTRONIC PROPERTIES CALCULATION (OPTIONAL) =====
    # Note: Electronic properties are added in build_core_workgraph, not here
    # This keeps core_workgraph focused on base calculations

    # ===== SLAB GENERATION =====
    # Use input slabs if provided, otherwise generate them
    if use_input_slabs:
        # User will provide pre-generated slabs after building
        # Skip slab generation entirely
        slab_namespace = None  # Will be set post-build
    elif miller_indices is None:
        # No miller_indices provided and no input_slabs - skip slab generation
        slab_namespace = None
    else:
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
    relaxed_slabs_output = {}
    slab_energies_output = {}
    slab_remote_output = {}
    unrelaxed_slab_energies_output = {}
    unrelaxed_slab_remote_output = {}
    relaxation_energies_output = {}

    if relax_slabs and slab_namespace is not None:
        # Use slab-specific parameters or fall back to bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        # OPTIONAL: SCF calculation on unrelaxed slabs (only if relaxation energy requested)
        if compute_relaxation_energy:
            scf_outputs = scf_slabs_scatter(
                slabs=slab_namespace,
                code=code,
                potential_family=potential_family,
                potential_mapping=slab_pot_map,
                parameters=slab_params,
                options=slab_opts,
                kpoints_spacing=slab_kpts,
                clean_workdir=clean_workdir,
            )
            unrelaxed_slab_energies_output = scf_outputs.energies
            unrelaxed_slab_remote_output = scf_outputs.remote_folders

        # ALWAYS: Relax slabs
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
        slab_remote_output = relaxation_outputs.remote_folders

        # OPTIONAL: Calculate relaxation energies (only if SCF was done)
        if compute_relaxation_energy:
            relax_energy_outputs = calculate_relaxation_energies_scatter(
                unrelaxed_energies=unrelaxed_slab_energies_output,
                relaxed_energies=slab_energies_output,
            )
            relaxation_energies_output = relax_energy_outputs.relaxation_energies

    # ===== THERMODYNAMICS CALCULATION (OPTIONAL) =====
    surface_energies_output = {}
    if compute_thermodynamics and relax_slabs and slab_namespace is not None:
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
    cleavage_energies_output = {}
    if compute_cleavage and relax_slabs and not use_input_slabs:
        # Only compute cleavage inside core_workgraph if slabs were generated/relaxed here
        # When use_input_slabs=True, cleavage will be added in build_core_workgraph
        cleavage_outputs = compute_cleavage_energies_scatter(
            slabs=relaxed_slabs_output,
            energies=slab_energies_output,
            bulk_structure=bulk_vasp.structure,
            bulk_energy=bulk_energy.result,
        )

        cleavage_energies_output = cleavage_outputs.cleavage_energies

    # ===== AIMD CALCULATION (OPTIONAL) =====
    aimd_outputs = {}
    if run_aimd and slab_namespace is not None:
        # Use AIMD-specific parameters or fall back to slab parameters
        aimd_params = aimd_parameters if aimd_parameters is not None else slab_params
        aimd_opts = aimd_options if aimd_options is not None else slab_opts
        aimd_pot_map = aimd_potential_mapping if aimd_potential_mapping is not None else slab_pot_map
        aimd_kpts = aimd_kpoints_spacing if aimd_kpoints_spacing is not None else slab_kpts

        # Run AIMD on all slabs in parallel
        aimd_results = aimd_slabs_scatter(
            slabs=slab_namespace,
            aimd_sequence=aimd_sequence,
            code=code,
            aimd_parameters=aimd_params,
            potential_family=potential_family,
            potential_mapping=aimd_pot_map,
            options=aimd_opts,
            kpoints_spacing=aimd_kpts,
            clean_workdir=clean_workdir,
        )

        aimd_outputs = aimd_results

    # Return all outputs
    # Note: Electronic properties outputs (bulk_bands, bulk_dos, bulk_electronic_properties_misc)
    # are added dynamically in build_core_workgraph, not returned here
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
        'slab_structures': slab_namespace if slab_namespace is not None else {},
        'relaxed_slabs': relaxed_slabs_output,
        'slab_energies': slab_energies_output,
        'unrelaxed_slab_energies': unrelaxed_slab_energies_output,
        'relaxation_energies': relaxation_energies_output,
        'surface_energies': surface_energies_output,
        'cleavage_energies': cleavage_energies_output,
        'slab_remote': slab_remote_output,
        'unrelaxed_slab_remote': unrelaxed_slab_remote_output,
        'aimd_results': aimd_outputs,
    }


def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,
    metal_name: str = None,
    oxygen_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    miller_indices: list = None,
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,  # Calculate relaxation energy by default
    input_slabs: dict = None,
    compute_cleavage: bool = True,
    restart_from_node: int = None,  # PK of previous workgraph to restart from
    compute_electronic_properties_bulk: bool = False,  # NEW: Enable DOS and bands
    bands_parameters: dict = None,  # NEW: INCAR parameters
    bands_options: dict = None,  # NEW: Scheduler options
    band_settings: dict = None,  # NEW: Band workflow settings
    run_aimd: bool = False,
    aimd_sequence: list = None,
    aimd_parameters: dict = None,
    aimd_options: dict = None,
    aimd_potential_mapping: dict = None,
    aimd_kpoints_spacing: float = None,
    name: str = 'FormationEnthalpy',
):
    """
    Build a centralized WorkGraph for bulk relaxation, formation enthalpy, and surface calculations.

    This function builds a flexible workflow that can handle:
    - Bulk structure relaxation (always performed)
    - Reference structure relaxation (if metal_name and oxygen_name provided)
    - Formation enthalpy calculation (if references provided)
    - Slab generation and relaxation (if miller_indices or input_slabs provided)
    - Relaxation energy calculation (controlled by compute_relaxation_energy flag, default: True)
    - Cleavage energy calculation (controlled by compute_cleavage flag, default: True)
    - Surface thermodynamics (controlled by compute_thermodynamics flag, default: True, requires references)
    - Electronic properties calculation (controlled by compute_electronic_properties_bulk flag, default: False)

    Args:
        structures_dir: Directory containing structure files
        bulk_name: Filename of bulk structure (e.g., 'ag3po4.cif')
        code_label: VASP code label in AiiDA. Default: 'VASP-VTST-6.4.3@bohr'
        potential_family: Potential family name. Default: 'PBE'
        bulk_potential_mapping: Element to potential mapping for bulk (e.g., {'Ag': 'Ag', 'P': 'P', 'O': 'O'})
        kpoints_spacing: K-points spacing in A^-1 * 2pi. Default: 0.3
        bulk_parameters: VASP INCAR parameters for bulk relaxation
        bulk_options: Scheduler options for bulk calculation
        clean_workdir: Clean work directory after completion. Default: False

        metal_name: Metal reference filename (e.g., 'Ag.cif'). Default: None
        oxygen_name: Oxygen reference filename (e.g., 'O2.cif'). Default: None
        metal_potential_mapping: Potential mapping for metal. Required if metal_name provided.
        metal_parameters: VASP parameters for metal. Required if metal_name provided.
        metal_options: Scheduler options for metal. Required if metal_name provided.
        oxygen_potential_mapping: Potential mapping for oxygen. Required if oxygen_name provided.
        oxygen_parameters: VASP parameters for oxygen. Required if oxygen_name provided.
        oxygen_options: Scheduler options for oxygen. Required if oxygen_name provided.
        nonmetal_name: Nonmetal reference filename for ternary oxides (e.g., 'P.cif'). Default: None
        nonmetal_potential_mapping: Potential mapping for nonmetal. Required if nonmetal_name provided.
        nonmetal_parameters: VASP parameters for nonmetal. Required if nonmetal_name provided.
        nonmetal_options: Scheduler options for nonmetal. Required if nonmetal_name provided.

        miller_indices: Miller indices for slab generation (e.g., [1, 0, 0]). Optional, default: None (skips slab generation if None and input_slabs not provided)
        min_slab_thickness: Minimum slab thickness in Angstroms. Default: 15.0
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms. Default: 15.0
        slab_parameters: VASP parameters for slab relaxation. Default: None (uses bulk_parameters)
        slab_options: Scheduler options for slabs. Default: None (uses bulk_options)
        slab_potential_mapping: Potential mapping for slabs. Default: None (uses bulk_potential_mapping)
        slab_kpoints_spacing: K-points spacing for slabs. Default: None (uses kpoints_spacing)

        relax_slabs: Relax generated slabs. Default: False
        compute_thermodynamics: Compute surface energies (requires metal_name and oxygen_name). Default: True
        thermodynamics_sampling: Grid resolution for chemical potential sampling. Default: 100
        compute_relaxation_energy: Compute relaxation energies (E_relaxed - E_unrelaxed). Default: True
        compute_cleavage: Compute cleavage energies. Default: True
        input_slabs: Pre-generated slab structures dict. Skips slab generation if provided. Default: None
        restart_from_node: PK of previous workgraph to restart from (extracts slabs and RemoteData). Default: None
        compute_electronic_properties_bulk: Compute DOS and bands for relaxed bulk. Default: False
        bands_parameters: Dict with 'scf', 'bands', 'dos' keys containing INCAR dicts, plus
                          optional 'scf_kpoints_distance'. Use get_electronic_properties_defaults().
                          Default: None
        bands_options: Scheduler options for bands calculation. Default: None (uses bulk_options)
        band_settings: Dict with band workflow settings. Use get_electronic_properties_defaults()
                       ['band_settings']. Default: None
        run_aimd: Run AIMD on generated/input slabs. Default: False
        aimd_sequence: List of AIMD stages [{'temperature': K, 'steps': N}, ...]. Default: None
        aimd_parameters: AIMD INCAR parameters (use get_aimd_defaults()). Default: None
        aimd_options: Scheduler options for AIMD. Default: None (uses slab_options)
        aimd_potential_mapping: Potential mapping for AIMD. Default: None (uses slab_potential_mapping)
        aimd_kpoints_spacing: K-points spacing for AIMD. Default: None (uses slab_kpoints_spacing)
        name: WorkGraph name. Default: 'FormationEnthalpy'

    Returns:
        WorkGraph instance ready to be submitted

    Examples:
        Bulk-only relaxation:
        >>> wg = build_core_workgraph(
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520, ...},
        ...     bulk_options={'resources': {'num_machines': 1}},
        ... )

        Formation enthalpy with all calculations:
        >>> wg = build_core_workgraph(
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     metal_name='Ag.cif',
        ...     oxygen_name='O2.cif',
        ...     nonmetal_name='P.cif',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     metal_potential_mapping={'Ag': 'Ag'},
        ...     oxygen_potential_mapping={'O': 'O'},
        ...     nonmetal_potential_mapping={'P': 'P'},
        ...     bulk_parameters={...},
        ...     metal_parameters={...},
        ...     oxygen_parameters={...},
        ...     nonmetal_parameters={...},
        ... )
    """
    # Handle restart from previous workgraph
    restart_folders = None
    restart_slabs = None
    if restart_from_node is not None:
        from teros.core.slabs import extract_restart_folders_from_node
        print(f"\n{'='*70}")
        print(f"RESTART MODE: Loading data from node {restart_from_node}")
        print(f"{'='*70}")

        # Load the previous workgraph node
        try:
            prev_node = orm.load_node(restart_from_node)

            # Extract restart folders (RemoteData)
            restart_folders = extract_restart_folders_from_node(restart_from_node)
            print(f"  ✓ Extracted restart folders: {list(restart_folders.keys())}")

            # Extract slab structures from previous run
            if hasattr(prev_node.outputs, 'slab_structures'):
                restart_slabs = {}
                for label in prev_node.outputs.slab_structures.keys():
                    restart_slabs[label] = prev_node.outputs.slab_structures[label]
                print(f"  ✓ Extracted slab structures: {list(restart_slabs.keys())}")
            else:
                print(f"  ! Warning: Previous node has no slab_structures, will generate new slabs")

            # Override input_slabs with slabs from previous run
            if restart_slabs:
                input_slabs = restart_slabs
                print(f"  → Using slabs from previous run")

            print(f"{'='*70}\n")

        except ValueError as e:
            print(f"  ✗ Error extracting restart data: {e}")
            print(f"  → Proceeding without restart\n")
            restart_folders = None
            restart_slabs = None

    # Special handling for input_slabs: stored nodes can't be passed through @task.graph
    use_input_slabs = input_slabs is not None and len(input_slabs) > 0

    # Automatically disable cleavage calculation when using manual slabs
    if use_input_slabs and compute_cleavage:
        print(f"\n{'='*70}")
        print(f"WARNING: Cleavage calculation disabled")
        print(f"  Reason: Manual slabs (input_slabs) provided")
        print(f"  Cleavage energies require complementary slab pairs from")
        print(f"  automatic slab generation, which is bypassed when using")
        print(f"  manually provided slabs.")
        print(f"{'='*70}\n")
        compute_cleavage = False

    # Build the workgraph (pass None for input_slabs to avoid serialization issues)
    # Note: parameters will be wrapped with {'incar': ...} inside the graph
    wg = core_workgraph.build(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping or {},
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters or {},
        bulk_options=bulk_options or {},
        clean_workdir=clean_workdir,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        metal_potential_mapping=metal_potential_mapping,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        oxygen_potential_mapping=oxygen_potential_mapping,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        nonmetal_name=nonmetal_name,
        nonmetal_potential_mapping=nonmetal_potential_mapping,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        miller_indices=miller_indices if not use_input_slabs else [0, 0, 1],  # Dummy value
        min_slab_thickness=min_slab_thickness if not use_input_slabs else 10.0,  # Dummy value
        min_vacuum_thickness=min_vacuum_thickness if not use_input_slabs else 10.0,  # Dummy value
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
        compute_relaxation_energy=compute_relaxation_energy,
        input_slabs=None,  # Always pass None to avoid serialization
        use_input_slabs=use_input_slabs,  # Pass the flag
        compute_cleavage=compute_cleavage,
        run_aimd=run_aimd,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        aimd_options=aimd_options,
        aimd_potential_mapping=aimd_potential_mapping,
        aimd_kpoints_spacing=aimd_kpoints_spacing,
        # Note: Electronic properties are handled manually below, not passed to core_workgraph
    )

    # If user provided slabs, manually add the relax_slabs_scatter task
    if use_input_slabs and relax_slabs:
        from teros.core.slabs import relax_slabs_scatter
        from teros.core.thermodynamics import identify_oxide_type, compute_surface_energies_scatter
        from aiida.orm import Code, load_code

        # Get code
        code = load_code(code_label)

        # Get parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        # Check if we have restart folders
        if restart_folders is not None:
            # RESTART MODE: Manually add VASP tasks with restart_folder
            print(f"  → Building workgraph with restart_folder for each slab")

            from aiida.plugins import WorkflowFactory
            from teros.core.slabs import extract_total_energy

            VaspWC = WorkflowFactory('vasp.v2.vasp')

            # Store task references for building output dictionaries
            vasp_tasks = {}
            energy_tasks = {}

            # Create VASP tasks for each slab with restart_folder
            for label, slab_struct in input_slabs.items():
                # Build VASP inputs
                vasp_inputs = {
                    'structure': slab_struct,
                    'code': code,
                    'parameters': {'incar': slab_params},
                    'options': slab_opts,
                    'potential_family': potential_family,
                    'potential_mapping': slab_pot_map,
                    'kpoints_spacing': slab_kpts,
                    'clean_workdir': clean_workdir,
                    'settings': orm.Dict(dict={'parser_settings': {'add_structure': True, 'add_trajectory': True}}),
                }

                # Add restart_folder if available for this slab
                if label in restart_folders:
                    vasp_inputs['restart_folder'] = restart_folders[label]
                    print(f"    → {label}: using restart_folder PK {restart_folders[label].pk}")

                # Add VASP task
                vasp_task = wg.add_task(
                    VaspWC,
                    name=f'VaspWorkChain_slab_{label}',
                    **vasp_inputs
                )
                vasp_tasks[label] = vasp_task

                # Add energy extraction task
                energy_task = wg.add_task(
                    extract_total_energy,
                    name=f'extract_energy_{label}',
                    energies=vasp_task.outputs.misc,
                )
                energy_tasks[label] = energy_task

            # Build output dictionaries that can be passed to thermodynamics/cleavage
            # These are socket dictionaries, not plain dictionaries
            relaxed_slabs_dict = {label: vasp_tasks[label].outputs.structure for label in input_slabs.keys()}
            slab_energies_dict = {label: energy_tasks[label].outputs.result for label in input_slabs.keys()}
            slab_remote_dict = {label: vasp_tasks[label].outputs.remote_folder for label in input_slabs.keys()}

            # Use a collector task to properly expose outputs in WorkGraph format
            from teros.core.slabs import collect_slab_outputs

            collector = wg.add_task(
                collect_slab_outputs,
                name='collect_slab_outputs_restart',
                structures=relaxed_slabs_dict,
                energies=slab_energies_dict,
                remote_folders=slab_remote_dict,
            )

            # Connect collector outputs to WorkGraph outputs
            wg.outputs.relaxed_slabs = collector.outputs.structures
            wg.outputs.slab_energies = collector.outputs.energies
            wg.outputs.slab_remote = collector.outputs.remote_folders
            wg.outputs.slab_structures = input_slabs

            print(f"  ✓ Created {len(vasp_tasks)} VASP tasks with restart capability")
            print(f"  ✓ All slab outputs connected via collector task")

        else:
            # NO RESTART: Use normal scf + relax + relaxation_energy calculation
            from teros.core.slabs import scf_slabs_scatter, calculate_relaxation_energies_scatter

            # Step 1: Add SCF task for unrelaxed slabs
            scf_task = wg.add_task(
                scf_slabs_scatter,
                name='scf_slabs_scatter',
                slabs=input_slabs,
                code=code,
                potential_family=potential_family,
                potential_mapping=slab_pot_map,
                parameters=slab_params,
                options=slab_opts,
                kpoints_spacing=slab_kpts,
                clean_workdir=clean_workdir,
            )

            # Step 2: Add relaxation task
            scatter_task_inputs = {
                'slabs': input_slabs,
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': slab_pot_map,
                'parameters': slab_params,
                'options': slab_opts,
                'kpoints_spacing': slab_kpts,
                'clean_workdir': clean_workdir,
            }

            scatter_task = wg.add_task(
                relax_slabs_scatter,
                name='relax_slabs_scatter',
                **scatter_task_inputs
            )

            # Step 3: Add relaxation energy calculation
            relax_energy_task = wg.add_task(
                calculate_relaxation_energies_scatter,
                name='calculate_relaxation_energies_scatter',
                unrelaxed_energies=scf_task.outputs.energies,
                relaxed_energies=scatter_task.outputs.energies,
            )

            # Connect slab outputs to graph outputs
            wg.outputs.relaxed_slabs = scatter_task.outputs.relaxed_structures
            wg.outputs.slab_energies = scatter_task.outputs.energies
            wg.outputs.slab_remote = scatter_task.outputs.remote_folders
            wg.outputs.unrelaxed_slab_energies = scf_task.outputs.energies
            wg.outputs.unrelaxed_slab_remote = scf_task.outputs.remote_folders
            wg.outputs.relaxation_energies = relax_energy_task.outputs.relaxation_energies
            wg.outputs.slab_structures = input_slabs

            # Use scatter task outputs for thermodynamics
            relaxed_slabs_dict = scatter_task.outputs.relaxed_structures
            slab_energies_dict = scatter_task.outputs.energies

        # Add thermodynamics calculation if requested
        if compute_thermodynamics:
            # Get tasks that were created in the core_workgraph
            bulk_vasp_task = wg.tasks['VaspWorkChain']  # Bulk relaxation task
            bulk_energy_task = wg.tasks['extract_total_energy']  # Bulk energy extraction
            formation_hf_task = wg.tasks['calculate_formation_enthalpy']

            # Add oxide type identification task with unique name
            oxide_type_task = wg.add_task(
                identify_oxide_type,
                name='identify_oxide_type_input_slabs',
                bulk_structure=bulk_vasp_task.outputs.structure,
            )

            # Use the output dictionaries created above (works for both restart and non-restart)
            # relaxed_slabs_dict and slab_energies_dict are defined in both branches

            # Add surface energies scatter task with unique name
            surface_energies_task = wg.add_task(
                compute_surface_energies_scatter,
                name='compute_surface_energies_input_slabs',
                slabs=relaxed_slabs_dict,
                energies=slab_energies_dict,
                bulk_structure=bulk_vasp_task.outputs.structure,
                bulk_energy=bulk_energy_task.outputs.result,
                reference_energies=formation_hf_task.outputs.result,
                formation_enthalpy=formation_hf_task.outputs.result,
                oxide_type=oxide_type_task.outputs.result,
                sampling=thermodynamics_sampling,
            )

            # Connect surface energies output
            wg.outputs.surface_energies = surface_energies_task.outputs.surface_energies
            print(f"  ✓ Thermodynamics calculation enabled (surface energies)")

        # Add cleavage energies calculation if requested
        # Note: Only add manually if we're in input_slabs mode (restart or manual slabs)
        # In normal mode, compute_cleavage is handled inside core_workgraph
        if compute_cleavage and use_input_slabs:
            from teros.core.cleavage import compute_cleavage_energies_scatter

            # Get tasks that were created in the core_workgraph
            bulk_vasp_task = wg.tasks['VaspWorkChain']  # Bulk relaxation task
            bulk_energy_task = wg.tasks['extract_total_energy']  # Bulk energy extraction

            # Use the output dictionaries created above (works for both restart and non-restart)
            # relaxed_slabs_dict and slab_energies_dict are defined in both branches

            # Add cleavage energies scatter task
            cleavage_task = wg.add_task(
                compute_cleavage_energies_scatter,
                name='compute_cleavage_energies_scatter_input_slabs',
                slabs=relaxed_slabs_dict,
                energies=slab_energies_dict,
                bulk_structure=bulk_vasp_task.outputs.structure,
                bulk_energy=bulk_energy_task.outputs.result,
            )

            # Connect cleavage energies output
            wg.outputs.cleavage_energies = cleavage_task.outputs.cleavage_energies
            print(f"  ✓ Cleavage energies calculation enabled")

    # NEW: Add electronic properties calculation if requested
    if compute_electronic_properties_bulk:
        from aiida.plugins import WorkflowFactory
        from aiida.orm import load_code

        # Get BandsWorkChain and wrap it as a task (same pattern as VaspWorkChain)
        BandsWorkChain = WorkflowFactory('vasp.v2.bands')
        BandsTask = task(BandsWorkChain)

        # Get bulk task from workgraph
        bulk_vasp_task = wg.tasks['VaspWorkChain']

        # Get code
        code = load_code(code_label)

        # Use bands-specific options or fall back to bulk options
        bands_opts = bands_options if bands_options is not None else bulk_options

        # Build inputs dictionary
        # BandsWorkChain expects inputs in namespaces: scf, bands, dos
        bands_inputs = {
            'structure': bulk_vasp_task.outputs.structure,  # Socket from bulk task
            'metadata': {
                'label': 'Bulk_Electronic_Properties',
                'description': 'DOS and bands calculation for relaxed bulk structure',
            }
        }

        # Add band_settings if provided
        if band_settings:
            bands_inputs['band_settings'] = band_settings

        # Build SCF namespace inputs (required)
        scf_inputs = {
            'code': code,
            'potential_family': potential_family,
            'potential_mapping': bulk_potential_mapping,
            'options': bands_opts,
            'clean_workdir': clean_workdir,
        }
        if bands_parameters and 'scf' in bands_parameters:
            scf_inputs['parameters'] = {'incar': bands_parameters['scf']}
        if bands_parameters and 'scf_kpoints_distance' in bands_parameters:
            scf_inputs['kpoints_spacing'] = bands_parameters['scf_kpoints_distance']
        
        bands_inputs['scf'] = scf_inputs

        # Build Bands namespace inputs (optional)
        if bands_parameters and 'bands' in bands_parameters:
            bands_inputs['bands'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': bulk_potential_mapping,
                'options': bands_opts,
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_parameters['bands']},
            }

        # Build DOS namespace inputs (optional)
        if bands_parameters and 'dos' in bands_parameters:
            bands_inputs['dos'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': bulk_potential_mapping,
                'options': bands_opts,
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_parameters['dos']},
            }

        # Add the bands task to the workgraph
        bands_task = wg.add_task(
            BandsTask,
            name='BandsWorkChain_bulk',
            **bands_inputs
        )

        # Connect all outputs from BandsWorkChain
        wg.outputs.bulk_bands = bands_task.outputs.band_structure
        wg.outputs.bulk_dos = bands_task.outputs.dos
        wg.outputs.bulk_primitive_structure = bands_task.outputs.primitive_structure
        wg.outputs.bulk_seekpath_parameters = bands_task.outputs.seekpath_parameters

        print(f"  ✓ Electronic properties calculation enabled (DOS and bands)")

    # Set the name
    wg.name = name

    return wg


def build_core_workgraph_with_map(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    oxygen_name: str,
    miller_indices: list = None,
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    metal_potential_mapping: dict = None,
    oxygen_potential_mapping: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    clean_workdir: bool = False,
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
    input_slabs: dict = None,
    compute_cleavage: bool = True,
    compute_electronic_properties_bulk: bool = False,  # NEW
    bands_parameters: dict = None,  # NEW
    bands_options: dict = None,  # NEW
    band_settings: dict = None,  # NEW
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
        oxygen_name=oxygen_name,
        nonmetal_name=nonmetal_name,
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
        input_slabs=input_slabs,
        compute_cleavage=compute_cleavage,
        compute_electronic_properties_bulk=compute_electronic_properties_bulk,  # NEW
        bands_parameters=bands_parameters,  # NEW
        bands_options=bands_options,  # NEW
        band_settings=band_settings,  # NEW
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
