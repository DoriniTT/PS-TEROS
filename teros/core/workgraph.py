"""
PS-TEROS Core WorkGraph

This module contains the core workflow for PS-TEROS calculations using
the pythonic scatter-gather pattern for parallel execution.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, dynamic, namespace, group
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
from teros.core.adsorption_energy import compute_adsorption_energies_scatter
from teros.core.workflow_presets import (
    resolve_preset,
    validate_preset_inputs,
    validate_flag_dependencies,
    check_old_style_api,
)

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


@task.calcfunction()
def create_dummy_float(value: float = 0.0) -> orm.Float:
    """Create a dummy Float node for workflows that don't need certain outputs."""
    return orm.Float(value)


@task.calcfunction()
def create_dummy_structure() -> orm.StructureData:
    """Create a dummy StructureData node for workflows that don't need certain outputs."""
    from ase import Atoms
    dummy = Atoms()
    return orm.StructureData(ase=dummy)

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
    # Slab electronic properties outputs
    'slab_bands', 'slab_dos', 'slab_primitive_structures', 'slab_seekpath_parameters',
    # Adsorption energy outputs
    'relaxed_complete_structures', 'separated_structures', 'substrate_energies', 'molecule_energies',
    'complete_energies', 'adsorption_energies',
])
def core_workgraph(
    structures_dir: str = None,
    bulk_name: str = None,
    code_label: str = None,
    bulk_code_label: str = None,
    metal_code_label: str = None,
    nonmetal_code_label: str = None,
    oxygen_code_label: str = None,
    slab_code_label: str = None,
    potential_family: str = None,
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
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    slab_relax_builder_inputs: dict = None,
    slab_scf_builder_inputs: dict = None,
    structure_specific_relax_builder_inputs: dict = None,
    structure_specific_scf_builder_inputs: dict = None,
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
    compute_electronic_properties_slabs: bool = False,  # NEW
    slab_electronic_properties: dict = None,  # NEW
    slab_bands_parameters: dict = None,  # NEW
    slab_bands_options: dict = None,  # NEW
    slab_band_settings: dict = None,  # NEW
    run_adsorption_energy: bool = False,  # NEW
    adsorption_structures: t.Annotated[dict, dynamic(orm.StructureData)] = None,  # NEW - Dynamic namespace of StructureData
    adsorption_formulas: t.Annotated[dict, dict] = None,  # NEW - Dict of {label: formula_string}
    adsorption_parameters: dict = None,  # NEW
    adsorption_options: dict = None,  # NEW
    adsorption_potential_mapping: dict = None,  # NEW
    adsorption_kpoints_spacing: float = None,  # NEW
    # Adsorption energy: Relaxation control (NEW)
    relax_before_adsorption: bool = False,
    adsorption_relax_builder_inputs: dict = None,
    adsorption_scf_builder_inputs: dict = None,
    # Adsorption energy: Atom fixing control (NEW)
    adsorption_fix_atoms: bool = False,
    adsorption_fix_type: str = None,
    adsorption_fix_thickness: float = 0.0,
    adsorption_fix_elements: list = None,
    # Adsorption energy: Structure-specific builder inputs (NEW)
    adsorption_structure_specific_relax_builder_inputs: dict = None,
    adsorption_structure_specific_scf_builder_inputs: dict = None,
    adsorption_structure_component_specific_scf_builder_inputs: dict = None,
    # Concurrency control
    max_concurrent_jobs: int = None,
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
    # Load the codes
    # Main code (fallback for all)
    code = orm.load_code(code_label)

    # Reference codes (use specific codes if provided, otherwise fall back to main code)
    bulk_code = orm.load_code(bulk_code_label) if bulk_code_label else code
    metal_code = orm.load_code(metal_code_label) if metal_code_label else code
    nonmetal_code = orm.load_code(nonmetal_code_label) if nonmetal_code_label else code
    oxygen_code = orm.load_code(oxygen_code_label) if oxygen_code_label else code
    slab_code = orm.load_code(slab_code_label) if slab_code_label else code

    # Get VASP workchain and wrap it as a task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Determine if we should compute formation enthalpy
    # (requires both metal_name and oxygen_name)
    compute_formation_enthalpy = (metal_name is not None and oxygen_name is not None)

    # ===== BULK RELAXATION (CONDITIONAL) =====
    # Only load and relax bulk if needed
    if structures_dir and bulk_name:
        bulk_struct = load_structure_from_file(f"{structures_dir}/{bulk_name}")

        bulk_vasp = VaspTask(
            structure=bulk_struct,
            code=bulk_code,
            parameters={'incar': bulk_parameters},
            options=bulk_options,
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=bulk_potential_mapping,
            clean_workdir=clean_workdir,
            settings=orm.Dict(dict=get_settings()),
        )

        bulk_energy = extract_total_energy(energies=bulk_vasp.misc)
    else:
        # No bulk structure provided - set to None
        bulk_vasp = None
        bulk_energy = None

    # ===== REFERENCE RELAXATIONS AND FORMATION ENTHALPY (CONDITIONAL) =====
    if compute_formation_enthalpy:
        # ===== METAL RELAXATION =====
        metal_struct = load_structure_from_file(f"{structures_dir}/{metal_name}")

        metal_vasp = VaspTask(
            structure=metal_struct,
            code=metal_code,
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
                code=nonmetal_code,
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
            code=oxygen_code,
            parameters={'incar': oxygen_parameters},
            options=oxygen_options,
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=oxygen_potential_mapping,
            clean_workdir=clean_workdir,
            settings=orm.Dict(dict=get_settings()),
        )

        oxygen_energy = extract_total_energy(energies=oxygen_vasp.misc)

        # ===== TASK DEPENDENCIES: Serial reference calculations =====
        # Chain: Bulk >> Metal >> Nonmetal >> Oxygen
        bulk_vasp >> metal_vasp
        if nonmetal_name is not None:
            metal_vasp >> nonmetal_vasp >> oxygen_vasp
        else:
            metal_vasp >> oxygen_vasp

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
        # Bulk-only mode or no bulk: Create dummy outputs with proper types
        if structures_dir and bulk_name:
            # We have bulk - use bulk structure as placeholder for reference structures
            metal_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
            metal_energy = type('obj', (object,), {'result': bulk_energy.result})()
            nonmetal_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
            nonmetal_energy = type('obj', (object,), {'result': bulk_energy.result})()
            oxygen_vasp = type('obj', (object,), {'structure': bulk_vasp.structure})()
            oxygen_energy = type('obj', (object,), {'result': bulk_energy.result})()
            formation_hf = type('obj', (object,), {'result': bulk_energy.result})()  # Placeholder
        else:
            # No bulk at all - set everything to None
            metal_vasp = None
            metal_energy = None
            nonmetal_vasp = None
            nonmetal_energy = None
            oxygen_vasp = None
            oxygen_energy = None
            formation_hf = None

    # ===== ELECTRONIC PROPERTIES CALCULATION (OPTIONAL) =====
    # Note: Electronic properties are added in build_core_workgraph, not here
    # This keeps core_workgraph focused on base calculations

    # ===== SLAB GENERATION =====
    # Use input slabs if provided, otherwise generate them
    if use_input_slabs:
        # User will provide pre-generated slabs after building
        # Skip slab generation entirely
        slab_namespace = None  # Will be set post-build
    elif miller_indices is None or not (structures_dir and bulk_name):
        # No miller_indices provided, no input_slabs, or no bulk - skip slab generation
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

    # ===== TASK DEPENDENCIES: Create reference group for stage separation =====
    # Group all reference calculations to chain to slab stage
    if compute_formation_enthalpy:
        if nonmetal_name is not None:
            references = group(bulk_vasp, metal_vasp, nonmetal_vasp, oxygen_vasp)
        else:
            references = group(bulk_vasp, metal_vasp, oxygen_vasp)
    else:
        # Bulk-only mode
        references = bulk_vasp if bulk_vasp is not None else None

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
                code=slab_code,  # Use slab_code
                potential_family=potential_family,
                potential_mapping=slab_pot_map,
                parameters=slab_params,
                options=slab_opts,
                kpoints_spacing=slab_kpts,
                clean_workdir=clean_workdir,
                max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
                scf_builder_inputs=slab_scf_builder_inputs,  # NEW
                structure_specific_scf_builder_inputs=structure_specific_scf_builder_inputs,  # NEW
            )
            unrelaxed_slab_energies_output = scf_outputs.energies
            unrelaxed_slab_remote_output = scf_outputs.remote_folders

            # ===== TASK DEPENDENCIES: References >> SCF =====
            if references is not None:
                references >> scf_outputs

        # ALWAYS: Relax slabs
        relaxation_outputs = relax_slabs_scatter(
            slabs=slab_namespace,
            code=slab_code,  # Use slab_code
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            parameters=slab_params,
            options=slab_opts,
            kpoints_spacing=slab_kpts,
            clean_workdir=clean_workdir,
            max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
            relax_builder_inputs=slab_relax_builder_inputs,  # NEW
            structure_specific_relax_builder_inputs=structure_specific_relax_builder_inputs,  # NEW
        )

        # ===== TASK DEPENDENCIES: SCF >> Relaxation or References >> Relaxation =====
        if compute_relaxation_energy:
            # Chain: SCF >> Relaxation
            scf_outputs >> relaxation_outputs
        else:
            # Chain: References >> Relaxation (no SCF stage)
            if references is not None:
                references >> relaxation_outputs

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
    if compute_thermodynamics and relax_slabs and slab_namespace is not None and (structures_dir and bulk_name):
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
    if compute_cleavage and relax_slabs and not use_input_slabs and (structures_dir and bulk_name):
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
    # Note: AIMD tasks are added manually in build_core_workgraph using wg.add_task()
    # Sequential dependencies cannot be created inside @task.graph
    aimd_outputs = {}
    # ===== SLAB ELECTRONIC PROPERTIES CALCULATION (OPTIONAL) =====
    # Initialize output dictionaries for slab electronic properties
    slab_bands_output = {}
    slab_dos_output = {}
    slab_primitive_structures_output = {}
    slab_seekpath_parameters_output = {}
    
    # Add electronic properties for auto-generated slabs (same pattern as cleavage)
    if compute_electronic_properties_slabs and relax_slabs and slab_electronic_properties and not use_input_slabs:
        from teros.core.slabs import calculate_electronic_properties_slabs_scatter
        from aiida.orm import load_code
        
        # Get code
        code = load_code(code_label)
        
        # Get default parameters
        default_params = slab_bands_parameters if slab_bands_parameters else {}
        default_opts = slab_bands_options if slab_bands_options else bulk_options
        default_settings = slab_band_settings if slab_band_settings else {}
        
        # Get slab parameters
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        
        # Calculate electronic properties for selected slabs
        slab_elec_outputs = calculate_electronic_properties_slabs_scatter(
            slabs=relaxed_slabs_output,
            slab_electronic_properties=slab_electronic_properties,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            clean_workdir=clean_workdir,
            default_bands_parameters=default_params,
            default_bands_options=default_opts,
            default_band_settings=default_settings,
            max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
        )
        
        # Collect outputs
        slab_bands_output = slab_elec_outputs.slab_bands
        slab_dos_output = slab_elec_outputs.slab_dos
        slab_primitive_structures_output = slab_elec_outputs.slab_primitive_structures
        slab_seekpath_parameters_output = slab_elec_outputs.slab_seekpath_parameters

    # Return all outputs
    # Note: Electronic properties outputs (bulk_bands, bulk_dos, bulk_electronic_properties_misc)
    # are added dynamically in build_core_workgraph, not returned here
    # For workflows without bulk, return dummy nodes to satisfy @task.graph output types
    dummy_float = create_dummy_float(0.0)
    dummy_struct = create_dummy_structure()
    
    # Determine if we have bulk based on input parameters (not future values)
    has_bulk = structures_dir and bulk_name
    
    return {
        'bulk_energy': bulk_energy.result if has_bulk else dummy_float,
        'bulk_structure': bulk_vasp.structure if has_bulk else dummy_struct,
        'metal_energy': metal_energy.result if compute_formation_enthalpy else dummy_float,
        'metal_structure': metal_vasp.structure if compute_formation_enthalpy else dummy_struct,
        'nonmetal_energy': nonmetal_energy.result if compute_formation_enthalpy else dummy_float,
        'nonmetal_structure': nonmetal_vasp.structure if compute_formation_enthalpy else dummy_struct,
        'oxygen_energy': oxygen_energy.result if compute_formation_enthalpy else dummy_float,
        'oxygen_structure': oxygen_vasp.structure if compute_formation_enthalpy else dummy_struct,
        'formation_enthalpy': formation_hf.result if compute_formation_enthalpy else dummy_float,
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
        'slab_bands': slab_bands_output,
        'slab_dos': slab_dos_output,
        'slab_primitive_structures': slab_primitive_structures_output,
        'slab_seekpath_parameters': slab_seekpath_parameters_output,
    }


def build_core_workgraph(
    structures_dir: str = None,
    bulk_name: str = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    calculator: str = 'vasp',  # NEW: 'vasp' or 'cp2k'
    aimd_code_label: str = None,  # NEW: Optional separate code for AIMD
    bulk_code_label: str = None,  # NEW: Optional separate code for bulk calculations
    metal_code_label: str = None,  # NEW: Optional separate code for metal reference
    nonmetal_code_label: str = None,  # NEW: Optional separate code for nonmetal reference
    oxygen_code_label: str = None,  # NEW: Optional separate code for oxygen reference
    slab_code_label: str = None,  # NEW: Optional separate code for slab calculations
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
    # NEW: Slab builder inputs (full control over VASP builder)
    slab_relax_builder_inputs: dict = None,  # Builder inputs for slab relaxation
    slab_scf_builder_inputs: dict = None,  # Builder inputs for slab SCF (unrelaxed)
    structure_specific_relax_builder_inputs: dict = None,  # Per-slab overrides for relaxation
    structure_specific_scf_builder_inputs: dict = None,  # Per-slab overrides for SCF
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    workflow_preset: str = None,  # NEW: Workflow preset name
    relax_slabs: bool = None,  # CHANGED: Now defaults to None (preset controls)
    compute_thermodynamics: bool = None,  # CHANGED: Now defaults to None
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = None,  # CHANGED: Now defaults to None
    input_slabs: dict = None,
    compute_cleavage: bool = None,  # CHANGED: Now defaults to None
    restart_from_node: int = None,  # PK of previous workgraph to restart from
    compute_electronic_properties_bulk: bool = None,  # CHANGED: Now defaults to None
    bands_parameters: dict = None,  # NEW: INCAR parameters
    bands_options: dict = None,  # NEW: Scheduler options
    band_settings: dict = None,  # NEW: Band workflow settings
    run_aimd: bool = None,  # CHANGED: Now defaults to None
    aimd_sequence: list = None,
    aimd_parameters: dict = None,
    aimd_options: dict = None,
    aimd_potential_mapping: dict = None,
    aimd_kpoints_spacing: float = None,
    # CP2K-specific parameters (NEW):
    basis_content: str = None,
    pseudo_content: str = None,
    cp2k_kind_section: list = None,
    # Fixed atoms configuration (NEW):
    fix_atoms: bool = False,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list = None,
    fix_components: str = "XYZ",
    compute_electronic_properties_slabs: bool = None,  # CHANGED: Now defaults to None
    slab_electronic_properties: dict = None,  # NEW: Per-slab parameter overrides
    slab_bands_parameters: dict = None,  # NEW: Default slab parameters
    slab_bands_options: dict = None,  # NEW: Default slab scheduler options
    slab_band_settings: dict = None,  # NEW: Default slab band settings
    run_adsorption_energy: bool = None,  # CHANGED: Now defaults to None
    adsorption_structures: dict = None,  # Dict of {label: StructureData} with substrate+adsorbate
    adsorption_formulas: dict = None,  # Dict of {label: formula_string} for adsorbates (e.g., {"site1": "OH"})
    adsorption_parameters: dict = None,  # VASP parameters for adsorption calculations
    adsorption_options: dict = None,  # Scheduler options for adsorption calculations
    adsorption_potential_mapping: dict = None,  # Element to potential mapping
    adsorption_kpoints_spacing: float = None,  # K-points spacing for adsorption
    # Adsorption energy: Relaxation control (NEW)
    relax_before_adsorption: bool = False,
    adsorption_relax_builder_inputs: dict = None,
    adsorption_scf_builder_inputs: dict = None,
    # Adsorption energy: Atom fixing control (NEW)
    adsorption_fix_atoms: bool = False,
    adsorption_fix_type: str = None,
    adsorption_fix_thickness: float = 0.0,
    adsorption_fix_elements: list = None,
    # Adsorption energy: Structure-specific builder inputs (NEW)
    adsorption_structure_specific_relax_builder_inputs: dict = None,
    adsorption_structure_specific_scf_builder_inputs: dict = None,
    adsorption_structure_component_specific_scf_builder_inputs: dict = None,
    max_concurrent_jobs: int = 4,  # NEW: Limit concurrent VASP calculations (None = unlimited)
    name: str = 'FormationEnthalpy',
):
    """
    Build a centralized WorkGraph for bulk relaxation, formation enthalpy, and surface calculations.

    This function builds a flexible workflow that can handle:
    - Bulk structure relaxation (always performed)
    - Reference structure relaxation (if metal_name and oxygen_name provided)
    - Formation enthalpy calculation (if references provided)
    - Slab generation and relaxation (if miller_indices or input_slabs provided)
    - Relaxation energy calculation (controlled by workflow preset or compute_relaxation_energy flag)
    - Cleavage energy calculation (controlled by workflow preset or compute_cleavage flag)
    - Surface thermodynamics (controlled by workflow preset or compute_thermodynamics flag, requires references)
    - Electronic properties calculation (controlled by workflow preset or compute_electronic_properties_bulk flag)
    - AIMD simulations (controlled by workflow preset or run_aimd flag)

    **NEW: Workflow Preset System**
    
    Instead of manually setting 7+ boolean flags, you can use workflow presets for common use cases:
    
    - 'surface_thermodynamics' (default): Complete surface thermodynamics workflow
    - 'surface_thermodynamics_unrelaxed': Quick unrelaxed surface energy screening
    - 'cleavage_only': Cleavage energy calculations only
    - 'relaxation_energy_only': Relaxation energy calculations only
    - 'bulk_only': Bulk relaxation only (no surfaces)
    - 'formation_enthalpy_only': Formation enthalpy without surfaces
    - 'electronic_structure_bulk_only': Electronic properties (DOS/bands) for bulk
    - 'aimd_only': AIMD simulation on slabs
    - 'comprehensive': Complete analysis (all features enabled)
    
    Use list_workflow_presets() to see all available presets with descriptions.
    
    Presets can be overridden by explicitly setting individual flags:
        wg = build_core_workgraph(
            workflow_preset='surface_thermodynamics',
            compute_cleavage=False,  # Override preset default
            ...
        )

    Args:
        workflow_preset: Name of workflow preset to use. If None, defaults to 'surface_thermodynamics'.
                        Use list_workflow_presets() to see available presets.
        workflow_preset: Name of workflow preset to use. If None, defaults to 'surface_thermodynamics'.
                        Use list_workflow_presets() to see available presets.
        structures_dir: Directory containing structure files
        bulk_name: Filename of bulk structure (e.g., 'ag3po4.cif')
        code_label: VASP code label in AiiDA. Default: 'VASP-VTST-6.4.3@bohr'
                   Used for slab calculations and as fallback for references.
        bulk_code_label: Optional separate VASP code for bulk calculations.
                        Default: None (uses code_label)
        metal_code_label: Optional separate VASP code for metal reference.
                         Default: None (uses code_label)
        nonmetal_code_label: Optional separate VASP code for nonmetal reference.
                            Default: None (uses code_label)
        oxygen_code_label: Optional separate VASP code for oxygen reference.
                          Default: None (uses code_label)
        slab_code_label: Optional separate VASP code for slab calculations.
                        Default: None (uses code_label)
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
        slab_relax_builder_inputs: Dict of builder parameters for slab relaxation (full control).
                                  Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                                  If provided, overrides slab_parameters, slab_options, slab_kpoints_spacing.
                                  Default: None
        slab_scf_builder_inputs: Dict of builder parameters for slab SCF (unrelaxed energy, full control).
                                Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                                If provided, overrides slab_parameters, slab_options, slab_kpoints_spacing.
                                Default: None

        relax_slabs: Relax generated slabs. Default: None (controlled by preset)
        compute_thermodynamics: Compute surface energies (requires metal_name and oxygen_name). Default: None (controlled by preset)
        thermodynamics_sampling: Grid resolution for chemical potential sampling. Default: 100
        compute_relaxation_energy: Compute relaxation energies (E_relaxed - E_unrelaxed). Default: None (controlled by preset)
        compute_cleavage: Compute cleavage energies. Default: None (controlled by preset)
        input_slabs: Pre-generated slab structures dict. Skips slab generation if provided. Default: None
        restart_from_node: PK of previous workgraph to restart from (extracts slabs and RemoteData). Default: None
        compute_electronic_properties_bulk: Compute DOS and bands for relaxed bulk. Default: None (controlled by preset)
        bands_parameters: Dict with 'scf', 'bands', 'dos' keys containing INCAR dicts, plus
                          optional 'scf_kpoints_distance'. Use get_electronic_properties_defaults().
                          Default: None
        bands_options: Scheduler options for bands calculation. Default: None (uses bulk_options)
        band_settings: Dict with band workflow settings. Use get_electronic_properties_defaults()
                       ['band_settings']. Default: None
        run_aimd: Run AIMD on generated/input slabs. Default: None (controlled by preset)
        aimd_sequence: List of AIMD stages [{'temperature': K, 'steps': N}, ...]. Default: None
        aimd_parameters: AIMD INCAR parameters (use get_aimd_defaults()). Default: None
        aimd_options: Scheduler options for AIMD. Default: None (uses slab_options)
        aimd_potential_mapping: Potential mapping for AIMD. Default: None (uses slab_potential_mapping)
        aimd_kpoints_spacing: K-points spacing for AIMD. Default: None (uses slab_kpoints_spacing)
        compute_electronic_properties_slabs: Compute DOS and bands for selected slabs. Default: None (controlled by preset)
        slab_electronic_properties: Dict mapping slab labels to parameter overrides.
                                    Format: {'term_0': {'bands_parameters': ..., 'bands_options': ..., 'band_settings': ...}}
                                    If a parameter key is missing, falls back to slab_bands_* defaults.
                                    Default: None
        slab_bands_parameters: Default parameters for all slabs. Use get_slab_electronic_properties_defaults().
                               Default: None
        slab_bands_options: Default scheduler options for slab electronic properties. Default: None (uses bulk_options)
        slab_band_settings: Default band settings for slabs. Use get_slab_electronic_properties_defaults()
                           ['band_settings']. Default: None

        basis_content: CP2K basis set content string. Required for CP2K calculations. Default: None
        pseudo_content: CP2K pseudopotential content string. Required for CP2K calculations. Default: None
        cp2k_kind_section: List of CP2K KIND section parameters. Default: None

        fix_atoms: Enable fixed atom constraints. Default: False
        fix_type: Type of constraint ('bottom', 'top', 'elements', 'manual'). Default: None
        fix_thickness: Thickness in Angstroms for 'bottom'/'top' fixing. Default: 0.0
        fix_elements: List of element symbols to fix (for 'elements' type). Default: None
        fix_components: Cartesian components to constrain ('X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ'). Default: 'XYZ'

        run_adsorption_energy: Enable adsorption energy calculations. Default: None (controlled by preset)
        adsorption_structures: Dict mapping labels to StructureData objects containing substrate+adsorbate systems.
                              Format: {'site1': structure1, 'site2': structure2, ...}
                              Each structure should contain both the substrate (e.g., slab) and adsorbate.
                              Default: None
        adsorption_formulas: Dict mapping labels to adsorbate chemical formulas.
                            Format: {'site1': 'OH', 'site2': 'OOH', ...}
                            The formula is used to identify and separate the adsorbate from the substrate.
                            Must match the actual adsorbate composition in adsorption_structures.
                            Default: None
        adsorption_parameters: VASP INCAR parameters for adsorption energy calculations.
                              Applied to all three calculations (complete, substrate, molecule).
                              Use get_adsorption_defaults() for recommended settings.
                              Default: None
        adsorption_options: Scheduler options for adsorption calculations (walltime, nodes, etc.).
                           Applied to all three calculations (complete, substrate, molecule).
                           Default: None
        adsorption_potential_mapping: Element to potential mapping for adsorption calculations.
                                     Format: {'O': 'O', 'H': 'H', ...}
                                     Should include all elements in substrate and adsorbate.
                                     Default: None (uses bulk_potential_mapping)
        adsorption_kpoints_spacing: K-points spacing in A^-1 * 2pi for adsorption calculations.
                                   Applies to all three calculations (complete, substrate, molecule).
                                   Default: None (uses kpoints_spacing)
        relax_before_adsorption: If True, relax the complete substrate+adsorbate structure before SCF calculation.
                                When enabled, the adsorption energy workflow becomes:
                                1. Relax complete system (if relax_before_adsorption=True)
                                2. SCF on relaxed/unrelaxed complete system
                                3. Separate and relax substrate
                                4. Separate and relax molecule
                                5. Calculate E_ads = E_complete - E_substrate - E_molecule
                                Default: False
        adsorption_relax_builder_inputs: Dict of builder parameters for relaxing the complete substrate+adsorbate.
                                        Format: {'parameters': {...}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                                        Only used if relax_before_adsorption=True.
                                        If not provided, uses adsorption_parameters/options/kpoints_spacing.
                                        Default: None
        adsorption_scf_builder_inputs: Dict of builder parameters for SCF calculation on complete system.
                                      Format: {'parameters': {...}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
                                      If not provided, uses adsorption_parameters/options/kpoints_spacing.
                                      Default: None

        name: WorkGraph name. Default: 'FormationEnthalpy'

    Returns:
        WorkGraph: Instance ready to be submitted. The workflow provides outputs including:
            - relaxed_bulk: Relaxed bulk structure (StructureData)
            - bulk_energy: Energy of relaxed bulk (Float)
            - formation_enthalpy: Formation enthalpy if references provided (Float)
            - slabs: Generated/input slab structures (dict)
            - relaxed_slabs: Relaxed slab structures if relax_slabs=True (dict)
            - surface_energies: Surface energies if compute_thermodynamics=True (dict)
            - cleavage_energies: Cleavage energies if compute_cleavage=True (dict)
            - relaxation_energies: Slab relaxation energies if compute_relaxation_energy=True (dict)
            - dos, bands: Electronic properties if compute_electronic_properties_bulk=True
            - relaxed_complete_structures: Relaxed substrate+adsorbate structures if run_adsorption_energy=True
                                          and relax_before_adsorption=True (dict). Only populated when relaxation
                                          is enabled for adsorption calculations.
            - adsorption_energies: Adsorption energies if run_adsorption_energy=True (dict)
            - substrate_energies: Substrate energies from adsorption calculations (dict)
            - molecule_energies: Molecule energies from adsorption calculations (dict)
            - complete_energies: Complete system energies from adsorption calculations (dict)

    Examples:
        Using workflow preset (recommended):
        >>> wg = build_core_workgraph(
        ...     workflow_preset='surface_thermodynamics',
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     metal_name='Ag.cif',
        ...     oxygen_name='O2.cif',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     bulk_parameters={...},
        ...     bulk_options={...},
        ... )
        
        Overriding preset defaults:
        >>> wg = build_core_workgraph(
        ...     workflow_preset='surface_thermodynamics',
        ...     compute_cleavage=False,  # Override
        ...     ...
        ... )
        
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

        Adsorption energy with relaxation (recommended for accurate energies):
        >>> from aiida.orm import load_node
        >>> # Load pre-generated substrate+adsorbate structures
        >>> oh_structure = load_node(12345).outputs.structure  # OH on Ag(111)
        >>> ooh_structure = load_node(12346).outputs.structure  # OOH on Ag(111)
        >>>
        >>> wg = build_core_workgraph(
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6, ...},
        ...     bulk_options={'resources': {'num_machines': 2}},
        ...     # Enable adsorption energy calculations
        ...     run_adsorption_energy=True,
        ...     adsorption_structures={
        ...         'oh_site': oh_structure,
        ...         'ooh_site': ooh_structure,
        ...     },
        ...     adsorption_formulas={
        ...         'oh_site': 'OH',
        ...         'ooh_site': 'OOH',
        ...     },
        ...     # Relax substrate+adsorbate before calculating adsorption energy
        ...     relax_before_adsorption=True,
        ...     adsorption_relax_builder_inputs={
        ...         'parameters': {'IBRION': 2, 'NSW': 100, 'EDIFF': 1e-6, 'EDIFFG': -0.02},
        ...         'options': {'resources': {'num_machines': 2}, 'walltime': '24:00:00'},
        ...         'kpoints_spacing': 0.2,
        ...     },
        ...     # SCF parameters for relaxed structure
        ...     adsorption_scf_builder_inputs={
        ...         'parameters': {'NSW': 0, 'EDIFF': 1e-7, 'PREC': 'Accurate'},
        ...         'options': {'resources': {'num_machines': 1}, 'walltime': '6:00:00'},
        ...         'kpoints_spacing': 0.15,
        ...     },
        ...     # Shared settings for substrate and molecule calculations
        ...     adsorption_potential_mapping={'Ag': 'Ag', 'O': 'O', 'H': 'H'},
        ... )
        >>> # Access relaxed structures
        >>> relaxed_oh = wg.outputs.relaxed_complete_structures['oh_site']
        >>> relaxed_ooh = wg.outputs.relaxed_complete_structures['ooh_site']
        >>> # Access adsorption energies
        >>> e_ads_oh = wg.outputs.adsorption_energies['oh_site']
        >>> e_ads_ooh = wg.outputs.adsorption_energies['ooh_site']

        Adsorption energy without relaxation (faster, less accurate):
        >>> wg = build_core_workgraph(
        ...     structures_dir='structures',
        ...     bulk_name='ag3po4.cif',
        ...     bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     bulk_parameters={...},
        ...     bulk_options={...},
        ...     # Enable adsorption energy without relaxation
        ...     run_adsorption_energy=True,
        ...     adsorption_structures={'oh_site': oh_structure},
        ...     adsorption_formulas={'oh_site': 'OH'},
        ...     relax_before_adsorption=False,  # Skip relaxation (default)
        ...     # Only SCF parameters needed
        ...     adsorption_parameters={'NSW': 0, 'EDIFF': 1e-7, ...},
        ...     adsorption_options={'resources': {'num_machines': 1}},
        ...     adsorption_kpoints_spacing=0.15,
        ... )
    """
    # ========================================================================
    # WORKFLOW PRESET RESOLUTION
    # ========================================================================

    # Check for deprecated old-style API usage
    check_old_style_api(
        workflow_preset,
        relax_slabs,
        compute_thermodynamics,
        compute_cleavage,
        compute_relaxation_energy,
        compute_electronic_properties_bulk,
        compute_electronic_properties_slabs,
        run_aimd,
        run_adsorption_energy,
    )

    # Resolve preset and apply user overrides
    resolved_preset_name, resolved_flags = resolve_preset(
        workflow_preset,
        relax_slabs,
        compute_thermodynamics,
        compute_cleavage,
        compute_relaxation_energy,
        compute_electronic_properties_bulk,
        compute_electronic_properties_slabs,
        run_aimd,
        run_adsorption_energy,
    )
    
    # Extract resolved flags
    relax_slabs = resolved_flags['relax_slabs']
    compute_thermodynamics = resolved_flags['compute_thermodynamics']
    compute_cleavage = resolved_flags['compute_cleavage']
    compute_relaxation_energy = resolved_flags['compute_relaxation_energy']
    compute_electronic_properties_bulk = resolved_flags['compute_electronic_properties_bulk']
    compute_electronic_properties_slabs = resolved_flags['compute_electronic_properties_slabs']
    run_aimd = resolved_flags['run_aimd']
    run_adsorption_energy = resolved_flags['run_adsorption_energy']
    
    # Validate preset requirements
    validation_errors = validate_preset_inputs(
        resolved_preset_name,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        nonmetal_name=nonmetal_name,
        miller_indices=miller_indices,
        bands_parameters=bands_parameters,
        bands_options=bands_options,
        band_settings=band_settings,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        aimd_options=aimd_options,
        aimd_potential_mapping=aimd_potential_mapping,
        aimd_kpoints_spacing=aimd_kpoints_spacing,
        slab_bands_parameters=slab_bands_parameters,
        slab_bands_options=slab_bands_options,
        slab_band_settings=slab_band_settings,
        slab_electronic_properties=slab_electronic_properties,
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
    )
    
    if validation_errors:
        error_msg = "\n".join(validation_errors)
        raise ValueError(f"Preset validation failed:\n{error_msg}")
    
    # Check flag dependencies and emit warnings
    dependency_warnings = validate_flag_dependencies(
        resolved_flags,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        miller_indices=miller_indices,
        input_slabs=input_slabs,
        bands_parameters=bands_parameters,
        slab_bands_parameters=slab_bands_parameters,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
    )
    
    for warning in dependency_warnings:
        if warning.startswith("ERROR:"):
            raise ValueError(warning)
        print(f"\n⚠️  {warning}\n")
    
    # ========================================================================
    # VALIDATE REQUIRED INPUTS
    # ========================================================================
    
    # Determine if bulk structure is needed
    compute_formation_enthalpy = (metal_name is not None and oxygen_name is not None)
    needs_bulk = (
        compute_formation_enthalpy or 
        compute_thermodynamics or 
        compute_electronic_properties_bulk or
        (miller_indices is not None and input_slabs is None)  # Need bulk for slab generation
    )
    
    # Only require structures_dir and bulk_name if we actually need the bulk structure
    if needs_bulk:
        if structures_dir is None or bulk_name is None:
            raise ValueError(
                "structures_dir and bulk_name are required for workflows that need bulk structure.\n"
                f"Current preset '{resolved_preset_name}' requires bulk structure for:\n"
                + (f"  - Formation enthalpy calculation\n" if compute_formation_enthalpy else "")
                + (f"  - Surface thermodynamics\n" if compute_thermodynamics else "")
                + (f"  - Electronic properties (bulk)\n" if compute_electronic_properties_bulk else "")
                + (f"  - Slab generation (miller_indices provided)\n" if (miller_indices and not input_slabs) else "")
                + "\nEither provide structures_dir and bulk_name, or use input_slabs with aimd_only preset."
            )
    
    # Print workflow configuration
    print(f"\n{'='*70}")
    print(f"WORKFLOW CONFIGURATION")
    print(f"{'='*70}")
    print(f"Preset: {resolved_preset_name}")
    print(f"\nActive Components:")
    for flag_name, flag_value in resolved_flags.items():
        status = "✓" if flag_value else "✗"
        print(f"  {status} {flag_name}: {flag_value}")
    print(f"{'='*70}\n")
    
    # ========================================================================
    # RESTART HANDLING
    # ========================================================================
    
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

    # ========================================================================
    # CP2K-SPECIFIC SETUP
    # ========================================================================

    basis_file = None
    pseudo_file = None
    if calculator == 'cp2k':
        import io
        from teros.core.builders.aimd_builder_cp2k import (
            get_basis_molopt_content,
            get_gth_potentials_content,
        )

        # Use provided content or defaults
        basis_str = basis_content if basis_content else get_basis_molopt_content()
        pseudo_str = pseudo_content if pseudo_content else get_gth_potentials_content()

        # Create SinglefileData nodes
        basis_file = orm.SinglefileData(
            io.BytesIO(basis_str.encode('utf-8')),
            filename='BASIS_MOLOPT'
        )
        pseudo_file = orm.SinglefileData(
            io.BytesIO(pseudo_str.encode('utf-8')),
            filename='GTH_POTENTIALS'
        )

        print(f"\n  ✓ Created CP2K basis and pseudopotential files")

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
        bulk_code_label=bulk_code_label,
        metal_code_label=metal_code_label,
        nonmetal_code_label=nonmetal_code_label,
        oxygen_code_label=oxygen_code_label,
        slab_code_label=slab_code_label,
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
        slab_relax_builder_inputs=slab_relax_builder_inputs,
        slab_scf_builder_inputs=slab_scf_builder_inputs,
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
        compute_electronic_properties_slabs=compute_electronic_properties_slabs,  # NEW
        slab_electronic_properties=slab_electronic_properties,  # NEW
        slab_bands_parameters=slab_bands_parameters,  # NEW
        slab_bands_options=slab_bands_options,  # NEW
        slab_band_settings=slab_band_settings,  # NEW
        run_adsorption_energy=run_adsorption_energy,  # NEW
        adsorption_structures=adsorption_structures,  # NEW
        adsorption_formulas=adsorption_formulas,  # NEW
        adsorption_parameters=adsorption_parameters,  # NEW
        adsorption_options=adsorption_options,  # NEW
        adsorption_potential_mapping=adsorption_potential_mapping,  # NEW
        adsorption_kpoints_spacing=adsorption_kpoints_spacing,  # NEW
        relax_before_adsorption=relax_before_adsorption,  # NEW
        adsorption_relax_builder_inputs=adsorption_relax_builder_inputs,  # NEW
        adsorption_scf_builder_inputs=adsorption_scf_builder_inputs,  # NEW
        adsorption_fix_atoms=adsorption_fix_atoms,  # NEW
        adsorption_fix_type=adsorption_fix_type,  # NEW
        adsorption_fix_thickness=adsorption_fix_thickness,  # NEW
        adsorption_fix_elements=adsorption_fix_elements,  # NEW
        adsorption_structure_specific_relax_builder_inputs=adsorption_structure_specific_relax_builder_inputs,  # NEW
        adsorption_structure_specific_scf_builder_inputs=adsorption_structure_specific_scf_builder_inputs,  # NEW
        adsorption_structure_component_specific_scf_builder_inputs=adsorption_structure_component_specific_scf_builder_inputs,  # NEW
        max_concurrent_jobs=max_concurrent_jobs,  # NEW
        # Note: Electronic properties are handled manually below, not passed to core_workgraph
    )

    # If user provided slabs, manually add the relax_slabs_scatter task
    if use_input_slabs and relax_slabs:
        from teros.core.slabs import relax_slabs_scatter
        from teros.core.thermodynamics import identify_oxide_type, compute_surface_energies_scatter
        from aiida.orm import Code, load_code

        # Get code (use slab_code_label if provided, otherwise code_label)
        code = load_code(slab_code_label if slab_code_label else code_label)

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
                code=slab_code,  # Use slab_code
                potential_family=potential_family,
                potential_mapping=slab_pot_map,
                parameters=slab_params,
                options=slab_opts,
                kpoints_spacing=slab_kpts,
                clean_workdir=clean_workdir,
                max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
                scf_builder_inputs=slab_scf_builder_inputs,  # NEW
                structure_specific_scf_builder_inputs=structure_specific_scf_builder_inputs,  # NEW
            )

            # Step 2: Add relaxation task
            scatter_task_inputs = {
                'slabs': input_slabs,
                'code': slab_code,  # Use slab_code
                'potential_family': potential_family,
                'potential_mapping': slab_pot_map,
                'parameters': slab_params,
                'options': slab_opts,
                'kpoints_spacing': slab_kpts,
                'clean_workdir': clean_workdir,
                'max_number_jobs': orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
                'relax_builder_inputs': slab_relax_builder_inputs,  # NEW
                'structure_specific_relax_builder_inputs': structure_specific_relax_builder_inputs,  # NEW
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

    # NEW: Add slab electronic properties calculation if requested
    # Currently only supported for input_slabs mode (including restart)
    if compute_electronic_properties_slabs and relax_slabs and slab_electronic_properties and use_input_slabs:
        from teros.core.slabs import calculate_electronic_properties_slabs_scatter

        # Get default parameters (per-slab overrides handled in scatter function)
        default_params = slab_bands_parameters if slab_bands_parameters else {}
        default_opts = slab_bands_options if slab_bands_options else bulk_options
        default_settings = slab_band_settings if slab_band_settings else {}

        # Get code
        code = load_code(code_label)

        # Get slab parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping

        # Determine which slabs output to use
        if restart_folders is not None:
            # Restart mode: use collector outputs
            relaxed_slabs_source = collector.outputs.structures
        else:
            # Non-restart input_slabs: use scatter task outputs
            relaxed_slabs_source = scatter_task.outputs.relaxed_structures

        # Add electronic properties task
        slab_elec_task = wg.add_task(
            calculate_electronic_properties_slabs_scatter,
            name='calculate_electronic_properties_slabs',
            slabs=relaxed_slabs_source,
            slab_electronic_properties=slab_electronic_properties,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            clean_workdir=clean_workdir,
            default_bands_parameters=default_params,
            default_bands_options=default_opts,
            default_band_settings=default_settings,
            max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
        )

        # Connect outputs
        wg.outputs.slab_bands = slab_elec_task.outputs.slab_bands
        wg.outputs.slab_dos = slab_elec_task.outputs.slab_dos
        wg.outputs.slab_primitive_structures = slab_elec_task.outputs.slab_primitive_structures
        wg.outputs.slab_seekpath_parameters = slab_elec_task.outputs.slab_seekpath_parameters

        print(f"  ✓ Slab electronic properties calculation enabled for {len(slab_electronic_properties)} slabs")

    # Add electronic properties calculation if requested
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

    # ===== AIMD CALCULATION (OPTIONAL) =====
    # Check if AIMD should be added
    should_add_aimd = run_aimd and aimd_sequence is not None and (input_slabs is not None or miller_indices is not None)

    if should_add_aimd:
        print("\n  → Adding AIMD stages to workflow")
        print(f"     Calculator: {calculator}")
        print(f"     Number of stages: {len(aimd_sequence)}")

        # Import appropriate scatter function based on calculator
        if calculator == 'vasp':
            from teros.core.aimd_functions import aimd_single_stage_scatter
            aimd_scatter_func = aimd_single_stage_scatter
        elif calculator == 'cp2k':
            from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k
            aimd_scatter_func = aimd_single_stage_scatter_cp2k
        else:
            raise ValueError(f"Unknown calculator: {calculator}")

        from aiida.orm import Code, load_code

        # Load code - use aimd_code_label if provided, otherwise use code_label
        aimd_code_to_use = aimd_code_label if aimd_code_label is not None else code_label
        code = load_code(aimd_code_to_use)
        print(f"     Using code: {aimd_code_to_use}")

        # Get parameters - use AIMD-specific or fall back to slab/bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        aimd_params = aimd_parameters if aimd_parameters is not None else slab_params
        aimd_opts = aimd_options if aimd_options is not None else slab_opts
        aimd_pot_map = aimd_potential_mapping if aimd_potential_mapping is not None else slab_pot_map
        aimd_kpts = aimd_kpoints_spacing if aimd_kpoints_spacing is not None else slab_kpts

        # Handle fixed atoms if requested
        fixed_atoms_lists = {}

        if fix_atoms and fix_type is not None:
            print(f"\n  → Preparing fixed atoms constraints")
            print(f"     Type: {fix_type}")
            print(f"     Thickness: {fix_thickness} Å")
            print(f"     Elements: {fix_elements if fix_elements else 'all'}")

            from teros.core.fixed_atoms import get_fixed_atoms_list

            # For input_slabs, calculate fixed atoms now
            if input_slabs is not None:
                for label, slab_struct in input_slabs.items():
                    fixed_list = get_fixed_atoms_list(
                        slab_struct,
                        fix_type=fix_type,
                        fix_thickness=fix_thickness,
                        fix_elements=fix_elements,
                    )
                    fixed_atoms_lists[label] = fixed_list
                    print(f"     {label}: {len(fixed_list)} atoms fixed")

        # Determine which slabs to use as initial structures
        if relax_slabs and 'relax_slabs_scatter' in wg.tasks:
            relax_task = wg.tasks['relax_slabs_scatter']
            initial_slabs_source = relax_task.outputs.relaxed_structures
            print(f"     Using relaxed slabs from relax_slabs_scatter task")
        elif relax_slabs and 'collect_slab_outputs_restart' in wg.tasks:
            collector_task = wg.tasks['collect_slab_outputs_restart']
            initial_slabs_source = collector_task.outputs.structures
            print(f"     Using relaxed slabs from restart collector")
        elif input_slabs is not None:
            initial_slabs_source = input_slabs
            print(f"     Using unrelaxed input slabs as initial structures")
        else:
            gen_task = wg.tasks['generate_slab_structures']
            initial_slabs_source = gen_task.outputs.slabs
            print(f"     Using unrelaxed generated slabs as initial structures")

        # Sequential AIMD stages
        current_structures = initial_slabs_source
        current_remotes = {}
        stage_tasks = []

        for stage_idx, stage_config in enumerate(aimd_sequence):
            # Get temperature for stage naming (support both old and new format temporarily)
            stage_temp = stage_config.get('TEBEG', stage_config.get('temperature', 0))
            stage_name = f"aimd_stage_{stage_idx:02d}_{stage_temp}K"

            # Print stage info (support both old and new format temporarily)
            stage_steps = stage_config.get('NSW', stage_config.get('steps', 0))
            print(f"     Stage {stage_idx}: {stage_temp}K × {stage_steps} steps")

            # Build stage inputs based on calculator
            if calculator == 'vasp':
                stage_inputs = {
                    'slabs': current_structures,
                    'stage_config': stage_config,
                    'code': code,
                    'base_aimd_parameters': aimd_params,
                    'structure_aimd_overrides': None,
                    'potential_family': potential_family,
                    'potential_mapping': aimd_pot_map,
                    'options': aimd_opts,
                    'kpoints_spacing': aimd_kpts,
                    'clean_workdir': clean_workdir,
                    'restart_folders': current_remotes,
                    'max_number_jobs': orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
                }
            elif calculator == 'cp2k':
                stage_inputs = {
                    'slabs': current_structures,
                    'temperature': stage_config['temperature'],
                    'steps': stage_config['steps'],
                    'code': code,
                    'base_aimd_parameters': aimd_params,
                    'structure_aimd_overrides': None,
                    'basis_file': basis_file,
                    'pseudo_file': pseudo_file,
                    'options': aimd_opts,  # Pass scheduler options directly
                    'clean_workdir': clean_workdir,
                    'restart_folders': current_remotes,
                    'max_number_jobs': orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
                }

                # Add fixed atoms for CP2K
                if fixed_atoms_lists:
                    # Pre-computed fixed atoms (for input_slabs)
                    stage_inputs['fixed_atoms_lists'] = fixed_atoms_lists
                    stage_inputs['fix_components'] = fix_components
                elif fix_atoms and fix_type is not None:
                    # Dynamic fixed atoms calculation (for auto-generated slabs)
                    stage_inputs['fix_type'] = fix_type
                    stage_inputs['fix_thickness'] = fix_thickness if fix_thickness else 0.0
                    stage_inputs['fix_elements'] = fix_elements
                    stage_inputs['fix_components'] = fix_components

            # Add task
            stage_task = wg.add_task(
                aimd_scatter_func,
                name=stage_name,
                **stage_inputs
            )

            stage_tasks.append(stage_task)

            # Wire outputs to next stage inputs
            current_structures = stage_task.outputs.structures
            current_remotes = stage_task.outputs.remote_folders

        print(f"  ✓ AIMD calculation enabled ({len(aimd_sequence)} sequential stages)")
        print(f"     Access AIMD outputs via: wg.tasks['aimd_stage_XX_XXXK'].outputs")

    # ===== ADSORPTION ENERGY CALCULATION (OPTIONAL) =====
    # Check if adsorption energy should be computed
    should_add_adsorption = run_adsorption_energy and adsorption_structures is not None and adsorption_formulas is not None

    if should_add_adsorption:
        print("\n  → Adding adsorption energy calculation")
        print(f"     Number of structures: {len(adsorption_structures)}")
        print(f"     Structure keys: {list(adsorption_structures.keys())}")
        print(f"     Adsorbate formulas: {adsorption_formulas}")
        if relax_before_adsorption:
            print(f"     Relaxation: ENABLED (relax → separate → SCF)")
        else:
            print(f"     Relaxation: DISABLED (separate → SCF)")

        from aiida.orm import Code, load_code

        # Load code
        code = load_code(code_label)

        # Get parameters - use adsorption-specific or fall back to slab/bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        ads_params = adsorption_parameters if adsorption_parameters is not None else slab_params
        ads_opts = adsorption_options if adsorption_options is not None else slab_opts
        ads_pot_map = adsorption_potential_mapping if adsorption_potential_mapping is not None else slab_pot_map
        ads_kpts = adsorption_kpoints_spacing if adsorption_kpoints_spacing is not None else slab_kpts

        # Extract relax and SCF parameters from builder_inputs or use old-style
        # For backward compatibility, extract INCAR params if only builder_inputs provided
        relax_params = None
        scf_params = None

        if adsorption_relax_builder_inputs is not None:
            # New style: extract INCAR from builder_inputs (for backward compat)
            relax_params = adsorption_relax_builder_inputs.get('parameters', {})
            if 'incar' in relax_params:
                relax_params = relax_params['incar']
        elif ads_params is not None:
            # Old style: use ads_params directly
            relax_params = dict(ads_params)

        if adsorption_scf_builder_inputs is not None:
            # New style: extract INCAR from builder_inputs (for backward compat)
            scf_params = adsorption_scf_builder_inputs.get('parameters', {})
            if 'incar' in scf_params:
                scf_params = scf_params['incar']
        elif ads_params is not None:
            # Old style: use ads_params directly
            scf_params = dict(ads_params)

        # Add adsorption energy scatter task (with NEW builder inputs support)
        adsorption_task = wg.add_task(
            compute_adsorption_energies_scatter,
            name='compute_adsorption_energies_scatter',
            structures=adsorption_structures,
            adsorbate_formulas=adsorption_formulas,
            code=code,
            potential_family=potential_family,
            potential_mapping=ads_pot_map,

            # Relaxation parameters (old-style for backward compat)
            relax_before_adsorption=relax_before_adsorption,
            relax_parameters=relax_params,
            scf_parameters=scf_params,

            # Common settings (old-style for backward compat)
            options=ads_opts,
            kpoints_spacing=ads_kpts,
            clean_workdir=clean_workdir,

            # NEW: Builder inputs (full control)
            relax_builder_inputs=adsorption_relax_builder_inputs,
            scf_builder_inputs=adsorption_scf_builder_inputs,
            structure_specific_relax_builder_inputs=adsorption_structure_specific_relax_builder_inputs,
            structure_specific_scf_builder_inputs=adsorption_structure_specific_scf_builder_inputs,
            structure_component_specific_scf_builder_inputs=adsorption_structure_component_specific_scf_builder_inputs,

            # Atom fixing parameters
            fix_atoms=adsorption_fix_atoms,
            fix_type=adsorption_fix_type,
            fix_thickness=adsorption_fix_thickness,
            fix_elements=adsorption_fix_elements,

            # Concurrency control
            max_number_jobs=orm.Int(max_concurrent_jobs) if max_concurrent_jobs is not None else None,
        )

        # Connect outputs
        wg.outputs.relaxed_complete_structures = adsorption_task.outputs.relaxed_complete_structures  # NEW
        wg.outputs.separated_structures = adsorption_task.outputs.separated_structures
        wg.outputs.substrate_energies = adsorption_task.outputs.substrate_energies
        wg.outputs.molecule_energies = adsorption_task.outputs.molecule_energies
        wg.outputs.complete_energies = adsorption_task.outputs.complete_energies
        wg.outputs.adsorption_energies = adsorption_task.outputs.adsorption_energies

        print(f"  ✓ Adsorption energy calculation enabled")
        if relax_before_adsorption:
            print(f"     Access relaxed structures via: wg.outputs.relaxed_complete_structures")
        print(f"     Access adsorption energies via: wg.outputs.adsorption_energies")

    # Set the name
    wg.name = name

    # CONCURRENCY CONTROL: Limit how many VASP calculations run simultaneously
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg

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
