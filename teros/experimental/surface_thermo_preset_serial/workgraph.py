"""Main workgraph for serial surface thermodynamics preset."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph
from ase.io import read

from teros.core.slabs import generate_slab_structures, extract_total_energy
from teros.core.hf import calculate_formation_enthalpy
from teros.core.thermodynamics import identify_oxide_type

from .slab_operations import (
    build_scf_slabs_nodes,
    build_relax_slabs_nodes,
    build_energy_extraction_nodes,
    build_relaxation_energy_nodes,
)
from .thermodynamics_operations import (
    build_surface_energy_nodes,
    select_surface_energy_by_oxide_type,
)
from .utils import (
    prepare_vasp_parameters,
    create_default_bulk_parameters,
    create_default_slab_parameters,
    create_default_scf_parameters,
)

VaspWorkChain = WorkflowFactory('vasp.vasp')


@task.graph(outputs=[
    'bulk_energy', 'bulk_structure',
    'metal_energy', 'oxygen_energy', 'nonmetal_energy',
    'metal_structure', 'oxygen_structure', 'nonmetal_structure',
    'formation_enthalpy', 'reference_energies',
    'slab_structures', 'relaxed_slabs',
    'slab_energies', 'unrelaxed_slab_energies', 'relaxation_energies',
    'surface_energies', 'oxide_type',
])
def surface_thermodynamics_serial_workgraph(
    structures_dir: str = None,
    bulk_name: str = None,
    code_label: str = None,
    potential_family: str = None,
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,
    # Reference materials
    metal_name: str = None,
    oxygen_name: str = None,
    nonmetal_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    # Slab parameters
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
    # Control flags
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,
    input_slabs: dict = None,
):
    """
    Serial surface thermodynamics workflow with flat graph structure.

    All VASP calculation nodes are added directly to the main graph,
    allowing max_concurrent_jobs to control concurrent execution.

    Args:
        structures_dir: Directory containing structure files
        bulk_name: Bulk structure filename
        code_label: VASP code label
        potential_family: Potential family name
        bulk_potential_mapping: Element to potential mapping for bulk
        kpoints_spacing: K-points spacing for bulk
        bulk_parameters: VASP parameters for bulk
        bulk_options: Computer options for bulk
        clean_workdir: Whether to clean working directory
        metal_name: Metal reference structure filename
        oxygen_name: Oxygen reference structure filename
        nonmetal_name: Optional nonmetal reference structure filename
        metal_potential_mapping: Element to potential mapping for metal
        metal_parameters: VASP parameters for metal
        metal_options: Computer options for metal
        oxygen_potential_mapping: Element to potential mapping for oxygen
        oxygen_parameters: VASP parameters for oxygen
        oxygen_options: Computer options for oxygen
        nonmetal_potential_mapping: Element to potential mapping for nonmetal
        nonmetal_parameters: VASP parameters for nonmetal
        nonmetal_options: Computer options for nonmetal
        miller_indices: List of Miller indices to generate slabs
        min_slab_thickness: Minimum slab thickness (Angstroms)
        min_vacuum_thickness: Minimum vacuum thickness (Angstroms)
        slab_parameters: VASP parameters for slabs
        slab_options: Computer options for slabs
        slab_potential_mapping: Element to potential mapping for slabs
        slab_kpoints_spacing: K-points spacing for slabs
        lll_reduce: Whether to LLL reduce slab
        center_slab: Whether to center slab in cell
        symmetrize: Whether to symmetrize slab
        primitive: Whether to use primitive slab
        in_unit_planes: Whether Miller indices are in unit planes
        max_normal_search: Maximum normal search range
        relax_slabs: Whether to relax slabs
        compute_thermodynamics: Whether to compute thermodynamics
        thermodynamics_sampling: Number of chemical potential sampling points
        compute_relaxation_energy: Whether to compute relaxation energy
        input_slabs: Optional pre-generated slabs dict

    Returns:
        Dictionary of outputs
    """
    # Create WorkGraph instance
    wg = WorkGraph(name="surface_thermodynamics_serial")

    # Load code
    code = orm.load_code(code_label)

    # Set default parameters
    if bulk_parameters is None:
        bulk_parameters = create_default_bulk_parameters()
    if slab_parameters is None:
        slab_parameters = create_default_slab_parameters()
    if slab_kpoints_spacing is None:
        slab_kpoints_spacing = kpoints_spacing

    # =========================================================================
    # PHASE 1: Bulk calculations
    # =========================================================================

    # Load bulk structure
    bulk_filepath = f"{structures_dir}/{bulk_name}"
    bulk_structure = orm.StructureData(ase=read(bulk_filepath))

    # Add bulk relaxation node
    bulk_vasp_params = prepare_vasp_parameters(
        base_parameters=bulk_parameters,
        code=code,
        potential_family=potential_family,
        potential_mapping=bulk_potential_mapping,
        kpoints_spacing=kpoints_spacing,
        options=bulk_options,
        clean_workdir=clean_workdir,
    )

    bulk_node = wg.add_task(
        VaspWorkChain,
        name="bulk_relax",
        structure=bulk_structure,
        **bulk_vasp_params,
    )

    # Extract bulk energy
    bulk_energy_node = wg.add_task(
        extract_total_energy,
        name="extract_bulk_energy",
        misc=bulk_node.outputs.misc,
    )

    # =========================================================================
    # PHASE 2: Reference calculations (for thermodynamics)
    # =========================================================================

    reference_nodes = {}
    reference_energy_nodes = {}

    if compute_thermodynamics:
        # Metal reference
        if metal_name:
            metal_filepath = f"{structures_dir}/{metal_name}"
            metal_structure = orm.StructureData(ase=read(metal_filepath))

            metal_vasp_params = prepare_vasp_parameters(
                base_parameters=metal_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=metal_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=metal_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['metal'] = wg.add_task(
                VaspWorkChain,
                name="metal_relax",
                structure=metal_structure,
                **metal_vasp_params,
            )

            reference_energy_nodes['metal'] = wg.add_task(
                extract_total_energy,
                name="extract_metal_energy",
                misc=reference_nodes['metal'].outputs.misc,
            )

        # Oxygen reference
        if oxygen_name:
            oxygen_filepath = f"{structures_dir}/{oxygen_name}"
            oxygen_structure = orm.StructureData(ase=read(oxygen_filepath))

            oxygen_vasp_params = prepare_vasp_parameters(
                base_parameters=oxygen_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=oxygen_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=oxygen_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['oxygen'] = wg.add_task(
                VaspWorkChain,
                name="oxygen_relax",
                structure=oxygen_structure,
                **oxygen_vasp_params,
            )

            reference_energy_nodes['oxygen'] = wg.add_task(
                extract_total_energy,
                name="extract_oxygen_energy",
                misc=reference_nodes['oxygen'].outputs.misc,
            )

        # Optional nonmetal reference
        if nonmetal_name:
            nonmetal_filepath = f"{structures_dir}/{nonmetal_name}"
            nonmetal_structure = orm.StructureData(ase=read(nonmetal_filepath))

            nonmetal_vasp_params = prepare_vasp_parameters(
                base_parameters=nonmetal_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=nonmetal_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=nonmetal_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['nonmetal'] = wg.add_task(
                VaspWorkChain,
                name="nonmetal_relax",
                structure=nonmetal_structure,
                **nonmetal_vasp_params,
            )

            reference_energy_nodes['nonmetal'] = wg.add_task(
                extract_total_energy,
                name="extract_nonmetal_energy",
                misc=reference_nodes['nonmetal'].outputs.misc,
            )

    # =========================================================================
    # Prepare reference energies dict
    # =========================================================================

    if compute_thermodynamics:
        # Build reference energies dict
        reference_energies_node = wg.add_task(
            build_reference_energies_dict,
            name="build_reference_energies",
            metal_energy=reference_energy_nodes.get('metal').outputs.result if 'metal' in reference_energy_nodes else None,
            oxygen_energy=reference_energy_nodes.get('oxygen').outputs.result if 'oxygen' in reference_energy_nodes else None,
            nonmetal_energy=reference_energy_nodes.get('nonmetal').outputs.result if 'nonmetal' in reference_energy_nodes else None,
            metal_structure=reference_nodes.get('metal').outputs.structure if 'metal' in reference_nodes else None,
            oxygen_structure=reference_nodes.get('oxygen').outputs.structure if 'oxygen' in reference_nodes else None,
            nonmetal_structure=reference_nodes.get('nonmetal').outputs.structure if 'nonmetal' in reference_nodes else None,
        )

        # Calculate formation enthalpy
        formation_enthalpy_node = wg.add_task(
            calculate_formation_enthalpy,
            name="formation_enthalpy",
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            reference_energies=reference_energies_node.outputs.result,
        )

        # Identify oxide type
        oxide_type_node = wg.add_task(
            identify_oxide_type,
            name="identify_oxide_type",
            bulk_structure=bulk_node.outputs.structure,
        )

    # =========================================================================
    # PHASE 3: Slab generation
    # =========================================================================

    if input_slabs:
        # Use provided slabs
        slab_structures = input_slabs
    else:
        # Generate slabs
        slab_gen_node = wg.add_task(
            generate_slab_structures,
            name="generate_slabs",
            bulk_structure=bulk_node.outputs.structure,
            miller_indices=orm.List(list=miller_indices or [[1,0,0], [1,1,0], [1,1,1]]),
            min_slab_thickness=orm.Float(min_slab_thickness),
            min_vacuum_thickness=orm.Float(min_vacuum_thickness),
            lll_reduce=orm.Bool(lll_reduce),
            center_slab=orm.Bool(center_slab),
            symmetrize=orm.Bool(symmetrize),
            primitive=orm.Bool(primitive),
            in_unit_planes=orm.Bool(in_unit_planes),
            max_normal_search=orm.Int(max_normal_search) if max_normal_search else None,
        )
        slab_structures = slab_gen_node.outputs.slabs

    # For direct node addition, we need to handle slabs as a known dict
    # In practice, we'll need to wait for slab generation or use input_slabs
    # For this implementation, let's assume input_slabs is provided as a dict

    if not input_slabs:
        raise NotImplementedError(
            "Dynamic slab generation not yet supported in serial mode. "
            "Please provide input_slabs as a dictionary."
        )

    # =========================================================================
    # PHASE 4: Slab calculations
    # =========================================================================

    # Build SCF (unrelaxed) nodes
    scf_params = create_default_scf_parameters()
    scf_nodes = build_scf_slabs_nodes(
        wg=wg,
        slabs=input_slabs,
        code=code,
        potential_family=potential_family,
        potential_mapping=slab_potential_mapping or bulk_potential_mapping,
        kpoints_spacing=slab_kpoints_spacing,
        parameters=scf_params,
        options=slab_options or bulk_options,
        clean_workdir=clean_workdir,
    )

    # Extract unrelaxed energies
    unrelaxed_energy_nodes = build_energy_extraction_nodes(
        wg=wg,
        vasp_nodes=scf_nodes,
        node_type="scf",
    )

    # Build relaxation nodes (if requested)
    if relax_slabs:
        relax_nodes = build_relax_slabs_nodes(
            wg=wg,
            slabs=input_slabs,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_potential_mapping or bulk_potential_mapping,
            kpoints_spacing=slab_kpoints_spacing,
            parameters=slab_parameters,
            options=slab_options or bulk_options,
            clean_workdir=clean_workdir,
        )

        # Extract relaxed energies
        relaxed_energy_nodes = build_energy_extraction_nodes(
            wg=wg,
            vasp_nodes=relax_nodes,
            node_type="relaxed",
        )

        # Calculate relaxation energies (if requested)
        if compute_relaxation_energy:
            relaxation_energy_nodes = build_relaxation_energy_nodes(
                wg=wg,
                unrelaxed_energies=unrelaxed_energy_nodes,
                relaxed_energies=relaxed_energy_nodes,
            )

    # =========================================================================
    # PHASE 5: Thermodynamics calculations
    # =========================================================================

    if compute_thermodynamics:
        # Use relaxed energies if available, otherwise SCF energies
        energy_nodes_for_thermo = relaxed_energy_nodes if relax_slabs else unrelaxed_energy_nodes

        # Build surface energy nodes
        surface_energy_nodes = build_surface_energy_nodes(
            wg=wg,
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            slab_structures=input_slabs,
            slab_energies=energy_nodes_for_thermo,
            reference_energies=reference_energies_node.outputs.result,
            formation_enthalpy=formation_enthalpy_node.outputs.result,
            oxide_type=oxide_type_node.outputs.result,
            sampling=thermodynamics_sampling,
        )

        # Select appropriate surface energy based on oxide type
        final_surface_energies = {}
        for slab_id, surf_eng_nodes in surface_energy_nodes.items():
            final_surface_energies[slab_id] = wg.add_task(
                select_surface_energy_by_oxide_type,
                name=f"select_surface_energy_{slab_id}",
                oxide_type=oxide_type_node.outputs.result,
                binary_result=surf_eng_nodes['binary'].outputs.result,
                ternary_result=surf_eng_nodes['ternary'].outputs.result,
            )

    # =========================================================================
    # Collect outputs
    # =========================================================================

    outputs = {
        'bulk_energy': bulk_energy_node.outputs.result,
        'bulk_structure': bulk_node.outputs.structure,
        'slab_structures': input_slabs,
    }

    if compute_thermodynamics:
        outputs['metal_energy'] = reference_energy_nodes.get('metal').outputs.result if 'metal' in reference_energy_nodes else None
        outputs['oxygen_energy'] = reference_energy_nodes.get('oxygen').outputs.result if 'oxygen' in reference_energy_nodes else None
        outputs['nonmetal_energy'] = reference_energy_nodes.get('nonmetal').outputs.result if 'nonmetal' in reference_energy_nodes else None
        outputs['metal_structure'] = reference_nodes.get('metal').outputs.structure if 'metal' in reference_nodes else None
        outputs['oxygen_structure'] = reference_nodes.get('oxygen').outputs.structure if 'oxygen' in reference_nodes else None
        outputs['nonmetal_structure'] = reference_nodes.get('nonmetal').outputs.structure if 'nonmetal' in reference_nodes else None
        outputs['formation_enthalpy'] = formation_enthalpy_node.outputs.result
        outputs['reference_energies'] = reference_energies_node.outputs.result
        outputs['oxide_type'] = oxide_type_node.outputs.result

        # Surface energies
        outputs['surface_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in final_surface_energies.items()
        }

    if relax_slabs:
        outputs['relaxed_slabs'] = {
            slab_id: node.outputs.structure
            for slab_id, node in relax_nodes.items()
        }
        outputs['slab_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in relaxed_energy_nodes.items()
        }

    outputs['unrelaxed_slab_energies'] = {
        slab_id: node.outputs.result
        for slab_id, node in unrelaxed_energy_nodes.items()
    }

    if compute_relaxation_energy and relax_slabs:
        outputs['relaxation_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in relaxation_energy_nodes.items()
        }

    return outputs


@task.calcfunction
def build_reference_energies_dict(
    metal_energy: orm.Float = None,
    oxygen_energy: orm.Float = None,
    nonmetal_energy: orm.Float = None,
    metal_structure: orm.StructureData = None,
    oxygen_structure: orm.StructureData = None,
    nonmetal_structure: orm.StructureData = None,
) -> orm.Dict:
    """
    Build reference energies dictionary for thermodynamics calculations.

    Args:
        metal_energy: Total energy of metal reference
        oxygen_energy: Total energy of oxygen reference
        nonmetal_energy: Total energy of nonmetal reference
        metal_structure: Metal reference structure
        oxygen_structure: Oxygen reference structure
        nonmetal_structure: Nonmetal reference structure

    Returns:
        Dictionary with '*_energy_per_atom' keys
    """
    ref_dict = {}

    if metal_energy and metal_structure:
        metal_ase = metal_structure.get_ase()
        num_metal_atoms = len(metal_ase)
        ref_dict['metal_energy_per_atom'] = metal_energy.value / num_metal_atoms

    if oxygen_energy and oxygen_structure:
        oxygen_ase = oxygen_structure.get_ase()
        num_oxygen_atoms = len(oxygen_ase)
        ref_dict['oxygen_energy_per_atom'] = oxygen_energy.value / num_oxygen_atoms

    if nonmetal_energy and nonmetal_structure:
        nonmetal_ase = nonmetal_structure.get_ase()
        num_nonmetal_atoms = len(nonmetal_ase)
        ref_dict['nonmetal_energy_per_atom'] = nonmetal_energy.value / num_nonmetal_atoms

    return orm.Dict(dict=ref_dict)
