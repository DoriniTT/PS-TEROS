"""Main workgraph for serial surface thermodynamics preset."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, namespace, dynamic
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


@task.calcfunction
def generate_multiple_miller_slabs(
    bulk_structure: orm.StructureData,
    miller_indices_list: orm.List,
    min_slab_thickness: orm.Float,
    min_vacuum_thickness: orm.Float,
    lll_reduce: orm.Bool,
    center_slab: orm.Bool,
    symmetrize: orm.Bool,
    primitive: orm.Bool,
    in_unit_planes: orm.Bool,
    max_normal_search: orm.Int,
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """
    Generate slab structures for multiple Miller indices from a bulk structure.

    This function generates slabs for each Miller index in the provided list,
    naming them as slab_{h}{k}{l}_term_{n}.

    Args:
        bulk_structure: Relaxed bulk structure
        miller_indices_list: List of Miller indices, e.g., [[1,0,0], [1,1,0]]
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        lll_reduce: Whether to reduce with LLL algorithm
        center_slab: Whether to center slab in cell
        symmetrize: Whether to generate symmetric terminations
        primitive: Whether to use primitive cell
        in_unit_planes: Whether thickness in unit planes
        max_normal_search: Max search for surface normal

    Returns:
        Dict with 'slabs' key containing all generated slabs
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    adaptor = AseAtomsAdaptor()
    ase_structure = bulk_structure.get_ase()
    pymatgen_structure = adaptor.get_structure(ase_structure)

    if primitive.value:
        analyzer = SpacegroupAnalyzer(pymatgen_structure)
        pymatgen_structure = analyzer.get_primitive_standard_structure()

    all_slabs = {}

    for miller in miller_indices_list.get_list():
        miller_str = ''.join(map(str, miller))

        generator = SlabGenerator(
            pymatgen_structure,
            miller,
            min_slab_thickness.value,
            min_vacuum_thickness.value,
            center_slab=center_slab.value,
            in_unit_planes=in_unit_planes.value if in_unit_planes else False,
            max_normal_search=max_normal_search.value if max_normal_search else None,
            lll_reduce=lll_reduce.value,
        )

        slabs = generator.get_slabs(symmetrize=symmetrize.value)

        for index, slab in enumerate(slabs):
            orthogonal_slab = slab.get_orthogonal_c_slab()
            ase_slab = adaptor.get_atoms(orthogonal_slab)
            slab_id = f"slab_{miller_str}_term_{index}"
            all_slabs[slab_id] = orm.StructureData(ase=ase_slab)

    return {'slabs': all_slabs}


def generate_slabs_from_relaxed_bulk_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,
    bulk_builder: dict = None,
    # Slab generation parameters
    miller_indices: list = None,
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
):
    """
    Preliminary workflow: Relax bulk structure and generate slabs.

    This is Stage 1 of the two-stage workflow. It:
    1. Relaxes the bulk structure
    2. Generates slabs from the relaxed bulk structure
    3. Returns the generated slabs as outputs

    The returned slabs can then be used as input_slabs for the main
    surface_thermodynamics_serial_workgraph.

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
        bulk_builder: Full builder dictionary for bulk (optional)
        miller_indices: List of Miller indices to generate slabs
        min_slab_thickness: Minimum slab thickness (Angstroms)
        min_vacuum_thickness: Minimum vacuum thickness (Angstroms)
        lll_reduce: Whether to LLL reduce slab
        center_slab: Whether to center slab in cell
        symmetrize: Whether to symmetrize slab
        primitive: Whether to use primitive slab
        in_unit_planes: Whether Miller indices are in unit planes
        max_normal_search: Maximum normal search range

    Returns:
        WorkGraph that outputs the generated slabs
    """
    from pathlib import Path
    from ase.io import read

    wg = WorkGraph("generate_slabs_from_relaxed_bulk")

    # Load bulk structure
    structure_path = Path(structures_dir) / bulk_name
    bulk_ase = read(str(structure_path))
    bulk_structure = orm.StructureData(ase=bulk_ase)

    # Phase 1: Bulk relaxation
    if bulk_builder:
        bulk_node = wg.add_task(
            VaspWorkChain,
            name="bulk_relax",
            structure=bulk_structure,
            **bulk_builder,
        )
    else:
        # Load code
        code = orm.load_code(code_label)

        # Prepare parameters
        bulk_params = create_default_bulk_parameters()
        bulk_vasp_params = prepare_vasp_parameters(
            base_parameters=bulk_parameters or bulk_params,
            code=code,
            potential_family=potential_family,
            potential_mapping=bulk_potential_mapping or {},
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

    # Phase 2: Generate slabs from relaxed bulk
    slab_gen_node = wg.add_task(
        generate_multiple_miller_slabs,
        name="generate_slabs",
        bulk_structure=bulk_node.outputs.structure,  # Relaxed bulk
        miller_indices_list=orm.List(miller_indices),
        min_slab_thickness=orm.Float(min_slab_thickness),
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        symmetrize=orm.Bool(symmetrize),
        primitive=orm.Bool(primitive),
        in_unit_planes=orm.Bool(in_unit_planes),
        max_normal_search=orm.Int(max_normal_search) if max_normal_search else orm.Int(1),
    )

    # Set workgraph outputs
    wg.outputs = {
        'slabs': slab_gen_node.outputs.slabs,
        'relaxed_bulk': bulk_node.outputs.structure,
    }

    return wg


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
    # Builder dictionaries for full control (optional)
    bulk_builder: dict = None,
    metal_builder: dict = None,
    oxygen_builder: dict = None,
    nonmetal_builder: dict = None,
    slab_scf_builders: dict = None,  # Dict of {slab_id: builder}
    slab_relax_builders: dict = None,  # Dict of {slab_id: builder}
):
    """
    Serial surface thermodynamics workflow with flat graph structure.

    All VASP calculation nodes are added directly to the main graph,
    allowing max_concurrent_jobs to control concurrent execution.

    When decorated with @task.graph, the workgraph is created automatically by the decorator.
    This function should define the workflow logic and return outputs.

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
        input_slabs: REQUIRED. Pre-generated slabs as dict of slab_id -> StructureData.
                     Dynamic generation not yet supported in serial mode.

    Returns:
        Dictionary of outputs including:
        - slab_structures: The input_slabs dict (passed through)
        - bulk_energy, bulk_structure: Bulk calculation results
        - relaxed_slabs, slab_energies: Slab calculation results (if relax_slabs=True)
        - unrelaxed_slab_energies: SCF energies for all slabs
        - surface_energies: Thermodynamic surface energies (if compute_thermodynamics=True)
        - formation_enthalpy, oxide_type, reference_energies: Thermodynamic properties
    """
    # Create WorkGraph instance
    wg = WorkGraph(name="surface_thermodynamics_serial")

    # Load code with validation
    try:
        code = orm.load_code(code_label)
    except Exception as e:
        raise ValueError(
            f"Code '{code_label}' not found or not accessible. "
            f"Use 'verdi code list' to see available codes. "
            f"Error: {str(e)}"
        )

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
    if bulk_builder:
        # Use provided builder for full control
        bulk_node = wg.add_task(
            VaspWorkChain,
            name="bulk_relax",
            structure=bulk_structure,
            **bulk_builder,
        )
    else:
        # Use standard parameter preparation
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
        energies=bulk_node.outputs.misc,
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

            if metal_builder:
                # Use provided builder
                reference_nodes['metal'] = wg.add_task(
                    VaspWorkChain,
                    name="metal_relax",
                    structure=metal_structure,
                    **metal_builder,
                )
            else:
                # Use standard parameter preparation
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
                energies=reference_nodes['metal'].outputs.misc,
            )

        # Oxygen reference
        if oxygen_name:
            oxygen_filepath = f"{structures_dir}/{oxygen_name}"
            oxygen_structure = orm.StructureData(ase=read(oxygen_filepath))

            if oxygen_builder:
                # Use provided builder
                reference_nodes['oxygen'] = wg.add_task(
                    VaspWorkChain,
                    name="oxygen_relax",
                    structure=oxygen_structure,
                    **oxygen_builder,
                )
            else:
                # Use standard parameter preparation
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
                energies=reference_nodes['oxygen'].outputs.misc,
            )

        # Optional nonmetal reference
        if nonmetal_name:
            nonmetal_filepath = f"{structures_dir}/{nonmetal_name}"
            nonmetal_structure = orm.StructureData(ase=read(nonmetal_filepath))

            if nonmetal_builder:
                # Use provided builder
                reference_nodes['nonmetal'] = wg.add_task(
                    VaspWorkChain,
                    name="nonmetal_relax",
                    structure=nonmetal_structure,
                    **nonmetal_builder,
                )
            else:
                # Use standard parameter preparation
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
                energies=reference_nodes['nonmetal'].outputs.misc,
            )

    # =========================================================================
    # Calculate formation enthalpy
    # =========================================================================

    if compute_thermodynamics:
        # For binary oxides, use metal as dummy for nonmetal if not provided
        nonmetal_structure_out = (
            reference_nodes['nonmetal'].outputs.structure if 'nonmetal' in reference_nodes
            else reference_nodes['metal'].outputs.structure
        )
        nonmetal_energy_out = (
            reference_energy_nodes['nonmetal'].outputs.result if 'nonmetal' in reference_energy_nodes
            else reference_energy_nodes['metal'].outputs.result
        )

        formation_enthalpy_node = wg.add_task(
            calculate_formation_enthalpy,
            name="formation_enthalpy",
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            metal_structure=reference_nodes['metal'].outputs.structure,
            metal_energy=reference_energy_nodes['metal'].outputs.result,
            oxygen_structure=reference_nodes['oxygen'].outputs.structure,
            oxygen_energy=reference_energy_nodes['oxygen'].outputs.result,
            nonmetal_structure=nonmetal_structure_out,
            nonmetal_energy=nonmetal_energy_out,
        )

        # Build reference energies dictionary
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

        # Identify oxide type
        oxide_type_node = wg.add_task(
            identify_oxide_type,
            name="identify_oxide_type",
            bulk_structure=bulk_node.outputs.structure,
        )

    # =========================================================================
    # PHASE 3: Slab generation from RELAXED bulk
    # =========================================================================

    # Three modes: pre-provided slabs, generate from relaxed bulk, or generate from input bulk
    slab_gen_node = None
    expected_slab_ids = []

    if input_slabs:
        # Mode 1: Use pre-provided slabs (no slab generation needed)
        expected_slab_ids = list(input_slabs.keys())

    elif miller_indices:
        # Mode 2: Generate slabs from RELAXED bulk structure
        # This waits for bulk relaxation to complete, then generates slabs at runtime

        slab_gen_node = wg.add_task(
            generate_multiple_miller_slabs,
            name="generate_slabs_from_relaxed_bulk",
            bulk_structure=bulk_node.outputs.structure,  # Relaxed bulk!
            miller_indices_list=orm.List(miller_indices),
            min_slab_thickness=orm.Float(min_slab_thickness),
            min_vacuum_thickness=orm.Float(min_vacuum_thickness),
            lll_reduce=orm.Bool(lll_reduce),
            center_slab=orm.Bool(center_slab),
            symmetrize=orm.Bool(symmetrize),
            primitive=orm.Bool(primitive),
            in_unit_planes=orm.Bool(in_unit_planes),
            max_normal_search=orm.Int(max_normal_search) if max_normal_search else orm.Int(1),
        )

        # Pre-allocate expected slab IDs based on miller_indices
        # Assume max 20 terminations per Miller index (conservative estimate)
        # Can be increased if needed for complex structures
        max_terminations_per_miller = 20
        for miller in miller_indices:
            miller_str = ''.join(map(str, miller))
            for term_idx in range(max_terminations_per_miller):
                expected_slab_ids.append(f"slab_{miller_str}_term_{term_idx}")
    else:
        raise ValueError(
            "Either input_slabs or miller_indices must be provided. "
            "Provide input_slabs as a dict of {slab_id: StructureData}, "
            "OR provide miller_indices as a list of tuples like [(1,0,0), (1,1,0)]"
        )

    # =========================================================================
    # PHASE 4: Slab calculations (SCF and Relaxation)
    # =========================================================================

    scf_nodes = {}
    relax_nodes = {}
    unrelaxed_energy_nodes = {}
    relaxed_energy_nodes = {}

    # Create VASP tasks for each expected slab
    scf_params = create_default_scf_parameters()

    for slab_id in expected_slab_ids:
        # Determine slab structure source
        if input_slabs:
            # Mode 1: Use pre-provided slab
            slab_structure = input_slabs[slab_id]
        else:
            # Mode 2: Use dynamically generated slab from calcfunction output
            # The output namespace is: slab_gen_node.outputs.slabs.<slab_id>
            slab_structure = getattr(slab_gen_node.outputs.slabs, slab_id)

        # SCF calculation
        if slab_scf_builders and slab_id in slab_scf_builders:
            scf_nodes[slab_id] = wg.add_task(
                VaspWorkChain,
                name=f"scf_slab_{slab_id}",
                structure=slab_structure,
                **slab_scf_builders[slab_id],
            )
        else:
            scf_vasp_params = prepare_vasp_parameters(
                base_parameters=scf_params,
                code=code,
                potential_family=potential_family,
                potential_mapping=slab_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=slab_kpoints_spacing,
                options=slab_options or bulk_options,
                clean_workdir=clean_workdir,
            )
            scf_nodes[slab_id] = wg.add_task(
                VaspWorkChain,
                name=f"scf_slab_{slab_id}",
                structure=slab_structure,
                **scf_vasp_params,
            )

        # Extract SCF energy
        unrelaxed_energy_nodes[slab_id] = wg.add_task(
            extract_total_energy,
            name=f"extract_scf_energy_{slab_id}",
            energies=scf_nodes[slab_id].outputs.misc,
        )

        # Relaxation calculation (if requested)
        if relax_slabs:
            if slab_relax_builders and slab_id in slab_relax_builders:
                relax_nodes[slab_id] = wg.add_task(
                    VaspWorkChain,
                    name=f"relax_slab_{slab_id}",
                    structure=slab_structure,
                    **slab_relax_builders[slab_id],
                )
            else:
                relax_vasp_params = prepare_vasp_parameters(
                    base_parameters=slab_parameters,
                    code=code,
                    potential_family=potential_family,
                    potential_mapping=slab_potential_mapping or bulk_potential_mapping,
                    kpoints_spacing=slab_kpoints_spacing,
                    options=slab_options or bulk_options,
                    clean_workdir=clean_workdir,
                )
                relax_nodes[slab_id] = wg.add_task(
                    VaspWorkChain,
                    name=f"relax_slab_{slab_id}",
                    structure=slab_structure,
                    **relax_vasp_params,
                )

            # Extract relaxation energy
            relaxed_energy_nodes[slab_id] = wg.add_task(
                extract_total_energy,
                name=f"extract_relaxed_energy_{slab_id}",
                energies=relax_nodes[slab_id].outputs.misc,
            )

    # Calculate relaxation energies (if requested)
    relaxation_energy_nodes = {}
    if relax_slabs and compute_relaxation_energy:
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

        # Build slab structures dict for thermodynamics
        slab_structures_for_thermo = {}
        for slab_id in expected_slab_ids:
            if input_slabs:
                slab_structures_for_thermo[slab_id] = input_slabs[slab_id]
            else:
                slab_structures_for_thermo[slab_id] = getattr(slab_gen_node.outputs.slabs, slab_id)

        # Build surface energy nodes
        surface_energy_nodes = build_surface_energy_nodes(
            wg=wg,
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            slab_structures=slab_structures_for_thermo,
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
    }

    # Note: slab_structures are inputs, not outputs, so we don't return them
    # The input_slabs dict is available from the workgraph's execution context

    if compute_thermodynamics:
        # Add reference energies (only if they exist)
        if 'metal' in reference_energy_nodes:
            outputs['metal_energy'] = reference_energy_nodes['metal'].outputs.result
        if 'oxygen' in reference_energy_nodes:
            outputs['oxygen_energy'] = reference_energy_nodes['oxygen'].outputs.result
        if 'nonmetal' in reference_energy_nodes:
            outputs['nonmetal_energy'] = reference_energy_nodes['nonmetal'].outputs.result

        # Add reference structures (only if they exist)
        if 'metal' in reference_nodes:
            outputs['metal_structure'] = reference_nodes['metal'].outputs.structure
        if 'oxygen' in reference_nodes:
            outputs['oxygen_structure'] = reference_nodes['oxygen'].outputs.structure
        if 'nonmetal' in reference_nodes:
            outputs['nonmetal_structure'] = reference_nodes['nonmetal'].outputs.structure

        # Add thermodynamics results
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

    # Set workgraph outputs
    wg.outputs = outputs

    return wg


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
