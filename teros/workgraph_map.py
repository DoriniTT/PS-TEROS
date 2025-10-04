"""
Helper module for building TEROS workgraphs with VASP.

This module provides convenience functions to build complete TEROS workgraphs
using VASP as the DFT code, automatically loading structures from files and
creating properly configured builders.
"""
from pathlib import Path
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph
from teros.core.workgraph import create_teros_workgraph


def build_core_workgraph_with_map(
    structures_dir,
    bulk_name,
    metal_name,
    nonmetal_name,
    oxygen_name,
    code_label,
    potential_family,
    bulk_potential_mapping,
    metal_potential_mapping,
    nonmetal_potential_mapping,
    oxygen_potential_mapping,
    kpoints_spacing,
    bulk_parameters,
    bulk_options,
    metal_parameters,
    metal_options,
    nonmetal_parameters,
    nonmetal_options,
    oxygen_parameters,
    oxygen_options,
    clean_workdir=True,
    # Slab generation parameters
    miller_indices=None,
    min_slab_thickness=None,
    min_vacuum_thickness=None,
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
    in_unit_planes=False,
    max_normal_search=None,
    # Slab relaxation parameters
    relax_slabs=True,
    slab_parameters=None,
    slab_options=None,
    slab_potential_mapping=None,
    slab_kpoints_spacing=None,
    name="TEROS_Workflow",
):
    """
    Build a complete TEROS workgraph for VASP calculations.
    
    This function loads structure files, creates VASP builders with proper
    pseudopotentials and parameters, and assembles them into a complete
    TEROS workgraph.
    
    :param structures_dir: Directory containing structure files
    :param bulk_name: Filename of bulk structure CIF
    :param metal_name: Filename of metal reference structure CIF
    :param nonmetal_name: Filename of non-metal reference structure CIF
    :param oxygen_name: Filename of O2 molecule structure CIF
    :param code_label: AiiDA code label for VASP
    :param potential_family: VASP potential family name
    :param bulk_potential_mapping: Dict mapping element symbols to potential names for bulk
    :param metal_potential_mapping: Dict mapping element symbols to potential names for metal
    :param nonmetal_potential_mapping: Dict mapping element symbols to potential names for nonmetal
    :param oxygen_potential_mapping: Dict mapping element symbols to potential names for oxygen
    :param kpoints_spacing: K-points spacing for bulk and references
    :param bulk_parameters: Dict of VASP INCAR parameters for bulk
    :param bulk_options: Dict of scheduler options for bulk
    :param metal_parameters: Dict of VASP INCAR parameters for metal
    :param metal_options: Dict of scheduler options for metal
    :param nonmetal_parameters: Dict of VASP INCAR parameters for nonmetal
    :param nonmetal_options: Dict of scheduler options for nonmetal
    :param oxygen_parameters: Dict of VASP INCAR parameters for oxygen
    :param oxygen_options: Dict of scheduler options for oxygen
    :param clean_workdir: Whether to clean working directories after completion
    :param miller_indices: Miller indices for slab generation
    :param min_slab_thickness: Minimum slab thickness in Angstroms
    :param min_vacuum_thickness: Minimum vacuum thickness in Angstroms
    :param lll_reduce: Whether to apply LLL reduction
    :param center_slab: Whether to center the slab
    :param symmetrize: Whether to generate symmetric terminations
    :param primitive: Whether to use primitive cell
    :param in_unit_planes: Whether to use unit planes
    :param max_normal_search: Maximum normal search parameter
    :param relax_slabs: Whether to relax slabs (not implemented in this version)
    :param slab_parameters: Dict of VASP INCAR parameters for slabs
    :param slab_options: Dict of scheduler options for slabs
    :param slab_potential_mapping: Dict mapping element symbols to potential names for slabs
    :param slab_kpoints_spacing: K-points spacing for slabs
    :param name: Name for the workgraph
    :return: Configured WorkGraph instance
    """
    from aiida_vasp.utils.aiida_utils import get_data_class
    
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.vasp')
    
    # Load structures
    structures_path = Path(structures_dir)
    bulk_structure = orm.StructureData(ase=read(structures_path / bulk_name))
    metal_structure = orm.StructureData(ase=read(structures_path / metal_name))
    nonmetal_structure = orm.StructureData(ase=read(structures_path / nonmetal_name))
    oxygen_structure = orm.StructureData(ase=read(structures_path / oxygen_name))
    
    # Get code
    code = orm.load_code(code_label)
    
    # Helper function to create builder
    def create_vasp_builder(structure, parameters, options, potential_mapping, kpoints_spacing_val):
        """Create a VASP builder with given parameters."""
        builder = VaspWorkChain.get_builder()
        builder.code = code
        builder.structure = structure
        builder.parameters = orm.Dict(dict={'incar': parameters})  # VASP expects parameters under 'incar' key
        builder.options = orm.Dict(dict=options)
        builder.potential_family = orm.Str(potential_family)
        builder.potential_mapping = orm.Dict(dict=potential_mapping)
        
        # Set k-points
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints_mesh_from_density(kpoints_spacing_val)
        builder.kpoints = kpoints
        
        # Clean workdir option
        builder.clean_workdir = orm.Bool(clean_workdir)
        
        return builder
    
    # Create builders for bulk and references
    builder_bulk = create_vasp_builder(
        bulk_structure, bulk_parameters, bulk_options,
        bulk_potential_mapping, kpoints_spacing
    )
    
    # Determine element names from structure
    bulk_elements = set(bulk_structure.get_ase().get_chemical_symbols())
    metals = [e for e in bulk_elements if e != 'O' and e != list(nonmetal_potential_mapping.keys())[0]]
    nonmetal_elem = list(nonmetal_potential_mapping.keys())[0]
    
    # Create reference builders
    reference_builders = {}
    
    # Metal reference
    if metals:
        metal_elem = metals[0]
        reference_builders[metal_elem] = create_vasp_builder(
            metal_structure, metal_parameters, metal_options,
            metal_potential_mapping, kpoints_spacing
        )
    
    # Non-metal reference
    reference_builders[nonmetal_elem] = create_vasp_builder(
        nonmetal_structure, nonmetal_parameters, nonmetal_options,
        nonmetal_potential_mapping, kpoints_spacing
    )
    
    # Oxygen reference
    reference_builders['O'] = create_vasp_builder(
        oxygen_structure, oxygen_parameters, oxygen_options,
        oxygen_potential_mapping, kpoints_spacing
    )
    
    # Create slab builder if slab relaxation is enabled
    if slab_parameters and slab_options and slab_potential_mapping:
        builder_slab = create_vasp_builder(
            bulk_structure,  # Will be replaced with slab structure
            slab_parameters, slab_options,
            slab_potential_mapping, slab_kpoints_spacing
        )
    else:
        builder_slab = None
    
    # Convert parameters to AiiDA types
    miller_indices_aiida = orm.List(list=miller_indices) if miller_indices else orm.List(list=[1, 0, 0])
    min_slab_thickness_aiida = orm.Float(min_slab_thickness) if min_slab_thickness else orm.Float(10.0)
    min_vacuum_thickness_aiida = orm.Float(min_vacuum_thickness) if min_vacuum_thickness else orm.Float(15.0)
    lll_reduce_aiida = orm.Bool(lll_reduce)
    center_slab_aiida = orm.Bool(center_slab)
    symmetrize_aiida = orm.Bool(symmetrize)
    primitive_aiida = orm.Bool(primitive)
    in_unit_planes_aiida = orm.Bool(in_unit_planes)
    max_normal_search_aiida = orm.Int(max_normal_search) if max_normal_search else None
    
    # Create the workgraph
    wg = create_teros_workgraph.build(
        dft_workchain=VaspWorkChain,
        builder_bulk=builder_bulk,
        builder_slab=builder_slab,
        reference_builders=reference_builders,
        workgraph_name=name,
        code="VASP",
        miller_indices=miller_indices_aiida,
        min_slab_thickness=min_slab_thickness_aiida,
        min_vacuum_thickness=min_vacuum_thickness_aiida,
        lll_reduce=lll_reduce_aiida,
        center_slab=center_slab_aiida,
        symmetrize=symmetrize_aiida,
        primitive=primitive_aiida,
        in_unit_planes=in_unit_planes_aiida,
        max_normal_search=max_normal_search_aiida,
    )
    
    return wg


# Import read function
from ase.io import read
