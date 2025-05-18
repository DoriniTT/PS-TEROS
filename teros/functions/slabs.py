from aiida.orm import StructureData, Int
from aiida_workgraph import task
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from aiida.orm import StructureData, List, Float, Bool, Int

@task.calcfunction(outputs=[{"name": "structures"}])
def get_slabs(
    relaxed_structure: StructureData,
    miller_indices: List,
    min_slab_thickness: Float,
    min_vacuum_thickness: Float,
    lll_reduce: Bool,
    center_slab: Bool,
    symmetrize: Bool,
    primitive: Bool,
    in_unit_planes: Bool,
    max_normal_search: Int = None
):
    """
    Generate slab structures from a bulk crystal structure.

    Args:
        relaxed_structure: AiiDA StructureData of the bulk crystal
        miller_indices: AiiDA List for Miller indices (e.g., [1, 0, 0])
        min_slab_thickness: AiiDA Float for minimum slab thickness in Angstroms
        min_vacuum_thickness: AiiDA Float for minimum vacuum thickness in Angstroms
        lll_reduce: AiiDA Bool whether to reduce the cell using LLL algorithm
        center_slab: AiiDA Bool whether to center the slab
        symmetrize: AiiDA Bool whether to generate symmetrically distinct terminations
        primitive: AiiDA Bool whether to use the primitive cell
        in_unit_planes: AiiDA Bool whether to restrict to unit planes
        max_normal_search: AiiDA Int for max normal search (optional)

    Returns:
        Dictionary with 'structures' key containing different slab terminations
    """
    # --- Helper functions for structure conversion ---
    def get_pymatgen_structure(structure):
        """Convert AiiDA StructureData to pymatgen Structure."""
        return structure.get_pymatgen()

    def get_aiida_structure(structure):
        """Convert pymatgen Structure to AiiDA StructureData."""
        ase_atoms = AseAtomsAdaptor().get_atoms(structure)
        return StructureData(ase=ase_atoms)

    # --- Slab generation parameters ---
    # Convert AiiDA types to Python types
    py_miller_indices = tuple(miller_indices.get_list())
    py_min_slab_thickness = min_slab_thickness.value
    py_min_vacuum_thickness = min_vacuum_thickness.value
    py_lll_reduce = lll_reduce.value
    py_center_slab = center_slab.value
    py_symmetrize = symmetrize.value
    py_primitive = primitive.value
    py_in_unit_planes = in_unit_planes.value
    py_max_normal_search = max_normal_search.value if max_normal_search is not None else None

    # --- Convert input structure to pymatgen format ---
    bulk_structure = get_pymatgen_structure(relaxed_structure)

    # --- Optionally reduce to primitive cell for cleaner slabs ---
    if py_primitive:
        analyzer = SpacegroupAnalyzer(bulk_structure)
        bulk_structure = analyzer.get_primitive_standard_structure()

    # --- Create the slab generator object ---
    slab_gen = SlabGenerator(
        bulk_structure,
        py_miller_indices,
        py_min_slab_thickness,
        py_min_vacuum_thickness,
        lll_reduce=py_lll_reduce,
        center_slab=py_center_slab,
        max_normal_search=py_max_normal_search,
        in_unit_planes=py_in_unit_planes
    )

    # --- Generate all possible slabs for the given orientation ---
    slabs = slab_gen.get_slabs(symmetrize=py_symmetrize)

    # --- Convert slabs to orthogonal cells and then to AiiDA structures ---
    aiida_slabs = {}
    for i, slab in enumerate(slabs):
        ortho_slab = slab.get_orthogonal_c_slab()  # Convert to orthogonal cell along c-axis
        super_slab = ortho_slab.make_supercell((1, 1, 1))  # No expansion, but could be changed
        aiida_slabs[f"s_{i}"] = get_aiida_structure(super_slab)  # Convert to AiiDA StructureData

    # --- Return the dictionary of slab structures ---
    return {'structures': aiida_slabs}