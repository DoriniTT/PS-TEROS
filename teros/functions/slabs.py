from aiida.orm import StructureData
from aiida_workgraph import task
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from aiida.orm import StructureData

@task.calcfunction(outputs=[{"name": "structures"}])
def get_slabs(relaxed_structure):
    """
    Generate slab structures from a bulk crystal structure.

    Args:
        relaxed_structure: AiiDA StructureData of the bulk crystal

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
    miller_indices = (1, 0, 0)           # Miller index for surface orientation (e.g., (100) surface)
    min_slab_thickness = 10.0            # Minimum thickness of the slab in Angstroms
    min_vacuum_thickness = 15.0          # Minimum thickness of vacuum region in Angstroms
    lll_reduce = True                    # Whether to reduce the cell using the LLL algorithm for better numerical stability
    center_slab = True                   # Whether to center the slab in the simulation cell
    symmetrize = True                    # Whether to generate symmetrically distinct slab terminations
    primitive = True                     # Whether to use the primitive cell for slab generation
    max_normal_search = None             # Use default search for surface normal (None lets pymatgen decide)
    in_unit_planes = False               # Whether to restrict to unit planes (False allows more general slabs)

    # --- Convert input structure to pymatgen format ---
    bulk_structure = get_pymatgen_structure(relaxed_structure)

    # --- Optionally reduce to primitive cell for cleaner slabs ---
    if primitive:
        analyzer = SpacegroupAnalyzer(bulk_structure)
        bulk_structure = analyzer.get_primitive_standard_structure()

    # --- Create the slab generator object ---
    slab_gen = SlabGenerator(
        bulk_structure,
        miller_indices,
        min_slab_thickness,
        min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        max_normal_search=max_normal_search,
        in_unit_planes=in_unit_planes
    )

    # --- Generate all possible slabs for the given orientation ---
    slabs = [slab_gen.get_slabs(symmetrize=symmetrize)[0]]

    # --- Convert slabs to orthogonal cells and then to AiiDA structures ---
    aiida_slabs = {}
    for i, slab in enumerate(slabs):
        ortho_slab = slab.get_orthogonal_c_slab()  # Convert to orthogonal cell along c-axis
        super_slab = ortho_slab.make_supercell((1, 1, 1))  # No expansion, but could be changed
        aiida_slabs[f"s_{i}"] = get_aiida_structure(super_slab)  # Convert to AiiDA StructureData

    # --- Return the dictionary of slab structures ---
    return {'structures': aiida_slabs}