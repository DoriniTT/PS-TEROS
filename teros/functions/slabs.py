from aiida.orm import StructureData, Int
from aiida_workgraph import task
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from aiida.orm import StructureData, List, Float, Bool, Int
import csv

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
    Generate slab structures from a bulk crystal structure using Pymatgen's SlabGenerator.

    This calcfunction wraps Pymatgen's slab generation capabilities to produce various
    surface terminations for a given bulk material and Miller index.

    :param relaxed_structure: AiiDA ``StructureData`` node of the bulk crystal.
    :type relaxed_structure: aiida.orm.StructureData
    :param miller_indices: AiiDA ``List`` of the Miller indices for slab generation (e.g., ``List(list=[1,0,0])``).
    :type miller_indices: aiida.orm.List
    :param min_slab_thickness: AiiDA ``Float`` specifying the minimum slab thickness in Angstroms.
    :type min_slab_thickness: aiida.orm.Float
    :param min_vacuum_thickness: AiiDA ``Float`` specifying the minimum vacuum thickness in Angstroms.
    :type min_vacuum_thickness: aiida.orm.Float
    :param lll_reduce: AiiDA ``Bool`` indicating whether to reduce the cell using the LLL algorithm
                       before slab generation to get a more orthogonal cell. Default: False.
    :type lll_reduce: aiida.orm.Bool
    :param center_slab: AiiDA ``Bool`` indicating whether to center the slab in the c direction of the cell. Default: True.
    :type center_slab: aiida.orm.Bool
    :param symmetrize: AiiDA ``Bool`` indicating whether to generate symmetrically distinct terminations.
                       If False, generates all unique terminations. Default: False.
    :type symmetrize: aiida.orm.Bool
    :param primitive: AiiDA ``Bool`` indicating whether to find the primitive cell of the bulk structure
                      before generating slabs. Default: True.
    :type primitive: aiida.orm.Bool
    :param in_unit_planes: AiiDA ``Bool`` indicating whether to restrict Miller indices to unit planes. Default: False.
    :type in_unit_planes: aiida.orm.Bool
    :param max_normal_search: An optional AiiDA ``Int``. Pymatgen parameter: max normal search for finding Miller indices.
                              Corresponds to ``max_sites`` in Pymatgen's ``SlabGenerator.get_slabs()`` if ``symmetrize=False``.
                              Default: None (uses Pymatgen's default).
    :type max_normal_search: aiida.orm.Int, optional

    :return: A dictionary where the key 'structures' maps to another dictionary.
             This inner dictionary has slab identifiers (e.g., "s_0", "s_1") as keys
             and AiiDA ``StructureData`` nodes of the generated slabs as values.
             Example: ``{'structures': {'s_0': <StructureData>, 's_1': <StructureData>}}``
    :rtype: dict[str, dict[str, aiida.orm.StructureData]]
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

    # --- Helper function to analyze stoichiometric deviations ---
    def analyze_defects(slab_structure, bulk_structure):
        """
        Analyze stoichiometric deviations in a slab compared to the bulk structure.
        
        :param slab_structure: pymatgen Structure of the slab
        :param bulk_structure: pymatgen Structure of the bulk
        :return: dict with element symbols as keys and excess amounts as values
        """
        # Get bulk composition and normalize to per formula unit
        bulk_composition = bulk_structure.composition
        # Get the reduced composition and factor
        bulk_reduced, factor = bulk_composition.get_reduced_composition_and_factor()
        
        # Get slab composition
        slab_composition = slab_structure.composition
        
        # Calculate excess for each element in the bulk
        defects = {}
        for element in bulk_reduced.elements:
            bulk_amount = bulk_reduced[element]
            # Scale to slab size (using the same element to determine scaling)
            slab_amount = slab_composition[element] if element in slab_composition else 0
            # Calculate excess relative to the reduced bulk formula
            excess = slab_amount - bulk_amount
            defects[element.name] = excess
            
        # Check for elements in slab that are not in bulk
        for element in slab_composition.elements:
            if element not in bulk_reduced.elements:
                defects[element.name] = slab_composition[element]
                
        return defects

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
    defect_data = []
    
    # Get all elements in the bulk structure for consistent CSV columns
    bulk_elements = [element.name for element in bulk_structure.composition.elements]
    
    for i, slab in enumerate(slabs):
        ortho_slab = slab.get_orthogonal_c_slab()  # Convert to orthogonal cell along c-axis
        super_slab = ortho_slab.make_supercell((1, 1, 1))  # No expansion, but could be changed
        aiida_slabs[f"s_{i}"] = get_aiida_structure(super_slab)  # Convert to AiiDA StructureData
        
        # Analyze defects for this slab
        defects = analyze_defects(super_slab, bulk_structure)
        defect_entry = {"termination": f"s_{i}"}
        defect_entry.update(defects)
        defect_data.append(defect_entry)

    # --- Write defect data to CSV file ---
    # Ensure all elements from bulk are included in CSV columns
    csv_columns = ["termination"] + bulk_elements
    
    # Add any additional elements found in slabs but not in bulk
    all_elements = set(bulk_elements)
    for entry in defect_data:
        all_elements.update(entry.keys())
    all_elements.discard("termination")
    
    # Reorder columns: termination first, then elements in alphabetical order
    csv_columns = ["termination"] + sorted(list(all_elements))
    
    with open("DEFECT_TYPES", "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for entry in defect_data:
            # Ensure all columns are present in each row
            row = {col: entry.get(col, 0) for col in csv_columns}
            writer.writerow(row)

    # --- Return the dictionary of slab structures ---
    return {'structures': aiida_slabs}
