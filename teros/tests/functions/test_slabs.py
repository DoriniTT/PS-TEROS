import pytest
from aiida.orm import StructureData, List, Float, Bool, Int
from pymatgen.core import Lattice, Structure
from teros.functions.slabs import get_slabs

@pytest.fixture
def bulk_si_structure():
    """
    Returns an AiiDA StructureData object for bulk Silicon.
    """
    # Create a pymatgen Structure for Silicon using spacegroup
    structure = Structure.from_spacegroup("Fd-3m", Lattice.cubic(5.43), ["Si"], [[0, 0, 0]])
    
    # Convert to AiiDA StructureData
    structure_data = StructureData(pymatgen_structure=structure)
    return structure_data

# Placeholder for actual tests that will use the fixture
def test_example_usage_of_fixture(bulk_si_structure):
    assert isinstance(bulk_si_structure, StructureData)
    # The Fd-3m spacegroup for Si results in a conventional cell with 8 atoms.
    assert len(bulk_si_structure.sites) == 8 
    assert bulk_si_structure.get_formula() == "Si8" # Or "Si" if AiiDA normalizes it. AiiDA/Pymatgen might simplify this.
    # Let's check the formula more robustly
    assert "Si" in bulk_si_structure.get_formula()

    # To be more precise for Fd-3m Si conventional cell:
    # A conventional cell of Silicon (diamond structure, Fd-3m) has 8 atoms.
    # The primitive cell has 2 atoms.
    # Structure.from_spacegroup by default creates a conventional cell.
    # So, we expect 8 atoms in the structure.
    
    # Let's get the pymatgen structure back to check the number of sites if AiiDA simplified it
    pmg_structure = bulk_si_structure.get_pymatgen()
    assert len(pmg_structure.sites) == 8
    assert pmg_structure.formula == "Si8"


def test_get_slabs_basic(bulk_si_structure):
    """
    Test the get_slabs function with basic inputs for Si(100).
    """
    inputs = {
        "relaxed_structure": bulk_si_structure,
        "miller_indices": List(list=[1,0,0]),
        "min_slab_thickness": Float(10.0),
        "min_vacuum_thickness": Float(10.0),
        "lll_reduce": Bool(False),
        "center_slab": Bool(True),
        "symmetrize": Bool(False), # Get all unique terminations
        "primitive": Bool(True),   # Use primitive cell of bulk Si
        "in_unit_planes": Bool(False),
        # max_normal_search is intentionally omitted to use the default
    }

    # Call the calcfunction
    # For calcfunctions, the direct call returns the output dictionary
    output = get_slabs(**inputs)

    # Perform assertions
    assert isinstance(output, dict), "Output should be a dictionary."
    assert "structures" in output, "Key 'structures' missing in output."
    assert isinstance(output["structures"], dict), "'structures' should be a dictionary."
    assert len(output["structures"]) > 0, "No slab structures were generated."

    for slab_name, slab_structure in output["structures"].items():
        assert isinstance(slab_structure, StructureData), \
            f"Slab '{slab_name}' is not an AiiDA StructureData object."
        # Further checks on slab properties can be added here, e.g.,
        # - periodicity (should be 2D periodic, 3rd cell vector non-periodic)
        # - minimum thickness and vacuum (approximate, due to discrete atomic layers)
        # - number of sites (should be reasonable for Si(100) slab)
        # - Miller indices if stored in attributes or extras

        # Check for vacuum along z by looking at pbc
        assert slab_structure.pbc == (True, True, False), \
            f"Slab '{slab_name}' should have 2D periodicity (True, True, False)."

        # Check that the slab is indeed derived from Silicon
        assert "Si" in slab_structure.get_formula(), \
            f"Slab '{slab_name}' formula {slab_structure.get_formula()} does not contain Si."


def test_get_slabs_symmetrize(bulk_si_structure):
    """
    Test the get_slabs function with symmetrize=True and symmetrize=False
    and compare the number of generated slabs for Si(111).
    """
    base_inputs = {
        "relaxed_structure": bulk_si_structure,
        "miller_indices": List(list=[1,1,1]), # Miller (1,1,1)
        "min_slab_thickness": Float(10.0),
        "min_vacuum_thickness": Float(10.0),
        "lll_reduce": Bool(False),
        "center_slab": Bool(True),
        "primitive": Bool(True),
        "in_unit_planes": Bool(False),
    }

    # First call: symmetrize=True
    inputs_symm_true = base_inputs.copy()
    inputs_symm_true["symmetrize"] = Bool(True)
    output_symm_true = get_slabs(**inputs_symm_true)
    num_slabs_symm_true = len(output_symm_true["structures"])

    # Second call: symmetrize=False
    inputs_symm_false = base_inputs.copy()
    inputs_symm_false["symmetrize"] = Bool(False)
    output_symm_false = get_slabs(**inputs_symm_false)
    num_slabs_symm_false = len(output_symm_false["structures"])

    # Perform assertions
    assert num_slabs_symm_true > 0, "No slabs generated with symmetrize=True for Si(111)"
    assert num_slabs_symm_false > 0, "No slabs generated with symmetrize=False for Si(111)"
    
    assert num_slabs_symm_true <= num_slabs_symm_false, \
        f"Number of slabs with symmetrize=True ({num_slabs_symm_true}) " \
        f"should be less than or equal to symmetrize=False ({num_slabs_symm_false})."


def test_get_slabs_thickness(bulk_si_structure):
    """
    Test that generated slabs meet the minimum slab and vacuum thickness requirements.
    """
    min_slab_val = 15.0
    min_vacuum_val = 20.0
    tolerance = 1.0  # Angstroms

    inputs = {
        "relaxed_structure": bulk_si_structure,
        "miller_indices": List(list=[1,0,0]),
        "min_slab_thickness": Float(min_slab_val),
        "min_vacuum_thickness": Float(min_vacuum_val),
        "lll_reduce": Bool(False),
        "center_slab": Bool(True),
        "symmetrize": Bool(False),
        "primitive": Bool(True),
        "in_unit_planes": Bool(False),
    }

    output = get_slabs(**inputs)

    assert "structures" in output and len(output["structures"]) > 0, \
        "No slabs generated for thickness test Si(100)"

    for slab_name, slab_sd in output["structures"].items():
        slab_pmg = slab_sd.get_pymatgen()

        # Check 1: Total cell height (c-lattice vector)
        # This is what Pymatgen's SlabGenerator tries to achieve for slab_thickness + vacuum_thickness
        total_cell_c = slab_pmg.lattice.c
        expected_total_height = min_slab_val + min_vacuum_val
        assert total_cell_c >= expected_total_height - tolerance, \
            (f"Slab '{slab_name}': Cell c-vector length {total_cell_c:.2f} Å is less than "
             f"min_slab ({min_slab_val} Å) + min_vacuum ({min_vacuum_val} Å) - tolerance ({tolerance} Å) = {expected_total_height - tolerance:.2f} Å.")

        # Check 2: Detailed check of actual atomic slab thickness and vacuum
        if not slab_pmg.sites: # Skip if no atoms in the structure
            continue

        z_coords = slab_pmg.cart_coords[:, 2]
        actual_slab_extent = max(z_coords) - min(z_coords)

        # Assert actual slab thickness
        assert actual_slab_extent >= min_slab_val - tolerance, \
            (f"Slab '{slab_name}': Actual atomic slab extent {actual_slab_extent:.2f} Å is less than "
             f"min_slab_thickness ({min_slab_val} Å) - tolerance ({tolerance} Å) = {min_slab_val - tolerance:.2f} Å.")

        # Assert actual vacuum thickness
        # The vacuum is on both sides if centered, or one side otherwise.
        # The total vacuum in the cell is total_cell_c - actual_slab_extent.
        actual_vacuum_gap = total_cell_c - actual_slab_extent
        assert actual_vacuum_gap >= min_vacuum_val - tolerance, \
            (f"Slab '{slab_name}': Actual vacuum gap {actual_vacuum_gap:.2f} Å is less than "
             f"min_vacuum_thickness ({min_vacuum_val} Å) - tolerance ({tolerance} Å) = {min_vacuum_val - tolerance:.2f} Å.")


def test_get_slabs_primitive(bulk_si_structure):
    """
    Test the get_slabs function with primitive=True and primitive=False
    and compare the number of atoms in the generated slabs for Si(100).
    """
    base_inputs = {
        "relaxed_structure": bulk_si_structure, # Conventional cell fixture
        "miller_indices": List(list=[1,0,0]),
        "min_slab_thickness": Float(10.0),
        "min_vacuum_thickness": Float(10.0),
        "lll_reduce": Bool(False),
        "center_slab": Bool(True),
        "symmetrize": Bool(False),
        "in_unit_planes": Bool(False),
    }

    # First call: primitive=True
    inputs_primitive_true = base_inputs.copy()
    inputs_primitive_true["primitive"] = Bool(True)
    output_primitive_true = get_slabs(**inputs_primitive_true)

    assert "structures" in output_primitive_true and len(output_primitive_true["structures"]) > 0, \
        "No slabs generated with primitive=True for Si(100)"
    slab_primitive = list(output_primitive_true["structures"].values())[0] # Get the first slab
    num_atoms_primitive = len(slab_primitive.sites)

    # Second call: primitive=False
    inputs_primitive_false = base_inputs.copy()
    inputs_primitive_false["primitive"] = Bool(False)
    output_primitive_false = get_slabs(**inputs_primitive_false)
    
    assert "structures" in output_primitive_false and len(output_primitive_false["structures"]) > 0, \
        "No slabs generated with primitive=False for Si(100)"
    slab_conventional = list(output_primitive_false["structures"].values())[0] # Get the first slab
    num_atoms_conventional = len(slab_conventional.sites)

    # Perform assertions
    # For Si (100) and a reasonable thickness, slabs from primitive bulk should have fewer atoms
    # than slabs from conventional bulk, because the cross-sectional area of the slab unit cell
    # derived from the primitive bulk cell is smaller.
    assert num_atoms_primitive < num_atoms_conventional, \
        (f"Number of atoms in slab from primitive bulk ({num_atoms_primitive}) "
         f"is not less than from conventional bulk ({num_atoms_conventional}) for Si(100). "
         f"This might occur for very thin slabs or specific Miller indices/thickness combinations. "
         f"Primitive slab cell: {slab_primitive.cell}, Conventional slab cell: {slab_conventional.cell}")
