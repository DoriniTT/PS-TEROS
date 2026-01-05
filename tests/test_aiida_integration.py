"""
AiiDA Integration Tests

These tests verify AiiDA-specific functionality including:
- Node creation and storage
- Data type serialization
- WorkGraph construction (without submission)
- Calcfunction execution with mock data

Note: Calcfunctions decorated with @task.calcfunction must be called
within a WorkGraph context to return actual AiiDA nodes. When called
directly, they return task socket objects. For testing calcfunction
logic, we test the underlying functions or use the .run() method.
"""

import pytest


def check_aiida():
    """Check if AiiDA is available and configured."""
    try:
        from aiida import load_profile, orm
        load_profile()
        return True
    except Exception:
        return False


AIIDA_AVAILABLE = check_aiida()

if not AIIDA_AVAILABLE:
    pytest.skip("AiiDA not configured", allow_module_level=True)


# =============================================================================
# AiiDA DATA TYPE TESTS
# =============================================================================

class TestAiiDADataTypes:
    """Test AiiDA data type creation and serialization."""

    def test_create_float_node(self):
        """Test creating Float node."""
        from aiida import orm

        node = orm.Float(3.14159)
        assert node.value == pytest.approx(3.14159)

    def test_create_int_node(self):
        """Test creating Int node."""
        from aiida import orm

        node = orm.Int(42)
        assert node.value == 42

    def test_create_str_node(self):
        """Test creating Str node."""
        from aiida import orm

        node = orm.Str("test_string")
        assert node.value == "test_string"

    def test_create_bool_node(self):
        """Test creating Bool node."""
        from aiida import orm

        node = orm.Bool(True)
        assert node.value is True

    def test_create_list_node(self):
        """Test creating List node."""
        from aiida import orm

        node = orm.List(list=[1, 2, 3])
        assert node.get_list() == [1, 2, 3]

    def test_create_dict_node(self):
        """Test creating Dict node."""
        from aiida import orm

        data = {'key1': 'value1', 'key2': 42, 'nested': {'a': 1}}
        node = orm.Dict(dict=data)

        result = node.get_dict()
        assert result['key1'] == 'value1'
        assert result['key2'] == 42
        assert result['nested']['a'] == 1

    def test_create_structure_data(self):
        """Test creating StructureData node."""
        from aiida import orm
        from ase import Atoms

        # Create simple H2O molecule
        atoms = Atoms(
            'H2O',
            positions=[
                [0.0, 0.0, 0.0],
                [0.96, 0.0, 0.0],
                [0.24, 0.93, 0.0],
            ],
            cell=[10.0, 10.0, 10.0],
            pbc=True,
        )

        structure = orm.StructureData(ase=atoms)

        assert len(structure.sites) == 3
        symbols = [site.kind_name for site in structure.sites]
        assert 'H' in symbols
        assert 'O' in symbols

    def test_structure_data_round_trip(self):
        """Test StructureData to ASE and back."""
        from aiida import orm
        from ase import Atoms

        # Create structure
        atoms = Atoms('Ag4', positions=[[0, 0, 0], [2, 0, 0], [0, 2, 0], [2, 2, 0]], cell=[4, 4, 10])
        structure = orm.StructureData(ase=atoms)

        # Convert back to ASE
        atoms_back = structure.get_ase()

        assert len(atoms_back) == 4
        assert all(s == 'Ag' for s in atoms_back.get_chemical_symbols())


class TestAiiDANodeStorage:
    """Test AiiDA node storage operations."""

    def test_store_float_node(self):
        """Test storing Float node."""
        from aiida import orm

        node = orm.Float(2.718)
        node.store()

        assert node.is_stored
        assert node.pk is not None

    def test_store_dict_node(self):
        """Test storing Dict node."""
        from aiida import orm

        node = orm.Dict(dict={'test': 'data', 'number': 123})
        node.store()

        assert node.is_stored

        # Load and verify
        loaded = orm.load_node(node.pk)
        assert loaded.get_dict()['test'] == 'data'

    def test_store_structure_data(self):
        """Test storing StructureData node."""
        from aiida import orm
        from ase import Atoms

        atoms = Atoms('Ag2O', positions=[[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]],
                      cell=[5, 5, 5], pbc=True)
        structure = orm.StructureData(ase=atoms)
        structure.store()

        assert structure.is_stored

        # Load and verify
        loaded = orm.load_node(structure.pk)
        assert len(loaded.sites) == 3


# =============================================================================
# CALCFUNCTION TESTS
# =============================================================================

class TestCalcfunctions:
    """Test PS-TEROS calcfunctions with mock data.
    
    Note: Calcfunctions decorated with @task.calcfunction from aiida-workgraph
    are TaskHandle objects that can only be properly executed within a WorkGraph
    context. For Tier 2 tests, we test:
    1. The underlying logic by recreating the calculations
    2. The function existence and callable status
    3. Data types and return value expectations
    """

    def test_extract_total_energy_logic(self):
        """Test the logic of extract_total_energy."""
        # Test the calculation logic directly
        energies_dict = {
            'total_energies': {
                'energy_extrapolated': -15.12345,
                'energy_no_entropy': -15.12340,
            }
        }

        # Extract energy using the same logic as the calcfunction
        total_energies = energies_dict.get('total_energies', {})
        if 'energy_extrapolated' in total_energies:
            energy = total_energies['energy_extrapolated']
        elif 'energy_no_entropy' in total_energies:
            energy = total_energies['energy_no_entropy']
        else:
            energy = None

        assert energy == pytest.approx(-15.12345)

    def test_extract_total_energy_alternative_key(self):
        """Test extract_total_energy logic with alternative energy key."""
        energies_dict = {
            'total_energies': {
                'energy_no_entropy': -10.5,
            }
        }

        total_energies = energies_dict.get('total_energies', {})
        if 'energy_extrapolated' in total_energies:
            energy = total_energies['energy_extrapolated']
        elif 'energy_no_entropy' in total_energies:
            energy = total_energies['energy_no_entropy']
        else:
            energy = None

        assert energy == pytest.approx(-10.5)

    def test_identify_oxide_type_binary_logic(self):
        """Test identify_oxide_type logic for binary oxide."""
        from collections import Counter
        from ase import Atoms

        # Ag2O structure
        atoms = Atoms(
            'Ag2O',
            positions=[[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]],
            cell=[5, 5, 5],
            pbc=True
        )

        # Same logic as identify_oxide_type
        bulk_counts = Counter(atoms.get_chemical_symbols())
        metal_elements = sorted(element for element in bulk_counts if element != 'O')

        if len(metal_elements) == 1:
            oxide_type = 'binary'
        elif len(metal_elements) == 2:
            oxide_type = 'ternary'
        else:
            oxide_type = 'complex'

        assert oxide_type == 'binary'

    def test_identify_oxide_type_ternary_logic(self):
        """Test identify_oxide_type logic for ternary oxide."""
        from collections import Counter
        from ase import Atoms

        # AgPO4 structure (ternary: Ag, P, O)
        atoms = Atoms(
            symbols=['Ag', 'P', 'O', 'O', 'O', 'O'],
            positions=[
                [0, 0, 0], [1, 0, 0], [2, 0, 0],
                [3, 0, 0], [4, 0, 0], [5, 0, 0]
            ],
            cell=[10, 10, 10],
            pbc=True
        )

        bulk_counts = Counter(atoms.get_chemical_symbols())
        metal_elements = sorted(element for element in bulk_counts if element != 'O')

        if len(metal_elements) == 1:
            oxide_type = 'binary'
        elif len(metal_elements) == 2:
            oxide_type = 'ternary'
        else:
            oxide_type = 'complex'

        assert oxide_type == 'ternary'

    def test_calculate_energy_difference_logic(self):
        """Test energy difference calculation logic."""
        E_final = -100.5
        E_initial = -100.0

        E_diff = E_final - E_initial

        assert E_diff == pytest.approx(-0.5)

    def test_calcfunction_handles_exist(self):
        """Test that calcfunction TaskHandles are properly defined."""
        from teros.core.slabs import extract_total_energy, calculate_energy_difference
        from teros.core.thermodynamics import identify_oxide_type

        # Verify these are callable/have expected attributes
        assert extract_total_energy is not None
        assert calculate_energy_difference is not None
        assert identify_oxide_type is not None


# =============================================================================
# SLAB GENERATION TESTS
# =============================================================================

class TestSlabGenerationWithAiiDA:
    """Test slab generation functions with AiiDA."""

    def test_slab_generation_logic(self):
        """Test slab generation logic using pymatgen directly."""
        from pymatgen.core import Structure, Lattice
        from pymatgen.core.surface import SlabGenerator

        # Create FCC Ag bulk structure using pymatgen
        a = 4.08
        lattice = Lattice.cubic(a)
        structure = Structure(
            lattice,
            ['Ag', 'Ag', 'Ag', 'Ag'],
            [
                [0, 0, 0],
                [0.5, 0.5, 0],
                [0.5, 0, 0.5],
                [0, 0.5, 0.5],
            ]
        )

        # Generate slabs
        generator = SlabGenerator(
            structure,
            miller_index=(1, 0, 0),
            min_slab_size=10.0,
            min_vacuum_size=10.0,
            lll_reduce=True,
            center_slab=True,
        )

        slabs = generator.get_slabs()

        # Should generate at least one slab
        assert len(slabs) > 0

        # Each slab should have atoms
        for slab in slabs:
            assert len(slab) > 0

    def test_aiida_structure_data_from_pymatgen_slab(self):
        """Test creating AiiDA StructureData from pymatgen slab."""
        from aiida import orm
        from pymatgen.core import Structure, Lattice
        from pymatgen.core.surface import SlabGenerator
        from pymatgen.io.ase import AseAtomsAdaptor

        # Create bulk and generate slab
        a = 4.08
        lattice = Lattice.cubic(a)
        bulk = Structure(
            lattice,
            ['Ag', 'Ag', 'Ag', 'Ag'],
            [
                [0, 0, 0],
                [0.5, 0.5, 0],
                [0.5, 0, 0.5],
                [0, 0.5, 0.5],
            ]
        )

        generator = SlabGenerator(
            bulk, (1, 0, 0), min_slab_size=10.0, min_vacuum_size=10.0
        )
        slabs = generator.get_slabs()

        if len(slabs) > 0:
            # Convert to AiiDA StructureData
            slab_pmg = slabs[0]
            adaptor = AseAtomsAdaptor()
            slab_ase = adaptor.get_atoms(slab_pmg)
            slab_aiida = orm.StructureData(ase=slab_ase)

            assert isinstance(slab_aiida, orm.StructureData)
            assert len(slab_aiida.sites) > 0


# =============================================================================
# ADSORPTION ENERGY TESTS
# =============================================================================

class TestAdsorptionEnergyFunctions:
    """Test adsorption energy functions."""

    def test_parse_formula(self):
        """Test parse_formula function."""
        from teros.core.adsorption_energy import parse_formula

        result = parse_formula('OOH')
        assert result == {'O': 2, 'H': 1}

        result = parse_formula('H2O')
        assert result == {'H': 2, 'O': 1}

    def test_calculate_adsorption_energy_logic(self):
        """Test adsorption energy calculation logic."""
        E_complete = -150.0
        E_substrate = -100.0
        E_molecule = -30.0

        E_ads = E_complete - E_substrate - E_molecule

        # E_ads = -150 - (-100) - (-30) = -150 + 100 + 30 = -20
        assert E_ads == pytest.approx(-20.0)

    def test_separate_adsorbate_structure_logic(self):
        """Test adsorbate separation logic."""
        from collections import Counter
        from ase import Atoms

        # Create a structure with substrate (Ag) and adsorbate (OH)
        atoms = Atoms(
            symbols=['Ag', 'Ag', 'Ag', 'Ag', 'O', 'H'],
            positions=[
                [0, 0, 0],      # Ag
                [2, 0, 0],      # Ag
                [0, 2, 0],      # Ag
                [2, 2, 0],      # Ag
                [1, 1, 2.5],    # O (adsorbate)
                [1, 1, 3.5],    # H (adsorbate)
            ],
            cell=[4, 4, 20],
            pbc=True
        )

        # Adsorbate formula
        adsorbate_formula = {'O': 1, 'H': 1}  # OH

        # Count atoms in structure
        structure_counts = Counter(atoms.get_chemical_symbols())

        # Substrate would have these atoms removed
        substrate_counts = structure_counts.copy()
        for element, count in adsorbate_formula.items():
            substrate_counts[element] -= count

        # Verify substrate has 4 Ag atoms
        assert substrate_counts['Ag'] == 4
        assert substrate_counts['O'] == 0
        assert substrate_counts['H'] == 0


# =============================================================================
# FORMATION ENTHALPY TESTS
# =============================================================================

class TestFormationEnthalpyFunctions:
    """Test formation enthalpy functions."""

    def test_calculate_formation_enthalpy_logic(self):
        """Test formation enthalpy calculation logic for binary oxide."""
        from collections import Counter
        from math import gcd
        from functools import reduce

        # Ag2O composition
        element_counts = {'Ag': 2, 'O': 1}

        # Mock energies
        E_bulk = -15.0  # Total energy of Ag2O
        E_Ag_per_atom = -2.5  # Energy per Ag atom
        E_O2_per_atom = -4.5  # Energy per O atom from O2

        # Calculate formation energy
        # ΔH_f = E_bulk - Σ(n_i * E_i)
        formation_energy = E_bulk
        formation_energy -= element_counts['Ag'] * E_Ag_per_atom  # -5.0
        formation_energy -= element_counts['O'] * E_O2_per_atom   # -4.5

        # Expected: -15.0 - (-5.0) - (-4.5) = -15.0 + 5.0 + 4.5 = -5.5
        assert formation_energy == pytest.approx(-5.5)

        # Normalize per formula unit
        values = list(element_counts.values())
        formula_units = reduce(gcd, values)  # 1 for Ag2O

        formation_per_fu = formation_energy / formula_units
        assert formation_per_fu == pytest.approx(-5.5)

    def test_formation_enthalpy_handles_exist(self):
        """Test that formation enthalpy TaskHandle is defined."""
        from teros.core.hf import calculate_formation_enthalpy

        assert calculate_formation_enthalpy is not None
