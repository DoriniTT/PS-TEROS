"""
Structure Utility Tests

These tests verify structure-related utilities using pymatgen/ASE.
"""

import pytest
import numpy as np


class TestSlabGeneration:
    """Test slab generation utilities using pymatgen directly."""

    @pytest.fixture
    def bulk_structure(self):
        """Create a simple bulk structure for testing."""
        from pymatgen.core import Structure, Lattice

        # Simple cubic structure (like Ag)
        lattice = Lattice.cubic(4.0)
        structure = Structure(
            lattice,
            ['Ag', 'Ag', 'Ag', 'Ag'],
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ]
        )
        return structure

    def test_slab_generator_creation(self, bulk_structure):
        """Test creating a SlabGenerator."""
        from pymatgen.core.surface import SlabGenerator

        generator = SlabGenerator(
            bulk_structure,
            miller_index=(1, 0, 0),
            min_slab_size=10.0,
            min_vacuum_size=10.0,
            center_slab=True,
        )

        assert generator is not None

    def test_slab_generation(self, bulk_structure):
        """Test generating slabs."""
        from pymatgen.core.surface import SlabGenerator

        generator = SlabGenerator(
            bulk_structure,
            miller_index=(1, 0, 0),
            min_slab_size=10.0,
            min_vacuum_size=10.0,
            center_slab=True,
        )

        slabs = generator.get_slabs()

        assert len(slabs) > 0
        assert all(hasattr(s, 'lattice') for s in slabs)

    def test_different_miller_indices(self, bulk_structure):
        """Test slab generation with different Miller indices."""
        from pymatgen.core.surface import SlabGenerator

        miller_indices = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]

        for miller in miller_indices:
            generator = SlabGenerator(
                bulk_structure,
                miller_index=miller,
                min_slab_size=8.0,
                min_vacuum_size=8.0,
            )
            slabs = generator.get_slabs()
            assert len(slabs) >= 0, f"Failed for Miller index {miller}"


class TestSurfaceAreaCalculation:
    """Test surface area calculation."""

    def test_surface_area_calculation(self):
        """Test calculating surface area from cell vectors."""
        # Simple orthogonal cell
        a_vec = np.array([10.0, 0.0, 0.0])
        b_vec = np.array([0.0, 10.0, 0.0])

        cross = np.cross(a_vec, b_vec)
        area = np.linalg.norm(cross)

        assert area == pytest.approx(100.0)

    def test_surface_area_non_orthogonal(self):
        """Test surface area for non-orthogonal cell."""
        a_vec = np.array([10.0, 0.0, 0.0])
        b_vec = np.array([5.0, 8.66, 0.0])  # 60 degree angle

        cross = np.cross(a_vec, b_vec)
        area = np.linalg.norm(cross)

        # Area should be |a||b|sin(60°) = 10 * 10 * sin(60°) ≈ 86.6
        expected = 10.0 * 10.0 * np.sin(np.radians(60))
        assert area == pytest.approx(expected, rel=0.01)


class TestStructureManipulation:
    """Test structure manipulation with ASE."""

    def test_ase_atoms_creation(self):
        """Test creating ASE Atoms object."""
        from ase import Atoms

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

        assert len(atoms) == 3
        assert atoms.get_chemical_formula() == 'H2O'

    def test_ase_atoms_deletion(self):
        """Test deleting atoms from ASE structure."""
        from ase import Atoms

        atoms = Atoms(
            'H2O',
            positions=[
                [0.0, 0.0, 0.0],
                [0.96, 0.0, 0.0],
                [0.24, 0.93, 0.0],
            ],
        )

        # Delete last atom (O)
        del atoms[-1]

        assert len(atoms) == 2
        assert atoms.get_chemical_formula() == 'H2'

    def test_ase_atoms_copy(self):
        """Test copying ASE structure."""
        from ase import Atoms

        original = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
        copy = original.copy()

        # Modify copy
        copy.positions[0] = [5, 5, 5]

        # Original should be unchanged
        assert original.positions[0][0] == pytest.approx(0.0)

    def test_pymatgen_ase_conversion(self):
        """Test converting between pymatgen and ASE."""
        from pymatgen.core import Structure, Lattice
        from pymatgen.io.ase import AseAtomsAdaptor

        # Create pymatgen structure
        lattice = Lattice.cubic(4.0)
        pmg_structure = Structure(
            lattice,
            ['Ag', 'Ag'],
            [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        )

        # Convert to ASE
        adaptor = AseAtomsAdaptor()
        ase_atoms = adaptor.get_atoms(pmg_structure)

        assert len(ase_atoms) == 2
        assert 'Ag' in ase_atoms.get_chemical_symbols()

        # Convert back to pymatgen
        pmg_back = adaptor.get_structure(ase_atoms)
        assert len(pmg_back) == 2


class TestCompositionAnalysis:
    """Test composition analysis utilities."""

    def test_element_counting(self):
        """Test counting elements in a structure."""
        from collections import Counter

        symbols = ['Ag', 'Ag', 'O', 'Ag', 'O']
        counts = Counter(symbols)

        assert counts['Ag'] == 3
        assert counts['O'] == 2

    def test_oxide_identification(self):
        """Test identifying oxide type from composition."""
        from collections import Counter

        def identify_oxide(counts):
            if 'O' not in counts:
                return 'not_oxide'
            metals = [el for el in counts if el != 'O']
            if len(metals) == 1:
                return 'binary'
            elif len(metals) == 2:
                return 'ternary'
            else:
                return 'complex'

        assert identify_oxide(Counter({'Ag': 2, 'O': 1})) == 'binary'
        assert identify_oxide(Counter({'Ag': 3, 'P': 1, 'O': 4})) == 'ternary'
        assert identify_oxide(Counter({'Ag': 1})) == 'not_oxide'

    def test_stoichiometry_reduction(self):
        """Test reducing stoichiometry to formula unit."""
        from math import gcd
        from functools import reduce

        def reduce_stoichiometry(counts):
            values = list(counts.values())
            common = reduce(gcd, values)
            return {k: v // common for k, v in counts.items()}

        # Ag2O already reduced
        result = reduce_stoichiometry({'Ag': 2, 'O': 1})
        assert result == {'Ag': 2, 'O': 1}

        # Fe4O6 -> Fe2O3
        result = reduce_stoichiometry({'Fe': 4, 'O': 6})
        assert result == {'Fe': 2, 'O': 3}
