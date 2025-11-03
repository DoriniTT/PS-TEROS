"""Tests for AIMD tasks."""
from aiida import orm, load_profile
from ase.build import bulk
from teros.core.aimd.tasks import create_supercell

# Load AiiDA profile for tests
load_profile('presto')


def test_create_supercell_basic():
    """Test supercell creation."""
    from ase import Atoms
    # Create simple structure
    atoms = bulk('Al', 'fcc', a=4.0)
    structure = orm.StructureData(ase=atoms)

    # Create 2x2x1 supercell
    # Pass ASE Atoms (simulating how AiiDA unpacks StructureData)
    # Function returns ASE Atoms (pythonjob converts to StructureData)
    supercell_ase = create_supercell._callable(structure.get_ase(), [2, 2, 1])

    # Check result is ASE Atoms
    assert isinstance(supercell_ase, Atoms)
    assert len(supercell_ase) == len(atoms) * 4  # 2*2*1

    # Convert to StructureData for cell comparison
    supercell = orm.StructureData(ase=supercell_ase)
    original_cell = structure.cell
    super_cell = supercell.cell
    assert abs(super_cell[0][0] - 2 * original_cell[0][0]) < 0.01
    assert abs(super_cell[1][1] - 2 * original_cell[1][1]) < 0.01
    assert abs(super_cell[2][2] - original_cell[2][2]) < 0.01


def test_create_supercell_3x3x2():
    """Test larger supercell."""
    from ase import Atoms
    atoms = bulk('Fe', 'bcc', a=2.87)
    structure = orm.StructureData(ase=atoms)

    supercell_ase = create_supercell._callable(structure.get_ase(), [3, 3, 2])

    assert isinstance(supercell_ase, Atoms)
    assert len(supercell_ase) == len(atoms) * 18  # 3*3*2
