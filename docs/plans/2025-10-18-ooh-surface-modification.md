# OOH Surface Modification Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Create a script to modify LaMnO3 surface structure by adding an OOH radical to the most exposed O atom coordinated with the topmost surface Mn atom.

**Architecture:** Hybrid approach using Pymatgen for structure analysis and coordination detection, ASE for geometry manipulation and structure output. The script loads a CIF slab, identifies the surface Mn (highest z), finds its coordinated O atoms, selects the most exposed O (highest z), and adds O and H atoms to form OOH radical with proper bond lengths.

**Tech Stack:** Pymatgen (structure analysis), ASE (geometry manipulation), Python 3.x

---

## Task 1: Create structure loader and validator

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/structure_loader.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_structure_loader.py`

**Step 1: Write the failing test**

Create the test file:

```python
# test_structure_loader.py
import pytest
from pathlib import Path
from structure_loader import load_structure, validate_structure


def test_load_structure_from_cif():
    """Test loading CIF file returns valid Structure object"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    assert structure is not None
    assert len(structure) > 0  # Has atoms
    assert structure.lattice is not None


def test_validate_structure_has_required_elements():
    """Test validation checks for Mn and O elements"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    result = validate_structure(structure)

    assert result["valid"] is True
    assert "Mn" in result["elements"]
    assert "O" in result["elements"]
    assert result["num_mn"] > 0
    assert result["num_o"] > 0


def test_load_structure_invalid_path():
    """Test loading non-existent file raises error"""
    with pytest.raises(FileNotFoundError):
        load_structure("nonexistent.cif")
```

**Step 2: Run test to verify it fails**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy/teros/experimental/adsorption_energy/structures/lamno3
source ~/envs/aiida/bin/activate && python -m pytest test_structure_loader.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'structure_loader'"

**Step 3: Write minimal implementation**

Create the implementation file:

```python
# structure_loader.py
from pathlib import Path
from pymatgen.core import Structure


def load_structure(cif_path):
    """
    Load a structure from a CIF file.

    Args:
        cif_path: Path to CIF file (str or Path)

    Returns:
        pymatgen.core.Structure object

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    cif_path = Path(cif_path)

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    structure = Structure.from_file(str(cif_path))
    return structure


def validate_structure(structure):
    """
    Validate that structure has required elements.

    Args:
        structure: pymatgen.core.Structure object

    Returns:
        dict with keys: valid (bool), elements (list), num_mn (int), num_o (int)
    """
    elements = [str(site.specie) for site in structure]
    element_set = set(elements)

    num_mn = elements.count("Mn")
    num_o = elements.count("O")

    valid = "Mn" in element_set and "O" in element_set

    return {
        "valid": valid,
        "elements": list(element_set),
        "num_mn": num_mn,
        "num_o": num_o
    }
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_structure_loader.py -v
```

Expected: PASS (all 3 tests pass)

**Step 5: Commit**

```bash
git add test_structure_loader.py structure_loader.py
git commit -m "feat: add structure loader with validation for LaMnO3 surface"
```

---

## Task 2: Identify surface Mn atom

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/surface_detector.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_surface_detector.py`

**Step 1: Write the failing test**

```python
# test_surface_detector.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn


def test_find_surface_mn_returns_topmost():
    """Test finding Mn with highest z-coordinate"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    mn_index, mn_position = find_surface_mn(structure)

    assert mn_index is not None
    assert mn_index >= 0
    assert len(mn_position) == 3  # x, y, z coordinates

    # Verify it's actually a Mn atom
    assert str(structure[mn_index].specie) == "Mn"

    # Verify it's the highest z among all Mn atoms
    mn_sites = [i for i, site in enumerate(structure) if str(site.specie) == "Mn"]
    for idx in mn_sites:
        assert structure[mn_index].coords[2] >= structure[idx].coords[2]


def test_find_surface_mn_returns_cartesian_coords():
    """Test that returned position is in Cartesian coordinates"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    mn_index, mn_position = find_surface_mn(structure)

    # Cartesian coordinates should match structure's cartesian coords
    import numpy as np
    expected_coords = structure[mn_index].coords
    assert np.allclose(mn_position, expected_coords, atol=1e-6)
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_surface_detector.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'surface_detector'"

**Step 3: Write minimal implementation**

```python
# surface_detector.py
import numpy as np


def find_surface_mn(structure):
    """
    Find the topmost Mn atom in the structure (highest z-coordinate).

    Args:
        structure: pymatgen.core.Structure object

    Returns:
        tuple: (mn_index, mn_position)
            mn_index: int, index of surface Mn in structure
            mn_position: np.array, Cartesian coordinates [x, y, z]
    """
    # Find all Mn atoms
    mn_sites = [(i, site) for i, site in enumerate(structure) if str(site.specie) == "Mn"]

    if not mn_sites:
        raise ValueError("No Mn atoms found in structure")

    # Find Mn with highest z-coordinate (Cartesian)
    surface_mn_idx = None
    max_z = -np.inf

    for idx, site in mn_sites:
        z_coord = site.coords[2]  # Cartesian z
        if z_coord > max_z:
            max_z = z_coord
            surface_mn_idx = idx

    surface_mn_position = structure[surface_mn_idx].coords

    return surface_mn_idx, surface_mn_position
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_surface_detector.py -v
```

Expected: PASS (all 2 tests pass)

**Step 5: Commit**

```bash
git add test_surface_detector.py surface_detector.py
git commit -m "feat: add surface Mn detector for topmost Mn identification"
```

---

## Task 3: Find coordinated O atoms

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/coordination_finder.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_coordination_finder.py`

**Step 1: Write the failing test**

```python
# test_coordination_finder.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens


def test_find_coordinated_oxygens_returns_list():
    """Test finding O atoms coordinated to surface Mn"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    assert isinstance(o_indices, list)
    assert len(o_indices) > 0  # Should have at least some coordinated O

    # Verify all are actually O atoms
    for idx in o_indices:
        assert str(structure[idx].specie) == "O"


def test_coordinated_oxygens_within_distance():
    """Test that all coordinated O are within reasonable Mn-O bond distance"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)

    import numpy as np
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        assert distance <= 2.5  # Within cutoff
        assert distance > 1.5  # Not unreasonably close


def test_coordinated_oxygens_sorted_by_distance():
    """Test that returned O atoms are sorted by distance to Mn"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    import numpy as np
    distances = []
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        distances.append(distance)

    # Check sorted
    assert distances == sorted(distances)
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_coordination_finder.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'coordination_finder'"

**Step 3: Write minimal implementation**

```python
# coordination_finder.py
import numpy as np


def find_coordinated_oxygens(structure, mn_index, cutoff=2.5):
    """
    Find O atoms coordinated to a specific Mn atom.

    Args:
        structure: pymatgen.core.Structure object
        mn_index: int, index of Mn atom in structure
        cutoff: float, maximum Mn-O distance in Angstroms (default 2.5)

    Returns:
        list: indices of O atoms coordinated to Mn, sorted by distance
    """
    mn_position = structure[mn_index].coords

    # Find all O atoms within cutoff distance
    coordinated_o = []

    for i, site in enumerate(structure):
        if str(site.specie) == "O":
            o_position = site.coords
            distance = np.linalg.norm(mn_position - o_position)

            if distance <= cutoff:
                coordinated_o.append((i, distance))

    # Sort by distance
    coordinated_o.sort(key=lambda x: x[1])

    # Return only indices
    return [idx for idx, dist in coordinated_o]
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_coordination_finder.py -v
```

Expected: PASS (all 3 tests pass)

**Step 5: Commit**

```bash
git add test_coordination_finder.py coordination_finder.py
git commit -m "feat: add coordination finder for Mn-O bonds"
```

---

## Task 4: Select most exposed O atom

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/oxygen_selector.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_oxygen_selector.py`

**Step 1: Write the failing test**

```python
# test_oxygen_selector.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen


def test_select_most_exposed_oxygen_returns_highest_z():
    """Test selecting O with highest z-coordinate from coordinated O atoms"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)

    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    assert exposed_o_idx in o_indices
    assert len(exposed_o_pos) == 3

    # Verify it's the highest z among coordinated O
    for idx in o_indices:
        assert exposed_o_pos[2] >= structure[idx].coords[2]


def test_select_most_exposed_oxygen_is_actually_oxygen():
    """Test that selected atom is actually O"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)

    exposed_o_idx, _ = select_most_exposed_oxygen(structure, o_indices)

    assert str(structure[exposed_o_idx].specie) == "O"


def test_select_most_exposed_oxygen_empty_list_raises_error():
    """Test that empty O list raises appropriate error"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    with pytest.raises(ValueError, match="No oxygen atoms provided"):
        select_most_exposed_oxygen(structure, [])
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_oxygen_selector.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'oxygen_selector'"

**Step 3: Write minimal implementation**

```python
# oxygen_selector.py
import numpy as np


def select_most_exposed_oxygen(structure, o_indices):
    """
    Select the most exposed O atom (highest z-coordinate) from a list of O indices.

    Args:
        structure: pymatgen.core.Structure object
        o_indices: list of int, indices of O atoms to consider

    Returns:
        tuple: (exposed_o_index, exposed_o_position)
            exposed_o_index: int, index of most exposed O
            exposed_o_position: np.array, Cartesian coordinates [x, y, z]

    Raises:
        ValueError: If o_indices is empty
    """
    if not o_indices:
        raise ValueError("No oxygen atoms provided")

    # Find O with highest z-coordinate
    max_z = -np.inf
    exposed_o_idx = None

    for idx in o_indices:
        z_coord = structure[idx].coords[2]
        if z_coord > max_z:
            max_z = z_coord
            exposed_o_idx = idx

    exposed_o_position = structure[exposed_o_idx].coords

    return exposed_o_idx, exposed_o_position
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_oxygen_selector.py -v
```

Expected: PASS (all 3 tests pass)

**Step 5: Commit**

```bash
git add test_oxygen_selector.py oxygen_selector.py
git commit -m "feat: add oxygen selector for most exposed O atom"
```

---

## Task 5: Construct OOH radical geometry

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/ooh_constructor.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_ooh_constructor.py`

**Step 1: Write the failing test**

```python
# test_ooh_constructor.py
import pytest
import numpy as np
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical


def test_construct_ooh_adds_two_atoms():
    """Test that OOH construction adds exactly 2 atoms (O and H)"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    original_num_atoms = len(structure)

    # Convert to ASE for modification
    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    assert len(modified_atoms) == original_num_atoms + 2


def test_construct_ooh_adds_correct_species():
    """Test that added atoms are O and H"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    # Check that new atoms are O and H
    new_atom_symbols = modified_atoms.get_chemical_symbols()[original_num_atoms:]
    assert 'O' in new_atom_symbols
    assert 'H' in new_atom_symbols


def test_construct_ooh_bond_lengths():
    """Test that O-O and O-H bond lengths are reasonable"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    positions = modified_atoms.get_positions()

    # Get positions of base O, new O, and H
    base_o_pos = exposed_o_pos
    new_o_pos = positions[original_num_atoms]  # First new atom (O)
    h_pos = positions[original_num_atoms + 1]  # Second new atom (H)

    # Check O-O bond length (~1.45 Å for OOH)
    oo_distance = np.linalg.norm(new_o_pos - base_o_pos)
    assert 1.3 < oo_distance < 1.6, f"O-O bond length {oo_distance:.2f} Å out of range"

    # Check O-H bond length (~0.97 Å for O-H)
    oh_distance = np.linalg.norm(h_pos - new_o_pos)
    assert 0.85 < oh_distance < 1.1, f"O-H bond length {oh_distance:.2f} Å out of range"


def test_construct_ooh_orientation():
    """Test that OOH points away from surface (positive z direction)"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    positions = modified_atoms.get_positions()
    base_o_z = exposed_o_pos[2]
    new_o_z = positions[original_num_atoms][2]
    h_z = positions[original_num_atoms + 1][2]

    # New O and H should be above base O
    assert new_o_z > base_o_z
    assert h_z > new_o_z
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_ooh_constructor.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'ooh_constructor'"

**Step 3: Write minimal implementation**

```python
# ooh_constructor.py
import numpy as np
from ase import Atoms


def construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos,
                         oo_bond_length=1.45, oh_bond_length=0.97):
    """
    Add O and H atoms to form OOH radical on the exposed O atom.

    Args:
        ase_atoms: ASE Atoms object (will be modified)
        exposed_o_idx: int, index of base O atom
        exposed_o_pos: np.array, position of base O atom [x, y, z]
        oo_bond_length: float, O-O bond length in Angstroms (default 1.45)
        oh_bond_length: float, O-H bond length in Angstroms (default 0.97)

    Returns:
        ASE Atoms object with OOH radical added
    """
    # Create a copy to avoid modifying original
    modified_atoms = ase_atoms.copy()

    # Calculate positions for new O and H atoms
    # Place them along positive z direction (away from surface)

    # New O atom position: base_O + (0, 0, oo_bond_length)
    new_o_position = exposed_o_pos + np.array([0.0, 0.0, oo_bond_length])

    # H atom position: new_O + (0, 0, oh_bond_length)
    h_position = new_o_position + np.array([0.0, 0.0, oh_bond_length])

    # Add new O atom
    from ase import Atom
    modified_atoms.append(Atom('O', position=new_o_position))

    # Add H atom
    modified_atoms.append(Atom('H', position=h_position))

    return modified_atoms
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_ooh_constructor.py -v
```

Expected: PASS (all 4 tests pass)

**Step 5: Commit**

```bash
git add test_ooh_constructor.py ooh_constructor.py
git commit -m "feat: add OOH radical constructor with proper bond geometry"
```

---

## Task 6: Export modified structure

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/structure_exporter.py`
- Create test: `teros/experimental/adsorption_energy/structures/lamno3/test_structure_exporter.py`

**Step 1: Write the failing test**

```python
# test_structure_exporter.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical
from structure_exporter import export_structure_to_cif
from pymatgen.io.ase import AseAtomsAdaptor


def test_export_structure_creates_file():
    """Test that export creates a CIF file"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"

    # Clean up if exists
    if output_path.exists():
        output_path.unlink()

    export_structure_to_cif(modified_atoms, output_path)

    assert output_path.exists()

    # Clean up
    output_path.unlink()


def test_export_structure_readable():
    """Test that exported CIF can be loaded back"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"
    export_structure_to_cif(modified_atoms, output_path)

    # Load back and verify
    loaded_structure = load_structure(output_path)
    assert len(loaded_structure) == original_num_atoms + 2

    # Clean up
    output_path.unlink()


def test_export_structure_contains_h():
    """Test that exported structure contains H atom"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"
    export_structure_to_cif(modified_atoms, output_path)

    loaded_structure = load_structure(output_path)
    elements = [str(site.specie) for site in loaded_structure]
    assert 'H' in elements

    # Clean up
    output_path.unlink()
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_structure_exporter.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'structure_exporter'"

**Step 3: Write minimal implementation**

```python
# structure_exporter.py
from pathlib import Path
from ase.io import write


def export_structure_to_cif(ase_atoms, output_path):
    """
    Export ASE Atoms object to CIF file.

    Args:
        ase_atoms: ASE Atoms object
        output_path: Path to output CIF file (str or Path)
    """
    output_path = Path(output_path)

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write CIF file
    write(str(output_path), ase_atoms, format='cif')
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_structure_exporter.py -v
```

Expected: PASS (all 3 tests pass)

**Step 5: Commit**

```bash
git add test_structure_exporter.py structure_exporter.py
git commit -m "feat: add structure exporter for CIF output"
```

---

## Task 7: Create main script with CLI

**Files:**
- Create: `teros/experimental/adsorption_energy/structures/lamno3/add_ooh_to_surface.py`
- Test: Run the main script on actual LaMnO3_100_A4_surface.cif

**Step 1: Write integration script**

```python
# add_ooh_to_surface.py
#!/usr/bin/env python
"""
Add OOH radical to LaMnO3 surface structure.

This script:
1. Loads a surface slab structure (CIF)
2. Identifies the topmost Mn atom (surface Mn)
3. Finds O atoms coordinated to this Mn
4. Selects the most exposed O (highest z)
5. Adds O and H to form OOH radical
6. Exports modified structure to new CIF file
"""

import sys
from pathlib import Path
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from structure_loader import load_structure, validate_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical
from structure_exporter import export_structure_to_cif


def main():
    """Main execution function"""

    # Configuration
    input_cif = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    output_cif = Path(__file__).parent / "LaMnO3_100_A4_surface_OOH.cif"

    print("=" * 60)
    print("OOH Surface Modification Script")
    print("=" * 60)
    print()

    # Step 1: Load and validate structure
    print(f"Loading structure from: {input_cif}")
    structure = load_structure(input_cif)

    validation = validate_structure(structure)
    if not validation["valid"]:
        print("ERROR: Structure validation failed")
        print(f"  Missing required elements")
        sys.exit(1)

    print(f"✓ Structure loaded successfully")
    print(f"  Total atoms: {len(structure)}")
    print(f"  Elements: {', '.join(validation['elements'])}")
    print(f"  Mn atoms: {validation['num_mn']}")
    print(f"  O atoms: {validation['num_o']}")
    print()

    # Step 2: Find surface Mn
    print("Finding surface Mn atom...")
    mn_index, mn_position = find_surface_mn(structure)
    print(f"✓ Surface Mn identified")
    print(f"  Index: {mn_index}")
    print(f"  Position (Cartesian): [{mn_position[0]:.4f}, {mn_position[1]:.4f}, {mn_position[2]:.4f}]")
    print()

    # Step 3: Find coordinated O atoms
    print("Finding coordinated O atoms...")
    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)
    print(f"✓ Found {len(o_indices)} coordinated O atoms")
    for i, o_idx in enumerate(o_indices, 1):
        o_pos = structure[o_idx].coords
        distance = np.linalg.norm(mn_position - o_pos)
        print(f"  O{i}: index={o_idx}, distance={distance:.3f} Å")
    print()

    # Step 4: Select most exposed O
    print("Selecting most exposed O atom...")
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)
    print(f"✓ Most exposed O selected")
    print(f"  Index: {exposed_o_idx}")
    print(f"  Position (Cartesian): [{exposed_o_pos[0]:.4f}, {exposed_o_pos[1]:.4f}, {exposed_o_pos[2]:.4f}]")
    print(f"  Z-coordinate: {exposed_o_pos[2]:.4f} Å")
    print()

    # Step 5: Convert to ASE and construct OOH
    print("Constructing OOH radical...")
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    print(f"✓ OOH radical constructed")
    print(f"  Original atoms: {len(ase_atoms)}")
    print(f"  Modified atoms: {len(modified_atoms)}")
    print(f"  Added: 1 O + 1 H")

    # Get positions of new atoms
    positions = modified_atoms.get_positions()
    new_o_pos = positions[len(ase_atoms)]
    h_pos = positions[len(ase_atoms) + 1]

    oo_distance = np.linalg.norm(new_o_pos - exposed_o_pos)
    oh_distance = np.linalg.norm(h_pos - new_o_pos)

    print(f"  O-O bond length: {oo_distance:.3f} Å")
    print(f"  O-H bond length: {oh_distance:.3f} Å")
    print()

    # Step 6: Export modified structure
    print(f"Exporting modified structure to: {output_cif}")
    export_structure_to_cif(modified_atoms, output_cif)
    print(f"✓ Structure exported successfully")
    print()

    print("=" * 60)
    print("OOH modification completed successfully!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2: Run the integration script**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy/teros/experimental/adsorption_energy/structures/lamno3
source ~/envs/aiida/bin/activate && python add_ooh_to_surface.py
```

Expected output:
```
============================================================
OOH Surface Modification Script
============================================================

Loading structure from: LaMnO3_100_A4_surface.cif
✓ Structure loaded successfully
  Total atoms: 37
  Elements: La, Mn, O
  Mn atoms: 7
  O atoms: 22

Finding surface Mn atom...
✓ Surface Mn identified
  Index: [number]
  Position (Cartesian): [x, y, z]

Finding coordinated O atoms...
✓ Found [n] coordinated O atoms
  O1: index=[num], distance=[dist] Å
  ...

Selecting most exposed O atom...
✓ Most exposed O selected
  Index: [number]
  Position (Cartesian): [x, y, z]

Constructing OOH radical...
✓ OOH radical constructed
  Original atoms: 37
  Modified atoms: 39
  Added: 1 O + 1 H
  O-O bond length: 1.450 Å
  O-H bond length: 0.970 Å

Exporting modified structure to: LaMnO3_100_A4_surface_OOH.cif
✓ Structure exported successfully

============================================================
OOH modification completed successfully!
============================================================
```

**Step 3: Verify output file exists and is valid**

Run:
```bash
ls -lh LaMnO3_100_A4_surface_OOH.cif
source ~/envs/aiida/bin/activate && python -c "from structure_loader import load_structure; s = load_structure('LaMnO3_100_A4_surface_OOH.cif'); print(f'Atoms: {len(s)}'); print(f'Elements: {set([str(site.specie) for site in s])}')"
```

Expected:
- File exists with reasonable size
- Contains 39 atoms (37 original + 2 new)
- Elements include H

**Step 4: Run all tests to ensure integration works**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest test_*.py -v
```

Expected: ALL tests pass

**Step 5: Commit main script**

```bash
git add add_ooh_to_surface.py
git commit -m "feat: add main script for OOH surface modification with full workflow"
```

---

## Task 8: Clean up and final verification

**Step 1: Run full test suite**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy/teros/experimental/adsorption_energy/structures/lamno3
source ~/envs/aiida/bin/activate && python -m pytest test_*.py -v --tb=short
```

Expected: All tests pass

**Step 2: Run main script again to generate final output**

Run:
```bash
source ~/envs/aiida/bin/activate && python add_ooh_to_surface.py
```

Expected: Successful execution with output CIF file

**Step 3: Visual inspection (optional)**

If ASE GUI available:
```bash
source ~/envs/aiida/bin/activate && ase gui LaMnO3_100_A4_surface_OOH.cif
```

**Step 4: Create documentation**

Create `README.md` in the same directory documenting usage:

```markdown
# OOH Surface Modification Tool

## Overview

This tool modifies a LaMnO3 surface slab structure by adding an OOH radical to the most exposed oxygen atom coordinated with the topmost surface Mn atom.

## Usage

```bash
python add_ooh_to_surface.py
```

## Input

- `LaMnO3_100_A4_surface.cif`: Input surface slab structure

## Output

- `LaMnO3_100_A4_surface_OOH.cif`: Modified structure with OOH radical

## Methodology

1. Identify topmost Mn atom (highest z-coordinate)
2. Find O atoms coordinated to this Mn (within 2.5 Å)
3. Select most exposed O (highest z-coordinate)
4. Add OOH radical:
   - O-O bond: 1.45 Å
   - O-H bond: 0.97 Å
   - Orientation: perpendicular to surface (+z direction)

## Dependencies

- pymatgen
- ASE
- numpy

## Testing

Run all tests:
```bash
pytest test_*.py -v
```
```

**Step 5: Final commit**

```bash
git add README.md
git commit -m "docs: add README for OOH surface modification tool"
```

---

## Completion Checklist

- [ ] All 8 tasks completed
- [ ] All tests passing
- [ ] Main script executes successfully
- [ ] Output CIF file generated and validated
- [ ] Documentation created
- [ ] Code committed with descriptive messages
- [ ] No debug print statements left in code
- [ ] All imports are used
- [ ] File paths are absolute and correct

---

## Notes

- Bond lengths used (typical for OOH radical):
  - O-O: 1.45 Å
  - O-H: 0.97 Å

- Mn-O coordination cutoff: 2.5 Å (typical for perovskites)

- Orientation: OOH placed perpendicular to surface (+z direction) for simplicity. For more sophisticated placement, could use surface normal or Mn-O bond direction.

- All tests use the actual `LaMnO3_100_A4_surface.cif` file in the directory for integration testing.
