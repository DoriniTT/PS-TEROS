# Adsorption Energy Module Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Create a core module for calculating adsorption energies from substrate+adsorbate structures using AiiDA-WorkGraph and VASP.

**Architecture:** Module follows the thermodynamics/cleavage pattern with calcfunctions for structure separation and energy calculation, plus task.graph scatter-gather for parallel VASP relaxations. Uses pymatgen StructureGraph for robust adsorbate identification via connectivity analysis.

**Tech Stack:** AiiDA-WorkGraph, pymatgen, ASE, VASP (via aiida-vasp)

---

## Task 1: Create adsorbate separation calcfunction with connectivity analysis

**Files:**
- Create: `teros/core/adsorption_energy.py`
- Create test: `teros/core/test_adsorption_energy.py`

**Step 1: Write the failing test**

Create the test file:

```python
# teros/core/test_adsorption_energy.py
import pytest
import numpy as np
from aiida import orm
from pymatgen.core import Structure, Lattice
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from adsorption_energy import separate_adsorbate_structure


def create_test_surface_with_oh():
    """Create a simple Ag(111) surface with OH adsorbate for testing."""
    # Create 2x2 Ag slab
    lattice = Lattice.from_parameters(a=5.8, b=5.8, c=20.0,
                                     alpha=90, beta=90, gamma=90)

    # Ag atoms in bottom layer
    ag_positions = [
        [0.0, 0.0, 10.0],
        [2.9, 0.0, 10.0],
        [0.0, 2.9, 10.0],
        [2.9, 2.9, 10.0],
    ]

    # OH adsorbate on top (bonded O-H cluster)
    oh_positions = [
        [1.45, 1.45, 12.5],  # O atom above surface
        [1.45, 1.45, 13.5],  # H atom above O
    ]

    species = ['Ag'] * 4 + ['O', 'H']
    positions = ag_positions + oh_positions

    structure = Structure(lattice, species, positions, coords_are_cartesian=True)

    return orm.StructureData(pymatgen=structure)


def test_separate_adsorbate_returns_three_structures():
    """Test that separation returns substrate, molecule, and complete structures."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = separate_adsorbate_structure(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    assert 'substrate' in result
    assert 'molecule' in result
    assert 'complete' in result

    assert isinstance(result['substrate'], orm.StructureData)
    assert isinstance(result['molecule'], orm.StructureData)
    assert isinstance(result['complete'], orm.StructureData)


def test_separate_adsorbate_correct_atom_counts():
    """Test that separated structures have correct number of atoms."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = separate_adsorbate_structure(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    # Original: 4 Ag + 2 (O,H) = 6 atoms
    # Substrate: 4 Ag
    # Molecule: 1 O + 1 H = 2 atoms
    # Complete: 6 atoms

    substrate_ase = result['substrate'].get_ase()
    molecule_ase = result['molecule'].get_ase()
    complete_ase = result['complete'].get_ase()

    assert len(substrate_ase) == 4
    assert len(molecule_ase) == 2
    assert len(complete_ase) == 6


def test_separate_adsorbate_correct_species():
    """Test that molecule contains only adsorbate species."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = separate_adsorbate_structure(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    molecule_ase = result['molecule'].get_ase()
    symbols = molecule_ase.get_chemical_symbols()

    assert 'O' in symbols
    assert 'H' in symbols
    assert 'Ag' not in symbols


def test_separate_adsorbate_same_cell():
    """Test that all three structures have the same cell."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = separate_adsorbate_structure(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    complete_ase = result['complete'].get_ase()
    substrate_ase = result['substrate'].get_ase()
    molecule_ase = result['molecule'].get_ase()

    complete_cell = complete_ase.get_cell()
    substrate_cell = substrate_ase.get_cell()
    molecule_cell = molecule_ase.get_cell()

    assert np.allclose(complete_cell, substrate_cell)
    assert np.allclose(complete_cell, molecule_cell)


def test_separate_adsorbate_invalid_formula():
    """Test that invalid formula raises error."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('')  # Empty formula

    with pytest.raises(ValueError, match="Adsorbate formula cannot be empty"):
        separate_adsorbate_structure(
            structure=complete,
            adsorbate_formula=adsorbate_formula
        )


def test_separate_adsorbate_no_match():
    """Test that non-existent adsorbate raises error."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OOH')  # Not in structure

    with pytest.raises(ValueError, match="Could not find adsorbate"):
        separate_adsorbate_structure(
            structure=complete,
            adsorbate_formula=adsorbate_formula
        )
```

**Step 2: Run test to verify it fails**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'adsorption_energy'"

**Step 3: Write minimal implementation**

Create the implementation file:

```python
# teros/core/adsorption_energy.py
"""
Adsorption Energy Calculations Module

This module provides functions to calculate adsorption energies from
substrate+adsorbate structures using AiiDA-WorkGraph and VASP.
"""

from __future__ import annotations

import typing as t
from collections import Counter

import numpy as np
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Composition


def parse_formula(formula_str: str) -> dict[str, int]:
    """
    Parse chemical formula string into element counts.

    Args:
        formula_str: Chemical formula (e.g., "OOH", "H2O")

    Returns:
        Dictionary mapping element symbols to counts

    Example:
        >>> parse_formula("OOH")
        {'O': 2, 'H': 1}
    """
    comp = Composition(formula_str)
    return {str(el): int(count) for el, count in comp.items()}


@task.calcfunction
def separate_adsorbate_structure(
    structure: orm.StructureData,
    adsorbate_formula: orm.Str,
) -> t.Annotated[dict, namespace(
    substrate=orm.StructureData,
    molecule=orm.StructureData,
    complete=orm.StructureData
)]:
    """
    Separate substrate+adsorbate structure into three components.

    Uses connectivity analysis (pymatgen StructureGraph) to identify bonded
    clusters and match the adsorbate by chemical formula.

    Args:
        structure: Complete structure (substrate + adsorbate)
        adsorbate_formula: Chemical formula of adsorbate (e.g., "OOH", "OH")

    Returns:
        Dictionary with three StructureData nodes:
        - substrate: Structure with adsorbate removed, original cell preserved
        - molecule: Adsorbate in original cell at original position
        - complete: Original structure (for provenance)

    Raises:
        ValueError: If formula is invalid, adsorbate not found, or multiple matches
    """
    from pymatgen.io.ase import AseAtomsAdaptor

    # Validation: check formula is not empty
    formula_str = adsorbate_formula.value.strip()
    if not formula_str:
        raise ValueError("Adsorbate formula cannot be empty")

    # Parse adsorbate formula
    try:
        target_composition = Counter(parse_formula(formula_str))
    except Exception as e:
        raise ValueError(f"Invalid chemical formula '{formula_str}': {e}")

    # Get pymatgen structure
    pmg_structure = structure.get_pymatgen()

    # Validate structure has enough atoms
    total_target_atoms = sum(target_composition.values())
    if len(pmg_structure) < total_target_atoms:
        raise ValueError(
            f"Structure has only {len(pmg_structure)} atoms, "
            f"cannot contain adsorbate {formula_str} ({total_target_atoms} atoms)"
        )

    # Build connectivity graph using CrystalNN
    # CrystalNN works for both crystals and molecules
    sg = StructureGraph.with_local_env_strategy(pmg_structure, CrystalNN())

    # Get connected components (bonded clusters)
    connected_components = sg.get_subgraphs_as_molecules()

    # Find cluster matching adsorbate formula
    matching_clusters = []
    for i, molecule in enumerate(connected_components):
        # Get composition of this cluster
        cluster_composition = Counter(
            {str(el): int(count) for el, count in molecule.composition.items()}
        )

        if cluster_composition == target_composition:
            matching_clusters.append((i, molecule))

    # Validation: check we found exactly one match
    if len(matching_clusters) == 0:
        available_formulas = [
            mol.composition.reduced_formula
            for mol in connected_components
        ]
        raise ValueError(
            f"Could not find adsorbate '{formula_str}' in structure. "
            f"Found connected clusters with formulas: {available_formulas}"
        )

    if len(matching_clusters) > 1:
        raise ValueError(
            f"Found {len(matching_clusters)} clusters matching '{formula_str}'. "
            f"Structure is ambiguous - please provide structure with unique adsorbate."
        )

    # Get the matched adsorbate
    cluster_idx, adsorbate_molecule = matching_clusters[0]

    # Get indices of adsorbate atoms in original structure
    # The molecule object contains the sites, we need to map back to indices
    adsorbate_indices = set()

    for mol_site in adsorbate_molecule:
        # Find matching site in original structure
        for i, struct_site in enumerate(pmg_structure):
            # Check if same element and same position
            if (str(struct_site.specie) == str(mol_site.specie) and
                np.allclose(struct_site.coords, mol_site.coords, atol=1e-3)):
                adsorbate_indices.add(i)
                break

    # Validation: check substrate is not too small
    substrate_indices = [i for i in range(len(pmg_structure)) if i not in adsorbate_indices]
    if len(substrate_indices) < 3:
        raise ValueError(
            f"Substrate only has {len(substrate_indices)} atoms after removing adsorbate. "
            f"Structure may be invalid."
        )

    # Create three structures using ASE for manipulation
    adaptor = AseAtomsAdaptor()
    ase_structure = adaptor.get_atoms(pmg_structure)

    # 1. Substrate: remove adsorbate atoms, keep cell
    substrate_ase = ase_structure.copy()
    del substrate_ase[[i for i in adsorbate_indices]]

    # 2. Molecule: keep only adsorbate atoms, keep cell and position
    molecule_ase = ase_structure.copy()
    del molecule_ase[substrate_indices]

    # 3. Complete: return as-is for provenance
    complete_ase = ase_structure.copy()

    # Convert back to AiiDA StructureData
    return {
        'substrate': orm.StructureData(ase=substrate_ase),
        'molecule': orm.StructureData(ase=molecule_ase),
        'complete': orm.StructureData(ase=complete_ase),
    }
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_returns_three_structures -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_correct_atom_counts -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_correct_species -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_same_cell -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_invalid_formula -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_separate_adsorbate_no_match -v
```

Expected: All 6 tests PASS

**Step 5: Commit**

```bash
git add teros/core/adsorption_energy.py teros/core/test_adsorption_energy.py
git commit -m "feat: add adsorbate separation with connectivity analysis"
```

---

## Task 2: Create adsorption energy calculation calcfunction

**Files:**
- Modify: `teros/core/adsorption_energy.py`
- Modify: `teros/core/test_adsorption_energy.py`

**Step 1: Write the failing test**

Add to test file:

```python
# Add to teros/core/test_adsorption_energy.py

from adsorption_energy import calculate_adsorption_energy


def test_calculate_adsorption_energy_correct_formula():
    """Test that adsorption energy uses correct formula."""
    # E_ads = E_complete - E_substrate - E_molecule

    E_complete = orm.Float(-100.0)
    E_substrate = orm.Float(-80.0)
    E_molecule = orm.Float(-15.0)

    E_ads = calculate_adsorption_energy(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    # E_ads = -100 - (-80) - (-15) = -100 + 80 + 15 = -5.0
    expected = -5.0
    assert abs(E_ads.value - expected) < 1e-9


def test_calculate_adsorption_energy_positive():
    """Test case where adsorption is endothermic (positive E_ads)."""
    E_complete = orm.Float(-50.0)
    E_substrate = orm.Float(-40.0)
    E_molecule = orm.Float(-15.0)

    E_ads = calculate_adsorption_energy(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    # E_ads = -50 - (-40) - (-15) = -50 + 40 + 15 = 5.0 (endothermic)
    expected = 5.0
    assert abs(E_ads.value - expected) < 1e-9


def test_calculate_adsorption_energy_returns_float():
    """Test that result is Float node."""
    E_complete = orm.Float(-100.0)
    E_substrate = orm.Float(-80.0)
    E_molecule = orm.Float(-15.0)

    E_ads = calculate_adsorption_energy(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    assert isinstance(E_ads, orm.Float)
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_calculate_adsorption_energy_correct_formula -v
```

Expected: FAIL with "ImportError: cannot import name 'calculate_adsorption_energy'"

**Step 3: Write minimal implementation**

Add to implementation file:

```python
# Add to teros/core/adsorption_energy.py


@task.calcfunction
def calculate_adsorption_energy(
    E_complete: orm.Float,
    E_substrate: orm.Float,
    E_molecule: orm.Float,
) -> orm.Float:
    """
    Calculate adsorption energy from component energies.

    Uses the formula:
        E_ads = E_complete - E_substrate - E_molecule

    A negative adsorption energy indicates exothermic (favorable) adsorption.
    A positive adsorption energy indicates endothermic (unfavorable) adsorption.

    Args:
        E_complete: Total energy of substrate+adsorbate system (eV)
        E_substrate: Total energy of bare substrate (eV)
        E_molecule: Total energy of isolated adsorbate molecule (eV)

    Returns:
        Adsorption energy in eV
    """
    E_ads = E_complete.value - E_substrate.value - E_molecule.value
    return orm.Float(E_ads)
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_calculate_adsorption_energy_correct_formula -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_calculate_adsorption_energy_positive -v
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_calculate_adsorption_energy_returns_float -v
```

Expected: All 3 tests PASS

**Step 5: Commit**

```bash
git add teros/core/adsorption_energy.py teros/core/test_adsorption_energy.py
git commit -m "feat: add adsorption energy calculation function"
```

---

## Task 3: Create scatter-gather workflow for parallel processing

**Files:**
- Modify: `teros/core/adsorption_energy.py`

**Step 1: Write the failing test**

Add to test file:

```python
# Add to teros/core/test_adsorption_energy.py

from adsorption_energy import compute_adsorption_energies_scatter


def test_compute_adsorption_energies_scatter_signature():
    """Test that scatter function exists with correct signature."""
    # This test just checks the function exists and accepts right parameters
    # Full integration test will be in examples folder

    import inspect

    sig = inspect.signature(compute_adsorption_energies_scatter)
    params = list(sig.parameters.keys())

    # Check required parameters exist
    assert 'structures' in params
    assert 'adsorbate_formulas' in params
    assert 'code' in params
    assert 'potential_family' in params
    assert 'potential_mapping' in params
    assert 'parameters' in params
    assert 'options' in params
```

**Step 2: Run test to verify it fails**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_compute_adsorption_energies_scatter_signature -v
```

Expected: FAIL with "ImportError: cannot import name 'compute_adsorption_energies_scatter'"

**Step 3: Write minimal implementation**

Add to implementation file:

```python
# Add to teros/core/adsorption_energy.py

# Import at top of file
from .slabs import extract_total_energy


@task.graph
def compute_adsorption_energies_scatter(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    adsorbate_formulas: t.Annotated[dict[str, str], dict],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
)]:
    """
    Scatter-gather workflow for calculating adsorption energies.

    This workflow:
    1. Separates all substrate+adsorbate structures in parallel
    2. Relaxes all systems (3N VASP jobs) in parallel
    3. Calculates adsorption energies for each system

    Args:
        structures: Dynamic namespace of complete structures
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'system1': 'OOH', 'system2': 'OH'}
        code: AiiDA code for VASP
        potential_family: Pseudopotential family name
        potential_mapping: Element to potential mapping
        parameters: VASP INCAR parameters
        options: Scheduler options for VASP calculations
        kpoints_spacing: K-points spacing in Angstrom^-1 (optional)
        clean_workdir: Whether to clean remote working directories

    Returns:
        Dictionary with namespaces:
        - separated_structures: Dict of separated systems for each input
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies of complete systems
        - adsorption_energies: Final adsorption energies
    """
    # Validate inputs
    if structures.keys() != adsorbate_formulas.keys():
        raise ValueError(
            f"Mismatch between structures and adsorbate_formulas keys. "
            f"Structures: {list(structures.keys())}, "
            f"Formulas: {list(adsorbate_formulas.keys())}"
        )

    # Get VASP workflow
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    vasp_task_cls = task(vasp_wc)

    # Output dictionaries
    separated_dict: dict[str, dict] = {}
    substrate_energies: dict[str, orm.Float] = {}
    molecule_energies: dict[str, orm.Float] = {}
    complete_energies: dict[str, orm.Float] = {}
    adsorption_energies: dict[str, orm.Float] = {}

    # Phase 1: Separate structures (parallel)
    for key, structure in structures.items():
        adsorbate_str = adsorbate_formulas[key]

        separated = separate_adsorbate_structure(
            structure=structure,
            adsorbate_formula=orm.Str(adsorbate_str)
        ).result

        separated_dict[key] = separated

    # Phase 2: VASP relaxations (parallel, 3N jobs)
    for key, separated in separated_dict.items():
        # Helper function to create VASP inputs
        def create_vasp_inputs(struct: orm.StructureData) -> dict:
            inputs: dict[str, t.Any] = {
                'structure': struct,
                'code': code,
                'parameters': {'incar': dict(parameters)},
                'options': dict(options),
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'clean_workdir': clean_workdir,
                'settings': orm.Dict(dict={
                    'parser_settings': {
                        'add_trajectory': True,
                        'add_structure': True,
                        'add_kpoints': True,
                    }
                }),
            }
            if kpoints_spacing is not None:
                inputs['kpoints_spacing'] = kpoints_spacing
            return inputs

        # Substrate relaxation
        substrate_calc = vasp_task_cls(**create_vasp_inputs(separated['substrate']))
        substrate_energies[key] = extract_total_energy(energies=substrate_calc.misc).result

        # Molecule relaxation
        molecule_calc = vasp_task_cls(**create_vasp_inputs(separated['molecule']))
        molecule_energies[key] = extract_total_energy(energies=molecule_calc.misc).result

        # Complete system relaxation
        complete_calc = vasp_task_cls(**create_vasp_inputs(separated['complete']))
        complete_energies[key] = extract_total_energy(energies=complete_calc.misc).result

    # Phase 3: Calculate adsorption energies (parallel)
    for key in structures.keys():
        E_ads = calculate_adsorption_energy(
            E_complete=complete_energies[key],
            E_substrate=substrate_energies[key],
            E_molecule=molecule_energies[key],
        ).result

        adsorption_energies[key] = E_ads

    # Return all results
    return {
        'separated_structures': separated_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': complete_energies,
        'adsorption_energies': adsorption_energies,
    }
```

**Step 4: Run test to verify it passes**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py::test_compute_adsorption_energies_scatter_signature -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add teros/core/adsorption_energy.py teros/core/test_adsorption_energy.py
git commit -m "feat: add scatter-gather workflow for parallel adsorption energy calculations"
```

---

## Task 4: Update core module exports

**Files:**
- Modify: `teros/core/__init__.py`

**Step 1: Add exports to __init__.py**

Add the new functions to the exports:

```python
# Modify teros/core/__init__.py

# Add to imports section:
from .adsorption_energy import (
    separate_adsorbate_structure,
    calculate_adsorption_energy,
    compute_adsorption_energies_scatter,
)

# Add to __all__ list:
__all__ = [
    'get_structure_from_file',
    'prepare_vasp_inputs',
    'calculate_formation_enthalpy',
    'generate_slab_structures',
    'relax_slabs_scatter',
    'extract_total_energy',
    'identify_oxide_type',
    'calculate_surface_energy_ternary',
    'calculate_surface_energy_binary',
    'compute_surface_energies_scatter',
    'calculate_cleavage_energy',
    'compute_cleavage_energies_scatter',
    'list_workflow_presets',
    'get_preset_config',
    'get_preset_summary',
    'WORKFLOW_PRESETS',
    'DEFAULT_PRESET',
    'separate_adsorbate_structure',
    'calculate_adsorption_energy',
    'compute_adsorption_energies_scatter',
]
```

**Step 2: Test imports work**

Run:
```bash
source ~/envs/aiida/bin/activate && python -c "from teros.core import separate_adsorbate_structure, calculate_adsorption_energy, compute_adsorption_energies_scatter; print('Imports successful')"
```

Expected: "Imports successful"

**Step 3: Commit**

```bash
git add teros/core/__init__.py
git commit -m "feat: export adsorption energy functions from core module"
```

---

## Task 5: Create example workflow

**Files:**
- Create: `examples/adsorption_energy/test_oh_ag111/run_adsorption_energy.py`
- Create: `examples/adsorption_energy/test_oh_ag111/README.md`

**Step 1: Create directory structure**

Run:
```bash
mkdir -p examples/adsorption_energy/test_oh_ag111/structures
```

**Step 2: Create example workflow script**

```python
# examples/adsorption_energy/test_oh_ag111/run_adsorption_energy.py
#!/usr/bin/env python
"""
Example: Calculate adsorption energy of OH on Ag(111) surface.

This example demonstrates the adsorption energy module workflow:
1. Load substrate+adsorbate structures
2. Separate into substrate, molecule, complete
3. Relax all three systems with VASP
4. Calculate adsorption energy
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from aiida_workgraph import WorkGraph

# Import PSTEROS core functions
from teros.core import (
    get_structure_from_file,
    compute_adsorption_energies_scatter,
)


def main():
    """Run adsorption energy calculation workflow."""
    load_profile()

    print("=" * 70)
    print("Adsorption Energy Calculation: OH on Ag(111)")
    print("=" * 70)
    print()

    # Configuration
    base_path = Path(__file__).parent

    # Load VASP code (adjust for your setup)
    code_label = 'vasp@localhost'  # Modify to match your VASP code
    try:
        code = orm.load_code(code_label)
    except:
        print(f"ERROR: Could not load VASP code '{code_label}'")
        print("Please create a VASP code in AiiDA first:")
        print(f"  verdi code create core.code.installed")
        sys.exit(1)

    # Load structures
    # For this example, we'll create a simple test structure programmatically
    # In real use, you'd load from CIF files using get_structure_from_file()

    from pymatgen.core import Structure, Lattice
    from ase import Atoms

    # Create simple Ag slab with OH
    lattice = Lattice.from_parameters(a=5.8, b=5.8, c=20.0,
                                     alpha=90, beta=90, gamma=90)

    ag_positions = [[0.0, 0.0, 10.0], [2.9, 0.0, 10.0],
                    [0.0, 2.9, 10.0], [2.9, 2.9, 10.0]]
    oh_positions = [[1.45, 1.45, 12.5], [1.45, 1.45, 13.5]]

    species = ['Ag'] * 4 + ['O', 'H']
    positions = ag_positions + oh_positions

    structure = Structure(lattice, species, positions, coords_are_cartesian=True)

    # Prepare input data
    structures = {
        'site1': orm.StructureData(pymatgen=structure),
    }

    adsorbate_formulas = {
        'site1': 'OH',
    }

    # VASP parameters (minimal for testing, adjust for production)
    parameters = {
        'ENCUT': 400,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'NSW': 50,
        'ISIF': 2,  # Relax ions only, keep cell fixed
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options
    options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 4,
        },
        'max_wallclock_seconds': 3600,
        'queue_name': 'debug',  # Adjust for your cluster
    }

    # Potential mapping
    potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
        'H': 'H',
    }

    print("Setting up WorkGraph...")

    # Create WorkGraph
    wg = WorkGraph('adsorption_energy_oh_ag111')

    # Add adsorption energy calculation task
    ads_task = wg.add_task(
        compute_adsorption_energies_scatter,
        name='compute_adsorption_energies',
        structures=structures,
        adsorbate_formulas=adsorbate_formulas,
        code=code,
        potential_family='PBE',  # Adjust for your setup
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    print("Submitting WorkGraph to AiiDA daemon...")
    wg.submit(wait=False)

    print()
    print("=" * 70)
    print(f"WorkGraph submitted! PK: {wg.pk}")
    print()
    print("Monitor progress with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")
    print()
    print("After completion, check results with:")
    print(f"  verdi process show {wg.pk}")
    print("=" * 70)

    return wg.pk


if __name__ == '__main__':
    pk = main()
    sys.exit(0)
```

**Step 3: Create README documentation**

```markdown
# Adsorption Energy Example: OH on Ag(111)

## Overview

This example demonstrates how to calculate adsorption energies using the PSTEROS adsorption energy module.

## Workflow

1. **Input**: Substrate+adsorbate structure (Ag slab with OH)
2. **Separation**: Automatically separates into:
   - Bare Ag substrate
   - Isolated OH molecule
   - Complete Ag+OH system
3. **Relaxation**: VASP relaxes all three structures in parallel
4. **Calculation**: Computes E_ads = E_complete - E_substrate - E_molecule

## Usage

### Prerequisites

1. Configure AiiDA VASP code:
   ```bash
   verdi code create core.code.installed
   ```

2. Set up pseudopotentials:
   ```bash
   verdi data core.upf uploadfamily --path /path/to/potentials --name PBE
   ```

### Run the workflow

```bash
source ~/envs/aiida/bin/activate
cd examples/adsorption_energy/test_oh_ag111
python run_adsorption_energy.py
```

### Monitor progress

```bash
verdi process list
verdi process show <PK>
verdi process report <PK>
```

## Expected Results

- **Substrate energy**: Total energy of bare Ag(111) slab
- **Molecule energy**: Total energy of isolated OH in same cell
- **Complete energy**: Total energy of Ag(111)+OH system
- **Adsorption energy**: E_ads (negative = favorable adsorption)

## Notes

- Adjust VASP parameters for production calculations
- Use appropriate k-point sampling for your system
- Check convergence with respect to slab thickness and vacuum
- This example uses minimal parameters for fast testing
```

**Step 4: Test the example can be imported**

Run:
```bash
source ~/envs/aiida/bin/activate && python -c "import sys; sys.path.insert(0, 'examples/adsorption_energy/test_oh_ag111'); import run_adsorption_energy; print('Example script loads successfully')"
```

Expected: "Example script loads successfully"

**Step 5: Commit**

```bash
git add examples/adsorption_energy/
git commit -m "docs: add example workflow for OH/Ag(111) adsorption energy"
```

---

## Task 6: Run full test suite

**Step 1: Run all unit tests**

Run:
```bash
source ~/envs/aiida/bin/activate && python -m pytest teros/core/test_adsorption_energy.py -v
```

Expected: All tests PASS

**Step 2: Clear Python cache**

Run:
```bash
find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true
```

**Step 3: Restart AiiDA daemon**

Run:
```bash
verdi daemon restart
```

Expected: Daemon restarts successfully

**Step 4: Verify imports work**

Run:
```bash
source ~/envs/aiida/bin/activate && python -c "
from teros.core import (
    separate_adsorbate_structure,
    calculate_adsorption_energy,
    compute_adsorption_energies_scatter
)
print('All adsorption energy functions imported successfully')
"
```

Expected: "All adsorption energy functions imported successfully"

---

## Task 7: Create integration documentation

**Files:**
- Create: `docs/adsorption_energy_module.md`

**Step 1: Create module documentation**

```markdown
# Adsorption Energy Module

## Overview

The adsorption energy module provides tools to calculate adsorption energies from substrate+adsorbate structures using AiiDA-WorkGraph and VASP.

## Key Features

- **Automatic structure separation** using connectivity analysis (pymatgen StructureGraph)
- **Parallel VASP calculations** for substrate, molecule, and complete system
- **Consistent cell handling** - all three systems use the same cell
- **Integration with existing PSTEROS modules** - follows same patterns as thermodynamics/cleavage

## Core Functions

### `separate_adsorbate_structure()`

Separates substrate+adsorbate into three components using connectivity analysis.

**Parameters:**
- `structure` (StructureData): Complete substrate+adsorbate structure
- `adsorbate_formula` (Str): Chemical formula of adsorbate (e.g., "OOH", "OH")

**Returns:**
- Dictionary with `substrate`, `molecule`, and `complete` StructureData nodes

**Example:**
```python
from teros.core import separate_adsorbate_structure
from aiida import orm

result = separate_adsorbate_structure(
    structure=complete_structure,
    adsorbate_formula=orm.Str('OH')
)

substrate = result['substrate']
molecule = result['molecule']
```

### `calculate_adsorption_energy()`

Calculates adsorption energy from component energies.

**Formula:** E_ads = E_complete - E_substrate - E_molecule

**Parameters:**
- `E_complete` (Float): Energy of substrate+adsorbate system
- `E_substrate` (Float): Energy of bare substrate
- `E_molecule` (Float): Energy of isolated molecule

**Returns:**
- Float: Adsorption energy in eV (negative = favorable)

**Example:**
```python
from teros.core import calculate_adsorption_energy
from aiida import orm

E_ads = calculate_adsorption_energy(
    E_complete=orm.Float(-100.0),
    E_substrate=orm.Float(-80.0),
    E_molecule=orm.Float(-15.0)
)
print(f"Adsorption energy: {E_ads.value:.3f} eV")
```

### `compute_adsorption_energies_scatter()`

Complete workflow: separation + VASP relaxations + energy calculation.

**Parameters:**
- `structures`: Dictionary of complete structures
- `adsorbate_formulas`: Dictionary mapping structure keys to adsorbate formulas
- `code`: VASP code
- `potential_family`: Pseudopotential family
- `potential_mapping`: Element to potential mapping
- `parameters`: VASP INCAR parameters
- `options`: Scheduler options
- `kpoints_spacing`: K-point spacing (optional)
- `clean_workdir`: Clean remote directories (default: True)

**Returns:**
Dictionary with:
- `separated_structures`: Separated systems for each input
- `substrate_energies`: Substrate energies
- `molecule_energies`: Molecule energies
- `complete_energies`: Complete system energies
- `adsorption_energies`: Final E_ads values

**Example:**
```python
from teros.core import compute_adsorption_energies_scatter

results = compute_adsorption_energies_scatter(
    structures={'site1': structure1, 'site2': structure2},
    adsorbate_formulas={'site1': 'OOH', 'site2': 'OH'},
    code=vasp_code,
    potential_family='PBE',
    potential_mapping={'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    parameters={'ENCUT': 520, 'EDIFF': 1e-6, ...},
    options={'resources': {'num_machines': 1}, ...},
)

for key, E_ads in results['adsorption_energies'].items():
    print(f"{key}: {E_ads.value:.3f} eV")
```

## Adsorbate Identification

The module uses **connectivity analysis** to identify adsorbates:

1. Builds bonding graph using pymatgen's StructureGraph with CrystalNN
2. Finds all connected molecular clusters
3. Matches cluster composition to adsorbate formula
4. Extracts adsorbate atoms while preserving substrate

**Advantages:**
- Robust to different geometries
- No geometric assumptions needed
- Handles complex adsorbates
- Works with periodic boundaries

**Requirements:**
- Adsorbate must be a single bonded cluster
- Adsorbate formula must be unique in structure
- No covalent bonds between adsorbate and substrate

## Best Practices

1. **Structure preparation:**
   - Ensure adsorbate is distinguishable from substrate
   - Use sufficient vacuum (>10 Å) to prevent periodic interactions
   - Pre-relax substrate to minimize reconstruction

2. **VASP parameters:**
   - Use consistent INCAR settings for all three calculations
   - ISIF=2 (relax ions only, keep cell fixed)
   - Converge k-points and ENCUT
   - Use appropriate ISMEAR for metallic/insulating systems

3. **Cell handling:**
   - All three systems use the same cell for consistency
   - Avoids basis set superposition error
   - Molecule stays at original position (no centering)

4. **Testing:**
   - Compare with literature values for validation
   - Check sensitivity to vacuum thickness
   - Verify substrate convergence (slab thickness, k-points)

## Integration with Experimental Tools

The module can use structures generated by experimental tools:

```python
# Use experimental surface builder to create substrate+adsorbate
from teros.experimental.adsorption_energy.lamno3.surface_builder import add_ooh_to_surface

# Generate structure with OOH
complete_structure = add_ooh_to_surface(substrate_cif)

# Calculate adsorption energy
results = compute_adsorption_energies_scatter(
    structures={'ooh_site': complete_structure},
    adsorbate_formulas={'ooh_site': 'OOH'},
    ...
)
```

## See Also

- Example: `examples/adsorption_energy/test_oh_ag111/`
- Related modules: `teros.core.thermodynamics`, `teros.core.cleavage`
- Experimental tools: `teros.experimental.adsorption_energy/`
```

**Step 2: Commit documentation**

```bash
git add docs/adsorption_energy_module.md
git commit -m "docs: add comprehensive module documentation for adsorption energy"
```

---

## Completion Checklist

- [ ] Task 1: Adsorbate separation with connectivity analysis
- [ ] Task 2: Adsorption energy calculation function
- [ ] Task 3: Scatter-gather workflow implementation
- [ ] Task 4: Core module exports updated
- [ ] Task 5: Example workflow created
- [ ] Task 6: All tests passing
- [ ] Task 7: Documentation complete

## Validation Steps

After completing all tasks:

1. **Unit tests pass:**
   ```bash
   pytest teros/core/test_adsorption_energy.py -v
   ```

2. **Imports work:**
   ```bash
   python -c "from teros.core import separate_adsorbate_structure, calculate_adsorption_energy, compute_adsorption_energies_scatter"
   ```

3. **Daemon runs without errors:**
   ```bash
   verdi daemon restart && verdi daemon status
   ```

4. **Example can be imported:**
   ```bash
   cd examples/adsorption_energy/test_oh_ag111 && python -c "import run_adsorption_energy"
   ```

5. **Integration test** (optional, requires VASP):
   ```bash
   cd examples/adsorption_energy/test_oh_ag111
   python run_adsorption_energy.py
   verdi process list  # Check it submitted successfully
   ```

## Notes

- The module follows the same patterns as `thermodynamics.py` and `cleavage.py`
- All structures maintain the same cell to avoid basis set errors
- Connectivity analysis is more robust than geometric heuristics
- The scatter-gather pattern enables efficient parallel processing (3N VASP jobs)
- Integration testing requires configured VASP code and pseudopotentials
