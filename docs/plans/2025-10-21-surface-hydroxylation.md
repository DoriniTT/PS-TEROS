# Surface Hydroxylation Module Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create PS-TEROS module that generates hydroxylated/oxygen-deficient surface variants and relaxes them in parallel with controlled concurrency.

**Architecture:** Nested WorkGraph pattern - main WorkGraph orchestrates CalcFunction for structure generation (wrapping surface_modes.py), child WorkGraph for semaphore-limited parallel VASP relaxations, and CalcFunction for result collection.

**Tech Stack:** AiiDA-WorkGraph, aiida-vasp (vasp.v2.relax), ASE (structure manipulation), NumPy

---

## Task 1: Create Module Structure

**Files:**
- Create: `teros/core/surface_hydroxylation/__init__.py`
- Create: `teros/core/surface_hydroxylation/workgraph.py`
- Create: `teros/core/surface_hydroxylation/tasks.py`
- Create: `teros/core/surface_hydroxylation/relaxations.py`
- Create: `teros/core/surface_hydroxylation/utils.py`

**Step 1: Create module directory**

```bash
mkdir -p teros/core/surface_hydroxylation
```

**Step 2: Create __init__.py with exports**

File: `teros/core/surface_hydroxylation/__init__.py`

```python
"""Surface hydroxylation module for PS-TEROS."""

from .workgraph import SurfaceHydroxylationWorkGraph

__all__ = ['SurfaceHydroxylationWorkGraph']
```

**Step 3: Create placeholder files**

Create empty files:
- `teros/core/surface_hydroxylation/workgraph.py`
- `teros/core/surface_hydroxylation/tasks.py`
- `teros/core/surface_hydroxylation/relaxations.py`
- `teros/core/surface_hydroxylation/utils.py`

Each should have a module docstring:

```python
"""[Component description]."""
```

**Step 4: Verify directory structure**

```bash
ls -la teros/core/surface_hydroxylation/
```

Expected: 5 files (__init__.py, workgraph.py, tasks.py, relaxations.py, utils.py)

**Step 5: Commit**

```bash
git add teros/core/surface_hydroxylation/
git commit -m "feat: create surface_hydroxylation module structure"
```

---

## Task 2: Move and Refactor surface_modes.py

**Files:**
- Move: `teros/experimental/vacancies_hydroxilation/surface_modes.py` → `teros/core/surface_hydroxylation/surface_modes.py`
- Modify: `teros/core/surface_hydroxylation/__init__.py`

**Step 1: Copy surface_modes.py to new location**

```bash
cp teros/experimental/vacancies_hydroxilation/surface_modes.py teros/core/surface_hydroxylation/surface_modes.py
```

**Step 2: Verify imports work**

```bash
/home/thiagotd/envs/aiida/bin/python -c "from teros.core.surface_hydroxylation.surface_modes import SurfaceModifier; print('Import successful')"
```

Expected: "Import successful"

**Step 3: Update __init__.py to export SurfaceModifier**

File: `teros/core/surface_hydroxylation/__init__.py`

```python
"""Surface hydroxylation module for PS-TEROS."""

from .workgraph import SurfaceHydroxylationWorkGraph
from .surface_modes import SurfaceModifier

__all__ = ['SurfaceHydroxylationWorkGraph', 'SurfaceModifier']
```

**Step 4: Test import from module**

```bash
/home/thiagotd/envs/aiida/bin/python -c "from teros.core.surface_hydroxylation import SurfaceModifier; print('Export successful')"
```

Expected: "Export successful"

**Step 5: Commit**

```bash
git add teros/core/surface_hydroxylation/surface_modes.py teros/core/surface_hydroxylation/__init__.py
git commit -m "feat: move surface_modes.py to core module"
```

---

## Task 3: Implement Utils (Structure Conversion)

**Files:**
- Modify: `teros/core/surface_hydroxylation/utils.py`

**Step 1: Write test for AiiDA → ASE conversion**

Create test file: `tests/core/surface_hydroxylation/test_utils.py`

```python
"""Tests for surface_hydroxylation utilities."""

import pytest
from aiida.orm import StructureData
from ase import Atoms
from teros.core.surface_hydroxylation.utils import aiida_to_ase, ase_to_aiida


def test_aiida_to_ase_conversion():
    """Test StructureData → ASE Atoms conversion."""
    # Create simple AiiDA structure (2-atom H2 molecule)
    structure = StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    structure.append_atom(position=(5.0, 5.0, 5.0), symbols='H')
    structure.append_atom(position=(5.5, 5.0, 5.0), symbols='H')

    # Convert to ASE
    atoms = aiida_to_ase(structure)

    # Verify
    assert isinstance(atoms, Atoms)
    assert len(atoms) == 2
    assert atoms.get_chemical_symbols() == ['H', 'H']
```

**Step 2: Run test to verify it fails**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_utils.py::test_aiida_to_ase_conversion -v
```

Expected: FAIL with "cannot import name 'aiida_to_ase'"

**Step 3: Implement aiida_to_ase**

File: `teros/core/surface_hydroxylation/utils.py`

```python
"""Utility functions for surface_hydroxylation module."""

from aiida.orm import StructureData
from ase import Atoms


def aiida_to_ase(structure: StructureData) -> Atoms:
    """
    Convert AiiDA StructureData to ASE Atoms.

    Args:
        structure: AiiDA StructureData object

    Returns:
        ASE Atoms object
    """
    return structure.get_ase()


def ase_to_aiida(atoms: Atoms) -> StructureData:
    """
    Convert ASE Atoms to AiiDA StructureData.

    Args:
        atoms: ASE Atoms object

    Returns:
        AiiDA StructureData object
    """
    return StructureData(ase=atoms)
```

**Step 4: Run test to verify it passes**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_utils.py::test_aiida_to_ase_conversion -v
```

Expected: PASS

**Step 5: Write test for ASE → AiiDA conversion**

Add to `tests/core/surface_hydroxylation/test_utils.py`:

```python
def test_ase_to_aiida_conversion():
    """Test ASE Atoms → StructureData conversion."""
    # Create ASE Atoms
    atoms = Atoms('H2', positions=[[0, 0, 0], [0.74, 0, 0]], cell=[10, 10, 10], pbc=True)

    # Convert to AiiDA
    structure = ase_to_aiida(atoms)

    # Verify
    assert isinstance(structure, StructureData)
    assert len(structure.sites) == 2
```

**Step 6: Run test to verify it passes**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_utils.py::test_ase_to_aiida_conversion -v
```

Expected: PASS

**Step 7: Commit**

```bash
git add teros/core/surface_hydroxylation/utils.py tests/core/surface_hydroxylation/test_utils.py
git commit -m "feat: add structure conversion utilities"
```

---

## Task 4: Implement generate_structures CalcFunction

**Files:**
- Modify: `teros/core/surface_hydroxylation/tasks.py`

**Step 1: Write test for generate_structures**

Create: `tests/core/surface_hydroxylation/test_tasks.py`

```python
"""Tests for surface_hydroxylation task functions."""

import pytest
from aiida.orm import StructureData, Dict
from aiida.engine import run
from ase.build import fcc111
from teros.core.surface_hydroxylation.tasks import generate_structures


def test_generate_structures_hydrogen_mode():
    """Test structure generation in hydrogen mode."""
    # Create simple test slab (2x2 Pt(111) with O adlayer)
    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
    # Add O atoms on top
    slab.extend(fcc111('O', size=(2, 2, 1), vacuum=0.0))
    slab.center(vacuum=10.0, axis=2)

    structure = StructureData(ase=slab)

    params = Dict({
        'mode': 'hydrogen',
        'species': 'O',
        'z_window': 0.5,
        'which_surface': 'top',
        'oh_dist': 0.98,
        'include_empty': False,
        'deduplicate_by_coverage': True,
        'coverage_bins': 3
    })

    # Run calcfunction
    result = run(generate_structures, structure=structure, params=params)

    # Verify outputs
    assert 'manifest' in result
    assert 'structures' in result
    assert isinstance(result['manifest'], Dict)
    assert isinstance(result['structures'], list)
    assert len(result['structures']) > 0
    assert all(isinstance(s, StructureData) for s in result['structures'])
```

**Step 2: Run test to verify it fails**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_tasks.py::test_generate_structures_hydrogen_mode -v
```

Expected: FAIL with "cannot import name 'generate_structures'"

**Step 3: Implement generate_structures**

File: `teros/core/surface_hydroxylation/tasks.py`

```python
"""Task functions (CalcFunctions) for surface_hydroxylation module."""

import tempfile
from pathlib import Path
import json

from aiida.engine import calcfunction
from aiida.orm import StructureData, Dict, List

from .utils import aiida_to_ase, ase_to_aiida
from .surface_modes import SurfaceModifier


@calcfunction
def generate_structures(structure: StructureData, params: Dict) -> dict:
    """
    Generate surface variants using surface_modes.py.

    Args:
        structure: Input relaxed slab structure
        params: Dict with surface_modes parameters:
            - mode: str ('vacancies'/'hydrogen'/'combine')
            - species: str (default 'O')
            - z_window: float (default 0.5)
            - which_surface: str ('top'/'bottom'/'both')
            - oh_dist: float (default 0.98)
            - include_empty: bool (default False)
            - supercell: list[int] or None
            - deduplicate_by_coverage: bool
            - coverage_bins: int or None

    Returns:
        dict with:
            - manifest: Dict (parsed manifest from surface_modes)
            - structures: List (generated StructureData variants)
    """
    # Convert AiiDA → ASE
    atoms = aiida_to_ase(structure)

    # Extract parameters with defaults
    p = params.get_dict()
    mode = p.get('mode', 'hydrogen')
    species = p.get('species', 'O')
    z_window = p.get('z_window', 0.5)
    which_surface = p.get('which_surface', 'top')
    oh_dist = p.get('oh_dist', 0.98)
    include_empty = p.get('include_empty', False)
    supercell = p.get('supercell', None)
    deduplicate = p.get('deduplicate_by_coverage', False)
    coverage_bins = p.get('coverage_bins', None)

    # Convert supercell to tuple if provided
    if supercell is not None:
        supercell = tuple(supercell)

    # Create temporary directory for outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir)

        # Create SurfaceModifier instance
        sm = SurfaceModifier(
            atoms=atoms,
            species=species,
            z_window=z_window,
            which_surface=which_surface,
            oh_dist=oh_dist,
            include_empty=include_empty,
            outdir=outdir,
            fmt='vasp',
            supercell=supercell,
            deduplicate_by_coverage=deduplicate,
            coverage_bins=coverage_bins
        )

        # Run appropriate mode
        if mode == 'vacancies':
            manifest = sm.run_vacancies()
        elif mode == 'hydrogen':
            manifest = sm.run_hydrogen()
        elif mode == 'combine':
            manifest = sm.run_combine()
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Read generated structures
        structures = []
        for variant in manifest['variants']:
            filepath = Path(variant['file'])
            from ase.io import read
            variant_atoms = read(filepath.as_posix())
            structures.append(ase_to_aiida(variant_atoms))

    # Return manifest and structures
    return {
        'manifest': Dict(dict=manifest),
        'structures': List(list=structures)
    }
```

**Step 4: Run test to verify it passes**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_tasks.py::test_generate_structures_hydrogen_mode -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add teros/core/surface_hydroxylation/tasks.py tests/core/surface_hydroxylation/test_tasks.py
git commit -m "feat: implement generate_structures calcfunction"
```

---

## Task 5: Implement collect_results CalcFunction

**Files:**
- Modify: `teros/core/surface_hydroxylation/tasks.py`

**Step 1: Write test for collect_results**

Add to `tests/core/surface_hydroxylation/test_tasks.py`:

```python
from aiida.orm import Int, Float, Namespace
from teros.core.surface_hydroxylation.tasks import collect_results


def test_collect_results_with_mixed_success_failure():
    """Test result collection with some successes and failures."""
    # Create mock relaxation outputs
    from aiida.orm import Node

    # Mock successful relaxation
    success_node = Node()
    success_node.store()
    success_node.set_extra('exit_status', 0)
    success_node.set_extra('total_energy', -123.45)

    # Mock failed relaxation
    fail_node = Node()
    fail_node.store()
    fail_node.set_extra('exit_status', 400)
    fail_node.set_extra('error_message', 'Convergence failed')

    # Create namespace
    relaxations = Namespace({
        'relax_0': success_node,
        'relax_1': fail_node
    })

    # Create manifest
    manifest = Dict({
        'variants': [
            {'name': 'oh_001_0.5', 'OH_coverage': 0.5},
            {'name': 'oh_002_1.0', 'OH_coverage': 1.0}
        ]
    })

    # Run collect_results
    result = run(collect_results, relaxations=relaxations, manifest=manifest)

    # Verify
    assert 'successful_relaxations' in result
    assert 'failed_relaxations' in result
    assert 'statistics' in result

    stats = result['statistics'].get_dict()
    assert stats['total'] == 2
    assert stats['succeeded'] == 1
    assert stats['failed'] == 1
```

**Step 2: Run test to verify it fails**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_tasks.py::test_collect_results_with_mixed_success_failure -v
```

Expected: FAIL with "cannot import name 'collect_results'"

**Step 3: Implement collect_results**

Add to `teros/core/surface_hydroxylation/tasks.py`:

```python
@calcfunction
def collect_results(relaxations, manifest: Dict) -> dict:
    """
    Collect and organize relaxation results.

    Args:
        relaxations: Namespace {relax_0, relax_1, ..., relax_N} from RelaxationsWorkGraph
        manifest: Original manifest Dict from generate_structures

    Returns:
        dict with:
            - successful_relaxations: List of Dicts with structure, energy, coverage, metadata
            - failed_relaxations: List of Dicts with name, coverage, error info
            - statistics: Dict with total, succeeded, failed counts
    """
    manifest_dict = manifest.get_dict()
    variants = manifest_dict['variants']

    successful = []
    failed = []

    # Iterate through relaxations by index
    for idx, variant in enumerate(variants):
        relax_key = f'relax_{idx}'

        # Get relaxation output
        if relax_key not in relaxations:
            # Relaxation not found (shouldn't happen)
            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage'),
                'error_message': 'Relaxation output not found'
            })
            continue

        relax_output = relaxations[relax_key]

        # Check exit status
        exit_status = relax_output.exit_status if hasattr(relax_output, 'exit_status') else relax_output.get_extra('exit_status', 1)

        if exit_status == 0:
            # Success - extract structure and energy
            structure = relax_output.outputs.structure if hasattr(relax_output, 'outputs') else None
            energy = relax_output.outputs.total_energy if hasattr(relax_output, 'outputs') else Float(relax_output.get_extra('total_energy'))

            successful.append({
                'name': variant['name'],
                'structure': structure,
                'energy': energy,
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage'),
                'metadata': Dict(dict=variant)
            })
        else:
            # Failure - record error
            error_msg = relax_output.get_extra('error_message', 'Unknown error') if hasattr(relax_output, 'get_extra') else 'Unknown error'

            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage'),
                'exit_status': Int(exit_status),
                'error_message': str(error_msg)
            })

    # Calculate statistics
    statistics = Dict(dict={
        'total': len(variants),
        'succeeded': len(successful),
        'failed': len(failed)
    })

    return {
        'successful_relaxations': List(list=successful),
        'failed_relaxations': List(list=failed),
        'statistics': statistics
    }
```

**Step 4: Run test to verify it passes**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/test_tasks.py::test_collect_results_with_mixed_success_failure -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add teros/core/surface_hydroxylation/tasks.py tests/core/surface_hydroxylation/test_tasks.py
git commit -m "feat: implement collect_results calcfunction"
```

---

## Task 6: Implement RelaxationsWorkGraph with Semaphore

**Files:**
- Modify: `teros/core/surface_hydroxylation/relaxations.py`

**Step 1: Study existing scatter-gather VASP patterns**

Read existing PS-TEROS code to understand pattern:

```bash
grep -r "scatter.*relax" teros/core/ | head -5
```

Look for existing WorkGraph examples with semaphore usage.

**Step 2: Write skeleton RelaxationsWorkGraph**

File: `teros/core/surface_hydroxylation/relaxations.py`

```python
"""Child WorkGraph for parallel VASP relaxations."""

from aiida_workgraph import WorkGraph
from aiida.orm import List, Dict, Int


class RelaxationsWorkGraph(WorkGraph):
    """WorkGraph for parallel VASP relaxations with semaphore limiting."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setup_inputs()
        self.setup_tasks()

    def setup_inputs(self):
        """Define input sockets."""
        self.add_input('structures', List, required=True)
        self.add_input('builder_config', Dict, required=True)
        self.add_input('max_parallel', Int, required=True)

    def setup_tasks(self):
        """Create VASP relaxation tasks with semaphore."""
        # Get inputs
        structures = self.inputs.structures.value
        builder_config = self.inputs.builder_config.value
        max_parallel = self.inputs.max_parallel.value

        # Create semaphore
        semaphore = self.add_semaphore('relax_semaphore', max_parallel)

        # Create relaxation task for each structure
        for idx, structure in enumerate(structures):
            task = self.add_task(
                'aiida.workflows:vasp.v2.relax',
                name=f'relax_{idx}'
            )

            # Set structure input
            task.set_input('structure', structure)

            # Apply builder config
            self._apply_builder_config(task, builder_config)

            # Apply semaphore
            task.waiting_on = semaphore

            # Add to outputs
            self.add_output(f'relax_{idx}', task.outputs.output)

    def _apply_builder_config(self, task, config):
        """Apply builder configuration to VASP task."""
        config_dict = config.get_dict() if hasattr(config, 'get_dict') else config

        # Apply each config key to task
        for key, value in config_dict.items():
            if hasattr(task, key):
                setattr(task, key, value)
            else:
                # Try setting as input
                try:
                    task.set_input(key, value)
                except Exception:
                    pass  # Skip unknown keys
```

**Step 3: Verify imports work**

```bash
/home/thiagotd/envs/aiida/bin/python -c "from teros.core.surface_hydroxylation.relaxations import RelaxationsWorkGraph; print('Import successful')"
```

Expected: "Import successful"

**Step 4: Commit**

```bash
git add teros/core/surface_hydroxylation/relaxations.py
git commit -m "feat: implement RelaxationsWorkGraph with semaphore"
```

---

## Task 7: Implement Main SurfaceHydroxylationWorkGraph

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py`

**Step 1: Write skeleton main WorkGraph**

File: `teros/core/surface_hydroxylation/workgraph.py`

```python
"""Main WorkGraph for surface hydroxylation workflow."""

from aiida_workgraph import WorkGraph
from aiida.orm import StructureData, Dict, Int

from .tasks import generate_structures, collect_results
from .relaxations import RelaxationsWorkGraph


class SurfaceHydroxylationWorkGraph(WorkGraph):
    """
    Main workflow for surface hydroxylation studies.

    Inputs:
        structure: StructureData - input relaxed slab
        surface_params: Dict - parameters for surface_modes.py
        builder_config: Dict - complete VASP builder configuration
        max_parallel_jobs: Int - max concurrent VASP jobs

    Outputs:
        successful_relaxations: List of successful relaxation results
        failed_relaxations: List of failed relaxation results
        statistics: Dict with total, succeeded, failed counts
        manifest: Dict with surface_modes manifest
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setup_workflow()

    def setup_workflow(self):
        """Set up the workflow structure."""
        # Input sockets
        self.add_input('structure', StructureData, required=True)
        self.add_input('surface_params', Dict, required=True)
        self.add_input('builder_config', Dict, required=True)
        self.add_input('max_parallel_jobs', Int, required=True)

        # Task 1: Generate structures
        gen_task = self.add_task(
            generate_structures,
            name='generate_structures',
            structure=self.inputs.structure,
            params=self.inputs.surface_params
        )

        # Task 2: Relaxations WorkGraph
        relax_wg = self.add_task(
            RelaxationsWorkGraph,
            name='relaxations',
            structures=gen_task.outputs.structures,
            builder_config=self.inputs.builder_config,
            max_parallel=self.inputs.max_parallel_jobs
        )

        # Task 3: Collect results
        collect_task = self.add_task(
            collect_results,
            name='collect_results',
            relaxations=relax_wg.outputs,
            manifest=gen_task.outputs.manifest
        )

        # Output sockets
        self.add_output('successful_relaxations', collect_task.outputs.successful_relaxations)
        self.add_output('failed_relaxations', collect_task.outputs.failed_relaxations)
        self.add_output('statistics', collect_task.outputs.statistics)
        self.add_output('manifest', gen_task.outputs.manifest)
```

**Step 2: Verify imports work**

```bash
/home/thiagotd/envs/aiida/bin/python -c "from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph; print('Import successful')"
```

Expected: "Import successful"

**Step 3: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: implement main SurfaceHydroxylationWorkGraph"
```

---

## Task 8: Create Test Example

**Files:**
- Create: `examples/surface_hydroxylation/run_test.py`
- Create: `examples/surface_hydroxylation/README.md`

**Step 1: Create example directory**

```bash
mkdir -p examples/surface_hydroxylation
```

**Step 2: Create test run script**

File: `examples/surface_hydroxylation/run_test.py`

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
Test script for surface_hydroxylation module.

Uses small Pt(111) slab with O adlayer for fast testing.
"""

from aiida import orm, load_profile
from aiida_workgraph import WorkGraph
from ase.build import fcc111
from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph

# Load AiiDA profile
load_profile('psteros')

# Create test slab (2x2 Pt(111) with O adlayer)
print("Creating test structure...")
slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
slab.extend(fcc111('O', size=(2, 2, 1), vacuum=0.0))
slab.center(vacuum=10.0, axis=2)

structure = orm.StructureData(ase=slab)
print(f"Test slab: {len(slab)} atoms")

# Surface generation parameters
surface_params = orm.Dict({
    'mode': 'hydrogen',
    'species': 'O',
    'z_window': 0.5,
    'which_surface': 'top',
    'oh_dist': 0.98,
    'include_empty': False,
    'deduplicate_by_coverage': True,
    'coverage_bins': 3
})

# VASP builder configuration (lightweight for testing)
builder_config = orm.Dict({
    'metadata': {
        'options': {
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4},
            'max_wallclock_seconds': 3600,
            'account': 'your_account',
            'queue_name': 'debug'
        }
    },
    'kpoints_distance': orm.Float(0.5),  # Coarse for testing
    'parameters': orm.Dict({
        'EDIFF': 1e-4,
        'ENCUT': 400,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 10,  # Very few steps for testing
        'IBRION': 2,
        'ISIF': 2
    }),
    'relax': orm.Dict({
        'positions': True,
        'shape': False,
        'volume': False
    })
})

# Parallelization limit
max_parallel = orm.Int(2)

print("\nCreating workflow...")
wg = SurfaceHydroxylationWorkGraph()
wg.name = 'surface_hydroxylation_test'

# Set inputs
wg.inputs.structure = structure
wg.inputs.surface_params = surface_params
wg.inputs.builder_config = builder_config
wg.inputs.max_parallel_jobs = max_parallel

# Submit workflow
print("Submitting workflow...")
result = wg.submit()
print(f"Workflow submitted: PK = {result.pk}")
print(f"\nMonitor with: verdi process show {result.pk}")
print(f"Wait and check: sleep 60 && verdi process show {result.pk}")
```

**Step 3: Create README documentation**

File: `examples/surface_hydroxylation/README.md`

```markdown
# Surface Hydroxylation Module - Test Example

This example demonstrates the surface_hydroxylation module with a small test system.

## Test System

- **Structure:** 2×2 Pt(111) slab with O adlayer (~20 atoms)
- **Mode:** Hydrogen (hydroxylation only)
- **Coverage bins:** 3 (generates ~3-5 structures)
- **Max parallel:** 2 VASP jobs at once

## Running the Test

1. Ensure AiiDA daemon is running:
   ```bash
   verdi daemon status
   verdi daemon start  # if not running
   ```

2. Run the test script:
   ```bash
   /home/thiagotd/envs/aiida/bin/python run_test.py
   ```

3. Monitor progress:
   ```bash
   verdi process show <PK>
   verdi process report <PK>
   ```

4. Wait for completion (~5-10 minutes with lightweight settings):
   ```bash
   sleep 300
   verdi process show <PK>
   ```

## Expected Results

**Success criteria:**
- Main workflow exits with status `[0]`
- `outputs.statistics` shows total structures generated
- `outputs.successful_relaxations` contains relaxed structures and energies
- All VASP relaxation tasks visible in process tree

**Check outputs:**
```python
from aiida import orm
node = orm.load_node(<PK>)

# Check statistics
stats = node.outputs.statistics.get_dict()
print(f"Total: {stats['total']}, Succeeded: {stats['succeeded']}, Failed: {stats['failed']}")

# Check successful relaxations
for relax in node.outputs.successful_relaxations:
    print(f"{relax['name']}: {relax['energy'].value} eV")
```

## Troubleshooting

**No structures generated:**
- Check surface_params (z_window might be too small)
- Verify O atoms are present in input slab

**All relaxations fail:**
- Check VASP configuration in builder_config
- Verify compute resources are correct
- Check VASP pseudopotentials are available

**Semaphore not limiting parallelization:**
- Check RelaxationsWorkGraph implementation
- Verify max_parallel_jobs is being used correctly
```

**Step 4: Make script executable**

```bash
chmod +x examples/surface_hydroxylation/run_test.py
```

**Step 5: Commit**

```bash
git add examples/surface_hydroxylation/
git commit -m "feat: add test example for surface_hydroxylation"
```

---

## Task 9: Integration Testing

**Step 1: Clear Python cache**

```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

**Step 2: Restart AiiDA daemon**

```bash
verdi daemon restart
```

**Step 3: Run test example**

```bash
cd examples/surface_hydroxylation
/home/thiagotd/envs/aiida/bin/python run_test.py
```

Expected output:
```
Creating test structure...
Test slab: 20 atoms

Creating workflow...
Submitting workflow...
Workflow submitted: PK = XXXX

Monitor with: verdi process show XXXX
```

**Step 4: Monitor workflow progress**

```bash
verdi process show <PK>
```

Expected: Workflow running, tasks visible

**Step 5: Wait for completion**

```bash
sleep 300  # Wait 5 minutes
verdi process show <PK>
```

Expected: Exit status [0], all tasks completed

**Step 6: Verify outputs**

```bash
verdi process report <PK>
```

Expected: No errors, successful completion message

**Step 7: Document test results**

If test passes, commit verification:

```bash
git add -A
git commit -m "test: verify surface_hydroxylation module integration"
```

---

## Task 10: Production Testing

**Files:**
- Create: `examples/surface_hydroxylation/run_production.py`

**Step 1: Create production example**

File: `examples/surface_hydroxylation/run_production.py`

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
Production example for surface_hydroxylation module.

Uses realistic perovskite oxide surface with production VASP settings.
"""

from aiida import orm, load_profile
from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph

# Load AiiDA profile
load_profile('psteros')

# Load production slab structure
# TODO: Replace with actual relaxed surface from surface_thermodynamics
structure = orm.load_node(<PK_of_relaxed_slab>)

print(f"Production slab: {len(structure.sites)} atoms")

# Surface generation parameters - combined mode
surface_params = orm.Dict({
    'mode': 'combine',  # Vacancies + hydroxylation
    'species': 'O',
    'z_window': 0.5,
    'which_surface': 'top',
    'oh_dist': 0.98,
    'include_empty': False,
    'deduplicate_by_coverage': True,
    'coverage_bins': 10  # More bins for better coverage sampling
})

# Production VASP settings
builder_config = orm.Dict({
    'metadata': {
        'options': {
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
            'max_wallclock_seconds': 3600 * 10,  # 10 hours
            'account': 'your_production_account',
            'queue_name': 'normal'
        }
    },
    'kpoints_distance': orm.Float(0.3),  # Converged k-points
    'parameters': orm.Dict({
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 200,
        'IBRION': 2,
        'ISIF': 2,
        'LREAL': False,
        'PREC': 'Accurate'
    }),
    'relax': orm.Dict({
        'positions': True,
        'shape': False,
        'volume': False
    })
})

# Limit parallel jobs based on cluster capacity
max_parallel = orm.Int(5)

print("\nCreating production workflow...")
wg = SurfaceHydroxylationWorkGraph()
wg.name = 'surface_hydroxylation_production'

wg.inputs.structure = structure
wg.inputs.surface_params = surface_params
wg.inputs.builder_config = builder_config
wg.inputs.max_parallel_jobs = max_parallel

print("Submitting production workflow...")
result = wg.submit()
print(f"Workflow submitted: PK = {result.pk}")
print(f"\nThis will take several hours. Monitor with:")
print(f"  verdi process show {result.pk}")
```

**Step 2: Add production notes to README**

Add to `examples/surface_hydroxylation/README.md`:

```markdown
## Production Usage

See `run_production.py` for production example with:
- Realistic surface structure (~100 atoms)
- Combined mode (vacancies + hydroxylation)
- 10 coverage bins (~15-20 structures)
- Production VASP settings
- 5 parallel jobs

Expected runtime: Several hours depending on structure size and settings.
```

**Step 3: Commit**

```bash
git add examples/surface_hydroxylation/run_production.py examples/surface_hydroxylation/README.md
git commit -m "docs: add production example and usage notes"
```

---

## Final Verification

**Step 1: Verify all tests pass**

```bash
/home/thiagotd/envs/aiida/bin/python -m pytest tests/core/surface_hydroxylation/ -v
```

Expected: All tests PASS

**Step 2: Verify imports work**

```bash
/home/thiagotd/envs/aiida/bin/python -c "from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph; print('Module ready')"
```

Expected: "Module ready"

**Step 3: Check git status**

```bash
git status
```

Expected: All changes committed, working directory clean

**Step 4: Review commit history**

```bash
git log --oneline | head -15
```

Expected: ~10-12 commits with clear messages

---

## Completion Checklist

- [ ] Module structure created (`teros/core/surface_hydroxylation/`)
- [ ] `surface_modes.py` moved and refactored
- [ ] Utils implemented (structure conversion)
- [ ] `generate_structures` calcfunction implemented
- [ ] `collect_results` calcfunction implemented
- [ ] `RelaxationsWorkGraph` with semaphore implemented
- [ ] Main `SurfaceHydroxylationWorkGraph` implemented
- [ ] Test example created and verified
- [ ] Production example documented
- [ ] All tests passing
- [ ] Documentation complete

**When all items checked:** Module is complete and ready for merge to develop branch.
