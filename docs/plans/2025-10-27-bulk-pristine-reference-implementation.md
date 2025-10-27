# Bulk and Pristine Slab Reference Calculations Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add automatic bulk and pristine slab reference calculations to the surface hydroxylation workflow to enable surface energy analysis per Section S2.

**Architecture:** Add two sequential VaspTask instances (bulk relaxation with ISIF=3, pristine slab relaxation with ISIF=2) before existing structure generation in SurfaceHydroxylationWorkGraph. Extend input parameters to accept bulk structure (CIF/PK/StructureData) and bulk_builder_inputs. Add four new outputs: bulk_structure, bulk_energy, pristine_structure, pristine_energy.

**Tech Stack:** AiiDA, AiiDA-WorkGraph, VASP, ASE (for CIF reading)

**Design Document:** `/home/thiagotd/git/fosfato/calculos/hydroxylation/docs/plans/2025-10-27-bulk-pristine-reference-calculations-design.md`

---

## Task 1: Add New Parameters to SurfaceHydroxylationWorkGraph

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:22-32`

**Step 1: Add bulk parameters to function signature**

Modify the `@task.graph` decorator outputs and function signature:

```python
@task.graph(outputs=[
    'manifest',
    'structures',
    'energies',
    'bulk_structure',        # NEW
    'bulk_energy',           # NEW
    'pristine_structure',    # NEW
    'pristine_energy',       # NEW
])
def SurfaceHydroxylationWorkGraph(
    structure: orm.StructureData,
    surface_params: dict,
    code: orm.InstalledCode,
    builder_inputs: dict,
    bulk_structure: orm.StructureData,       # NEW - Required
    bulk_builder_inputs: dict,               # NEW - Required
    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    structure_specific_builder_inputs: dict = None,
) -> dict:
```

**Step 2: Update docstring**

Add to the docstring Args section (after `builder_inputs`):

```python
    """
    ...existing docstring...

    Args:
        ...existing parameters...
        bulk_structure: Bulk crystal structure for reference calculations (StructureData).
                       Required for surface energy calculations per Section S2.
        bulk_builder_inputs: VASP parameters for bulk relaxation (dict). Must include
                            ISIF=3 for cell relaxation. Same format as builder_inputs.
        ...remaining parameters...
```

**Step 3: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add bulk_structure and bulk_builder_inputs parameters to SurfaceHydroxylationWorkGraph"
```

---

## Task 2: Add Bulk Relaxation Task

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:133-145`

**Step 1: Add bulk relaxation task implementation**

Insert after line 136 (after `surface_params = orm.Dict(dict=surface_params)`):

```python
    # =========================================================================
    # Task 0: Bulk Relaxation (NEW)
    # =========================================================================
    from aiida_workgraph.tasks import VaspWorkChain as VaspTask

    # Prepare settings with parser configuration
    bulk_settings = bulk_builder_inputs.get('settings', {})
    if not bulk_settings:
        bulk_settings = {
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }
    if not isinstance(bulk_settings, orm.Dict):
        bulk_settings = orm.Dict(dict=bulk_settings)

    # Create bulk relaxation task
    bulk_vasp = VaspTask(
        structure=bulk_structure,
        code=code,
        parameters=orm.Dict(dict=bulk_builder_inputs['parameters']),
        kpoints_spacing=orm.Float(bulk_builder_inputs.get('kpoints_spacing', 0.5)),
        potential_family=orm.Str(bulk_builder_inputs['potential_family']),
        potential_mapping=orm.Dict(dict=bulk_builder_inputs.get('potential_mapping', {})),
        options=orm.Dict(dict=bulk_builder_inputs['options']),
        settings=bulk_settings,
        clean_workdir=orm.Bool(bulk_builder_inputs.get('clean_workdir', False)),
    )

    # Extract energy from bulk relaxation
    from .utils import extract_total_energy_from_misc
    bulk_energy = extract_total_energy_from_misc(misc=bulk_vasp.outputs.misc)
```

**Step 2: Check if extract_total_energy_from_misc exists**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/bulk-pristine-reference
grep -n "extract_total_energy" teros/core/surface_hydroxylation/utils.py
```

Expected: Function exists or needs to be created

**Step 3: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add bulk relaxation task to workflow"
```

---

## Task 3: Add Pristine Slab Relaxation Task

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py` (after bulk task)

**Step 1: Add pristine slab relaxation task**

Insert immediately after the bulk relaxation task:

```python
    # =========================================================================
    # Task 0.5: Pristine Slab Relaxation (NEW)
    # =========================================================================

    # Prepare settings for pristine slab
    pristine_settings = builder_inputs.get('settings', {})
    if not pristine_settings:
        pristine_settings = {
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }
    if not isinstance(pristine_settings, orm.Dict):
        pristine_settings = orm.Dict(dict=pristine_settings)

    # Create pristine slab relaxation task
    pristine_vasp = VaspTask(
        structure=structure,  # Use input slab structure
        code=code,
        parameters=orm.Dict(dict=builder_inputs['parameters']),
        kpoints_spacing=orm.Float(builder_inputs.get('kpoints_spacing', 0.5)),
        potential_family=orm.Str(builder_inputs['potential_family']),
        potential_mapping=orm.Dict(dict=builder_inputs.get('potential_mapping', {})),
        options=orm.Dict(dict=builder_inputs['options']),
        settings=pristine_settings,
        clean_workdir=orm.Bool(builder_inputs.get('clean_workdir', False)),
    )

    pristine_energy = extract_total_energy_from_misc(misc=pristine_vasp.outputs.misc)
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add pristine slab relaxation task to workflow"
```

---

## Task 4: Update Return Statement with New Outputs

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:160-167`

**Step 1: Add new outputs to return dict**

Modify the return statement at the end of `SurfaceHydroxylationWorkGraph`:

```python
    # =========================================================================
    # Return outputs including reference calculations
    # =========================================================================
    return {
        'manifest': gen_outputs.manifest,
        'structures': relax_outputs.structures,
        'energies': relax_outputs.energies,
        'bulk_structure': bulk_vasp.outputs.structure,      # NEW
        'bulk_energy': bulk_energy,                         # NEW
        'pristine_structure': pristine_vasp.outputs.structure,  # NEW
        'pristine_energy': pristine_energy,                 # NEW
    }
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add bulk and pristine outputs to workflow return"
```

---

## Task 5: Add Bulk Parameters to build_surface_hydroxylation_workgraph

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:169-180`

**Step 1: Add bulk parameters to builder function**

Modify function signature:

```python
def build_surface_hydroxylation_workgraph(
    structure: orm.StructureData = None,
    structure_pk: int = None,
    surface_params: dict = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    builder_inputs: dict = None,

    # NEW BULK PARAMETERS
    bulk_structure: orm.StructureData = None,    # NEW
    bulk_structure_pk: int = None,               # NEW
    bulk_cif_path: str = None,                   # NEW
    bulk_builder_inputs: dict = None,            # NEW

    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    structure_specific_builder_inputs: dict = None,
    name: str = 'SurfaceHydroxylation',
) -> WorkGraph:
```

**Step 2: Update docstring**

Add to Args section:

```python
        bulk_structure: Bulk crystal structure (StructureData).
                       Provide one of: bulk_structure, bulk_structure_pk, or bulk_cif_path.
        bulk_structure_pk: PK of bulk StructureData node (int).
        bulk_cif_path: Path to CIF file for bulk structure (str).
                      Example: 'ag3po4.cif'
        bulk_builder_inputs: VASP parameters for bulk relaxation (dict).
                            Must include ISIF=3 for cell relaxation.
                            If None, uses sensible defaults.
```

**Step 3: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add bulk parameters to build_surface_hydroxylation_workgraph"
```

---

## Task 6: Add Bulk Structure Validation Logic

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:344-392`

**Step 1: Add validation after code loading (after line 343)**

Insert after the code loading section:

```python
    # ========================================================================
    # BULK STRUCTURE VALIDATION (NEW)
    # ========================================================================

    # Validate exactly one bulk input is provided
    bulk_inputs_provided = sum([
        bulk_structure is not None,
        bulk_structure_pk is not None,
        bulk_cif_path is not None,
    ])

    if bulk_inputs_provided == 0:
        raise ValueError(
            "Bulk structure is required for surface energy calculations.\n"
            "Provide one of: bulk_structure, bulk_structure_pk, or bulk_cif_path.\n\n"
            "Examples:\n"
            "  bulk_cif_path='ag3po4.cif'\n"
            "  bulk_structure_pk=1234\n"
            "  bulk_structure=orm.load_node(1234)"
        )

    if bulk_inputs_provided > 1:
        raise ValueError(
            "Provide exactly ONE bulk structure input.\n"
            f"Found: bulk_structure={bulk_structure is not None}, "
            f"bulk_structure_pk={bulk_structure_pk is not None}, "
            f"bulk_cif_path={bulk_cif_path is not None}"
        )

    # Load/convert bulk structure to StructureData
    if bulk_cif_path is not None:
        try:
            from ase.io import read
            from pathlib import Path
            cif_path = Path(bulk_cif_path)
            if not cif_path.exists():
                raise FileNotFoundError(f"CIF file not found: {bulk_cif_path}")
            atoms = read(str(cif_path))
            bulk_structure = orm.StructureData(ase=atoms)
            print(f"   ✓ Loaded bulk structure from CIF: {bulk_cif_path}")
            print(f"   Formula: {atoms.get_chemical_formula()}")
        except Exception as e:
            raise ValueError(
                f"Could not read CIF file: {bulk_cif_path}\n"
                f"Error: {e}\n"
                f"Check file exists and is valid CIF format."
            )

    elif bulk_structure_pk is not None:
        try:
            bulk_structure = orm.load_node(bulk_structure_pk)
            if not isinstance(bulk_structure, orm.StructureData):
                raise ValueError(
                    f"Node {bulk_structure_pk} is not a StructureData.\n"
                    f"Found: {type(bulk_structure).__name__}"
                )
            print(f"   ✓ Loaded bulk structure from PK: {bulk_structure_pk}")
        except Exception as e:
            raise ValueError(
                f"Could not load bulk structure from PK {bulk_structure_pk}.\n"
                f"Error: {e}\n"
                f"Use 'verdi data core.structure list' to see available structures."
            )

    # If bulk_structure provided directly, validate type
    elif bulk_structure is not None:
        if not isinstance(bulk_structure, orm.StructureData):
            raise ValueError(
                f"bulk_structure must be StructureData, got {type(bulk_structure).__name__}"
            )
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add bulk structure validation and loading logic"
```

---

## Task 7: Add Default bulk_builder_inputs

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py` (after bulk validation)

**Step 1: Add default bulk builder inputs**

Insert after bulk structure validation:

```python
    # Set default bulk_builder_inputs if not provided
    if bulk_builder_inputs is None:
        print("   ⚠ Using default bulk_builder_inputs (ISIF=3, ENCUT=500)")
        # Default bulk relaxation parameters (based on default_ag3po4_builders.py)
        bulk_builder_inputs = {
            'parameters': {
                'incar': {
                    'PREC': 'Accurate',
                    'ENCUT': 500,
                    'EDIFF': 1e-6,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ALGO': 'Fast',
                    'LREAL': False,
                    'NELM': 200,
                    'LWAVE': False,
                    'LCHARG': False,
                    'ISIF': 3,        # Cell + ionic relaxation for bulk
                    'NSW': 500,
                    'IBRION': 2,
                    'EDIFFG': -0.01,  # Tighter convergence for bulk
                }
            },
            'kpoints_spacing': 0.3,  # Denser k-points for bulk
            'potential_family': 'PBE',
            'potential_mapping': {},
            'options': {
                'resources': {
                    'num_machines': 1,
                    'num_mpiprocs_per_machine': 16,
                },
                'queue_name': 'normal',
                'max_wallclock_seconds': 3600 * 4,
            },
            'clean_workdir': False,
        }
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add default bulk_builder_inputs with ISIF=3"
```

---

## Task 8: Update WorkGraph Build Call

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:405-421`

**Step 1: Pass bulk parameters to WorkGraph.build**

Modify the `SurfaceHydroxylationWorkGraph.build()` call:

```python
    # Build the WorkGraph using the @task.graph function
    wg = SurfaceHydroxylationWorkGraph.build(
        structure=structure,
        surface_params=surface_params,
        code=code,
        builder_inputs=builder_inputs,
        bulk_structure=bulk_structure,                  # NEW
        bulk_builder_inputs=bulk_builder_inputs,        # NEW
        max_parallel_jobs=max_parallel_jobs,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        fix_elements=fix_elements,
        structure_specific_builder_inputs=structure_specific_builder_inputs,
    )
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: pass bulk parameters to WorkGraph build"
```

---

## Task 9: Update organize_hydroxylation_results Helper

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:424-494`

**Step 1: Add reference_data to return dict**

Modify the return statement in `organize_hydroxylation_results()`:

```python
def organize_hydroxylation_results(workflow_node):
    """
    Organize hydroxylation workflow results into successful/failed/statistics.

    ...existing docstring...

    Returns:
        dict with:
            - successful_relaxations: list of dicts with successful results
            - failed_relaxations: list of dicts with failed results
            - statistics: dict with total/succeeded/failed counts
            - reference_data: dict with bulk and pristine reference calculations (NEW)
    """
    # Get outputs
    manifest = workflow_node.outputs.manifest.get_dict()
    structures = workflow_node.outputs.structures
    energies = workflow_node.outputs.energies

    # NEW: Extract reference data
    reference_data = {
        'bulk_structure_pk': workflow_node.outputs.bulk_structure.pk,
        'bulk_energy': workflow_node.outputs.bulk_energy.value,
        'pristine_structure_pk': workflow_node.outputs.pristine_structure.pk,
        'pristine_energy': workflow_node.outputs.pristine_energy.value,
    }

    variants = manifest['variants']
    successful = []
    failed = []

    # ... existing variant processing code ...

    return {
        'successful_relaxations': successful,
        'failed_relaxations': failed,
        'statistics': {
            'total': len(variants),
            'succeeded': len(successful),
            'failed': len(failed)
        },
        'reference_data': reference_data,  # NEW
    }
```

**Step 2: Commit**

```bash
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add reference_data to organize_hydroxylation_results"
```

---

## Task 10: Check/Add extract_total_energy_from_misc Utility

**Files:**
- Check: `teros/core/surface_hydroxylation/utils.py`
- Possibly create new function

**Step 1: Check if function exists**

```bash
grep -n "extract_total_energy" teros/core/surface_hydroxylation/utils.py
```

**Step 2a: If function exists, verify it works with misc output**

Read the function and verify it extracts from VaspWorkChain.outputs.misc

**Step 2b: If function doesn't exist, create it**

Add to `teros/core/surface_hydroxylation/utils.py`:

```python
from aiida import orm
from aiida_workgraph import task

@task()
def extract_total_energy_from_misc(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP WorkChain misc output.

    Args:
        misc: Dict node containing various VASP outputs including total_energies

    Returns:
        Float node with final total energy in eV
    """
    misc_dict = misc.get_dict()

    # Get total_energies array from misc
    if 'total_energies' not in misc_dict:
        raise ValueError("total_energies not found in misc output")

    total_energies = misc_dict['total_energies']

    if not total_energies:
        raise ValueError("total_energies array is empty")

    # Return final energy (last ionic step)
    final_energy = float(total_energies[-1])
    return orm.Float(final_energy)
```

**Step 3: Commit (if created new function)**

```bash
git add teros/core/surface_hydroxylation/utils.py
git commit -m "feat: add extract_total_energy_from_misc utility function"
```

---

## Task 11: Create Test Example Folder

**Files:**
- Create: `examples/hydroxylation_with_bulk_reference/`

**Step 1: Create example directory structure**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/bulk-pristine-reference
mkdir -p examples/hydroxylation_with_bulk_reference
cd examples/hydroxylation_with_bulk_reference
```

**Step 2: Copy ag3po4.cif from test location**

```bash
cp /home/thiagotd/git/fosfato/calculos/hydroxylation/ag3po4.cif .
```

**Step 3: Commit**

```bash
git add examples/hydroxylation_with_bulk_reference/ag3po4.cif
git commit -m "test: add ag3po4.cif for bulk reference example"
```

---

## Task 12: Create run_hydroxylation_v2.py Example Script

**Files:**
- Create: `examples/hydroxylation_with_bulk_reference/run_hydroxylation_v2.py`

**Step 1: Create example script**

Create file with complete example based on run_hydroxylation_v1.py but with bulk additions:

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
Surface Hydroxylation v2 - Ag3PO4 with Bulk and Pristine Reference Calculations

This example demonstrates the updated hydroxylation workflow with automatic
bulk and pristine slab relaxations for surface energy calculations (Section S2).

NEW in v2:
- Automatic bulk relaxation from CIF file (ISIF=3 for cell relaxation)
- Automatic pristine slab relaxation (reference γ₀)
- Complete reference data for surface thermodynamics analysis
- Four new outputs: bulk_structure, bulk_energy, pristine_structure, pristine_energy

Usage:
    source ~/envs/aiida/bin/activate
    python run_hydroxylation_v2.py
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

def main():
    """Run hydroxylation workflow with bulk and pristine reference calculations."""

    print("\n" + "="*70)
    print("SURFACE HYDROXYLATION v2 - Ag3PO4 with Reference Calculations")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   ✓ Daemon is running")

    # Load slab structure (use your actual slab structure PK or file)
    print("\n2. Loading slab structure...")
    # MODIFY THIS: Use your actual slab structure
    structure_file = Path("st2_ag3po4_110_12A.vasp")

    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        structure = orm.StructureData(ase=atoms)
        print(f"   ✓ Loaded from file: {structure_file}")
    else:
        # Or load from database
        structure_pk = 1234  # MODIFY THIS
        structure = orm.load_node(structure_pk)
        print(f"   ✓ Loaded from PK: {structure_pk}")

    print(f"   Composition: {structure.get_composition()}")
    print(f"   Structure PK: {structure.pk}")

    # Surface modification parameters
    print("\n3. Surface modification parameters:")
    surface_params = {
        'mode': 'combine',               # Hydroxylation + vacancies
        'species': 'O',                  # Target oxygen atoms
        'z_window': 0.5,                 # Surface detection window (Å)
        'which_surface': 'top',          # Modify top surface only
        'oh_dist': 0.98,                 # O-H bond distance (Å)
        'include_empty': False,          # Not needed - pristine calculated automatically
        'deduplicate_by_coverage': True, # Enable deduplication
        'coverage_bins': 3,              # Sample 3 coverages for quick test
    }

    print(f"   Mode: {surface_params['mode']}")
    print(f"   Coverage bins: {surface_params['coverage_bins']}")
    print(f"   Expected structures: ~3-5")

    # Slab VASP parameters (ISIF=2, no cell relaxation)
    print("\n4. Slab VASP Configuration (ISIF=2):")
    code_label = 'VASPGAM-6.5.0@lovelace-parexp'

    slab_builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': True,
                'LCHARG': True,
                'NCORE': 3,
                'KPAR': 4,
                'ISIF': 2,              # Ionic relaxation only
                'NSW': 500,
                'IBRION': 2,
                'EDIFFG': -0.1,
            }
        },
        'kpoints_spacing': 1.0,
        'potential_family': 'PBE',
        'potential_mapping': {
            'Ag': 'Ag',
            'P': 'P',
            'O': 'O',
            'H': 'H',
        },
        'options': {
            'resources': {
                'num_machines': 2,
                'num_cores_per_machine': 48,
            },
            'queue_name': 'parexp',
        },
        'clean_workdir': False,
    }

    print(f"   VASP code: {code_label}")
    print(f"   ISIF: {slab_builder_inputs['parameters']['incar']['ISIF']} (ionic only)")
    print(f"   ENCUT: {slab_builder_inputs['parameters']['incar']['ENCUT']} eV")

    # NEW: Bulk VASP parameters (ISIF=3, cell relaxation)
    print("\n5. Bulk VASP Configuration (ISIF=3):")

    bulk_builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-6,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Fast',
                'LREAL': False,
                'LWAVE': False,
                'LCHARG': False,
                'ISIF': 3,              # Cell + ionic relaxation
                'NSW': 500,
                'IBRION': 2,
                'EDIFFG': -0.01,        # Tighter convergence
            }
        },
        'kpoints_spacing': 0.3,  # Denser k-points
        'potential_family': 'PBE',
        'potential_mapping': {
            'Ag': 'Ag',
            'P': 'P',
            'O': 'O',
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 48,
            },
            'queue_name': 'parexp',
        },
        'clean_workdir': False,
    }

    print(f"   ISIF: {bulk_builder_inputs['parameters']['incar']['ISIF']} (cell + ionic)")
    print(f"   ENCUT: {bulk_builder_inputs['parameters']['incar']['ENCUT']} eV")
    print(f"   K-points: {bulk_builder_inputs['kpoints_spacing']} Å⁻¹ (denser than slab)")

    # Parallelization and fixing
    max_parallel = 5
    fix_type = 'bottom'
    fix_thickness = 5.0

    print(f"\n6. Workflow Configuration:")
    print(f"   Batch size: {max_parallel} structures")
    print(f"   Fix type: {fix_type}")
    print(f"   Fix thickness: {fix_thickness} Å")

    # Build workflow with bulk reference
    print("\n7. Building workflow with bulk and pristine reference...")

    wg = build_surface_hydroxylation_workgraph(
        structure=structure,
        surface_params=surface_params,
        code_label=code_label,
        builder_inputs=slab_builder_inputs,

        # NEW: Bulk structure and parameters
        bulk_cif_path='ag3po4.cif',              # Option 1: CIF file
        # bulk_structure_pk=1234,                # Option 2: StructureData PK
        bulk_builder_inputs=bulk_builder_inputs,

        max_parallel_jobs=max_parallel,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        name='Ag3PO4_Hydroxylation_v2_with_references',
    )

    print("   ✓ Workflow built successfully")

    # Submit
    print("\n8. Submitting to AiiDA daemon...")
    result = wg.submit()
    pk = result.pk

    # Write PK to file
    pk_file = Path("workflow_pk.txt")
    with open(pk_file, 'w') as f:
        f.write(f"{pk}\n")
    print(f"   ✓ Workflow PK written to: {pk_file}")

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkflow PK: {pk}")

    print(f"\nMonitor with:")
    print(f"  verdi process show {pk}")
    print(f"  verdi process report {pk}")

    print(f"\nExpected workflow steps:")
    print(f"  1. Bulk relaxation (ISIF=3, cell optimization)")
    print(f"  2. Pristine slab relaxation (reference γ₀)")
    print(f"  3. Generate structure variants (~{surface_params['coverage_bins']})")
    print(f"  4. Relax all variants (batch={max_parallel})")

    print(f"\nNEW outputs in v2:")
    print(f"  - bulk_structure, bulk_energy")
    print(f"  - pristine_structure, pristine_energy")
    print(f"  - Plus existing: manifest, structures, energies")

    print(f"\nAfter completion, analyze results:")
    print(f"  python -c \"")
    print(f"from aiida import orm")
    print(f"from teros.core.surface_hydroxylation import organize_hydroxylation_results")
    print(f"node = orm.load_node({pk})")
    print(f"results = organize_hydroxylation_results(node)")
    print(f"ref = results['reference_data']")
    print(f"print('Bulk energy:', ref['bulk_energy'], 'eV')")
    print(f"print('Pristine energy:', ref['pristine_energy'], 'eV')")
    print(f"print('Successful variants:', len(results['successful_relaxations']))")
    print(f"  \"")

    print(f"\n{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

**Step 2: Make script executable**

```bash
chmod +x examples/hydroxylation_with_bulk_reference/run_hydroxylation_v2.py
```

**Step 3: Commit**

```bash
git add examples/hydroxylation_with_bulk_reference/run_hydroxylation_v2.py
git commit -m "test: add run_hydroxylation_v2.py example with bulk reference"
```

---

## Task 13: Create README for Example

**Files:**
- Create: `examples/hydroxylation_with_bulk_reference/README.md`

**Step 1: Create README**

```markdown
# Hydroxylation Workflow with Bulk and Pristine Reference Calculations

This example demonstrates the **v2 hydroxylation workflow** with automatic bulk and pristine slab reference calculations for surface energy analysis.

## What's New in v2

- **Automatic bulk relaxation** from CIF file (ISIF=3 for cell optimization)
- **Automatic pristine slab relaxation** for reference surface energy (γ₀)
- **Four new outputs**: `bulk_structure`, `bulk_energy`, `pristine_structure`, `pristine_energy`
- Complete reference data for surface thermodynamics (Section S2)

## Requirements

- Relaxed slab structure (StructureData or VASP file)
- Bulk crystal structure (CIF file, StructureData, or PK)
- VASP code configured in AiiDA
- AiiDA daemon running

## Files

- `run_hydroxylation_v2.py`: Main workflow script with bulk reference
- `ag3po4.cif`: Bulk Ag₃PO₄ structure for reference calculations
- `README.md`: This file

## Usage

### 1. Prepare Structures

Ensure you have:
- Your relaxed slab structure (modify PK in script)
- Bulk structure CIF file (provided: `ag3po4.cif`)

### 2. Configure Parameters

Edit `run_hydroxylation_v2.py`:

```python
# Slab structure (choose one)
structure_pk = 1234  # Your slab PK
# OR
structure_file = Path("your_slab.vasp")

# Bulk structure (choose one)
bulk_cif_path='ag3po4.cif'
# OR
bulk_structure_pk=5678
# OR
bulk_structure=orm.load_node(5678)
```

### 3. Run Workflow

```bash
source ~/envs/aiida/bin/activate
cd /home/thiagotd/git/PS-TEROS/examples/hydroxylation_with_bulk_reference
python run_hydroxylation_v2.py
```

### 4. Monitor Progress

```bash
# Get workflow PK from output or workflow_pk.txt
verdi process show <PK>
verdi process report <PK>

# Watch live updates
watch -n 30 verdi process show <PK>
```

## Expected Workflow Steps

1. **Bulk relaxation** (ISIF=3): Cell + ionic optimization of ag3po4.cif
2. **Pristine slab relaxation** (ISIF=2): Ionic relaxation of unmodified slab
3. **Generate variants**: Create hydroxylated/vacancy structures
4. **Relax variants**: Parallel VASP relaxations

## Outputs

After completion, access results:

```python
from aiida import orm
from teros.core.surface_hydroxylation import organize_hydroxylation_results

node = orm.load_node(<workflow_pk>)
results = organize_hydroxylation_results(node)

# Reference calculations
ref = results['reference_data']
print(f"Bulk energy: {ref['bulk_energy']} eV")
print(f"Pristine energy: {ref['pristine_energy']} eV")
print(f"Bulk PK: {ref['bulk_structure_pk']}")
print(f"Pristine PK: {ref['pristine_structure_pk']}")

# Variant results
for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV (coverage={r['coverage']:.2f})")
```

## Surface Energy Calculations

The reference data enables surface free energy calculations per Section S2:

- **E_bulk**: Use `ref['bulk_energy']`
- **γ₀**: Calculate from `ref['pristine_energy']` using equations 4-10
- **γ_modified**: Calculate for each variant using bulk and pristine references

## Comparison to v1

**v1 (run_hydroxylation_v1.py)**:
- Only relaxes modified surface structures
- No bulk or pristine reference calculations
- Cannot calculate surface energies

**v2 (run_hydroxylation_v2.py)**:
- Automatic bulk relaxation (ISIF=3)
- Automatic pristine slab relaxation
- Complete reference data for thermodynamics
- Ready for surface phase diagram generation

## Troubleshooting

**Error: "Bulk structure is required"**
- Provide one of: `bulk_cif_path`, `bulk_structure_pk`, or `bulk_structure`

**Error: "Could not read CIF file"**
- Check file exists: `ls ag3po4.cif`
- Check file format is valid CIF

**Bulk relaxation fails**
- Check `bulk_builder_inputs` has ISIF=3
- Verify VASP parameters are suitable for your system
- Check scheduler resources are adequate

## Related Documentation

- Design document: `/home/thiagotd/git/fosfato/calculos/hydroxylation/docs/plans/2025-10-27-bulk-pristine-reference-calculations-design.md`
- Section S2: `surface_energy_calc_procedure.tex`
- Main module: `teros/core/surface_hydroxylation/workgraph.py`
```

**Step 2: Commit**

```bash
git add examples/hydroxylation_with_bulk_reference/README.md
git commit -m "docs: add README for hydroxylation v2 example"
```

---

## Task 14: Test Import and Syntax

**Files:**
- Test: `teros/core/surface_hydroxylation/workgraph.py`

**Step 1: Clear Python cache**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/bulk-pristine-reference
find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete 2>/dev/null
```

**Step 2: Test import**

```bash
source ~/envs/aiida/bin/activate
python -c "from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph; print('✓ Import successful')"
```

Expected: "✓ Import successful" (no errors)

**Step 3: If errors, fix syntax/imports**

Check error messages and fix:
- Missing imports
- Syntax errors
- Typos

**Step 4: Commit any fixes**

```bash
git add <fixed_files>
git commit -m "fix: resolve import/syntax errors"
```

---

## Task 15: Restart AiiDA Daemon

**Files:**
- None (daemon restart)

**Step 1: Restart daemon to load new code**

```bash
verdi daemon restart
```

**Step 2: Wait for daemon to start**

```bash
sleep 5
verdi daemon status
```

Expected: Daemon running with workers

**Step 3: Verify profile**

```bash
verdi profile list
verdi status
```

Expected: All services running

---

## Task 16: Dry-Run Test (Optional but Recommended)

**Files:**
- Test: Build workflow without submission

**Step 1: Create test script**

Create `examples/hydroxylation_with_bulk_reference/test_build.py`:

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
Test building the v2 workflow without submission.
"""

from aiida import orm, load_profile
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

load_profile(profile='presto')

# Use minimal test structures
structure = orm.load_node(1234)  # Your test slab PK

try:
    wg = build_surface_hydroxylation_workgraph(
        structure=structure,
        surface_params={'mode': 'hydrogen', 'coverage_bins': 1},
        code_label='VASPGAM-6.5.0@lovelace-parexp',
        builder_inputs={'parameters': {'incar': {'ENCUT': 400}}},
        bulk_cif_path='ag3po4.cif',
        bulk_builder_inputs={'parameters': {'incar': {'ISIF': 3, 'ENCUT': 400}}},
        max_parallel_jobs=1,
    )
    print("✓ Workflow built successfully")
    print(f"✓ WorkGraph: {wg.name}")
    print(f"✓ Tasks: {len(wg.tasks)}")
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
```

**Step 2: Run test**

```bash
python test_build.py
```

Expected: "✓ Workflow built successfully"

**Step 3: If successful, remove test script (optional)**

```bash
rm test_build.py
```

---

## Task 17: Final Commit and Summary

**Files:**
- None

**Step 1: Check git status**

```bash
git status
git log --oneline -10
```

**Step 2: Verify all changes committed**

Ensure no uncommitted changes remain

**Step 3: Create summary commit message**

If any documentation updates needed:

```bash
git commit --allow-empty -m "feat: complete bulk and pristine reference calculations implementation

Summary of changes:
- Added bulk_structure and bulk_builder_inputs parameters
- Implemented bulk relaxation task (ISIF=3)
- Implemented pristine slab relaxation task
- Added four new outputs: bulk_structure/energy, pristine_structure/energy
- Added CIF file loading support
- Extended organize_hydroxylation_results with reference_data
- Created v2 example with ag3po4.cif
- Full documentation in examples/hydroxylation_with_bulk_reference/

Enables surface energy calculations per Section S2.
Ready for production testing."
```

---

## Testing Plan

After implementation, test the workflow:

### Test 1: CIF File Input

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/bulk-pristine-reference/examples/hydroxylation_with_bulk_reference
python run_hydroxylation_v2.py
```

**Success criteria:**
- Workflow submits successfully
- Bulk relaxation starts (check with `verdi process show <pk>`)
- Pristine relaxation starts after bulk
- All tasks complete with exit status 0
- Four new outputs present

### Test 2: StructureData PK Input

Modify `run_hydroxylation_v2.py`:

```python
# Replace bulk_cif_path with:
bulk_structure_pk=<existing_bulk_pk>,
```

Run and verify same success criteria.

### Test 3: Validation Error Handling

Test error messages:

```python
# Missing bulk structure
wg = build_surface_hydroxylation_workgraph(
    structure=structure,
    # No bulk inputs
)
# Expected: ValueError with helpful message

# Multiple bulk inputs
wg = build_surface_hydroxylation_workgraph(
    structure=structure,
    bulk_cif_path='ag3po4.cif',
    bulk_structure_pk=1234,  # Both provided
)
# Expected: ValueError about multiple inputs
```

### Test 4: Results Organization

After workflow completes:

```python
from teros.core.surface_hydroxylation import organize_hydroxylation_results
results = organize_hydroxylation_results(node)

assert 'reference_data' in results
assert 'bulk_energy' in results['reference_data']
assert 'pristine_energy' in results['reference_data']
print("✓ Reference data present")
```

---

## Merge Checklist

Before merging to develop:

- [ ] All commits follow conventional format
- [ ] No uncommitted changes
- [ ] Import test passes
- [ ] At least one full workflow test completed successfully
- [ ] All four new outputs (bulk/pristine structure/energy) verified
- [ ] organize_hydroxylation_results returns reference_data
- [ ] Example README is clear and accurate
- [ ] No broken imports or syntax errors
- [ ] Daemon restarted after changes

---

## Post-Merge: Update CLAUDE.md

After merging to develop, update the main CLAUDE.md:

```markdown
## Recent Updates

### Hydroxylation Module v2 (2025-10-27)

The surface hydroxylation module now includes automatic bulk and pristine slab reference calculations:

- Bulk relaxation with ISIF=3 (cell optimization)
- Pristine slab relaxation for γ₀ reference
- New parameters: bulk_structure, bulk_cif_path, bulk_structure_pk, bulk_builder_inputs
- New outputs: bulk_structure, bulk_energy, pristine_structure, pristine_energy
- Enables surface energy calculations per Section S2

See: `examples/hydroxylation_with_bulk_reference/` for usage.
```

---

**END OF IMPLEMENTATION PLAN**
