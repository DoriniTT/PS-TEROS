# Surface Hydroxylation Module - Examples

This directory contains examples demonstrating the `surface_hydroxylation` module for generating and relaxing hydroxylated/oxygen-deficient surface variants.

## Overview

The surface_hydroxylation module provides:
- **Structure Generation**: Create surface variants with different OH coverages and O vacancies
- **Parallel Relaxation**: Relax multiple structures with batch control
- **Result Collection**: Organize results by coverage and success status

**Key Features:**
- Builder function pattern for easy workflow creation
- Direct VASP INCAR parameter configuration
- Post-processing helper for result organization
- Complete provenance tracking with AiiDA

## Quick Start

```python
from teros.core.surface_hydroxylation import (
    build_surface_hydroxylation_workgraph,
    organize_hydroxylation_results,
)

# Build workflow
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,                    # Your relaxed slab
    surface_params={'mode': 'hydrogen'},  # Hydroxylation
    code_label='VASP-6.4.1@cluster',
    vasp_config={
        'parameters': {'ENCUT': 520, 'ISIF': 2, 'EDIFFG': -0.02},
        'kpoints_spacing': 0.3,
        'potential_family': 'PBE',
    },
    options={'resources': {'num_machines': 1}},
    max_parallel_jobs=3,
)

# Submit
result = wg.submit()

# After completion: organize results
from aiida import orm
node = orm.load_node(result.pk)
results = organize_hydroxylation_results(node)
print(results['statistics'])  # Total, succeeded, failed
```

See sections below for detailed examples and parameter descriptions.

## Examples

### 1. Structure Generation Test (Recommended for Testing)

**File:** `test_structure_generation.py`

Tests ONLY the structure generation component without submitting VASP jobs.

**Features:**
- Fast execution (< 1 minute)
- Tests `generate_structures` CalcFunction
- Verifies hydroxylation mode works correctly
- Checks coverage-based deduplication
- Good for development and debugging

**Test System:**
- 2×2 Pt(111) slab with O adlayer (~20 atoms)
- Hydrogen mode (hydroxylation only)
- 3 coverage bins (generates ~3-5 structures)

**Running:**
```bash
# Ensure AiiDA daemon is running
verdi daemon status

# Run test
/home/thiagotd/envs/aiida/bin/python test_structure_generation.py
```

**Expected Output:**
```
======================================================================
SURFACE HYDROXYLATION - STRUCTURE GENERATION TEST
======================================================================

1. Loading AiiDA profile...
   ✓ Profile loaded: psteros

2. Creating test structure...
   - Pt slab: 16 atoms
   - Total atoms: 20
   - O atoms: 4
   ✓ Structure created

3. Setting up surface generation parameters...
   Parameters:
   - Mode: hydrogen (hydroxylation)
   - Species: O
   - Coverage bins: 3
   - Deduplication: enabled

4. Running structure generation...
   (This is a CalcFunction - will be stored in provenance)

5. Analyzing results...
   ✓ Manifest created
   - Total variants generated: 3

   Variant details:
   [0] oh_001_0.25
       - Coverage: 0.25
       - Atoms: 21
       - PK: 12345
   [1] oh_002_0.50
       - Coverage: 0.50
       - Atoms: 22
       - PK: 12346
   [2] oh_003_0.75
       - Coverage: 0.75
       - Atoms: 23
       - PK: 12347

6. Verification...
   ✓ All 3 structures created successfully
   ✓ Coverage range: 0.25 - 0.75

======================================================================
STRUCTURE GENERATION TEST COMPLETED SUCCESSFULLY
======================================================================
```

**What it Tests:**
- ✓ Module imports work correctly
- ✓ Structure conversion (AiiDA ↔ ASE)
- ✓ SurfaceModifier integration
- ✓ CalcFunction provenance storage
- ✓ Coverage-based deduplication
- ✓ Multiple structure output handling

### 2. Full Workflow Test (Production-Like)

**File:** `test_full_workflow.py`

Tests the COMPLETE workflow including VASP relaxations.

**WARNING:** This submits actual VASP jobs! Only run with proper VASP setup.

**Features:**
- Complete workflow test
- Includes parallel VASP relaxations
- Tests semaphore-based concurrency control
- Production-like validation

**Test System:**
- Same 2×2 Pt(111) with O adlayer
- Lightweight VASP parameters (NSW=10, ENCUT=400)
- Max 2 parallel jobs
- Expected runtime: ~10-20 minutes

**Prerequisites:**
- VASP code configured in AiiDA
- Pseudopotentials installed
- Compute resources available

**Running:**
```bash
# Check VASP setup
verdi code show VASP-VTST-6.4.3@bohr

# Check daemon
verdi daemon status

# Run test
/home/thiagotd/envs/aiida/bin/python test_full_workflow.py
```

**Monitoring:**
```bash
# Get workflow PK from output
verdi process show <PK>

# Monitor progress (refresh every 30s)
watch -n 30 verdi process show <PK>

# Check detailed report
verdi process report <PK>
```

**Expected Results:**

**Success criteria:**
- Main workflow exits with `[0]`
- All tasks complete without errors
- `outputs.manifest` shows total structures generated
- `outputs.structures` and `outputs.energies` contain raw namespace outputs

**Checking outputs (using post-processing helper):**
```python
from aiida import orm
from teros.core.surface_hydroxylation import organize_hydroxylation_results

# Load workflow node
node = orm.load_node(<PK>)

# Organize results using helper function
results = organize_hydroxylation_results(node)

# Check statistics
print(f"Total: {results['statistics']['total']}")
print(f"Succeeded: {results['statistics']['succeeded']}")
print(f"Failed: {results['statistics']['failed']}")

# Check successful relaxations
for result in results['successful_relaxations']:
    name = result['name']
    energy = result['energy']
    coverage = result['coverage']
    structure_pk = result['structure_pk']

    print(f"{name}:")
    print(f"  Coverage: {coverage:.2f}")
    print(f"  Energy: {energy:.6f} eV")
    print(f"  Structure PK: {structure_pk}")

    # Load relaxed structure
    relaxed = orm.load_node(structure_pk)
```

**Note:** The workflow returns raw namespace outputs (manifest, structures, energies).
Use `organize_hydroxylation_results()` helper function to organize results into
successful/failed categories with statistics.

### 3. Production Example

**File:** `run_production.py`

Production-ready example for real research calculations.

**Features:**
- Complete workflow with production VASP settings
- Combined mode (vacancies + hydroxylation)
- 10 coverage bins (~15-20 structures)
- Parallel execution with controlled concurrency
- Full provenance tracking

**System Requirements:**
- Realistic perovskite oxide surface (~100+ atoms)
- Production VASP settings (converged ENCUT, k-points)
- Cluster compute resources configured
- VASP pseudopotentials installed

**Setup:**
1. Edit `run_production.py` and replace placeholders:
   - Set structure PK from surface_thermodynamics output
   - Update VASP code label for your cluster
   - Update compute resources (account, queue, cores)
   - Adjust VASP parameters for your system

2. Review parameters:
```python
# Surface modification
surface_params = {
    'mode': 'combine',           # Both vacancies and hydroxylation
    'coverage_bins': 10,         # Good coverage sampling
    'deduplicate_by_coverage': True
}

# VASP configuration (direct INCAR parameters)
vasp_config = {
    'parameters': {              # Direct INCAR parameters
        'PREC': 'Accurate',
        'ENCUT': 520,            # Converged cutoff
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'ISIF': 2,               # Relax positions only
        'NSW': 200,              # Max ionic steps
        'IBRION': 2,             # CG relaxation
        'EDIFFG': -0.02,         # Force convergence (eV/Å)
    },
    'kpoints_spacing': 0.3,      # Converged k-points (Å⁻¹)
    'potential_family': 'PBE',
    'clean_workdir': False,
}

# Scheduler options
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 40,
    },
    'queue_name': 'par40',
    'max_wallclock_seconds': 3600 * 10,
}

# Parallelization
max_parallel_jobs = 5            # Adjust for cluster capacity
```

3. Submit using builder function:
```bash
/home/thiagotd/envs/aiida/bin/python run_production.py
```

The script uses `build_surface_hydroxylation_workgraph()` builder function which
handles code loading and workflow creation automatically.

**Expected Runtime:**
- Structure generation: < 1 minute
- VASP relaxations: Several hours (depends on system size)
- Total: ~6-12 hours for typical perovskite surface

**Monitoring:**
```bash
# Check overall workflow status
verdi process show <WORKFLOW_PK>

# Monitor running VASP jobs
verdi process list -a -p 1  # Last 24 hours

# Check specific relaxation
verdi process report <RELAX_PK>
```

**Post-Processing:**
```python
from aiida import orm
from teros.core.surface_hydroxylation import organize_hydroxylation_results

# Load completed workflow
wf = orm.load_node(<WORKFLOW_PK>)

# Organize results using helper function
results = organize_hydroxylation_results(wf)

# Get statistics
stats = results['statistics']
print(f"Generated {stats['total']} structures")
print(f"Successfully relaxed {stats['succeeded']}")
print(f"Failed {stats['failed']}")

# Analyze successful relaxations
successful = results['successful_relaxations']

# Find lowest energy configuration
results_sorted = sorted(successful, key=lambda r: r['energy'])
best = results_sorted[0]

print(f"\nLowest energy configuration:")
print(f"  Name: {best['name']}")
print(f"  Coverage: {best['coverage']:.2f}")
print(f"  Energy: {best['energy']:.6f} eV")
print(f"  Structure PK: {best['structure_pk']}")

# Load structure for further analysis
best_structure = orm.load_node(best['structure_pk'])
```

**Important:** The workflow returns raw namespace outputs. Always use
`organize_hydroxylation_results()` helper to process results after completion.

**Important Notes:**
- Review VASP settings for your specific system
- Test with small system first (`test_full_workflow.py`)
- Monitor disk space (VASP outputs can be large)
- Plan for cluster downtime (workflows can be restarted)

## Module Architecture

### Builder Function Pattern

The module follows PS-TEROS builder function pattern with two-tier design:

**High-level (User-facing):**
```python
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,              # or structure=orm.load_node(1234)
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    vasp_config=vasp_config,
    options=options,
    max_parallel_jobs=3,
    name='SurfaceHydroxylation',
)
result = wg.submit()
```

**Low-level (@task.graph):**
```python
@task.graph(outputs=['manifest', 'structures', 'energies'])
def SurfaceHydroxylationWorkGraph(...):
    # Actual workflow implementation
    pass
```

### Workflow Structure

```
build_surface_hydroxylation_workgraph (builder function)
  │
  ├─ Input validation (structure, code, parameters)
  ├─ Set defaults for optional parameters
  └─ Build and return WorkGraph
        │
        └─> SurfaceHydroxylationWorkGraph (@task.graph)
              │
              ├─> generate_structures (CalcFunction)
              │     ├─ Converts AiiDA → ASE
              │     ├─ Runs SurfaceModifier
              │     ├─ Generates structure variants
              │     └─ Returns namespace: manifest + structures dict
              │
              ├─> relax_slabs_with_semaphore (child @task.graph)
              │     ├─ Scatters: Creates VASP task for each structure
              │     ├─ Batch control: max_parallel_jobs limit
              │     ├─ Gathers: Collects structures, energies
              │     └─ Returns namespace: structures dict + energies dict
              │
              └─> Returns raw namespace outputs
                    ├─ manifest (Dict)
                    ├─ structures (namespace dict)
                    └─ energies (namespace dict)

organize_hydroxylation_results (post-processing helper)
  │
  ├─ Input: Completed workflow node
  ├─ Extracts: manifest, structures, energies namespaces
  ├─ Processes: Matches results, separates success/failure
  └─ Returns: {successful_relaxations, failed_relaxations, statistics}
```

### Key Components

**build_surface_hydroxylation_workgraph** (`teros/core/surface_hydroxylation/workgraph.py`)
- Builder function (user-facing API)
- Validates inputs, sets defaults
- Loads structure and code
- Returns ready-to-submit WorkGraph

**SurfaceHydroxylationWorkGraph** (`teros/core/surface_hydroxylation/workgraph.py`)
- Main workflow (@task.graph decorator)
- Orchestrates structure generation and relaxation
- Returns raw namespace outputs (manifest, structures, energies)

**generate_structures** (`teros/core/surface_hydroxylation/tasks.py`)
- CalcFunction wrapping SurfaceModifier
- Input: Single structure + parameters
- Output: Namespace with manifest + structures dict (indexed by '0', '1', ...)

**relax_slabs_with_semaphore** (`teros/core/surface_hydroxylation/relaxations.py`)
- Child @task.graph for parallel relaxations
- Uses VASP workchain (`vasp.v2.vasp`)
- Batch control via max_parallel parameter
- Returns: Namespace with structures dict + energies dict (indexed)

**organize_hydroxylation_results** (`teros/core/surface_hydroxylation/workgraph.py`)
- Post-processing helper (runs after workflow completion)
- Input: Completed workflow node
- Extracts and organizes raw namespace outputs
- Returns: Dict with successful_relaxations, failed_relaxations, statistics
- **Must be used** to get organized results from workflow

## Parameters

### surface_params (Dict)

Parameters for structure generation:

```python
surface_params = {
    'mode': 'hydrogen',              # 'hydrogen', 'vacancies', or 'combine'
    'species': 'O',                  # Target species for modification
    'z_window': 0.5,                 # Z-coordinate window for surface detection (Å)
    'which_surface': 'top',          # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,                 # O-H bond distance for hydroxylation (Å)
    'include_empty': False,          # Include structure with no modifications
    'supercell': None,               # Supercell expansion [nx, ny, nz] or None
    'deduplicate_by_coverage': True, # Enable coverage-based deduplication
    'coverage_bins': 5,              # Number of coverage bins for sampling
}
```

**Mode options:**
- `hydrogen`: Add OH groups to O sites (hydroxylation)
- `vacancies`: Remove O atoms (oxygen vacancies)
- `combine`: Both hydroxylation and vacancies

**Coverage deduplication:**
- Reduces thousands of structures to ~N bins worth
- Each bin samples representative configurations
- Typical: `coverage_bins=5` gives ~5-10 structures

### vasp_config (Dict)

VASP configuration with direct INCAR parameters:

```python
vasp_config = {
    # Direct INCAR parameters (passed to vasp.v2.vasp workflow)
    'parameters': {
        'PREC': 'Accurate',       # Precision level
        'ENCUT': 520,             # Plane-wave cutoff (eV)
        'EDIFF': 1e-6,            # Electronic convergence (eV)
        'ISMEAR': 0,              # Smearing method (Gaussian)
        'SIGMA': 0.05,            # Smearing width (eV)
        'ALGO': 'Normal',         # Electronic minimization algorithm
        'LREAL': False,           # Real-space projection
        'NELM': 100,              # Max electronic steps
        'LWAVE': False,           # Write WAVECAR
        'LCHARG': False,          # Write CHGCAR
        # Relaxation parameters
        'ISIF': 2,                # Relax positions only (2), or positions+cell (3)
        'NSW': 200,               # Max ionic steps
        'IBRION': 2,              # Ionic relaxation (2=CG, 1=RMM-DIIS)
        'EDIFFG': -0.02,          # Force convergence (eV/Å, negative = force criterion)
    },
    # K-points
    'kpoints_spacing': 0.3,       # K-points spacing (Å⁻¹)

    # Pseudopotentials
    'potential_family': 'PBE',    # Potential family name
    'potential_mapping': {},      # Optional element-specific mapping

    # Cleanup
    'clean_workdir': False,       # Keep calculation files
}
```

### options (Dict)

Scheduler options (separate from vasp_config):

```python
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 40,  # Or num_cores_per_machine
    },
    'queue_name': 'par40',
    'max_wallclock_seconds': 3600 * 10,
    # Optional:
    # 'account': 'my_account',
    # 'prepend_text': 'module load vasp/6.4.2',
    # 'custom_scheduler_commands': '#SBATCH --constraint=haswell',
}
```

**Note:** The `code_label` is passed separately to the builder function, not in vasp_config.

### max_parallel_jobs (Int)

Maximum number of structures to process in this workflow run (batch control).

**Guidelines:**
- Use batch approach: processes first N structures only
- Typical: 1-2 for testing, 5-10 for production batches
- Each job uses resources specified in `options`
- Increase this value in subsequent runs to process more structures

**Example workflow:**
```python
# First run: test with 2 structures
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=2)

# After verification: process more structures
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=10)
```

## Usage Patterns

### Basic Usage (Recommended)

```python
from aiida import orm
from teros.core.surface_hydroxylation import (
    build_surface_hydroxylation_workgraph,
    organize_hydroxylation_results,
)

# 1. Define parameters
surface_params = {
    'mode': 'hydrogen',
    'species': 'O',
    'coverage_bins': 5,
}

vasp_config = {
    'parameters': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISIF': 2,
        'NSW': 200,
        'IBRION': 2,
        'EDIFFG': -0.02,
    },
    'kpoints_spacing': 0.3,
    'potential_family': 'PBE',
}

options = {
    'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
    'queue_name': 'normal',
}

# 2. Build and submit workflow
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    vasp_config=vasp_config,
    options=options,
    max_parallel_jobs=3,
)
result = wg.submit()
pk = result.pk

# 3. Wait for completion
# Monitor: verdi process show <pk>

# 4. Process results after completion
node = orm.load_node(pk)
results = organize_hydroxylation_results(node)

print(f"Total: {results['statistics']['total']}")
print(f"Succeeded: {results['statistics']['succeeded']}")
print(f"Failed: {results['statistics']['failed']}")

# 5. Analyze successful relaxations
for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV (coverage={r['coverage']:.2f})")
```

### Advanced: Accessing Raw Outputs

If you need direct access to the workflow outputs:

```python
node = orm.load_node(pk)

# Raw namespace outputs
manifest = node.outputs.manifest.get_dict()
structures = node.outputs.structures  # Namespace dict
energies = node.outputs.energies      # Namespace dict

# Access specific structure
structure_0 = structures['0']
energy_0 = energies['0'].value

# For most use cases, use organize_hydroxylation_results() instead
```

### Incremental Processing (Batch Approach)

Process structures in batches to test before full production run:

```python
# Step 1: Test with first 2 structures
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    vasp_config=vasp_config,
    options=options,
    max_parallel_jobs=2,  # Process only first 2
)
result1 = wg.submit()

# Step 2: After verification, process more
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    vasp_config=vasp_config,
    options=options,
    max_parallel_jobs=10,  # Process more structures
)
result2 = wg.submit()
```

## Troubleshooting

### Structure Generation Issues

**No structures generated:**
- Check `z_window` - might be too restrictive
- Verify target species exists in surface region
- Try `include_empty=True` to get baseline structure

**Too many structures:**
- Enable `deduplicate_by_coverage=True`
- Reduce `coverage_bins` (e.g., 3-5 for testing)

**Wrong species modified:**
- Verify `species` parameter matches structure
- Check `which_surface` setting (top/bottom/both)

### VASP Relaxation Issues

**All relaxations fail:**
- Verify VASP code: `verdi code show <code_label>`
- Check pseudopotentials: `verdi data core.upf listfamilies`
- Validate compute resources in `options`
- Check VASP parameters in `vasp_config['parameters']`

**Some relaxations fail:**
- Check individual job reports: `verdi process report <relax_PK>`
- Common: Convergence issues with extreme configurations
- Use `organize_hydroxylation_results()` to see failed_relaxations list
- Review failed configurations to identify problematic structures

**Batch control not working:**
- Verify `max_parallel_jobs` parameter
- Check workflow outputs to confirm number of structures processed
- Note: Batch approach processes first N structures only

**Jobs stuck in queue:**
- Check queue status: `squeue -u $USER`
- Verify walltime is sufficient in `options`
- Check cluster load and queue limits

**Cannot access results:**
- Use `organize_hydroxylation_results()` helper function
- Do NOT try to access namespace dict outputs directly in CalcFunctions
- Workflow returns raw outputs: manifest, structures, energies

### Provenance Issues

**CalcFunction errors:**
- Clear cache: `find . -type d -name __pycache__ -exec rm -rf {} +`
- Restart daemon: `verdi daemon restart`
- Check imports work: `python -c "from teros.core.surface_hydroxylation import *"`

**Result node missing:**
- Ensure CalcFunction completed: `verdi process show <PK>`
- Check exit status is `[0]`
- Verify outputs exist: `verdi node show <PK>`

## Development Notes

### Testing Strategy

1. **Start with structure generation test** (this is fast)
2. Verify manifest and structures look correct
3. Only then run full workflow with VASP

### Production Usage

For real research:
- Use realistic surface structure (~100-200 atoms)
- Set `coverage_bins=10` for better sampling
- Use production VASP parameters (converged ENCUT, k-points)
- Set appropriate `max_parallel_jobs` for cluster

### Integration with PS-TEROS

This module integrates with:
- Surface thermodynamics workflow (provides input slabs)
- Electronic structure analysis (uses output structures)
- Custom workflow builder (can be added as step)

## Important Architectural Decisions

### Why Builder Function Pattern?

The module uses a two-tier design following PS-TEROS patterns:
- **Builder function** (`build_surface_hydroxylation_workgraph`): User-facing API
- **@task.graph function** (`SurfaceHydroxylationWorkGraph`): Internal implementation

This separation provides:
- Input validation and helpful error messages
- Sensible defaults for optional parameters
- Code loading handled automatically
- Consistent API across PS-TEROS modules

### Why Post-Processing Helper?

The workflow returns **raw namespace outputs** (manifest, structures, energies) because:
- **WorkGraph limitation**: Namespace dicts containing AiiDA nodes cannot be passed to CalcFunctions due to JSON serialization constraints
- **Solution**: Return raw namespaces + provide Python helper for post-processing

**DO NOT** try to organize results inside the workflow - use `organize_hydroxylation_results()` after completion.

### Why Direct INCAR Parameters?

VASP configuration uses direct INCAR parameters (`vasp_config['parameters']`) because:
- **Simplicity**: Matches standard `vasp.v2.vasp` workflow inputs
- **Clarity**: No confusing translation layer from nested relax/base structure
- **Flexibility**: Full control over all VASP INCAR tags

This module uses `vasp.v2.vasp` workflow (generic VASP), NOT `vasp.v2.relax` (specialized relaxation).

## References

**Module location:** `teros/core/surface_hydroxylation/`

**Key files:**
- `workgraph.py` - Builder function, main workflow, post-processing helper
- `tasks.py` - CalcFunctions (generate_structures)
- `relaxations.py` - Parallel relaxation @task.graph
- `surface_modes.py` - Structure modification logic (SurfaceModifier)
- `utils.py` - Helper functions (AiiDA ↔ ASE conversion)
- `__init__.py` - Module exports

**Exported API:**
```python
from teros.core.surface_hydroxylation import (
    build_surface_hydroxylation_workgraph,  # Builder function (use this)
    organize_hydroxylation_results,         # Post-processing helper (use this)
    SurfaceHydroxylationWorkGraph,          # Low-level @task.graph
    SurfaceModifier,                        # Structure modification class
)
```

**Documentation:**
- This README: Usage examples and troubleshooting
- Implementation plan: `docs/plans/2025-10-21-surface-hydroxylation.md`
- Module tests: `tests/core/surface_hydroxylation/`
- PS-TEROS module builder skill: `~/.claude/skills/psteros-module-builder/`

## Support

For issues or questions:
1. Check this README troubleshooting section
2. Review module tests for examples
3. Check AiiDA provenance for error details
4. Contact PS-TEROS development team
