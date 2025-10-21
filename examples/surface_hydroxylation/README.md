# Surface Hydroxylation Module - Examples

This directory contains examples demonstrating the `surface_hydroxylation` module for generating and relaxing hydroxylated/oxygen-deficient surface variants.

## Overview

The surface_hydroxylation module provides:
- **Structure Generation**: Create surface variants with different OH coverages and O vacancies
- **Parallel Relaxation**: Relax multiple structures with controlled concurrency
- **Result Collection**: Organize results by coverage and success status

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
- `outputs.statistics` shows total structures generated
- `outputs.successful_relaxations` contains relaxed structures and energies

**Checking outputs:**
```python
from aiida import orm

# Load workflow node
node = orm.load_node(<PK>)

# Check statistics
stats = node.outputs.statistics.get_dict()
print(f"Total: {stats['total']}")
print(f"Succeeded: {stats['succeeded']}")
print(f"Failed: {stats['failed']}")

# Check successful relaxations
results = node.outputs.successful_relaxations.get_dict()['results']
for result in results:
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

## Module Architecture

### Workflow Structure

```
SurfaceHydroxylationWorkGraph (main workflow)
  │
  ├─> generate_structures (CalcFunction)
  │     ├─ Converts AiiDA → ASE
  │     ├─ Runs SurfaceModifier
  │     ├─ Generates structure variants
  │     └─ Returns manifest + structures
  │
  ├─> relax_slabs_with_semaphore (child WorkGraph)
  │     ├─ Scatters: Creates VASP task for each structure
  │     ├─ Limits: Max parallel via semaphore
  │     ├─ Gathers: Collects structures, energies, exit_statuses
  │     └─ Returns indexed results
  │
  └─> collect_results (CalcFunction)
        ├─ Matches structures to manifest
        ├─ Separates successful vs failed
        ├─ Calculates statistics
        └─ Returns organized results
```

### Key Components

**generate_structures** (`teros/core/surface_hydroxylation/tasks.py`)
- CalcFunction wrapping SurfaceModifier
- Input: Single structure + parameters
- Output: Manifest + multiple structures (structure_0, structure_1, ...)

**relax_slabs_with_semaphore** (`teros/core/surface_hydroxylation/relaxations.py`)
- Child WorkGraph for parallel relaxations
- Uses VASP workchain (`vasp.v2.vasp`)
- Returns: structures, energies, exit_statuses, errors (all indexed)

**collect_results** (`teros/core/surface_hydroxylation/tasks.py`)
- CalcFunction for result organization
- Groups by success/failure
- Stores structure PKs (not full structures) for efficiency

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

### builder_config (Dict)

VASP builder configuration:

```python
builder_config = {
    'code': orm.load_code('VASP-VTST-6.4.3@bohr'),
    'parameters': {
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 200,
        'IBRION': 2,
        'ISIF': 2,
    },
    'potential_family': 'PBE',
    'potential_mapping': {},  # Optional element-specific mapping
    'kpoints_spacing': 0.3,   # K-points spacing (Å⁻¹)
    'options': {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 40},
        'queue_name': 'par40',
        'max_wallclock_seconds': 3600 * 10,
    },
    'clean_workdir': False,
    'settings': None,  # Optional parser settings
}
```

### max_parallel_jobs (Int)

Maximum number of concurrent VASP relaxations.

**Guidelines:**
- Consider cluster capacity and queue limits
- Typical: 2-5 for testing, 10-20 for production
- Each job uses resources specified in `builder_config.options`

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

**Some relaxations fail:**
- Check individual job reports: `verdi process report <relax_PK>`
- Common: Convergence issues with extreme configurations
- Check failed_relaxations output for patterns

**Semaphore not limiting:**
- Note: Current implementation relies on WorkGraph's natural concurrency
- For explicit limits, consider external job management

**Jobs stuck in queue:**
- Check queue status: `squeue -u $USER`
- Verify walltime is sufficient
- Check cluster load

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

## References

**Module location:** `teros/core/surface_hydroxylation/`

**Key files:**
- `workgraph.py` - Main workflow definition
- `tasks.py` - CalcFunctions (generate_structures, collect_results)
- `relaxations.py` - Parallel relaxation WorkGraph
- `surface_modes.py` - Structure modification logic
- `utils.py` - Helper functions

**Documentation:**
- Implementation plan: `docs/plans/2025-10-21-surface-hydroxylation.md`
- Module tests: `tests/core/surface_hydroxylation/`

## Support

For issues or questions:
1. Check this README troubleshooting section
2. Review module tests for examples
3. Check AiiDA provenance for error details
4. Contact PS-TEROS development team
