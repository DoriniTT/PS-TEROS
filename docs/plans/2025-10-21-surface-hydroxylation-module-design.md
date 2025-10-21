# Surface Hydroxylation Module Design

**Date:** 2025-10-21
**Author:** Design session with Claude
**Status:** Design approved, ready for implementation

## Overview

The `surface_hydroxylation` module enables systematic exploration of hydroxylated and oxygen-deficient surface structures for computational surface chemistry studies. It generates surface variants with controlled hydroxylation (OH groups) and oxygen vacancies at different coverages, then relaxes all structures in parallel using VASP.

## Scope and Purpose

**In scope:**
- Take a relaxed surface slab as input
- Generate surface variants (hydroxylation + vacancies) using surface_modes.py
- Relax all variants with VASP using controlled parallelization
- Store relaxed structures and total energies in AiiDA database

**Out of scope (for now):**
- Thermodynamic analysis (phase diagrams, formation energies)
- Pourbaix/electrochemical analysis
- These will be added in future iterations

## Requirements Summary

| Requirement | Decision |
|------------|----------|
| **Module location** | `teros/core/surface_hydroxylation/` |
| **Module name** | `surface_hydroxylation` |
| **Script integration** | CalcFunction wrapper around surface_modes.py |
| **Parallelization** | WorkGraph awaitable/semaphore mechanism |
| **VASP configuration** | Single shared builder config for all structures |
| **Error handling** | Continue with partial results, report failures separately |
| **Output format** | Namespace with successful/failed relaxations + metadata |

## Architecture

### Design Pattern: Nested WorkGraph

```
SurfaceHydroxylationWorkGraph (Main)
│
├─► generate_structures (CalcFunction)
│   - Wraps surface_modes.py
│   - Returns: manifest Dict + List[StructureData]
│
├─► RelaxationsWorkGraph (Child WorkGraph)
│   - Creates N VaspRelaxWorkChain tasks
│   - Applies semaphore for max_parallel_jobs
│   - Returns: namespace {relax_0, relax_1, ..., relax_N}
│
└─► collect_results (CalcFunction)
    - Processes relaxation outputs
    - Returns: namespace {successful, failed, statistics}
```

**Rationale:** Nested WorkGraph pattern matches existing PS-TEROS conventions (AIMD, surface_thermodynamics) and provides good modularity while keeping components reusable.

## Module Structure

```
teros/core/surface_hydroxylation/
├── __init__.py
├── workgraph.py              # SurfaceHydroxylationWorkGraph
├── tasks.py                  # CalcFunctions (generate_structures, collect_results)
├── relaxations.py            # RelaxationsWorkGraph with semaphore
├── surface_modes.py          # Moved from teros/experimental/
└── utils.py                  # Helper functions
```

## Component Specifications

### 1. generate_structures (CalcFunction)

**Location:** `tasks.py`

**Signature:**
```python
@calcfunction
def generate_structures(structure: StructureData,
                       params: Dict) -> dict:
    """
    Wraps surface_modes.py to generate surface variants.

    Args:
        structure: Input relaxed slab structure
        params: Dict with surface_modes.py parameters:
            - mode: str ('vacancies'/'hydrogen'/'combine')
            - species: str (default 'O')
            - z_window: float (default 0.5)
            - which_surface: str ('top'/'bottom'/'both')
            - oh_dist: float (default 0.98)
            - include_empty: bool (default False)
            - supercell: List[int] or None ([nx, ny, nz])
            - deduplicate_by_coverage: bool
            - coverage_bins: int or None

    Returns:
        {
            'manifest': Dict (parsed manifest.json),
            'structures': List[StructureData] (all variants)
        }
    """
```

**Implementation notes:**
- Convert AiiDA StructureData → ASE Atoms
- Import and call SurfaceModifier from surface_modes.py
- Parse manifest.json output
- Convert generated structures → List[StructureData]
- Return both manifest and structures for downstream tasks

### 2. RelaxationsWorkGraph (Child WorkGraph)

**Location:** `relaxations.py`

**Signature:**
```python
class RelaxationsWorkGraph(WorkGraph):
    """Parallel VASP relaxations with semaphore limiting."""

    def setup(self,
              structures: List[StructureData],
              builder_config: Dict,
              max_parallel: int):
        """
        Creates N VaspRelaxWorkChain tasks with semaphore.

        Args:
            structures: List of surface variants to relax
            builder_config: Complete Dict to populate builder
                (metadata, kpoints_distance, parameters, relax, ...)
            max_parallel: Max concurrent jobs (semaphore limit)
        """
```

**Implementation notes:**
- Create semaphore context with `max_parallel` limit
- For each structure:
  - Create `VaspRelaxWorkChain` task with name `relax_{idx}`
  - Populate builder from `builder_config` Dict (direct pass-through)
  - Set structure input
  - Apply semaphore constraint: `task.waiting_on = semaphore`
  - Add to namespace output
- Semaphore ensures only `max_parallel` jobs run concurrently
- When one finishes, next waiting job starts automatically

### 3. collect_results (CalcFunction)

**Location:** `tasks.py`

**Signature:**
```python
@calcfunction
def collect_results(relaxations: namespace,
                   manifest: Dict) -> dict:
    """
    Collects and organizes relaxation results.

    Args:
        relaxations: Namespace {relax_0, relax_1, ..., relax_N}
        manifest: Original manifest from generate_structures

    Returns:
        namespace with:
            - successful_relaxations: List[Dict]
                [{structure, energy, coverage, name, metadata}, ...]
            - failed_relaxations: List[Dict]
                [{name, coverage, exit_status, error_message}, ...]
            - statistics: Dict {total, succeeded, failed}
    """
```

**Implementation notes:**
- Iterate through relaxations namespace by index
- Check `exit_status` for each relaxation
- **Success (exit_status == 0):**
  - Extract `structure` and `total_energy`
  - Match with manifest metadata using index
  - Add to `successful_relaxations` list
- **Failure (exit_status != 0):**
  - Extract error information
  - Add to `failed_relaxations` list
- Calculate statistics (total, success count, failure count)
- Return organized namespace

### 4. SurfaceHydroxylationWorkGraph (Main)

**Location:** `workgraph.py`

**User Interface:**
```python
from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph

wg = SurfaceHydroxylationWorkGraph()

# Inputs:
wg.input.structure = input_slab_structure      # StructureData
wg.input.builder_config = builder_dict         # Dict (complete builder)
wg.input.surface_params = params_dict          # Dict (surface_modes params)
wg.input.max_parallel_jobs = Int(5)           # Max concurrent jobs

# Outputs:
wg.outputs.successful_relaxations   # List[Dict] with structure, energy, coverage
wg.outputs.failed_relaxations       # List[Dict] with coverage, error info
wg.outputs.statistics              # Dict: {total, succeeded, failed}
wg.outputs.manifest                # Dict: original manifest from surface_modes
```

## Input Specifications

### builder_config Structure

User provides complete VASP builder configuration as Dict with full control:

```python
builder_dict = Dict({
    'metadata': {
        'options': {
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
            'max_wallclock_seconds': 3600 * 10,
            'account': 'project_account',
            'queue_name': 'normal'
        }
    },
    'kpoints_distance': 0.3,
    'parameters': {  # INCAR settings (direct pass-through)
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 200,
        'IBRION': 2,
        'ISIF': 2,
        # ... all other INCAR tags
    },
    'relax': {
        'positions': True,
        'shape': False,
        'volume': False
    },
    # ... any other builder attributes
})
```

**Rationale:** Direct pass-through gives users full control over VASP settings without abstraction layers. Same config applies to all generated structures (single shared configuration requirement).

### surface_params Structure

```python
params_dict = Dict({
    'mode': 'combine',              # 'vacancies', 'hydrogen', 'combine'
    'species': 'O',                 # Target species for modification
    'z_window': 0.5,                # Å window to identify surface atoms
    'which_surface': 'top',         # 'top', 'bottom', 'both'
    'oh_dist': 0.98,                # O–H distance in Å
    'include_empty': False,         # Include unmodified baseline
    'supercell': None,              # [nx, ny, nz] or None
    'deduplicate_by_coverage': True,
    'coverage_bins': 10             # Number of coverage bins for sampling
})
```

## Output Specifications

### successful_relaxations

List of Dicts, one per successful relaxation:
```python
{
    'name': 'combo_001_1.5_0.75',           # From manifest
    'structure': StructureData,             # Relaxed structure
    'energy': Float,                        # Total energy (eV)
    'coverage': Tuple[float, float],        # (vac_coverage, oh_coverage)
    'metadata': Dict                        # Full variant metadata from manifest
}
```

### failed_relaxations

List of Dicts, one per failed relaxation:
```python
{
    'name': 'combo_015_2.0_1.5',
    'coverage': Tuple[float, float],
    'exit_status': Int,                     # VASP exit status
    'error_message': str                    # Error description
}
```

### statistics

```python
{
    'total': 15,
    'succeeded': 13,
    'failed': 2
}
```

## Parallelization Strategy

**Implementation:** WorkGraph awaitable/semaphore mechanism

**How it works:**
1. User sets `max_parallel_jobs = 5` (example)
2. `RelaxationsWorkGraph` creates 15 VaspRelaxWorkChain tasks (example)
3. Semaphore limits concurrent execution to 5
4. First 5 tasks start immediately
5. When one completes, next waiting task starts
6. Process continues until all tasks complete

**Advantages:**
- Native WorkGraph pattern
- Optimal resource utilization (no waiting for slowest job in batch)
- Clean implementation
- Easy to monitor and debug

## Error Handling

**Strategy:** Continue with partial results

- If relaxation succeeds (exit_status == 0): Extract structure + energy
- If relaxation fails (exit_status != 0): Record in `failed_relaxations`
- Workflow always completes successfully
- User decides how to handle partial datasets

**Rationale:** Better to get 13/15 structures than lose everything. Failures are explicitly tracked and reported.

## Integration Details

### surface_modes.py Integration

- Move from `teros/experimental/vacancies_hydroxilation/` to `teros/core/surface_hydroxylation/`
- Import as Python module (not subprocess call)
- Use `SurfaceModifier` class directly in `generate_structures` calcfunction
- Leverage existing functionality: `run_vacancies()`, `run_hydrogen()`, `run_combine()`

### Provenance Tracking

All steps tracked in AiiDA graph:
```
Input slab (StructureData)
    ↓
generate_structures (CalcFunction)
    ↓
List of variant structures (List[StructureData])
    ↓
VaspRelaxWorkChain × N (parallel, semaphore-limited)
    ↓
Relaxed structures + energies
    ↓
collect_results (CalcFunction)
    ↓
Final outputs (namespace)
```

Full lineage queryable through AiiDA database.

## Testing Strategy

**Test location:** `/home/thiagotd/git/PS-TEROS/examples/surface_hydroxylation/`

**Test plan:**
1. Use small test slab (e.g., 2×2 perovskite (001) surface, ~40 atoms)
2. Configure surface_params to generate ~10 structures:
   - `mode = 'combine'`
   - `deduplicate_by_coverage = True`
   - `coverage_bins = 5`
3. Set `max_parallel_jobs = 3` to test semaphore
4. Configure lightweight VASP settings (fast relaxation for testing)
5. Launch workflow
6. Wait for completion (sleep 60 or check with `verdi process show`)
7. **Success criteria:**
   - Main workflow node returns `exit_status = [0]`
   - `outputs.statistics` shows expected total count
   - `outputs.successful_relaxations` contains structures and energies
   - All relaxation tasks visible in AiiDA graph

**Production validation:**
After testing, run with production parameters:
- Realistic surface slab (~100 atoms)
- Full coverage range (15-20 structures)
- Production VASP settings
- Verify results are physically reasonable

## Future Extensions (Out of Current Scope)

Once core module is stable, add:

1. **Thermodynamic analysis:**
   - Formation energies vs coverage
   - Phase diagrams (coverage vs chemical potential)
   - Grand canonical ensemble analysis

2. **Pourbaix/electrochemical analysis:**
   - Potential-dependent stability
   - pH-potential phase diagrams
   - Surface speciation vs applied potential

3. **Additional surface modifications:**
   - Substitutional doping
   - Adsorbed molecules (H₂O, CO₂, etc.)
   - Multi-species hydroxylation (F-OH, Cl-OH)

4. **Advanced sampling:**
   - Cluster expansion for coverage interpolation
   - Machine learning surrogate models
   - Active learning for structure selection

## Dependencies

- AiiDA-core
- AiiDA-WorkGraph
- aiida-vasp (vasp.v2.relax plugin)
- ASE (for structure manipulation in surface_modes.py)
- NumPy (for surface_modes.py calculations)

## Implementation Checklist

- [ ] Create module directory structure
- [ ] Move and refactor surface_modes.py
- [ ] Implement generate_structures calcfunction
- [ ] Implement RelaxationsWorkGraph with semaphore
- [ ] Implement collect_results calcfunction
- [ ] Implement main SurfaceHydroxylationWorkGraph
- [ ] Create test example in examples/surface_hydroxylation/
- [ ] Write usage documentation
- [ ] Test with small system
- [ ] Validate with production parameters
- [ ] Create pull request to develop branch

---

**Design approved:** Ready for implementation
