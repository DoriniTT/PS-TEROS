# AIMD Standalone Module Design

**Date:** 2025-11-03
**Status:** Design Complete
**Location:** `teros/core/aimd/`

## Purpose

Create a standalone AIMD module that gives users complete control over multi-stage molecular dynamics calculations on multiple structures. The module accepts pre-existing structures (from hydroxylation, adsorption, or any source), runs sequential AIMD stages with full parameter control, and returns organized results.

## Requirements

- Accept StructureData nodes or PKs directly (no bulk workflow coupling)
- Run sequential AIMD stages on all structures (stage N restarts from stage N-1)
- Apply supercell transformation optionally per structure before AIMD
- Control builder parameters at three levels: global, per-structure, per-stage
- Override any parameter for specific (structure, stage) combinations (matrix control)
- Limit concurrent VASP calculations with `max_concurrent_jobs`
- Return organized results with energies, trajectories, and structures per stage

## Module Structure

```
teros/core/aimd/
├── __init__.py           # Exports main functions
├── workgraph.py          # Main entry: build_aimd_workgraph()
├── tasks.py              # WorkGraph tasks for supercell and AIMD orchestration
├── utils.py              # Result organization, validation, merging
└── README.md             # Module documentation
```

## Architecture

### Reuse Existing Code

The module imports and orchestrates existing AIMD functions:

- `prepare_aimd_parameters()` from `teros.core.aimd` - inject temperature/steps into INCAR
- `aimd_single_stage_scatter()` from `teros.core.aimd` - parallel AIMD on multiple structures
- `extract_total_energy()` from `teros.core.slabs` - extract energies from outputs

### WorkGraph Flow

1. **Input preparation** (workgraph.py)
   - Load structures (resolve PKs if needed)
   - Validate `aimd_stages` format: `[{'temperature': K, 'steps': N}, ...]`
   - Apply supercells if `supercell_specs` provided

2. **Stage loop** (sequential)
   - For each stage:
     - Merge builder inputs: global → structure override → stage override → matrix override
     - Call `aimd_single_stage_scatter()` with merged parameters
     - Chain remote_folders: stage N output → stage N+1 restart input
     - Pass `max_concurrent_jobs` to control parallelism within stage

3. **Output collection** (utils.py)
   - Organize results by structure and stage
   - Return nested namespace with structures, energies, trajectories, remote_folders

### Concurrency Model

- **Within stage**: all structures run in parallel (limited by `max_concurrent_jobs`)
- **Across stages**: strictly sequential (stage 2 waits for stage 1)

## API Design

### Main Function

```python
def build_aimd_workgraph(
    structures: dict[str, orm.StructureData | int],
    aimd_stages: list[dict],
    code_label: str,
    builder_inputs: dict,

    supercell_specs: dict[str, list[int]] = None,
    structure_overrides: dict[str, dict] = None,
    stage_overrides: dict[int, dict] = None,
    matrix_overrides: dict[tuple, dict] = None,

    max_concurrent_jobs: int = None,
    name: str = 'AIMDWorkGraph',
) -> WorkGraph
```

### Parameters

- **structures**: `{name: StructureData or PK}` - input structures for AIMD
- **aimd_stages**: `[{'temperature': 300, 'steps': 100}, ...]` - sequential stages
- **code_label**: `'VASP6.5.0@cluster02'` - VASP code label
- **builder_inputs**: default builder config (parameters, kpoints_spacing, options, etc.)
- **supercell_specs**: `{structure_name: [nx, ny, nz]}` - optional supercell per structure
- **structure_overrides**: `{structure_name: builder_inputs}` - override per structure
- **stage_overrides**: `{stage_idx: builder_inputs}` - override per stage (0-indexed)
- **matrix_overrides**: `{(structure_name, stage_idx): builder_inputs}` - specific combinations
- **max_concurrent_jobs**: limit parallel VASP calculations (None = unlimited)
- **name**: WorkGraph name

### Builder Override Priority

```
matrix_overrides > stage_overrides > structure_overrides > builder_inputs
```

Higher priority wins. Merging uses deep merge (nested dicts combined, not replaced).

### Example Usage

```python
# Simple: same settings everywhere
wg = build_aimd_workgraph(
    structures={'slab1': structure1, 'slab2': structure2},
    aimd_stages=[
        {'temperature': 300, 'steps': 100},   # Equilibration
        {'temperature': 300, 'steps': 500},   # Production
    ],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={
        'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400, ...}},
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'potential_mapping': {},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    },
    max_concurrent_jobs=4,
)

# Advanced: full control
wg = build_aimd_workgraph(
    structures={'slab1': pk1, 'slab2': pk2},
    aimd_stages=[...],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={...},

    # Supercell for slab1 only
    supercell_specs={'slab1': [2, 2, 1]},

    # slab2 needs more resources everywhere
    structure_overrides={
        'slab2': {'options': {'resources': {'num_cores_per_machine': 48}}}
    },

    # Stage 1 (production) needs tighter convergence
    stage_overrides={
        1: {'parameters': {'incar': {'EDIFF': 1e-7}}}
    },

    # slab1 + stage 1 needs special algorithm
    matrix_overrides={
        ('slab1', 1): {'parameters': {'incar': {'ALGO': 'All'}}}
    },

    max_concurrent_jobs=2,
)
```

## Outputs

### WorkGraph Namespace

```python
{
    'results': {
        'slab1': {
            0: {
                'structure': StructureData,
                'energy': Float,
                'remote_folder': RemoteData,
                'trajectory': TrajectoryData,
            },
            1: {...},
        },
        'slab2': {
            0: {...},
            1: {...},
        },
    },
    'supercells': {  # Only if supercell_specs provided
        'slab1': StructureData,
    }
}
```

### Result Organization Helper

```python
def organize_aimd_results(node: orm.WorkGraphNode) -> dict:
    """
    Organize AIMD results into user-friendly format.

    Returns:
        {
            'summary': {
                'total_structures': int,
                'total_stages': int,
                'successful': int,
                'failed': int,
            },
            'results': {
                structure_name: [
                    {
                        'stage': 0,
                        'temperature': 300,
                        'steps': 100,
                        'energy': float,
                        'structure_pk': int,
                        'trajectory_pk': int,
                    },
                    {...},
                ]
            },
            'failed_calculations': [
                {'structure': 'slab1', 'stage': 0, 'error': '...'},
            ],
        }
    """
```

## Validation & Error Handling

### Input Validation

- **aimd_stages**: each dict must contain `temperature` and `steps` keys
- **supercell_specs**: values must be 3-element lists `[nx, ny, nz]` with positive integers
- **override keys**: structure names must exist in `structures`, stage indices in range
- **structures**: PKs must resolve to StructureData nodes

### Error Handling

- **Invalid inputs**: raise `ValueError` with clear message
- **VASP calculation fails**: continue other structures, report in `failed_calculations`
- **Supercell creation fails**: raise error immediately (cannot proceed)
- **Missing override keys**: raise `KeyError` with helpful message

### Validation Functions

```python
def validate_stage_sequence(stages: list[dict]) -> None:
    """Raise ValueError if stages missing required keys."""

def validate_supercell_spec(spec: list[int]) -> None:
    """Raise ValueError if spec not [nx, ny, nz] with positive integers."""

def merge_builder_inputs(base: dict, override: dict) -> dict:
    """Deep merge override into base, return new dict."""

def create_supercell_with_pymatgen(
    structure: orm.StructureData,
    spec: list[int]
) -> orm.StructureData:
    """Create supercell using pymatgen, return new StructureData."""
```

## Implementation Details

### Supercell Creation (tasks.py)

```python
@task
def create_supercell(structure: orm.StructureData, spec: list[int]) -> orm.StructureData:
    """
    Create supercell using pymatgen.

    Args:
        structure: Input structure
        spec: [nx, ny, nz] supercell dimensions

    Returns:
        Supercell StructureData
    """
    from pymatgen.core import Structure

    # Convert to pymatgen
    pmg_struct = structure.get_pymatgen()

    # Create supercell
    pmg_supercell = pmg_struct * spec

    # Convert back to AiiDA
    supercell = orm.StructureData(pymatgen=pmg_supercell)
    return supercell
```

### Builder Merging (utils.py)

```python
def merge_builder_inputs(base: dict, override: dict) -> dict:
    """
    Deep merge override into base.

    Nested dicts are recursively merged.
    Non-dict values in override replace base values.

    Example:
        base = {'parameters': {'incar': {'ENCUT': 400, 'PREC': 'Normal'}}}
        override = {'parameters': {'incar': {'ENCUT': 500}}}
        result = {'parameters': {'incar': {'ENCUT': 500, 'PREC': 'Normal'}}}
    """
    from copy import deepcopy

    result = deepcopy(base)

    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge_builder_inputs(result[key], value)
        else:
            result[key] = deepcopy(value)

    return result
```

### Stage Orchestration (workgraph.py)

```python
def build_aimd_workgraph(...):
    """Build AIMD workgraph with sequential stages."""

    wg = WorkGraph(name=name)

    # Set max_concurrent_jobs on workgraph
    if max_concurrent_jobs:
        wg.max_number_jobs = max_concurrent_jobs

    # 1. Load and prepare structures
    prepared_structures = {}
    supercell_outputs = {}

    for struct_name, struct_input in structures.items():
        # Load structure if PK
        if isinstance(struct_input, int):
            struct = orm.load_node(struct_input)
        else:
            struct = struct_input

        # Create supercell if requested
        if supercell_specs and struct_name in supercell_specs:
            sc_task = wg.add_task(
                create_supercell,
                structure=struct,
                spec=supercell_specs[struct_name],
                name=f'create_supercell_{struct_name}',
            )
            prepared_structures[struct_name] = sc_task.outputs.result
            supercell_outputs[struct_name] = sc_task.outputs.result
        else:
            prepared_structures[struct_name] = struct

    # 2. Run sequential AIMD stages
    stage_outputs = {}
    previous_remote_folders = {}

    for stage_idx, stage_config in enumerate(aimd_stages):
        # Build builder inputs for each structure in this stage
        stage_builders = {}
        for struct_name in prepared_structures:
            builder = merge_all_overrides(
                struct_name, stage_idx,
                builder_inputs, structure_overrides,
                stage_overrides, matrix_overrides,
            )
            stage_builders[struct_name] = builder

        # Create AIMD task for this stage
        aimd_task = wg.add_task(
            aimd_single_stage_scatter,
            slabs=prepared_structures,
            temperature=stage_config['temperature'],
            steps=stage_config['steps'],
            code=orm.load_code(code_label),
            aimd_parameters=stage_builders,  # Per-structure builders
            restart_folders=previous_remote_folders if stage_idx > 0 else {},
            max_number_jobs=max_concurrent_jobs,
            name=f'stage_{stage_idx}_aimd',
        )

        # Store outputs
        stage_outputs[stage_idx] = aimd_task.outputs
        prepared_structures = aimd_task.outputs.structures
        previous_remote_folders = aimd_task.outputs.remote_folders

    # 3. Organize outputs
    wg.add_output('results', stage_outputs)
    if supercell_outputs:
        wg.add_output('supercells', supercell_outputs)

    return wg
```

## Testing Strategy

Create example script: `examples/vasp/step_XX_aimd_standalone.py`

Test cases:
1. Single structure, single stage (minimal)
2. Multiple structures, single stage (parallelism)
3. Single structure, multiple stages (restart chaining)
4. Multiple structures, multiple stages (full workflow)
5. With supercell transformation
6. With structure overrides
7. With stage overrides
8. With matrix overrides (full control)
9. With `max_concurrent_jobs` limiting

Each test validates:
- WorkGraph completes with exit_status=0
- Output namespace contains expected structures/energies
- Stage N restarts from stage N-1 (check CONTCAR usage)
- Overrides applied correctly (check INCAR values)
- Supercells created with correct dimensions

## Migration Path

Current `workgraph.py` AIMD integration remains unchanged. The new module provides an alternative API for users who want:
- Direct structure input (not bulk+slab workflow)
- Full builder control
- Simpler API focused on AIMD only

Users can choose:
- **Option 1**: `build_core_workgraph(workflow_preset='aimd_only', ...)` - full workflow
- **Option 2**: `build_aimd_workgraph(structures={...}, ...)` - standalone AIMD

Both use the same underlying `aimd_single_stage_scatter()` function.

## Documentation

Create `teros/core/aimd/README.md` explaining:
- Module purpose and design
- API reference with examples
- Relationship to `workgraph.py` AIMD integration
- Parameter override system with priority rules
- Common use cases and patterns

## Success Criteria

Module implementation succeeds when:
1. Example script completes with all exit_status=0
2. Results show correct stage progression (energies, structures)
3. Overrides apply correctly at all levels
4. Supercells have correct dimensions
5. `max_concurrent_jobs` limits parallelism
6. Code mirrors hydroxylation module structure
7. Documentation clear and complete
