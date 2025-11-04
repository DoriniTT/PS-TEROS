# AIMD Standalone Module

## Overview

The AIMD standalone module (`teros.core.aimd`) provides a simplified, flexible interface for running Ab Initio Molecular Dynamics (AIMD) simulations on pre-existing structures, independent of the main bulk+slab workflow. This module is designed for users who need fine-grained control over AIMD parameters and want to run MD simulations on structures from any source.

**Key capabilities:**
- Direct structure input (StructureData nodes or PKs)
- Sequential multi-stage AIMD with automatic restart chaining
- Per-structure parameter customization via three-level override system
- Optional supercell transformations before AIMD
- Concurrency control for parallel calculations
- Complete AiiDA provenance tracking

## Scientific Background

### Ab Initio Molecular Dynamics

AIMD combines molecular dynamics simulation with electronic structure calculations, computing forces on atoms from first-principles quantum mechanics (DFT) at each timestep. This enables:

- **Temperature effects**: Study materials at finite temperature
- **Structural dynamics**: Observe atomic motion, phase transitions, diffusion
- **Free energy sampling**: Calculate thermodynamic properties via MD trajectories
- **Reaction pathways**: Discover reaction mechanisms through thermal fluctuations

### Born-Oppenheimer Molecular Dynamics

PS-TEROS uses Born-Oppenheimer MD (VASP's IBRION=0), where:
1. Electronic structure is converged at each ionic step
2. Forces on atoms are computed from converged electronic state
3. Newton's equations of motion are integrated to update positions
4. Temperature is controlled via thermostat (Nosé-Hoover, MDALGO=2)

**Key VASP parameters for AIMD:**
- `IBRION = 0`: Switch to MD mode
- `MDALGO = 2`: Nosé-Hoover thermostat (canonical ensemble, NVT)
- `POTIM`: MD timestep in femtoseconds (typically 1-3 fs)
- `TEBEG/TEEND`: Temperature in Kelvin
- `NSW`: Number of MD steps
- `SMASS`: Nosé mass (0.0 = automatic)

## Module Architecture

### Design Philosophy

The standalone AIMD module follows these principles:

1. **Simplicity**: Minimal required parameters, sensible defaults
2. **Flexibility**: Three-level override system for per-structure customization
3. **Composability**: Works with any AiiDA StructureData source
4. **Transparency**: All parameters explicit, no hidden defaults
5. **Scalability**: Handles 1 to 100+ structures efficiently

### Workflow Structure

```
Input Structures
    ↓
Optional Supercell Creation (parallel)
    ↓
Stage 0: Equilibration (all structures in parallel)
    ↓
Stage 1: Production (restarts from Stage 0)
    ↓
Stage N: Extended sampling (restarts from Stage N-1)
    ↓
Output: Trajectories, final structures, energies
```

**Key behaviors:**
- Within each stage: structures run **in parallel** (limited by `max_concurrent_jobs`)
- Across stages: strictly **sequential** (Stage N+1 waits for Stage N)
- Restart chaining: automatic `remote_folder` linking between stages

## Quick Start

### Basic Example: Single Structure, Two Stages

```python
from teros.core.aimd import build_aimd_workgraph
from aiida import orm, load_profile
from ase.io import read

load_profile('presto')

# Load structure
structure = orm.StructureData(ase=read('my_slab.cif'))

# Define AIMD stages using VASP-native parameter names
aimd_stages = [
    {'TEBEG': 300, 'NSW': 100},   # Equilibration
    {'TEBEG': 300, 'NSW': 500},   # Production
]

# Builder inputs (VASP parameters)
builder_inputs = {
    'parameters': {
        'incar': {
            'PREC': 'Normal',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 0,      # MD mode
            'MDALGO': 2,      # Nosé-Hoover
            'POTIM': 2.0,     # 2 fs timestep
        }
    },
    'kpoints_spacing': 0.5,
    'potential_family': 'PBE',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    },
    'clean_workdir': False,
}

# Build and submit
wg = build_aimd_workgraph(
    structures={'my_slab': structure},
    aimd_stages=aimd_stages,
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=builder_inputs,
    name='MyAIMD',
)

wg.submit(wait=False)
print(f"WorkGraph PK: {wg.pk}")
```

## Override System

The override system enables per-structure INCAR parameter customization across three levels:

### Override Levels

1. **Structure-level**: Apply to all stages of specific structures
2. **Stage-level**: Apply to all structures in specific stages
3. **Matrix-level**: Apply to specific (structure, stage) combinations

**Priority order**: matrix > stage > structure > base

### Supported Parameters

The override system supports **INCAR parameters only**:
- `ENCUT`, `PREC`, `EDIFF`, `EDIFFG`
- `ALGO`, `ISIF`, `IBRION`
- `NCORE`, `KPAR`, `LREAL`
- Any other VASP INCAR tag

### Not Supported

These builder inputs remain **uniform** across all structures:
- `kpoints_spacing` - K-points grid density
- `options` - Scheduler settings (cores, walltime, etc.)
- `potential_mapping` - Element to pseudopotential mapping
- `potential_family` - Pseudopotential family

**Workaround**: Use separate `build_aimd_workgraph()` calls for structures needing different kpoints/options.

### Example: All Three Override Levels

```python
wg = build_aimd_workgraph(
    structures={'slab1': s1, 'slab2': s2},
    aimd_stages=[
        {'TEBEG': 300, 'NSW': 100},   # Stage 0
        {'TEBEG': 300, 'NSW': 500},   # Stage 1
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs={
        'parameters': {'incar': {'ENCUT': 400, 'PREC': 'Normal'}},
        # ... other parameters
    },

    # Structure-level: slab2 needs higher cutoff everywhere
    structure_overrides={
        'slab2': {'parameters': {'incar': {'ENCUT': 500}}}
    },

    # Stage-level: Stage 1 (production) needs tighter convergence
    stage_overrides={
        1: {'parameters': {'incar': {'EDIFF': 1e-7}}}
    },

    # Matrix-level: slab1 + stage 1 needs special algorithm
    matrix_overrides={
        ('slab1', 1): {'parameters': {'incar': {'ALGO': 'All'}}}
    },
)
```

**Result:**
- `slab1, stage 0`: ENCUT=400, PREC=Normal (base)
- `slab1, stage 1`: ENCUT=400, PREC=Normal, EDIFF=1e-7, ALGO=All (stage + matrix)
- `slab2, stage 0`: ENCUT=500, PREC=Normal (structure override)
- `slab2, stage 1`: ENCUT=500, PREC=Normal, EDIFF=1e-7 (structure + stage)

### Override Format

All override parameters must follow this nested structure:

```python
{
    'parameters': {
        'incar': {
            'INCAR_KEY': value,
            # ... more INCAR parameters
        }
    }
}
```

**Important**: Only the `'parameters'/'incar'` path is processed. Other keys (like `'kpoints_spacing'`, `'options'`) are ignored.

## Common Use Cases

### Use Case 1: Temperature Series

Run the same structure at multiple temperatures:

```python
# Option 1: Separate workflows (if different ENCUT needed)
for temp in [300, 400, 500, 600]:
    wg = build_aimd_workgraph(
        structures={'slab': structure},
        aimd_stages=[{'temperature': temp, 'steps': 200}],
        code_label='VASP-6.5.1@cluster02',
        builder_inputs=base_inputs,
        name=f'AIMD_{temp}K',
    )
    wg.submit(wait=False)

# Option 2: Stage-based (if all parameters identical)
wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[
        {'temperature': 300, 'steps': 200},
        {'temperature': 400, 'steps': 200},
        {'temperature': 500, 'steps': 200},
        {'temperature': 600, 'steps': 200},
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=base_inputs,
)
```

### Use Case 2: Multiple Structures with Different Convergence

```python
wg = build_aimd_workgraph(
    structures={
        'AgO_100': structure1,
        'AgO_110': structure2,
        'AgO_111': structure3,
    },
    aimd_stages=[
        {'temperature': 300, 'steps': 50},    # Equilibration
        {'temperature': 300, 'steps': 200},   # Production
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs={
        'parameters': {'incar': {'ENCUT': 400, 'EDIFF': 1e-5}},
        # ... other parameters
    },

    # 110 surface needs tighter convergence
    structure_overrides={
        'AgO_110': {'parameters': {'incar': {'ENCUT': 500, 'EDIFF': 1e-6}}}
    },

    max_concurrent_jobs=2,  # Only 2 VASP jobs at once
)
```

### Use Case 3: Equilibration + Production with Different Precision

```python
wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[
        {'temperature': 300, 'steps': 100},   # Equilibration
        {'temperature': 300, 'steps': 500},   # Production
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs={
        'parameters': {'incar': {'ENCUT': 400, 'PREC': 'Normal'}},
        # ... other parameters
    },

    # Production stage needs higher accuracy
    stage_overrides={
        1: {'parameters': {'incar': {'PREC': 'Accurate', 'EDIFF': 1e-7}}}
    },
)
```

### Use Case 4: Supercell Creation

```python
wg = build_aimd_workgraph(
    structures={'small_slab': structure},
    aimd_stages=[{'temperature': 300, 'steps': 200}],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=base_inputs,

    # Create 3x3x1 supercell before AIMD
    supercell_specs={'small_slab': [3, 3, 1]},
)
```

## API Reference

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

#### Parameters

**structures** : `dict[str, orm.StructureData | int]`
- Input structures for AIMD
- Format: `{name: StructureData_node}` or `{name: PK}`
- Example: `{'slab1': structure, 'slab2': 12345}`

**aimd_stages** : `list[dict]`
- Sequential AIMD stages using VASP-native parameter names
- Required: `TEBEG` (initial temperature K), `NSW` (MD steps)
- Optional: `TEEND` (final temperature, defaults to TEBEG), `POTIM` (timestep fs), `MDALGO` (thermostat), `SMASS` (Nosé mass)
- Example: `[{'TEBEG': 300, 'NSW': 100}, {'TEBEG': 400, 'NSW': 200}]`
- Temperature annealing: `[{'TEBEG': 300, 'TEEND': 500, 'NSW': 200}]`

**code_label** : `str`
- VASP code label from AiiDA
- Format: `'CodeName@ComputerName'`
- Example: `'VASP-6.5.1@cluster02'`

**builder_inputs** : `dict`
- Default VASP builder configuration
- Required keys:
  - `parameters`: `{'incar': {...}}` - VASP INCAR settings
  - `kpoints_spacing`: `float` - K-points spacing
  - `potential_family`: `str` - Pseudopotential family
  - `potential_mapping`: `dict` - Element to potential mapping
  - `options`: `dict` - Scheduler options
  - `clean_workdir`: `bool` - Whether to clean work directory

**supercell_specs** : `dict[str, list[int]]`, optional
- Supercell transformations to apply before AIMD
- Format: `{structure_name: [nx, ny, nz]}`
- Example: `{'slab1': [2, 2, 1]}` creates 2×2×1 supercell

**structure_overrides** : `dict[str, dict]`, optional
- Per-structure builder overrides
- Only `'parameters'/'incar'` keys are applied
- Format: `{structure_name: {'parameters': {'incar': {...}}}}`

**stage_overrides** : `dict[int, dict]`, optional
- Per-stage builder overrides (0-indexed)
- Only `'parameters'/'incar'` keys are applied
- Format: `{stage_idx: {'parameters': {'incar': {...}}}}`

**matrix_overrides** : `dict[tuple, dict]`, optional
- Per-(structure, stage) builder overrides
- Only `'parameters'/'incar'` keys are applied
- Format: `{(structure_name, stage_idx): {'parameters': {'incar': {...}}}}`

**max_concurrent_jobs** : `int`, optional
- Maximum number of concurrent VASP calculations
- Default: `None` (unlimited within WorkGraph limits)

**name** : `str`, optional
- WorkGraph name for identification
- Default: `'AIMDWorkGraph'`

#### Returns

**WorkGraph**
- AiiDA WorkGraph ready to submit
- Submit with: `wg.submit(wait=False)`

### Output Structure

Results are accessible through the WorkGraph node after completion:

```python
wg_node = orm.load_node(wg.pk)

# Access individual stage outputs
stage_0_structures = wg_node.outputs.stage_0_structures
stage_0_energies = wg_node.outputs.stage_0_energies
stage_1_structures = wg_node.outputs.stage_1_structures

# Access supercells (if created)
supercell_struct1 = wg_node.outputs.supercell_struct1
```

## Monitoring and Verification

### Check Workflow Status

```bash
# Overall status
verdi process show <WorkGraph_PK>

# Specific stage
verdi process show <stage_PK>

# Logs and error messages
verdi process report <PK>

# Real-time monitoring
verdi process watch <PK>
```

### Verify INCAR Overrides Applied

```bash
# 1. Get WorkGraph PK from submission output
verdi process show <WorkGraph_PK>

# 2. Find stage_0_aimd task PK
verdi process show <WorkGraph_PK> | grep stage_0_aimd

# 3. Get VaspWorkChain PKs from stage task
verdi process show <stage_0_PK>

# 4. Get VaspCalculation PK from VaspWorkChain
verdi process show <VaspWorkChain_PK>

# 5. Check INCAR for calculation
verdi calcjob inputcat <VaspCalculation_PK> INCAR | grep -E 'ENCUT|PREC|ALGO'
```

**Example verification for override system:**

```bash
# Expected for ag1, stage 0: ENCUT=400, PREC=Normal, ALGO=Normal
verdi calcjob inputcat 124520 INCAR | grep -E 'ENCUT|PREC|ALGO'
# Output:
# ENCUT = 400
# PREC = Normal
# ALGO = Normal

# Expected for ag2, stage 0: ENCUT=500 (structure override)
verdi calcjob inputcat 124523 INCAR | grep -E 'ENCUT|PREC|ALGO'
# Output:
# ENCUT = 500
# PREC = Normal
# ALGO = Normal
```

## Best Practices

### AIMD Parameter Selection

**Timestep (POTIM)**:
- Light elements (H, C, O): 0.5-1.5 fs
- Medium elements (Si, Al, Ag): 1.5-2.5 fs
- Heavy elements (Pb, Au, W): 2.5-4.0 fs
- Rule of thumb: ~1/20 of fastest vibrational period

**Number of steps (NSW)**:
- Equilibration: 50-200 steps (0.1-0.5 ps)
- Production: 500-5000 steps (1-10 ps)
- Diffusion studies: 10,000+ steps (20+ ps)

**Convergence (EDIFF)**:
- Equilibration: 1e-5 eV (faster, less accurate)
- Production: 1e-6 to 1e-7 eV (slower, more accurate)

**k-points**:
- Use Γ-point only for large supercells (>200 atoms)
- For smaller cells: kpoints_spacing=0.4-0.5 Å⁻¹

### Multi-Stage Strategy

**Recommended stage sequence:**

1. **Coarse equilibration** (50-100 steps, PREC=Normal, EDIFF=1e-5)
   - Fast thermalization
   - Removes initial strain

2. **Fine equilibration** (100-200 steps, PREC=Accurate, EDIFF=1e-6)
   - Better energy conservation
   - Stabilize temperature

3. **Production** (500-5000 steps, PREC=Accurate, EDIFF=1e-6)
   - Collect statistics
   - Calculate properties

### Resource Management

**For multiple structures:**

```python
# Limit parallel jobs to avoid cluster overload
wg = build_aimd_workgraph(
    structures={f'struct_{i}': struct for i in range(20)},
    aimd_stages=[{'temperature': 300, 'steps': 200}],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=base_inputs,
    max_concurrent_jobs=5,  # Only 5 VASP jobs at once
)
```

**For long trajectories:**

```python
# Break into stages to checkpoint progress
wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[
        {'temperature': 300, 'steps': 500},   # Checkpoint every 500 steps
        {'temperature': 300, 'steps': 500},
        {'temperature': 300, 'steps': 500},
        {'temperature': 300, 'steps': 500},
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=base_inputs,
)
```

## Troubleshooting

### Issue: AIMD not equilibrating

**Symptoms**: Temperature fluctuates wildly, energy not conserved

**Solutions**:
- Decrease `POTIM` (timestep too large)
- Increase `NSW` (needs more steps to equilibrate)
- Check `SMASS` (try 0.0 for automatic or -1 for microcanonical)
- Use `stage_overrides` to increase equilibration steps

### Issue: Different structures need different parameters

**Solution**: Use override system

```python
structure_overrides={
    'high_symmetry_slab': {
        'parameters': {'incar': {'ENCUT': 500, 'EDIFF': 1e-7}}
    },
    'defect_slab': {
        'parameters': {'incar': {'ALGO': 'All', 'NELM': 200}}
    },
}
```

### Issue: Out of memory errors

**Solutions**:
- Increase `options['resources']['num_cores_per_machine']`
- Use `NCORE` or `KPAR` in INCAR to improve parallelization
- Reduce `ENCUT` for equilibration stages
- Use `LREAL = Auto` for large systems

### Issue: Calculations not starting

**Check**:
1. Code exists: `verdi code show VASP-6.5.1@cluster02`
2. Daemon running: `verdi daemon status`
3. Structures valid: `verdi data core.structure show <PK>`
4. No errors: `verdi process report <PK>`

## Testing

Unit tests are located in `teros/core/aimd/test_*.py`:

```bash
# Run all AIMD tests
pytest teros/core/aimd/test_*.py -v

# Run override system tests only
pytest teros/core/aimd/test_overrides.py -v

# Run specific test
pytest teros/core/aimd/test_overrides.py::test_structure_overrides -v
```

Full workflow test:

```bash
python examples/vasp/step_19_aimd_with_overrides.py
```

## Module Structure

```
teros/core/aimd/
├── __init__.py           # Exports: build_aimd_workgraph, organize_aimd_results
├── workgraph.py          # Main entry: build_aimd_workgraph()
├── tasks.py              # WorkGraph tasks: create_supercell()
├── utils.py              # Validation and merging utilities
├── test_utils.py         # Unit tests for utils
├── test_tasks.py         # Unit tests for tasks
├── test_overrides.py     # Unit tests for override system
└── README.md             # Module-level documentation
```

## Relationship to Main Workflow

This module is **independent** from `teros.core.workgraph.build_core_workgraph()` but reuses the underlying `aimd_single_stage_scatter()` function.

**Choose standalone AIMD when:**
- You have pre-existing structures (from any source)
- You want simple AIMD-only workflows
- You need full control over all parameters
- You need per-structure parameter customization

**Choose main workflow when:**
- You want full bulk → slab → relaxation → AIMD pipeline
- You want workflow presets (`'aimd_only'`, `'surface_energy'`, etc.)
- You want automatic cleavage energy calculations

## Examples

### Example 1: Basic AIMD

See: `examples/vasp/step_18_aimd_standalone.py`

Single structure, two-stage AIMD with supercell creation.

### Example 2: Override System

See: `examples/vasp/step_19_aimd_with_overrides.py`

Demonstrates all three override levels with two structures and verification instructions.

## Implementation Notes

- Reuses `aimd_single_stage_scatter()` from `teros.core.aimd_functions`
- Supercells created with pymatgen via ASE adapter
- WorkGraph handles task orchestration and restart chaining
- All AiiDA nodes stored in provenance graph
- Override system uses shallow dictionary merging with priority-based application

## See Also

- Module README: `teros/core/aimd/README.md`
- Design document: `docs/plans/2025-11-03-aimd-override-system.md`
- Implementation plan: `docs/plans/2025-11-03-aimd-override-implementation.md`
- Example scripts: `examples/vasp/step_18_aimd_standalone.py`, `step_19_aimd_with_overrides.py`
- Main workflow guide: `docs/WORKFLOW_SYSTEM_EXPLAINED.md`
