# TODO: Explorer as Central Orchestration Hub

## Vision

Transform the `explorer` module into the **central orchestration hub** for all PS-TEROS workflows. Instead of just chaining VASP calculations, users could compose complex multi-physics workflows by combining different modules like **LEGO blocks**.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         EXPLORER: Central Hub                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Stage 1        Stage 2        Stage 3        Stage 4        Stage 5      │
│  ┌────────┐    ┌────────┐    ┌────────┐    ┌────────┐    ┌────────┐       │
│  │  VASP  │───►│  VASP  │───►│  VASP  │───►│  DOS   │───►│  AIMD  │       │
│  │ relax  │    │ relax  │    │supercell│    │ module │    │ module │       │
│  └────────┘    └────────┘    └────────┘    └────────┘    └────────┘       │
│      │              │              │              │              │          │
│      ▼              ▼              ▼              ▼              ▼          │
│   outputs        outputs        outputs        outputs        outputs      │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Motivation

Currently, to run a complex workflow like:
1. Relax a slab (rough → fine)
2. Create supercell and relax
3. Calculate DOS
4. Run AIMD simulation

Users must manually orchestrate separate workflows and pass data between them. With the "LEGO" approach, this becomes a single, unified workflow with automatic data passing.

## Example: The Dream API

```python
from teros.core.explorer import quick_workflow

stages = [
    # Stage 1-4: VASP relaxations (current functionality)
    {
        'name': 'relax_1x1_rough',
        'type': 'vasp',  # NEW: stage type
        'incar': incar_rough,
        'restart': None,
        'kpoints_spacing': 0.06,
        'fix_type': 'center',
        'fix_thickness': 7.0,
    },
    {
        'name': 'relax_1x1_fine',
        'type': 'vasp',
        'incar': incar_fine,
        'restart': None,
        'kpoints_spacing': 0.03,
        'fix_type': 'center',
        'fix_thickness': 7.0,
    },
    {
        'name': 'relax_2x2',
        'type': 'vasp',
        'supercell': [2, 2, 1],
        'incar': incar_fine,
        'restart': None,
        'kpoints': [6, 6, 1],
        'fix_type': 'center',
        'fix_thickness': 7.0,
    },

    # Stage 5: DOS calculation (NEW!)
    {
        'name': 'dos_calculation',
        'type': 'dos',  # Uses quick_dos internally
        'structure_from': 'relax_2x2',  # Use structure from previous stage
        'scf_incar': {'encut': 700, 'ediff': 1e-6, 'ismear': 0},
        'dos_incar': {'nedos': 3000, 'lorbit': 11, 'ismear': -5},
        'dos_kpoints_spacing': 0.02,
        'retrieve': ['DOSCAR'],
    },

    # Stage 6: AIMD simulation (NEW!)
    {
        'name': 'aimd_equilibration',
        'type': 'aimd',  # Uses AIMD module internally
        'structure_from': 'relax_2x2',
        'temperature': 300,  # Kelvin
        'timestep': 1.0,     # fs
        'nsteps': 5000,
        'ensemble': 'NVT',
        'thermostat': 'nose-hoover',
    },

    # Stage 7: Production AIMD (NEW!)
    {
        'name': 'aimd_production',
        'type': 'aimd',
        'structure_from': 'aimd_equilibration',
        'restart': 'aimd_equilibration',  # Continue from equilibration
        'temperature': 300,
        'nsteps': 50000,
    },

    # Stage 8: Surface thermodynamics (FUTURE)
    {
        'name': 'surface_energy',
        'type': 'thermodynamics',
        'structure_from': 'relax_2x2',
        'bulk_structure': bulk_structure,
        'reference_energies': {...},
    },
]

result = quick_workflow(
    structure=initial_structure,
    stages=stages,
    code_label='VASP-6.5.1@cluster',
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options=options,
)
```

## Module Integration Roadmap

### Phase 1: DOS Integration (Priority: HIGH)

**Goal:** Add `type: 'dos'` stages that use `quick_dos` internally.

**Implementation:**
```python
# In quick_vasp_sequential (or new quick_workflow):

if stage['type'] == 'dos':
    # Get structure from previous stage
    structure = stage_tasks[stage['structure_from']]['vasp'].outputs.structure

    # Create DOS task using BandsWorkChain
    dos_task = wg.add_task(
        BandsTask,
        name=f'dos_{stage_name}',
        structure=structure,
        code=code,
        # ... DOS-specific inputs
    )

    # Expose DOS outputs
    setattr(wg.outputs, f'{stage_name}_dos', dos_task.outputs.dos)
    setattr(wg.outputs, f'{stage_name}_energy', ...)
```

**Stage Configuration for DOS:**
| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `type` | Yes | - | `'dos'` |
| `name` | Yes | - | Unique stage identifier |
| `structure_from` | Yes | - | Stage name to get structure from |
| `scf_incar` | Yes | - | INCAR for SCF step |
| `dos_incar` | Yes | - | INCAR for DOS step |
| `kpoints_spacing` | No | base | K-points for SCF |
| `dos_kpoints_spacing` | No | 80% of SCF | K-points for DOS |
| `retrieve` | No | `[]` | Files to retrieve |

**New Result Functions:**
```python
get_dos_stage_results(result, 'dos_calculation')
# Returns: {'energy': ..., 'dos': ArrayData, 'projectors': ...}
```

---

### Phase 2: AIMD Integration (Priority: MEDIUM)

**Goal:** Add `type: 'aimd'` stages that use the AIMD module.

**Stage Configuration for AIMD:**
| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `type` | Yes | - | `'aimd'` |
| `name` | Yes | - | Unique stage identifier |
| `structure_from` | Yes | - | Stage name for initial structure |
| `restart` | No | `None` | Stage name for AIMD restart |
| `temperature` | Yes | - | Temperature (K) |
| `timestep` | No | 1.0 | Timestep (fs) |
| `nsteps` | Yes | - | Number of MD steps |
| `ensemble` | No | `'NVT'` | `'NVE'`, `'NVT'`, `'NPT'` |
| `thermostat` | No | `'nose-hoover'` | Thermostat type |
| `incar_overrides` | No | `{}` | Additional INCAR settings |

**Outputs:**
- `{stage_name}_trajectory` - Trajectory data
- `{stage_name}_final_structure` - Final structure
- `{stage_name}_temperature_history` - T vs time
- `{stage_name}_energy_history` - E vs time

---

### Phase 3: Surface Thermodynamics (Priority: LOW)

**Goal:** Add `type: 'thermodynamics'` stages.

**Stage Configuration:**
| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `type` | Yes | - | `'thermodynamics'` |
| `slab_structure_from` | Yes | - | Stage with relaxed slab |
| `bulk_structure` | Yes | - | Bulk structure (StructureData) |
| `bulk_energy` | No | - | Pre-computed bulk energy |
| `reference_energies` | Yes | - | O2, metal references |
| `temperature_range` | No | [300, 1000] | T range for phase diagram |

**Outputs:**
- `{stage_name}_surface_energy` - γ(T, Δμ)
- `{stage_name}_phase_diagram` - Stability regions

---

### Phase 4: Additional Modules (Priority: FUTURE)

| Module | Stage Type | Description |
|--------|------------|-------------|
| Adsorption Energy | `'adsorption'` | Calculate E_ads for adsorbates |
| Electronic Structure | `'bands'` | Band structure calculation |
| Phonons | `'phonons'` | Phonon dispersion (if supported) |
| NEB | `'neb'` | Nudged elastic band for barriers |
| Bader Analysis | `'bader'` | Charge analysis |

---

## Architecture Changes

### Current Architecture
```
quick_vasp_sequential()
    └── Only handles 'vasp' type stages
```

### Proposed Architecture
```
quick_workflow()  # New unified function
    ├── StageHandler (abstract base)
    │   ├── VaspStageHandler      # Current VASP stages
    │   ├── DosStageHandler       # DOS via BandsWorkChain
    │   ├── AimdStageHandler      # AIMD module
    │   ├── ThermodynamicsHandler # Surface energy
    │   └── ... (extensible)
    │
    └── Stage validation per type
```

### Key Design Decisions

1. **Backward Compatibility**
   - `quick_vasp_sequential()` remains unchanged
   - New `quick_workflow()` function for multi-module stages
   - Stages without `type` field default to `'vasp'`

2. **Stage Handler Pattern**
   ```python
   class StageHandler(ABC):
       @abstractmethod
       def validate_config(self, stage: dict) -> None:
           """Validate stage configuration."""
           pass

       @abstractmethod
       def create_tasks(self, wg: WorkGraph, stage: dict, context: dict) -> dict:
           """Create WorkGraph tasks for this stage."""
           pass

       @abstractmethod
       def get_outputs(self) -> dict:
           """Return output sockets to expose."""
           pass
   ```

3. **Context Passing**
   ```python
   context = {
       'stage_tasks': {},           # name -> task references
       'current_structure': ...,    # Current structure socket
       'current_remote': ...,       # Current remote folder
       'code': ...,                 # VASP code
       'potential_family': ...,
       'potential_mapping': ...,
       'options': ...,
   }
   ```

4. **Output Naming Convention**
   ```
   {stage_name}_{output_type}

   Examples:
   - relax_2x2_energy
   - relax_2x2_structure
   - dos_calculation_dos
   - dos_calculation_projectors
   - aimd_production_trajectory
   - aimd_production_final_structure
   ```

---

## Implementation Plan

### Step 1: Refactor Current Code (Preparation)
- [ ] Extract VASP stage handling into `VaspStageHandler` class
- [ ] Create `StageHandler` abstract base class
- [ ] Add `type` field support (default: `'vasp'`)
- [ ] Ensure backward compatibility

### Step 2: Implement DOS Integration
- [ ] Create `DosStageHandler` class
- [ ] Add DOS stage validation
- [ ] Wire BandsWorkChain as task
- [ ] Add DOS-specific output extraction
- [ ] Add `get_dos_stage_results()` function
- [ ] Write tests
- [ ] Update documentation

### Step 3: Implement AIMD Integration
- [ ] Create `AimdStageHandler` class
- [ ] Add AIMD stage validation
- [ ] Wire AIMD WorkGraph as task
- [ ] Handle AIMD restart between stages
- [ ] Add trajectory output extraction
- [ ] Write tests
- [ ] Update documentation

### Step 4: Create `quick_workflow()` Function
- [ ] New entry point that supports all stage types
- [ ] Registry of stage handlers
- [ ] Unified validation
- [ ] Unified result extraction

---

## Open Questions

1. **How to handle different codes per stage?**
   - Some stages might need different VASP versions
   - AIMD might use CP2K instead of VASP
   - Solution: Allow `code_label` override per stage?

2. **How to handle different computational resources?**
   - DOS needs fewer resources than AIMD
   - Solution: Allow `options` override per stage?

3. **How to handle stage dependencies beyond structure?**
   - Some modules need multiple inputs (thermodynamics needs bulk + slab)
   - Solution: Allow `inputs_from` dict mapping input names to stages?

4. **How to handle conditional stages?**
   - "Only run AIMD if relaxation converged"
   - Solution: Add `condition` field with simple expressions?

5. **How to visualize complex workflows?**
   - ASCII art in terminal?
   - Export to graphviz/mermaid?

---

## Example Use Cases

### Use Case 1: Surface Characterization Pipeline
```python
stages = [
    {'type': 'vasp', 'name': 'relax', ...},
    {'type': 'dos', 'name': 'dos', 'structure_from': 'relax', ...},
    {'type': 'bands', 'name': 'bands', 'structure_from': 'relax', ...},
    {'type': 'bader', 'name': 'charges', 'structure_from': 'relax', ...},
]
```

### Use Case 2: Phase Diagram Calculation
```python
stages = [
    {'type': 'vasp', 'name': 'bulk_relax', ...},
    {'type': 'vasp', 'name': 'slab_relax', ...},
    {'type': 'thermodynamics', 'name': 'phase_diagram',
     'bulk_from': 'bulk_relax', 'slab_from': 'slab_relax', ...},
]
```

### Use Case 3: AIMD with Analysis
```python
stages = [
    {'type': 'vasp', 'name': 'relax', ...},
    {'type': 'aimd', 'name': 'equilibrate', 'structure_from': 'relax', ...},
    {'type': 'aimd', 'name': 'production', 'restart': 'equilibrate', ...},
    {'type': 'analysis', 'name': 'rdf', 'trajectory_from': 'production', ...},
]
```

### Use Case 4: Defect Formation Energy
```python
stages = [
    {'type': 'vasp', 'name': 'pristine', 'structure': pristine, ...},
    {'type': 'vasp', 'name': 'vacancy', 'structure': vacancy, ...},
    {'type': 'defect_energy', 'name': 'formation',
     'pristine_from': 'pristine', 'defect_from': 'vacancy', ...},
]
```

---

## Notes

- Start with DOS integration as proof of concept
- Keep the architecture flexible for future modules
- Maintain the "simple API" philosophy of explorer
- Each stage type should have comprehensive validation
- Error messages should clearly indicate which stage failed and why

---

## References

- Current `quick_vasp_sequential()` implementation: `workgraph.py`
- DOS module: `quick_dos()` in `workgraph.py`
- AIMD module: `teros/core/aimd/`
- Thermodynamics: `teros/core/thermodynamics.py`
- Adsorption energy: `teros/core/adsorption_energy.py`

---

*Created: 2025-01-27*
*Status: PLANNING*
*Priority: HIGH for DOS, MEDIUM for AIMD, LOW for others*
