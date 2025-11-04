# MLFF Module Design - Questions for You

Based on your PS-TEROS architecture, I need your input on several design decisions. I've organized these by priority and linked them to your existing code patterns.

---

## 🔴 HIGH PRIORITY: Architecture Decisions

### Q1: Module Location & Organization

Looking at your existing code structure:
```
teros/core/
├── aimd_functions.py         # AIMD tasks
├── aimd_cp2k.py              # CP2K AIMD
└── aimd/                     # Standalone AIMD module
    ├── __init__.py
    ├── workgraph.py
    ├── tasks.py
    └── utils.py
```

**Question 1a**: Where should MLFF code live?

**Option A: Extend existing AIMD module**
```
teros/core/aimd/
├── workgraph.py          # Add build_mlff_workgraph() here
├── tasks.py              # Keep existing
├── utils.py              # Add MLFF validation here
└── mlff_utils.py         # NEW: MLFF-specific utilities
```

**Option B: Create separate MLFF module**
```
teros/core/mlff/          # NEW module
├── __init__.py
├── workgraph.py          # build_mlff_workgraph()
├── tasks.py              # MLFF-specific tasks
├── utils.py              # MLFF validation
└── ml_file_handlers.py   # ML_AB/ML_FFN management
```

**Option C: Integration layer**
```
teros/core/aimd/
├── workgraph.py          # Existing AIMD
├── mlff.py               # NEW: MLFF wrapper (calls workgraph.py)
└── utils.py              # Shared utilities
```

**Your preference?** A, B, or C? Why?

---

### Q2: API Design Philosophy

Your current AIMD module uses this pattern:
```python
wg = build_aimd_workgraph(
    structures={'slab1': s1},
    aimd_stages=[{'TEBEG': 300, 'NSW': 100}, ...],  # Simple list
    builder_inputs={...},                            # Shared config
)
```

**Question 2a**: Should MLFF use the same pattern (just add ML tags)?

**Option A: Keep it simple (minimal API change)**
```python
# User adds ML tags to aimd_stages
mlff_stages = [
    {'TEBEG': 300, 'NSW': 200, 'ML_LMLFF': True, 'ML_ISTART': 0},
    {'TEBEG': 300, 'NSW': 5000, 'ML_LMLFF': True, 'ML_ISTART': 2},
]

wg = build_aimd_workgraph(  # Same function!
    structures={'slab1': s1},
    aimd_stages=mlff_stages,
    builder_inputs={...},
)
```

**Option B: Explicit MLFF function**
```python
# Dedicated MLFF function with guided config
wg = build_mlff_workgraph(  # New function
    structures={'slab1': s1},
    training_config={'temperature': 300, 'steps': 200},
    production_config={'temperature': 300, 'steps': 5000},
    builder_inputs={...},
)
```

**Option C: Mode flag**
```python
# Use mode to switch behavior
wg = build_aimd_workgraph(
    structures={'slab1': s1},
    aimd_stages=[...],
    mode='mlff',  # NEW flag
    mlff_settings={...},
    builder_inputs={...},
)
```

**Your preference?** A, B, or C?

**Follow-up**: Who are your primary users?
- [ ] You only (expert mode is fine)
- [ ] Your research group (moderate guidance needed)
- [ ] External users (need lots of validation/docs)

---

### Q3: Override System Interaction

Your AIMD module has a sophisticated override system:
```python
build_aimd_workgraph(
    structure_overrides={slab_name: {...}},     # Per-structure
    stage_overrides={stage_idx: {...}},         # Per-stage
    matrix_overrides={(slab, stage): {...}},    # Per-combination
)
```

**Question 3a**: Should MLFF stages support the same overrides?

**Scenario**: Different structures need different ML hyperparameters
```python
# Structure 1: Default ML settings
# Structure 2: Needs stricter ML convergence
structure_overrides = {
    'slab2': {
        'parameters': {
            'incar': {
                'ML_WTOTEN': 1.0,  # Higher energy weight
                'ML_WTIFOR': 2.0,  # Higher force weight
            }
        }
    }
}
```

**Question 3b**: Should ML_ISTART be overridable per structure?

**Scenario**: Train ML model on slab1, use it for slab2
```python
matrix_overrides = {
    ('slab1', 0): {'ML_ISTART': 0},  # slab1 trains
    ('slab2', 0): {'ML_ISTART': 2},  # slab2 uses pre-trained model from slab1
}
```

**Is this a use case you need?** YES / NO / MAYBE

---

## 🟡 MEDIUM PRIORITY: Workflow Behavior

### Q4: Multi-Structure MLFF Strategy

Your AIMD module runs structures **in parallel** within each stage:
```python
structures = {'slab1': s1, 'slab2': s2, 'slab3': s3}
# Stage 0: All 3 slabs run simultaneously (limited by max_concurrent_jobs)
# Stage 1: All 3 slabs run simultaneously (restart from Stage 0)
```

**Question 4a**: For MLFF, should each structure train its own ML model?

**Scenario**:
```python
structures = {
    'Ag_100': ag_100_structure,
    'Ag_110': ag_110_structure,
    'Ag_111': ag_111_structure,
}
```

**Independent models (current AIMD pattern)**:
- Each slab → separate VASP calc → separate ML_AB/ML_FFN
- 3 different ML models (one per surface)
- Simple to implement (already supported!)

**Shared model (requires new logic)**:
- All slabs → train one combined ML model
- Need to merge ML_AB files across slabs
- More complex, but single model might generalize better

**Your use case:**
- [ ] Independent models (surfaces are different, need separate models)
- [ ] Shared model (surfaces are similar, one model for all)
- [ ] Don't know yet

---

### Q5: Stage Sequencing

Your AIMD module chains stages sequentially:
```
Stage 0 (all slabs) → Stage 1 (all slabs) → Stage 2 (all slabs)
         ↓                     ↓                     ↓
    remote_folders        remote_folders        remote_folders
```

**Question 5a**: Should MLFF stages be strictly sequential (like AIMD)?

**Scenario**: Separate training and production phases
```
Training stage (200 steps, all slabs)
         ↓
    Wait for ALL slabs to finish training
         ↓
Production stage (5000 steps, all slabs)
```

**Alternative**: Allow per-slab progression
```
Slab1: Train → Done → Production (start immediately)
Slab2: Train → (still running...)
Slab3: Train → Done → Production (start immediately)
```

**Your preference:**
- [ ] Sequential (wait for all) - matches AIMD pattern
- [ ] Per-slab progression - faster but more complex
- [ ] Don't care / let me test both

---

### Q6: File Management Strategy

**Question 6a**: How important is provenance tracking for ML files?

Your system currently tracks:
- WAVECAR/CHGCAR via `remote_folder` (RemoteData node)
- All calculations in AiiDA provenance graph

For MLFF, we need to track:
- ML_AB (training data) - can be large (100s of MB)
- ML_FFN (neural network) - moderate size (10-50 MB)
- ML_LOGFILE (training log) - small

**Options:**

**Option 1: Implicit (rely on remote_folder)**
```python
# Assume VASP finds ML files in restart directory
# Pro: No code changes
# Con: No explicit tracking in AiiDA
```

**Option 2: Explicit (separate RemoteData nodes)**
```python
# Explicitly copy ML_AB/ML_FFN to next stage
# Pro: Clear provenance
# Con: More complex
```

**Option 3: Store as AiiDA data nodes**
```python
# Store ML_AB/ML_FFN as FolderData or SinglefileData
# Pro: Full provenance, can retrieve later
# Con: Large database size
```

**Your preference:**
- [ ] Implicit (if it works, keep it simple)
- [ ] Explicit (clear what's happening)
- [ ] Store in AiiDA (full tracking, even if large)

---

### Q7: Validation & Error Handling

Your AIMD module has validation:
```python
validate_stage_sequence(aimd_stages)  # Checks TEBEG/NSW present
validate_supercell_spec(spec)         # Checks [nx, ny, nz] format
```

**Question 7a**: What level of MLFF validation do you want?

**Level 1: Minimal (warnings only)**
```python
if 'ML_ISTART' in stage_config and stage_config['ML_ISTART'] > 0:
    if not restart_folders:
        print("Warning: ML_ISTART=1/2 but no restart folder. ML files may be missing.")
```

**Level 2: Moderate (block obvious errors)**
```python
if 'ML_ISTART' in stage_config:
    if stage_config['ML_ISTART'] == 0 and restart_folders:
        raise ValueError("ML_ISTART=0 starts fresh training, but restart_folder provided. Remove restart or use ML_ISTART=1")
```

**Level 3: Strict (validate entire MLFF workflow)**
```python
validate_mlff_stage_sequence(mlff_stages)
# - Check ML_ISTART progression (0 → 1 → 2, or 0 → 2)
# - Validate ML hyperparameters
# - Check ENCUT reduction makes sense
# - Verify file dependencies
```

**Your preference:**
- [ ] Level 1 (trust me, I know VASP)
- [ ] Level 2 (catch obvious mistakes)
- [ ] Level 3 (prevent all invalid configs)

**Question 7b**: If MLFF training fails, what should happen?

**Scenario**: Stage 0 completes but ML model quality is poor

- [ ] Continue anyway (user's responsibility to check)
- [ ] Warn user but continue
- [ ] Abort workflow
- [ ] Retry with different hyperparameters (auto-tuning)

---

## 🟢 LOW PRIORITY: Advanced Features

### Q8: Builder Integration

Your builders provide defaults:
```python
from teros.core.builders.aimd_builder import get_aimd_defaults

aimd_params = get_aimd_defaults(
    energy_cutoff=400,
    timestep=1.0,
    mdalgo=2,  # Nosé-Hoover
)
```

**Question 8a**: Should we create a MLFF builder?

```python
from teros.core.builders.mlff_builder import get_mlff_defaults

mlff_params = get_mlff_defaults(
    energy_cutoff=400,
    timestep=1.0,
    ml_algorithm='default',  # or 'mace', 'nequip', etc.
    training_steps=200,
    production_steps=5000,
)

# Returns complete mlff_stages list + builder_inputs
```

**Would this be useful?** YES / NO / LATER

---

### Q9: Testing Strategy

Your test scripts are in `examples/vasp/`:
```
step_18_aimd_standalone.py
step_19_aimd_with_overrides.py
```

**Question 9a**: What MLFF examples do you need?

- [ ] `step_XX_mlff_basic.py` - Simple training → production
- [ ] `step_XX_mlff_full.py` - Training → refinement → production
- [ ] `step_XX_mlff_multi_structure.py` - Multiple slabs with MLFF
- [ ] `step_XX_mlff_overrides.py` - Different ML params per structure
- [ ] All of the above

**Question 9b**: Do you want unit tests (pytest) or just example scripts?

Your current AIMD has both:
```
teros/core/aimd/test_utils.py
teros/core/aimd/test_overrides.py
examples/vasp/step_18_aimd_standalone.py
```

- [ ] Unit tests (for production code)
- [ ] Example scripts only (faster development)
- [ ] Both (most robust)

---

### Q10: Documentation & Workflow Presets

Your workflow might benefit from MLFF presets (like your main workflow has):

**Question 10a**: Would presets be useful?

```python
# Instead of manually configuring stages
wg = build_mlff_workgraph(
    structures={...},
    preset='basic',  # Trains for 200 steps, runs 5000 production
    # OR
    preset='full',   # Trains 200, refines 500, runs 5000 production
    # OR
    preset='annealing',  # Temperature ramping with MLFF
    # OR
    preset='custom', # Full manual control
    mlff_stages=[...],  # Only if preset='custom'
)
```

**Would you use presets?** YES / NO / MAYBE

---

## 📋 Summary: Please Answer

To design the MLFF module that fits YOUR needs, please answer:

### Architecture (Q1-Q3):
1. **Module location**: A (extend AIMD), B (separate MLFF), or C (wrapper)?
2. **API style**: A (simple), B (explicit), or C (mode flag)?
3. **Do you need override system** for MLFF parameters?

### Workflow (Q4-Q7):
4. **Multi-structure**: Independent models or shared models?
5. **Stage sequencing**: Sequential (wait for all) or per-slab?
6. **File management**: Implicit, explicit, or store in AiiDA?
7. **Validation level**: Minimal, moderate, or strict?

### Features (Q8-Q10):
8. **MLFF builder**: Useful now, later, or not needed?
9. **Testing**: Unit tests, examples, or both?
10. **Presets**: Would you use them?

### Use Case Context:
- **What's your primary goal with MLFF?**
  - [ ] Research/exploration
  - [ ] Production calculations
  - [ ] Method development
  - [ ] Teaching/learning

- **Timeline?**
  - [ ] Need it working this week (quick prototype)
  - [ ] 1-2 weeks (robust implementation)
  - [ ] 1+ month (full-featured)

- **System sizes you typically run?**
  - [ ] Small (< 50 atoms)
  - [ ] Medium (50-200 atoms)
  - [ ] Large (> 200 atoms)

- **Typical simulation length?**
  - [ ] Short (< 500 steps)
  - [ ] Medium (500-5000 steps)
  - [ ] Long (> 5000 steps)

---

## 🚀 After You Answer

Based on your responses, I'll:

1. **Design the module architecture** (where code lives, API design)
2. **Implement the chosen approach** (Option A, B, or C)
3. **Create test scripts** matching your use cases
4. **Write documentation** at the appropriate level
5. **Set up validation** (minimal, moderate, or strict)

---

## Quick Reaction Mode

If you want to get started FAST, here's my recommendation:

**Phase 1 (This week):**
- Test if MLFF works with current code (Option A, simple API)
- Run `step_XX_mlff_test.py` in the mlff_module folder
- See if restart automatically includes ML files

**Phase 2 (If Phase 1 works):**
- Document MLFF usage
- Create 2-3 example scripts
- You're done! 🎉

**Phase 2 (If Phase 1 fails):**
- Implement explicit ML file management
- Add basic validation
- Create MLFF helper function (Option B or C)

**Sound good?** Let me know which questions you want to tackle first, or if you want to just run the test and see what happens!
