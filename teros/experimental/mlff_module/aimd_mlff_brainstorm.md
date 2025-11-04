# AIMD with Machine Learning Force Fields (MLFF) - Brainstorming Document

**Date:** 2025-11-04
**Status:** Initial Brainstorming Phase
**Related:** Standalone AIMD Module (`teros/core/aimd/`)

---

## Overview

This document explores how to extend the PS-TEROS standalone AIMD module to support VASP's Machine Learning Force Field (MLFF) calculations. MLFF combines ab initio accuracy with reduced computational cost by training on-the-fly neural network potentials during MD simulations.

---

## Current AIMD Implementation Summary

### Architecture
- **Sequential stages**: Each stage runs AIMD at specified temperature/steps
- **Automatic restart chaining**: Stage N+1 continues from Stage N
- **Parallel execution**: All structures run in parallel within each stage
- **Simple configuration**: `aimd_stages = [{'TEBEG': T, 'NSW': N}, ...]`

### Key Components
1. `build_aimd_workgraph()` - Main orchestrator
2. `aimd_single_stage_scatter()` - Runs one stage across all structures
3. `prepare_aimd_parameters()` - Merges base + stage parameters
4. Override system - Per-structure/stage/matrix INCAR customization

---

## VASP MLFF Workflow Overview

### What Makes MLFF Different?

MLFF is NOT just adding a few INCAR tags. It's a **multi-phase workflow** with distinct stages:

### Phase 1: Training (Ab Initio)
- Run standard DFT-MD to collect training data
- VASP computes forces/energies from first principles
- **Key INCAR tags:**
  - `ML_LMLFF = .TRUE.` - Enable MLFF
  - `ML_ISTART = 0` - Start training from scratch
  - `NSW = 50-500` - Training steps (depends on system)
  - All standard DFT parameters (ENCUT, PREC, etc.)

### Phase 2: Refinement (Mixed Ab Initio + ML)
- Continue training with on-the-fly learning
- VASP uses ML predictions but refines with DFT periodically
- **Key INCAR tags:**
  - `ML_LMLFF = .TRUE.`
  - `ML_ISTART = 1` - Continue from existing ML_AB/ML_FFN
  - `ML_MODE = refit` or `train` - How to update the model
  - `ML_WTOTEN/ML_WTIFOR/ML_WTSIF` - Weight parameters for energy/forces/stress
  - Possibly reduced `NSW` compared to training

### Phase 3: Production (Pure ML)
- Run long MD simulations using trained ML potential only
- No DFT calculations (much faster!)
- **Key INCAR tags:**
  - `ML_LMLFF = .TRUE.`
  - `ML_ISTART = 2` - Read existing ML_AB/ML_FFN, no training
  - `NSW = 10000+` - Can run much longer simulations
  - Can reduce DFT parameters (ENCUT, KPOINTS) since ML model is used

### Critical Files
- **ML_AB** - Training data (positions, forces, energies)
- **ML_FFN** - Neural network weights
- These must be carried between stages via WAVECAR/CHGCAR restart mechanism

---

## Design Challenges

### 1. **File Management**
- **Challenge**: ML_AB and ML_FFN files must persist across stages
- **Current system**: Uses `remote_folder` for WAVECAR/CHGCAR restart
- **Question**: Does VASP automatically read ML_AB/ML_FFN from restart folder?
  - If YES: Existing restart mechanism might work!
  - If NO: Need explicit file copying between stages

### 2. **Stage Configuration**
- **Challenge**: Each MLFF phase has different INCAR requirements
- **Current system**: `aimd_stages = [{'TEBEG': T, 'NSW': N}, ...]`
- **MLFF needs:**
  ```python
  mlff_stages = [
      # Training phase
      {'TEBEG': 300, 'NSW': 200, 'ML_ISTART': 0, 'ML_LMLFF': True},

      # Refinement phase
      {'TEBEG': 300, 'NSW': 500, 'ML_ISTART': 1, 'ML_MODE': 'refit'},

      # Production phase (pure ML)
      {'TEBEG': 300, 'NSW': 5000, 'ML_ISTART': 2, 'ENCUT': 200},  # Reduced ENCUT!
  ]
  ```

### 3. **Validation Requirements**
- **Challenge**: MLFF has complex interdependencies between parameters
- **Examples:**
  - `ML_ISTART = 1 or 2` requires ML_AB/ML_FFN from previous stage
  - `ML_ISTART = 0` should NOT have ML files (clean start)
  - `ML_ISTART = 2` should have reduced DFT computational parameters
  - Certain ML parameters only valid with specific `ML_ISTART` values

### 4. **Error Handling**
- **Challenge**: MLFF training can fail (insufficient data, bad hyperparameters)
- **Current system**: No error recovery within workflow
- **Possible solutions:**
  - Validate ML model quality between stages
  - Automatic parameter adjustment on failure
  - Fallback to pure DFT if ML training fails

---

## Proposed Design Options

### Option A: Extend Existing Stage System (Conservative)

**Idea**: Keep current architecture, just allow more INCAR parameters in `aimd_stages`

```python
mlff_stages = [
    # Stage 0: Training
    {
        'TEBEG': 300,
        'NSW': 200,
        'ML_LMLFF': True,
        'ML_ISTART': 0,
        # ... other ML params
    },
    # Stage 1: Refinement
    {
        'TEBEG': 300,
        'NSW': 500,
        'ML_ISTART': 1,
        'ML_MODE': 'refit',
    },
    # Stage 2: Production
    {
        'TEBEG': 300,
        'NSW': 5000,
        'ML_ISTART': 2,
        'ENCUT': 200,  # Can reduce DFT parameters!
    },
]

wg = build_aimd_workgraph(
    structures={'slab1': s1},
    aimd_stages=mlff_stages,  # Just works with MLFF parameters!
    code_label='VASP6.5.0@cluster02',
    builder_inputs=base_config,
)
```

**Pros:**
- Minimal code changes
- Users already understand the stage concept
- Existing restart mechanism might handle ML_AB/ML_FFN automatically

**Cons:**
- No validation of MLFF-specific logic (ML_ISTART sequence, file dependencies)
- No guardrails against invalid MLFF configurations
- Users must know VASP MLFF details

**Implementation:**
1. Test if `remote_folder` restart automatically includes ML files
2. If yes: NO CODE CHANGES needed! Just document MLFF usage
3. If no: Modify restart to explicitly copy ML_AB/ML_FFN

---

### Option B: MLFF-Aware Stage System (Moderate)

**Idea**: Add MLFF-specific stage types with validation

```python
from teros.core.aimd import build_mlff_workgraph

mlff_config = {
    'training': {
        'temperature': 300,
        'steps': 200,
        'ml_params': {
            'ML_ISTART': 0,
            # Auto-set: ML_LMLFF=True
        }
    },
    'refinement': {
        'temperature': 300,
        'steps': 500,
        'ml_params': {
            'ML_ISTART': 1,
            'ML_MODE': 'refit',
        }
    },
    'production': {
        'temperature': 300,
        'steps': 5000,
        'ml_params': {
            'ML_ISTART': 2,
        },
        'reduced_dft': True,  # Flag to reduce ENCUT, KPOINTS automatically
    },
}

wg = build_mlff_workgraph(
    structures={'slab1': s1},
    mlff_config=mlff_config,
    code_label='VASP6.5.0@cluster02',
    builder_inputs=base_config,
)
```

**Pros:**
- Explicit MLFF workflow structure
- Can validate ML_ISTART sequence
- Can auto-reduce DFT parameters in production phase
- Clear separation from standard AIMD

**Cons:**
- New function to maintain
- Less flexible (what if user wants custom workflow?)
- More complex API

**Implementation:**
1. Create `build_mlff_workgraph()` wrapper
2. Validate MLFF stage sequence internally
3. Call existing `build_aimd_workgraph()` with expanded stages
4. Add MLFF-specific file management if needed

---

### Option C: Hybrid Approach (Flexible)

**Idea**: Allow both simple and advanced MLFF workflows

```python
# Simple: Let users specify raw MLFF stages (Option A)
wg = build_aimd_workgraph(
    aimd_stages=[
        {'TEBEG': 300, 'NSW': 200, 'ML_LMLFF': True, 'ML_ISTART': 0},
        {'TEBEG': 300, 'NSW': 5000, 'ML_ISTART': 2},
    ],
    # ... standard parameters
)

# Advanced: Use MLFF preset with validation (Option B)
wg = build_aimd_workgraph(
    aimd_stages='mlff_preset',  # Special string triggers preset
    mlff_config={
        'training_steps': 200,
        'refinement_steps': 500,
        'production_steps': 5000,
        'temperature': 300,
        'ml_algorithm': 'mace',  # or 'nequip'
    },
    # ... standard parameters
)
```

**Pros:**
- Flexibility for advanced users (Option A)
- Convenience for beginners (Option B)
- Single API surface

**Cons:**
- More complex implementation
- Harder to document

---

## Recommended Approach

### Phase 1: Experiment (Current Phase)
**Goal**: Determine if existing system already supports MLFF

**Steps:**
1. Create test script in `examples/vasp/step_XX_mlff_test.py`
2. Use current `build_aimd_workgraph()` with MLFF INCAR tags
3. Check if:
   - ML_AB/ML_FFN files are automatically included in restart
   - MLFF stages run successfully
   - Production stage (ML_ISTART=2) can read trained model

**If successful → Option A (document only, no code changes!)**
**If fails → Continue to Phase 2**

### Phase 2: Minimal Implementation
**Goal**: Make MLFF work with minimal changes

**Steps:**
1. If file management is the issue:
   - Modify `aimd_single_stage_scatter()` to detect MLFF mode
   - Explicitly copy ML_AB/ML_FFN in restart mechanism
2. Add basic validation:
   - Warn if ML_ISTART=1/2 but no restart folder
   - Warn if ML_ISTART=0 with restart folder (unexpected training start)

### Phase 3: Advanced Features (Optional)
**Goal**: Add convenience and safety

**Steps:**
1. Create `build_mlff_workgraph()` helper function (Option B)
2. Add MLFF presets:
   - `'mlff_basic'`: training → production
   - `'mlff_full'`: training → refinement → production
   - `'mlff_annealing'`: temperature ramping with MLFF
3. Add MLFF-specific validation:
   - Check ML_ISTART sequence
   - Validate ML hyperparameters
   - Auto-adjust DFT parameters in production stage

---

## Technical Questions to Resolve

### Q1: File Management
**Question**: Does VASP read ML_AB/ML_FFN from the same directory as WAVECAR/CHGCAR when restarting?

**How to test:**
```python
# Stage 0: Train MLFF
stage_0_incar = {
    'ML_LMLFF': True,
    'ML_ISTART': 0,
    'NSW': 100,
}

# Stage 1: Continue from stage 0 (standard restart)
stage_1_incar = {
    'ML_LMLFF': True,
    'ML_ISTART': 1,  # Should read ML_AB/ML_FFN from restart folder
    'NSW': 100,
}

# Run both stages with existing restart mechanism
# Check if Stage 1 successfully reads ML files
```

**Expected result**: If VASP finds ML_AB/ML_FFN in the restart directory automatically, no code changes needed!

---

### Q2: ENCUT/KPOINTS Reduction in Production
**Question**: How much can we reduce DFT parameters in ML_ISTART=2 mode?

**Considerations:**
- In pure ML mode (ML_ISTART=2), forces come from ML model, not DFT
- But VASP still needs ENCUT/KPOINTS for initial wavefunction?
- Or can we set minimal values?

**How to test:**
```python
# Production stage with minimal DFT
production_incar = {
    'ML_LMLFF': True,
    'ML_ISTART': 2,
    'NSW': 1000,
    'ENCUT': 200,  # Much lower than training (e.g., 400)
    # Or even remove ENCUT entirely?
}
```

**Research**: Check VASP wiki / examples for recommended production settings

---

### Q3: Error Handling
**Question**: How to detect if MLFF training failed?

**Possible indicators:**
- Check ML_LOGFILE for errors
- Monitor ML prediction accuracy (if VASP outputs this)
- Compare ML forces vs DFT forces in refinement stage

**Implementation idea:**
```python
@task
def validate_mlff_training(remote_folder):
    """
    Parse ML_LOGFILE from VASP output.
    Raise error if training failed or accuracy too low.
    """
    # Read ML_LOGFILE
    # Check for error messages
    # Parse accuracy metrics
    # Return success/failure flag
```

---

### Q4: Multi-Structure MLFF
**Question**: Can we train one ML model for multiple structures simultaneously?

**Scenario:**
```python
structures = {
    'slab_100': structure1,
    'slab_110': structure2,
    'slab_111': structure3,
}
```

**Options:**
1. **Independent models**: Each structure trains its own ML_AB/ML_FFN
   - Current system supports this naturally
   - Each structure → separate VASP calculation → separate ML files
2. **Shared model**: Train single ML model on all structures combined
   - Requires merging ML_AB files across structures
   - Not sure if VASP supports this workflow?
   - More complex implementation

**Recommendation**: Start with independent models (Option 1)

---

## Example MLFF Workflow

### Scenario: Train MLFF for Ag surface AIMD

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
STEP XX: AIMD with Machine Learning Force Fields (MLFF)

Demonstrates MLFF training workflow:
1. Training phase: Collect ab initio data (200 steps)
2. Refinement phase: Continue training with refit (500 steps)
3. Production phase: Pure ML simulation (5000 steps)

Material: Ag (111) surface
"""

from aiida import load_profile, orm
from teros.core.aimd import build_aimd_workgraph
from ase.io import read

load_profile('presto')

# Load Ag structure
ag_structure = orm.StructureData(ase=read('structures/Ag.cif'))

# MLFF stages using VASP-native parameters
mlff_stages = [
    # Stage 0: Training (ab initio)
    {
        'TEBEG': 300,
        'TEEND': 300,
        'NSW': 200,  # Collect 200 training snapshots
        'ML_LMLFF': True,
        'ML_ISTART': 0,  # Start training from scratch
        # ML hyperparameters (optional, VASP has defaults)
        'ML_WTOTEN': 0.1,  # Weight for energy training
        'ML_WTIFOR': 1.0,  # Weight for force training
    },

    # Stage 1: Refinement (mixed ab initio + ML)
    {
        'TEBEG': 300,
        'TEEND': 300,
        'NSW': 500,
        'ML_LMLFF': True,
        'ML_ISTART': 1,  # Continue from previous ML_AB/ML_FFN
        'ML_MODE': 'refit',  # Refit existing model
    },

    # Stage 2: Production (pure ML)
    {
        'TEBEG': 300,
        'TEEND': 300,
        'NSW': 5000,  # Long simulation with ML forces
        'ML_LMLFF': True,
        'ML_ISTART': 2,  # Use trained model, no more training
        # Question: Can we reduce ENCUT here?
        'ENCUT': 200,  # Reduced from training value (400)
    },
]

# Base DFT parameters (used in training/refinement)
builder_inputs = {
    'parameters': {
        'incar': {
            'PREC': 'Normal',
            'ENCUT': 400,  # Full DFT in training
            'EDIFF': 1e-5,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 0,
            'MDALGO': 2,
            'POTIM': 1.0,  # 1 fs timestep
            'LWAVE': True,  # Need for restart
            'LCHARG': True,
        }
    },
    'kpoints_spacing': 0.5,
    'potential_family': 'PBE',
    'potential_mapping': {'Ag': 'Ag'},
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    },
    'clean_workdir': False,  # Keep ML files!
}

# Build MLFF workgraph (using existing function!)
wg = build_aimd_workgraph(
    structures={'ag_slab': ag_structure},
    aimd_stages=mlff_stages,
    code_label='VASP6.5.0@cluster02',
    builder_inputs=builder_inputs,

    # Optional: Create supercell for better sampling
    supercell_specs={'ag_slab': [2, 2, 1]},

    name='MLFF_Training_Ag',
)

print(f"Submitting MLFF workflow...")
wg.submit(wait=False)
print(f"WorkGraph PK: {wg.pk}")
print(f"\nExpected stages:")
print(f"  1. Training (200 steps, DFT): Collect training data")
print(f"  2. Refinement (500 steps, mixed): Improve ML model")
print(f"  3. Production (5000 steps, ML): Fast simulation")
print(f"\nMonitor with: verdi process status {wg.pk}")
```

---

## Implementation Checklist

### Phase 1: Testing
- [ ] Create `examples/vasp/step_XX_mlff_basic.py` test script
- [ ] Test 3-stage MLFF workflow (training → refinement → production)
- [ ] Verify ML_AB/ML_FFN files are preserved in restart
- [ ] Check if ML_ISTART=2 stage successfully uses trained model
- [ ] Document findings

### Phase 2: Minimal Changes (If Needed)
- [ ] Modify `aimd_single_stage_scatter()` to detect MLFF mode
- [ ] Ensure ML files copied in restart mechanism (if not automatic)
- [ ] Add basic warnings for common MLFF mistakes
- [ ] Update documentation with MLFF examples

### Phase 3: Advanced Features (Optional)
- [ ] Create `build_mlff_workgraph()` helper function
- [ ] Add MLFF validation utilities
- [ ] Implement MLFF presets
- [ ] Add ML model quality checking between stages
- [ ] Create analysis tools for ML_LOGFILE

---

## Open Questions for Discussion

1. **Scope**: Should MLFF support be:
   - Fully automatic (presets, validation) → More work, easier for users
   - Manual (just document INCAR tags) → Less work, requires VASP expertise

2. **Architecture**: Which option resonates most?
   - Option A: No code changes, just documentation
   - Option B: New `build_mlff_workgraph()` function
   - Option C: Hybrid approach (support both)

3. **Validation**: How strict should MLFF validation be?
   - No validation (trust user knows VASP)
   - Basic warnings (ML_ISTART sequence check)
   - Strict validation (block invalid configs)

4. **Error handling**: Should workflow auto-retry MLFF training failures?
   - Yes: More robust, but complex logic
   - No: User must debug manually

5. **Multi-structure**: Support shared ML models across structures?
   - Phase 1: Independent models only (simple)
   - Phase 2: Shared model option (complex, research needed)

---

## Next Steps

1. **Read VASP documentation** (you have the links!):
   - Understand file dependencies (ML_AB, ML_FFN, restart behavior)
   - Learn recommended ML hyperparameters
   - Check production stage DFT parameter requirements

2. **Create test script** (`step_XX_mlff_basic.py`):
   - Implement simplest possible MLFF workflow
   - Use existing `build_aimd_workgraph()` with MLFF INCAR tags
   - See what works and what breaks

3. **Decide on approach** based on test results:
   - If test works → Option A (document only)
   - If fails → Minimal fixes (Phase 2)
   - If users need help → Advanced features (Phase 3)

4. **Discuss priorities**:
   - What's the main use case? (Training once vs routine MLFF runs)
   - Who are the users? (VASP experts vs beginners)
   - Timeline? (Quick prototype vs robust production tool)

---

## References

- VASP MLFF Wiki: https://www.vasp.at/wiki/index.php/Machine_learning_force_field_calculations:_Basics
- VASP MLFF Best Practices: https://vasp.at/wiki/Best_practices_for_machine-learned_force_fields
- PS-TEROS AIMD Module: `/home/user/PS-TEROS/teros/core/aimd/`
- Example AIMD Script: `/home/user/PS-TEROS/examples/vasp/step_18_aimd_standalone.py`
