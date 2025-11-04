# VASP Machine Learning Force Fields - Technical Knowledge

**Source**: VASP Wiki (via web search 2025-11-04)

---

## Core Concepts

### What is MLFF?

VASP's Machine Learning Force Field (MLFF) allows you to:
- Generate ML models from ab initio calculations on-the-fly
- Use trained ML models for fast MD simulations
- Improve and refine existing ML models
- Apply ML force fields without further ab initio calculations

**Key advantage**: 10-100x faster MD after initial training, while maintaining near-DFT accuracy.

---

## Critical INCAR Tags

### Primary Control: ML_LMLFF

```
ML_LMLFF = .TRUE.   # Enable MLFF (REQUIRED!)
```

**CRITICAL**: Without `ML_LMLFF = .TRUE.`, all other ML_* tags are **completely ignored** and VASP runs regular DFT calculations!

---

### Operational Modes: ML_ISTART

`ML_ISTART` controls what VASP does with ML models:

#### ML_ISTART = 0 (Start Training from Scratch)
```
ML_LMLFF = .TRUE.
ML_ISTART = 0
ML_MODE = TRAIN
```

**Behavior:**
- Starts ML training from scratch
- Runs ab initio DFT frequently (no ML model exists yet)
- Generates training data → ML_ABN file
- Saves trained model → ML_FFN file
- Intensive at start, improves as model learns

**Use when**: Starting new ML model training

---

#### ML_ISTART = 1 (Continue Training with Existing Data)
```
ML_LMLFF = .TRUE.
ML_ISTART = 1
ML_MODE = TRAIN
```

**Behavior:**
- Reads pre-existing ML_AB file (training data)
- Continues training/learning from that data
- Adds new ab initio data during MD
- Updates ML_FFN model

**Use when**: Continuing training, expanding existing dataset

**CRITICAL**: Requires ML_AB file from previous run!

---

#### ML_ISTART = 2 (Production - Prediction Only)
```
ML_LMLFF = .TRUE.
ML_ISTART = 2
```

**Behavior:**
- Reads trained ML_FFN model
- **NO ab initio calculations** (pure ML predictions!)
- **NO learning** (model is fixed)
- Very fast MD simulations

**Use when**: Running production MD with trained model

**CRITICAL**: Requires ML_FFN file from training run!

---

#### ML_ISTART = 3 (Learning from Given Data Only)
```
ML_LMLFF = .TRUE.
ML_ISTART = 3
```

**Behavior:**
- Learns only from provided ML_AB file
- **No MD time steps**
- Generates new ML_FFN from existing data
- Single-step operation

**Use when**: Post-processing training data, parameter tuning

---

#### ML_ISTART = 4 (Refit Existing Data)
```
ML_LMLFF = .TRUE.
ML_ISTART = 4
ML_MODE = REFIT
```

**Behavior:**
- Refits force field based on existing ML_AB
- Only single step executed
- **No ab initio calculations**
- **Much faster** than retraining from scratch

**Use when**: Optimizing ML hyperparameters with same data

---

## Training Control: ML_MODE

```
ML_MODE = TRAIN    # Default: On-the-fly training
ML_MODE = REFIT    # Refit existing data (faster)
```

**ML_MODE = TRAIN**:
- Used with ML_ISTART = 0 or 1
- Samples structures via MD
- Builds dataset piece by piece
- Automatically assembles training data

**ML_MODE = REFIT**:
- Used with ML_ISTART = 4
- Refits ML model with different hyperparameters
- Uses existing ML_AB data
- No new ab initio calculations

---

## Critical Files

### ML_AB (Training Data)
**Content**: Ab initio positions, forces, energies, stress tensors

**Role**:
- Input for ML_ISTART = 1, 3, 4 (requires pre-existing data)
- Output as ML_ABN from ML_ISTART = 0, 1 (new training data)

**Workflow**:
```bash
# After training run (ML_ISTART=0):
VASP generates: ML_ABN

# To continue training:
cp ML_ABN ML_AB
# Then run with ML_ISTART=1
```

---

### ML_FFN (Force Field Model)
**Content**: Neural network weights and parameters

**Role**:
- Output from training (ML_ISTART = 0, 1, 3, 4)
- Input for production (ML_ISTART = 2)

**Workflow**:
```bash
# After training run:
VASP generates: ML_FFN

# For production run:
cp ML_FFN <production_directory>/ML_FFN
# Then run with ML_ISTART=2 (reads ML_FFN automatically)
```

---

### ML_LOGFILE (Training Log)
**Content**: Training progress, errors, statistics

**Role**: Monitor ML model quality during training

**What to check**:
- Prediction errors (should decrease over time)
- Energy/force RMSE (target: < 1 meV/atom for energy)
- When training has converged

---

## Recommended Workflow

### Phase 1: Initial Training (ML_ISTART = 0)

```python
# INCAR for training
training_incar = {
    # Enable MLFF
    'ML_LMLFF': True,
    'ML_ISTART': 0,        # Start from scratch
    'ML_MODE': 'TRAIN',    # Train on-the-fly

    # MD settings (sample configuration space)
    'IBRION': 0,           # MD
    'MDALGO': 2,           # Nosé-Hoover
    'TEBEG': 300,
    'TEEND': 300,
    'NSW': 200,            # 200 training snapshots (adjust!)
    'POTIM': 1.0,          # 1 fs timestep

    # DFT settings (full accuracy for training!)
    'PREC': 'Normal',
    'ENCUT': 400,          # Full cutoff
    'EDIFF': 1e-5,
    'ISMEAR': 0,
    'SIGMA': 0.05,

    # ML hyperparameters (optional, VASP has good defaults)
    'ML_WTOTEN': 0.1,      # Weight for energy training
    'ML_WTIFOR': 1.0,      # Weight for force training
    # 'ML_WTSIF': 0.0,     # Weight for stress (default 0)

    # Output control
    'LWAVE': True,         # For restart
    'LCHARG': True,
}

# After completion:
# - Check ML_LOGFILE for convergence
# - ML_ABN contains training data
# - ML_FFN contains trained model
```

**How many steps?**
- Small systems (< 50 atoms): 100-300 steps
- Medium systems (50-200 atoms): 200-500 steps
- Large systems (> 200 atoms): 500-1000 steps
- Check ML_LOGFILE to see when errors plateau

---

### Phase 2: Refinement (Optional, ML_ISTART = 1 or 4)

**Option A: Continue training with more data**
```python
refinement_incar = {
    'ML_LMLFF': True,
    'ML_ISTART': 1,        # Continue from existing data
    'ML_MODE': 'TRAIN',
    'NSW': 500,            # Add 500 more snapshots
    # ... same DFT/MD settings
}

# Before running:
# cp ML_ABN ML_AB  (use previous training data)
```

**Option B: Refit with different hyperparameters**
```python
refit_incar = {
    'ML_LMLFF': True,
    'ML_ISTART': 4,        # Refit existing data
    'ML_MODE': 'REFIT',
    # Adjust ML hyperparameters
    'ML_WTOTEN': 1.0,      # Try different weights
    'ML_WTIFOR': 2.0,
}

# Uses existing ML_AB, generates new ML_FFN
# Much faster than retraining!
```

---

### Phase 3: Production (ML_ISTART = 2)

```python
production_incar = {
    # Enable MLFF in prediction mode
    'ML_LMLFF': True,
    'ML_ISTART': 2,        # Prediction only!

    # MD settings (same as training)
    'IBRION': 0,
    'MDALGO': 2,
    'TEBEG': 300,
    'TEEND': 300,
    'NSW': 10000,          # Much longer! ML is fast
    'POTIM': 1.0,

    # DFT settings (can potentially reduce?)
    'PREC': 'Normal',
    'ENCUT': 200,          # Question: Can we reduce this?
    'EDIFF': 1e-5,

    # Output control
    'LWAVE': False,        # No longer need for ML runs
    'LCHARG': False,
}

# Before running:
# - Copy ML_FFN from training run to this directory
# - VASP reads ML_FFN automatically
# - NO ab initio calculations (forces come from ML!)
# - Should be ~10-100x faster than DFT
```

---

## Restart Strategy for PS-TEROS

### Question: How are ML files passed between stages?

**Option 1: Automatic (Hope)**
```
Training stage output: remote_folder contains ML_ABN + ML_FFN
                                ↓
Production stage input: remote_folder (VASP finds ML files?)
```

**Test this first!** If VASP automatically finds ML files in restart directory, no code changes needed.

---

**Option 2: Explicit Copying (If Option 1 fails)**
```python
# In aimd_single_stage_scatter() or similar:

if stage_uses_mlff_prediction(stage_config):  # ML_ISTART = 2
    # Need to ensure ML_FFN is available
    # Either:
    # A) Copy explicitly from previous stage remote_folder
    # B) Store ML_FFN as separate RemoteData node
    # C) Use VASP restart mechanism (if it handles ML files)

    # Implementation depends on test results!
```

---

## Best Practices from VASP Wiki

### Training Phase

1. **Temperature sampling**:
   - Cover the temperature range you plan to simulate
   - Example: If production will be 300K, train near 300K
   - For broader range, train at 300K and 500K

2. **Configuration sampling**:
   - Ensure diverse structures in training set
   - Don't train only on equilibrated structures
   - Include some "rattled" or high-temperature configs

3. **Convergence monitoring**:
   - Watch ML_LOGFILE
   - Energy RMSE should be < 1 meV/atom
   - Force RMSE should be < 50 meV/Å
   - If errors don't decrease, need more data or better hyperparameters

4. **DFT parameters during training**:
   - Use **full accuracy** (same as regular DFT)
   - `ENCUT`, `KPOINTS`, `PREC` should be production quality
   - ML model quality = quality of training data!

---

### Production Phase

1. **Use fixed force field**:
   - Finish sampling and training first
   - Use prediction-only mode (ML_ISTART=2)
   - Don't mix training and production

2. **Temperature extrapolation**:
   - ML models work best near training conditions
   - Avoid temperatures far from training range
   - If needed, retrain at target temperature

3. **System changes**:
   - ML model is specific to system composition
   - Different stoichiometry → retrain
   - Different elements → retrain
   - Same elements, different structure → might work (test carefully!)

4. **Long simulations**:
   - ML allows 10,000+ step runs (vs 1000 DFT)
   - Watch for drift or unphysical behavior
   - Validate against DFT periodically

---

## Important Considerations

### File Sizes
- **ML_AB**: Can be large (100s of MB for long training)
- **ML_FFN**: Moderate (10-50 MB typically)
- **Impact**: If storing in AiiDA database, consider size limits

### Transferability
- ML models are **NOT universally transferable**
- Model for Ag(100) ≠ Model for Ag(111) (usually)
- Model for Ag ≠ Model for Ag₂O
- Train separate models for different systems (usually)

### Combining Models
- VASP forum discusses "Combining MLFF Output"
- Possible to merge ML_AB files from multiple runs
- Can train single model on combined data
- Advanced feature, needs investigation

---

## Key Takeaways for PS-TEROS

1. **ML_LMLFF = .TRUE. is mandatory** - without it, nothing happens!

2. **ML_ISTART controls workflow**:
   - 0 = Start training
   - 1 = Continue training
   - 2 = Production (fast!)
   - 3 = Learn from data only
   - 4 = Refit

3. **File management is critical**:
   - ML_ABN → ML_AB (rename for continuation)
   - ML_FFN must be present for ML_ISTART = 2
   - Need to test if restart_folder includes these

4. **Training needs full DFT accuracy**:
   - Can't reduce ENCUT during training
   - Production phase might allow reduction (test!)

5. **Each structure likely needs own model**:
   - Current scatter pattern (parallel slabs) works naturally
   - Each slab trains its own ML_AB/ML_FFN
   - Shared models possible but more complex

---

## Next Steps for Implementation

1. **Test restart behavior**:
   ```bash
   # Run step_XX_mlff_test.py
   # Check if ML_FFN is found in production stage
   # Look for "reading ML_FFN" or "reading ML_AB" in OUTCAR
   ```

2. **Based on test results**:
   - **Success**: Document MLFF usage, done!
   - **Failure**: Implement explicit ML file copying

3. **Add validation** (if desired):
   ```python
   def validate_mlff_stage(stage_config, has_restart):
       if not stage_config.get('ML_LMLFF'):
           return  # Not an MLFF stage

       istart = stage_config.get('ML_ISTART', 0)

       if istart in [1, 2] and not has_restart:
           raise ValueError(f"ML_ISTART={istart} requires restart folder with ML files")

       if istart == 0 and has_restart:
           warn("ML_ISTART=0 starts fresh training but restart folder provided")
   ```

4. **Create examples**:
   - Basic: Train (ML_ISTART=0) → Produce (ML_ISTART=2)
   - Full: Train → Refine → Produce
   - Multi-slab: Multiple surfaces with independent ML models

---

## References

- VASP Wiki: Machine learning force field calculations: Basics
- VASP Wiki: Best practices for machine-learned force fields
- VASP Wiki: ML_ISTART, ML_LMLFF, ML_MODE tags
- VASP Tutorials: Part 2 - Machine learning force fields
- VASP Forum: Restarting MLFF retraining, Combining MLFF Output
