# AIMD vs MLFF: Visual Comparison

## Standard AIMD Workflow (Current Implementation)

```
┌─────────────────────────────────────────────────────────────┐
│                    Standard AIMD                            │
└─────────────────────────────────────────────────────────────┘

Stage 0: Equilibration (DFT forces)
├── Input: Structure, TEBEG=300K, NSW=50
├── Process: DFT calculates forces at EVERY step
├── Output: Equilibrated structure + remote_folder
└── Time: ~1 hour (50 DFT calculations)

        ↓ (restart with WAVECAR/CHGCAR)

Stage 1: Production (DFT forces)
├── Input: Stage 0 output, TEBEG=300K, NSW=500
├── Process: DFT calculates forces at EVERY step
├── Output: MD trajectory
└── Time: ~10 hours (500 DFT calculations)

Total time: ~11 hours for 550 MD steps
All steps use full DFT → EXPENSIVE but ACCURATE
```

---

## MLFF Workflow (Proposed Implementation)

```
┌─────────────────────────────────────────────────────────────┐
│           AIMD with Machine Learning (MLFF)                 │
└─────────────────────────────────────────────────────────────┘

Stage 0: Training (DFT forces + ML learning)
├── Input: Structure, TEBEG=300K, NSW=200, ML_ISTART=0
├── Process:
│   ├── DFT calculates forces (ab initio)
│   ├── VASP trains neural network on-the-fly
│   └── Saves training data (ML_AB) and model (ML_FFN)
├── Output: Trained ML model + remote_folder
└── Time: ~4 hours (200 DFT calculations)

        ↓ (restart with WAVECAR/CHGCAR + ML_AB + ML_FFN)

Stage 1: Refinement [OPTIONAL] (Mixed DFT + ML)
├── Input: Stage 0 output, NSW=500, ML_ISTART=1
├── Process:
│   ├── ML predicts forces (fast)
│   ├── DFT refines periodically (expensive)
│   └── Model continues learning
├── Output: Improved ML model
└── Time: ~2 hours (mostly ML, some DFT)

        ↓ (restart with updated ML model)

Stage 2: Production (Pure ML forces)
├── Input: Trained model, NSW=5000, ML_ISTART=2
├── Process:
│   ├── ML predicts forces (NO DFT!)
│   ├── Can reduce ENCUT (DFT not needed)
│   └── Run very long simulations
├── Output: Long MD trajectory
└── Time: ~1 hour (5000 ML predictions)

Total time: ~7 hours for 5700 MD steps
   Training: 200 steps DFT (expensive)
   Production: 5000 steps ML (cheap!)

RESULT: 10x more MD steps in LESS time!
```

---

## Side-by-Side Comparison

### Computational Cost

| Workflow | Total Steps | DFT Steps | ML Steps | Time  | Cost     |
|----------|-------------|-----------|----------|-------|----------|
| Standard | 550         | 550       | 0        | 11 hr | High     |
| MLFF     | 5700        | 200       | 5500     | 7 hr  | Moderate |

**MLFF advantage**: 10x more sampling for LESS computational cost!

---

### File Management

#### Standard AIMD
```
Stage 0 outputs:
├── OUTCAR (forces, energies)
├── CONTCAR (final structure)
├── WAVECAR (wavefunctions) ← restart
└── CHGCAR (charge density) ← restart

Stage 1 reads:
├── WAVECAR (from Stage 0)
└── CHGCAR (from Stage 0)
```

#### MLFF
```
Stage 0 (Training) outputs:
├── OUTCAR
├── CONTCAR
├── WAVECAR ← restart
├── CHGCAR ← restart
├── ML_AB ← training data (NEW!)
└── ML_FFN ← neural network (NEW!)

Stage 1 (Production) reads:
├── WAVECAR (from Stage 0)
├── CHGCAR (from Stage 0)
├── ML_AB ← CRITICAL!
└── ML_FFN ← CRITICAL!

Question: Does restart_folder include ML files automatically?
```

---

### INCAR Configuration

#### Standard AIMD
```python
stage_0 = {
    'TEBEG': 300,
    'TEEND': 300,
    'NSW': 50,
    # Standard DFT
    'ENCUT': 400,
    'PREC': 'Normal',
    # MD settings
    'IBRION': 0,
    'MDALGO': 2,
}
```

#### MLFF Training
```python
stage_0_mlff = {
    'TEBEG': 300,
    'TEEND': 300,
    'NSW': 200,
    # Standard DFT (for training)
    'ENCUT': 400,
    'PREC': 'Normal',
    # MD settings
    'IBRION': 0,
    'MDALGO': 2,
    # MLFF training (NEW!)
    'ML_LMLFF': True,   # Enable MLFF
    'ML_ISTART': 0,     # Start training
    'ML_WTOTEN': 0.1,   # Energy weight
    'ML_WTIFOR': 1.0,   # Force weight
}
```

#### MLFF Production
```python
stage_1_mlff = {
    'TEBEG': 300,
    'TEEND': 300,
    'NSW': 5000,  # Much longer!
    # Reduced DFT (not used for forces)
    'ENCUT': 200,  # Can reduce?
    # MD settings
    'IBRION': 0,
    'MDALGO': 2,
    # MLFF production (NEW!)
    'ML_LMLFF': True,
    'ML_ISTART': 2,  # Use model, no training
}
```

---

## Code Usage Comparison

### Standard AIMD (Current)

```python
from teros.core.aimd import build_aimd_workgraph

wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[
        {'TEBEG': 300, 'NSW': 50},   # Equilibration
        {'TEBEG': 300, 'NSW': 500},  # Production
    ],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={
        'parameters': {'incar': {...}},
        # ... other settings
    },
)
```

### MLFF (Option A: No Code Changes)

```python
from teros.core.aimd import build_aimd_workgraph

wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[
        # Training: Add ML parameters
        {
            'TEBEG': 300,
            'NSW': 200,
            'ML_LMLFF': True,  # NEW
            'ML_ISTART': 0,    # NEW
        },
        # Production: Use ML model
        {
            'TEBEG': 300,
            'NSW': 5000,       # Much longer!
            'ML_LMLFF': True,  # NEW
            'ML_ISTART': 2,    # NEW
            'ENCUT': 200,      # Reduced
        },
    ],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={
        'parameters': {'incar': {...}},
        # ... other settings
    },
)

# Same function, just different INCAR tags!
# Works IF restart includes ML files automatically
```

### MLFF (Option B: New Function)

```python
from teros.core.aimd import build_mlff_workgraph

wg = build_mlff_workgraph(
    structures={'slab': structure},
    mlff_config={
        'training': {
            'temperature': 300,
            'steps': 200,
        },
        'production': {
            'temperature': 300,
            'steps': 5000,
        },
    },
    code_label='VASP6.5.0@cluster02',
    builder_inputs={...},
)

# Dedicated MLFF function
# Validates ML parameters automatically
# Handles ML file management explicitly
```

---

## Implementation Complexity

### Option A: Extend Current System
```
Code changes: None (just documentation)
Testing effort: Low (add ML INCAR tags, test restart)
Maintenance: None
User experience: Expert users only (must know VASP MLFF)

Risk: Restart might not include ML files → need Option B
```

### Option B: MLFF-Aware Function
```
Code changes: Moderate (new function + validation)
Testing effort: Moderate (validate file handling, error cases)
Maintenance: Moderate (new code to maintain)
User experience: Guided (presets, validation, error messages)

Risk: More code complexity, but safer for users
```

### Option C: Hybrid
```
Code changes: Moderate (support both approaches)
Testing effort: High (test both paths)
Maintenance: High (two interfaces)
User experience: Flexible (experts use A, beginners use B)

Risk: API complexity, but maximum flexibility
```

---

## Performance Comparison (Example: 100-atom system)

### Standard AIMD: 1000 steps

```
Time per DFT step: ~2 minutes
Total time: 2000 minutes (~33 hours)
Computational cost: High (1000 DFT calculations)
Sampling quality: Good (1000 snapshots)
```

### MLFF: 200 training + 10000 production

```
Training phase:
  - 200 DFT steps: 400 minutes (~7 hours)

Production phase:
  - 10000 ML steps: ~60 minutes (~1 hour)

Total time: 460 minutes (~8 hours)
Computational cost: Moderate (only 200 DFT calculations)
Sampling quality: Excellent (10200 snapshots)

RESULT: 10x more sampling in 1/4 the time!
```

---

## When to Use Each Approach

### Use Standard AIMD when:
- [ ] System is small (< 50 atoms)
- [ ] Need very few MD steps (< 100)
- [ ] Require absolute DFT accuracy at every step
- [ ] No plans for extensive sampling

### Use MLFF when:
- [x] System is medium-large (> 50 atoms)
- [x] Need long simulations (> 1000 steps)
- [x] Want to explore phase space extensively
- [x] Computational budget is limited
- [x] ML accuracy (1 meV/atom) is acceptable

---

## Key Decision: Which Option?

### Choose Option A if:
- You're comfortable with VASP MLFF
- Want to test MLFF quickly
- Don't need validation/guardrails
- Existing restart mechanism works

### Choose Option B if:
- Need to support novice users
- Want robust validation
- Restart doesn't include ML files
- Building production workflows

### Choose Option C if:
- Have diverse user base
- Want maximum flexibility
- Can invest in development time
- Need both expert and beginner modes

---

## Recommended Testing Path

1. **Week 1: Experiment**
   ```bash
   # Test if Option A works out of the box
   python teros/experimental/step_XX_mlff_test.py
   ```

2. **Week 2: Evaluate**
   - Did ML files propagate? → Option A
   - Files missing? → Option B
   - Decide on scope based on use cases

3. **Week 3+: Implement**
   - Option A: Document only (fast!)
   - Option B: Implement new function (moderate)
   - Option C: Hybrid approach (slow)

---

## Summary

| Aspect              | Standard AIMD | MLFF                    |
|---------------------|---------------|-------------------------|
| **Complexity**      | Simple        | Multi-stage             |
| **Speed per step**  | Slow (DFT)    | Fast (ML) after training|
| **Total time**      | High          | Lower for long runs     |
| **Accuracy**        | Exact (DFT)   | ~1 meV/atom (ML)        |
| **Sampling**        | Limited       | Extensive               |
| **Use case**        | Quick checks  | Production sampling     |
| **Learning curve**  | Low           | Moderate                |
| **File management** | Simple        | Complex (ML files)      |
| **Code changes**    | None          | TBD (test first!)       |

**Bottom line**: MLFF enables 10-100x more MD sampling at similar or lower cost, but requires careful workflow design.
