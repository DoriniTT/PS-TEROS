# MLFF Module - Development Folder

**Purpose**: Brainstorming and experimental implementation of VASP Machine Learning Force Field (MLFF) support for PS-TEROS

**Status**: 🚧 Experimental / Planning Phase

---

## 📁 Contents

### 1. **QUESTIONS_FOR_YOU.md** (START HERE!)
   **What**: Detailed questionnaire about your needs and preferences
   **Why**: Your answers will guide the implementation approach
   **Action**: Please read and answer the questions!

### 2. **VASP_MLFF_KNOWLEDGE.md**
   **What**: Technical reference about VASP MLFF (from VASP Wiki)
   **Content**:
   - ML_ISTART modes (0, 1, 2, 3, 4)
   - ML_LMLFF and ML_MODE tags
   - File management (ML_AB, ML_FFN)
   - Best practices from VASP
   - Recommended workflows

### 3. **MLFF_COMPARISON.md**
   **What**: Visual comparison of Standard AIMD vs MLFF workflows
   **Content**:
   - Side-by-side workflow diagrams
   - Computational cost analysis
   - File management differences
   - Code usage examples
   - Performance benchmarks

### 4. **MLFF_DISCUSSION.md**
   **What**: Quick summary and discussion guide
   **Content**:
   - Three implementation options (A, B, C)
   - Critical technical questions
   - Next actions and timeline
   - Decision points

### 5. **aimd_mlff_brainstorm.md**
   **What**: Comprehensive technical brainstorming
   **Content**:
   - Design challenges
   - Proposed design options
   - Implementation strategies
   - Open questions
   - Full analysis

### 6. **step_XX_mlff_test.py** (RUNNABLE TEST!)
   **What**: Experimental test script
   **Purpose**: Test if MLFF works with existing AIMD code
   **Workflow**:
   - Stage 0: Training (ML_ISTART=0, 20 steps)
   - Stage 1: Production (ML_ISTART=2, 50 steps)
   **Action**: Review, enable submission, and run!

---

## 🚀 Quick Start

### Option 1: Answer Questions First (Recommended)
```bash
# Read the questionnaire
cat QUESTIONS_FOR_YOU.md

# Think about:
# - Where should code live?
# - What API style do you prefer?
# - How much validation do you need?
# - What are your use cases?
```

### Option 2: Test First, Decide Later (Fast)
```bash
# Run the test script to see if MLFF works out of the box
cd /home/user/PS-TEROS/teros/experimental/mlff_module

# Review the test script
cat step_XX_mlff_test.py

# Edit to enable submission (uncomment wg.submit() line)
# Then run:
source ~/envs/aiida/bin/activate
python step_XX_mlff_test.py

# Monitor:
verdi process status <PK>

# Critical check: Did Stage 1 find ML_FFN from Stage 0?
# Look for "reading ML_FFN" in Stage 1 OUTCAR
```

---

## 📊 Decision Tree

```
                    START HERE
                         |
                         ↓
          ┌──────────────────────────┐
          │  Run step_XX_mlff_test.py │
          └──────────────┬─────────────┘
                         |
                         ↓
              Does ML_ISTART=2 work?
                    /        \
                 YES          NO
                  |            |
                  ↓            ↓
        ┌─────────────┐    ┌──────────────────┐
        │  OPTION A   │    │  OPTION B or C   │
        │  Document   │    │  Implement       │
        │  only!      │    │  ML file mgmt    │
        └─────────────┘    └──────────────────┘
                  |            |
                  ↓            ↓
        ┌─────────────────────────────┐
        │  Answer QUESTIONS_FOR_YOU   │
        │  to finalize design         │
        └─────────────────────────────┘
                  |
                  ↓
          Implementation Phase
```

---

## 🎯 Implementation Options Summary

### Option A: Minimal (No Code Changes)
- **What**: Use existing `build_aimd_workgraph()`, just add ML INCAR tags
- **Pros**: Zero code changes, test immediately
- **Cons**: No validation, expert users only
- **When**: If test shows ML files work with restart mechanism

### Option B: MLFF-Aware Function
- **What**: Create `build_mlff_workgraph()` with validation
- **Pros**: Guided workflow, safe for beginners
- **Cons**: More code, less flexible
- **When**: If restart doesn't handle ML files, or need validation

### Option C: Hybrid
- **What**: Support both simple (A) and guided (B) approaches
- **Pros**: Maximum flexibility
- **Cons**: More complex
- **When**: Diverse user base with different expertise levels

---

## 🔍 Key Technical Questions

### Q1: File Management (CRITICAL!)
**Does restart_folder automatically include ML_AB/ML_FFN?**
- If YES → Option A (easy!)
- If NO → Option B/C (need explicit file handling)

### Q2: ENCUT Reduction
**Can we reduce ENCUT in production phase (ML_ISTART=2)?**
- Test script uses `ENCUT=200` in Stage 1 (vs 400 in Stage 0)
- If works → Document it
- If fails → Keep ENCUT constant

### Q3: Multi-Structure Strategy
**Should each slab train its own ML model?**
- Independent models: Already supported (current scatter pattern)
- Shared model: Requires merging ML_AB files (complex)

---

## 📚 Reading Order

1. **First Time**: Read this README
2. **Understand MLFF**: Read `VASP_MLFF_KNOWLEDGE.md`
3. **See Differences**: Read `MLFF_COMPARISON.md`
4. **Quick Context**: Read `MLFF_DISCUSSION.md`
5. **Design Input**: Answer `QUESTIONS_FOR_YOU.md`
6. **Deep Dive**: Read `aimd_mlff_brainstorm.md`
7. **Test It**: Run `step_XX_mlff_test.py`

---

## 🔬 Testing Strategy

### Phase 1: Feasibility Test (This Week)
```bash
# Goal: Does MLFF work with current code?
python step_XX_mlff_test.py

# Success criteria:
# ✓ Stage 0 completes with ML_ABN + ML_FFN generated
# ✓ Stage 1 finds and uses ML_FFN
# ✓ Stage 1 runs much faster than Stage 0
# ✓ Reduced ENCUT works (or document why not)
```

### Phase 2: Design Decision (After Test)
```bash
# Based on test results, choose Option A, B, or C
# Answer QUESTIONS_FOR_YOU.md
# Finalize architecture
```

### Phase 3: Implementation (TBD)
```bash
# Implement chosen option
# Create unit tests (if needed)
# Create example scripts
# Document usage
```

---

## 📋 Checklist

### Before Implementation:
- [ ] Run `step_XX_mlff_test.py`
- [ ] Check if ML files work with restart
- [ ] Read all documentation
- [ ] Answer `QUESTIONS_FOR_YOU.md`
- [ ] Decide on Option A, B, or C

### During Implementation:
- [ ] Create/modify code in chosen location
- [ ] Add validation (if needed)
- [ ] Handle ML file management (if needed)
- [ ] Write unit tests (if desired)
- [ ] Create example scripts

### After Implementation:
- [ ] Test with real systems
- [ ] Document usage
- [ ] Create tutorial
- [ ] Move to `teros/core/` (graduate from experimental!)

---

## 🤝 Next Steps

### For You:
1. **Read documentation** in this folder
2. **Answer QUESTIONS_FOR_YOU.md**
3. **Run step_XX_mlff_test.py** (optional but recommended)
4. **Share findings** (did ML files work? any errors?)

### For Us (After Your Input):
1. Choose implementation approach
2. Design module architecture
3. Implement chosen option
4. Create tests and examples
5. Document everything

---

## 📞 Questions?

If you have questions about:
- **VASP MLFF specifics**: Check `VASP_MLFF_KNOWLEDGE.md`
- **Design options**: Check `MLFF_COMPARISON.md` and `MLFF_DISCUSSION.md`
- **Implementation details**: Check `aimd_mlff_brainstorm.md`
- **How to proceed**: Check `QUESTIONS_FOR_YOU.md`

---

## 🎉 Goal

**Enable fast, accurate molecular dynamics simulations using VASP's machine learning force fields, integrated seamlessly with PS-TEROS's standalone AIMD module.**

**Target**: 10-100x more MD sampling at similar or lower computational cost!
