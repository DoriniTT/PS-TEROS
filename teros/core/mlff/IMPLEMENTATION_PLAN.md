# MLFF Implementation Plan

Step-by-step development of MLFF module for PS-TEROS.

---

## Phase 1: Basic Workflow ✅ (Current)

**Goal**: Simple training → production workflow

**Completed**:
- [x] Module structure (`teros/core/mlff/`)
- [x] `build_mlff_workgraph()` function
- [x] Two-stage workflow (ML_ISTART=0 → ML_ISTART=2)
- [x] Si bulk test script
- [x] Reuses existing AIMD infrastructure

**Test**: `python examples/vasp/step_20_mlff_si.py`

**Critical check**: Does Stage 1 find ML_FFN from Stage 0?

---

## Phase 2: Validation & Refinement

**Goals**:
- Add ML_ISTART=1 refinement stage support
- Basic validation (check ML_LMLFF, ML_ISTART values)
- Error messages for common mistakes

**Tasks**:
- [ ] Add optional refinement stage
- [ ] Validate MLFF parameters
- [ ] Check file dependencies
- [ ] Warn on invalid configs

---

## Phase 3: Advanced Features

**Goals**:
- Multi-stage temperature annealing
- ML model quality checking
- Per-structure ML parameters
- MLFF builder with defaults

**Tasks**:
- [ ] Parse ML_LOGFILE for quality metrics
- [ ] Support arbitrary stage sequences
- [ ] Add structure-specific ML settings
- [ ] Create `get_mlff_defaults()` builder

---

## Phase 4: Production Ready

**Goals**:
- Comprehensive testing
- Documentation
- Examples for common use cases

**Tasks**:
- [ ] Unit tests
- [ ] Multiple example scripts
- [ ] Best practices guide
- [ ] Performance benchmarks

---

## Current Status

**Working**: Basic two-stage MLFF workflow
**Next**: Test Phase 1, then move to Phase 2
**Timeline**: 1+ month, incremental development
