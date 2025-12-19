# PS-TEROS CI/CD Roadmap

This document outlines the continuous integration and continuous deployment strategy for PS-TEROS.

**Last Updated**: 2025-11-26

---

## Overview

PS-TEROS uses a 2-tier CI/CD strategy:

| Tier | Name | Trigger | Duration | Resources |
|------|------|---------|----------|-----------|
| 1+2 | Unit & Integration Tests | Every PR/Push | ~5-10 min | GitHub Actions (PostgreSQL + RabbitMQ) |
| 3 | End-to-End Tests | Weekly/Release | Hours | Self-hosted runner |

---

## Combined Tests (Tier 1 + Tier 2) ✅ IMPLEMENTED

**Status**: Active  
**Trigger**: Every push and pull request  
**Duration**: ~5-10 minutes  

### What It Tests

1. **Syntax & Import Validation**
   - All Python files can be imported without errors
   - No syntax errors in codebase

2. **Pure Python Function Tests**
   - `parse_formula()` - Chemical formula parsing
   - `deep_merge_dicts()` - Dictionary merging utility
   - Workflow preset validation
   - Miller index handling
   - Stoichiometry calculations

3. **AiiDA Data Type Tests**
   - Node creation (Float, Int, Str, Bool, List, Dict)
   - StructureData creation from ASE
   - Node storage and retrieval

4. **Calcfunction Logic Tests**
   - Energy extraction logic
   - Oxide type identification
   - Surface area calculations
   - Adsorption energy formulas
   - Formation enthalpy calculations

5. **Task Function Existence**
   - Verify TaskHandle objects are defined
   - Check scatter functions exist

6. **Code Style**
   - Flake8 linting (warnings only, non-blocking)

### Files

```
.github/workflows/ci.yml          # Main CI workflow
tests/
├── __init__.py
├── conftest.py                   # Pytest fixtures
├── test_imports.py               # Import validation
├── test_pure_functions.py        # Pure Python tests
├── test_workflow_presets.py      # Preset validation
├── test_structure_utils.py       # Structure utilities
├── test_aiida_integration.py     # AiiDA data types & logic
└── test_workgraph_construction.py # WorkGraph validation
```

### GitHub Actions Configuration

Uses PostgreSQL and RabbitMQ service containers:
- `postgres:14` - Database backend
- `rabbitmq:3-management` - Message broker
- Tests run on Python 3.9, 3.10, 3.11

### Running Locally

```bash
# Install test dependencies
pip install pytest pytest-cov flake8

# Ensure AiiDA is configured
verdi status

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=teros --cov-report=html

# Run linting
flake8 teros/ --max-line-length=120 --ignore=E501,W503
```

---

## Tier 3: End-to-End Tests 🔜 PLANNED

**Status**: Planned  
**Trigger**: Weekly schedule or manual (release)  
**Duration**: Hours (depends on calculations)  

### What It Will Test

1. **Real VASP Calculations**
   - Bulk relaxation on cluster02
   - Slab generation and relaxation
   - Formation enthalpy calculation

2. **Full Workflow Execution**
   - Complete `surface_thermodynamics` preset
   - Multi-stage AIMD simulations
   - Adsorption energy workflow

3. **Output Validation**
   - Compare against reference values
   - Energy differences within tolerance
   - Structure RMSD checks

### Implementation Notes

- Requires self-hosted GitHub Actions runner
- Runner must have access to cluster02
- AiiDA daemon must be running
- Reference outputs stored in `tests/references/`

### Self-Hosted Runner Setup (Future)

```bash
# On the machine with AiiDA access:
mkdir actions-runner && cd actions-runner
curl -o actions-runner-linux-x64.tar.gz -L https://github.com/actions/runner/releases/download/v2.x.x/actions-runner-linux-x64-2.x.x.tar.gz
tar xzf ./actions-runner-linux-x64.tar.gz
./config.sh --url https://github.com/YOUR_ORG/PS-TEROS --token YOUR_TOKEN
./run.sh
```

---

## Reference Values

For Tier 3 tests, we need reference outputs from known-good calculations:

| System | Property | Reference Value | Tolerance |
|--------|----------|-----------------|-----------|
| Ag2O bulk | Energy/atom | TBD | ±0.001 eV |
| Ag2O (100) | Surface energy | TBD | ±0.01 J/m² |
| Ag2O | Formation enthalpy | TBD | ±0.01 eV/f.u. |

**Action Item**: Run reference calculations and store outputs.

---

## Future Enhancements

### Automated Changelog
- Parse commit messages for changelog generation
- Use conventional commits format

### Version Bumping
- Automatic version bump on release
- Semantic versioning (major.minor.patch)

### Documentation CI
- Build Sphinx docs on PR
- Deploy to GitHub Pages on merge

### Code Coverage Tracking
- Upload coverage to Codecov
- Coverage badges in README

---

## Workflow Files Reference

### Main CI Workflow (`.github/workflows/ci.yml`)

Runs on every push/PR:
- Python 3.9, 3.10, 3.11
- Install dependencies
- Run pytest
- Run flake8

### Release Workflow (Future)

Runs on tag push:
- Build package
- Upload to PyPI
- Create GitHub release

---

## Troubleshooting

### Common CI Failures

| Error | Cause | Solution |
|-------|-------|----------|
| Import error | Missing dependency | Add to `setup.py` |
| AiiDA not configured | Missing profile | Use mock/skip decorator |
| Timeout | Long-running test | Move to Tier 2/3 |
| Flake8 error | Style violation | Fix or add to ignore list |

### Skipping Tests

```python
import pytest

@pytest.mark.skipif(not AIIDA_AVAILABLE, reason="AiiDA not configured")
def test_requires_aiida():
    ...

@pytest.mark.slow
def test_long_running():
    ...
```

Run without slow tests:
```bash
pytest -m "not slow"
```

---

## Contributing

When adding new features:

1. Add unit tests in `tests/test_*.py`
2. Ensure all Tier 1 tests pass locally
3. Update this roadmap if adding new test categories
