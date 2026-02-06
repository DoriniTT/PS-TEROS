# PS-TEROS Testing Guide

This document provides comprehensive information about the PS-TEROS testing framework, which uses a three-tier strategy for efficient development and validation.

## Table of Contents

1. [Overview](#overview)
2. [Tier1 - Pure Python Tests](#tier1---pure-python-tests)
3. [Tier2 - AiiDA Integration Tests](#tier2---aiida-integration-tests)
4. [Tier3 - End-to-End Tests](#tier3---end-to-end-tests)
5. [Quick Start](#quick-start)
6. [Testing Workflows](#testing-workflows)
7. [Missing Tests for Remaining Bricks](#missing-tests-for-remaining-bricks)
8. [Troubleshooting](#troubleshooting)

## Overview

The testing framework consists of **three tiers**, each optimizing for a specific combination of speed and validation depth:

| Tier | Purpose | Duration | Needs VASP | Needs AiiDA | Use Case |
|------|---------|----------|-----------|------------|----------|
| **Tier1** | Pure Python logic validation | seconds | ❌ | ❌ | Parsing, validation, algorithms |
| **Tier2** | AiiDA WorkGraph construction | 5-10 min | ❌ | ✅ | Calcfunctions, task creation |
| **Tier3** | Result extraction validation | seconds | ✅ | ✅ | Data processing, schemas |

### Tier Testing Coverage (Updated 2026-02-06)

- **~120 tier1 tests** passing (pure Python validation)
- **45 tier2 tests** passing (lego brick integration)
- **30 tier3 tests** passing (pre-computed VASP scenarios)
- **6 reference calculations** generated and stored in lego_reference_pks.json
- All tests properly isolated with pytest markers

### Brick Test Coverage Status

| Brick | Tier1 | Tier2 | Tier3 | Reference PKs | Status |
|-------|-------|-------|-------|---------------|--------|
| vasp | ✅ | ✅ | ✅ (8 tests) | 2 (relax, scf) | Complete |
| dos | ✅ | ✅ | ✅ (5 tests) | 1 (scf+dos) | Complete |
| batch | ✅ | ✅ | ✅ (5 tests) | 1 (encut scan) | Complete |
| aimd | ✅ | ✅ | ✅ (5 tests) | 1 (short MD) | Complete |
| sequential | ✅ | ✅ | ✅ (7 tests) | 1 (relax→scf) | Complete |
| bader | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| convergence | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| thickness | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| hubbard_response | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| hubbard_analysis | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| qe | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| cp2k | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| generate_neb_images | ✅ | ❌ | ❌ | 0 | **Needs tests** |
| neb | ✅ | ❌ | ❌ | 0 | **Needs tests** |

## Tier1 - Pure Python Tests

### What It Tests

Pure Python logic without any external dependencies:
- File parsing (VASP, ASE formats)
- Data structure validation
- Parameter merging/transformation
- Connection validation rules
- Utility functions

### What It Doesn't Need

- AiiDA profile
- VASP installation
- Any cluster connectivity
- Any compute resources

### Key Test Files

| File | Coverage |
|------|----------|
| `tests/test_lego_connections.py` | Port declarations, connection validation |
| `tests/test_lego_bricks.py` | Brick configuration validation |
| `tests/test_lego_concurrent_and_outputs.py` | Output naming, sequential workflows |
| `tests/test_pure_functions.py` | General utility functions |

### Running Tier1 Tests

```bash
# All tier1 tests
pytest tests/ -m tier1 -v

# Specific test file
pytest tests/test_lego_connections.py -m tier1 -v

# Specific test class
pytest tests/test_lego_connections.py::TestPortTypeRegistry -m tier1 -v

# Specific test
pytest tests/test_lego_connections.py::TestPortTypeRegistry::test_all_outputs_recognized -m tier1 -v
```

### Example Tier1 Test

```python
@pytest.mark.tier1
class TestPortTypeRegistry:
    """Test that all port types are registered and recognized."""
    
    def test_all_output_types_recognized(self):
        from teros.core.lego.bricks.connections import PORT_TYPES, ALL_PORTS
        
        for brick_name, ports in ALL_PORTS.items():
            for output_name, output_spec in ports['outputs'].items():
                assert output_spec['type'] in PORT_TYPES
```

## Tier2 - AiiDA Integration Tests

### What It Tests

WorkGraph construction and calcfunction behavior with real AiiDA nodes:
- Calcfunction execution with mock data
- WorkGraph task creation and wiring
- Output namespace management
- Real AiiDA node handling (Dict, StructureData, TrajectoryData, etc.)
- Result extraction logic

### Requirements

- AiiDA profile configured and active
- Python packages installed (no VASP code needed)
- ~5-10 minutes runtime

### Key Test Files

| File | Coverage |
|------|----------|
| `tests/test_lego_vasp_integration.py` | VASP brick calcfunctions, energy extraction |
| `tests/test_lego_dos_integration.py` | DOS validation, SCF+DOS workflow |
| `tests/test_lego_batch_integration.py` | Batch stage validation and creation |
| `tests/test_lego_aimd_integration.py` | AIMD calcfunctions, trajectory handling |
| `tests/test_lego_sequential_integration.py` | Multi-stage pipeline validation |
| `tests/test_lego_aimd_trajectory_concatenation.py` | Trajectory concatenation logic |

### Running Tier2 Tests

```bash
# All tier2 tests
pytest tests/ -m tier2 -v

# Lego module tier2 tests
pytest tests/test_lego_*_integration.py -m tier2 -v

# Specific test class
pytest tests/test_lego_vasp_integration.py::TestVaspExtractEnergy -m tier2 -v
```

### Setup Required

```bash
# Activate environment and AiiDA
source ~/envs/aiida/bin/activate
verdi profile set-default presto
verdi daemon restart  # Ensure daemon is running
```

### Example Tier2 Test

```python
@pytest.mark.tier2
@pytest.mark.requires_aiida
class TestVaspExtractEnergy:
    """Test extract_energy calcfunction with mock misc Dict."""
    
    def test_extract_energy_from_total_energies(self):
        from aiida import orm
        from aiida_workgraph import WorkGraph
        from teros.core.lego.tasks import extract_energy
        
        misc = orm.Dict(dict={
            'total_energies': {
                'energy_extrapolated': -10.12345,
            }
        })
        
        wg = WorkGraph(name='test_extract_energy')
        wg.add_task(extract_energy, name='extract', misc=misc)
        wg.run()
        
        assert wg.tasks['extract'].state == 'FINISHED'
        result = wg.tasks['extract'].outputs.result.value
        assert abs(float(result) - (-10.12345)) < 1e-8
```

## Tier3 - End-to-End Tests

### What It Tests

Result extraction and validation against pre-computed VASP calculations:
- Energy extraction accuracy
- Structure output presence/absence
- Trajectory dimensions and coordinate systems
- Misc dict schemas
- Result consistency (relax vs. SCF energy agreement)
- File retrieval

### Requirements

- AiiDA profile with reference PKs loaded
- Pre-computed VASP calculations stored in `tests/fixtures/lego_reference_pks.json`
- Each calculation runs once; results reused across tests

### Reference PKs

Stored in `tests/fixtures/lego_reference_pks.json`:

```json
{
    "vasp": {
        "relax_si": {
            "pk": 12345,
            "code_label": "VASP-6.5.1@localwork",
            "description": "Si diamond relaxation"
        },
        "scf_sno2": {
            "pk": 12346,
            "code_label": "VASP-6.5.1@localwork",
            "description": "SnO2 static SCF"
        }
    },
    "dos": {
        "dos_sno2": {
            "pk": 12347,
            "code_label": "VASP-6.5.1@localwork"
        }
    },
    "batch": {
        "batch_si_encut": {
            "pk": 12348,
            "code_label": "VASP-6.5.1@localwork",
            "stage_names": ["encut_scan"],
            "stage_types": {"encut_scan": "batch"},
            "stage_namespaces": {…}
        }
    },
    "aimd": {
        "aimd_si": {
            "pk": 12349,
            "code_label": "VASP-6.5.1@localwork",
            "stage_names": ["short_md"],
            "stage_types": {"short_md": "aimd"},
            "stage_namespaces": {…}
        }
    },
    "sequential": {
        "relax_then_scf_si": {
            "pk": 12350,
            "code_label": "VASP-6.5.1@localwork",
            "stage_names": ["relax", "scf"],
            "stage_types": {"relax": "vasp", "scf": "vasp"},
            "stage_namespaces": {…}
        }
    }
}
```

### Key Test Files

| File | Scenarios Tested |
|------|------------------|
| `tests/test_lego_vasp_integration.py` | Si relax, SnO2 SCF |
| `tests/test_lego_dos_integration.py` | SnO2 DOS |
| `tests/test_lego_batch_integration.py` | Si ENCUT scan (3 values) |
| `tests/test_lego_aimd_integration.py` | Si AIMD (5 steps) |
| `tests/test_lego_sequential_integration.py` | Si relax → SCF pipeline |

### Running Tier3 Tests

```bash
# Check reference PK status
python tests/generate_lego_references.py --status

# All tier3 tests (skip if PKs missing)
pytest tests/ -m tier3 -v

# Lego module tier3 tests
pytest tests/test_lego_*_integration.py -m tier3 -v

# Specific test class
pytest tests/test_lego_vasp_integration.py::TestVaspRelaxResultExtraction -m tier3 -v
```

### Example Tier3 Test

```python
@pytest.mark.tier3
@pytest.mark.requires_aiida
class TestVaspRelaxResultExtraction:
    """Validate result extraction from completed VASP relaxation."""
    
    def test_get_results_returns_energy(self, reference_pks):
        from tests.conftest import load_node_or_skip
        from teros.core.lego.results import get_results
        
        pk = reference_pks['vasp']['relax_si']['pk']
        load_node_or_skip(pk)  # Skip if node missing or PK null
        
        result = get_results(pk)
        assert result['energy'] is not None
        assert result['energy'] < 0  # Si energy should be negative
```

## Generating Tier3 References

### Initial Setup

The reference generator creates lightweight VASP calculations for all brick types:

```bash
# Generate all references
python tests/generate_lego_references.py --code VASP-6.5.1@localwork

# Generate specific brick only
python tests/generate_lego_references.py --code VASP-6.5.1@localwork --brick vasp
```

### Available Bricks for Reference Generation

- `vasp` - VASP brick (Si relax + SnO2 SCF)
- `dos` - DOS brick (SnO2 DOS)
- `batch` - Batch brick (Si ENCUT scan)
- `aimd` - AIMD brick (Si MD)
- `sequential` - Sequential brick (Si relax → SCF)

### Monitoring Progress

```bash
# Check status of all submitted calculations
python tests/generate_lego_references.py --status

# Wait for all to finish (polls every 15s)
python tests/generate_lego_references.py --wait

# Wait with custom poll interval
python tests/generate_lego_references.py --wait --poll 30
```

### Reference Specifications

Each reference is designed to be lightweight:

| Scenario | Parameters | Purpose |
|----------|-----------|---------|
| Si relax | ENCUT=300, NSW=20, ISIF=3 | Test structure output + energy |
| SnO2 SCF | ENCUT=300, NSW=0 | Test no structure (NSW=0) + energy |
| SnO2 DOS | ENCUT=300, ISMEAR=-5 | Test dual-stage DOS workflow |
| Si ENCUT batch | 3 values: 200, 250, 300 | Test parameter variation |
| Si AIMD | TEBEG=300K, NSW=5 | Test trajectory output |
| Si relax→SCF | 2 stages: relax + static | Test sequential workflow |

## Quick Start

### One-Time Setup

```bash
# 1. Activate environment
source ~/envs/aiida/bin/activate
verdi profile set-default presto

# 2. Start AiiDA daemon
verdi daemon restart

# 3. Generate tier3 references (one time)
python tests/generate_lego_references.py --code VASP-6.5.1@localwork
python tests/generate_lego_references.py --wait
```

### Daily Testing Workflow

```bash
# After code changes, always restart daemon
verdi daemon restart

# 1. Fast tier1 validation (seconds)
pytest tests/test_lego_connections.py -m tier1 -v

# 2. Deeper tier2 integration (minutes)
pytest tests/test_lego_vasp_integration.py -m tier2 -v

# 3. Final tier3 validation (seconds, reuses pre-computed data)
pytest tests/test_lego_vasp_integration.py -m tier3 -v
```

## Testing Workflows

### Scenario: Fixing a Tier2 Test Failure

```bash
# 1. Run failing test with verbose output
pytest tests/test_lego_vasp_integration.py::TestVaspExtractEnergy::test_extract_energy_from_total_energies -m tier2 -vv

# 2. Edit code and restart daemon
verdi daemon restart

# 3. Re-run test
pytest tests/test_lego_vasp_integration.py::TestVaspExtractEnergy::test_extract_energy_from_total_energies -m tier2 -vv

# 4. If fixed, run full tile2 suite
pytest tests/test_lego_*_integration.py -m tier2 -v
```

### Scenario: Adding a New Lego Brick

```bash
# Add unit tests for configuration validation
# File: tests/test_lego_bricks.py
# Add: @pytest.mark.tier1 class TestMyBrickValidateStage

# Add integration tests for calcfunctions
# File: tests/test_lego_<mybrick>_integration.py
# Add: @pytest.mark.tier2 class TestMyBrickFunctionality

# Test tier1
pytest tests/test_lego_bricks.py::TestMyBrickValidateStage -m tier1 -v

# Test tier2
pytest tests/test_lego_<mybrick>_integration.py -m tier2 -v

# Optional: Add tier3 tests after generating references
```

### Scenario: Debugging a Result Extraction Issue

```bash
# 1. Find the pre-computed PK
python tests/generate_lego_references.py --status

# 2. Load and inspect manually
python
>>> from aiida import orm, load_profile
>>> load_profile()
>>> wg = orm.load_node(12345)
>>> wg.outputs
# Browse outputs

# 3. Run the affected tier3 test with debugging
pytest tests/test_lego_vasp_integration.py::TestVaspRelaxResultExtraction::test_get_results_returns_energy -m tier3 -vv

# 4. Fix extraction logic in teros/core/lego/results.py
# 5. Restart daemon and re-test
verdi daemon restart
pytest tests/test_lego_vasp_integration.py::TestVaspRelaxResultExtraction -m tier3 -v
```

## Missing Tests for Remaining Bricks

The following bricks need tier2 and tier3 tests implemented. This section provides a roadmap for adding comprehensive test coverage.

### Bader Brick

**What to test:**

**Tier2:**
- `test_bader_validation.py`:
  - Validate `charge_from` parameter (must reference VASP stage with `laechg=True`)
  - Validate INCAR requirements are present
  - Test error handling when CHGCAR/AECCAR files missing

**Tier3:**
- `test_lego_bader_integration.py`:
  - Test `get_stage_results()` extracts bader charges correctly
  - Validate output schema: `{'charges': {...}, 'pk': int, 'stage': str, 'type': 'bader'}`
  - Test per-atom charge assignment matches expected format

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - VASP relax with laechg=True
# - Bader stage using charge_from='relax'
```

### Convergence Brick

**What to test:**

**Tier2:**
- `test_convergence_validation.py`:
  - Validate `convergence_type` parameter ('encut' or 'kpoints')
  - Validate `values` parameter is a list with ≥2 elements
  - Test error handling when `structure_from` invalid

**Tier3:**
- `test_lego_convergence_integration.py`:
  - Test `get_stage_results()` extracts convergence data
  - Validate output schema: `{'convergence_results': {...}, 'energies': [...]}`
  - Test energy array length matches values array length
  - Test convergence threshold detection

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - VASP relax for bulk structure
# - Convergence stage: ENCUT scan [200, 250, 300]
```

### Thickness Brick

**What to test:**

**Tier2:**
- `test_thickness_validation.py`:
  - Validate `miller_indices` parameter (list of 3 ints)
  - Validate `layer_counts` parameter (list of ≥2 ints)
  - Test error handling for missing `structure_from`/`energy_from` or `bulk_incar`
  - Test standalone vs. chained input modes

**Tier3:**
- `test_lego_thickness_integration.py`:
  - Test `get_stage_results()` extracts thickness convergence data
  - Validate output schema: `{'convergence_results': {...}, 'surface_energies': {...}}`
  - Test per-thickness energy extraction
  - Test convergence threshold detection

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - VASP bulk relax (e.g., SnO2)
# - Thickness stage: [3, 5] layers on (110) surface
```

### Hubbard Response/Analysis Bricks

**What to test:**

**Tier2:**
- `test_hubbard_validation.py`:
  - **Response brick:**
    - Validate `ground_state_from` references a VASP stage with `lorbit=11`, `lwave=True`, `ldau=False`
    - Validate `target_species` is in structure
    - Validate `potential_values` list excludes 0.0
    - Validate `ldaul` is 2 or 3
  - **Analysis brick:**
    - Validate `response_from` references a hubbard_response stage
    - Validate `target_species` matches response stage

**Tier3:**
- `test_lego_hubbard_integration.py`:
  - **Response tests:**
    - Test `get_stage_results()` extracts response data for each potential
    - Validate output schema: `{'responses': {...}, 'occupations': {...}}`
    - Test occupations array matches number of potentials
  - **Analysis tests:**
    - Test `get_stage_results()` extracts U value from linear regression
    - Validate output schema: `{'U_value': float, 'chi_slope': float, 'chi0': float}`
    - Test U value is positive and physically reasonable (0-10 eV range)

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - VASP relax for SnO2
# - Ground state SCF with ldau=False, lorbit=11, lwave=True, lcharg=True
# - Hubbard response: target_species='Sn', potentials=[-0.1, 0.1]
# - Hubbard analysis: from response stage
```

### QE Brick

**What to test:**

**Tier2:**
- `test_qe_validation.py`:
  - Validate `calculation` parameter ('scf', 'relax', 'vc-relax')
  - Validate `parameters` dict structure (CONTROL, SYSTEM, ELECTRONS sections)
  - Validate pseudopotential family exists
  - Test error handling for invalid structure

**Tier3:**
- `test_lego_qe_integration.py`:
  - Test `get_stage_results()` extracts energy correctly
  - Test relaxed structure extraction (for relax/vc-relax)
  - Validate output schema matches VASP brick pattern
  - Test energy is negative for converged calculations

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - QE SCF for Si primitive cell
# - QE relax for Si primitive cell
# Use lightweight parameters: ecutwfc=30, k=0.15
```

### CP2K Brick

**What to test:**

**Tier2:**
- `test_cp2k_validation.py`:
  - Validate `run_type` parameter ('ENERGY', 'GEO_OPT', 'CELL_OPT', 'MD')
  - Validate `parameters` dict structure
  - Validate basis set and pseudopotential availability
  - Test error handling for invalid structure

**Tier3:**
- `test_lego_cp2k_integration.py`:
  - Test `get_stage_results()` extracts energy correctly
  - Test relaxed structure extraction (for GEO_OPT/CELL_OPT)
  - Test trajectory extraction (for MD)
  - Validate output schema

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - CP2K ENERGY for H2O molecule
# - CP2K GEO_OPT for H2O molecule
# Use lightweight parameters: cutoff=300 Ry
```

### Generate NEB Images Brick

**What to test:**

**Tier2:**
- `test_neb_images_validation.py`:
  - Validate `initial_from` and `final_from` parameters
  - Validate `n_images` is positive integer
  - Validate `method` is 'idpp' or 'linear'
  - Test error handling for mismatched initial/final structures

**Tier3:**
- `test_lego_neb_images_integration.py`:
  - Test image generation from two relaxed endpoints
  - Validate output: dict with n_images StructureData entries
  - Test IDPP vs linear interpolation generates different paths
  - Test atom count preserved across all images

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - VASP relax for initial structure (e.g., O in FCC hollow site)
# - VASP relax for final structure (e.g., O in FCC bridge site)
# - Generate_neb_images: n_images=3, method='idpp'
```

### NEB Brick

**What to test:**

**Tier2:**
- `test_neb_validation.py`:
  - Validate `images_from` parameter references valid stage
  - Validate NEB-specific INCAR parameters (IMAGES, SPRING, IBRION=3)
  - Test error handling when images count mismatch

**Tier3:**
- `test_lego_neb_integration.py`:
  - Test `get_stage_results()` extracts barrier height
  - Test extraction of per-image energies
  - Validate output schema: `{'barrier': float, 'energies': [...], 'forces': [...]}`
  - Test convergence detection (max force < threshold)

**Reference generation:**
```bash
# Add to generate_lego_references.py:
# - Sequential workflow:
#   1. Relax initial structure
#   2. Relax final structure
#   3. Generate 3 intermediate images
#   4. NEB with images_from='generate_images'
# Use lightweight: ENCUT=200, NSW=10
```

## Implementation Checklist for Adding Brick Tests

When adding tests for a new brick, follow this procedure:

### 1. Create Tier2 Test File

```bash
# Create test file
touch tests/test_lego_<brick>_integration.py
```

Structure:
```python
import pytest
from aiida import orm


@pytest.mark.tier2
@pytest.mark.requires_aiida
class Test<Brick>Validation:
    """Test brick-specific configuration validation."""
    
    def test_stage_validation_succeeds(self):
        """Valid configuration should pass validate_stage()."""
        from teros.core.lego.bricks.<brick> import validate_stage
        
        stage = {
            'name': 'test_stage',
            'type': '<brick>',
            # ... valid config
        }
        stage_names = {'previous_stage'}
        
        # Should not raise
        validate_stage(stage, stage_names)
    
    def test_stage_validation_rejects_invalid(self):
        """Invalid configuration should raise ValueError."""
        from teros.core.lego.bricks.<brick> import validate_stage
        
        stage = {
            'name': 'test_stage',
            'type': '<brick>',
            # ... invalid config (e.g., missing required field)
        }
        stage_names = set()
        
        with pytest.raises(ValueError, match='expected error message'):
            validate_stage(stage, stage_names)


@pytest.mark.tier2
@pytest.mark.requires_aiida
class Test<Brick>CalcfunctionExecution:
    """Test calcfunction behavior with mock data."""
    
    def test_result_extraction_from_mock_data(self):
        """Test extraction logic with synthetic inputs."""
        # Create mock AiiDA nodes
        mock_output = orm.Dict(dict={'key': 'value'})
        
        # Test extraction function
        from teros.core.lego.bricks.<brick> import get_stage_results
        # ... test logic
```

### 2. Add Tier3 Tests to Same File

```python
@pytest.mark.tier3
@pytest.mark.requires_aiida
@pytest.mark.localwork
class Test<Brick>ResultExtraction:
    """Test result extraction from pre-computed calculations."""
    
    @pytest.fixture
    def reference_pks(self):
        """Load reference PKs from fixtures."""
        import json
        from pathlib import Path
        
        pks_file = Path(__file__).parent / 'fixtures' / 'lego_reference_pks.json'
        with open(pks_file) as f:
            return json.load(f)
    
    def test_returns_expected_output_type(self, reference_pks):
        """Brick should return correct output type."""
        from tests.conftest import build_sequential_result
        from teros.core.lego.results import get_stage_results
        
        scenario = reference_pks['<brick>']['<scenario_name>']
        seq_result = build_sequential_result(scenario)
        
        stage_name = seq_result['__stage_names__'][0]
        result = get_stage_results(seq_result, stage_name)
        
        assert result is not None
        assert 'pk' in result
        assert 'type' in result
        assert result['type'] == '<brick>'
        # ... type-specific assertions
```

### 3. Update Reference Generator

Edit `tests/generate_lego_references.py`:

```python
def generate_<brick>_references(code_label: str):
    """Generate <brick> reference calculations."""
    print('[<brick>]')
    
    # Load or create test structure
    structure = ...
    
    # Build stages
    stages = [
        {
            'name': '<stage_name>',
            'type': '<brick>',
            # ... configuration
        }
    ]
    
    # Submit workflow
    from teros.core.lego import quick_vasp_sequential
    result = quick_vasp_sequential(
        structure=structure,
        code_label=code_label,
        stages=stages,
        # ...
    )
    
    pk = result['__workgraph_pk__']
    print(f'  Submitted <brick> <scenario>: PK={pk}')
    
    # Save to JSON
    update_reference_pks('<brick>', '<scenario_name>', {
        'pk': pk,
        'code_label': code_label,
        'description': '...',
        'stage_names': [...],
        'stage_types': {...},
        'stage_namespaces': {...},
    })
```

Add to main brick list:
```python
BRICKS = ['vasp', 'dos', 'batch', 'aimd', 'sequential', '<brick>']
```

### 4. Generate and Verify References

```bash
# Generate references
python tests/generate_lego_references.py --code VASP-6.5.1@localwork --brick <brick>

# Wait for completion
python tests/generate_lego_references.py --wait

# Verify status
python tests/generate_lego_references.py --status
```

### 5. Run Tests

```bash
# Run tier2
pytest tests/test_lego_<brick>_integration.py -m tier2 -v

# Run tier3 (after references generated)
pytest tests/test_lego_<brick>_integration.py -m tier3 -v
```

### 6. Update Documentation

- Update AGENTS.md testing section brick coverage table
- Update this guide's brick coverage status table
- Add brick-specific testing notes if needed

## Troubleshooting

### Issue: Tier2 tests fail with "AiiDA not configured" skip

**Solution:**
```bash
source ~/envs/aiida/bin/activate
verdi profile set-default presto
verdi daemon restart
pytest tests/ -m tier2 -v
```

### Issue: Tier3 tests skip with "Reference PKs not found"

**Solution:**
```bash
# Generate references (one time setup)
python tests/generate_lego_references.py --code VASP-6.5.1@localwork
python tests/generate_lego_references.py --wait

# Then run tier3
pytest tests/ -m tier3 -v
```

### Issue: Code changes don't take effect in tests

**Solution:**
```bash
# Always restart daemon after code changes
verdi daemon restart

# Verify daemon is running
verdi daemon status

# Then re-run tests
pytest tests/ -m tier2 -v
```

### Issue: Some tests marked "slow" taking too long

**Solution:**
```bash
# Exclude slow tests during development
pytest tests/ -m "not slow" -v

# Run slow tests only when necessary
pytest tests/ -m slow -v
```

### Issue: Tier2 test failure: "ImportError: cannot import name ..."

**Solution:**
This typically means code changes weren't daemon-reloaded:
```bash
# Restart daemon and try again
verdi daemon restart
pytest tests/test_lego_*_integration.py -m tier2 -v
```

### Issue: Tier3 test "Node PK=12345 not found in database"

**Solution:**
The reference PK doesn't exist in your AiiDA database. Either:
1. Generate references: `python tests/generate_lego_references.py --code VASP-6.5.1@localwork`
2. Or load the PK from a dump: contact repository maintainers

## Test Count Summary (Updated 2026-02-06)

**Completed Bricks: 5 of 14**

- Tier1: **~120 tests** across all files (fast, all passing)
- Tier2: **45 tests** across 5 brick modules (all passing)
- Tier3: **30 tests** across 5 brick modules (all passing)
- **Reference PKs: 6** (vasp: 2, dos: 1, batch: 1, aimd: 1, sequential: 1)

**Remaining Work:** 9 bricks need tier2/tier3 tests (bader, convergence, thickness, hubbard_response, hubbard_analysis, qe, cp2k, generate_neb_images, neb)

## CI/CD Integration

The CI pipeline runs:

```bash
# Fast tier1 (always runs) - ~seconds
pytest tests/ -m tier1 -v

# Tier2 (always runs) - ~5-10 min
pytest tests/ -m tier2 -v

# Tier3 (skipped in CI) - requires artifact storage
# pytest tests/ -m tier3 -v  # Not in CI by default
```

## Contributing Tests

When adding new tests:

1. **Tier1** (pure Python) - Always add if testing logic, parsing, validation
2. **Tier2** (AiiDA integration) - Add when testing WorkGraph construction or calcfunctions
3. **Tier3** (end-to-end) - Add after generating reference PKs, for result validation

Mark tests appropriately:
```python
@pytest.mark.tier1  # Pure Python
@pytest.mark.requires_aiida  # Needs AiiDA
@pytest.mark.tier2  # AiiDA integration
@pytest.mark.tier3  # End-to-end
@pytest.mark.slow  # Long-running
```

## References

- [AGENTS.md](../AGENTS.md) - Developer quick reference
- [tests/conftest.py](../tests/conftest.py) - Pytest configuration and fixtures
- [tests/generate_lego_references.py](../tests/generate_lego_references.py) - Reference generator CLI
- [tests/fixtures/lego_reference_pks.json](../tests/fixtures/lego_reference_pks.json) - Reference storage
