# Bug Fixes: Adsorption Energy Workflow

This document summarizes **two critical bugs** fixed in the adsorption energy test scripts.

## Bug 1: Missing `convergence_mode` Parameter

### Issue

When running the adsorption energy test scripts, Phase 1 (relaxation) failed with:

```
AttributeError: 'AttributeDict' object has no attribute 'convergence_mode'. Did you mean: 'convergence_on'?
```

This happened in the `VaspRelaxWorkChain.analyze_convergence` method at line 510 of the aiida-vasp plugin.

### Root Cause

The `relax_settings` dictionary was missing the **required** `convergence_mode` parameter.

According to the aiida-vasp documentation (see `/home/thiagotd/.claude/skills/aiida-builder/working-with-aiida-vasp-workchains.md`), the `convergence_mode` parameter is required when `convergence_on` is `True`.

### Valid Values

- `'last'`: Check only the change in the last ionic step (recommended for most cases)
- `'inout'`: Check difference between input and output structures

### Fix Applied

Added `'convergence_mode': 'last'` to `relax_settings` in both scripts:

1. **`run_quick_test.py`** - Line 117
2. **`run_lamno3_oh_adsorption.py`** - Line 129

## Bug 2: Missing `kpoints_spacing` in SCF Builder Inputs

### Issue

After fixing Bug 1, Phase 3 (SCF calculations) failed with:

```
InputValidationError: Must supply either 'kpoints' or 'kpoints_spacing'
```

This happened in the `VaspWorkChain.init_inputs` method at line 567 of the aiida-vasp plugin.

### Root Cause

The `scf_builder_inputs` dictionary was missing the **required** `kpoints_spacing` parameter.

When using the builder inputs approach, **each builder inputs dict must be self-contained** with all required parameters including kpoints configuration. The top-level `kpoints_spacing` parameter is NOT automatically propagated to builder inputs.

### Fix Applied

Added `'kpoints_spacing': kpoints_spacing` to `scf_builder_inputs` in both scripts:

1. **`run_quick_test.py`** - Line 164
2. **`run_lamno3_oh_adsorption.py`** - Line 189

### Reference

See the documentation in `teros/core/workgraph.py` lines 808-914 for the correct format:

```python
adsorption_scf_builder_inputs={
    'parameters': {'incar': {...}},
    'options': {...},
    'potential_family': 'PBE',
    'potential_mapping': {...},
    'kpoints_spacing': 0.15,  # <-- REQUIRED!
}
```

## Complete Templates

### Required relax_settings Template

```python
relax_builder_inputs = {
    'relax_settings': {
        # Basic relaxation control
        'algo': 'cg',                           # 'cg', 'diis', or 'fire'
        'force_cutoff': 0.1,                    # eV/Å
        'steps': 500,                            # Maximum ionic steps

        # Degrees of freedom
        'positions': True,                       # Relax atomic positions
        'shape': False,                          # Relax cell shape
        'volume': False,                         # Relax cell volume

        # Convergence checking (workchain level)
        'convergence_on': True,                  # Enable convergence checks
        'convergence_absolute': False,           # Use relative convergence
        'convergence_max_iterations': 5,         # Max restart iterations
        'convergence_positions': 0.01,           # Angstrom
        'convergence_volume': 0.01,              # Angstrom^3
        'convergence_shape_lengths': 0.1,        # Angstrom
        'convergence_shape_angles': 0.1,         # Degrees
        'convergence_mode': 'last',              # *** REQUIRED *** 'last' or 'inout'

        # Control
        'perform': True,                         # Set False to skip relaxation
    },
    'vasp': {
        'parameters': {'incar': {...}},
        'options': {...},
        'potential_family': 'PBE',
        'potential_mapping': {...},
        'kpoints_spacing': 0.6,                  # REQUIRED for VaspRelaxWorkChain
        'clean_workdir': False,
    }
}
```

### Required scf_builder_inputs Template

```python
scf_builder_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ALGO': 'Normal',
        'ENCUT': 400,
        'EDIFF': 1e-5,
        # ... other INCAR parameters ...
    }},
    'options': {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 48},
        'queue_name': 'parexp',
    },
    'potential_family': 'PBE',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O', 'H': 'H'},
    'kpoints_spacing': 0.6,  # *** REQUIRED *** for VaspWorkChain
}
```

## Testing

After both fixes, run:

```bash
# Clear Python cache (recommended)
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

# Restart daemon
verdi daemon restart

# Verify syntax
python3 -m py_compile run_quick_test.py

# Run the test
source ~/envs/aiida/bin/activate
python run_quick_test.py

# Monitor
verdi process list -a -p1
verdi process show <PK>
```

## Status

✓ Bug 1 fixed: `convergence_mode` added to `relax_settings`
✓ Bug 2 fixed: `kpoints_spacing` added to `scf_builder_inputs`
✓ Both scripts syntax verified
✓ Ready for full workflow testing

The quick test should now run successfully through all 4 phases:
1. ✓ Phase 1: Relaxation (with convergence_mode)
2. ✓ Phase 2: Separation (automatic)
3. ✓ Phase 3: SCF calculations (with kpoints_spacing)
4. ✓ Phase 4: Adsorption energy calculation

Expected runtime: ~10-20 minutes for quick test, ~12-18 hours for LaMnO3.
