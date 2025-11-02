# Custom VASP Calculation Module Examples

This directory contains examples for the `custom_calculation` module.

## Examples

### test_single.py
Run a custom VASP calculation on a single structure.

Usage:
```bash
source ~/envs/aiida/bin/activate
python test_single.py
```

### test_multiple.py
Run custom VASP calculations on multiple structures with the same settings.

### test_different_settings.py
Run custom VASP calculations on multiple structures with different settings for each.

## Requirements

- AiiDA profile: `psteros`
- VASP code: `VASP-6.4.1@cluster02` (or modify code_label in scripts)
- AiiDA daemon running: `verdi daemon start`
