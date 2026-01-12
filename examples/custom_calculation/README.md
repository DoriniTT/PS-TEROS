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

### test_dos.py
Run DOS calculations using the `build_dos_calculation_workgraph` function.

This workflow:
1. Performs SCF calculation to get charge density (CHGCAR)
2. Performs non-SCF DOS calculation using the CHGCAR

Usage:
```bash
source ~/envs/aiida/bin/activate
python test_dos.py
```

To retrieve results after completion:
```python
from aiida import orm
from teros.core.custom_calculation import get_dos_results

wg = orm.load_node(<PK>)
results = get_dos_results(wg)
print(results['dos'])        # DOS data
print(results['projectors']) # Projected DOS data
```

## Requirements

- AiiDA profile: `psteros`
- VASP code: `VASP-6.4.1@cluster02` (or modify code_label in scripts)
- AiiDA daemon running: `verdi daemon start`
