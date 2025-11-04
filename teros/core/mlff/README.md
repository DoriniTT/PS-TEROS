# MLFF Module

Machine Learning Force Field workflows for VASP.

Based on VASP Liquid Si MLFF example.

## Quick Start

```python
from teros.core.mlff import build_mlff_workgraph
from aiida import orm, load_profile
from ase.io import read

load_profile('presto')
structure = orm.StructureData(ase=read('Si.cif'))

# Training only (simplest)
wg = build_mlff_workgraph(
    structures={'si': structure},
    training_steps=100,
    temperature=300,
    code_label='VASP6.5.0@cluster02',
    potential_family='PBE',
    potential_mapping={'Si': 'Si'},
)

# Training + Production
wg = build_mlff_workgraph(
    structures={'si': structure},
    training_steps=100,
    production_steps=1000,  # ML-driven MD
    temperature=300,
    code_label='VASP6.5.0@cluster02',
    potential_family='PBE',
    potential_mapping={'Si': 'Si'},
)

wg.submit(wait=False)
```

## Workflow

**Training (ML_ISTART=0)**:
- DFT+MD with on-the-fly ML learning
- VASP trains neural network automatically
- Generates: ML_AB (data), ML_FFN (model), ML_LOGFILE (log)

**Production (ML_ISTART=2, optional)**:
- Pure ML predictions (no DFT!)
- 10-100x faster than DFT
- Uses ML_FFN from training

## Key Features

- **Simple API**: Just specify steps and temperature
- **Good defaults**: INCAR settings from VASP examples
- **Optional production**: Add production_steps if needed
- **Reuses AIMD**: Built on existing PS-TEROS infrastructure

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| structures | Yes | - | Dict of structures |
| training_steps | Yes | - | NSW for training |
| temperature | Yes | - | Temperature in K |
| code_label | Yes | - | VASP code |
| potential_family | Yes | - | Potential family |
| potential_mapping | Yes | - | Element mapping |
| production_steps | No | None | NSW for production |
| kpoints_spacing | No | 0.5 | K-points spacing |
| encut | No | 400 | Energy cutoff (eV) |
| prec | No | 'Normal' | VASP precision |
| options | No | cluster02 | Scheduler options |

## Examples

Training only (start here):
```bash
python examples/vasp/step_20_mlff_si.py
```

## Files

- `workgraph.py` - Main function
- `README.md` - This file
- `IMPLEMENTATION_PLAN.md` - Development roadmap

## See Also

- VASP Wiki: Liquid Si - MLFF
- Planning: `teros/experimental/mlff_module/`
