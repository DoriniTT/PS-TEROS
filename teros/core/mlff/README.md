# MLFF Module

Machine Learning Force Field workflows for PS-TEROS using VASP.

## Quick Start

```python
from teros.core.mlff import build_mlff_workgraph
from aiida import orm, load_profile
from ase.io import read

load_profile('presto')

structure = orm.StructureData(ase=read('Si.cif'))

builder_inputs = {
    'parameters': {
        'incar': {
            'PREC': 'Normal',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 0,
            'MDALGO': 2,
            'POTIM': 1.0,
            'LWAVE': True,
            'LCHARG': True,
        }
    },
    'kpoints_spacing': 0.5,
    'potential_family': 'PBE',
    'potential_mapping': {'Si': 'Si'},
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    },
    'clean_workdir': False,
}

wg = build_mlff_workgraph(
    structures={'si_bulk': structure},
    training_steps=200,
    production_steps=5000,
    temperature=300,
    code_label='VASP6.5.0@cluster02',
    builder_inputs=builder_inputs,
)

wg.submit(wait=False)
```

## Workflow

```
Stage 0: Training (ML_ISTART=0)
  - DFT+ML learning for N steps
  - Generates ML_AB, ML_FFN
  ↓
Stage 1: Production (ML_ISTART=2)
  - Pure ML predictions
  - 10-100x faster than DFT
```

## Files

- `workgraph.py` - Main `build_mlff_workgraph()` function
- `__init__.py` - Module exports

## Implementation Status

**Phase 1 (Current)**:
- Basic two-stage workflow (training → production)
- Reuses existing AIMD infrastructure
- Simple API with minimal parameters

**Phase 2 (Future)**:
- Refinement stage support (ML_ISTART=1)
- Validation utilities
- ML model quality checking
- Multi-stage temperature annealing

## See Also

- Test script: `examples/vasp/step_XX_mlff_si.py`
- Planning docs: `teros/experimental/mlff_module/`
