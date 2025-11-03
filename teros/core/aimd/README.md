# AIMD Standalone Module

Run multi-stage AIMD calculations on pre-existing structures with full parameter control.

## Features

- Accept StructureData nodes or PKs directly
- Sequential AIMD stages with automatic restart chaining
- Optional supercell transformation per structure
- Full builder parameter control: global, per-structure, per-stage, and matrix overrides
- Concurrency control with max_concurrent_jobs

## Quick Start

```python
from aiida import orm, load_profile
from teros.core.aimd import build_aimd_workgraph

load_profile('presto')

wg = build_aimd_workgraph(
    structures={'slab1': structure1, 'slab2': structure2},
    aimd_stages=[
        {'temperature': 300, 'steps': 100},   # Equilibration
        {'temperature': 300, 'steps': 500},   # Production
    ],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={
        'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400, ...}},
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    },
    max_concurrent_jobs=4,
)

wg.submit()
```

See `examples/vasp/step_XX_aimd_standalone.py` for complete examples.
