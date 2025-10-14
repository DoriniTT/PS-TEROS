# AIMD (Ab Initio Molecular Dynamics) in PS-TEROS

## Overview

AIMD capability allows running sequential molecular dynamics simulations on slab structures as a parallel, independent analysis alongside traditional relaxation calculations.

## Features

- Sequential AIMD stages with automatic restart chaining
- Flexible temperature/timestep configuration
- Runs on all generated or input slabs
- Full provenance tracking of all stages
- Compatible with binary and ternary oxides

## Usage

```python
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_aimd_defaults

# Define AIMD parameters
aimd_parameters = get_aimd_defaults(
    energy_cutoff=400,
    timestep=1.0,  # 1 fs
)

# Define AIMD sequence
aimd_sequence = [
    {'temperature': 300, 'steps': 1000},
    {'temperature': 300, 'steps': 1000},
    {'temperature': 500, 'steps': 1000},
]

# Build workflow
wg = build_core_workgraph(
    # ... standard parameters ...
    run_aimd=True,
    aimd_sequence=aimd_sequence,
    aimd_parameters=aimd_parameters,
    aimd_options=slab_options,
)
```

## Output Structure

```
aimd_results/
  term_0/
    stage_00_300K_structure
    stage_00_300K_trajectory
    stage_00_300K_energy
    stage_00_300K_remote
    stage_00_300K_retrieved
    stage_01_300K_structure
    ...
    final_structure
    final_remote
    final_trajectory
  term_1/
    ...
```

## Examples

- Binary oxide: `examples/complete/complete_ag2o_aimd_example.py`
- Ternary oxide: `examples/complete/complete_ag3po4_aimd_example.py`

## VASP Parameters

AIMD uses NVT ensemble (Nosé-Hoover thermostat):
- `IBRION = 0` (MD mode)
- `MDALGO = 2` (NVT)
- `POTIM` (timestep in fs)
- `TEBEG/TEEND` (set automatically per stage)
- `NSW` (set automatically per stage)

## Notes

- Each stage restarts from previous using WAVECAR/CHGCAR
- Fixed atoms support coming in future update
- External restart capability coming in future update
