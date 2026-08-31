# Surface Energy WorkGraph Integration Test

Integration test for automatic surface energy calculations in hydroxylation workflows.

## What This Tests

- Complete WorkGraph execution with surface energy calculations
- Task dependencies and data flow
- Correct integration of all 3 formation reactions
- Output structure and accessibility

## Test System

- Small Ag3PO4 slab (minimal atoms for fast execution)
- Pristine + 1 OH modification
- Mocked VASP calculations (pre-computed energies)

## Running

```bash
source ~/envs/aiida/bin/activate
python test_integration.py
```

Expected runtime: ~30 seconds (all mocked, no real VASP)

## Success Criteria

- WorkGraph completes with exit code 0
- All 3 reaction results accessible via outputs
- Surface energies are positive (J/m²)
- Formation energies are reasonable (eV)
