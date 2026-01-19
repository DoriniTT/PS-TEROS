# summary

kind: output

source:
```prose
output summary = session "Verify decorator standardization"
  context: { aimd_fix, surface_energy_fix }
```

---

## Verification Results

### 1. teros/core/aimd/tasks.py

| Check | Status |
|-------|--------|
| No `from aiida.engine import calcfunction` | PASS |
| Uses `@task.calcfunction` decorator | PASS |
| flake8 (max-line-length=120, ignore=E501,W503,E402,F401) | PASS (no errors) |

**Details:**
- Import: `from aiida_workgraph import task, dynamic, namespace` (line 5)
- Decorator: `@task.calcfunction` on `create_supercell` (line 8)
- Function renamed from `create_supercell_calcfunc` to `create_supercell`


### 2. teros/core/surface_hydroxylation/surface_energy.py

| Check | Status |
|-------|--------|
| No `from aiida.engine import calcfunction` | PASS |
| All 8 functions use `@task.calcfunction` | PASS |
| flake8 (max-line-length=120, ignore=E501,W503,E402,F401) | 1 warning (unrelated) |

**Details:**
- Import: `from aiida_workgraph import task` (line 9)
- Decorated functions:
  - `calc_delta_g_general_reaction1` (line 243)
  - `calc_delta_g_general_reaction2` (line 275)
  - `calc_delta_g_general_reaction3` (line 311)
  - `calc_delta_g_reaction1` (line 346)
  - `calc_delta_g_reaction2` (line 396)
  - `calc_delta_g_reaction3` (line 436)
  - `calc_gamma_s` (line 541)
  - `calc_gamma` (line 577)

**Flake8 note:** Line 693 has `F541 f-string is missing placeholders` - this is a pre-existing issue unrelated to the decorator changes.


## Summary

All decorator standardization changes were successful:

- Both files now use `from aiida_workgraph import task` instead of `from aiida.engine import calcfunction`
- All decorated functions use `@task.calcfunction` pattern
- No linting errors introduced by the changes
