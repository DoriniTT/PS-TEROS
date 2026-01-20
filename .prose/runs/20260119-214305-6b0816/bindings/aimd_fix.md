# aimd_fix

kind: let

source:
```prose
aimd_fix = session: fixer
  prompt: """Update teros/core/aimd/tasks.py..."""
```

---

Updated `teros/core/aimd/tasks.py` to use modern task decorators:

1. **Removed** `from aiida.engine import calcfunction` import
2. **Changed** `@calcfunction` to `@task.calcfunction` decorator
3. **Renamed** function from `create_supercell_calcfunc` to `create_supercell`
4. **Removed** the wrapper line `create_supercell = task(create_supercell_calcfunc)`
5. **Verified** `create_supercells_scatter` already uses `create_supercell` (no changes needed)

The file now follows the modern pattern where `@task.calcfunction` directly decorates the function, eliminating the need for a separate wrapper step.
