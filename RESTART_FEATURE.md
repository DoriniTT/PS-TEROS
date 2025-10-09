# Restart Feature - Phase 1: Expose RemoteData Nodes

## Overview
This document describes the implementation of Phase 1 of the restart feature, which exposes RemoteData nodes from slab relaxations to the main workgraph outputs.

## Changes Made

### 1. Modified `teros/core/slabs.py`

**Function:** `relax_slabs_scatter`

**Changes:**
- Added `remote_folders=dynamic(orm.RemoteData)` to the return type annotation
- Created a `remote_folders` dictionary to collect RemoteData nodes from each slab relaxation
- For each VASP relaxation task, now captures `relaxation.remote_folder` in addition to structure and energy
- Returns the `remote_folders` dictionary in the output

**Purpose:** Captures and returns the RemoteData nodes from each individual slab relaxation calculation.

### 2. Modified `teros/core/workgraph.py`

**Function:** `core_workgraph`

**Changes:**
- Added `'slab_remote'` to the outputs list in the `@task.graph` decorator
- Initialized `slab_remote_output = {}` variable
- When `relax_slabs=True`, extracts `remote_folders` from `relaxation_outputs`
- Returns `slab_remote_output` in the output dictionary
- Updated docstring to document the new `slab_remote` output

**Function:** `build_core_workgraph`

**Changes:**
- In the section handling `use_input_slabs=True`, added connection for remote folders:
  - `wg.outputs.slab_remote = scatter_task.outputs.remote_folders`

**Purpose:** Exposes the RemoteData nodes at the main workgraph level, making them accessible in the workflow outputs.

## How to Use

After running a PS-TEROS calculation with `relax_slabs=True`, you can now access the RemoteData nodes:

```bash
# View the main workgraph node
verdi process show <PK>

# The outputs section will now include:
# slab_remote
#     term_0          <PK>  RemoteData
#     term_1          <PK>  RemoteData
#     term_2          <PK>  RemoteData
#     ...
```

You can access individual RemoteData nodes:

```python
from aiida import orm

# Load the main workgraph node
wg = orm.load_node(18818)

# Access remote folders for each slab
remote_term_0 = wg.outputs.slab_remote.term_0
remote_term_1 = wg.outputs.slab_remote.term_1

# Get the remote path
print(remote_term_0.get_remote_path())
```

## Next Steps: Phase 2 - Restart Functionality

Phase 2 will implement the actual restart capability:

1. Add parameters to accept RemoteData nodes from a previous run
2. Modify the VASP workchain calls to use the RemoteData for restarting calculations
3. Implement logic to skip completed calculations and only restart failed/incomplete ones
4. Add validation to ensure RemoteData nodes are compatible with the current calculation

## Testing

To test the changes:

1. **Clear Python cache** (IMPORTANT - daemon uses cached .pyc files):
   ```bash
   find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
   find . -name "*.pyc" -delete 2>/dev/null
   ```

2. **Restart the AiiDA daemon:**
   ```bash
   verdi daemon restart
   ```

3. **Run a PS-TEROS calculation** with `relax_slabs=True`

4. **After completion, check the outputs:**
   ```bash
   verdi process show <PK>
   ```

5. **Verify that `slab_remote` appears** in the outputs with RemoteData nodes for each slab

### Example Output

```
Outputs             PK     Type
------------------  -----  -------------
slab_remote
    term_0          19858  RemoteData
    term_1          19865  RemoteData
```

### Verified Working

✓ Tested with Ag2O (100) slab relaxation (PK 19774)  
✓ RemoteData nodes successfully exposed for each slab termination  
✓ Remote paths accessible on compute cluster (bohr)

## Benefits

- **Restart Capability:** Enables restarting failed calculations from the last checkpoint
- **Debugging:** Allows access to calculation files on remote machines for troubleshooting
- **Analysis:** Facilitates post-processing analysis using files from the remote folder
- **Efficiency:** Sets foundation for implementing incremental calculations without re-running completed steps
