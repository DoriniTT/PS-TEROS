# AIMD Restart Folders Fix

## Problem
After fixing the initial YAML serialization error, a second issue appeared when trying to run stage 1 with restart from stage 0:

```
Error in task aimd_stage_01_300K: type `<class 'aiida.orm.nodes.data.remote.base.RemoteData'>` 
is not supported as it is not json-serializable
```

## Root Cause
The `restart_folders` parameter had no type annotation:
```python
restart_folders=None
```

Without the proper type annotation, WorkGraph didn't know that `restart_folders` should be treated as input links (AiiDA nodes). Instead, it tried to serialize the RemoteData objects as JSON-serializable attributes, which failed.

## Solution

### 1. Added proper type annotation with empty dict as default
Changed from:
```python
restart_folders=None
```

To:
```python
restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {}
```

**Key insight**: Using `{}` (empty dict) as default instead of `None` avoids the Union/Optional type issue that causes "unhashable type: 'dict'" errors in Python's typing system.

### 2. Updated restart_folders checks
Changed from:
```python
if restart_folders is not None and slab_label in restart_folders:
```

To:
```python
if restart_folders and slab_label in restart_folders:
```

This works because:
- Empty dict `{}` evaluates to `False`
- Dict with items evaluates to `True`
- Both cases are handled correctly

### 3. Updated workgraph.py to initialize with empty dict
Changed from:
```python
current_remotes = None
```

To:
```python
current_remotes = {}  # Empty dict for first stage (no restart)
```

## Files Modified
1. `teros/core/aimd.py` - Fixed type annotation and restart_folders checks
2. `teros/core/workgraph.py` - Initialize current_remotes with empty dict

## Test Results
✅ **Workflow PK 35124** - Both AIMD stages running successfully:
- **aimd_stage_00_300K** (PK 35348): Finished successfully
- **aimd_stage_01_300K** (PK 35432): Running with restart from stage 0
- No JSON serialization errors
- RemoteData properly passed as input links

## Key Takeaway
When using optional AiiDA node inputs in @task.graph functions:
- ✅ **DO**: Use empty dict `{}` as default
- ✅ **DO**: Add proper type annotation with `dynamic()`
- ❌ **DON'T**: Use `None` as default (causes Union/Optional issues)
- ❌ **DON'T**: Omit type annotation (causes JSON serialization errors)

## Verification Commands
```bash
# Check workflow status
verdi process show 35124

# Check AIMD stage 1
verdi process show 35432

# Verify restart folders are being used
verdi process show 35432 | grep restart_folders -A 5
```

## Pattern for Sequential Tasks with Restart Capability
```python
@task.graph
def sequential_task(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    # ... other parameters
):
    results = {}
    for label, structure in structures.items():
        inputs = {'structure': structure}
        
        # Add restart if available
        if restart_folders and label in restart_folders:
            inputs['restart_folder'] = restart_folders[label]
        
        # Run calculation
        task = SomeWorkChain(**inputs)
        results[label] = task.outputs
    
    return results
```
