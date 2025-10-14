# CP2K AIMD Integration Debugging Guide

## Problem Statement

When running `step_07a_aimd_autogenerate_slabs.py`, the workflow fails during the AIMD stage with:

```
TypeError: 'int' object is not iterable
```

The error occurs in AiiDA's serialization layer (`clean_value` function) when trying to store the CP2K task node.

## Understanding the Error

The error happens because AiiDA's node storage system expects certain data types to be wrapped in AiiDA ORM types (like `orm.Int`, `orm.Bool`, `orm.Dict`), but somewhere in our code we're passing plain Python types (like `int`, `bool`, `dict`) where AiiDA's serialization expects something iterable or an AiiDA node.

## Root Cause Analysis Strategy

### 1. Check the Error Traceback

The traceback shows:
```
File "/aiida/orm/implementation/utils.py", line 108, in clean_value
    for key, value in value.items():
TypeError: 'int' object is not iterable
```

This means `clean_value` is trying to iterate over something expecting a dictionary, but it's getting an integer instead.

### 2. Compare with Working Example

Look at `/home/thiagotd/projects/fosfato-h2o/calculos/hydrosurface/_110/st-2/calculo/wg_AIMD_CP2K_3D_single_slab_updated.py`:

**Working code (lines 504-505):**
```python
aimd_task = Cp2kTask(
    cp2k=cp2k_inputs,
    max_iterations=orm.Int(builder_inputs['max_iterations']),
    clean_workdir=orm.Bool(builder_inputs['clean_workdir']),
)
```

**Broken code (teros/core/aimd_cp2k.py lines 118-122):**
```python
aimd_task = Cp2kTask(
    cp2k=cp2k_inputs,
    max_iterations=orm.Int(3),
    clean_workdir=orm.Bool(clean_workdir),
)
```

Both wrap the parameters correctly, so the issue is NOT here.

### 3. Check the CP2K Input Structure

The problem is likely in the `cp2k_inputs` dictionary or the `aimd_parameters` dictionary that gets passed to it.

**Current structure (aimd_cp2k.py lines 93-111):**
```python
cp2k_inputs = {
    'structure': slab_structure,
    'parameters': orm.Dict(dict=stage_params),
    'code': code,
    'metadata': {
        'options': dict(options),  # ← POTENTIAL ISSUE
    },
    'file': {
        'basis': basis_file,
        'pseudo': pseudo_file,
    },
    'settings': orm.Dict(dict={...}),
}
```

**Working example structure:**
```python
cp2k_inputs = {
    'structure': prev_structure,
    'parameters': orm.Dict(dict=params_with_fixed),
    'code': builder_inputs['code'],
    'metadata': builder_inputs['metadata'],  # ← DIFFERENT!
    'file': builder_inputs['file'],
    'settings': orm.Dict(dict=builder_inputs['settings']),
}
```

**Key difference:** The working example passes `metadata` directly, while our code nests it under `'options'`.

### 4. Check How Options Are Passed

Look at workgraph.py line 1426:
```python
'options': aimd_opts,
```

Then in aimd_cp2k.py line 98:
```python
'metadata': {
    'options': dict(options),
}
```

This creates: `metadata: {options: {resources: {...}, queue_name: '...'}}`

But the working example has:
```python
metadata = {
    'options': {
        'resources': resources,
        'queue_name': QUEUE_NAME,
    }
}
```

So `builder_inputs['metadata']` already contains the full metadata structure.

## The Actual Problem

**HYPOTHESIS:** The `options` parameter passed to `aimd_single_stage_scatter_cp2k` is already the full metadata dictionary, not just the options sub-dictionary.

Looking at the function signature (aimd_cp2k.py line 26):
```python
options: t.Annotated[dict, Cp2kBaseWorkChain.spec().inputs['metadata']],
```

This annotation says `options` should be the FULL metadata structure, not just the scheduler options!

## Solution

### Option 1: Rename Parameter (Clearest)

Change the parameter name from `options` to `metadata` to match its actual content:

**In aimd_cp2k.py:**
```python
def aimd_single_stage_scatter_cp2k(
    slabs: ...,
    temperature: float,
    steps: int,
    code: orm.Code,
    aimd_parameters: dict,
    basis_file: orm.SinglefileData,
    pseudo_file: orm.SinglefileData,
    metadata: t.Annotated[dict, Cp2kBaseWorkChain.spec().inputs['metadata']],  # ← RENAMED
    clean_workdir: bool,
    ...
):
    ...
    cp2k_inputs = {
        'structure': slab_structure,
        'parameters': orm.Dict(dict=stage_params),
        'code': code,
        'metadata': metadata,  # ← PASS DIRECTLY
        'file': {
            'basis': basis_file,
            'pseudo': pseudo_file,
        },
        'settings': orm.Dict(dict={...}),
    }
```

**In workgraph.py line 1426:**
```python
stage_inputs = {
    ...
    'metadata': {  # ← BUILD FULL METADATA STRUCTURE
        'options': aimd_opts,
    },
    ...
}
```

### Option 2: Keep Current Name, Fix Usage

If keeping the parameter name as `options`:

**In aimd_cp2k.py line 97-99:**
```python
'metadata': {
    'options': options,  # ← REMOVE dict() call
},
```

But the annotation is misleading since it says the parameter IS metadata, not just options.

## Debugging Steps to Confirm

1. **Add debug output before task creation:**

```python
# In aimd_cp2k.py, before line 118
print(f"DEBUG: Creating CP2K task for {slab_label}")
print(f"  options type: {type(options)}")
print(f"  options content: {options}")
print(f"  stage_params type: {type(stage_params)}")
```

2. **Check what gets passed from workgraph.py:**

```python
# In workgraph.py, before line 1437
print(f"DEBUG: AIMD stage {stage_idx} inputs:")
print(f"  options: {aimd_opts}")
```

3. **Run the test and check output:**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/aiida/bin/activate
verdi daemon restart
python examples/cp2k/step_07a_aimd_autogenerate_slabs.py
```

4. **Check the daemon log if no debug output:**

```bash
tail -100 ~/.aiida/daemon/log/daemon-40000.log | grep -A 20 "DEBUG:"
```

## Additional Checks

### Check for Unwrapped Integers in CP2K Parameters

Look at `teros/core/builders/aimd_builder_cp2k.py`:

```python
def get_aimd_defaults_cp2k(...):
    return {
        "GLOBAL": {
            "RUN_TYPE": "MD",
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "MD": {
                "ENSEMBLE": "NVT",
                "TIMESTEP": timestep,  # float
                "THERMOSTAT": {
                    "TYPE": thermostat,  # string
                    "REGION": "GLOBAL",
                },
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "DFT": {
                ...
                "MGRID": {
                    "CUTOFF": cutoff,      # float
                    "REL_CUTOFF": rel_cutoff,  # float
                    "NGRIDS": 4            # ← PLAIN INTEGER!
                },
                ...
                "SCF": {
                    ...
                    "EPS_SCF": eps_scf,    # float
                    "MAX_SCF": max_scf,    # ← PLAIN INTEGER!
                    ...
                },
                ...
            },
        },
    }
```

These integers (`NGRIDS`, `MAX_SCF`) are in the nested dictionary that gets passed to `orm.Dict(dict=stage_params)`.

**This is OK** because they're inside a Python dict that gets wrapped in `orm.Dict()`. The serialization happens correctly.

## Most Likely Issue

The problem is the metadata/options nesting. The parameter is annotated as `metadata` but named `options`, and we're double-nesting it.

## Quick Fix to Test

Edit `teros/core/aimd_cp2k.py` line 97-99:

```python
# OLD (BROKEN):
'metadata': {
    'options': dict(options),
},

# NEW (FIX):
'metadata': options,
```

This assumes `options` parameter already contains the full metadata structure with `{'options': {...}}` inside it.

## Verification

After applying the fix:

1. Clear Python cache:
```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

2. Restart daemon:
```bash
verdi daemon restart
```

3. Run test:
```bash
python examples/cp2k/step_07a_aimd_autogenerate_slabs.py
```

4. Wait 30 seconds and check status:
```bash
verdi process status <PK>
```

If the error persists, check the daemon log for the actual line that's failing and work backward from there.

## Alternative: Check cp2k.base Workchain Spec

If still stuck, examine the actual CP2K workchain specification:

```python
from aiida.plugins import WorkflowFactory
Cp2kWorkChain = WorkflowFactory('cp2k.base')
print(Cp2kWorkChain.spec().inputs)
```

This will show exactly what structure it expects for the inputs.
