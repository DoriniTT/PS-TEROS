# IMPORTANT NOTES: User-Provided Slabs Feature

## Critical Requirements

### ⚠️ **MUST Store Structures Before Passing**

When using user-provided slabs, you **MUST** call `.store()` on each `StructureData` node before adding it to the `input_slabs` dictionary:

```python
# ❌ WRONG - Will cause errors
input_slabs = {
    "term_0": orm.StructureData(ase=read("slab.cif")),
}

# ✅ CORRECT - Store first
structure = orm.StructureData(ase=read("slab.cif"))
structure.store()  # This is required!
input_slabs = {
    "term_0": structure,
}
```

### Why is this required?

The implementation uses post-build manual task injection to avoid provenance cycle issues. Stored structures can be safely passed as task inputs without creating serialization problems.

## How It Works Internally

1. **Build Phase**: WorkGraph is built with `use_input_slabs=True` flag
   - Slab generation is skipped (`slab_namespace = None`)
   - Graph built without scatter task for slabs

2. **Post-Build Phase**: After graph is built
   - `relax_slabs_scatter` task is manually added to the WorkGraph
   - User slabs (stored nodes) are passed directly as task inputs
   - Task outputs are wired to graph outputs

3. **Execution Phase**: Workflow runs normally
   - Scatter task receives stored slab structures
   - Creates parallel VASP calculations for each slab
   - No slab generation task is called

## Verification

To verify your slabs are being used correctly:

```bash
# After submission, check the workflow
verdi process show <PK>

# Look for:
# 1. NO "generate_slab_structures" in Called tasks
# 2. "relax_slabs_scatter" with your slab PKs as inputs
#3. Check the scatter task shows your exact slab PKs
verdi process show <scatter_PK>
```

## Common Issues

### Issue 1: "Cycle in graph" error
**Cause**: Trying to pass unstored structures  
**Solution**: Call `.store()` on each structure before adding to dict

### Issue 2: Wrong slabs being relaxed
**Cause**: Implementation error in earlier versions  
**Solution**: Update to latest version, verify no "generate_slab_structures" task

### Issue 3: Serialization error
**Cause**: Passing dict of nodes through build parameters  
**Solution**: This is handled internally now, ensure you're using latest code

## Best Practices

1. **Always store structures immediately after creation**:
   ```python
   structure = orm.StructureData(ase=atoms)
   structure.store()
   ```

2. **Use descriptive keys**:
   ```python
   input_slabs = {
       "oxygen_terminated": structure1,
       "metal_terminated": structure2,
   }
   ```

3. **Validate structures before storing**:
   ```python
   atoms = read("slab.cif")
   # Check formula, cell, positions, etc.
   print(atoms.get_chemical_formula())
   print(atoms.cell)
   # Then create and store
   structure = orm.StructureData(ase=atoms)
   structure.store()
   ```

4. **Keep track of PKs**:
   ```python
   for key, structure in input_slabs.items():
       print(f"{key}: PK {structure.pk}")
   ```

## Implementation Architecture

This feature differs from automatic generation in its architecture:

**Automatic Generation**:
- Slab generation task is part of `@task.graph`
- Runs during graph execution
- Output connects to scatter task automatically

**User-Provided**:
- No generation task created
- Scatter task added post-build
- Manual wiring of inputs/outputs
- Avoids provenance cycles with stored nodes

Both produce identical outputs and user experience, but internal execution differs.

## Version Information

- Feature added: PS-TEROS v0.2.0
- Current implementation: Post-build injection with `use_input_slabs` flag
- Tested with: AiiDA 2.x, AiiDA-WorkGraph latest

## See Also

- `IMPLEMENTATION_SUMMARY.md` - Technical details
- `README_INPUT_SLABS.md` - User guide
- `QUICKSTART.md` - Quick start
- `slabs_input_relax.py` - Working example
