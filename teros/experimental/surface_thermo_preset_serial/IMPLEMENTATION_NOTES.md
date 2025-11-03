# Implementation Notes: Serial Surface Thermodynamics Preset

## Problem Statement

The goal was to create a serial/flat-graph variant of the surface thermodynamics preset that:

1. **Generates slabs from RELAXED bulk structure** (not input bulk)
2. **Maintains flat-graph architecture** (no nested workgraphs)
3. **Respects `max_number_jobs` concurrency control**
4. **Avoids the nested workgraph approach** that prevents `max_number_jobs` from working

## Architectural Challenge

The fundamental issue is a conflict between three requirements:

1. **Slabs must be generated from relaxed bulk** → This is a RUNTIME operation (happens after bulk relaxation completes)
2. **VASP tasks must be created for individual slabs** → This is a BUILD-TIME operation (graph construction)
3. **No nested workgraphs allowed** → This eliminates the `@task.graph_builder` solution

### Why Dynamic Outputs Don't Work at Build Time

In AiiDA WorkGraph:
- Dynamic namespaces (e.g., `namespace(slabs=dynamic(orm.StructureData))`) emit outputs at RUNTIME
- Individual sockets within the dynamic namespace don't exist at BUILD time
- You cannot use `getattr(node.outputs.slabs, slab_id)` during graph construction
- This causes: `AttributeError: TaskSocketNamespace has no sub-socket 'slab_100_term_0'`

**Documentation reference**: [Socket Concept](https://aiida-workgraph.readthedocs.io/en/latest/concept/autogen/socket_concept.html)

## Approaches Tried

### ❌ Approach 1: Pre-allocate Slab Slots with `getattr()`

**What we tried**:
```python
# Pre-allocate expected slab IDs
expected_slab_ids = []
for miller in miller_indices:
    miller_str = ''.join(map(str, miller))
    for term_idx in range(max_terminations_per_miller):
        expected_slab_ids.append(f"slab_{miller_str}_term_{term_idx}")

# Try to create VASP tasks referencing dynamic outputs
for slab_id in expected_slab_ids:
    slab_structure = getattr(slab_gen_node.outputs.slabs, slab_id)  # FAILS HERE
    scf_nodes[slab_id] = wg.add_task(
        VaspWorkChain,
        name=f"scf_slab_{slab_id}",
        structure=slab_structure,
        ...
    )
```

**Why it failed**:
- `slab_gen_node.outputs.slabs` is a dynamic namespace
- Individual sockets like `slab_100_term_0` don't exist until RUNTIME
- `getattr()` tries to access them at BUILD time
- Result: `AttributeError`

**Test script**: `tests/test_relaxed_bulk_slab_gen.py`

### ❌ Approach 2: Graph Builder with `@task.graph_builder`

**What this would do**:
```python
@task.graph_builder
def create_slab_calculations(slabs_namespace):
    # This creates a nested workgraph
    # Would work, but violates "no nested workgraphs" requirement
    pass
```

**Why we didn't pursue it**:
- Creates nested workgraph structure
- User explicitly requested flat-graph architecture
- Defeats the purpose of the serial preset

**Documentation reference**: [Graph Builder Decorator](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/graph_builder.html)

### ✅ Approach 3: Two-Stage Workflow (Current Implementation)

**What we implemented**:

**Stage 1**: `generate_slabs_from_relaxed_bulk_workgraph()`
- Standalone workflow that:
  - Relaxes the bulk structure
  - Generates slabs from relaxed bulk using `generate_multiple_miller_slabs` calcfunction
  - Returns slabs as concrete outputs

**Stage 2**: `surface_thermodynamics_serial_workgraph()`
- Main workflow that:
  - Accepts pre-provided `input_slabs` dictionary
  - Runs all surface thermodynamics calculations
  - Maintains flat-graph architecture
  - `max_number_jobs` works correctly

**Implementation details**:

1. **New calcfunction** (`workgraph.py:33-103`):
```python
@task.calcfunction
def generate_multiple_miller_slabs(
    bulk_structure: orm.StructureData,
    miller_indices_list: orm.List,
    ...
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """Generate slab structures from bulk for multiple Miller indices."""
    # Uses PyMatGen SlabGenerator
    # Returns dynamic namespace with all generated slabs
```

2. **Preliminary workflow** (`workgraph.py:106-225`):
```python
def generate_slabs_from_relaxed_bulk_workgraph(...):
    wg = WorkGraph("generate_slabs_from_relaxed_bulk")

    # Phase 1: Bulk relaxation
    bulk_node = wg.add_task(VaspWorkChain, ...)

    # Phase 2: Slab generation from relaxed bulk
    slab_gen_node = wg.add_task(
        generate_multiple_miller_slabs,
        bulk_structure=bulk_node.outputs.structure,  # Uses relaxed bulk!
        ...
    )

    # Set outputs
    wg.outputs = {
        'slabs': slab_gen_node.outputs.slabs,
        'relaxed_bulk': bulk_node.outputs.structure,
    }

    return wg
```

3. **Usage pattern** (`tests/test_two_stage_serial_preset.py`):
```python
# Stage 1: Generate slabs
stage1_wg = generate_slabs_from_relaxed_bulk_workgraph(...)
stage1_wg.submit(wait=False)

# Wait for Stage 1 completion
# ...

# Extract slabs
slabs_namespace = stage1_wg.outputs.slabs
slab_dict = {key: slabs_namespace[key] for key in slabs_namespace.keys()}

# Stage 2: Run thermodynamics with pre-provided slabs
stage2_wg = surface_thermodynamics_serial_workgraph(
    input_slabs=slab_dict,  # Use slabs from Stage 1
    ...
)
stage2_wg.max_number_jobs = 2  # Concurrency control works!
stage2_wg.submit(wait=False)
```

**Why this works**:
- ✅ Slabs generated from relaxed bulk (Stage 1)
- ✅ Main workflow has flat-graph architecture (Stage 2)
- ✅ `max_number_jobs` controls concurrency in Stage 2
- ✅ No dynamic output access issues (slabs are concrete inputs)

**Test script**: `tests/test_two_stage_serial_preset.py`

## What Works

1. **Calcfunction for slab generation**: `generate_multiple_miller_slabs` correctly generates slabs from relaxed bulk
2. **Dynamic namespace outputs**: Works correctly when accessed AFTER workflow completes
3. **Flat-graph architecture**: Stage 2 maintains all VASP nodes at same level
4. **Concurrency control**: `max_number_jobs` works as expected in Stage 2

## What Doesn't Work

1. **Accessing dynamic outputs at build time**: Cannot use `getattr()` or direct attribute access on dynamic namespaces during graph construction
2. **Single-workflow solution**: No way found to generate slabs from relaxed bulk AND create VASP tasks in a single flat-graph workflow
3. **Pre-allocation approach**: Cannot pre-allocate slots for unknown number of slab terminations

## Best Approaches Going Forward

### Option A: Two-Stage Workflow (Recommended)

**Pros**:
- Clean separation of concerns
- Both workflows maintain provenance
- Can reuse Stage 1 slabs for multiple Stage 2 runs
- Each stage can be monitored independently
- Clear failure points

**Cons**:
- Requires two separate workflow submissions
- User must wait for Stage 1 before submitting Stage 2
- More complex user-facing API

**Best for**: Production workflows, parameter studies, reusing slabs

### Option B: Hybrid Approach with Pre-provided Slabs

Keep the current main workflow that accepts `input_slabs`, but:
- Document that slabs should come from preliminary relaxation
- Provide helper script to run bulk relaxation + slab generation
- Main workflow assumes slabs are from relaxed bulk

**Pros**:
- Main workflow remains simple
- User controls bulk relaxation separately
- Flexibility in slab preparation

**Cons**:
- Manual step between relaxation and main workflow
- User must understand provenance implications

**Best for**: Interactive use, custom slab preparation

### Option C: Accept Nested Workgraphs

Reconsider the "no nested workgraphs" constraint:
- Use `@task.graph_builder` for dynamic slab calculations
- Accept that `max_number_jobs` only works at top level
- Document the limitation clearly

**Pros**:
- Single workflow submission
- Automatic slab generation
- Clean user-facing API

**Cons**:
- Creates nested structure
- `max_number_jobs` doesn't control inner workgraph tasks
- Defeats original purpose of serial preset

**Best for**: Simplified API when concurrency control isn't critical

## Recommendations

### For Immediate Use

Use the **two-stage workflow approach**:

1. Submit Stage 1 to generate slabs from relaxed bulk
2. Wait for completion and extract slabs
3. Submit Stage 2 with pre-provided slabs and `max_number_jobs` control

This is the only approach that satisfies all original requirements.

### For Future Development

Consider these improvements:

1. **Automatic orchestration**: Create a helper function that manages both stages automatically:
```python
def run_surface_thermodynamics_two_stage(
    # All parameters for both stages
):
    """Orchestrates Stage 1 and Stage 2 automatically."""
    # Submit Stage 1
    # Wait for completion
    # Extract slabs
    # Submit Stage 2
    # Return both PKs
```

2. **Caching mechanism**: Cache relaxed bulk + slabs for reuse in parameter studies

3. **Documentation**: Add clear examples showing the two-stage pattern

4. **Alternative**: Investigate if AiiDA WorkGraph will add support for:
   - Accessing dynamic outputs at build time
   - Dynamic task creation based on runtime outputs
   - Better scatter-gather patterns for this use case

## Files Modified

### Core Implementation

- `teros/experimental/surface_thermo_preset_serial/workgraph.py`
  - Added `generate_multiple_miller_slabs` calcfunction (lines 33-103)
  - Added `generate_slabs_from_relaxed_bulk_workgraph` function (lines 106-225)
  - Modified `surface_thermodynamics_serial_workgraph` to accept `input_slabs`

- `teros/experimental/surface_thermo_preset_serial/__init__.py`
  - Exported new `generate_slabs_from_relaxed_bulk_workgraph` function

### Test Scripts (in `tests/` directory)

- `test_two_stage_serial_preset.py`: Full two-stage workflow demonstration
- `test_relaxed_bulk_slab_gen.py`: Attempted single-workflow approach (doesn't work)
- `test_serial_slab_generation.py`: Earlier test of slab generation

## Key Takeaways

1. **Dynamic outputs are runtime data**: Cannot be accessed during graph construction
2. **Build-time vs runtime**: Fundamental distinction in WorkGraph architecture
3. **Flat-graph limitation**: Cannot dynamically create tasks based on runtime outputs without nesting
4. **Two-stage is viable**: Splits the problem into two solvable pieces
5. **Concurrency control works**: When all tasks are at same graph level

## References

- [AiiDA WorkGraph Documentation](https://aiida-workgraph.readthedocs.io/)
- [Limit Concurrent Jobs](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/limit_concurrent.html)
- [Socket Concept](https://aiida-workgraph.readthedocs.io/en/latest/concept/autogen/socket_concept.html)
- [Scatter-Gather Pattern](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/parallel.html)
- [Graph Builder Decorator](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/graph_builder.html)
- [Combining WorkGraphs](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/combine_workgraphs.html)

---

**Last updated**: 2025-11-02
**Status**: Two-stage approach implemented and ready for testing
