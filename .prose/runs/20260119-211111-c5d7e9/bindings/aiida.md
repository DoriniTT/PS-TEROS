# aiida

kind: let

source:
```prose
aiida = session: integration
  prompt: "Review PS-TEROS AiiDA integration..."
```

---

# AiiDA/WorkGraph Integration Review: PS-TEROS

**Overall Grade: A-** (92/100)

---

## Key Findings

### 1. Task Decorator Usage (Grade: A+)

**Proper @task.calcfunction usage:**
- Correctly used for data transformations
- All return AiiDA data types
- Properly annotated with namespace for dynamic outputs

**Correct @task.graph pattern:**
- Used for composite workflows and scatter-gather
- Proper type hints with `t.Annotated[dict, namespace(...)]`

**Minor Issue:**
- Some files mix old `@calcfunction` + `task()` wrapper with new `@task.calcfunction`
- Standardize to `@task.calcfunction` throughout

### 2. Output Socket Patterns (Grade: A-)

**Correct pattern used:**
```python
return {
    'bulk_energy': bulk_energy.result,
    'bulk_structure': bulk_vasp.structure,
    'formation_enthalpy': formation_hf.result,
}
```

Outputs correctly appear in `verdi process show`.

### 3. Scatter-Gather Implementation (Grade: A+)

**Textbook implementation:**
```python
@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
)]:
    for label, structure in slabs.items():
        relaxation = relax_task_cls(**relax_inputs)
        relaxed[label] = relaxation.structure
```

**Strengths:**
- Dynamic type annotations properly used
- Independent task creation for parallel execution
- Per-structure parameter overrides via deep_merge_dicts

### 4. VaspWorkChain Integration (Grade: A)

**Builder pattern:**
- Backward compatible (old-style + new builder_inputs)
- Deep copy prevents mutation bugs
- Consistent settings application

### 5. Provenance Tracking (Grade: A+)

- Complete provenance for formation enthalpy, surface energy, cleavage
- Links to bulk, slab, and reference calculations
- Chemical potential grids stored
- Oxide type identification in provenance

### 6. Error Handling (Grade: B+)

**Strengths:**
- Early validation with informative messages
- Physics-aware validation (stoichiometry checks)

**Gaps:**
- No retry logic for VASP convergence failures
- Missing validation that `structure` output exists (NSW=0 case)

---

## Notable Innovations

### 1. Placeholder Pattern
```python
if compute_formation_enthalpy:
    # actual calculations
else:
    metal_vasp = TaskOutputPlaceholder(structure=bulk_vasp.structure)
```
Satisfies type requirements without dummy VASP calculations.

### 2. Three-Level Override Hierarchy
1. Global defaults
2. Structure-specific overrides
3. Always-override (structure, code)

### 3. Adsorbate Separation with Graph Analysis
- CrystalNN + distance fallback for peroxide bonds
- Robust handling of both covalent and weak bonds

---

## Priority Recommendations

### High Priority
1. Standardize task decorators (@task.calcfunction everywhere)
2. Add VASP failure handlers with retry logic

### Medium Priority
3. Output socket consistency check for NSW=0
4. Document scatter-gather pattern for users

### Low Priority
5. Type hint consistency (orm.Code vs load_code())
