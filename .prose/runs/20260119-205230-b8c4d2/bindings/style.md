# style

kind: let

source:
```prose
style = session: reviewer
  prompt: "Review for code style and readability"
```

---

# PS-TEROS Code Style and Readability Review

**Date:** 2026-01-19
**Reviewer:** Code Analysis System
**Scope:** Full codebase analysis focusing on naming conventions, organization, documentation, and consistency

## Executive Summary

PS-TEROS demonstrates **strong overall code quality** with well-organized architecture, comprehensive documentation, and consistent patterns. The codebase shows evidence of careful design, particularly in the WorkGraph-based workflow orchestration and the modular structure. However, there are opportunities for improvement in naming consistency, documentation completeness, and code organization in specific areas.

**Overall Rating: 8.5/10**

---

## 1. Naming Conventions

### ✅ Strengths

1. **Module Names**: Excellent use of descriptive, purpose-driven names
   - `thermodynamics.py` - Clear domain indication
   - `formation_enthalpy.py` → `hf.py` - Concise abbreviation for internal use
   - `adsorption_energy.py` - Self-documenting
   - `workflow_presets.py` - Intent is immediately clear

2. **Function Names**: Generally follow PEP 8 snake_case with clear verb patterns
   ```python
   # Good examples:
   calculate_formation_enthalpy()
   generate_slab_structures()
   extract_total_energy()
   deep_merge_dicts()
   identify_oxide_type()
   ```

3. **Constant Names**: Proper UPPER_SNAKE_CASE with semantic meaning
   ```python
   EV_PER_ANGSTROM2_TO_J_PER_M2 = 16.0217663
   STOICHIOMETRY_RTOL = 1e-3
   DEFAULT_KPOINTS_SPACING = 0.03
   ```

4. **Class Names**: Appropriate PascalCase with clear purpose
   ```python
   TaskOutputPlaceholder
   EnergyOutputPlaceholder
   FormationEnthalpyPlaceholder
   JanafDatabase
   WorkGraphAnalysis
   ```

### ⚠️ Areas for Improvement

1. **Inconsistent Abbreviation Usage**
   - `hf.py` (formation enthalpy) - abbreviated
   - `thermodynamics.py` - full name
   - `adsorption_energy.py` - full name

   **Recommendation:** Either use full names consistently or establish abbreviation guidelines (e.g., abbreviate only when >3 words)

2. **Prefix/Suffix Inconsistency**
   ```python
   # Mixed patterns:
   build_core_workgraph()           # build_ prefix
   get_ag2o_defaults()              # get_ prefix
   calculate_formation_enthalpy()   # calculate_ prefix
   compute_surface_energies_scatter() # compute_ prefix
   ```

   **Recommendation:** Standardize verb prefixes:
   - `build_*` for WorkGraph constructors
   - `calculate_*` for @task.calcfunction operations
   - `compute_*` for high-level computational workflows
   - `get_*` for data retrieval/accessor functions

3. **Variable Naming in Long Functions**
   - Some parameter-heavy functions (e.g., `core_workgraph` with 50+ parameters) could benefit from parameter grouping
   - Variable names like `ref_data` and `formation_data` are generic

4. **Magic String Keys**
   ```python
   # Found in multiple files:
   misc_dict['total_energies']['energy_extrapolated']
   bulk_counts[element_M]
   ```

   **Recommendation:** Define string constants for frequently used dictionary keys

---

## 2. Code Organization

### ✅ Strengths

1. **Excellent Module Structure**: Clear separation of concerns
   ```
   teros/core/
   ├── workgraph.py          # Main orchestration (1992 lines - appropriate for entry point)
   ├── thermodynamics.py     # Domain-specific calculations
   ├── slabs.py              # Slab operations
   ├── hf.py                 # Formation enthalpy
   ├── cleavage.py           # Cleavage energy
   ├── constants.py          # Centralized constants
   ├── utils.py              # Shared utilities
   └── convergence/          # Submodule pattern
       ├── __init__.py
       ├── workgraph.py
       ├── tasks.py
       └── slabs.py
   ```

2. **Consistent Submodule Pattern**: Well-defined structure for complex features
   - Each submodule has: `__init__.py`, `workgraph.py`, `tasks.py`, `utils.py`
   - Clean export hierarchy through `__init__.py`
   - Examples: `fukui/`, `convergence/`, `surface_energy/`, `aimd/`

3. **Import Organization**: Generally follows PEP 8
   ```python
   # Good example from thermodynamics.py:
   from __future__ import annotations  # Future compatibility

   import typing as t                   # stdlib
   from collections import Counter
   from functools import reduce
   from math import gcd

   import numpy as np                   # third-party
   from aiida import orm
   from aiida_workgraph import dynamic, namespace, task
   ```

4. **Separation of Concerns**:
   - Pure calculations in `@task.calcfunction`
   - Workflow logic in `@task.graph`
   - Configuration in `builders/`
   - Utilities in `utils.py`
   - Physical constants in `constants.py`

### ⚠️ Areas for Improvement

1. **File Length Concerns**
   - `workgraph.py`: 1992 lines - very long but acceptable as main entry point
   - `adsorption_energy.py`: 1647 lines - could be split into multiple files
   - `slabs.py`: 848 lines - consider extracting helper functions

   **Recommendation:** For files >1000 lines:
   - Extract helper functions to separate `_helpers.py` or `_internal.py`
   - Consider splitting by feature (e.g., `slabs.py` → `slabs/generation.py`, `slabs/relaxation.py`, `slabs/analysis.py`)

2. **Duplicate Functionality**
   ```python
   # Found in multiple files:
   def get_settings():  # Appears in workgraph.py and surface_energy/workgraph.py
       return {'parser_settings': {...}}

   def load_structure_from_file(filepath: str):  # Multiple locations
       atoms = read(filepath)
       return orm.StructureData(ase=atoms)
   ```

   **Recommendation:** Consolidate duplicate utilities in `utils.py`

3. **Internal vs Public API Clarity**
   - Some helper functions lack `_` prefix despite being internal
   - `__all__` exports are comprehensive but could be organized by category

   **Recommendation:**
   - Use `_` prefix for internal functions: `_prepare_fukui_inputs()` (good example)
   - Group `__all__` exports by category with comments

4. **Experimental Code Location**
   ```
   teros/experimental/
   ├── surface_thermo_preset_serial/
   ├── max_jobs_investigation/
   └── zone_approach/
   ```

   **Recommendation:** Add clear deprecation notices or promotion path to stable API

---

## 3. Documentation Quality

### ✅ Strengths

1. **Comprehensive Module Docstrings**
   ```python
   # Excellent example from thermodynamics.py:
   """Ab initio atomistic thermodynamics for oxide surfaces.

   This module implements surface energy calculations as a function of chemical potential
   for both binary and ternary oxide slabs, following the ab initio atomistic thermodynamics
   framework.
   """
   ```

2. **Detailed Function Docstrings**: Google-style format with examples
   ```python
   def deep_merge_dicts(base: dict, override: dict) -> dict:
       """
       Deep merge override dict into base dict.

       Args:
           base: Base dictionary
           override: Override dictionary (values take precedence)

       Returns:
           Merged dictionary (new dict, inputs are not modified)

       Example:
           >>> base = {'a': 1, 'b': {'c': 2, 'd': 3}}
           >>> override = {'b': {'c': 99}, 'e': 5}
           >>> result = deep_merge_dicts(base, override)
           >>> result
           {'a': 1, 'b': {'c': 99, 'd': 3}, 'e': 5}
       """
   ```

3. **Physical Context in Documentation**
   - Constants include derivations and references
   - Thermodynamic equations documented with mathematical notation
   - VASP parameter explanations with usage guidance

4. **Excellent CLAUDE.md Documentation**
   - Clear module structure overview
   - Practical examples with working code
   - Troubleshooting tables
   - Architecture patterns explained

### ⚠️ Areas for Improvement

1. **Missing Type Hints in Some Functions**
   ```python
   # Inconsistent typing:
   def get_settings():  # Missing return type hint
       return {...}

   # Should be:
   def get_settings() -> dict:
       return {...}
   ```

2. **Incomplete Parameter Documentation**
   - Some functions with 20+ parameters have partial docstrings
   - `core_workgraph()` function has 50+ parameters but docstring is truncated

   **Recommendation:** Ensure every parameter is documented, even with brief descriptions

3. **Missing Docstrings for Internal Functions**
   ```python
   # Found in several files:
   def _prepare_convergence_inputs(builder_inputs: dict, code: orm.InstalledCode) -> dict:
       # Good implementation but missing docstring
   ```

4. **Inconsistent Example Quality**
   - Some functions have detailed examples
   - Others lack examples despite complexity
   - Test files exist but not always linked in docstrings

   **Recommendation:** Add "See Also" sections pointing to tests and examples

5. **TODO Comments**
   - Only 1 TODO found: `teros/core/aimd/utils.py:128: # TODO: Implementation in next task`
   - Generally good, but should track in issue system instead

---

## 4. Consistency Patterns

### ✅ Strengths

1. **Consistent WorkGraph Pattern**
   ```python
   @task.graph(outputs=['result1', 'result2'])
   def my_workflow(param: orm.StructureData) -> dict:
       """Clear pattern used throughout codebase"""
       task1 = VaspTask(...)
       result = process_task(output=task1.outputs.misc)
       return {'result1': result.outputs.value}
   ```

2. **Uniform Error Handling**
   ```python
   # Consistent validation pattern:
   if 'O' not in bulk_counts:
       raise ValueError('Structure contains no oxygen; not an oxide.')
   ```

3. **Predictable Function Signatures for Similar Operations**
   - All scatter-gather functions follow same pattern
   - All `@task.calcfunction` decorators have consistent return types
   - All WorkGraph builders accept similar parameter structures

4. **Consistent Use of AiiDA Types**
   ```python
   # Good practice throughout:
   def calculate_X(
       structure: orm.StructureData,
       energy: orm.Float,
       params: orm.Dict,
   ) -> orm.Dict:
   ```

### ⚠️ Areas for Improvement

1. **Mixed Logging Approaches**
   - Some modules use `logger = get_logger(__name__)`
   - Others use `logger = logging.getLogger(__name__)`
   - Some files have `logger` defined but never used
   - Print statements found in testing module (acceptable for test output)

2. **Inconsistent Parameter Ordering**
   ```python
   # Different parameter order across functions:
   function1(structure, code, parameters, options, ...)
   function2(code, structure, options, parameters, ...)
   ```

   **Recommendation:** Standardize parameter order:
   1. Required AiiDA nodes (structure, code)
   2. Configuration dicts (parameters, options)
   3. Optional overrides
   4. Flags/booleans

3. **Builder Pattern Variations**
   ```python
   # Multiple approaches found:
   # 1. builder_inputs dict (new style)
   # 2. parameters/options dicts (old style)
   # 3. Hybrid approach with priority logic
   ```

   **Recommendation:** Document the migration path and deprecation timeline for old style

4. **Scatter-Gather Naming**
   ```python
   # Inconsistent suffixes:
   relax_slabs_scatter()           # _scatter suffix
   compute_surface_energies_scatter()  # _scatter suffix
   scf_slabs_scatter()             # _scatter suffix

   # vs non-scatter equivalents (single operations)
   calculate_formation_enthalpy()  # no suffix needed
   ```

   **Good practice**: The `_scatter` suffix is clear, but consider documenting this pattern

---

## 5. Readability Assessment

### ✅ Strengths

1. **Clear Variable Names in Complex Logic**
   ```python
   # Example from thermodynamics.py:
   element_M = metal_elements[0]
   element_N_ref = metal_elements[1]
   element_O = 'O'

   bulk_energy_per_fu = bulk_energy.value / formula_units_in_bulk
   ```

2. **Well-Structured Complex Functions**
   - Logical sections with comments
   - Early validation and error handling
   - Clear return value construction

3. **Descriptive Comments for Non-Obvious Logic**
   ```python
   # Normalize to per-atom energies
   energy_per_atom = [e / n_atoms for e in energies] if energies else []

   # Reference energy (highest cutoff, assumed most accurate)
   ref_energy = energy_per_atom[-1] if energy_per_atom else 0.0
   ```

4. **Type Annotations Improve Readability**
   ```python
   def extract_restart_folders_from_node(node_pk: int) -> dict[str, orm.RemoteData]:
       """Clear return type makes usage obvious"""
   ```

### ⚠️ Areas for Improvement

1. **Long Parameter Lists Reduce Readability**
   ```python
   # core_workgraph() has 50+ parameters
   # Consider parameter objects or builder pattern:

   @dataclass
   class WorkGraphConfig:
       bulk_config: BulkConfig
       slab_config: SlabConfig
       computational_config: ComputationalConfig
   ```

2. **Nested Dictionary Access**
   ```python
   # Hard to read:
   misc_dict['total_energies']['energy_extrapolated']

   # Could be improved with helper functions:
   get_extrapolated_energy(misc_dict)
   ```

3. **Complex List Comprehensions**
   ```python
   # Found in convergence module - could be multi-line:
   cutoff_values = [item.get('cutoff', item.get('encut', 0)) for item in data['data']]
   ```

4. **Magic Numbers in Calculations**
   ```python
   # Some calculations have unexplained constants
   # Good: EV_PER_ANGSTROM2_TO_J_PER_M2 (defined in constants.py)
   # Could improve: Document formula derivations inline for complex calculations
   ```

---

## 6. Specific Module Analysis

### workgraph.py (Core Orchestrator)
- **Lines:** 1992
- **Complexity:** High (justified as main entry point)
- **Documentation:** Good module docstring, function docstrings need completion
- **Recommendation:** Consider extracting builder logic to separate module

### utils.py (Utilities)
- **Lines:** 461
- **Organization:** Excellent with clear section comments
- **Documentation:** Comprehensive with examples
- **Recommendation:** Continue this as model for other modules

### thermodynamics.py
- **Lines:** 589
- **Organization:** Clear separation of binary/ternary logic
- **Documentation:** Excellent mathematical descriptions
- **Recommendation:** Exemplary module structure

### adsorption_energy.py
- **Lines:** 1647
- **Complexity:** High, many helper functions
- **Recommendation:** Split into submodule:
  - `adsorption_energy/core.py`
  - `adsorption_energy/builders.py`
  - `adsorption_energy/separation.py`

### constants.py
- **Lines:** 111
- **Organization:** Perfect - categorical grouping with derivations
- **Documentation:** Excellent with references
- **Recommendation:** Use as template for other constant definitions

---

## 7. Priority Recommendations

### High Priority (Impact: Major)

1. **Standardize verb prefixes for functions**
   - `build_*` → WorkGraph constructors
   - `calculate_*` → Single calculations (@task.calcfunction)
   - `compute_*` → Complex workflows
   - `get_*` → Data access

2. **Complete type hints throughout codebase**
   - Add return type hints to all functions
   - Use `from __future__ import annotations` consistently

3. **Document parameter priorities in builder functions**
   - Clear documentation of old-style vs new-style API
   - Migration guide for deprecated patterns

### Medium Priority (Impact: Moderate)

4. **Consolidate duplicate utility functions**
   - Single `load_structure_from_file()` in utils.py
   - Single `get_settings()` function

5. **Add docstrings to all internal functions**
   - Even `_internal()` functions should explain purpose
   - Link to related public API

6. **Organize __all__ exports by category**
   ```python
   __all__ = [
       # Core workflow functions
       "build_core_workgraph",

       # Thermodynamic calculations
       "calculate_formation_enthalpy",
       "calculate_surface_energy_binary",

       # Utilities
       "deep_merge_dicts",
       ...
   ]
   ```

### Low Priority (Impact: Minor)

7. **Refactor long files (>1000 lines)**
   - Extract helpers to separate modules
   - Maintain backward compatibility

8. **Improve magic string handling**
   - Define constants for common dictionary keys
   - Use dataclasses for complex return types

9. **Enhance example documentation**
   - Link examples in docstrings
   - Add "See Also" sections

---

## 8. Code Style Compliance

### PEP 8 Adherence: 8/10
- ✅ Import organization generally correct
- ✅ Naming conventions mostly followed
- ⚠️ Some lines exceed 120 characters (configured in flake8)
- ⚠️ Occasional inconsistent spacing around operators

### Type Hints (PEP 484): 7/10
- ✅ Most functions have parameter type hints
- ⚠️ Return type hints missing in ~20% of functions
- ⚠️ Some complex types could use TypeAlias

### Documentation (PEP 257): 9/10
- ✅ Excellent docstring coverage
- ✅ Google-style format consistent
- ⚠️ Some internal functions lack docstrings

---

## 9. Positive Patterns to Maintain

1. **Comprehensive CLAUDE.md Documentation**
   - Excellent developer onboarding resource
   - Clear architecture explanations
   - Practical troubleshooting guides

2. **Constants Module Design**
   - Physical constants with derivations
   - Clear categorization
   - Referenced CODATA standards

3. **Placeholder Pattern for Optional Outputs**
   ```python
   class TaskOutputPlaceholder(NamedTuple):
       structure: Any
       misc: Optional[Any] = None
   ```
   - Type-safe alternative to dynamic objects
   - Clear intent and documentation

4. **Deep Merge Utility**
   - Essential for parameter overrides
   - Well-documented with examples
   - Prevents common parameter bugs

5. **Preset System Architecture**
   - Three-tier approach (presets → flags → dependencies)
   - Flexible while providing sensible defaults
   - Well-documented use cases

---

## 10. Conclusion

PS-TEROS demonstrates **mature software engineering practices** with well-thought-out architecture, comprehensive documentation, and consistent design patterns. The codebase is generally maintainable and readable.

### Key Strengths
- Excellent module organization with clear separation of concerns
- Comprehensive documentation with practical examples
- Consistent use of AiiDA and WorkGraph patterns
- Strong constant management and physical documentation

### Key Opportunities
- Standardize function naming prefixes
- Complete type hints coverage
- Refactor longest files (>1000 lines)
- Consolidate duplicate utilities

### Overall Assessment
The codebase is **production-ready** with minor improvements needed. The documentation quality and architectural consistency are particularly commendable. With the recommended refinements, this would be an exemplary scientific computing codebase.

**Recommended Actions:**
1. Create style guide document based on existing patterns
2. Add pre-commit hooks for type checking and PEP 8
3. Schedule refactoring sprint for files >1000 lines
4. Document migration path from old-style to new-style API

---

**Review Completed:** 2026-01-19
**Files Analyzed:** 100+ Python files in teros/core
**Lines Reviewed:** ~15,000+ lines of code
