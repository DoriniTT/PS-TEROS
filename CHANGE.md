# Changelog

## [v0.2.0] - User-Provided Slab Structures Feature

### New Feature: Input Slabs Support
- Added ability to provide pre-generated slab structures as input to PS-TEROS workflows
- Gives users full flexibility to use custom slab structures instead of automatic generation
- Maintains full backward compatibility with existing code

### API Changes
- Added `input_slabs` parameter to core workflow functions:
  - `core_workgraph()` in `teros/core/workgraph.py`
  - `build_core_workgraph()`
  - `build_core_workgraph_with_map()`
- Made slab generation parameters optional when `input_slabs` is provided:
  - `miller_indices`, `min_slab_thickness`, `min_vacuum_thickness` no longer required
  - Generation parameters ignored when user provides slabs

### Implementation Details
- Modified `teros/core/workgraph.py`:
  - Added conditional logic to use provided slabs or generate them automatically
  - Updated function signatures to accept `input_slabs` parameter
  - Enhanced docstrings with new parameter documentation
  - Added validation to ensure generation parameters are provided when needed

### New Examples
- `examples/slabs/slabs_input_relax.py`: Complete working example using user-provided slabs
- `examples/slabs/compare_modes.py`: Comparison demonstration of both modes
- `examples/slabs/QUICKSTART.md`: 5-minute quick start guide
- `examples/slabs/input_structures/`: Directory for user slab structure files

### New Documentation
- `docs/USER_PROVIDED_SLABS.md`: Technical documentation and API changes
- `examples/slabs/README_INPUT_SLABS.md`: Comprehensive user guide
- `examples/slabs/input_structures/README.md`: Instructions for input files

### Benefits
- **Flexibility**: Use any slab generation method or tool
- **Control**: Exact control over surface structures and terminations
- **Reproducibility**: Use exact structures from literature or previous work
- **Efficiency**: Skip generation when structures are already available
- **Compatibility**: Works with all ASE-supported file formats (CIF, POSCAR, xyz, etc.)

### Use Cases
- Reproducing published slab structures
- Custom surface reconstructions
- Adding adsorbates or defects to surfaces
- Using manually edited or specialized slab configurations
- Integration with external slab generation tools
- Testing specific surface terminations

### Backward Compatibility
✅ Fully backward compatible - all existing scripts work without modification

---

## [v0.1.2] - Enthalpy of Formation Bug Fix & Oxygen Chemical Potential Limits

### Bug Fixes
- Fixed calculation of enthalpy of formation to correct numerical errors in previous versions.
- Corrected the limits for the oxygen chemical potential to ensure physically meaningful bounds are enforced.

---

## [v0.1.1] - Added DEFECT_TYPES File Generation

### New Feature
- Added functionality to analyze stoichiometric deviations in surface slab terminations
- Generated DEFECT_TYPES file documents excess/deficit elements in each termination
- Enables subsequent calculations with charge compensations for specific terminations

### Implementation
- Modified `get_slabs` function in `functions/slabs.py`:
  - Added `analyze_defects` helper function to compare slab/bulk compositions
  - Calculates element excess relative to bulk stoichiometry for each termination
  - Outputs results to DEFECT_TYPES with columns: termination, [elements]

### File Format
DEFECT_TYPES contains:
- termination: Slab identifier (s_0, s_1, etc.)
- One column per element in bulk structure
