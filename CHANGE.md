# Changelog

## [v1.0.0] - Major Update: AiiDA/WorkGraph Modernization & New Features

### Breaking Changes
- Updated to latest versions of AiiDA and AiiDA-WorkGraph
- Modernized entire codebase for compatibility with current AiiDA ecosystem

### New Features

#### Manual Termination Input
- Added ability to manually input pre-generated slab structures via `input_slabs` parameter
- Bypasses automatic termination generation when custom slabs are provided
- Enables precise control over surface structures and terminations

#### Cleavage Energy Module
- New module `teros/core/cleavage.py` for calculating cleavage energies
- Computes energy required to split crystals into complementary surfaces
- Automatic pairing of complementary terminations following pymatgen convention
- Parallel computation using scatter-gather pattern
- Formula: `Ec(i,j) = 1/(2A) * (E_i^slab + E_j^slab - n*E_bulk)`
- Returns energies in both eV/Ų and J/m²
- Reference: See `examples/cleavage/*.md` for implementation details

#### Default Builders
- Added `teros.default_builders` module for simplified workflow setup
- Pre-configured parameter sets for common oxide systems
- Reduces workflow script complexity by ~200 lines
- Easy parameter override capability
- Examples in `examples/default_builders/`

#### Restart Functionality
- New restart mode using `restart_from_node` parameter
- Accepts PK of previous PS-TEROS calculation
- Reuses RemoteData from incomplete slab relaxations as restart points
- Automatically extracts structures and remote folders from previous run
- Enables continuation of failed or unconverged calculations
- VASP reads WAVECAR and CONTCAR from previous calculation
- Reference: See `examples/restart/*.md` for usage patterns

### Implementation Details
- Modified `teros/core/workgraph.py`: Added restart logic and cleavage integration
- Modified `teros/core/slabs.py`: Added restart folder extraction and collection functions
- Updated all builder functions to support new parameters
- Enhanced documentation throughout core modules

### Backward Compatibility
⚠️ Requires update to latest AiiDA and AiiDA-WorkGraph versions

---

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
