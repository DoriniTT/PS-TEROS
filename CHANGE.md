# Changelog

## [v0.2.0] - 2025-11-04 - Major Feature Release

### Overview

Version 0.2.0 represents a major advancement in PS-TEROS functionality with four flagship features: standalone AIMD module with override system, universal concurrency control, workflow preset system, and user-provided slab structures. This release dramatically improves workflow flexibility, resource management, and ease of use.

---

### New Feature: AIMD Standalone Module with Override System

**Module:** `teros.core.aimd`

**Overview:** Complete standalone AIMD (Ab Initio Molecular Dynamics) module enabling flexible MD simulations on any structures with per-structure, per-stage, and per-combination parameter customization.

**Key Capabilities:**
- **Direct structure input:** Use any AiiDA StructureData or PK
- **Sequential multi-stage AIMD:** Automatic restart chaining between stages (equilibration → production)
- **Three-level override system:** Structure-level, stage-level, and matrix-level INCAR parameter customization
- **Supercell transformations:** Optional supercell creation before AIMD
- **Concurrency control:** Limit parallel VASP calculations with `max_concurrent_jobs`
- **Complete provenance:** Full AiiDA tracking of all calculations

**Three-Level Override System:**
1. **Structure-level:** Apply to all stages of specific structures
2. **Stage-level:** Apply to all structures in specific stages
3. **Matrix-level:** Apply to specific (structure, stage) combinations
- **Priority order:** matrix > stage > structure > base

**API:**
```python
from teros.core.aimd import build_aimd_workgraph

wg = build_aimd_workgraph(
    structures={'slab1': structure1, 'slab2': structure2},
    aimd_stages=[
        {'temperature': 300, 'steps': 100},   # Equilibration
        {'temperature': 300, 'steps': 500},   # Production
    ],
    code_label='VASP-6.5.1@cluster02',
    builder_inputs=base_inputs,

    # Override system
    structure_overrides={'slab2': {'parameters': {'incar': {'ENCUT': 500}}}},
    stage_overrides={1: {'parameters': {'incar': {'EDIFF': 1e-7}}}},
    matrix_overrides={('slab1', 1): {'parameters': {'incar': {'ALGO': 'All'}}}},

    max_concurrent_jobs=2,
)
```

**New Files:**
- `teros/core/aimd/` - Complete module package
- `docs/aimd_standalone_module.md` - Comprehensive user guide (19KB)
- `examples/vasp/step_18_aimd_standalone.py` - Basic example
- `examples/vasp/step_19_aimd_with_overrides.py` - Override system demo

**Use Cases:**
- Temperature series (100K → 300K → 500K)
- Different convergence requirements per structure
- Multi-phase equilibration and production runs
- Custom MD protocols with stage-specific precision

**Documentation:** `docs/aimd_standalone_module.md`

---

### New Feature: Universal Concurrency Control

**Feature:** `max_concurrent_jobs` parameter for all PS-TEROS modules

**Version:** v2.2.0 (part of v0.2.0 release)

**Overview:** Extended concurrency control to 100% of PS-TEROS workflow modules, providing unified resource management across all calculation types.

**Complete Module Coverage:**
1. ✅ Core workflows (`build_core_workgraph`, `build_core_workgraph_with_map`)
2. ✅ AIMD modules (VASP and CP2K)
3. ✅ Surface hydroxylation workflows
4. ✅ Custom calculations
5. ✅ All scatter functions (slabs, electronic properties, adsorption)

**Key Implementation:**
- **Nested workgraph support** via `get_current_graph()` API
- **Automatic propagation** through all nesting levels
- **Unified API** across all modules
- **Production ready** with verified behavior

**Usage:**
```python
# Core workflow
wg = build_core_workgraph(
    max_concurrent_jobs=4,  # Max 4 VASP jobs at once
    # ... parameters
)

# AIMD workflow
wg = build_aimd_workgraph(
    max_concurrent_jobs=2,  # Serial-like execution
    # ... parameters
)

# Hydroxylation workflow
wg = build_surface_hydroxylation_workgraph(
    max_parallel=10,         # Process first 10 structures
    max_concurrent_jobs=3,   # Run 3 at a time
    # ... parameters
)
```

**Benefits:**
- **Resource control:** Match cluster capacity (e.g., 24 cores / 6 cores per job = 4 max)
- **Serial execution:** Set to 1 for debugging or minimal resources
- **Queue optimization:** Control job submission rate
- **Cluster-agnostic:** Works with any scheduler/queue system

**Common Values:**
- `1` - Serial mode (debugging, minimal resources)
- `4` - Default (balanced, small clusters)
- `8+` - Higher concurrency (medium/large clusters)
- `None` - Unlimited (queue-managed systems)

**Files Modified:**
- `teros/core/aimd.py` - VASP AIMD
- `teros/core/aimd_cp2k.py` - CP2K AIMD
- `teros/core/surface_hydroxylation/relaxations.py` - Hydroxylation
- `teros/core/custom_calculation/workgraph.py` - Custom calculations
- `teros/core/workgraph.py` - Propagation to AIMD stages

**Documentation:**
- `docs/CONCURRENCY_CONTROL.md` - Complete feature guide
- `docs/MAX_CONCURRENT_JOBS_ALL_MODULES.md` - Implementation summary

---

### New Feature: Workflow Preset System

**Overview:** Three-tier workflow configuration system providing high-level convenience presets with fine-grained override capability.

**Three Tiers:**
1. **Named workflow presets:** One parameter activates entire workflows
2. **Independent component flags:** Override preset defaults as needed
3. **Automatic dependency resolution:** Smart validation and auto-enabling

**11 Available Presets:**
1. `surface_thermodynamics` - Complete surface thermodynamics (default)
2. `surface_thermodynamics_unrelaxed` - Quick screening with unrelaxed slabs
3. `cleavage_only` - Cleavage energy calculations
4. `relaxation_energy_only` - Relaxation energy analysis
5. `bulk_only` - Bulk optimization only
6. `formation_enthalpy_only` - Formation enthalpy calculation
7. `electronic_structure_bulk_only` - DOS/bands for bulk
8. `electronic_structure_slabs_only` - DOS/bands for slabs
9. `electronic_structure_bulk_and_slabs` - DOS/bands for both
10. `aimd_only` - Molecular dynamics on slabs
11. `comprehensive` - Everything enabled

**API Changes:**
```python
# New simplified API
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # One line!
    structures_dir='structures',
    bulk_name='ag2o.cif',
    # ... rest (all flags set automatically)
)

# With overrides
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=False,  # Override: disable cleavage
    # ... rest
)
```

**Key Updates to `surface_thermodynamics` Preset:**
- Cleavage energies now **optional** (disabled by default)
- Relaxation energies now **optional** (disabled by default)
- Migration: Add `compute_cleavage=True` and/or `compute_relaxation_energy=True` to enable

**New Helper Functions:**
- `list_workflow_presets()` - List all available presets
- `get_preset_summary(name)` - Detailed preset configuration
- `resolve_preset(name, **overrides)` - Resolve final configuration

**Benefits:**
- **Simplicity:** One parameter replaces 7+ boolean flags
- **Validation:** Automatic parameter requirement checking
- **Flexibility:** Override any default as needed
- **Safety:** Clear error messages for missing requirements

**Documentation:**
- `docs/WORKFLOW_PRESETS_GUIDE.md` - Complete user guide (15KB)
- `docs/WORKFLOW_PRESETS_EXAMPLES.md` - Example workflows (19KB)
- `docs/WORKFLOW_SYSTEM_EXPLAINED.md` - Architecture overview (10KB)
- `docs/WORKFLOW_MIGRATION_GUIDE.md` - Migration from old API (11KB)

**Backward Compatibility:** ✅ Fully backward compatible

---

### New Feature: Surface Hydroxylation Module

**Module:** `teros.core.surface_hydroxylation`

**Overview:** Automated generation and relaxation of surface variants with different hydroxylation coverages and oxygen vacancy configurations for studying surface chemistry and catalytic activity.

**Key Capabilities:**
- **Hydroxylation:** Add OH groups to surface oxygen atoms
- **Oxygen vacancies:** Remove oxygen atoms to create defects
- **Combined mode:** Both hydroxylation and vacancies
- **Coverage-based deduplication:** Reduces thousands of structures to ~10 representative configurations
- **Batch VASP relaxation:** Parallel relaxation with concurrency control
- **Complete provenance:** Full AiiDA tracking of all calculations

**Scientific Background:**

Surface hydroxylation occurs when water dissociatively adsorbs on oxide surfaces, converting surface oxygen (O) to hydroxyl groups (OH). This process is fundamental in:
- Catalysis (different activity when hydroxylated)
- Water splitting (photocatalytic intermediates)
- Surface stability (stabilizes unstable terminations)
- Electronic properties (modifies work function)

**Coverage-Based Deduplication Algorithm:**

For surfaces with many oxygen atoms, the number of possible configurations grows exponentially:
- 10 oxygen atoms → 1,024 combinations
- 20 oxygen atoms → 1,048,576 combinations

The deduplication algorithm:
1. Generates all unique configurations
2. Groups by coverage (OH/nm² or vacancies/nm²)
3. Bins into N equal-width coverage ranges
4. Samples representative structures from each bin

**Result:** Reduces thousands to ~N structures spanning 0-100% coverage.

**API:**
```python
from teros.core.surface_hydroxylation import (
    build_surface_hydroxylation_workgraph,
    organize_hydroxylation_results,
)

# Define surface parameters
surface_params = {
    'mode': 'hydrogen',              # 'hydrogen', 'vacancies', or 'combine'
    'species': 'O',                  # Target oxygen atoms
    'z_window': 0.5,                 # Surface detection window (Å)
    'which_surface': 'top',          # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,                 # O-H bond distance (Å)
    'coverage_bins': 5,              # Number of coverage bins
    'deduplicate_by_coverage': True, # Enable deduplication
}

# Build and submit workflow
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,               # Input slab PK
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,   # VASP configuration
    max_concurrent_jobs=3,           # Concurrency control
    fix_type='bottom',               # Optional: fix bottom atoms
    fix_thickness=5.0,               # Optional: fix bottom 5Å
)
result = wg.submit()

# After completion, organize results
node = orm.load_node(result.pk)
results = organize_hydroxylation_results(node)
```

**Three Modes:**

1. **Hydrogen mode:** Add OH groups to surface oxygen atoms
   - Controls: `oh_dist` (O-H bond distance), `coverage_bins`

2. **Vacancies mode:** Remove surface oxygen atoms
   - Creates defects and active sites

3. **Combine mode:** Both hydroxylation and vacancies
   - Studies complex surface chemistry scenarios

**Advanced Features:**

- **Selective dynamics:** Fix bottom atoms during relaxation
  - `fix_type='bottom'` - Fix bottom layer
  - `fix_thickness=5.0` - Fix bottom 5Å

- **Per-structure builder overrides:**
  ```python
  structure_specific_builder_inputs={
      0: {'parameters': {'incar': {'ALGO': 'Normal'}}},
      2: {'kpoints_spacing': 0.2},
  }
  ```

- **Batch parallelization:** `max_parallel` limits number of structures to process
- **Concurrency control:** `max_concurrent_jobs` limits simultaneous VASP jobs

**Output Organization:**

Results use descriptive naming: `{index}_{variant_name}`
- `0_oh_000_3_7572` - First structure, 3.76 OH/nm²
- `1_oh_001_7_5145` - Second structure, 7.51 OH/nm²
- `2_vac_000_5_2341` - Third structure, 5.23 vacancies/nm²

**New Files:**
- `teros/core/surface_hydroxylation/` - Complete module package
- `docs/surface_hydroxylation.md` - Comprehensive user guide (30KB)
- `examples/hydroxylation/` - Example workflows

**Use Cases:**
- Catalysis studies (hydroxylated vs pristine surfaces)
- Water splitting intermediates
- Defect chemistry (oxygen vacancy effects)
- Surface reconstruction analysis
- Coverage-dependent property studies

**Integration:**

Works seamlessly with other PS-TEROS modules:
- Generate slabs with main workflow
- Apply hydroxylation/vacancies
- Calculate surface energies with thermodynamics module
- Run AIMD on hydroxylated surfaces

**Documentation:** `docs/surface_hydroxylation.md`

---

### New Feature: User-Provided Slab Structures

**Overview:** Ability to provide pre-generated slab structures as input to PS-TEROS workflows.

**Benefits:**
- **Flexibility:** Use any slab generation method or tool
- **Control:** Exact control over surface structures and terminations
- **Reproducibility:** Use exact structures from literature or previous work
- **Efficiency:** Skip generation when structures already available
- **Compatibility:** Works with all ASE-supported formats (CIF, POSCAR, xyz, etc.)

**API Changes:**
```python
# Load pre-generated slabs
slabs = {
    'term_0': orm.StructureData(ase=read('slab_100.cif')),
    'term_1': orm.StructureData(ase=read('slab_110.cif')),
}

# Build workflow with user slabs
wg = build_core_workgraph(
    input_slabs=slabs,  # Bypasses automatic generation
    structures_dir='structures',
    bulk_name='ag2o.cif',
    # miller_indices NOT required when using input_slabs
    # ... rest of parameters
)
```

**Modified Functions:**
- `core_workgraph()` - Added `input_slabs` parameter
- `build_core_workgraph()` - Added `input_slabs` parameter
- `build_core_workgraph_with_map()` - Added `input_slabs` parameter

**Use Cases:**
- Reproducing published slab structures
- Custom surface reconstructions
- Adding adsorbates or defects to surfaces
- Using manually edited structures
- Integration with external tools
- Testing specific terminations

**New Files:**
- `examples/slabs/slabs_input_relax.py` - Complete example
- `examples/slabs/compare_modes.py` - Comparison demo
- `examples/slabs/QUICKSTART.md` - 5-minute quick start
- `docs/USER_PROVIDED_SLABS.md` - Technical documentation (6KB)

**Backward Compatibility:** ✅ Fully backward compatible

---

### All Features Summary

**Version 0.2.0 includes:**

1. ✅ **AIMD Standalone Module** - Flexible MD with override system
2. ✅ **Universal Concurrency Control** - 100% module coverage
3. ✅ **Workflow Preset System** - 11 presets, simplified API
4. ✅ **Surface Hydroxylation Module** - Automated OH/vacancy generation with deduplication
5. ✅ **User-Provided Slabs** - Custom structure input
6. ✅ **Enhanced Documentation** - 9 new/updated docs (>110KB)

**Files Modified:** 10 core modules + 6 examples
**Documentation Added/Updated:** 9 files
**Backward Compatibility:** ✅ 100% - All existing scripts work unchanged

---

### Migration Guide

**From v0.1.x to v0.2.0:**

No breaking changes. All v0.1.x scripts work unchanged.

**Recommended updates for new features:**

```python
# OLD (still works)
wg = build_core_workgraph(
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=True,
    # ... 7+ boolean flags
)

# NEW (recommended)
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # Cleaner!
    compute_cleavage=True,                     # Enable optional features
    compute_relaxation_energy=True,
    max_concurrent_jobs=4,                     # Resource control
)
```

**To use new AIMD module:**
```python
from teros.core.aimd import build_aimd_workgraph

wg = build_aimd_workgraph(
    structures={'slab': structure},
    aimd_stages=[{'temperature': 300, 'steps': 500}],
    code_label='VASP-6.5.1@cluster',
    builder_inputs=inputs,
)
```

---

### Documentation

**New Documentation (9 files, >110KB):**
1. `docs/aimd_standalone_module.md` - AIMD module guide (19KB)
2. `docs/surface_hydroxylation.md` - Hydroxylation module guide (30KB)
3. `docs/CONCURRENCY_CONTROL.md` - Concurrency control guide (7KB)
4. `docs/MAX_CONCURRENT_JOBS_ALL_MODULES.md` - Implementation summary (10KB)
5. `docs/WORKFLOW_PRESETS_GUIDE.md` - Preset system guide (15KB)
6. `docs/WORKFLOW_PRESETS_EXAMPLES.md` - Preset examples (19KB)
7. `docs/WORKFLOW_SYSTEM_EXPLAINED.md` - Architecture overview (10KB)
8. `docs/WORKFLOW_MIGRATION_GUIDE.md` - Migration guide (11KB)
9. `docs/USER_PROVIDED_SLABS.md` - User slabs guide (6KB)

**Updated Documentation:**
- `docs/README.md` - Updated structure documentation

---

### Testing

**New Examples:**
- `examples/vasp/step_18_aimd_standalone.py` - Basic AIMD
- `examples/vasp/step_19_aimd_with_overrides.py` - Override system
- `examples/hydroxylation/` - Hydroxylation workflow examples
- `examples/slabs/slabs_input_relax.py` - User-provided slabs
- `examples/slabs/compare_modes.py` - Generation vs input comparison

**Verification:**
- All modules tested with `max_concurrent_jobs`
- Override system verified with multiple structures
- Preset system validated across all 11 presets
- Hydroxylation module tested with all three modes
- Coverage-based deduplication verified with large surfaces
- User slab input tested with CIF, POSCAR formats

---

### Known Issues

None. All features production-ready.

---

### Contributors

- Thiago T. Dorini - All features

---

## [Unreleased] - Workflow Preset System Updates

### Updated: Workflow Presets (October 2025)

**Summary:** Improved workflow preset system with optional components and new electronic structure presets.

**Modified Presets:**

1. **`surface_thermodynamics` preset:**
   - Cleavage energies now **optional** (disabled by default)
   - Relaxation energies now **optional** (disabled by default)
   - **Migration:** Add `compute_cleavage=True` and/or `compute_relaxation_energy=True` if you need these features

2. **`aimd_only` preset:**
   - Slab relaxation now **disabled** by default
   - AIMD runs on unrelaxed slabs (generated from bulk)
   - **Migration:** Add `relax_slabs=True` if you need slabs relaxed before AIMD

**New Presets:**

3. **`electronic_structure_slabs_only`:**
   - Electronic properties (DOS and band structure) for slabs only
   - Requires: `slab_bands_parameters`, `slab_band_settings`

4. **`electronic_structure_bulk_and_slabs`:**
   - Electronic properties for both bulk and slabs
   - Requires: `bands_parameters`, `band_settings`, `slab_bands_parameters`, `slab_band_settings`

**New Documentation:**

- `docs/WORKFLOW_SYSTEM_EXPLAINED.md` - Comprehensive three-tier system guide
- `docs/WORKFLOW_MIGRATION_GUIDE.md` - Migration guide for existing scripts
- Updated `docs/WORKFLOW_PRESETS_GUIDE.md` with accurate preset definitions
- Updated `docs/WORKFLOW_PRESETS_EXAMPLES.md` with new examples

**Backward Compatibility:**

✅ **Fully backward compatible**
- All existing presets still work
- Old scripts will run unchanged (may need explicit flags for optional features)
- No breaking changes to API

**See also:**
- [Workflow Updates Summary](examples/step_by_step/WORKFLOW_UPDATES_SUMMARY.md)
- [Migration Guide](docs/WORKFLOW_MIGRATION_GUIDE.md)

---

## [Unreleased] - Relaxation Energy Module

### New Feature: Relaxation Energy Calculation

**Overview**: Added optional calculation of relaxation energies for slab terminations, quantifying the energetic stabilization from atomic relaxation at surfaces.

**Formula**: `E_relax = E_relaxed - E_unrelaxed`

**Implementation**:
- Three separate WorkGraphs created when `compute_relaxation_energy=True`:
  1. `scf_slabs_scatter` - SCF calculations on unrelaxed slabs (NSW=0, IBRION=-1)
  2. `relax_slabs_scatter` - Slab relaxation (existing, unchanged)
  3. `calculate_relaxation_energies_scatter` - Energy difference calculation
- All workgraphs use scatter-gather pattern for parallel execution
- Fully optional: controlled by `compute_relaxation_energy` parameter (default: False)
- **Backward compatible**: Default behavior unchanged

### API Changes

**New Parameter**:
```python
build_core_workgraph(
    # ... existing parameters ...
    compute_relaxation_energy=False,  # NEW: Enable relaxation energy calculation
)
```

**New Outputs** (when `compute_relaxation_energy=True`):
- `unrelaxed_slab_energies`: Total energies from SCF calculations (NSW=0, IBRION=-1)
- `unrelaxed_slab_remote`: RemoteData nodes from SCF calculations
- `relaxation_energies`: E_relaxed - E_unrelaxed for each slab termination

### Modified Files

**Core Implementation**:
- `teros/core/slabs.py`:
  - Added `scf_slabs_scatter()` - SCF workgraph for unrelaxed slabs
  - Added `calculate_relaxation_energies_scatter()` - Energy calculation workgraph
  - Added `calculate_energy_difference()` - Helper calcfunction
  - Added `scf_relax_and_calculate_relaxation_energy()` - Combined function (alternative approach)

- `teros/core/workgraph.py`:
  - Added `compute_relaxation_energy` parameter to `core_workgraph()`
  - Added conditional logic to create SCF and relaxation energy workgraphs
  - Updated `build_core_workgraph()` to pass new parameter
  - Enhanced docstrings with new parameter documentation

### New Files

**Examples**:
- `examples/slabs/ag2o_100_relaxation_energy.py`: Complete working example
- `examples/slabs/monitor_relaxation_energy.py`: Monitoring script
- `examples/slabs/README_RELAXATION_ENERGY.md`: Example documentation

**Documentation**:
- `docs/RELAXATION_ENERGY.md`: Comprehensive user guide and API reference
- `QUICKSTART_RELAXATION_ENERGY.md`: Quick start guide
- `RELAXATION_ENERGY_IMPLEMENTATION.md`: Technical implementation details
- `IMPLEMENTATION_SUCCESS.md`: Verification and test results

### Usage Example

```python
wg = build_core_workgraph(
    structures_dir="/path/to/structures",
    bulk_name="oxide.cif",
    metal_name="metal.cif",
    oxygen_name="O2.cif",
    # ... other parameters ...
    relax_slabs=True,                    # Enable slab relaxation
    compute_relaxation_energy=True,      # Enable relaxation energy (OPTIONAL)
)
```
### Backward Compatibility

✅ **Fully backward compatible**
- Default `compute_relaxation_energy=False` preserves existing behavior
- No changes required to existing scripts
- SCF calculations only performed when explicitly requested

---

## [Unreleased] - Electronic Properties Module (DOS & Band Structure)

### New Feature: Electronic Structure Calculations

**Overview**: Added comprehensive electronic properties calculations for both bulk and slab structures, enabling density of states (DOS) and band structure analysis using the vasp.v2.bands workchain with seekpath integration.

**Capabilities**:
- **Bulk electronic properties**: DOS and band structure for relaxed bulk structures
- **Slab electronic properties**: DOS and band structure for selected slab terminations
- **Automatic high-symmetry path generation**: Using seekpath-aiida or pymatgen
- **Three-stage workflow**: SCF → Band structure → DOS
- **Material-agnostic**: Works for any system (binary, ternary, metals, semiconductors, insulators)
- **Configurable k-point sampling**: Independent control for SCF, bands, and DOS meshes

### Implementation

**Core Module**:
- `teros/core/builders/electronic_properties_builder.py`:
  - `get_electronic_properties_defaults()`: Builder for bulk structures
  - `get_slab_electronic_properties_defaults()`: Builder for slab structures (denser k-point sampling)
  - Returns complete parameter sets for vasp.v2.bands workchain

**Calculation Stages**:
1. **SCF stage**: Self-consistent calculation with LWAVE=True, LCHARG=True
2. **Band structure stage**: Non-self-consistent (ICHARG=11) along high-symmetry paths
3. **DOS stage**: Non-self-consistent with dense k-point mesh

### API Changes

**New Parameters for Bulk Electronic Properties**:
```python
build_core_workgraph(
    # ... existing parameters ...
    compute_electronic_properties_bulk=False,  # Enable bulk DOS/bands
    bands_parameters=None,                     # INCAR parameters (from builder)
    band_settings=None,                        # Band workflow settings
    bands_options=None,                        # Scheduler options
)
```

**New Parameters for Slab Electronic Properties**:
```python
build_core_workgraph(
    # ... existing parameters ...
    compute_electronic_properties_slabs=False,      # Enable slab DOS/bands
    slab_electronic_properties=None,                # Per-slab configuration dict
    slab_bands_parameters=None,                     # Default INCAR parameters
    slab_band_settings=None,                        # Default band settings
    slab_bands_options=None,                        # Default scheduler options
)
```

**New Outputs**:

For bulk (when `compute_electronic_properties_bulk=True`):
- `bulk_bands`: Band structure along high-symmetry paths
- `bulk_dos`: Density of states
- `bulk_primitive_structure`: Primitive cell used for band calculation
- `bulk_seekpath_parameters`: Seekpath symmetry information and k-point paths

For slabs (when `compute_electronic_properties_slabs=True`):
- `slab_bands`: Dictionary of band structures (one per selected slab)
- `slab_dos`: Dictionary of DOS data (one per selected slab)
- `slab_primitive_structures`: Dictionary of primitive cells
- `slab_seekpath_parameters`: Dictionary of seekpath info

### Configuration Examples

**Bulk Electronic Properties**:
```python
from teros.core.builders import get_electronic_properties_defaults

# Get default configuration
ep_defaults = get_electronic_properties_defaults(
    energy_cutoff=520,
    electronic_convergence=1e-5,
    ncore=4,
    ispin=2,  # Spin-polarized
    kpoints_mesh_density=0.3,     # SCF mesh
    band_kpoints_distance=0.2,    # Band path density
    dos_kpoints_distance=0.2,     # DOS mesh density
    line_density=0.2,             # Points along paths
    nedos=2000,                   # DOS grid points
    band_mode="seekpath-aiida",   # Auto high-symmetry paths
)

wg = build_core_workgraph(
    compute_electronic_properties_bulk=True,
    bands_parameters=ep_defaults,
    band_settings=ep_defaults['band_settings'],
    bands_options=bulk_options,
    # ... other parameters ...
)
```

**Slab Electronic Properties**:
```python
from teros.core.builders import get_slab_electronic_properties_defaults

# Get slab-tuned configuration (denser k-points for 2D systems)
slab_ep_defaults = get_slab_electronic_properties_defaults(
    energy_cutoff=520,
    kpoints_mesh_density=0.25,    # Denser than bulk
    band_kpoints_distance=0.15,   # Denser path sampling
    line_density=0.15,            # More points along paths
)

# Define which slabs to calculate
slab_electronic_properties = {
    'term_0': {
        'bands_parameters': slab_ep_defaults,
        'bands_options': slab_options,
        'band_settings': slab_ep_defaults['band_settings'],
    },
    'term_1': {
        'bands_parameters': slab_ep_defaults,
        'bands_options': slab_options,
        'band_settings': slab_ep_defaults['band_settings'],
    },
}

wg = build_core_workgraph(
    compute_electronic_properties_slabs=True,
    slab_electronic_properties=slab_electronic_properties,
    slab_bands_parameters=slab_ep_defaults,
    slab_band_settings=slab_ep_defaults['band_settings'],
    slab_bands_options=slab_options,
    # ... other parameters ...
)
```

### Modified Files

**Core Implementation**:
- `teros/core/workgraph.py`:
  - Added electronic properties integration for both bulk and slabs
  - Added conditional logic to create vasp.v2.bands workgraphs
  - Updated `build_core_workgraph()` with new parameters

### New Files

**Core Modules**:
- `teros/core/builders/electronic_properties_builder.py`: Material-agnostic parameter builders

**Examples**:
- `examples/electronic_properties/bulk_dos_bands_ag2o.py`: Bulk DOS/bands example
- `examples/electronic_properties/slab_electronic_properties_example.py`: Slab DOS/bands example
- `examples/complete/complete_ag2o_example.py`: Updated with electronic properties
- `examples/complete/complete_ag3po4_example.py`: Updated with both bulk and slab electronic properties

### Key Features

**Flexibility**:
- Independent control over SCF, band, and DOS k-point densities
- Per-slab configuration for heterogeneous terminations
- Material-agnostic builders work for any system

**Integration**:
- Seamless integration with vasp.v2.bands workchain
- Automatic primitive cell detection and standardization
- Seekpath integration for high-symmetry paths

**Performance**:
- Three-stage workflow with WAVECAR/CHGCAR reuse
- Parallel calculation across multiple slabs
- Configurable computational resources per stage

### Backward Compatibility

✅ **Fully backward compatible**
- Default `compute_electronic_properties_bulk=False` and `compute_electronic_properties_slabs=False`
- No changes required to existing scripts
- Electronic properties calculations only performed when explicitly requested

---

## [Unreleased] - AIMD Module (Ab Initio Molecular Dynamics)

### New Feature: AIMD on Slab Structures

**Overview**: Added comprehensive ab initio molecular dynamics (AIMD) module for running sequential MD simulations on slab structures with automatic restart chaining between stages.

**Capabilities**:
- **Parallel AIMD**: Run AIMD on multiple slab terminations simultaneously
- **Sequential stages**: Chain multiple AIMD stages (equilibration → production)
- **Automatic restart**: Each stage uses output of previous stage as restart point
- **Isothermal runs**: TEBEG=TEEND for constant temperature MD
- **Flexible staging**: Independent control of temperature and timesteps per stage
- **Trajectory output**: Full MD trajectory and energies from each stage

### Implementation

**Core Module**:
- `teros/core/aimd.py`:
  - `aimd_single_stage_scatter()`: Run single AIMD stage on all slabs in parallel
  - `prepare_aimd_parameters()`: Helper to inject stage-specific parameters (temperature, NSW)
  - `get_settings()`: Parser settings for trajectory and structure output

**Workflow Pattern**:
```
Stage 1: Equilibration (300K, 100 steps) on all slabs in parallel
    ↓ (automatic restart via remote_folders)
Stage 2: Production (300K, 500 steps) on all slabs in parallel
    ↓
Final: Trajectories, structures, energies for all slabs
```

### API

**Function Signature**:
```python
aimd_single_stage_scatter(
    slabs: dict[str, StructureData],          # Slabs to run AIMD on
    temperature: float,                        # Target temperature (K)
    steps: int,                                # Number of MD steps (NSW)
    code: Code,                                # VASP code
    aimd_parameters: dict,                     # Base AIMD INCAR (IBRION=0, MDALGO, etc.)
    potential_family: str,                     # Potential family
    potential_mapping: dict,                   # Element mapping
    options: dict,                             # Scheduler options
    kpoints_spacing: float,                    # K-point spacing
    clean_workdir: bool,                       # Clean work directory
    restart_folders: dict[str, RemoteData] = {},  # Optional restart from previous stage
) -> dict
```

**Returns**:
```python
{
    'structures': {slab_label: StructureData},    # Final structures
    'remote_folders': {slab_label: RemoteData},   # For next stage restart
    'energies': {slab_label: Float},              # Total energies
}
```

### Usage Example

**Single-Stage AIMD**:
```python
from teros.core.aimd import aimd_single_stage_scatter

# Base AIMD parameters (constant for all stages)
aimd_parameters = {
    'IBRION': 0,       # Molecular dynamics
    'MDALGO': 2,       # Nosé-Hoover thermostat
    'POTIM': 2.0,      # Timestep (fs)
    'SMASS': 3.0,      # Nosé mass
    'PREC': 'Normal',
    'ENCUT': 400,
    'ISMEAR': 0,
    'SIGMA': 0.1,
    'LWAVE': True,
    'LCHARG': True,
}

# Run AIMD stage
stage1 = aimd_single_stage_scatter(
    slabs=relaxed_slabs,         # From slab relaxation
    temperature=300,              # 300 K
    steps=100,                    # 100 MD steps
    code=code,
    aimd_parameters=aimd_parameters,
    potential_family='PBE',
    potential_mapping={'Ag': 'Ag', 'O': 'O'},
    options=options,
    kpoints_spacing=0.5,
    clean_workdir=False,
)
```

**Multi-Stage Sequential AIMD**:
```python
# Stage 1: Equilibration
stage1 = aimd_single_stage_scatter(
    slabs=relaxed_slabs,
    temperature=300,
    steps=100,                   # Equilibration: 100 steps
    # ... other parameters ...
)

# Stage 2: Production (automatic restart from Stage 1)
stage2 = aimd_single_stage_scatter(
    slabs=stage1['structures'],          # Use output structures
    temperature=300,
    steps=500,                           # Production: 500 steps
    restart_folders=stage1['remote_folders'],  # Restart from Stage 1 WAVECAR
    # ... other parameters ...
)

# Stage 3: Extended production (restart from Stage 2)
stage3 = aimd_single_stage_scatter(
    slabs=stage2['structures'],
    temperature=300,
    steps=1000,
    restart_folders=stage2['remote_folders'],
    # ... other parameters ...
)
```

**Temperature Ramping**:
```python
# Ramp from 100K to 300K
temps = [100, 150, 200, 250, 300]
stages = []
current_slabs = initial_slabs
current_folders = {}

for temp in temps:
    stage = aimd_single_stage_scatter(
        slabs=current_slabs,
        temperature=temp,
        steps=50,                        # Short equilibration at each temp
        restart_folders=current_folders,
        # ... other parameters ...
    )
    current_slabs = stage['structures']
    current_folders = stage['remote_folders']
    stages.append(stage)
```

### Modified Files

**Core Implementation**:
- `teros/core/aimd.py`: New module with AIMD functionality

### New Files

**Examples**:
- `examples/complete/complete_ag2o_aimd_example.py`: AIMD workflow for Ag₂O slabs
- `examples/complete/complete_ag3po4_aimd_example.py`: AIMD workflow for Ag₃PO₄ slabs

### Key Features

**Parallel Execution**:
- Scatter-gather pattern runs AIMD on all slabs simultaneously
- Each slab runs independently with identical MD parameters

**Restart Chaining**:
- Automatic restart via `restart_folders` parameter
- VASP reads WAVECAR, CHGCAR from previous stage
- Seamless multi-stage workflows

**Flexibility**:
- Independent control of temperature and steps per stage
- Easy temperature ramping and multi-phase protocols
- Compatible with any MDALGO (Nosé-Hoover, Langevin, etc.)

**Output**:
- Full MD trajectory for each slab (via parser_settings)
- Final structures for analysis or further stages
- Total energies for convergence monitoring

### Typical Workflow

1. **Slab relaxation** → Get equilibrium structures
2. **AIMD equilibration** → 100-200 steps at target temperature
3. **AIMD production** → 500-1000+ steps with restart from equilibration
4. **Analysis** → Extract trajectory, compute RDF, MSD, etc.

### Backward Compatibility

✅ **Fully backward compatible**
- AIMD is a separate module, does not affect existing workflows
- No changes to `build_core_workgraph()` required
- Users explicitly call `aimd_single_stage_scatter()` when needed

---

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
