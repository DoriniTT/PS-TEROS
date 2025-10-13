# Changelog

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
