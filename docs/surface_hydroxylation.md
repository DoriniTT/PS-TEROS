# Surface Hydroxylation Module

## Overview

The surface hydroxylation module (`teros.core.surface_hydroxylation`) enables automated generation and relaxation of surface variants with different hydroxylation coverages and oxygen vacancy configurations. This is essential for studying surface chemistry, catalytic activity, and surface stability under reactive conditions.

**Key capabilities:**
- Generate hydroxylated surface structures (OH groups)
- Generate oxygen-deficient surfaces (vacancies)
- Combined hydroxylation and vacancy modes
- Coverage-based deduplication (reduces thousands to ~10 representative structures)
- Batch VASP relaxation with configurable parallelization
- Complete AiiDA provenance tracking

## Scientific Background

### Hydroxylation

Surface hydroxylation occurs when water dissociatively adsorbs on oxide surfaces, converting surface oxygen atoms (O) to hydroxyl groups (OH). This process is fundamental in:

- **Catalysis**: Many oxide catalysts show drastically different activity when hydroxylated
- **Water splitting**: Hydroxylated surfaces are intermediates in photocatalytic water oxidation
- **Surface stability**: Hydroxylation can stabilize otherwise unstable surface terminations
- **Electronic properties**: OH groups modify work function and band alignment

The hydroxylation coverage is calculated as:
```
Coverage (OH/nm²) = N_OH / A_surface
```

Where:
- `N_OH`: Number of OH groups added
- `A_surface`: Surface area in nm²

### Oxygen Vacancies

Oxygen vacancies are created by removing oxygen atoms from the surface. These defects:

- Introduce localized electronic states in the band gap
- Create active sites for catalytic reactions
- Modify electrical conductivity (especially in oxides)
- Affect stability and phase transitions

The vacancy coverage is calculated similarly:
```
Coverage (vac/nm²) = N_vac / A_surface
```

### Coverage-Based Deduplication

For a 4-oxygen surface, there are C(4,2) = 6 ways to add 2 OH groups. For larger surfaces, this grows exponentially:
- 10 oxygen atoms: 1,024 total combinations
- 20 oxygen atoms: 1,048,576 total combinations

**The coverage-based deduplication algorithm**:
1. Generates all unique configurations
2. Groups by coverage (OH/nm² or vac/nm²)
3. Bins into N equal-width coverage ranges
4. Samples representative structures from each bin

**Result**: Reduces thousands of structures to ~N representative configurations spanning 0-100% coverage.

## Module Architecture

### Builder Function Pattern

The module follows PS-TEROS two-tier design:

```python
# High-level (user-facing)
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=params,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,  # Complete VASP builder configuration
    max_parallel_jobs=3,
    fix_type='bottom',              # Optional: fix bottom atoms
    fix_thickness=5.0,              # Optional: fix bottom 5Å
    structure_specific_builder_inputs={  # Optional: per-structure overrides
        0: {'parameters': {'incar': {'ALGO': 'Normal'}}},
        2: {'kpoints_spacing': 0.2},
    }
)
```

**Benefits**:
- Input validation with clear error messages
- Sensible defaults for optional parameters
- Automatic code/structure loading
- Consistent API across PS-TEROS modules

### Workflow Components

```
SurfaceHydroxylationWorkGraph
  │
  ├─> generate_structures (CalcFunction)
  │   └─ Uses SurfaceModifier to create variants
  │
  ├─> relax_slabs_with_semaphore (@task.graph)
  │   └─ Batch VASP relaxations with vasp.v2.vasp
  │
  └─> Returns raw namespace outputs:
      - manifest (Dict): Metadata for all variants
      - structures (namespace): {idx_variantname: StructureData}
      - energies (namespace): {idx_variantname: Float}
```

### Output Naming Convention

Outputs use descriptive keys in the format: `{index}_{variant_name}`

**Examples**:
- `0_oh_000_3_7572`: First structure, 3.76 OH/nm²
- `1_oh_001_7_5145`: Second structure, 7.51 OH/nm²
- `2_vac_000_5_2341`: Third structure, 5.23 vacancies/nm²

**Note**: Dots are replaced with underscores for AiiDA link label compatibility.

## Usage

### Basic Example

```python
from aiida import orm
from teros.core.surface_hydroxylation import (
    build_surface_hydroxylation_workgraph,
    organize_hydroxylation_results,
)

# 1. Define surface parameters
surface_params = {
    'mode': 'hydrogen',              # 'hydrogen', 'vacancies', or 'combine'
    'species': 'O',                  # Target oxygen atoms
    'z_window': 0.5,                 # Surface detection window (Å)
    'which_surface': 'top',          # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,                 # O-H bond distance (Å)
    'coverage_bins': 5,              # Number of coverage bins
    'deduplicate_by_coverage': True, # Enable deduplication
}

# 2. Configure VASP builder (complete WorkChain parameters)
builder_inputs = {
    'parameters': {
        'incar': {  # Note: nested under 'incar' key
            'PREC': 'Accurate',
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'ISIF': 2,          # Relax positions only
            'NSW': 200,
            'IBRION': 2,
            'EDIFFG': -0.02,    # Force convergence (eV/Å)
            'LWAVE': False,
            'LCHARG': False,
        }
    },
    'kpoints_spacing': 0.3,
    'potential_family': 'PBE',
    'potential_mapping': {},  # Auto-detect potentials
    'options': {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'queue_name': 'par40',
        'max_wallclock_seconds': 3600 * 10,
    },
    'clean_workdir': False,
}

# 3. Build and submit workflow
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,
    max_parallel_jobs=3,
)
result = wg.submit()
pk = result.pk

# 5. After completion, organize results
node = orm.load_node(pk)
results = organize_hydroxylation_results(node)

print(f"Total: {results['statistics']['total']}")
print(f"Succeeded: {results['statistics']['succeeded']}")
print(f"Failed: {results['statistics']['failed']}")

# 6. Analyze successful relaxations
for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV (coverage={r['coverage']:.2f})")
```

### Hydroxylation Mode

Study OH coverage from 0-100%:

```python
surface_params = {
    'mode': 'hydrogen',
    'species': 'O',
    'coverage_bins': 10,  # 10 representative coverages
}
```

**Use cases**:
- Water dissociation studies
- Proton transfer mechanisms
- Catalytic cycle intermediates

### Vacancy Mode

Study oxygen-deficient surfaces:

```python
surface_params = {
    'mode': 'vacancies',
    'species': 'O',
    'coverage_bins': 8,
}
```

**Use cases**:
- Defect formation energies
- Electronic structure of reduced oxides
- Oxygen evolution reaction (OER) mechanisms

### Combined Mode

Study both hydroxylation and vacancies:

```python
surface_params = {
    'mode': 'combine',
    'species': 'O',
    'coverage_bins': 10,
}
```

**Use cases**:
- Comprehensive surface phase diagrams
- Competitive adsorption studies
- Realistic atmospheric conditions

## Parameters Reference

### surface_params (Dict)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | str | Required | 'hydrogen', 'vacancies', or 'combine' |
| `species` | str | 'O' | Target atomic species for modification |
| `z_window` | float | 0.5 | Z-coordinate window for surface detection (Å) |
| `which_surface` | str | 'top' | 'top', 'bottom', or 'both' |
| `oh_dist` | float | 0.98 | O-H bond distance for hydroxylation (Å) |
| `include_empty` | bool | False | Include structure with no modifications |
| `supercell` | list or None | None | Supercell expansion [nx, ny, nz] |
| `deduplicate_by_coverage` | bool | True | Enable coverage-based deduplication |
| `coverage_bins` | int | 5 | Number of coverage bins for sampling |

### builder_inputs (Dict)

Complete VASP builder configuration for `vasp.v2.vasp` WorkChain.

**Required keys:**
```python
builder_inputs = {
    'parameters': {
        'incar': {  # IMPORTANT: INCAR tags nested under 'incar' key
            'PREC': 'Accurate',
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'ALGO': 'Normal',
            'LREAL': False,
            'LWAVE': False,
            'LCHARG': False,
            'ISIF': 2,           # Ionic relaxation mode
            'NSW': 200,          # Max ionic steps
            'IBRION': 2,         # Conjugate gradient
            'EDIFFG': -0.02,     # Force convergence
        }
    },
    'potential_family': 'PBE',  # Or 'PBE.54', 'LDA', etc.
    'options': {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'queue_name': 'par40',
        'max_wallclock_seconds': 3600 * 10,
    },
}
```

**Optional keys:**
- `kpoints_spacing`: Float (Å⁻¹, default: 0.5)
- `potential_mapping`: Dict (element→potential, default: {})
- `settings`: Dict (parser settings, default: add trajectory/structure/kpoints)
- `dynamics`: Dict (selective dynamics - auto-added if fix_type is set)
- `clean_workdir`: Bool (default: False)
- `max_iterations`: Int (max restart attempts, default: 5)
- `verbose`: Bool (detailed logging, default: False)

**Note**: The workflow automatically sets `code` and `structure` for each relaxation.

### max_parallel_jobs (Int)

Batch control: processes first N structures only.

**Strategy for large calculations**:
```python
# Step 1: Test with 2 structures
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=2)

# Step 2: After verification, process more
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=10)
```

## Post-Processing Results

### organize_hydroxylation_results()

**Why needed**: WorkGraph cannot pass namespace dicts containing AiiDA nodes to CalcFunctions due to JSON serialization limitations. Solution: return raw namespaces + provide Python helper.

```python
from teros.core.surface_hydroxylation import organize_hydroxylation_results

node = orm.load_node(workflow_pk)
results = organize_hydroxylation_results(node)

# Results structure:
{
    'successful_relaxations': [
        {
            'name': 'oh_000_3.7572',
            'structure_pk': 12345,
            'energy': -96.375793,
            'coverage': 3.7572,
            'metadata': {...}  # Full variant info
        },
        ...
    ],
    'failed_relaxations': [
        {
            'name': 'oh_005_15.0289',
            'coverage': 15.0289,
            'error_message': 'Relaxation failed - no structure/energy output'
        },
        ...
    ],
    'statistics': {
        'total': 10,
        'succeeded': 8,
        'failed': 2
    }
}
```

### Finding Optimal Coverage

```python
results = organize_hydroxylation_results(node)

# Sort by energy
sorted_results = sorted(
    results['successful_relaxations'],
    key=lambda r: r['energy']
)

# Lowest energy configuration
best = sorted_results[0]
print(f"Optimal coverage: {best['coverage']:.2f} OH/nm²")
print(f"Energy: {best['energy']:.6f} eV")

# Load structure for further analysis
best_structure = orm.load_node(best['structure_pk'])
```

### Coverage vs Energy Plots

```python
import matplotlib.pyplot as plt

coverages = [r['coverage'] for r in results['successful_relaxations']]
energies = [r['energy'] for r in results['successful_relaxations']]

plt.figure(figsize=(8, 6))
plt.scatter(coverages, energies, s=100, alpha=0.7)
plt.xlabel('OH Coverage (OH/nm²)')
plt.ylabel('Total Energy (eV)')
plt.title('Surface Energy vs Hydroxylation Coverage')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('coverage_vs_energy.png', dpi=300)
```

## Advanced Usage

### Custom Surface Detection

For surfaces with complex terminations:

```python
surface_params = {
    'mode': 'hydrogen',
    'species': 'O',
    'z_window': 1.0,        # Larger window for rough surfaces
    'which_surface': 'both', # Modify both surfaces
}
```

### Supercell Expansion

For low-coverage studies:

```python
surface_params = {
    'mode': 'hydrogen',
    'species': 'O',
    'supercell': [2, 2, 1],  # 2×2 expansion in xy plane
    'coverage_bins': 10,
}
```

### Selective Dynamics (Fixing Substrate Atoms)

**NEW FEATURE**: Fix substrate atoms during relaxation using VASP selective dynamics.

When studying surface hydroxylation, it's often desirable to fix the substrate atoms and only relax the adsorbates (OH groups). This is accomplished using three parameters:

```python
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,
    max_parallel_jobs=3,
    # Selective dynamics parameters
    fix_type='bottom',       # Where to fix: 'bottom', 'top', or 'center'
    fix_thickness=5.0,       # Fix region thickness in Angstroms
    fix_elements=['Ag', 'P'],  # Optional: only fix specific elements
)
```

**Parameters:**
- `fix_type`: `str` - Where to fix atoms (default: `None` = no fixing)
  - `'bottom'`: Fix bottom region of slab (most common)
  - `'top'`: Fix top region of slab
  - `'center'`: Fix central region of slab
  - `None`: No atoms fixed (all relax freely)

- `fix_thickness`: `float` - Thickness of fix region in Angstroms (default: 0.0)
  - Example: `5.0` = fix bottom 5Å of slab

- `fix_elements`: `list[str]` or `None` - Element symbols to fix (default: `None`)
  - `None`: Fix all elements in the specified region
  - `['Ag', 'P']`: Only fix Ag and P atoms in the region, O and H relax freely
  - Useful when you want adsorbates to relax but substrate to stay fixed

**How it works:**
1. Workflow calculates Z-coordinates of all atoms
2. Identifies atoms in the fix region based on `fix_type` and `fix_thickness`
3. Optionally filters by element symbols if `fix_elements` is provided
4. Creates VASP `positions_dof` array: `[True, True, True]` = relax, `[False, False, False]` = fix
5. Adds `dynamics` dict to VASP builder inputs automatically

**Common use cases:**

1. **Fix entire substrate, relax adsorbates:**
```python
fix_type='bottom'
fix_thickness=5.0
fix_elements=None  # Fix all atoms in bottom 5Å
```

2. **Fix only metal atoms, allow surface oxygen to relax:**
```python
fix_type='bottom'
fix_thickness=5.0
fix_elements=['Ag']  # Fix only Ag in bottom 5Å, O/P relax
```

3. **No fixing (default behavior):**
```python
fix_type=None  # All atoms relax freely
```

**Note**: The `dynamics` parameter in `builder_inputs` will be automatically overwritten if `fix_type` is set. To have full manual control over selective dynamics, set `fix_type=None` and provide `dynamics` directly in `builder_inputs`.

### Structure-Specific Builder Inputs

**NEW FEATURE**: Override VASP parameters for specific structures using deep merge.

Sometimes specific structures need different VASP settings (e.g., failed calculations, difficult convergence). Instead of rerunning the entire workflow, you can provide structure-specific overrides:

```python
# Default builder used for all structures
builder_inputs = {
    'parameters': {
        'incar': {
            'PREC': 'Normal',
            'ENCUT': 400,
            'ALGO': 'Fast',
            'EDIFFG': -0.02,
            # ... other INCAR tags
        }
    },
    'kpoints_spacing': 0.5,
    'potential_family': 'PBE',
    'options': {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
        'max_wallclock_seconds': 3600,
    },
}

# Override parameters for specific structures
structure_specific_builder_inputs = {
    0: {  # Structure 0: increase precision
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 600,
            }
        },
        'kpoints_spacing': 0.2,  # Denser k-points
    },
    2: {  # Structure 2: different algorithm + more resources
        'parameters': {
            'incar': {
                'ALGO': 'Normal',  # Slower but more stable
                'EDIFFG': -0.05,   # Looser convergence
            }
        },
        'options': {
            'resources': {'num_machines': 2},  # More resources
            'max_wallclock_seconds': 7200,     # More time
        },
    },
    5: {  # Structure 5: completely custom settings
        'parameters': {
            'incar': {
                'IBRION': 1,   # RMM-DIIS instead of CG
                'POTIM': 0.2,  # Smaller step size
            }
        },
    }
    # Structures 1, 3, 4, 6, ... use default builder_inputs
}

wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params=surface_params,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,  # Default for all
    structure_specific_builder_inputs=structure_specific_builder_inputs,  # Overrides
    max_parallel_jobs=10,
)
```

**How deep merging works:**

The workflow uses **deep merge** strategy, which means:
- Only specified parameters are overridden
- All other parameters are inherited from default `builder_inputs`
- Nested dicts are recursively merged

**Example:**

```python
# Default
builder_inputs = {
    'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400, 'EDIFF': 1e-5, 'ALGO': 'Fast'}},
    'kpoints_spacing': 0.5,
    'potential_family': 'PBE',
    'options': {'resources': {'num_machines': 1}},
}

# Structure 0 override
structure_specific_builder_inputs = {
    0: {
        'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 600}},
        'kpoints_spacing': 0.2,
    }
}

# Resulting builder for structure 0 (after deep merge):
{
    'parameters': {
        'incar': {
            'PREC': 'Accurate',  # Overridden
            'ENCUT': 600,        # Overridden
            'EDIFF': 1e-5,       # Inherited from default
            'ALGO': 'Fast',      # Inherited from default
        }
    },
    'kpoints_spacing': 0.2,          # Overridden
    'potential_family': 'PBE',       # Inherited from default
    'options': {'resources': {'num_machines': 1}},  # Inherited from default
}
```

**Common use cases:**

1. **Rerun only failed calculations with different parameters:**
```python
# After first run, identify failed structures (indices 2, 5, 7)
structure_specific_builder_inputs = {
    2: {'parameters': {'incar': {'ALGO': 'Normal', 'EDIFFG': -0.05}}},
    5: {'parameters': {'incar': {'IBRION': 1, 'NSW': 300}}},
    7: {'options': {'max_wallclock_seconds': 7200}},  # More time
}
# Set max_parallel_jobs=10 to include all structures, but only 2,5,7 will use custom settings
```

2. **High-precision calculation for specific coverage:**
```python
# Use high precision only for low-coverage structures (indices 0, 1)
structure_specific_builder_inputs = {
    0: {'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 600}}, 'kpoints_spacing': 0.2},
    1: {'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 600}}, 'kpoints_spacing': 0.2},
}
```

3. **Different algorithms for difficult structures:**
```python
# Structure 3 has convergence issues - try different algorithm
structure_specific_builder_inputs = {
    3: {
        'parameters': {
            'incar': {
                'ALGO': 'All',      # Try all algorithms sequentially
                'IBRION': 1,        # RMM-DIIS
                'POTIM': 0.2,       # Smaller ionic step
                'EDIFFG': -0.05,    # Looser convergence
            }
        },
        'options': {'max_wallclock_seconds': 10800},  # 3 hours instead of 1
    }
}
```

**Important notes:**
- Keys are **structure indices** (0, 1, 2, ...) corresponding to generation order
- Indices match the numeric prefix in output keys (e.g., `0_oh_000_3_7572` → index 0)
- Integer keys are automatically converted to strings for AiiDA serialization
- All structures not listed will use the default `builder_inputs`
- This is the **recommended pattern** for all PS-TEROS modules going forward

### Direct Namespace Access

For custom post-processing:

```python
node = orm.load_node(workflow_pk)

# Raw outputs
manifest = node.outputs.manifest.get_dict()
structures = node.outputs.structures  # Namespace dict
energies = node.outputs.energies      # Namespace dict

# Access specific structure
structure_key = '0_oh_000_3_7572'
structure = structures[structure_key]
energy = energies[structure_key].value

# Get variant metadata
variant = manifest['variants'][0]
print(f"Coverage: {variant['OH_coverage']}")
print(f"Decorated atoms: {variant['decorated_indices']}")
```

## Troubleshooting

### All Structures Failed

**Problem**: All relaxations fail immediately.

**Solutions**:
1. Verify VASP code: `verdi code show <code_label>`
2. Check pseudopotentials: `verdi data core.upf listfamilies`
3. Validate VASP parameters in `builder_inputs['parameters']['incar']`
4. Test single structure relaxation manually

### Some Structures Failed

**Problem**: High-coverage or low-coverage configurations fail to converge.

**Explanation**: Extreme configurations may be:
- Physically unstable (too many/few OH groups)
- Require different relaxation settings
- Have geometry optimization difficulties

**Solutions**:
```python
# Option 1: Use structure-specific builder inputs for failed structures
structure_specific_builder_inputs = {
    2: {  # Assuming structure 2 failed
        'parameters': {
            'incar': {
                'EDIFFG': -0.05,  # Looser (was -0.02)
                'NSW': 300,        # More steps
                'IBRION': 1,       # RMM-DIIS instead of CG
            }
        }
    }
}

# Option 2: Use looser convergence for all structures
builder_inputs = {
    'parameters': {
        'incar': {
            'EDIFFG': -0.05,  # Looser (was -0.02)
            'NSW': 300,        # More steps
            'IBRION': 1,       # RMM-DIIS instead of CG
            # ... other INCAR tags
        }
    },
    # ... other builder parameters
}
```

### No Structures Generated

**Problem**: `generate_structures` returns empty or very few structures.

**Solutions**:
1. Check `z_window` - might be too restrictive
2. Verify target species exists in surface region
3. Try `include_empty=True` to get baseline
4. Increase `coverage_bins` for more sampling

### Cache Issues

**Problem**: Changes to code not reflected in new runs.

**Solutions**:
```bash
# Clear Python cache
find . -type d -name __pycache__ -exec rm -rf {} +
find . -name "*.pyc" -delete

# Restart daemon
verdi daemon restart

# Use different parameters to bypass AiiDA cache
surface_params = {'coverage_bins': 6}  # Changed from 5
```

## Performance Considerations

### Computational Cost

For a typical perovskite oxide surface (~100 atoms):

| Coverage Bins | Structures Generated | VASP Time (each) | Total Time |
|---------------|---------------------|------------------|------------|
| 5 | ~5-8 | 2-4 hours | 10-32 hours |
| 10 | ~10-15 | 2-4 hours | 20-60 hours |
| 20 | ~20-30 | 2-4 hours | 40-120 hours |

**Recommendation**: Start with `coverage_bins=5`, analyze results, then increase if needed.

### Batch Strategy

```python
# Day 1: Test run
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=2)
# Verify results look reasonable

# Day 2: Partial production
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=5)

# Day 3: Full production
wg = build_surface_hydroxylation_workgraph(..., max_parallel_jobs=15)
```

### Memory Requirements

- **Structure generation**: Minimal (<1 GB RAM)
- **VASP relaxations**: Depends on system size and VASP settings
  - 100-atom slab: ~4-8 GB per job
  - 200-atom slab: ~8-16 GB per job

**Cluster configuration**:
```python
builder_inputs = {
    'parameters': {'incar': {...}},  # INCAR settings
    'potential_family': 'PBE',
    'kpoints_spacing': 0.3,
    'options': {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'max_wallclock_seconds': 3600 * 10,  # 10 hours
    },
}
```

## Integration with PS-TEROS Workflows

### Input from Surface Thermodynamics

```python
# Get relaxed slab from surface_thermodynamics workflow
slab_workflow = orm.load_node(12345)
relaxed_slab_pk = slab_workflow.outputs.relaxed_slabs['LaO_term'].pk

# Use as input for hydroxylation
wg = build_surface_hydroxylation_workgraph(
    structure_pk=relaxed_slab_pk,
    surface_params=surface_params,
    ...
)
```

### Output to Electronic Structure Analysis

```python
# Get best hydroxylated structure
results = organize_hydroxylation_results(node)
best = sorted(results['successful_relaxations'], key=lambda r: r['energy'])[0]

# Use in electronic structure workflow
from teros.core.electronic_properties import build_dos_workflow

dos_wg = build_dos_workflow(
    structure_pk=best['structure_pk'],
    ...
)
```

## API Reference

### build_surface_hydroxylation_workgraph()

```python
def build_surface_hydroxylation_workgraph(
    structure: orm.StructureData = None,
    structure_pk: int = None,
    surface_params: dict = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    builder_inputs: dict = None,
    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list[str] = None,
    structure_specific_builder_inputs: dict = None,
    name: str = 'SurfaceHydroxylation',
) -> WorkGraph
```

Build a WorkGraph for surface hydroxylation/vacancy calculations.

**Parameters:**
- `structure` (StructureData): Input relaxed slab structure (either this or `structure_pk`)
- `structure_pk` (int): PK of relaxed slab structure (either this or `structure`)
- `surface_params` (dict): Surface modification parameters (see Parameters Reference)
- `code_label` (str): Label of VASP code (default: 'VASP-VTST-6.4.3@bohr')
- `builder_inputs` (dict): Complete VASP builder configuration (see Parameters Reference)
- `max_parallel_jobs` (int): Number of structures to process (default: 2)
- `fix_type` (str): Where to fix atoms - 'bottom'/'top'/'center'/None (default: None)
- `fix_thickness` (float): Thickness of fix region in Angstroms (default: 0.0)
- `fix_elements` (list[str]): Optional element symbols to fix (default: None)
- `structure_specific_builder_inputs` (dict): Per-structure VASP overrides (default: None)
- `name` (str): Workflow name (default: 'SurfaceHydroxylation')

**Returns:**
- `WorkGraph`: Ready-to-submit workflow instance

**Example:**
```python
wg = build_surface_hydroxylation_workgraph(
    structure_pk=1234,
    surface_params={'mode': 'hydrogen', 'coverage_bins': 5},
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,
    max_parallel_jobs=3,
    fix_type='bottom',
    fix_thickness=5.0,
)
result = wg.submit()
```

### SurfaceHydroxylationWorkGraph()

```python
@task.graph(outputs=['manifest', 'structures', 'energies'])
def SurfaceHydroxylationWorkGraph(
    structure: orm.StructureData,
    surface_params: dict,
    code: orm.InstalledCode,
    builder_inputs: dict,
    max_parallel_jobs: int = 2,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list[str] = None,
    structure_specific_builder_inputs: dict = None,
) -> dict
```

Main workflow task graph for surface hydroxylation (lower-level API).

**Parameters:** Same as `build_surface_hydroxylation_workgraph()` but takes `code` object instead of `code_label`.

**Returns:**
- `dict` with outputs:
  - `manifest` (Dict): Metadata for all generated variants
  - `structures` (namespace): {idx_variantname: StructureData} for successful relaxations
  - `energies` (namespace): {idx_variantname: Float} for successful relaxations

**Note:** Most users should use `build_surface_hydroxylation_workgraph()` instead.

### organize_hydroxylation_results()

```python
def organize_hydroxylation_results(workflow_node) -> dict
```

Organize hydroxylation workflow results into structured format.

**Parameters:**
- `workflow_node`: Completed WorkGraph node (orm.WorkChainNode)

**Returns:**
- `dict` with:
  - `successful_relaxations` (list): List of dicts with successful results
    - Each contains: `name`, `structure_pk`, `energy`, `coverage`, `metadata`
  - `failed_relaxations` (list): List of dicts with failed results
    - Each contains: `name`, `coverage`, `error_message`
  - `statistics` (dict): Summary with `total`, `succeeded`, `failed` counts

**Example:**
```python
node = orm.load_node(workflow_pk)
results = organize_hydroxylation_results(node)

print(f"Succeeded: {results['statistics']['succeeded']}")
print(f"Failed: {results['statistics']['failed']}")

for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV")
```

### SurfaceModifier

```python
class SurfaceModifier:
    def __init__(
        self,
        atoms: Atoms,
        species: str = "O",
        z_window: float = 0.5,
        which_surface: str = "top",
        oh_dist: float = 0.98,
        include_empty: bool = False,
        outdir: Path | str = "outputs",
        fmt: str = "vasp",
        supercell: tuple[int, int, int] | None = None,
        deduplicate_by_coverage: bool = False,
        coverage_bins: int | None = None,
    )
```

Surface structure modification engine (used internally by `generate_structures` task).

**Methods:**
- `run_vacancies()`: Generate vacancy configurations
- `run_hydrogen()`: Generate hydroxylation configurations
- `run_combine()`: Generate combined vacancy + hydroxylation configurations
- `run_complete()`: Run all three modes

**Note:** This class is used internally. Most users interact with it through `surface_params` in the workflow functions.

## Citation

If you use this module in your research, please cite:

```bibtex
@software{psteros_surface_hydroxylation,
  title = {PS-TEROS Surface Hydroxylation Module},
  author = {{PS-TEROS Development Team}},
  year = {2025},
  url = {https://github.com/your-org/PS-TEROS},
}
```

## References

1. **Hydroxylation on oxide surfaces**: Carrasco et al., J. Chem. Phys. 2012
2. **Oxygen vacancy formation**: Ganduglia-Pirovano et al., Surf. Sci. Rep. 2007
3. **Coverage-based deduplication**: Original algorithm developed for PS-TEROS
4. **AiiDA-WorkGraph**: Uhrin et al., Computational Materials Science 2021

## Support

For issues or questions:
1. Check this documentation and README troubleshooting
2. Review example scripts in `examples/surface_hydroxylation/`
3. Check AiiDA provenance: `verdi process report <PK>`
4. Contact PS-TEROS development team
