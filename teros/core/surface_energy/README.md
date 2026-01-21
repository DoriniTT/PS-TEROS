# Metal and Intermetallic Surface Energy Module

This module computes surface energies for **elemental metals** and **stoichiometric intermetallics** using the simple thermodynamic formula:

$$\gamma = \frac{E_{\text{slab}} - N \cdot E_{\text{bulk}}^{\text{atom}}}{2A}$$

Where:
- $E_{\text{slab}}$: Total DFT energy of the slab
- $N$: Number of atoms in the slab
- $E_{\text{bulk}}^{\text{atom}}$: Bulk energy per atom
- $A$: Surface area (Å²)
- Factor of 2 accounts for two surfaces (top and bottom)

## Supported Materials

### Elemental Metals
Single-element metals like Au, Ag, Cu, Pt, Pd, Ni, Fe, etc.

### Stoichiometric Intermetallics
Binary or ternary intermetallic compounds where:
- **Stoichiometric surface**: The slab composition matches the bulk stoichiometry
- **Symmetric surface**: Top and bottom surfaces are equivalent

Examples: PdIn, AuCu, NiAl, Cu3Au, TiAl, Ni3Al, etc.

> **Note**: For non-stoichiometric or asymmetric intermetallic surfaces, chemical
> potential-dependent formulations are required (to be implemented separately).

## Key Differences from Oxide Thermodynamics

| Feature | Metal/Intermetallic Surface Energy | Oxide Thermodynamics |
|---------|-----------------------------------|----------------------|
| Chemical potential | None (stoichiometric) | Δμ_O, Δμ_Metal dependencies |
| Formula | Simple single-value γ | γ(Δμ) function |
| Reference calculations | None needed | O₂, metal references required |
| Complexity | Low | High |
| Applicability | Stoichiometric surfaces | Any composition |

## Module Structure

```
surface_energy/
├── __init__.py               # Module exports
├── surface_energy.py         # Surface energy calcfunctions
├── workgraph.py              # WorkGraph builder
├── wulff.py                  # Wulff shape construction
├── stoichiometric_finder.py  # Stoichiometric+symmetric finder (EXPERIMENTAL)
└── README.md                 # This file
```

## Key Functions

### `identify_compound_type`
Helper function that identifies whether a structure is an elemental metal or intermetallic.

**Returns:** Dict with:
- `compound_type`: 'elemental' or 'intermetallic'
- `elements`: List of element symbols
- `composition`: Dict mapping element to count
- `formula`: Chemical formula (e.g., 'Au', 'PdIn', 'Au3Cu')
- `stoichiometry`: Dict mapping element to ratio

### `calculate_metal_surface_energy`
AiiDA calcfunction that computes surface energy for a single slab.

**Inputs:**
- `bulk_structure`: Relaxed bulk StructureData
- `bulk_energy`: Total energy of bulk (Float)
- `slab_structure`: Relaxed slab StructureData
- `slab_energy`: Total energy of slab (Float)

**Returns:** Dict with:
- `gamma_eV_A2`: Surface energy in eV/Å²
- `gamma_J_m2`: Surface energy in J/m²
- `area_A2`: Surface area in Å²
- `N_slab`, `N_bulk`: Atom counts
- `compound_type`: 'elemental' or 'intermetallic'
- `formula`: Chemical formula
- `elements`: List of elements
- `composition`: Bulk composition dict
- `slab_composition`: Slab composition dict
- `is_stoichiometric`: Boolean flag

### `build_metal_surface_energy_workgraph`
Main entry point for users. Builds a complete workflow.

**Key Parameters:**
- `bulk_structure_path`: Path to CIF file
- `miller_indices`: List of orientations, e.g., `[[1,1,1], [1,0,0], [1,1,0]]`
- `code_label`: VASP code in AiiDA
- `potential_mapping`: e.g., `{'Au': 'Au'}` or `{'Pd': 'Pd', 'In': 'In'}`

### Wulff Shape Functions (NEW)

#### `build_wulff_shape`
AiiDA calcfunction that constructs the Wulff shape from surface energies. Automatically integrated into the workflow.

**Inputs:**
- `bulk_structure`: Bulk crystal structure (for lattice and symmetry)
- `surface_energies`: Output from gather_surface_energies

**Returns:** Dict with shape analysis (see Wulff Shape section below)

#### `get_symmetrically_equivalent_miller_indices`
Get all symmetrically equivalent Miller indices for a given orientation.

```python
get_symmetrically_equivalent_miller_indices(structure, (1, 1, 1))
# Returns: [(1,1,1), (-1,-1,-1), (1,1,-1), (1,-1,1), ...]
```

#### `expand_surface_energies_with_symmetry`
Expand calculated surface energies to all symmetry-equivalent orientations.

```python
expand_surface_energies_with_symmetry(structure, surface_energies_dict)
# Returns: {(1,1,1): 0.79, (-1,-1,-1): 0.79, ...}
```

#### `visualize_wulff_shape`
Create a 3D visualization of the Wulff shape (requires matplotlib).

```python
visualize_wulff_shape(bulk_structure, surface_energies, save_path='wulff.png')
```

#### `get_wulff_shape_summary`
Generate a human-readable summary of Wulff shape results.

```python
print(get_wulff_shape_summary(wg.outputs.wulff_shape))
```

## Workflow Architecture

```
WorkGraph<MetalSurfaceEnergy>
├── bulk_relax (VaspWorkChain)          # Runs ONCE
├── bulk_energy (extract_total_energy)
├── bulk_energy_per_atom
├── surface_hkl_111                      # Per-orientation sub-workflows
│   ├── generate_slab_structures
│   ├── relax_slabs_scatter
│   └── compute_metal_surface_energies_scatter
├── surface_hkl_100
├── surface_hkl_110
├── gather_surface_energies             # Consolidates all results
└── build_wulff_shape                   # Constructs equilibrium crystal shape (NEW)
```

## Usage Example - Elemental Metal (Au)

```python
from aiida import load_profile
from teros.core.surface_energy import build_metal_surface_energy_workgraph

load_profile('myprofile')

wg = build_metal_surface_energy_workgraph(
    bulk_structure_path='/path/to/au.cif',
    code_label='VASP-6.5.1@cluster',
    potential_family='PBE',
    potential_mapping={'Au': 'Au'},

    # Miller indices (all in ONE workflow)
    miller_indices=[[1,1,1], [1,0,0], [1,1,0]],

    # Slab parameters
    min_slab_thickness=20.0,
    min_vacuum_thickness=20.0,

    # VASP parameters
    bulk_parameters={
        'prec': 'Accurate',
        'encut': 500,
        'ismear': 1,
        'sigma': 0.2,
        'ibrion': 2,
        'isif': 3,  # Full relaxation for bulk
        'nsw': 100,
        'ediffg': -0.01,
    },

    name='Au_surface_energy',
)

wg.submit(wait=False)
print(f"Submitted: PK={wg.pk}")
```

## Usage Example - Intermetallic (PdIn)

```python
from aiida import load_profile
from teros.core.surface_energy import build_metal_surface_energy_workgraph

load_profile('myprofile')

wg = build_metal_surface_energy_workgraph(
    bulk_structure_path='/path/to/pdin.cif',
    code_label='VASP-6.5.1@cluster',
    potential_family='PBE',
    potential_mapping={'Pd': 'Pd', 'In': 'In'},  # Map both elements

    # Miller indices for intermetallic
    miller_indices=[[1,1,0], [1,0,0], [1,1,1]],

    # Slab parameters
    min_slab_thickness=20.0,
    min_vacuum_thickness=20.0,
    symmetrize=True,  # Important for stoichiometric surfaces!

    # VASP parameters
    bulk_parameters={
        'prec': 'Accurate',
        'encut': 500,
        'ismear': 1,
        'sigma': 0.2,
        'ibrion': 2,
        'isif': 3,
        'nsw': 100,
        'ediffg': -0.01,
    },

    name='PdIn_surface_energy',
)

wg.submit(wait=False)
print(f"Submitted: PK={wg.pk}")
```

## Output Structure

After completion, `verdi process show <PK>` displays:

```
Outputs                  PK    Type
-----------------------  ----  -------------
bulk_energy              XXXX  Float
bulk_energy_per_atom     XXXX  Float
bulk_structure           XXXX  StructureData
surface_energies         XXXX  Dict          <- All results consolidated
wulff_shape              XXXX  Dict          <- Wulff shape analysis (NEW)
relaxed_slabs_hkl_111
    term_0               XXXX  StructureData
relaxed_slabs_hkl_100
    term_0               XXXX  StructureData
...
```

The `surface_energies` Dict contains:

### For Elemental Metals:
```json
{
    "hkl_111": {
        "term_0": {
            "gamma_eV_A2": 0.048,
            "gamma_J_m2": 0.77,
            "area_A2": 21.5,
            "N_slab": 8,
            "compound_type": "elemental",
            "formula": "Au",
            "elements": ["Au"],
            "is_stoichiometric": true
        }
    }
}
```

### For Intermetallics:
```json
{
    "hkl_110": {
        "term_0": {
            "gamma_eV_A2": 0.082,
            "gamma_J_m2": 1.31,
            "area_A2": 15.2,
            "N_slab": 16,
            "compound_type": "intermetallic",
            "formula": "InPd",
            "elements": ["In", "Pd"],
            "composition": {"In": 4, "Pd": 4},
            "slab_composition": {"In": 8, "Pd": 8},
            "is_stoichiometric": true
        }
    }
}
```

## Expected Surface Energies (PBE)

### Elemental Metals
| Metal | (111) | (100) | (110) |
|-------|-------|-------|-------|
| Au    | ~0.8  | ~0.9  | ~1.3  |
| Ag    | ~0.6  | ~0.7  | ~0.9  |
| Cu    | ~1.3  | ~1.5  | ~1.6  |
| Pt    | ~1.5  | ~1.9  | ~2.0  |

*Values in J/m². Actual results depend on computational parameters.*

### Intermetallics
Surface energies for intermetallics vary widely depending on:
- Crystal structure (B2, L10, L12, etc.)
- Element combination
- Surface orientation and termination

Typical values range from 0.5 to 2.5 J/m².

## VASP Parameter Recommendations

For metals and intermetallics, use:
- `ISMEAR = 1` (Methfessel-Paxton) or `ISMEAR = 2`
- `SIGMA = 0.1-0.2` eV
- `ENCUT >= 1.3 × ENMAX` from POTCAR
- Dense k-points: spacing <= 0.03 Å⁻¹

For slabs:
- `ISIF = 2` (relax ions only, fix cell)
- Sufficient vacuum (>= 15 Å)
- Consider dipole corrections if asymmetric
- Use `symmetrize=True` in slab generation for stoichiometric surfaces

## Important Considerations for Intermetallics

1. **Stoichiometry**: Ensure the slab generation produces stoichiometric slabs by using `symmetrize=True`

2. **Termination**: For non-centrosymmetric intermetallics, different terminations may have different compositions. The `is_stoichiometric` flag in the output indicates whether each termination preserves bulk stoichiometry.

3. **Validation**: Always check the `slab_composition` in the output to verify the surface composition matches expectations.

4. **Future Extension**: For non-stoichiometric surfaces, chemical potential-dependent surface energies (γ(Δμ)) will be implemented in a separate module.

---

## Wulff Shape Construction (NEW)

The Wulff construction determines the equilibrium crystal shape that minimizes total surface energy for a fixed volume. This module automatically constructs the Wulff shape from calculated surface energies.

### Theory

The Wulff shape is the convex hull where each facet's distance from the center is proportional to its surface energy. Low-energy facets dominate the equilibrium shape.

### Symmetry Expansion

For cubic crystals (FCC, BCC), symmetry dramatically reduces required DFT calculations:

| Miller Family | Calculated | Symmetry Equivalents | Total Facets |
|---------------|------------|---------------------|--------------|
| {111} | 1 | 8 | 8 |
| {110} | 1 | 12 | 12 |
| {100} | 1 | 6 | 6 |
| **Total** | **3** | - | **26** |

Calculating just 3 unique orientations provides 26 facets for accurate Wulff shape construction.

### Wulff Shape Output

The `wulff_shape` Dict contains:

```json
{
    "wulff_shape_valid": true,
    "shape_factor": 5.34,
    "anisotropy": 0.072,
    "weighted_surface_energy": 0.82,
    "dominant_facet": "(-1, -1, -1)",
    "dominant_facet_fraction": 0.806,
    "facet_fractions": {
        "(-1, -1, -1)": 0.101,
        "(1, 1, 1)": 0.101,
        ...
    },
    "miller_energy_dict": {
        "(1, 1, 1)": 0.79,
        "(1, 0, 0)": 0.94,
        "(1, 1, 0)": 1.33
    },
    "expanded_miller_energy_dict": {
        "(1, 1, 1)": 0.79,
        "(-1, -1, -1)": 0.79,
        ...
    },
    "n_calculated_orientations": 3,
    "n_expanded_orientations": 26,
    "total_surface_area": 14.2,
    "volume": 1.0
}
```

| Field | Description |
|-------|-------------|
| `wulff_shape_valid` | Whether construction succeeded |
| `shape_factor` | 1.0 = sphere, >1 = anisotropic |
| `anisotropy` | Degree of shape anisotropy |
| `weighted_surface_energy` | Area-weighted average (J/m²) |
| `dominant_facet` | Miller index with largest area |
| `dominant_facet_fraction` | Fraction of total surface |
| `facet_fractions` | All facets and their area fractions |
| `miller_energy_dict` | Original calculated energies |
| `expanded_miller_energy_dict` | After symmetry expansion |

### Using Wulff Shape Results

```python
from aiida import orm
from teros.core.surface_energy import (
    get_wulff_shape_summary,
    visualize_wulff_shape
)

# Load completed workflow
wg = orm.load_node(12345)

# Print summary
print(get_wulff_shape_summary(wg.outputs.wulff_shape))
# Output:
# ============================================================
# WULFF SHAPE ANALYSIS
# ============================================================
#
# Orientations:
#   Calculated: 3
#   After symmetry expansion: 26
#
# Shape Properties:
#   Shape factor: 5.3415 (1.0 = sphere)
#   Anisotropy: 0.0724
#   Weighted surface energy: 0.8191 J/m²
#
# Dominant Facet:
#   (-1, -1, -1) (80.6% of surface)
# ...

# Visualize (requires matplotlib)
ax = visualize_wulff_shape(
    wg.outputs.bulk_structure,
    wg.outputs.surface_energies,
    color_set='PuBu',
    direction=(1, 1, 1),
    save_path='au_wulff.png'
)
```

### Limitations

1. **Single-digit Miller indices**: Key parsing assumes indices 0-9. High-index surfaces (≥10) not supported.

2. **Negative surface energies**: Surfaces with γ < 0 are skipped with a warning (indicate unstable surfaces or calculation errors).

3. **Minimum orientations**: For a valid 3D Wulff shape, at least 3 non-coplanar Miller families should be calculated.

### Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| `wulff_shape_valid: False` | No valid surface energies | Check DFT calculations completed |
| Negative energy warning | Unstable surface or error | Check VASP output, increase slab thickness |
| Flat Wulff shape | Low anisotropy | Normal for some materials |
| Missing facets | Not expanded from calculated | Add more Miller indices |

---

## Stoichiometric+Symmetric Surface Finder (EXPERIMENTAL)

This experimental feature helps find surfaces that are BOTH stoichiometric AND symmetric, which is challenging for intermetallics and oxides but valuable because it enables using the simple surface energy formula without chemical potential dependencies.

### Why This Matters

| Surface Type | Simple Formula | Thermodynamics |
|--------------|----------------|----------------|
| Stoichiometric + Symmetric | γ = (E_slab - N·E_bulk/atom) / 2A | Not needed |
| Non-stoichiometric OR asymmetric | Invalid | γ(Δμ) required |

For simple metals (Au, Ag, Cu), nearly all surfaces meet both criteria. For complex materials, finding valid surfaces requires systematic searching.

### Search Strategies

The finder uses multiple strategies in order:

1. **filter_first**: Generate all slabs with `symmetrize=False`, filter post-hoc (fastest)
2. **symmetrize_check**: Generate with `symmetrize=True`, verify stoichiometry preserved
3. **thickness_scan**: Scan thickness 10-30Å to find valid terminations
4. **bond_preservation**: Preserve chemical units (PO4, SO4) using bonds parameter

### Basic Usage

```python
from teros.core.surface_energy import find_stoichiometric_symmetric_slabs

# Find valid slabs for a Miller index
results = find_stoichiometric_symmetric_slabs(
    structure,
    miller_index=(1, 1, 0),
    min_slab_thickness=15.0,
)

if results:
    best_slab = results[0].slab
    print(f"Found valid slab: {best_slab.composition}")
```

### Feasibility Analysis

Before running expensive DFT calculations, analyze which Miller indices have valid surfaces:

```python
from teros.core.surface_energy import analyze_miller_feasibility, get_feasibility_summary

reports = analyze_miller_feasibility(
    structure,
    miller_indices=[(1,0,0), (1,1,0), (1,1,1)],
)

print(get_feasibility_summary(reports))

# Or check individually
for miller, report in reports.items():
    if report.has_valid_surfaces:
        print(f"{miller}: OK - use simple formula")
    else:
        print(f"{miller}: NEEDS THERMODYNAMICS")
```

### WorkGraph Integration

Enable stoichiometric filtering in the main workflow:

```python
wg = build_metal_surface_energy_workgraph(
    bulk_structure_path='/path/to/pdin.cif',
    miller_indices=[[1,1,0], [1,0,0]],

    # EXPERIMENTAL: Only generate stoichiometric+symmetric slabs
    require_stoichiometric_symmetric=True,
    stoichiometric_strategies=['filter_first', 'thickness_scan'],
    stoichiometric_max_thickness=30.0,

    code_label='VASP-6.5.1@cluster',
    potential_mapping={'Pd': 'Pd', 'In': 'In'},
    ...
)
```

### Bond Preservation for Oxides

For materials with polyhedral units (phosphates, sulfates, etc.):

```python
results = find_stoichiometric_symmetric_slabs(
    ag3po4_structure,
    miller_index=(1, 0, 0),
    strategies=['bond_preservation', 'thickness_scan'],
    bonds={('P', 'O'): 1.9},  # Preserve PO4 tetrahedra
)

# Or auto-detect bonds
results = find_stoichiometric_symmetric_slabs(
    structure,
    miller_index=(1, 0, 0),
    auto_detect_bonds=True,
)
```

### Error Handling

When no valid surface is found, a descriptive error is raised:

```python
from teros.core.surface_energy import (
    find_stoichiometric_symmetric_slabs,
    NoStoichiometricSymmetricSurfaceError,
)

try:
    results = find_stoichiometric_symmetric_slabs(structure, (1, 1, 1))
except NoStoichiometricSymmetricSurfaceError as e:
    print(e)
    # Output includes:
    # - Miller index
    # - Strategies tried
    # - Statistics (terminations checked, stoichiometric only, symmetric only)
    # - Recommendation to use thermodynamics approach
```

### Data Classes

#### `SlabSearchResult`
```python
@dataclass
class SlabSearchResult:
    slab: Optional[Slab]           # PyMatGen Slab or None
    is_stoichiometric: bool        # Matches bulk composition
    is_symmetric: bool             # Equivalent top/bottom
    strategy_used: str             # Which strategy found it
    thickness_angstrom: float      # Slab thickness
    termination_index: int         # Termination number
    miller_index: tuple            # Miller index
    bonds_broken: int = 0          # For bond_preservation
    warnings: list[str] = []       # Any warnings
```

#### `MillerFeasibilityReport`
```python
@dataclass
class MillerFeasibilityReport:
    miller_index: tuple            # Miller index analyzed
    has_valid_surfaces: bool       # Any valid surfaces found?
    n_terminations_checked: int    # Total terminations
    n_stoichiometric: int          # Stoichiometric only
    n_symmetric: int               # Symmetric only
    n_both: int                    # Both criteria met
    best_thickness: float          # Recommended thickness
    recommended_strategy: str      # Best strategy
    notes: list[str]               # Recommendations
```

### Wulff Shape Filtering

Filter Wulff shape to only include stoichiometric surfaces:

```python
from teros.core.surface_energy import build_wulff_shape

wulff_result = build_wulff_shape(
    bulk_structure,
    surface_energies,
    only_stoichiometric_symmetric=orm.Bool(True),  # Filter non-stoichiometric
)
```

### Module Structure

```
surface_energy/
├── __init__.py               # Module exports
├── surface_energy.py         # Surface energy calcfunctions
├── workgraph.py              # WorkGraph builder
├── wulff.py                  # Wulff shape construction
├── stoichiometric_finder.py  # Stoichiometric+symmetric finder (NEW)
└── README.md                 # This file
```

### Limitations

1. **Search completeness**: Not all possible thicknesses/terminations can be checked
2. **Bond detection**: Auto-detection may miss some bonds; provide explicit bonds dict for accuracy
3. **Complex structures**: Some materials may have no stoichiometric+symmetric surfaces for certain Miller indices

### When to Use Thermodynamics Instead

Use the thermodynamics module (`teros.core.thermodynamics`) when:
- The feasibility analysis shows no valid surfaces
- You need surface energies for non-stoichiometric terminations
- You want γ(Δμ) as a function of chemical potential

---

## Dependencies

- `pymatgen`: WulffShape, SpacegroupAnalyzer, SlabGenerator, structure conversion
- `aiida-workgraph`: Workflow orchestration
- `aiida-vasp`: VASP calculations
- `matplotlib`: Visualization (optional)
- `numpy`: Numerical operations

## See Also

- `teros.core.thermodynamics` - Oxide surface thermodynamics with chemical potentials
- `teros.core.slabs` - Slab generation utilities
- [Pymatgen WulffShape Tutorial](https://matgenb.materialsvirtuallab.org/2017/04/03/Slab-generation-and-Wulff-shape.html)
