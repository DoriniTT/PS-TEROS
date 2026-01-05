# Metal Surface Energy Module

This module computes surface energies for **elemental metals** using the simple thermodynamic formula:

$$\gamma = \frac{E_{\text{slab}} - N \cdot E_{\text{bulk}}^{\text{atom}}}{2A}$$

Where:
- $E_{\text{slab}}$: Total DFT energy of the slab
- $N$: Number of atoms in the slab
- $E_{\text{bulk}}^{\text{atom}}$: Bulk energy per atom
- $A$: Surface area (Å²)
- Factor of 2 accounts for two surfaces (top and bottom)

## Key Differences from Oxide Thermodynamics

| Feature | Metal Surface Energy | Oxide Thermodynamics |
|---------|---------------------|----------------------|
| Chemical potential | None | Δμ_O, Δμ_Metal dependencies |
| Formula | Simple single-value γ | γ(Δμ) function |
| Reference calculations | None needed | O₂, metal references required |
| Complexity | Low | High |

## Module Structure

```
metal_surface_energy/
├── __init__.py              # Module exports
├── surface_energy.py        # Surface energy calcfunctions
├── workgraph.py             # WorkGraph builder
└── README.md                # This file
```

## Key Functions

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
- `area_A2`, `N_slab`, `N_bulk`, `element`

### `build_metal_surface_energy_workgraph`
Main entry point for users. Builds a complete workflow.

**Key Parameters:**
- `bulk_structure_path`: Path to CIF file
- `miller_indices`: List of orientations, e.g., `[[1,1,1], [1,0,0], [1,1,0]]`
- `code_label`: VASP code in AiiDA
- `potential_mapping`: e.g., `{'Au': 'Au'}`

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
└── gather_surface_energies             # Consolidates all results
```

## Usage Example

```python
from aiida import load_profile
from teros.core.metal_surface_energy import build_metal_surface_energy_workgraph

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

## Output Structure

After completion, `verdi process show <PK>` displays:

```
Outputs                  PK    Type
-----------------------  ----  -------------
bulk_energy              XXXX  Float
bulk_energy_per_atom     XXXX  Float
bulk_structure           XXXX  StructureData
surface_energies         XXXX  Dict          ← All results consolidated
relaxed_slabs_hkl_111
    term_0               XXXX  StructureData
relaxed_slabs_hkl_100
    term_0               XXXX  StructureData
...
```

The `surface_energies` Dict contains:
```json
{
    "hkl_111": {
        "term_0": {
            "gamma_eV_A2": 0.048,
            "gamma_J_m2": 0.77,
            "area_A2": 21.5,
            "N_slab": 8,
            "element": "Au"
        }
    },
    "hkl_100": { ... },
    "hkl_110": { ... }
}
```

## Expected Surface Energies (PBE)

| Metal | (111) | (100) | (110) |
|-------|-------|-------|-------|
| Au    | ~0.8  | ~0.9  | ~1.3  |
| Ag    | ~0.6  | ~0.7  | ~0.9  |
| Cu    | ~1.3  | ~1.5  | ~1.6  |
| Pt    | ~1.5  | ~1.9  | ~2.0  |

*Values in J/m². Actual results depend on computational parameters.*

## VASP Parameter Recommendations

For metals, use:
- `ISMEAR = 1` (Methfessel-Paxton) or `ISMEAR = 2`
- `SIGMA = 0.1-0.2` eV
- `ENCUT ≥ 1.3 × ENMAX` from POTCAR
- Dense k-points: spacing ≤ 0.03 Å⁻¹

For slabs:
- `ISIF = 2` (relax ions only, fix cell)
- Sufficient vacuum (≥15 Å)
- Consider dipole corrections if asymmetric
