# Metal Surface Energy Module

This module computes surface energies for **elemental metals** (Au, Ag, Cu, Pt, etc.) using a simple formula without chemical potential dependencies.

## Formula

$$\gamma = \frac{E_{\text{slab}} - N \cdot E_{\text{bulk}}^{\text{atom}}}{2A}$$

- $E_{\text{slab}}$: Total DFT energy of the slab (eV)
- $N$: Number of atoms in the slab
- $E_{\text{bulk}}^{\text{atom}}$: Bulk energy per atom (eV)
- $A$: Surface area (Å²)
- Factor of 2: Two surfaces per slab

## Quick Start

```python
from aiida import load_profile
from teros.core.metal_surface_energy import build_metal_surface_energy_workgraph

load_profile('myprofile')

wg = build_metal_surface_energy_workgraph(
    bulk_structure_path='/path/to/au.cif',
    code_label='VASP-6.5.1@cluster',
    potential_family='PBE',
    potential_mapping={'Au': 'Au'},
    
    # Multiple orientations in one workflow
    miller_indices=[[1,1,1], [1,0,0], [1,1,0]],
    
    # Slab parameters
    min_slab_thickness=20.0,
    min_vacuum_thickness=20.0,
    
    # VASP parameters
    bulk_parameters={
        'prec': 'Accurate',
        'encut': 500,
        'ismear': 1,      # Methfessel-Paxton for metals
        'sigma': 0.2,
        'ibrion': 2,
        'isif': 3,        # Full relaxation for bulk
        'nsw': 100,
        'ediffg': -0.01,
    },
    
    name='Au_surface_energy',
)

wg.submit(wait=False)
```

## Workflow Architecture

```
WorkGraph<MetalSurfaceEnergy>
├── bulk_relax (VaspWorkChain)     ← Runs ONCE, shared by all orientations
├── bulk_energy
├── bulk_energy_per_atom
├── surface_hkl_111                ← Per-orientation sub-workflows
│   ├── generate_slab_structures
│   ├── relax_slabs_scatter
│   └── compute_surface_energies
├── surface_hkl_100
├── surface_hkl_110
└── gather_surface_energies        ← Consolidates all results
```

## Key Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `bulk_structure_path` | Path to bulk CIF file | `/path/to/au.cif` |
| `miller_indices` | List of orientations | `[[1,1,1], [1,0,0]]` |
| `potential_mapping` | Element to potential | `{'Au': 'Au'}` |
| `min_slab_thickness` | Slab thickness (Å) | `20.0` |
| `min_vacuum_thickness` | Vacuum gap (Å) | `20.0` |
| `kpoints_spacing` | K-mesh spacing | `0.02` |

## VASP Recommendations for Metals

```python
bulk_parameters = {
    'ismear': 1,      # Methfessel-Paxton (NOT Gaussian for metals!)
    'sigma': 0.2,     # Smearing width
    'encut': 500,     # ≥1.3 × ENMAX
    'ediff': 1e-6,
    'ibrion': 2,
    'isif': 3,        # Full relaxation for bulk
    'nsw': 100,
}

slab_parameters = bulk_parameters.copy()
slab_parameters['isif'] = 2  # Fix cell, relax ions only
```

## Outputs

After completion:

```
Outputs                  PK    Type
-----------------------  ----  -------------
bulk_energy              XXXX  Float
bulk_structure           XXXX  StructureData
surface_energies         XXXX  Dict          ← All results
relaxed_slabs_hkl_111
    term_0               XXXX  StructureData
...
```

The `surface_energies` Dict:
```json
{
    "hkl_111": {
        "term_0": {
            "gamma_J_m2": 0.77,
            "gamma_eV_A2": 0.048,
            "area_A2": 21.5,
            "N_slab": 8
        }
    },
    "hkl_100": { ... },
    "hkl_110": { ... }
}
```

## Expected Values (PBE)

| Metal | (111) J/m² | (100) J/m² | (110) J/m² |
|-------|------------|------------|------------|
| Au    | ~0.8       | ~0.9       | ~1.3       |
| Ag    | ~0.6       | ~0.7       | ~0.9       |
| Cu    | ~1.3       | ~1.5       | ~1.6       |
| Pt    | ~1.5       | ~1.9       | ~2.0       |

## Comparison with Oxide Module

| Feature | Metal Module | Oxide Thermodynamics |
|---------|--------------|---------------------|
| Formula | γ = (E_slab - N·E_bulk/atom)/(2A) | γ(Δμ_O, Δμ_M) |
| Chemical potential | None | Required |
| Reference calcs | None | O₂, metal refs |
| Use case | Pure metals | Oxides (Ag₂O, WO₃, etc.) |

## Files

- `surface_energy.py` - Core calcfunctions
- `workgraph.py` - WorkGraph builder
- `__init__.py` - Module exports

## See Also

- [Workflow Presets Guide](WORKFLOW_PRESETS_GUIDE.md)
- [Adsorption Energy Module](adsorption_energy_module.md)
- Example: `examples/surface_thermo_gold/surface_thermo_metals.py`
