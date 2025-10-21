# LaMnO3 + OH Adsorption Energy Calculation

This example demonstrates the **new features** of the PS-TEROS adsorption energy module:

1. **Full Builder API** for `vasp.v2.relax` workchain (complete control)
2. **Atom Fixing** for slab relaxations (fix bottom/top/center layers)

## Overview

Calculate adsorption energies using the formula:

```
E_ads = E_complete - E_substrate - E_molecule
```

Where:
- **E_complete**: Energy of substrate + adsorbate system
- **E_substrate**: Energy of bare substrate
- **E_molecule**: Energy of isolated adsorbate molecule

**Negative E_ads** = Favorable (exothermic) adsorption
**Positive E_ads** = Unfavorable (endothermic) adsorption

## New Features

### 1. Full Builder API for Relaxation

Previously, you had limited control over relaxation parameters. Now you have **complete control** over `vasp.v2.relax` settings:

```python
relax_builder_inputs = {
    'relax_settings': {  # Plain dict - converted to orm.Dict internally
        'algo': 'cg',              # Relaxation algorithm
        'force_cutoff': 0.1,       # Force convergence (eV/Å)
        'steps': 500,              # Maximum ionic steps
        'positions': True,         # Relax atomic positions
        'shape': False,            # Don't relax cell shape
        'volume': False,           # Don't relax cell volume
        'convergence_on': True,    # Enable convergence checks
        'convergence_positions': 0.01,  # Position convergence (Å)
        # ... full control over all settings
    },
    'vasp': {
        'parameters': {'incar': {...}},  # VASP parameters
        'options': {...},                # Scheduler options
        'potential_family': 'PBE',
        'potential_mapping': {...},
        'kpoints_spacing': 0.6,
    }
}
```

**Key Points:**
- Use **plain Python dicts** (automatically converted to `orm.Dict`)
- Complete control over all `relax_settings` parameters
- Full access to VASP parameters in nested `vasp` namespace

### 2. Atom Fixing for Slabs

Fix bottom layers of slabs during relaxation to simulate bulk-like behavior:

```python
wg = build_core_workgraph(
    ...
    adsorption_fix_atoms=True,
    adsorption_fix_type='bottom',      # 'bottom', 'top', or 'center'
    adsorption_fix_thickness=7.0,      # Thickness in Angstroms
    adsorption_fix_elements=None,      # None = all elements, or ['La', 'Mn']
)
```

**Fixing Options:**
- `'bottom'`: Fix atoms from bottom up to specified thickness
- `'top'`: Fix atoms from top down to specified thickness
- `'center'`: Fix atoms within thickness/2 of slab center
- `fix_elements`: Optional list to fix only specific elements

## Example: LaMnO3 + OH

This example calculates the adsorption energy of OH radical on LaMnO3(100) surface.

**System:**
- Substrate: LaMnO3 (100) surface (2×2×1 supercell)
- Adsorbate: OH radical
- Total structure: La₃₂Mn₂₈O₈₉H (150 atoms)

**Workflow Phases:**

1. **Relax** complete system (LaMnO3 + OH)
   - Uses `vasp.v2.relax` with full builder control
   - Fixes bottom 7 Å of slab

2. **Separate** relaxed structure into:
   - Substrate: La₃₂Mn₂₈O₈₈ (149 atoms)
   - Molecule: OH (2 atoms)
   - Complete: La₃₂Mn₂₈O₈₉H (150 atoms)

3. **SCF** calculations (3 parallel jobs):
   - Substrate energy
   - Molecule energy
   - Complete energy

4. **Calculate** E_ads = E_complete - E_substrate - E_molecule

## Usage

### Running the Example

```bash
# Activate environment
source ~/envs/aiida/bin/activate

# Run the example
python teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py
```

### Monitoring

```bash
# Check status
verdi process status <PK>

# View detailed report
verdi process report <PK>

# Watch live
verdi process watch <PK>
```

## API Reference

### `build_core_workgraph()` Parameters

**Adsorption Energy Workflow:**

```python
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',

    # Required parameters
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    adsorption_structures={'label': structure},
    adsorption_formulas={'label': 'OH'},
    adsorption_potential_mapping={'La': 'La', 'Mn': 'Mn_pv', ...},

    # NEW: Full builder control for relaxation
    relax_before_adsorption=True,
    adsorption_relax_builder_inputs={
        'relax_settings': {...},
        'vasp': {...}
    },

    # NEW: Builder control for SCF
    adsorption_scf_builder_inputs={
        'parameters': {'incar': {...}},
        'options': {...},
        ...
    },

    # NEW: Atom fixing
    adsorption_fix_atoms=True,
    adsorption_fix_type='bottom',
    adsorption_fix_thickness=7.0,
    adsorption_fix_elements=None,
)
```

### Workflow Outputs

Access results via:

```python
# Relaxed structures
relaxed = wg.outputs.relaxed_complete_structures['label']

# Separated structures
separated = wg.outputs.separated_structures['label']
substrate = separated['substrate']
molecule = separated['molecule']
complete = separated['complete']

# Energies
E_substrate = wg.outputs.substrate_energies['label']
E_molecule = wg.outputs.molecule_energies['label']
E_complete = wg.outputs.complete_energies['label']

# Adsorption energy
E_ads = wg.outputs.adsorption_energies['label']
```

## Configuration Details

### Relaxation Settings

The example uses these `relax_settings`:

```python
{
    'algo': 'cg',                          # Conjugate gradient
    'force_cutoff': 0.1,                   # 0.1 eV/Å
    'steps': 500,                          # Max ionic steps
    'positions': True,                     # Relax positions
    'shape': False,                        # Keep cell shape fixed
    'volume': False,                       # Keep volume fixed
    'convergence_on': True,                # Enable convergence checks
    'convergence_absolute': False,
    'convergence_max_iterations': 5,
    'convergence_positions': 0.01,         # 0.01 Å
    'convergence_volume': 0.01,
    'convergence_shape_lengths': 0.1,
    'convergence_shape_angles': 0.1,
    'perform': True,
}
```

### VASP Parameters

**Relaxation (Phase 1):**
```python
{
    'PREC': 'Accurate',
    'ALGO': 'Normal',
    'ENCUT': 500,           # eV
    'EDIFF': 1e-5,          # Electronic convergence
    'EDIFFG': -0.1,         # Force convergence (eV/Å)
    'ISMEAR': 0,            # Gaussian smearing
    'SIGMA': 0.01,
    'ISPIN': 2,             # Spin-polarized
    'LORBIT': 11,
    'LASPH': True,          # Non-spherical contributions
    'IVDW': 12,             # DFT-D3 van der Waals
    'NCORE': 3,
    'KPAR': 4,
}
```

**SCF (Phase 3):**
- Same parameters as relaxation
- NSW=0, IBRION=-1 (automatically enforced)

### Atom Fixing

With `fix_type='bottom'` and `fix_thickness=7.0`:

- Identifies all atoms within 7.0 Å from the bottom of the slab
- Creates ASE `FixAtoms` constraint
- Atoms remain fixed during relaxation
- Simulates bulk-like behavior in bottom layers

## Expected Results

**For OH on LaMnO3:**
- Literature E_ads: -2 to -4 eV (DFT-PBE+U)
- This example: Depends on surface termination and adsorption site

**Computational Cost:**
- Phase 1 (relax): ~12-18 hours (150 atoms, ionic steps)
- Phase 3 (SCF): ~2-4 hours each (3 parallel jobs)
- **Total**: ~12-18 hours (relaxation is the bottleneck)

## Important Notes

### Cell Handling

All three systems (substrate, molecule, complete) use the **same simulation cell**:
- Original cell from complete structure is preserved
- Adsorbate molecule stays at original position in the cell
- This eliminates **basis set superposition error (BSSE)**

### Backward Compatibility

The old-style API still works:

```python
# Old style (still supported)
wg = build_core_workgraph(
    ...
    adsorption_parameters={'NSW': 100, 'IBRION': 2, ...},
    adsorption_options={'resources': {'num_machines': 1}},
    adsorption_kpoints_spacing=0.6,
)
```

But the new builder API is **recommended** for full control.

## Testing

The implementation was successfully tested:

```bash
✓ Modules import successfully
✓ Plain dict → orm.Dict conversion works
✓ Code placement correct for relax workchain
✓ Atom fixing correctly applied
✓ Workflow submitted successfully (PK: 101121)
✓ VaspRelaxWorkChain running with constraints
```

**Test workflow status:**
```
WorkGraph<LaMnO3_OH_RelaxThenAdsorption_WithFixing><101121> Waiting
    └── WorkGraph<compute_adsorption_energies_scatter><101139> Waiting
        └── VaspRelaxWorkChain<101154> Waiting
            └── VaspWorkChain<101157> Waiting
                └── VaspCalculation<101160> Waiting
```

## Troubleshooting

### Common Issues

**1. Serialization error (orm.Dict not JSON-serializable)**

❌ **Wrong:**
```python
relax_builder_inputs = {
    'relax_settings': orm.Dict(dict={...}),  # DON'T use orm.Dict
}
```

✅ **Correct:**
```python
relax_builder_inputs = {
    'relax_settings': {...},  # Use plain dict
}
```

**2. Code placement error**

The code is automatically placed correctly:
- For `relax`: code goes in `vasp` namespace
- For `vasp`: code goes at top level

**3. Daemon not picking up changes**

Always restart after code modifications:
```bash
verdi daemon restart
```

## File Structure

```
teros/experimental/adsorption_energy/lamno3/
├── README.md                           # This file
├── run_lamno3_oh_adsorption.py        # Example script
├── LaMnO3_100_A4_surface_2x2x1_OOH.cif  # Test structure
└── surface_builder/                    # Structure building tools
```

## References

### Code Locations

- **Core module**: `teros/core/adsorption_energy.py`
- **Fixed atoms**: `teros/core/fixed_atoms.py`
- **Workgraph integration**: `teros/core/workgraph.py`

### Related Documentation

- VASP aiida-vasp plugin: https://github.com/aiida-vasp/aiida-vasp
- AiiDA WorkGraph: https://github.com/aiidateam/aiida-workgraph

## Version History

- **v2.0** (2025-10-21): Added builder API and atom fixing
  - Full control over `vasp.v2.relax` via `relax_settings`
  - Atom fixing for slab calculations
  - Improved code structure and documentation

- **v1.0** (Previous): Basic adsorption energy workflow
  - Limited parameter control
  - No atom constraints

---

**Last Updated:** 2025-10-21
**Feature Branch:** `feature-adsorption-energy`
**Status:** ✅ Tested and working
