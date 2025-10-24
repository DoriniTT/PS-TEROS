# Adsorption Energy Module

## Overview

The adsorption energy module (`teros.core.adsorption_energy`) provides tools to calculate adsorption energies from substrate+adsorbate structures using AiiDA-WorkGraph and VASP. The module uses advanced connectivity analysis to automatically separate structures and performs parallel DFT calculations to obtain adsorption energies.

## ✅ Key Improvement: Handles Bonded Adsorbates

**New subgraph algorithm** allows automatic structure separation even when adsorbates are **bonded to the surface**:

- ✅ **Works with chemisorbed species** (strong covalent bonds to surface)
- ✅ **Works with physisorbed species** (weak van der Waals interactions)
- ✅ **No geometric assumptions** needed about bonding configuration
- ✅ **Tested with OH on Ag(111)** - both hollow (weak) and top (strong) sites

**How:** Builds connectivity subgraph with edges **only between adsorbate atoms**, ignoring bonds to substrate. This allows identification of bonded clusters like OH, OOH, COOH regardless of surface attachment.

## Workflow Modes

### Mode 1: Direct SCF (Default)

**When to use:** Input structures are pre-relaxed or well-optimized

**Workflow:**
```
Input structure → Separate → SCF (3 jobs) → E_ads
```

**Example:**
```python
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={'system': structure},
    adsorption_formulas={'system': 'OH'},
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    adsorption_options={'resources': {'num_machines': 1}},
)
```

**VASP jobs:** 3N (substrate SCF, molecule SCF, complete SCF for each structure)

---

### Mode 2: Relax + SCF (New)

**When to use:** Input structures need geometry optimization before energy calculation

**Workflow:**
```
Input structure → Relax → Separate relaxed → SCF (3 jobs) → E_ads
```

**Example (Simplified API - Direct INCAR Parameters):**
```python
# Relaxation parameters (NSW > 0 for relaxation)
relax_params = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'IBRION': 2,
    'NSW': 100,      # NSW > 0 = relaxation
    'ISIF': 2,
    'EDIFFG': -0.05,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'LWAVE': False,
    'LCHARG': False,
}

# SCF parameters (NSW=0, IBRION=-1 set automatically)
scf_params = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'NELM': 200,
    'LWAVE': False,
    'LCHARG': False,
}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',

    # Structures
    adsorption_structures={'system': structure},
    adsorption_formulas={'system': 'OH'},
    adsorption_potential_mapping={'Ag': 'Ag', 'O': 'O', 'H': 'H'},

    # Code
    code_label='VASP@cluster',
    potential_family='PBE',

    # Simplified API: Direct INCAR parameters
    relax_before_adsorption=True,
    adsorption_relax_builder_inputs={'parameters': {'incar': relax_params}},
    adsorption_scf_builder_inputs={'parameters': {'incar': scf_params}},

    # Scheduler
    adsorption_options={'resources': {'num_machines': 1}},
    adsorption_kpoints_spacing=0.3,
)

# Access relaxed structure
relaxed = wg.outputs.relaxed_complete_structures['system']
```

**VASP jobs:** N + 3N (relax complete for each, then substrate/molecule/complete SCF)

**Key advantages:**
- Geometrically correct: Separates RELAXED structure, not input structure
- Relaxes only complete system (efficient)
- SCF is fast for all three components
- Eliminates geometry-related errors in E_ads

---

### Mode 3: Structure-Specific Parameter Overrides (NEW)

**When to use:** Different structures need different VASP settings (e.g., failed calculations, difficult convergence, or precision requirements vary by system)

**NEW FEATURE**: Override parameters for specific structures using deep merge strategy. This allows you to:
- Rerun only failed calculations with custom settings
- Use tighter convergence for specific structures
- Adjust resources for computationally demanding systems
- Test different algorithms for problematic structures

**How it works:**
- Structure indices (0, 1, 2, ...) match the order in `adsorption_structures` dict
- Deep merge: Only specified parameters are overridden, others inherited from defaults
- Applies to both relaxation and SCF phases independently
- All three SCF calculations (substrate, molecule, complete) use the same overrides

**Example:**
```python
# Default parameters for all structures
relax_params = {
    'PREC': 'Normal',
    'ENCUT': 400,
    'ALGO': 'Fast',
    'EDIFFG': -0.02,
    'NSW': 100,
    'IBRION': 2,
    'ISIF': 2,
}

scf_params = {
    'PREC': 'Normal',
    'ENCUT': 400,
    'EDIFF': 1e-6,
}

# Override for structure 0 (LaNiO3): Use higher precision
structure_specific_relax = {
    0: {  # Structure index 0
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 600,
                'ALGO': 'Normal',
            }
        },
        'kpoints_spacing': 0.8,  # Denser k-points
    }
}

# Override for structure 1 (LaCoO3): Tighter convergence + more resources
structure_specific_scf = {
    1: {  # Structure index 1
        'parameters': {
            'incar': {
                'EDIFF': 1e-7,
                'PREC': 'Accurate',
            }
        },
        'kpoints_spacing': 0.15,  # Very dense k-points for final energy
        'options': {
            'resources': {
                'num_machines': 2,  # More resources
            }
        }
    }
}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={
        'lanio3_oh': structure0,  # Index 0
        'lacoo3_oh': structure1,  # Index 1
    },
    adsorption_formulas={
        'lanio3_oh': 'OH',
        'lacoo3_oh': 'OH',
    },
    code_label='VASP@cluster',
    potential_family='PBE',
    adsorption_potential_mapping={'La': 'La', 'Ni': 'Ni', 'Co': 'Co', 'O': 'O', 'H': 'H'},

    # Default parameters (all structures)
    relax_before_adsorption=True,
    adsorption_relax_builder_inputs={'parameters': {'incar': relax_params}},
    adsorption_scf_builder_inputs={'parameters': {'incar': scf_params}},
    adsorption_options={'resources': {'num_machines': 1}},
    adsorption_kpoints_spacing=1.0,

    # NEW: Structure-specific overrides
    adsorption_structure_specific_relax_builder_inputs=structure_specific_relax,
    adsorption_structure_specific_scf_builder_inputs=structure_specific_scf,
)
```

**Result after deep merge:**
- Structure 0 (lanio3_oh) relaxation: Uses PREC='Accurate', ENCUT=600, ALGO='Normal', kpoints=0.8, but inherits EDIFFG, NSW, IBRION, ISIF from defaults
- Structure 0 SCF: Uses default scf_params
- Structure 1 (lacoo3_oh) relaxation: Uses default relax_params
- Structure 1 SCF: Uses EDIFF=1e-7, PREC='Accurate', kpoints=0.15, num_machines=2, but inherits other parameters from defaults
- All other structures: Use defaults for both relax and SCF

**Common use cases:**

1. **Rerun failed calculations:**
```python
# After first run, identify failed structures (e.g., indices 2, 5)
structure_specific_relax = {
    2: {'parameters': {'incar': {'ALGO': 'All', 'EDIFFG': -0.05}}},
    5: {'parameters': {'incar': {'IBRION': 1, 'NSW': 200}}},
}
# All structures will be processed, but only 2 and 5 will use custom settings
```

2. **Different precision levels:**
```python
# Use high precision for critical structures
structure_specific_scf = {
    0: {'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 600}}, 'kpoints_spacing': 0.15},
    1: {'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 600}}, 'kpoints_spacing': 0.15},
}
```

3. **Algorithm testing:**
```python
# Test different algorithms for problematic structures
structure_specific_relax = {
    3: {
        'parameters': {
            'incar': {
                'ALGO': 'All',      # Try all algorithms
                'IBRION': 1,        # RMM-DIIS
                'POTIM': 0.2,       # Smaller ionic step
                'EDIFFG': -0.05,    # Looser convergence
            }
        },
        'options': {'max_wallclock_seconds': 10800},  # More time
    }
}
```

**Important notes:**
- Structure indices are determined by insertion order into `adsorption_structures` dict
- Indices are 0-based (first structure = 0, second = 1, etc.)
- Keys can be integers or strings ('0', '1', etc.) - both are accepted
- Deep merge means only specified parameters override defaults
- This is the **recommended pattern** for all PS-TEROS modules going forward

---

## API Reference

### New Parameters (build_core_workgraph)

**`relax_before_adsorption`** (bool, default=False)
- If True, relax complete structures before separation and SCF
- Uses only `vasp.v2.vasp` WorkflowFactory (not `vasp.v2.relax`)
- Relaxation controlled by NSW parameter (NSW > 0)

**`adsorption_relax_builder_inputs`** (dict, optional)
- **Simplified API**: Pass INCAR parameters directly wrapped in nested dict
- Structure: `{'parameters': {'incar': {INCAR_PARAMS}}}`
- Must include: NSW > 0 (for relaxation), IBRION, ISIF, EDIFFG
- Example: `{'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'EDIFFG': -0.05}}}`
- Note: Common parameters (code, potential_family, options) are set separately

**`adsorption_scf_builder_inputs`** (dict, optional)
- **Simplified API**: Pass INCAR parameters directly wrapped in nested dict
- Structure: `{'parameters': {'incar': {INCAR_PARAMS}}}`
- NSW=0 and IBRION=-1 are **automatically enforced** (don't need to specify)
- Example: `{'parameters': {'incar': {'ENCUT': 520, 'EDIFF': 1e-6, 'NELM': 200}}}`
- Note: Common parameters (code, potential_family, options) are set separately

**`adsorption_structure_specific_relax_builder_inputs`** (dict, optional) - NEW
- Per-structure overrides for relaxation parameters using deep merge
- Structure: Dict mapping structure indices (0, 1, 2, ...) to builder dicts
- Example: `{0: {'parameters': {'incar': {'ALGO': 'Normal'}}, 'kpoints_spacing': 0.2}}`
- Only specified parameters are overridden; others inherited from `adsorption_relax_builder_inputs`
- Applies to relaxation of complete (substrate+adsorbate) structure only
- Integer keys automatically converted to strings for AiiDA compatibility

**`adsorption_structure_specific_scf_builder_inputs`** (dict, optional) - NEW
- Per-structure overrides for SCF parameters using deep merge
- Structure: Dict mapping structure indices (0, 1, 2, ...) to builder dicts
- Example: `{1: {'parameters': {'incar': {'EDIFF': 1e-7}}, 'options': {'resources': {'num_machines': 2}}}}`
- Only specified parameters are overridden; others inherited from `adsorption_scf_builder_inputs`
- Applies to ALL THREE SCF calculations (substrate, molecule, complete) for that structure
- Integer keys automatically converted to strings for AiiDA compatibility

**Key Change from Previous Version:**
- Old: Complex builder structure with multiple namespaces
- New: Direct INCAR parameters via simple wrapped dict
- Plugin: Uses only `vasp.v2.vasp` (not `vasp.v2.relax`)
- Control: Relaxation vs SCF determined by NSW parameter
- **Latest**: Structure-specific builder inputs for fine-grained control (v2.1+)

### Backward Compatibility

The simplified API can coexist with the old parameter style:

**Old style (deprecated but still works):**
```python
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={...},
    adsorption_formulas={...},
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'NSW': 100},
    adsorption_options={...},
    adsorption_potential_mapping={...},
)
```

**New simplified API (recommended):**
```python
relax_params = {'PREC': 'Accurate', 'ENCUT': 520, 'NSW': 100, 'IBRION': 2, ...}
scf_params = {'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6, 'NELM': 200, ...}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={...},
    adsorption_formulas={...},
    adsorption_potential_mapping={...},

    relax_before_adsorption=True,
    adsorption_relax_builder_inputs={'parameters': {'incar': relax_params}},
    adsorption_scf_builder_inputs={'parameters': {'incar': scf_params}},

    adsorption_options={...},
)
```

Migration is recommended for new workflows. Old scripts continue to work.

---

## Testing

### Test Cases

**Test 1: Simplified API with relaxation**
- File: `examples/vasp/step_12_adsorption_energy.py`
- Tests: New simplified API with direct INCAR parameters and relaxation enabled
- Plugin: Uses only `vasp.v2.vasp` (not `vasp.v2.relax`)

**Test 2: Production example (perovskite oxides)**
- Location: User workflows (e.g., LaNiO3, LaMnO3 with OOH adsorption)
- Tests: Real-world cases with relaxation on complex oxide surfaces

### Validation Checklist

- [ ] Simplified API works correctly (relax_builder_inputs + scf_builder_inputs)
- [ ] Uses only vasp.v2.vasp plugin (verify in WorkGraph)
- [ ] Relaxation completes successfully (NSW > 0)
- [ ] Separation works on relaxed structures
- [ ] SCF calculations automatically have NSW=0, IBRION=-1
- [ ] E_ads values are physically reasonable
- [ ] Relaxed structures show expected geometry changes

### Key API Changes

**Before (complex builder):**
- Used `vasp.v2.relax` for relaxation
- Complex nested builder structure with multiple namespaces
- Required separate parameters for each namespace

**After (simplified):**
- Uses only `vasp.v2.vasp` for all calculations
- Direct INCAR parameters via `{'parameters': {'incar': {...}}}`
- Common settings (code, options) set separately
- NSW parameter controls relaxation vs SCF automatically

## Theory

The adsorption energy is calculated using the standard formula:

```
E_ads = E_complete - E_substrate - E_molecule
```

Where:
- `E_complete`: Total energy of the substrate+adsorbate system (eV)
- `E_substrate`: Total energy of the bare substrate (eV)
- `E_molecule`: Total energy of the isolated adsorbate molecule (eV)

**Interpretation:**
- **Negative E_ads**: Exothermic adsorption (favorable, stable)
- **Positive E_ads**: Endothermic adsorption (unfavorable, unstable)

## Key Features

- **✅ Automatic structure separation** using advanced connectivity analysis (pymatgen StructureGraph)
- **✅ Handles bonded adsorbates** - Works with both chemisorbed and physisorbed species
- **Parallel VASP calculations** for substrate, molecule, and complete system via scatter-gather
- **Consistent cell handling** - all three systems use the same cell to avoid basis set superposition error
- **Integration with existing PSTEROS modules** - follows same patterns as thermodynamics/cleavage
- **Robust adsorbate identification** - no geometric assumptions or distance thresholds needed
- **Workflow preset integration** - Use `workflow_preset='adsorption_energy'` in `build_core_workgraph`

## Core Functions

### `separate_adsorbate_structure()`

Separates substrate+adsorbate into three components using connectivity analysis.

**Purpose:** Identifies and extracts the adsorbate from a complete structure by analyzing chemical bonds.

**Parameters:**
- `structure` (StructureData): Complete substrate+adsorbate structure
- `adsorbate_formula` (Str): Chemical formula of adsorbate (e.g., "OOH", "OH")

**Returns:**
Dictionary with three StructureData nodes:
- `substrate`: Structure with adsorbate removed, original cell preserved
- `molecule`: Adsorbate in original cell at original position
- `complete`: Original structure (for provenance)

**Raises:**
- `ValueError`: If formula is invalid, adsorbate not found, or multiple matches

**Example:**
```python
from teros.core import separate_adsorbate_structure
from aiida import orm

result = separate_adsorbate_structure(
    structure=complete_structure,
    adsorbate_formula=orm.Str('OH')
)

substrate = result['substrate']  # Bare surface
molecule = result['molecule']    # Isolated OH
complete = result['complete']    # Original structure
```

### `calculate_adsorption_energy()`

Calculates adsorption energy from component energies.

**Purpose:** Applies the adsorption energy formula to get E_ads from relaxed energies.

**Formula:** E_ads = E_complete - E_substrate - E_molecule

**Parameters:**
- `E_complete` (Float): Energy of substrate+adsorbate system (eV)
- `E_substrate` (Float): Energy of bare substrate (eV)
- `E_molecule` (Float): Energy of isolated molecule (eV)

**Returns:**
- Float: Adsorption energy in eV (negative = favorable)

**Example:**
```python
from teros.core import calculate_adsorption_energy
from aiida import orm

E_ads = calculate_adsorption_energy(
    E_complete=orm.Float(-100.0),
    E_substrate=orm.Float(-80.0),
    E_molecule=orm.Float(-15.0)
)
print(f"Adsorption energy: {E_ads.value:.3f} eV")  # Output: -5.000 eV
```

### `compute_adsorption_energies_scatter()`

Complete scatter-gather workflow: separation + VASP relaxations + energy calculation.

**Purpose:** Automates the entire workflow for multiple adsorption sites in parallel.

**Workflow:**
1. Separates all substrate+adsorbate structures in parallel
2. Relaxes all systems (3N VASP jobs) in parallel
3. Calculates adsorption energies for each system

**Parameters:**
- `structures` (dict): Dictionary of complete StructureData nodes
- `adsorbate_formulas` (dict): Dictionary mapping structure keys to adsorbate formulas
  - Example: `{'site1': 'OOH', 'site2': 'OH'}`
- `code` (Code): AiiDA code for VASP
- `potential_family` (str): Pseudopotential family name
- `potential_mapping` (dict): Element to potential mapping
  - Example: `{'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'}`
- `relax_before_adsorption` (bool): If True, relax complete structures before SCF (default: False)
- `relax_parameters` (dict, optional): VASP INCAR parameters for relaxation (old-style API)
- `scf_parameters` (dict, optional): VASP INCAR parameters for SCF (old-style API)
- `relax_builder_inputs` (dict, optional): Complete VASP builder dict for relaxation (new-style API)
- `scf_builder_inputs` (dict, optional): Complete VASP builder dict for SCF (new-style API)
- `structure_specific_relax_builder_inputs` (dict, optional): Per-structure overrides for relaxation (NEW)
  - Format: `{idx: {'parameters': {'incar': {...}}, 'kpoints_spacing': 0.2, ...}}`
- `structure_specific_scf_builder_inputs` (dict, optional): Per-structure overrides for SCF (NEW)
  - Format: `{idx: {'parameters': {'incar': {...}}, 'options': {...}, ...}}`
- `options` (dict): Scheduler options for VASP calculations
- `kpoints_spacing` (float, optional): K-points spacing in Angstrom^-1
- `clean_workdir` (bool, optional): Whether to clean remote working directories (default: True)
- `fix_atoms` (bool): Enable atom fixing during relaxation (default: False)
- `fix_type` (str, optional): Type of fixing ('bottom', 'top', 'center')
- `fix_thickness` (float, optional): Thickness in Angstroms for fixing region (default: 0.0)
- `fix_elements` (list[str], optional): Element symbols to fix (default: None = all elements)

**Returns:**
Dictionary with namespaces:
- `separated_structures`: Dict of separated systems for each input
- `substrate_energies`: Energies of bare substrates (eV)
- `molecule_energies`: Energies of isolated molecules (eV)
- `complete_energies`: Energies of complete systems (eV)
- `adsorption_energies`: Final E_ads values (eV)

**Example:**
```python
from teros.core import compute_adsorption_energies_scatter
from aiida import orm

# Define structures and adsorbates
structures = {
    'site1': structure1,  # LaMnO3 + OOH
    'site2': structure2,  # LaMnO3 + OH
}

adsorbate_formulas = {
    'site1': 'OOH',
    'site2': 'OH',
}

# VASP parameters
parameters = {
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'NSW': 100,
    'ISIF': 2,  # Relax ions only, keep cell fixed
    'LWAVE': False,
    'LCHARG': False,
}

# Scheduler options
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 16,
    },
    'max_wallclock_seconds': 3600 * 12,
    'queue_name': 'regular',
}

# Run workflow
results = compute_adsorption_energies_scatter(
    structures=structures,
    adsorbate_formulas=adsorbate_formulas,
    code=vasp_code,
    potential_family='PBE',
    potential_mapping={'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    parameters=parameters,
    options=options,
    kpoints_spacing=0.25,
)

# Access results
for key, E_ads in results['adsorption_energies'].items():
    print(f"{key}: {E_ads.value:.3f} eV")
```

## Adsorbate Identification Method

The module uses **advanced connectivity analysis** to identify adsorbates robustly:

### Algorithm

1. **Identify candidate atoms** - Find all atoms matching elements in adsorbate formula
2. **Build bonding graph** using pymatgen's StructureGraph with CrystalNN for entire structure
3. **Create adsorbate subgraph** - Extract edges **only between candidate atoms**, ignoring bonds to substrate
4. **Find connected clusters** in subgraph - groups of candidate atoms bonded to each other
5. **Match cluster composition** to adsorbate formula
6. **Extract adsorbate atoms** while preserving substrate

### Advantages

- **Robust to different geometries** - no distance thresholds needed
- **No geometric assumptions** - works for any adsorbate orientation
- **Handles complex adsorbates** - multi-atom molecules like OOH, COOH
- **Works with periodic boundaries** - correct handling of surface slabs
- **✅ Works with bonded adsorbates** - Handles both weakly and strongly bonded species

### How It Handles Bonded Adsorbates

The key innovation is the **subgraph approach**:

1. **Traditional approach fails**: When OH bonds to Ag, the entire system (Ag + OH) becomes one connected cluster
2. **Subgraph solution**: Build a graph with edges **only between O and H atoms** (adsorbate elements)
3. **Result**: OH forms a separate connected component, even when bonded to surface

**Example:**
```
Structure: Ag₄ + OH (with O bonded to Ag)

Traditional connectivity:
  ├─ Cluster 1: Ag₄-O-H (all connected) ❌

Subgraph approach (edges only between O-H):
  ├─ Cluster 1: Ag₄ (no edges in subgraph)
  └─ Cluster 2: O-H (bonded together) ✅
```

### Requirements

- Adsorbate atoms must form a **single bonded cluster** (all adsorbate atoms connected to each other)
- Adsorbate formula must be **unique in structure** (no duplicate adsorbates)
- ✅ **Bonds to substrate are fine** - Algorithm ignores substrate-adsorbate bonds
- ✅ **Works for any bonding configuration** - From weak physisorption to strong chemisorption

### Supported Adsorbates

Examples of adsorbates that work well:
- **Simple radicals**: OH, O, H, OOH
- **Molecules**: H2O, CO, CO2, N2
- **Complex radicals**: COOH, CHO, CH3

## Cell Handling

All three structures (substrate, molecule, complete) use **the same simulation cell**:

### Why Same Cell?

1. **Consistency**: Ensures DFT calculations use identical basis sets
2. **Avoids BSSE**: Minimizes basis set superposition error
3. **Proper comparison**: Energy differences are meaningful
4. **No re-centering**: Molecule stays at original position relative to substrate

### Implications

- Molecule energy includes some vacuum penalty (acceptable for E_ads calculation)
- All three VASP calculations use identical k-point grids
- Results are directly comparable

## Best Practices

### 1. Structure Preparation

**Substrate:**
- Use sufficient vacuum (>10 Angstrom) to prevent periodic interactions
- Pre-relax substrate to minimize reconstruction
- Ensure slab thickness is converged

**Adsorbate:**
- Place adsorbate at desired adsorption site
- Ensure reasonable geometry (bond lengths ~1 Angstrom)
- Adsorbate should be distinguishable from substrate by composition

**Validation:**
- Check that adsorbate is identified correctly using connectivity
- Verify no spurious bonds between adsorbate and substrate atoms

### 2. VASP Parameters

**Convergence:**
- Use **consistent INCAR settings** for all three calculations
- Converge k-points with respect to adsorption energy
- Converge ENCUT (typically 1.3x ENMAX)
- Test vacuum thickness convergence

**Relaxation:**
- `ISIF=2`: Relax ions only, keep cell fixed (required)
- `IBRION=2`: Conjugate gradient (recommended)
- `NSW`: Sufficient steps for convergence (50-100)
- `EDIFF`: 1e-6 or tighter

**Electronic structure:**
- `ISMEAR=0` (Gaussian) for molecules and insulators
- `ISMEAR=1` (Methfessel-Paxton) for metals
- `SIGMA=0.05`: Small smearing
- Spin polarization if needed (`ISPIN=2`, `MAGMOM`)

**Performance:**
- `LWAVE=False`: Don't write wavefunctions (saves space)
- `LCHARG=False`: Don't write charge density (saves space)
- `NCORE`: Set to number of cores per node for speed

### 3. Validation

**Energy checks:**
- Verify all three calculations converged (check OUTCAR for NSW reached)
- Check that E_ads is reasonable (compare with literature)
- Test sensitivity to vacuum thickness
- Verify k-point convergence

**Structure checks:**
- Ensure substrate doesn't reconstruct significantly
- Check molecule geometry after relaxation
- Verify adsorbate position is reasonable

**Comparison:**
- Compare with experimental values if available
- Test against literature DFT values
- Cross-check with different functionals (PBE, PBE+U, hybrid)

### 4. Common Issues

**Problem: "Could not find adsorbate among candidate atom clusters"**
- Solution: Check adsorbate formula matches structure composition
- Solution: Verify adsorbate atoms are bonded **to each other** (not just to substrate)
- Solution: Ensure adsorbate elements are present in structure
- Note: Bonds to substrate are OK! Algorithm ignores them.

**Problem: "Found multiple clusters matching"**
- Solution: Ensure only one adsorbate per structure
- Solution: Use different adsorbate formulas for each site
- Solution: If multiple identical adsorbates needed, submit as separate structures

**Problem: E_ads is very positive (unfavorable)**
- Check substrate relaxation - may need more NSW steps
- Verify k-points are converged
- Check spin polarization settings
- Consider using DFT+U for correlated systems
- Verify adsorbate is actually in contact with surface

**Problem: "Not enough atoms to form adsorbate"**
- Check structure has required elements (e.g., O and H for OH)
- Verify formula string is correct (case-sensitive)
- Check structure file loaded correctly

## Integration with Experimental Tools

The module can use structures generated by experimental tools:

```python
# Example: Use experimental surface builder to create substrate+adsorbate
from teros.experimental.adsorption_energy.lamno3.add_ooh_to_surface import main as add_ooh

# Generate structure with OOH
output_cif = "lamno3_001_ooh.cif"
add_ooh("lamno3_001_slab.cif", output_cif)

# Load structure
from teros.core import get_structure_from_file
complete_structure = get_structure_from_file(filepath=output_cif)

# Calculate adsorption energy
results = compute_adsorption_energies_scatter(
    structures={'ooh_site': complete_structure},
    adsorbate_formulas={'ooh_site': 'OOH'},
    code=vasp_code,
    potential_family='PBE',
    potential_mapping={'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    parameters=parameters,
    options=options,
)
```

## Usage in WorkGraph

### Complete Example

```python
from aiida import orm, load_profile
from aiida_workgraph import WorkGraph
from teros.core import compute_adsorption_energies_scatter, get_structure_from_file

load_profile()

# Load structures
structures = {
    'site1': get_structure_from_file(filepath='lamno3_ooh.cif'),
    'site2': get_structure_from_file(filepath='lamno3_oh.cif'),
}

# Define adsorbates
adsorbate_formulas = {
    'site1': 'OOH',
    'site2': 'OH',
}

# Load VASP code
code = orm.load_code('vasp@localhost')

# Create WorkGraph
wg = WorkGraph('adsorption_energies')

# Add task
ads_task = wg.add_task(
    compute_adsorption_energies_scatter,
    name='compute_adsorption_energies',
    structures=structures,
    adsorbate_formulas=adsorbate_formulas,
    code=code,
    potential_family='PBE',
    potential_mapping={'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    parameters={
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'NSW': 100,
        'ISIF': 2,
    },
    options={
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 16},
        'max_wallclock_seconds': 3600 * 12,
    },
    kpoints_spacing=0.25,
)

# Submit
wg.submit(wait=False)
print(f"WorkGraph submitted: {wg.pk}")
```

### Access Results

```python
from aiida import load_node

# Load completed WorkGraph
wg = load_node(PK)

# Get adsorption energies
ads_energies = wg.outputs.compute_adsorption_energies.adsorption_energies

for key, E_ads in ads_energies.items():
    print(f"{key}: {E_ads.value:.3f} eV")

# Get individual energies
substrate_energies = wg.outputs.compute_adsorption_energies.substrate_energies
molecule_energies = wg.outputs.compute_adsorption_energies.molecule_energies
complete_energies = wg.outputs.compute_adsorption_energies.complete_energies

# Get separated structures
separated = wg.outputs.compute_adsorption_energies.separated_structures
substrate_struct = separated['site1']['substrate']
molecule_struct = separated['site1']['molecule']
```

## Output Data Structure

### Per-site Results

For each input structure, the workflow returns:

**Separated structures:**
- `substrate`: Bare surface (StructureData)
- `molecule`: Isolated adsorbate (StructureData)
- `complete`: Original structure (StructureData)

**Energies:**
- `substrate_energies[key]`: E_substrate in eV (Float)
- `molecule_energies[key]`: E_molecule in eV (Float)
- `complete_energies[key]`: E_complete in eV (Float)
- `adsorption_energies[key]`: E_ads in eV (Float)

## Example

See `examples/adsorption_energy/test_oh_ag111/` for a complete working example with OH on Ag(111).

The example demonstrates:
- Structure setup (Ag slab with OH adsorbate)
- VASP parameter configuration
- WorkGraph setup and submission
- Result extraction and analysis

## Comparison with Other Modules

### Similar to Cleavage Energy Module

Both modules follow the same design pattern:
- Scatter-gather workflow for parallel processing
- Separate calcfunctions for component calculations
- Integration with VASP via aiida-vasp
- Provenance tracking via AiiDA

### Differences

| Feature | Adsorption Energy | Cleavage Energy |
|---------|------------------|-----------------|
| Input | Single structure | Complementary slab pairs |
| Separation | Connectivity analysis | Index-based splitting |
| Systems | 3 per input | 2 per pair |
| Formula | E_complete - E_sub - E_mol | (E_i + E_j - n*E_bulk)/(2A) |
| Units | eV | eV/A^2 and J/m^2 |

## See Also

- **Example**: `examples/adsorption_energy/test_oh_ag111/`
- **Related modules**: `teros.core.thermodynamics`, `teros.core.cleavage`
- **Experimental tools**: `teros.experimental.adsorption_energy/`
- **Tests**: `teros/core/test_adsorption_energy.py`

## References

Standard adsorption energy formula:
```
E_ads = E_complete - E_substrate - E_molecule
```

This formula is widely used in computational catalysis and surface science literature.
