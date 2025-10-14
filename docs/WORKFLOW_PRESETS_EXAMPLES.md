# PS-TEROS Workflow Presets - Code Examples

This document provides complete, runnable examples for each workflow preset.

---

## Table of Contents

1. [Common Setup](#common-setup)
2. [Example 1: Surface Thermodynamics (Default)](#example-1-surface-thermodynamics)
3. [Example 2: Quick Screening (Unrelaxed)](#example-2-quick-screening)
4. [Example 3: Cleavage Energy Only](#example-3-cleavage-energy-only)
5. [Example 4: Relaxation Energy Only](#example-4-relaxation-energy-only)
6. [Example 5: Bulk Only](#example-5-bulk-only)
7. [Example 6: Formation Enthalpy Only](#example-6-formation-enthalpy-only)
8. [Example 7: Electronic Structure (Bulk)](#example-7-electronic-structure)
9. [Example 8: Electronic Structure (Slabs Only)](#example-8-electronic-structure-slabs-only)
10. [Example 9: Electronic Structure (Bulk and Slabs)](#example-9-electronic-structure-bulk-and-slabs)
11. [Example 10: AIMD Simulation](#example-10-aimd-simulation)
12. [Example 11: Comprehensive Analysis](#example-11-comprehensive-analysis)
13. [Example 12: Overriding Preset Defaults](#example-12-overriding-defaults)

---

## Common Setup

All examples use this common setup:

```python
from teros.core.workgraph import build_core_workgraph
from aiida.engine import submit

# Common parameters
STRUCTURES_DIR = '/path/to/structures'
CODE_LABEL = 'VASP-VTST-6.4.3@bohr'
POTENTIAL_FAMILY = 'PBE'

# Common VASP parameters
BULK_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'ALGO': 'Normal',
    'IBRION': 2,
    'NSW': 100,
    'ISIF': 3,
}

BULK_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 48,
    },
    'max_wallclock_seconds': 3600 * 10,
    'queue_name': 'regular',
}
```

---

## Example 1: Surface Thermodynamics

Complete surface thermodynamics workflow (default preset).

```python
"""
Example: Complete surface thermodynamics workflow for Ag3PO4
Calculates surface energies, cleavage energies, and relaxation energies
"""

# Build workflow with preset
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # Activates full workflow
    
    # Structure files
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    
    # Code configuration
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    kpoints_spacing=0.4,
    clean_workdir=False,
    
    # Bulk calculation
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # Reference calculations
    metal_potential_mapping={'Ag': 'Ag'},
    metal_parameters=BULK_PARAMS.copy(),
    metal_options=BULK_OPTIONS,
    
    oxygen_potential_mapping={'O': 'O'},
    oxygen_parameters=BULK_PARAMS.copy(),
    oxygen_options=BULK_OPTIONS,
    
    nonmetal_potential_mapping={'P': 'P'},
    nonmetal_parameters=BULK_PARAMS.copy(),
    nonmetal_options=BULK_OPTIONS,
    
    # Slab configuration
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,
    
    # Thermodynamics
    thermodynamics_sampling=100,
    
    name='Ag3PO4_100_Thermodynamics'
)

# Submit workflow
result = submit(wg)
print(f"Submitted WorkGraph: {result.pk}")
```

**What this does:**
1. Relaxes Ag3PO4 bulk structure
2. Relaxes reference structures (Ag, P, O2)
3. Calculates formation enthalpy
4. Generates (100) slabs
5. Relaxes all slab terminations
6. Calculates surface energies with chemical potential sampling

**Note:** Cleavage and relaxation energies are NOT computed by default. To include them, add:
```python
compute_cleavage=True,
compute_relaxation_energy=True,
```

---

## Example 2: Quick Screening

Quick surface energy screening with unrelaxed slabs.

```python
"""
Example: Quick screening of multiple Miller indices
Fast assessment without slab relaxation
"""

wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics_unrelaxed',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # References
    metal_potential_mapping={'Ag': 'Ag'},
    metal_parameters=BULK_PARAMS.copy(),
    metal_options=BULK_OPTIONS,
    
    oxygen_potential_mapping={'O': 'O'},
    oxygen_parameters=BULK_PARAMS.copy(),
    oxygen_options=BULK_OPTIONS,
    
    nonmetal_potential_mapping={'P': 'P'},
    nonmetal_parameters=BULK_PARAMS.copy(),
    nonmetal_options=BULK_OPTIONS,
    
    # Slab - testing (110) surface
    miller_indices=[1, 1, 0],
    min_slab_thickness=15.0,
    min_vacuum_thickness=15.0,
    
    thermodynamics_sampling=50,  # Coarser sampling for screening
    
    name='Ag3PO4_110_Screening'
)

result = submit(wg)
print(f"Submitted screening workflow: {result.pk}")
```

**Use case:** Quickly test multiple Miller indices to identify most stable surfaces before running expensive relaxations.

---

## Example 3: Cleavage Energy Only

Calculate cleavage energies for complementary slab pairs.

```python
"""
Example: Cleavage energy analysis for Ag2O (100)
No thermodynamics, just cleavage energies
"""

wg = build_core_workgraph(
    workflow_preset='cleavage_only',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag2o.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # Slab configuration
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,
    
    name='Ag2O_100_Cleavage'
)

result = submit(wg)
print(f"Submitted cleavage workflow: {result.pk}")
```

**What this does:**
1. Relaxes bulk structure
2. Generates slabs
3. Relaxes slabs
4. Calculates cleavage energies (no formation enthalpy or surface energies)

---

## Example 4: Relaxation Energy Only

Analyze surface reconstruction by calculating relaxation energies.

```python
"""
Example: Relaxation energy analysis
Measures how much slabs relax from bulk termination
"""

wg = build_core_workgraph(
    workflow_preset='relaxation_energy_only',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # Slab configuration
    miller_indices=[1, 1, 1],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,
    
    name='Ag3PO4_111_Relaxation'
)

result = submit(wg)
print(f"Submitted relaxation energy workflow: {result.pk}")
```

**Use case:** Understanding surface reconstruction and comparing relaxation magnitudes across terminations.

---

## Example 5: Bulk Only

Simple bulk structure optimization.

```python
"""
Example: Bulk optimization only
No surface calculations
"""

wg = build_core_workgraph(
    workflow_preset='bulk_only',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    name='Ag3PO4_Bulk'
)

result = submit(wg)
print(f"Submitted bulk-only workflow: {result.pk}")
```

**Use case:** Testing calculation parameters, initial structure validation, or when you only need bulk properties.

---

## Example 6: Formation Enthalpy Only

Calculate formation enthalpy without surfaces.

```python
"""
Example: Formation enthalpy calculation
Useful for phase diagrams and stability analysis
"""

wg = build_core_workgraph(
    workflow_preset='formation_enthalpy_only',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    # Bulk
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # References
    metal_potential_mapping={'Ag': 'Ag'},
    metal_parameters=BULK_PARAMS.copy(),
    metal_options=BULK_OPTIONS,
    
    oxygen_potential_mapping={'O': 'O'},
    oxygen_parameters=BULK_PARAMS.copy(),
    oxygen_options=BULK_OPTIONS,
    
    nonmetal_potential_mapping={'P': 'P'},
    nonmetal_parameters=BULK_PARAMS.copy(),
    nonmetal_options=BULK_OPTIONS,
    
    name='Ag3PO4_FormationEnthalpy'
)

result = submit(wg)
print(f"Submitted formation enthalpy workflow: {result.pk}")
```

---

## Example 7: Electronic Structure

Calculate electronic properties (DOS and band structure) for bulk.

```python
"""
Example: Electronic structure calculation
DOS and band structure for bulk Ag3PO4
"""

from teros.core.builders import get_electronic_properties_defaults

# Get default electronic properties parameters
elec_defaults = get_electronic_properties_defaults()

wg = build_core_workgraph(
    workflow_preset='electronic_structure_bulk_only',
    
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # Electronic properties parameters
    bands_parameters=elec_defaults['bands_parameters'],
    bands_options=BULK_OPTIONS,
    band_settings=elec_defaults['band_settings'],
    
    name='Ag3PO4_Electronic_Structure'
)

result = submit(wg)
print(f"Submitted electronic structure workflow: {result.pk}")
```

**What this calculates:**
- Band structure along high-symmetry k-points
- Total and projected density of states
- Band gap (if semiconductor/insulator)

---

## Example 8: Electronic Structure (Slabs Only)

Calculate electronic properties for slabs only.

```python
"""
Example: Electronic structure for Ag3PO4 (100) slabs
DOS and band structure for surface terminations
"""

from teros.core.builders import get_slab_electronic_properties_defaults

# Get default slab electronic properties parameters (denser k-points for 2D)
slab_elec_defaults = get_slab_electronic_properties_defaults()

wg = build_core_workgraph(
    workflow_preset='electronic_structure_slabs_only',

    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',

    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,

    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,

    # Slab configuration
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,

    # Slab electronic properties parameters
    slab_bands_parameters=slab_elec_defaults['bands_parameters'],
    slab_bands_options=BULK_OPTIONS,
    slab_band_settings=slab_elec_defaults['band_settings'],

    name='Ag3PO4_100_Slab_Electronic'
)

result = submit(wg)
print(f"Submitted slab electronic structure workflow: {result.pk}")
```

**What this calculates:**
- Band structure for each slab termination
- Total and projected density of states for each slab
- Band gap for each termination (if semiconductor/insulator)

**Use case:** Comparing electronic properties of different surface terminations.

---

## Example 9: Electronic Structure (Bulk and Slabs)

Calculate electronic properties for both bulk and slabs.

```python
"""
Example: Complete electronic structure analysis
DOS and bands for both bulk and slabs
"""

from teros.core.builders import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults
)

# Electronic properties for bulk and slabs
elec_defaults = get_electronic_properties_defaults()
slab_elec_defaults = get_slab_electronic_properties_defaults()

wg = build_core_workgraph(
    workflow_preset='electronic_structure_bulk_and_slabs',

    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',

    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,

    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,

    # Slab configuration
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,

    # Bulk electronic properties
    bands_parameters=elec_defaults['bands_parameters'],
    bands_options=BULK_OPTIONS,
    band_settings=elec_defaults['band_settings'],

    # Slab electronic properties
    slab_bands_parameters=slab_elec_defaults['bands_parameters'],
    slab_bands_options=BULK_OPTIONS,
    slab_band_settings=slab_elec_defaults['band_settings'],

    name='Ag3PO4_100_Complete_Electronic'
)

result = submit(wg)
print(f"Submitted complete electronic structure workflow: {result.pk}")
```

**What this calculates:**
- All bulk electronic properties (from Example 7)
- All slab electronic properties (from Example 8)

**Use case:** Comprehensive electronic structure comparison between bulk and surfaces.

---

## Example 10: AIMD Simulation

Run Ab Initio Molecular Dynamics on slabs.

```python
"""
Example: AIMD simulation on Ag3PO4 (100) surface
Multi-stage temperature equilibration and production
"""

from teros.core.builders import get_aimd_defaults

# AIMD sequence: equilibration → production
aimd_sequence = [
    {'temperature': 300, 'steps': 1000},  # Equilibration
    {'temperature': 300, 'steps': 5000},  # Production run
]

# Get default AIMD parameters
aimd_params = get_aimd_defaults()

wg = build_core_workgraph(
    workflow_preset='aimd_only',

    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',

    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,

    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,

    # Slab configuration
    miller_indices=[1, 0, 0],
    min_slab_thickness=15.0,
    min_vacuum_thickness=15.0,

    # AIMD configuration
    aimd_sequence=aimd_sequence,
    aimd_parameters=aimd_params,
    aimd_options=BULK_OPTIONS,
    aimd_kpoints_spacing=0.5,  # Coarser for AIMD

    name='Ag3PO4_100_AIMD'
)

result = submit(wg)
print(f"Submitted AIMD workflow: {result.pk}")
```

**Note:** By default, `aimd_only` runs AIMD on **unrelaxed slabs** (generated from bulk). If you need relaxed slabs before AIMD, add:
```python
relax_slabs=True,
slab_parameters=BULK_PARAMS.copy(),
slab_options=BULK_OPTIONS,
```

---

## Example 11: Comprehensive Analysis

Complete analysis with all features enabled.

```python
"""
Example: Comprehensive analysis
Everything: thermodynamics + electronic properties + AIMD
"""

from teros.core.builders import (
    get_electronic_properties_defaults,
    get_aimd_defaults,
    get_slab_electronic_properties_defaults
)

# Electronic properties
elec_defaults = get_electronic_properties_defaults()
slab_elec_defaults = get_slab_electronic_properties_defaults()

# AIMD
aimd_params = get_aimd_defaults()
aimd_sequence = [
    {'temperature': 300, 'steps': 500},
]

wg = build_core_workgraph(
    workflow_preset='comprehensive',
    
    # Structures
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    
    # Code
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    # Bulk
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    # References
    metal_potential_mapping={'Ag': 'Ag'},
    metal_parameters=BULK_PARAMS.copy(),
    metal_options=BULK_OPTIONS,
    
    oxygen_potential_mapping={'O': 'O'},
    oxygen_parameters=BULK_PARAMS.copy(),
    oxygen_options=BULK_OPTIONS,
    
    nonmetal_potential_mapping={'P': 'P'},
    nonmetal_parameters=BULK_PARAMS.copy(),
    nonmetal_options=BULK_OPTIONS,
    
    # Slabs
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,
    
    # Thermodynamics
    thermodynamics_sampling=100,
    
    # Bulk electronic properties
    bands_parameters=elec_defaults['bands_parameters'],
    bands_options=BULK_OPTIONS,
    band_settings=elec_defaults['band_settings'],
    
    # Slab electronic properties
    slab_bands_parameters=slab_elec_defaults['bands_parameters'],
    slab_bands_options=BULK_OPTIONS,
    slab_band_settings=slab_elec_defaults['band_settings'],
    
    # AIMD
    aimd_sequence=aimd_sequence,
    aimd_parameters=aimd_params,
    aimd_options=BULK_OPTIONS,
    
    name='Ag3PO4_100_Comprehensive'
)

result = submit(wg)
print(f"Submitted comprehensive workflow: {result.pk}")
```

**Warning:** This is expensive! Only use for final production runs.

---

## Example 12: Overriding Defaults

Override specific features from a preset.

```python
"""
Example: Custom configuration by overriding preset
Start with surface_thermodynamics but disable cleavage
"""

wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    
    # Override: disable cleavage energies
    compute_cleavage=False,
    
    # Override: disable relaxation energies
    compute_relaxation_energy=False,
    
    # Rest of parameters as usual
    structures_dir=STRUCTURES_DIR,
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    nonmetal_name='P.cif',
    
    code_label=CODE_LABEL,
    potential_family=POTENTIAL_FAMILY,
    
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    bulk_parameters=BULK_PARAMS,
    bulk_options=BULK_OPTIONS,
    
    metal_potential_mapping={'Ag': 'Ag'},
    metal_parameters=BULK_PARAMS.copy(),
    metal_options=BULK_OPTIONS,
    
    oxygen_potential_mapping={'O': 'O'},
    oxygen_parameters=BULK_PARAMS.copy(),
    oxygen_options=BULK_OPTIONS,
    
    nonmetal_potential_mapping={'P': 'P'},
    nonmetal_parameters=BULK_PARAMS.copy(),
    nonmetal_options=BULK_OPTIONS,
    
    miller_indices=[1, 0, 0],
    slab_parameters=BULK_PARAMS.copy(),
    slab_options=BULK_OPTIONS,
    
    name='Ag3PO4_100_Custom'
)

result = submit(wg)
print(f"Submitted custom workflow: {result.pk}")
```

**Result:** Surface energies calculated, but cleavage and relaxation energies skipped.

---

## Summary

These examples cover all workflow presets and common use cases. Key points:

1. **Use presets** for simplicity and consistency
2. **Override** when you need custom behavior
3. **Start simple** (bulk_only) and gradually add complexity
4. **Check requirements** - each preset has different parameter needs
5. **Use helper functions** - `get_electronic_properties_defaults()`, `get_aimd_defaults()`, etc.

For more details, see the [Workflow Presets User Guide](WORKFLOW_PRESETS_GUIDE.md).
