# PS-TEROS Step-by-Step Examples

This directory contains a series of progressive examples that test each component of the PS-TEROS workflow individually. Each step builds upon the previous ones, allowing you to understand and test the workflow chain systematically.

## Purpose

These examples are designed to:
1. **Test individual components** of the PS-TEROS workflow
2. **Debug issues** step-by-step rather than running everything at once
3. **Learn** how each calculation fits into the complete workflow
4. **Validate** your AiiDA/VASP setup progressively

## Prerequisites

1. **AiiDA setup:**
   ```bash
   verdi profile set-default psteros
   verdi status
   verdi daemon start  # if not running
   ```

2. **Structure files:**
   - Located in `../complete/structures/`
   - Need: `ag2o.cif`, `Ag.cif`, `O2.cif`

3. **VASP code:**
   - Code label: `VASP-VTST-6.4.3@bohr` (adjust in scripts if different)
   - Potential family: `PBE`

## Step-by-Step Progression

### Step 1: Bulk Relaxation Only
**File:** `step_01_bulk_only.py`  
**Preset:** `bulk_only`

**What it tests:**
- Basic AiiDA/VASP integration
- Bulk structure relaxation only

**Runtime:** ~30 minutes (depends on system)

**Expected outputs:**
- `bulk_energy`: Total energy of relaxed Ag2O
- `bulk_structure`: Relaxed structure

**Why start here:**
If this fails, check your AiiDA code/computer setup before proceeding.

```bash
python step_01_bulk_only.py
verdi process status <PK>
```

---

### Step 2: Formation Enthalpy
**File:** `step_02_formation_enthalpy.py`  
**Preset:** `formation_enthalpy_only`

**What it tests:**
- Multiple parallel relaxations (Ag2O, Ag, O2)
- Formation enthalpy calculation

**Runtime:** ~1 hour (parallel execution)

**Expected outputs:**
- `bulk_energy`, `metal_energy`, `oxygen_energy`
- `formation_enthalpy`: ΔHf = E(Ag2O) - 2*E(Ag) - 0.5*E(O2)

**Why this step:**
Tests that parallel calculations and formation enthalpy logic work correctly.

```bash
python step_02_formation_enthalpy.py
verdi process status <PK>
```

---

### Step 3: Slab Relaxation
**File:** `step_03_slabs_relaxation.py`  
**Preset:** `relaxation_energy_only`

**What it tests:**
- Slab generation from bulk
- Slab relaxation
- Relaxation energy calculation (E_relaxed - E_unrelaxed)

**Runtime:** ~2-4 hours (depends on number of terminations)

**Expected outputs:**
- `slab_structures`: Generated (111) terminations
- `unrelaxed_slab_energies`: SCF energies
- `slab_energies`: Relaxed energies
- `relaxation_energies`: ΔE for each slab

**Why this step:**
Tests the slab workflow without the complexity of thermodynamics.

```bash
python step_03_slabs_relaxation.py
verdi process status <PK>
```

---

### Step 4: Cleavage Energy
**File:** `step_04_cleavage_energy.py`  
**Preset:** `cleavage_only`

**What it tests:**
- Bulk and slab relaxations
- Cleavage energy calculation for complementary pairs

**Runtime:** ~2-4 hours

**Expected outputs:**
- `cleavage_energies`: E_cleave = E_slab1 + E_slab2 - E_bulk

**Why this step:**
Tests cleavage energy logic separately from thermodynamics.

```bash
python step_04_cleavage_energy.py
verdi process status <PK>
```

---

### Step 5: Surface Thermodynamics (Complete)
**File:** `step_05_surface_thermodynamics.py`  
**Preset:** `surface_thermodynamics` (default)

**What it tests:**
- **Core workflow:**
  - Bulk and reference relaxations
  - Formation enthalpy
  - Slab generation and relaxation
  - Surface energies with chemical potential sampling
- **Optional (not included by default):**
  - Cleavage energies (add: `compute_cleavage=True`)
  - Relaxation energies (add: `compute_relaxation_energy=True`)

**Runtime:** ~3-4 hours

**Expected outputs:**
- `bulk_energy`, `metal_energy`, `oxygen_energy`
- `formation_enthalpy`
- `slab_structures`, `slab_energies`
- `surface_energies`: γ(μ_O) for each termination

**Why this step:**
This is the main PS-TEROS workflow. If Steps 1-4 work, this should work too.

```bash
python step_05_surface_thermodynamics.py
verdi process status <PK>
```

---

### Step 6: Electronic Properties (Bulk)
**File:** `step_06_electronic_properties.py`  
**Preset:** `electronic_structure_bulk_only`

**What it tests:**
- Bulk relaxation
- Electronic structure (DOS and band structure) for bulk

**Runtime:** ~1-2 hours

**Expected outputs:**
- `bulk_bands`: Band structure
- `bulk_dos`: Density of states
- `bulk_primitive_structure`, `bulk_seekpath_parameters`

**Why this step:**
Tests electronic properties workflow independently.

```bash
python step_06_electronic_properties.py
verdi process status <PK>
```

---

### Step 7: AIMD Simulation
**File:** `step_07_aimd_simulation.py`  
**Preset:** `aimd_only`

**What it tests:**
- Bulk relaxation
- Slab generation (no relaxation)
- AIMD (Ab Initio Molecular Dynamics) simulation

**Runtime:** ~8-24 hours (VERY EXPENSIVE!)

**Expected outputs:**
- `aimd_results`: Trajectories, energies for each stage

**Why this step:**
Tests AIMD workflow. Only run if you need molecular dynamics.

**⚠️ WARNING:** AIMD is computationally expensive. Make sure you have sufficient resources.
**⚠️ NOTE:** Slabs are NOT relaxed before AIMD in the `aimd_only` preset.

```bash
python step_07_aimd_simulation.py
verdi process status <PK>
```

---

### Step 8: Electronic Properties (Slabs Only)
**File:** `step_08_electronic_structure_slabs.py`  
**Preset:** `electronic_structure_slabs_only`

**What it tests:**
- Bulk relaxation
- Slab generation and relaxation
- Electronic structure (DOS and bands) for slabs

**Runtime:** ~3-5 hours

**Expected outputs:**
- `slab_bands`: Band structure for each slab
- `slab_dos`: Density of states for each slab
- `slab_primitive_structures`, `slab_seekpath_parameters`

**Why this step:**
Tests slab electronic properties workflow.

```bash
python step_08_electronic_structure_slabs.py
verdi process status <PK>
```

---

### Step 9: Electronic Properties (Bulk and Slabs)
**File:** `step_09_electronic_structure_bulk_and_slabs.py`  
**Preset:** `electronic_structure_bulk_and_slabs`

**What it tests:**
- Complete electronic structure analysis for both bulk and slabs
- Comprehensive DOS and band structure calculations

**Runtime:** ~4-6 hours

**Expected outputs:**
- All bulk electronic properties (from Step 6)
- All slab electronic properties (from Step 8)

**Why this step:**
Tests complete electronic structure workflow for comparison studies.

```bash
python step_09_electronic_structure_bulk_and_slabs.py
verdi process status <PK>
```

---

### Step 10: Custom Workflow
**File:** `step_10_custom_workflow.py`  
**Preset:** None (custom flags)

**What it tests:**
- Building workflows with individual flag configuration
- Custom combinations of workflow components

**Example:** Surface thermodynamics + cleavage, but no relaxation energies

**Why this step:**
Demonstrates maximum flexibility for custom workflows.

```bash
python step_10_custom_workflow.py
verdi process status <PK>
```

---

### Step 11: Preset with Overrides
**File:** `step_11_preset_with_overrides.py`  
**Preset:** `surface_thermodynamics` + overrides

**What it tests:**
- Starting with a preset and adding optional components
- Combining preset convenience with custom control

**Example:** Start with `surface_thermodynamics`, add cleavage and relaxation energies

**Why this step:**
Shows how to customize presets without building from scratch.

```bash
python step_11_preset_with_overrides.py
verdi process status <PK>
```

---

### Step 12: Surface Hydroxylation and Vacancy Generation
**File:** `step_12_surface_hydroxylation.py`
**Module:** `teros.core.surface_hydroxylation`

**What it tests:**
- Generate surface variants with different OH coverages
- Generate oxygen-deficient surfaces (vacancies)
- Coverage-based deduplication algorithm (reduces thousands to ~10 structures)
- Batch VASP relaxation with controlled parallelization
- Result organization and analysis

**Runtime:** ~30 minutes - 2 hours (depends on coverage_bins and system size)

**Expected outputs:**
- `manifest`: Dict with metadata for all generated variants
- `structures`: Namespace dict with relaxed structures (descriptive keys: `idx_variantname`)
- `energies`: Namespace dict with energies for each structure

**Output naming example:**
- `0_oh_000_3_7572`: First structure, 3.76 OH/nm² coverage
- `1_oh_001_7_5145`: Second structure, 7.51 OH/nm² coverage

**Why this step:**
Tests the surface_hydroxylation module for studying:
- Hydroxylation effects on surface stability
- Oxygen vacancy formation energies
- Coverage-dependent properties
- Surface reactivity under different conditions

**Post-processing:**
```python
from aiida import orm
from teros.core.surface_hydroxylation import organize_hydroxylation_results

node = orm.load_node(<PK>)
results = organize_hydroxylation_results(node)
print('Statistics:', results['statistics'])
for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV (coverage={r['coverage']:.2f})")
```

**Key features:**
- Three modes: `'hydrogen'` (OH groups), `'vacancies'` (remove O), `'combine'` (both)
- Coverage-based deduplication to reduce combinatorial explosion
- Descriptive output names for easy identification
- Builder function pattern for simplified workflow creation

**Documentation:** `../../docs/surface_hydroxylation.md`

```bash
python step_12_surface_hydroxylation.py
verdi process status <PK>
```

---

## Workflow Diagram

```
Step 1: Bulk Only
   └─> Ag2O relaxation

Step 2: Formation Enthalpy
   ├─> Ag2O relaxation
   ├─> Ag relaxation
   ├─> O2 relaxation
   └─> ΔHf calculation

Step 3: Slab Relaxation
   ├─> Ag2O relaxation
   ├─> Generate slabs
   ├─> SCF on slabs
   ├─> Relax slabs
   └─> ΔE_relax

Step 4: Cleavage Energy
   ├─> Ag2O relaxation
   ├─> Generate + relax slabs
   └─> E_cleave

Step 5: Surface Thermodynamics
   ├─> Ag2O, Ag, O2 relaxations
   ├─> ΔHf calculation
   ├─> Generate + relax slabs
   └─> γ(μ_O)
   (Optional: E_cleave, ΔE_relax)

Step 6: Electronic Properties (Bulk)
   ├─> Ag2O relaxation
   └─> DOS + Bands

Step 7: AIMD
   ├─> Ag2O relaxation
   ├─> Generate slabs (no relax)
   └─> AIMD simulation

Step 8: Electronic Properties (Slabs)
   ├─> Ag2O relaxation
   ├─> Generate + relax slabs
   └─> DOS + Bands (slabs)

Step 9: Electronic Properties (Bulk + Slabs)
   ├─> Ag2O relaxation
   ├─> Generate + relax slabs
   └─> DOS + Bands (both)

Step 10: Custom Workflow
   └─> User-defined component combination

Step 11: Preset with Overrides
   └─> Preset + custom additions

Step 12: Surface Hydroxylation
   ├─> Relaxed slab input
   ├─> Generate hydroxylated/vacancy structures
   │   └─> Coverage-based deduplication
   ├─> Batch VASP relaxations
   └─> Organize results by coverage
```

---

## Monitoring Workflows

### Check status:
```bash
verdi process status <PK>
verdi process show <PK>
```

### List all running:
```bash
verdi process list
verdi process list -a -p1  # More detailed
```

### Watch daemon:
```bash
verdi daemon status
watch -n 5 'verdi process list | head -20'
```

---

## Troubleshooting

### If Step 1 fails:
- Check `verdi code show <code_label>`
- Check `verdi computer test <computer_name>`
- Verify structure file exists

### If Step 2 fails:
- Make sure Step 1 works first
- Check all structure files exist (Ag.cif, O2.cif)
- Verify potential mappings are correct

### If Step 3+ fails:
- Work backwards: which earlier step works?
- Check daemon logs: `verdi daemon logshow`
- Look at specific process report: `verdi process report <PK>`

### General debugging:
```bash
# Restart daemon
verdi daemon restart

# Check failed process
verdi process show <PK>
verdi process report <PK>

# View calculation outputs
verdi calcjob outputcat <PK>
```

---

## Tips

1. **Start from Step 1** - don't skip ahead
2. **Wait for completion** - don't submit Step 2 until Step 1 finishes
3. **Check outputs** - use `verdi process show <PK>` to verify outputs
4. **Save PKs** - write down the PK of each successful run
5. **Be patient** - VASP calculations take time

---

## Customization

To adapt these examples for your system:

1. **Change code label:**
   ```python
   code_label = 'YOUR-VASP-CODE@YOUR-COMPUTER'
   ```

2. **Adjust resources:**
   ```python
   'num_machines': 2,  # More machines
   'num_cores_per_machine': 48,  # Different core count
   ```

3. **Different material:**
   - Replace structure files
   - Update potential mappings
   - Adjust VASP parameters for your system

4. **Different surface:**
   ```python
   miller_indices=[1, 0, 0]  # (100) surface
   ```

---

## Next Steps

After completing these step-by-step examples:

1. Run the **complete example** in `../complete/complete_ag2o_example.py`
2. Try the **workflow presets** in `../workflow_presets/`
3. Apply PS-TEROS to **your own materials**

---

## Questions?

- See main documentation: `../../docs/`
- **Workflow system guide**: `README_WORKFLOWS.md` (NEW!)
- Workflow preset details: Check the presets in the code

For detailed workflow configuration options, see `README_WORKFLOWS.md` which explains:
- All available workflow presets
- How to use custom workflows
- How to override preset defaults
- Examples for each use case

---

**Good luck testing your PS-TEROS workflows step by step!** 🚀
