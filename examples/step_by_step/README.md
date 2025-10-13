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
- **EVERYTHING combined:**
  - Bulk and reference relaxations
  - Formation enthalpy
  - Slab generation and relaxation
  - Surface energies with chemical potential sampling
  - Cleavage energies
  - Relaxation energies

**Runtime:** ~4-6 hours

**Expected outputs:**
- All of the above
- `surface_energies`: γ(μ_O) for each termination

**Why this step:**
This is the main PS-TEROS workflow. If Steps 1-4 work, this should work too.

```bash
python step_05_surface_thermodynamics.py
verdi process status <PK>
```

---

### Step 6: Electronic Properties
**File:** `step_06_electronic_properties.py`  
**Preset:** `electronic_structure_bulk_only`

**What it tests:**
- Bulk relaxation
- Electronic structure (DOS and band structure)

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
- Bulk and slab relaxations
- AIMD (Ab Initio Molecular Dynamics) simulation

**Runtime:** ~8-24 hours (VERY EXPENSIVE!)

**Expected outputs:**
- `aimd_results`: Trajectories, energies for each stage

**Why this step:**
Tests AIMD workflow. Only run if you need molecular dynamics.

**⚠️ WARNING:** AIMD is computationally expensive. Make sure you have sufficient resources.

```bash
python step_07_aimd_simulation.py
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

Step 5: Surface Thermodynamics (COMPLETE)
   ├─> Ag2O, Ag, O2 relaxations
   ├─> ΔHf calculation
   ├─> Generate + relax slabs
   ├─> ΔE_relax
   ├─> E_cleave
   └─> γ(μ_O)

Step 6: Electronic Properties
   ├─> Ag2O relaxation
   └─> DOS + Bands

Step 7: AIMD
   ├─> Ag2O relaxation
   ├─> Generate + relax slabs
   └─> AIMD simulation
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
- Workflow preset guide: `../../docs/WORKFLOW_PRESETS_GUIDE.md`
- Examples: `../../docs/WORKFLOW_PRESETS_EXAMPLES.md`

---

**Good luck testing your PS-TEROS workflows step by step!** 🚀
