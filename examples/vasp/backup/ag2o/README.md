# Ag2O Complete Example - PS-TEROS WorkGraph

This example demonstrates the complete workflow for calculating surface properties of **Ag2O** (silver oxide) using the PS-TEROS package with AiiDA-WorkGraph.

---

## Overview

**Material:** Ag2O (binary oxide with cuprite structure)

This example showcases ALL features of PS-TEROS:

1. **Bulk and reference structure relaxations**
   - Ag2O bulk (cuprite structure)
   - Ag metal reference
   - O2 oxygen reference

2. **Formation enthalpy calculation**
   - ΔH_f[Ag2O] from elements

3. **Slab generation**
   - (111) Miller index surfaces
   - Multiple terminations with complementary pairs

4. **Slab relaxation with unrelaxed SCF**
   - Unrelaxed single-point calculations
   - Full structural relaxation
   - Relaxation energy: E_relaxed - E_unrelaxed

5. **Cleavage energy calculation**
   - For complementary termination pairs

6. **Surface thermodynamics**
   - Chemical potential sampling (μ_O)
   - Surface energy γ(μ_O) for each termination

---

## File Structure

```
examples/ag2o/
├── README.md                        # This file
├── complete_ag2o_example.py         # Main workflow script
└── structures/
    ├── ag2o.cif                     # Bulk Ag2O structure
    ├── Ag.cif                       # Silver metal reference
    └── O2.cif                       # Oxygen molecule reference
```

---

## Material Properties

### Ag2O (Silver Oxide)
- **Structure:** Cuprite (cubic)
- **Space group:** Pn-3m (224)
- **Lattice parameter:** a ≈ 4.75 Å
- **Formula units per cell:** Z = 2
- **Composition:** Binary oxide (Ag₂O)

### Surface
- **Miller index:** (111)
- **Min slab thickness:** 15.0 Å
- **Min vacuum thickness:** 15.0 Å

---

## Calculation Parameters

### Bulk Ag2O Relaxation
```python
ENCUT = 520 eV
ISMEAR = 0 (Gaussian smearing for semiconductors)
SIGMA = 0.05
IBRION = 2 (Conjugate gradient)
ISIF = 3 (Relax cell and atoms)
EDIFFG = -0.01 eV/Å
```

### Metal Reference (Ag)
```python
ENCUT = 520 eV
ISMEAR = 1 (Methfessel-Paxton for metals)
SIGMA = 0.2
IBRION = 2
ISIF = 3
EDIFFG = -0.01 eV/Å
```

### Oxygen Reference (O2)
```python
ENCUT = 520 eV
ISMEAR = 0
SIGMA = 0.05
IBRION = 2
ISIF = 2 (Relax atoms only)
LREAL = False (Exact for small systems)
EDIFFG = -0.01 eV/Å
```

### Slab Relaxation
```python
ENCUT = 520 eV
ISMEAR = 0
SIGMA = 0.1
IBRION = 2
ISIF = 2 (Relax atoms only, keep cell fixed)
EDIFFG = -0.1 eV/Å (slightly relaxed for surfaces)
LASPH = True (aspherical contributions)
```

### Computational Resources
- **Cores per calculation:** 40
- **Queue:** par40
- **Code:** VASP-VTST-6.4.3@bohr
- **Potential family:** PBE

---

## Running the Example

### Prerequisites

1. **AiiDA profile setup:**
   ```bash
   verdi profile set-default psteros
   verdi status
   ```

2. **Daemon running:**
   ```bash
   verdi daemon start
   ```

3. **Clear Python cache (if needed):**
   ```bash
   find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
   ```

### Execute the Workflow

```bash
source ~/envs/aiida/bin/activate && python complete_ag2o_example.py
```

The script will:
- Load the AiiDA profile
- Set up all calculation parameters
- Build the complete WorkGraph
- Submit it to the daemon
- Print monitoring instructions

---

## Monitoring the Workflow

### Check workflow status:
```bash
verdi process show <PK>
verdi process report <PK>
```

### Monitor all running processes:
```bash
verdi process list
verdi process list -a -p1  # More detailed
```

### Wait for completion:
The workflow may take **several hours** depending on:
- Number of slab terminations generated
- Convergence behavior
- Cluster load

---

## Expected Outputs

Upon successful completion (exit code `[0]`), the workflow produces:

### 1. Bulk and Reference Energies
- `bulk_energy`: Total energy of relaxed Ag2O (eV)
- `metal_energy`: Total energy of relaxed Ag (eV)
- `oxygen_energy`: Total energy of relaxed O2 (eV)

### 2. Formation Enthalpy
- `formation_enthalpy`: ΔH_f[Ag2O] in eV/formula unit

### 3. Slab Structures
- `slab_structures`: All generated (111) terminations with complementary pairs

### 4. Slab Energies
- `unrelaxed_slab_energies`: Single-point SCF energies (eV)
- `slab_energies`: Relaxed slab energies (eV)
- `relaxation_energies`: E_relaxed - E_unrelaxed for each slab (eV)

### 5. Cleavage Energies
- `cleavage_energies`: For complementary termination pairs (eV/Å²)

### 6. Surface Thermodynamics
- `surface_energies`: γ(μ_O) for each termination (eV/Å²)
- Chemical potential range from Ag-rich to O-rich conditions

---

## Workflow Structure

```
Ag2O_Complete_Workflow
│
├── Parallel Relaxations
│   ├── Bulk Ag2O relaxation
│   ├── Metal Ag relaxation
│   └── Oxygen O2 relaxation
│
├── Formation Enthalpy Calculation
│   └── ΔH_f = E[Ag2O] - 2*E[Ag] - 0.5*E[O2]
│
├── Slab Generation
│   └── (111) surfaces with multiple terminations
│
├── For Each Slab Termination:
│   ├── Unrelaxed SCF calculation
│   ├── Full relaxation
│   └── Relaxation energy calculation
│
├── Cleavage Energy Calculation
│   └── For complementary pairs
│
└── Surface Thermodynamics
    ├── Oxide type identification (binary oxide)
    ├── Chemical potential sampling (μ_O)
    └── Surface energy calculation: γ(μ_O)
```

---

## Calculation Flags

All major features are enabled in this example:

```python
relax_slabs = True                    # Enable slab relaxation
compute_relaxation_energy = True      # Calculate E_relax - E_unrelaxed
compute_cleavage = True               # Calculate cleavage energies
compute_thermodynamics = True         # Calculate surface energies
thermodynamics_sampling = 50          # Grid points for μ sampling
```

---

## Validation Criteria

The workflow is considered **successful** when:

1. **Main node returns exit code `[0]`**
   ```bash
   verdi process show <PK>
   # Look for: State: finished [0]
   ```

2. **All expected outputs are present:**
   - Bulk and reference energies
   - Formation enthalpy
   - Slab structures
   - Relaxation energies
   - Cleavage energies
   - Surface energies

3. **Physical reasonableness:**
   - Formation enthalpy is negative (stable compound)
   - Relaxation energies are negative (lowering energy)
   - Surface energies are positive
   - Cleavage energies follow expected trends

---

## Differences from Ag3PO4 Example

This Ag2O example differs from the Ag3PO4 example in several key ways:

1. **Simpler chemistry:**
   - Binary oxide (Ag-O) vs ternary oxide (Ag-P-O)
   - No nonmetal reference needed (only Ag and O2)

2. **Surface orientation:**
   - (111) surfaces for Ag2O
   - (100) surfaces for Ag3PO4

3. **Thermodynamics:**
   - 1D chemical potential space (μ_O only)
   - 2D chemical potential space (μ_O, μ_P) for Ag3PO4

4. **Formation enthalpy:**
   ```
   Ag2O: ΔH_f = E[Ag2O] - 2*E[Ag] - 0.5*E[O2]
   Ag3PO4: ΔH_f = E[Ag3PO4] - 3*E[Ag] - E[P] - 2*E[O2]
   ```

---

## Troubleshooting

### Workflow fails to submit
```bash
verdi status           # Check daemon status
verdi daemon restart   # Restart daemon
verdi code list        # Verify code is available
```

### Calculations exceed walltime
- Increase NSW (more ionic steps)
- Relax EDIFFG convergence criteria
- Use more computational resources

### No slabs generated
- Check Miller indices are valid for the structure
- Reduce min_slab_thickness
- Check bulk relaxation succeeded

### Missing outputs
```bash
verdi process report <PK>  # Check for errors
verdi process show <PK>    # Inspect workflow state
```

---

## References

- **PS-TEROS documentation:** See main repository README
- **AiiDA documentation:** https://aiida.readthedocs.io
- **AiiDA-WorkGraph:** https://aiida-workgraph.readthedocs.io

---

## Notes

- This example uses production-quality parameters suitable for publication
- Calculation times may vary based on cluster availability
- All energies are reported in eV
- All areas are reported in Ų
- The workflow uses automatic k-point generation with spacing = 0.4 Å⁻¹

---

**Last updated:** 2025-10-11
**PS-TEROS version:** feature-relax-energy branch
