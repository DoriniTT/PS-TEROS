# Quick Test: Dual Adsorption Energy (OH + H2O on Ag)

## Overview

This is a **minimal, fast test** for the adsorption energy module that tests **TWO structures in parallel**. It uses the same workflow as `run_lamno3_oh_adsorption.py` but with much smaller structures and reduced computational parameters for rapid testing.

**Key advantage**: Tests parallel execution with multiple adsorbates (OH and H2O) to validate the scatter/gather workflow pattern.

## Comparison

| Feature | LaMnO3 Example | Quick Test (Dual) |
|---------|---------------|-------------------|
| **Structures tested** | 1 | 2 (parallel) |
| **Total atoms** | 150 (La32Mn28O89H) | 6-7 (Ag4OH + Ag4O2H) |
| **ENCUT** | 400 eV | 250 eV |
| **K-points** | 0.6 Å⁻¹ | 1.0 Å⁻¹ |
| **Spin polarization** | Yes (ISPIN=2) | No (ISPIN=1) |
| **Max relax steps** | 500 | 20 |
| **Convergence iterations** | 5 | 2 |
| **VASP jobs** | 4 (1 relax + 3 SCF) | 8 (2 relax + 6 SCF) |
| **Expected runtime** | 12-18 hours | 10-20 minutes |

## Structures

### Structure 1: OH on Ag(111)
- **System**: Ag(111) surface + OH radical
- **Ag slab**: 2×2 single-layer (4 atoms)
- **Adsorbate**: OH (2 atoms)
- **Total**: 6 atoms
- **File**: `Ag_OH_quick_test.cif`

### Structure 2: H2O on Ag(111)
- **System**: Ag(111) surface + H2O radical
- **Ag slab**: 2×2 single-layer (4 atoms)
- **Adsorbate**: H2O (3 atoms)
- **Total**: 7 atoms
- **File**: `Ag_H2O_quick_test.cif`

**Cell (both)**: 5.78 × 5.78 × 25.0 Å with PBC

## Features Tested

Despite being minimal, this test validates ALL key features:

1. **Parallel adsorption energy workflow** (4 phases × 2 structures)
   - Phase 1: Relax complete systems (2 parallel)
   - Phase 2: Structure separation (automatic)
   - Phase 3: SCF calculations (6 parallel: 3 per structure)
   - Phase 4: Energy calculations (2 values)

2. **Full builder API**
   - `adsorption_relax_builder_inputs` with `relax_settings`
   - `adsorption_scf_builder_inputs`

3. **Scatter/gather pattern**
   - Multiple structures processed in parallel
   - Independent adsorbate formulas (OH vs H2O)
   - Parallel VASP job execution

3. **Atom fixing**
   - `adsorption_fix_atoms=True`
   - `adsorption_fix_type='bottom'`
   - `adsorption_fix_thickness=2.0`

## Usage

```bash
# Make sure you're in the psteros profile
verdi profile list

# If not, switch to it
verdi profile set-default psteros

# Check daemon is running
verdi status
verdi daemon start  # if needed

# Run the test
source ~/envs/aiida/bin/activate
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy/teros/experimental/adsorption_energy/quick_test
python run_quick_test.py
```

## Monitoring

```bash
# Check status
verdi process list -a -p1

# Show details (use PK from submission)
verdi process show <PK>
verdi process status <PK>
verdi process report <PK>
```

## Expected Timeline

1. **Submission**: Immediate
2. **Phase 1 (Relax)**: 5-15 minutes
   - Relaxes **2 systems in parallel**:
     - Ag4OH (6 atoms)
     - Ag4O2H (7 atoms)
   - Max 20 ionic steps each
3. **Phase 2 (Separate)**: < 1 minute
   - Splits each into substrate, molecule, complete:
     - OH site: Ag4, OH, Ag4OH
     - H2O site: Ag4, H2O, Ag4O2H
4. **Phase 3 (SCF)**: 2-5 minutes each
   - **6 parallel single-point calculations**:
     - 3 for OH site (substrate, molecule, complete)
     - 3 for H2O site (substrate, molecule, complete)
5. **Phase 4 (Calculate)**: < 1 minute
   - Computes E_ads for both sites

**Total**: ~10-20 minutes (vs 12-18 hours for LaMnO3)

## Expected Outputs

After successful completion, you should see **outputs for BOTH structures**:

### OH Site Outputs
```
1. relaxed_complete_structures['oh_site']: Relaxed Ag4OH
2. separated_structures['oh_site']:
   - substrate: Ag4 (4 atoms)
   - molecule: OH (2 atoms)
   - complete: Ag4OH (6 atoms)
3. substrate_energies['oh_site']: E(Ag4)
4. molecule_energies['oh_site']: E(OH)
5. complete_energies['oh_site']: E(Ag4+OH)
6. adsorption_energies['oh_site']: E_ads for OH (eV)
```

### H2O Site Outputs
```
1. relaxed_complete_structures['h2o_site']: Relaxed Ag4O2H
2. separated_structures['h2o_site']:
   - substrate: Ag4 (4 atoms)
   - molecule: H2O (3 atoms)
   - complete: Ag4O2H (7 atoms)
3. substrate_energies['h2o_site']: E(Ag4)
4. molecule_energies['h2o_site']: E(H2O)
5. complete_energies['h2o_site']: E(Ag4+H2O)
6. adsorption_energies['h2o_site']: E_ads for H2O (eV)
```

## Success Criteria

The test is successful when:

1. Main node returns `[0]` (finished successfully)
2. All 4 phases complete without errors **for BOTH structures**
3. Final output contains:
   - `adsorption_energies['oh_site']`
   - `adsorption_energies['h2o_site']`
4. Both E_ads values are reasonable negative values (favorable adsorption)
5. Typical values: E_ads(H2O) < E_ads(OH) (H2O binds more strongly)

## Notes

- **This is for TESTING only** - parameters are intentionally minimal
- For production calculations, use parameters from `run_lamno3_oh_adsorption.py`
- The small system size means energies may not be physically accurate
- Focus is on validating the **workflow**, not getting accurate results

## Troubleshooting

If the test fails:

1. Check daemon is running: `verdi daemon status`
2. Check VASP code exists: `verdi code list`
3. Check potential family: `verdi data core.potcar listfamilies`
4. Check process report: `verdi process report <PK>`

## Next Steps

After this test succeeds:

1. Try the full LaMnO3 example for production-quality calculations
2. Modify parameters for your specific system
3. Explore other workflow presets in PS-TEROS
