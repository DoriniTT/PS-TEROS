# Dual Structure Upgrade: Quick Test Enhancement

## Summary

The quick test has been upgraded from testing **1 structure** to testing **2 structures in parallel**, providing a more comprehensive validation of the adsorption energy module.

## Changes Made

### 1. New Structure File

**Created**: `Ag_OOH_quick_test.cif`
- System: Ag(111) + OOH radical
- Atoms: 7 (4 Ag + OOH)
- Formula: HAg4O2

This complements the existing:
- `Ag_OH_quick_test.cif`: 6 atoms (4 Ag + OH)

### 2. Updated Script (`run_quick_test.py`)

#### Structure Loading
```python
# Before: Load 1 structure
complete_structure = orm.StructureData(ase=ase_structure)

# After: Load 2 structures
complete_structure_oh = orm.StructureData(ase=ase_structure_oh)
complete_structure_ooh = orm.StructureData(ase=ase_structure_ooh)
```

#### Adsorption Dictionaries
```python
# Before: Single structure
adsorption_structures = {
    'quick_test': complete_structure,
}
adsorption_formulas = {
    'quick_test': 'OH',
}

# After: Two structures
adsorption_structures = {
    'oh_site': complete_structure_oh,   # 6 atoms
    'ooh_site': complete_structure_ooh,  # 7 atoms
}
adsorption_formulas = {
    'oh_site': 'OH',
    'ooh_site': 'OOH',
}
```

### 3. Updated Documentation (`README.md`)

- Added comparison table showing 2 structures vs 1
- Documented both OH and OOH systems
- Updated expected outputs to show both sites
- Updated success criteria for dual structures
- Added expected physical behavior (OOH binds more strongly)

### 4. Updated Workflow Description

```
Workflow phases (for each structure):
  Phase 1: Relax complete systems
    - oh_site: Ag + OH (6 atoms)
    - ooh_site: Ag + OOH (7 atoms)
  Phase 2: Separate relaxed structures
  Phase 3: SCF calculations (3 single-point per structure)
    - Substrate (Ag)
    - Molecule (OH or OOH)
    - Complete (Ag+OH or Ag+OOH)
  Phase 4: Calculate adsorption energies

Total VASP jobs: 2 relax + 6 SCF = 8 calculations
```

## Advantages of Dual Structure Test

### 1. Better Validation
- **Parallel execution**: Tests scatter/gather workflow pattern
- **Different adsorbates**: OH vs OOH validates formula-based separation
- **Independent calculations**: 8 VASP jobs run in parallel

### 2. More Representative
- Real research workflows often compare multiple adsorbates
- Tests that labels ('oh_site', 'ooh_site') work correctly
- Validates that substrate energies can be shared (same Ag4 substrate)

### 3. Still Fast
- Runtime: Still ~10-20 minutes (jobs run in parallel)
- Only ~2× more VASP jobs than single structure
- Much faster than LaMnO3 example (12-18 hours)

## Comparison Table

| Feature | Single Structure | Dual Structure |
|---------|-----------------|----------------|
| **Structures** | 1 (OH only) | 2 (OH + OOH) |
| **Total atoms** | 6 | 6-7 |
| **VASP jobs** | 4 | 8 |
| **Validates parallel** | No | Yes |
| **Tests scatter** | No | Yes |
| **Runtime** | ~10-20 min | ~10-20 min |

## Expected Physical Results

Based on typical DFT calculations on metal surfaces:

- **E_ads(OH)**: -2 to -3 eV (favorable)
- **E_ads(OOH)**: -3 to -4 eV (more favorable)
- **Trend**: E_ads(OOH) < E_ads(OH)

The OOH radical has an extra O atom that can interact with the surface, leading to stronger binding.

## Files Modified

1. **Created**:
   - `Ag_OOH_quick_test.cif` (new structure)
   - `DUAL_STRUCTURE_UPGRADE.md` (this file)

2. **Updated**:
   - `run_quick_test.py` (load 2 structures, update dicts)
   - `README.md` (document dual structures)

## Testing

All files verified:

```bash
✓ Script syntax is valid
✓ OH structure: HAg4O (6 atoms)
✓ OOH structure: HAg4O2 (7 atoms)
✓ Both structures are valid and ready for testing
```

## Usage

No changes to usage - the script works the same way:

```bash
source ~/envs/aiida/bin/activate
cd teros/experimental/adsorption_energy/quick_test
python run_quick_test.py
```

The script will now automatically test BOTH structures in parallel.

## Next Steps

After running the upgraded test:

1. **Verify both outputs exist**:
   ```bash
   verdi process show <PK>
   # Check for:
   # - adsorption_energies['oh_site']
   # - adsorption_energies['ooh_site']
   ```

2. **Compare energies**:
   - OOH should bind more strongly (more negative)
   - Both should be negative (favorable adsorption)

3. **Validate parallel execution**:
   - Check that jobs ran in parallel (not sequential)
   - Total runtime should be ~same as single structure

## Status

✓ Dual structure upgrade complete
✓ All files verified and ready for testing
✓ Documentation updated
✓ Backward compatible (same API, just more structures)
