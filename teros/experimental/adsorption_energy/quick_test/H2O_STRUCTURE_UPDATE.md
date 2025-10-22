# Structure Update: OOH → H2O

## Issue

The OOH (hydroperoxyl radical) structure was unstable and not relaxing properly in VASP. This is expected because:
- OOH is a radical with an unpaired electron
- Radicals can have complex electronic structure
- May require spin-polarized calculations
- Geometry optimization can be challenging

## Solution

Replaced OOH with **H2O (water)** molecule:
- H2O is a closed-shell, stable molecule
- Well-characterized on metal surfaces
- Easier to converge in DFT calculations
- Still provides good test of dual-structure workflow

## Changes Made

### 1. Created New Structure File

**New**: `Ag_H2O_quick_test.cif`
- Formula: H2Ag4O (7 atoms: 4 Ag + H2O)
- H2O geometry:
  - O-H bond length: 0.96 Å (correct)
  - H-O-H angle: 104.5° (correct water angle)
- O atom positioned 12.5 Å above Ag surface
- H atoms positioned symmetrically with correct angle

**Removed**: `Ag_OOH_quick_test.cif`

### 2. Updated Script (`run_quick_test.py`)

All references updated:
```python
# Variable names
structure_file_ooh → structure_file_h2o
ase_structure_ooh → ase_structure_h2o
complete_structure_ooh → complete_structure_h2o

# Dictionary keys
'ooh_site' → 'h2o_site'

# Formula strings
'OOH' → 'H2O'
```

### 3. Updated Documentation

- `README.md`: All OOH references → H2O
- Docstring in script updated
- Print statements updated

## Comparison: OOH vs H2O

| Feature | OOH (old) | H2O (new) |
|---------|-----------|-----------|
| **Stability** | Radical (unstable) | Closed-shell (stable) |
| **Electrons** | Unpaired electron | All paired |
| **Spin** | Requires spin-polarized | Works with ISPIN=1 |
| **DFT convergence** | Challenging | Easy |
| **Physical relevance** | Less common | Very common |
| **Atoms** | 3 (O-O-H) | 3 (H-O-H) |

## Why H2O is Better for Testing

1. **Stability**: Closed-shell molecule, no unpaired electrons
2. **Convergence**: Reliably converges with minimal parameters
3. **Relevance**: Water adsorption is a benchmark test in surface science
4. **Complements OH**: Both oxygen-containing, different bonding
5. **Well-studied**: Extensive literature on H2O/metal interfaces

## Physical Expectations

Based on typical DFT results on Ag surfaces:

- **E_ads(OH)**: -2 to -3 eV (strong binding via Ag-O bond)
- **E_ads(H2O)**: -0.3 to -0.8 eV (weaker physisorption)
- **Trend**: E_ads(OH) << E_ads(H2O) (OH binds much more strongly)

OH forms a chemical bond with the surface, while H2O typically physisorbs (weaker interaction).

## Files Modified

1. **Created**:
   - `Ag_H2O_quick_test.cif` (new stable structure)
   - `H2O_STRUCTURE_UPDATE.md` (this file)

2. **Removed**:
   - `Ag_OOH_quick_test.cif` (unstable structure)

3. **Updated**:
   - `run_quick_test.py` (all OOH → H2O references)
   - `README.md` (documentation updated)

## Verification

All files verified and ready:

```
✓ Script syntax is valid
✓ Structure verification:
  OH:  HAg4O (6 atoms)
  H2O: H2Ag4O (7 atoms)
✓ Both structures are valid and ready for testing
```

## Usage

No changes to usage - the script works the same way:

```bash
source ~/envs/aiida/bin/activate
cd teros/experimental/adsorption_energy/quick_test
python run_quick_test.py
```

The workflow will now test:
1. **oh_site**: Ag + OH (chemisorption, strong binding)
2. **h2o_site**: Ag + H2O (physisorption, weak binding)

This provides a better test of the workflow with two very different binding modes!

## Status

✓ OOH structure removed
✓ H2O structure created with correct geometry
✓ Script updated and syntax verified
✓ Documentation updated
✓ Ready for testing
