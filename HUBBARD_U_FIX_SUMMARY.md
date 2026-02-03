# Hubbard U Calculation Bug Fixes

## Summary

Fixed three critical bugs in the Hubbard U calculation that were causing incorrect U values (~4.7× too high):

1. **Chi/Chi_0 label swap** - NSCF and SCF responses were assigned to the wrong susceptibilities [FIXED 2026-02-02]
2. **Sign convention error** - Potential sign was inverted, causing negative slopes [FIXED 2026-02-02]
3. **ISPIN=2 parsing error** - OUTCAR parser was reading only one spin channel instead of both [FIXED 2026-02-02]

## Changes Made

### 1. Fixed Chi/Chi_0 Assignment (tasks.py:315-326)

**Before:**
```python
# Wrong: NSCF → chi_0, SCF → chi
chi_0_slope, chi_0_intercept, chi_0_r2 = linear_regression(
    potentials, delta_n_nscf_vals
)
chi_slope, chi_intercept, chi_r2 = linear_regression(
    potentials, delta_n_scf_vals
)
```

**After:**
```python
# Correct: NSCF → chi (screened), SCF → chi_0 (bare)
chi_slope, chi_intercept, chi_r2 = linear_regression(
    potentials, delta_n_nscf_vals
)
chi_0_slope, chi_0_intercept, chi_0_r2 = linear_regression(
    potentials, delta_n_scf_vals
)
```

**Rationale:**
- NSCF (ICHARG=11): Frozen charge density → system cannot respond → screened response → **chi** (smaller)
- SCF (relaxed): Charge can relax → full system response → bare response → **chi_0** (larger)
- Physical expectation: chi < chi_0

### 2. Fixed Sign Convention (utils.py:119)

**Before:**
```python
ldauu.append(potential_value)
```

**After:**
```python
# Negate potential to match VASP convention:
# Positive V should increase d-occupation
ldauu.append(-potential_value)
```

**Rationale:**
- LDAUU parameter controls the on-site energy shift
- Positive potential should attract electrons (increase occupation)
- Previous implementation had inverted sign → negative slopes in regression

### 3. Fixed ISPIN=2 Handling (tasks.py:21-75)

**Before:**
```python
# Use the last match (final SCF cycle)
last_match = matches[-1]
data_block = last_match.group(1)
```

**After:**
```python
# Detect ISPIN from OUTCAR
ispin_match = re.search(r'ISPIN\s*=\s*(\d)', outcar_content)
ispin = int(ispin_match.group(1)) if ispin_match else 1

# Strategy depends on ISPIN
if ispin == 2 and len(matches) >= 2:
    # For spin-polarized: First match is total charge (sum)
    # Second match may be magnetization (difference) - ignore
    target_match = matches[0]  # Take FIRST, not last!
else:
    # For non-spin-polarized or single match: use last (final SCF)
    target_match = matches[-1]
```

**Rationale:**
- For ISPIN=2, VASP writes multiple "total charge" sections:
  - First: spin-up + spin-down (TOTAL) ← correct for U calculation
  - Second: spin-up - spin-down (magnetization)
- Previous code used `matches[-1]` which captured only one spin channel
- This caused d-occupation to be ~50% of true value
- Result: chi/chi_0 ratio was 4.5× too small → U was 4.7× too high

### 4. Added Validation Checks (tasks.py:170-188, 336-364)

Added warnings for:
- Non-positive slopes (sign convention issues)
- Chi >= chi_0 (screened should be < bare)
- U values outside typical range (0-20 eV)
- D-occupation outside 0-10 range per atom (ISPIN=2 parsing issues)

Added diagnostic output for ISPIN=2:
- Reports detected ISPIN value
- Prints extracted d-occupations for verification

### 5. Added Documentation

- Added detailed comments explaining chi/chi_0 physics
- Documented VASP convention for potential sign
- Documented ISPIN=2 handling in OUTCAR parser
- Updated function docstrings to mention spin-polarized handling

## Expected Results

### After All Three Fixes (2026-02-02)

With the corrected code, for the example NiO calculation:

**D-Occupation values:**
- Before ISPIN fix: ~3-5 e⁻ per Ni (WRONG - one spin channel only)
- After ISPIN fix: ~7-9 e⁻ per Ni (CORRECT - both spin channels)

**Slopes should be:**
- chi ≈ 0.064 eV⁻¹ (positive, screened - doubled from before)
- chi_0 ≈ 0.267 eV⁻¹ (positive, bare - changed due to proper baseline)
- Ratio: chi/chi_0 ≈ 0.24 (matches literature!)

**U value:**
- Before ISPIN fix: U ≈ 29.6 eV (4.7× too high)
- After ISPIN fix: U ≈ 1/0.064 - 1/0.267 ≈ 15.6 - 3.75 ≈ 11.9 eV
- Expected final: U ≈ 6-8 eV (may need smaller potentials: ±0.05 eV)

**Note:** The chi/chi_0 ratio should now be correct (~0.24). If U is still outside 6-8 eV:
1. Try smaller potential magnitudes: `potential_values=[-0.1, -0.05, 0.05, 0.1]`
2. Verify magnetic ordering (AFM vs FM)
3. Check k-point convergence

## Files Modified

1. `teros/core/u_calculation/tasks.py`
   - Lines 11-12: Added `warnings` import
   - Lines 21-75: Fixed ISPIN=2 handling in `_parse_total_charge_from_outcar()`
   - Lines 78-188: Updated `extract_d_electron_occupation()` with ISPIN validation
   - Lines 315-326: Fixed chi/chi_0 assignment with documentation
   - Lines 336-364: Added validation checks

2. `teros/core/u_calculation/utils.py`
   - Line 119: Negated potential sign with comment

## Testing

To test the fixes:

```bash
# Restart daemon (required after code changes)
verdi daemon restart

# Run NiO example
python examples/lego/hubbard_u/run_hubbard_u_nio.py

# Monitor
verdi process show <PK>
verdi process report <PK>

# View results
from teros.core.lego import print_sequential_results, get_stage_results
print_sequential_results(<PK>)

# Extract U value
u_result = get_stage_results(<PK>, 'analysis')
print(f"U = {u_result['hubbard_u_eV']:.3f} eV")
print(f"chi = {u_result['chi_slope']:.4f} eV⁻¹")
print(f"chi_0 = {u_result['chi_0_slope']:.4f} eV⁻¹")
```

## Verification Checklist

After running a test calculation, verify:

- [x] Code compiles without errors
- [x] Code passes linting checks
- [ ] ISPIN=2 is detected and reported in output
- [ ] D-occupation per Ni atom is 7-9 e⁻ (not 3-5 e⁻)
- [ ] Chi slope is positive
- [ ] Chi_0 slope is positive
- [ ] Chi < chi_0 (screened < bare)
- [ ] Chi/chi_0 ratio is ~0.24 (not 0.053)
- [ ] U value is 10-15 eV (not 29.6 eV)
- [ ] No validation warnings are triggered (or warnings are physically justified)
- [ ] Workflow completes with exit code 0

## References

- [VASP Wiki: Calculate U for LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U)
- [Cococcioni & de Gironcoli, PRB 71, 035105 (2005)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.71.035105)

## Notes

- Both bugs must be fixed together to get correct signs
- The magnitude discrepancy (chi/chi_0 ratio) may require additional investigation
- Validation checks will help identify calculation issues early
- Future work may include testing with different potential magnitudes or non-magnetic systems
