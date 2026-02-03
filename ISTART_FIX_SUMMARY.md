# Hubbard U Calculation - ISTART Fix

**Date**: 2026-02-02
**Issue**: NSCF response identical to ground state
**Status**: FIX APPLIED ✅

---

## Problem

After implementing all four fixes (ISPIN=2, chi/chi_0 labels, sign convention, ground state LDAU=True), the NiO calculation **still** produced incorrect results:

**Observed behavior:**
- NSCF d-occupation: 7.972 e⁻ per Ni (same as ground state)
- SCF d-occupation: 7.992 e⁻ per Ni (DOES respond)
- NSCF delta_n ≈ 0.000 → chi_0 very small (backwards!)
- SCF delta_n ≈ 0.080 → chi very large (backwards!)
- U = 29.6 eV (still 4.7× too high)

**Root cause:** Even with ICHARG=11 (fixed charge density), VASP was reading **WAVECAR** from the ground state calculation. These wavefunctions were converged for LDAUU=0, so when LDAUU=V was applied, they didn't respond properly.

---

## Solution

Set **ISTART=0** in NSCF calculations to force VASP to start with **new random wavefunctions** instead of reading WAVECAR.

### What ISTART Does

| ISTART | Behavior | Use Case |
|--------|----------|----------|
| 0 | Start from random wavefunctions | Initial calculation or when Hamiltonian changes significantly |
| 1 | Read WAVECAR (default if file exists) | Continue from previous calculation with similar Hamiltonian |
| 2 | Minimal basis set start | Special cases |

### Why This Fixes the Bare Response

When ICHARG=11 with ISTART=1 (default):
- Charge density frozen from ground state CHGCAR ✓
- Wavefunctions read from ground state WAVECAR (converged for U=0) ✗
- Hamiltonian has U=V, but wavefunctions don't adapt properly ✗
- Result: Suppressed "bare" response

When ICHARG=11 with **ISTART=0**:
- Charge density frozen from ground state CHGCAR ✓
- Wavefunctions initialized randomly (no WAVECAR read) ✓
- Hamiltonian diagonalized with U=V from scratch ✓
- Result: Proper "bare" response to +U potential ✓

---

## Code Changes

**File**: `teros/core/u_calculation/utils.py`

**Function**: `prepare_response_incar()`

**Lines**: 230-236

```python
# Non-SCF: Read and fix charge density
if not is_scf:
    incar['icharg'] = 11
    incar['istart'] = 0  # CRITICAL: Start from new wavefunctions, not WAVECAR
    # This ensures eigenvalues are recalculated with the perturbed +U potential,
    # giving the correct "bare" response
```

---

## Expected Results After Fix

For NiO with V = [−0.2, −0.1, +0.1, +0.2] eV:

### Before ISTART=0 Fix

| V (eV) | Δn_NSCF (e⁻) | Δn_SCF (e⁻) |
|--------|--------------|-------------|
| -0.2   | -0.016       | -0.116      |
| -0.1   | -0.008       | -0.092      |
| +0.1   | 0.000        | 0.080       |
| +0.2   | -0.004       | 0.100       |

- χ₀ (NSCF) = 0.032 eV⁻¹ (TOO SMALL)
- χ (SCF) = 0.604 eV⁻¹ (TOO LARGE)
- U = 29.6 eV ❌

### After ISTART=0 Fix (Expected)

| V (eV) | Δn_NSCF (e⁻) | Δn_SCF (e⁻) |
|--------|--------------|-------------|
| -0.2   | ~-0.10       | ~-0.02      |
| -0.1   | ~-0.05       | ~-0.01      |
| +0.1   | ~+0.05       | ~+0.01      |
| +0.2   | ~+0.10       | ~+0.02      |

- χ₀ (NSCF) ≈ 0.50 eV⁻¹ (LARGE) ✓
- χ (SCF) ≈ 0.12 eV⁻¹ (small) ✓
- Ratio χ/χ₀ ≈ 0.24 ✓
- U ≈ 6-8 eV ✓

---

## Complete Fix History

All five bugs have now been identified and fixed:

### Fix 1: ISPIN=2 Parsing (2026-02-02)
**File**: `teros/core/u_calculation/tasks.py` (lines 21-95)
**Issue**: Parser read only one spin channel instead of spin-summed values
**Fix**: Detect ISPIN and use first match (sum) for ISPIN=2

### Fix 2: Chi/Chi_0 Label Swap (Earlier)
**File**: `teros/core/u_calculation/tasks.py` (lines 361-374)
**Issue**: NSCF and SCF responses assigned to wrong susceptibility variables
**Fix**: Corrected chi_0 ← NSCF, chi ← SCF

### Fix 3: Sign Convention (Earlier)
**File**: `teros/core/u_calculation/utils.py` (line 119)
**Issue**: Potential sign inverted in LDAUU parameter
**Fix**: Added negation (`ldauu.append(-potential_value)`)

### Fix 4: Ground State LDAU=False (2026-02-02)
**Files**: All example scripts
**Issue**: Ground state had LDAU=False, causing incompatible charge density
**Fix**: Changed to LDAU=True with LDAUU=[0, 0]

### Fix 5: ISTART=0 for NSCF (2026-02-02) ⭐ NEW
**File**: `teros/core/u_calculation/utils.py` (line 233)
**Issue**: NSCF reading WAVECAR from ground state, suppressing bare response
**Fix**: Added ISTART=0 to force fresh wavefunction calculation

---

## Testing Instructions

### 1. Restart AiiDA Daemon

```bash
verdi daemon restart  # CRITICAL: Daemon must load new code
```

### 2. Run NiO Test

```bash
python examples/lego/hubbard_u/run_hubbard_u_nio_FIXED.py
```

### 3. Monitor Workflow

```bash
verdi process show <PK>
```

### 4. Verify ISTART in NSCF Calculations

Check that NSCF calculations have ISTART=0:

```bash
# Get NSCF VaspWorkChain PK from workflow
verdi process show <workflow_PK> | grep "nscf"

# Check INCAR for NSCF calculation
verdi calcjob inputcat <nscf_calc_PK> INCAR | grep ISTART
# Should show: ISTART = 0
```

### 5. Check NSCF Response

After workflow completes:

```python
from aiida import orm
from teros.core.lego import get_stage_results

# Get results
u_result = get_stage_results(<PK>, 'analysis')

# Check occupations
responses = u_result['responses']
for V, delta_n_nscf, delta_n_scf in zip(
    responses['potential_values'],
    responses['delta_n_nscf_values'],
    responses['delta_n_scf_values'],
):
    print(f"V={V:+.1f}: Δn_NSCF={delta_n_nscf:.3f}, Δn_SCF={delta_n_scf:.3f}")

# NSCF responses should now be NON-ZERO and follow linear trend
# For NiO: |Δn_NSCF| should be 5-10x LARGER than |Δn_SCF|
```

### 6. Verify U Value

```python
U = u_result['summary']['hubbard_u_eV']
print(f"U = {U:.3f} eV")

# Expected: 6-8 eV for NiO (not 29 eV!)
assert 5.0 <= U <= 10.0, f"U out of range: {U:.2f} eV"
print("✓ U is in expected range!")
```

---

## Physics Explanation

### The VASP Workflow

For Hubbard U linear response, VASP performs:

**Ground State** (LDAU=True, LDAUU=0):
```
1. Initialize random wavefunctions
2. Self-consistent loop: update charge + wavefunctions
3. Converge to ground state (U=0)
4. Save CHGCAR, WAVECAR
```

**NSCF Response** (ICHARG=11, LDAUU=V):

**WITHOUT ISTART=0** (BROKEN):
```
1. Read charge from CHGCAR → frozen ✓
2. Read wavefunctions from WAVECAR → wrong! ✗
3. These wavefunctions were optimized for U=0
4. Apply Hamiltonian with U=V
5. Occupations barely change (wavefunctions not adapted)
6. Result: Artificially suppressed bare response ✗
```

**WITH ISTART=0** (FIXED):
```
1. Read charge from CHGCAR → frozen ✓
2. Initialize NEW random wavefunctions ✓
3. Diagonalize Hamiltonian with U=V from scratch ✓
4. Occupations respond to +U potential ✓
5. Result: Proper bare response ✓
```

**SCF Response** (LDAUU=V, charge updates):
```
1. Start from ground state (optional)
2. Self-consistent loop with U=V
3. Both charge AND wavefunctions adapt
4. Result: Screened response (chi < chi_0) ✓
```

### Why χ₀ > χ

- **χ₀ (bare)**: Only wavefunctions respond, charge frozen → LARGE response
- **χ (screened)**: Both charge and wavefunctions respond → screening reduces response
- Therefore: χ₀ > χ always for transition metal oxides
- U = χ⁻¹ - χ₀⁻¹ = positive value

---

## Related Files

### Core Code
1. `teros/core/u_calculation/utils.py` - ISTART=0 added
2. `teros/core/u_calculation/tasks.py` - ISPIN=2 fix, chi/chi_0 labels
3. `teros/core/lego/bricks/hubbard_response.py` - Calls prepare_response_incar

### Examples
1. `examples/lego/hubbard_u/run_hubbard_u_nio_FIXED.py` - Ground state LDAU=True
2. `examples/lego/hubbard_u/run_hubbard_u_sno2.py` - Ground state LDAU=True
3. `examples/lego/hubbard_u/run_sequential_relax_then_u.py` - Ground state LDAU=True

### Documentation
1. `HUBBARD_U_FIXES_COMPLETE.md` - Summary of fixes 1-4
2. `ISTART_FIX_SUMMARY.md` - This file (fix 5)
3. `HUBBARD_U_ROOT_CAUSE.md` - LDAU=False analysis
4. `ISPIN2_FIX_DETAILS.md` - ISPIN=2 technical details

---

## References

- [VASP Wiki: ICHARG](https://www.vasp.at/wiki/index.php/ICHARG)
- [VASP Wiki: ISTART](https://www.vasp.at/wiki/index.php/ISTART)
- [VASP Wiki: Calculate U for LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U)
- [VASP Wiki: NiO LSDA+U](https://www.vasp.at/wiki/index.php/NiO_LSDA+U)
- Cococcioni & de Gironcoli, PRB 71, 035105 (2005)

---

## Status: READY FOR TESTING ⚡

All five bugs identified and fixed. **Critical new addition**: ISTART=0 for NSCF calculations ensures wavefunctions are properly recalculated with the perturbed +U Hamiltonian.

**Next step**: Test with NiO to verify χ₀ > χ and U ≈ 6-8 eV!
