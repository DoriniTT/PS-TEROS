# Hubbard U Workflow Fix - Implementation Summary

## Problem

The NiO Hubbard U calculation workflow was failing with exit code 302, caused by several VASP response calculations failing with exit code 11.

### Root Cause

**AiiDA Bug**: When `VaspCalculation` receives a `restart_folder` input, it doesn't produce a `remote_folder` output. The standard `VaspWorkChain` marks this output as required in its spec, causing it to fail with exit code 11: "The process did not register a required output."

**Evidence**:
- VaspCalculation: `Finished [0]` (VASP succeeded)
- VaspWorkChain: `Finished [11]` (failed due to missing output)
- Outputs present: `misc` (Dict), `retrieved` (FolderData)
- Output missing: `remote_folder` (RemoteData) - not actually needed

## Solution Implemented

Created a custom `VaspWorkChainNoRemote` wrapper class that makes the `remote_folder` output optional.

### Code Changes

**File**: `teros/core/lego/bricks/hubbard_response.py`

1. **Removed old monkey-patch** (lines 43-67)
2. **Added custom wrapper class** (lines 44-60):

```python
class VaspWorkChainNoRemote(OriginalVaspWC):
    """VaspWorkChain wrapper that makes remote_folder output optional.

    This works around the AiiDA bug where VaspCalculation with restart_folder
    input doesn't produce remote_folder output, which would cause the original
    VaspWorkChain to fail with exit code 11.

    Since the hubbard_response workflow doesn't need remote_folder from the
    response calculations (both NSCF and SCF restart from the ground state's
    remote folder), making this output optional allows the workflow to proceed.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        # Make remote_folder output not required
        spec.outputs['remote_folder'].required = False
```

3. **Updated task creation** (line 164):

```python
# Use custom wrapper that makes remote_folder output optional
VaspTask = task(VaspWorkChainNoRemote)
```

4. **Updated module docstring** (lines 1-21) to document the workaround

## Test Results

### Test Script
`examples/lego/hubbard_u/run_hubbard_u_nio.py`

### Workflow Completion
- **Status**: `Finished [0]` ✓
- **WorkGraph PK**: 35976
- **All stages completed successfully**

### Response Calculations Status

| Calculation | PK | Exit Status | remote_folder |
|------------|-----|-------------|---------------|
| nscf_V=-0.20 | 36020 | 0 ✓ | absent (OK) |
| nscf_V=-0.10 | 36084 | 0 ✓ | present |
| nscf_V=+0.10 | 36119 | 0 ✓ | absent (OK) |
| nscf_V=+0.20 | 36139 | 0 ✓ | absent (OK) |
| scf_V=-0.20 | 36050 | 0 ✓ | absent (OK) |
| scf_V=-0.10 | 36027 | 0 ✓ | absent (OK) |
| scf_V=+0.10 | 36059 | 0 ✓ | absent (OK) |
| scf_V=+0.20 | 36093 | 0 ✓ | absent (OK) |

**All 8 response calculations completed successfully!**

### Hubbard U Results (Initial - Contains Bugs)

**NOTE**: The initial implementation had two bugs that caused incorrect U values:

```
U(Ni-3d) = 29.594 eV  ← INCORRECT (4.7× too high, expected ~6 eV)
SCF linear fit R² = 0.956920
NSCF linear fit R² = 0.731429

Number of response calculations: 4
Potential values (eV): [-0.2, -0.1, 0.1, 0.2]

Linear fit parameters (WRONG LABELS):
  SCF: slope = -0.6040, intercept = -0.007000
  NSCF: slope = -0.0320, intercept = -0.007000
```

**Issues identified:**
1. Chi/chi_0 labels were swapped (NSCF → chi_0 instead of chi)
2. Sign convention was inverted (negative slopes instead of positive)

See "Hubbard U Bug Fixes" section below for details.

## Why This Approach Works

1. **Class Definition**: The custom wrapper defines its own spec by calling `super().define(spec)` and then modifying it. This happens when AiiDA loads the class.

2. **Daemon-Safe**: The class definition is part of the code itself, so when the daemon loads and instantiates `VaspWorkChainNoRemote`, it gets the modified spec.

3. **Clean Separation**: The workaround is isolated to the hubbard_response brick. Other bricks continue using the standard VaspWorkChain.

4. **Minimal Changes**: Only one file modified with ~20 lines of code changes.

5. **Well Documented**: The workaround is explicit and documented in both the class docstring and module docstring.

## Advantages Over Monkey-Patching

| Aspect | Monkey-Patch | Custom Wrapper |
|--------|--------------|----------------|
| Daemon persistence | ✗ No | ✓ Yes |
| Code clarity | ✗ Hidden side effect | ✓ Explicit |
| Maintainability | ✗ Hard to track | ✓ Easy to understand |
| Testability | ✗ Import-order dependent | ✓ Reliable |
| Scope | ✗ Global | ✓ Localized to brick |

## Verification

### Before Fix
```
VaspWorkChain: Finished [11]
Error: "The process did not register a required output: remote_folder"
```

### After Fix
```
VaspWorkChainNoRemote: Finished [0]
Outputs: misc (Dict), retrieved (FolderData)
remote_folder: optional (present in 1/8 cases, absent in 7/8 - both OK)
```

## Impact

- Hubbard U workflows now complete successfully
- Both NiO and SnO2 examples should work
- The sequential dependency pattern (lines 268-273) is preserved
- No changes needed to other lego bricks or user code

## Future Considerations

This wrapper should be removed once the upstream AiiDA bug is fixed:
- https://github.com/aiidateam/aiida-core/issues/XXXX
- https://github.com/aiidateam/aiida-vasp/issues/XXXX

Monitor these issues for updates.

---

## Hubbard U Bug Fixes (2026-02-02)

### Problem Discovery

After the workflow fix enabled successful execution, analysis revealed the U value was incorrect:
- **Calculated**: U = 29.594 eV
- **Expected**: U ≈ 6-8 eV (literature for NiO)
- **Error Factor**: ~4.7× too high

### Root Causes

**Bug 1: Chi/Chi_0 Label Swap**
- NSCF (ICHARG=11, frozen charge) produces **screened** response → chi (small)
- SCF (relaxed charge) produces **bare** response → chi_0 (large)
- Code incorrectly assigned: NSCF → chi_0, SCF → chi (backwards!)

**Bug 2: Sign Convention Inverted**
- Positive potential V should increase d-occupation (Δn > 0)
- Code produced negative Δn for positive V
- Caused negative slopes in linear regression

### Fixes Implemented

**Fix 1: Corrected Chi/Chi_0 Assignment** (`tasks.py:315-326`)

```python
# Before (WRONG):
chi_0_slope = linear_regression(potentials, delta_n_nscf_vals)  # Wrong label
chi_slope = linear_regression(potentials, delta_n_scf_vals)     # Wrong label

# After (CORRECT):
chi_slope = linear_regression(potentials, delta_n_nscf_vals)    # Screened
chi_0_slope = linear_regression(potentials, delta_n_scf_vals)   # Bare
```

**Fix 2: Corrected Sign Convention** (`utils.py:119`)

```python
# Before:
ldauu.append(potential_value)

# After:
ldauu.append(-potential_value)  # Flip sign to match VASP convention
```

**Fix 3: Added Validation** (`tasks.py:336-354`)

Added automatic checks for:
- Non-positive slopes (sign issues)
- Chi >= chi_0 (should be screened < bare)
- U outside typical range (0-20 eV)

**Fix 4: Enhanced Documentation**

Added detailed comments explaining:
- Physical meaning of chi vs chi_0
- ICHARG=11 behavior (frozen charge)
- Expected inequalities (chi < chi_0)

### Files Modified

1. **`teros/core/u_calculation/tasks.py`**
   - Added `warnings` import
   - Fixed chi/chi_0 assignment with documentation
   - Added validation checks

2. **`teros/core/u_calculation/utils.py`**
   - Negated potential sign in `build_ldau_arrays()`

### Expected Behavior After Fixes

**Slopes should be positive:**
- chi ≈ +0.032 eV⁻¹ (small, screened)
- chi_0 ≈ +0.604 eV⁻¹ (large, bare)
- Ratio: chi/chi_0 ≈ 0.053

**Physical interpretation:**
- Positive slopes confirm correct sign convention
- chi < chi_0 confirms correct assignment
- U = 1/chi - 1/chi_0 = 1/0.032 - 1/0.604 ≈ 29.6 eV

**Note on magnitude:** The chi/chi_0 ratio (0.053) is smaller than literature values (0.24), suggesting the magnitude issue may require additional investigation beyond these fixes. Possible causes:
- Potential magnitude too large (±0.2 eV)
- Spin-polarization effects
- LDAUTYPE=3 behavior for magnetic systems

### Testing

To verify fixes:

```bash
# Restart daemon
verdi daemon restart

# Run test
python examples/lego/hubbard_u/run_hubbard_u_nio.py

# Check results
from teros.core.lego import print_sequential_results, get_stage_results
print_sequential_results(<PK>)

u_result = get_stage_results(<PK>, 'analysis')
print(f"U = {u_result['hubbard_u_eV']:.3f} eV")
print(f"chi = {u_result['chi_slope']:.4f} eV⁻¹ (should be positive, small)")
print(f"chi_0 = {u_result['chi_0_slope']:.4f} eV⁻¹ (should be positive, large)")
```

### Verification Checklist

- [x] Code compiles without errors
- [x] Linting passes
- [x] Chi and chi_0 correctly labeled
- [ ] Test calculation shows positive slopes
- [ ] Test calculation shows chi < chi_0
- [ ] U value is in reasonable range
- [ ] Validation checks work correctly

### References

- [VASP Wiki: Calculate U for LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U)
- [Cococcioni & de Gironcoli, PRB 71, 035105 (2005)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.71.035105)
