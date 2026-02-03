# Testing Guide: ISPIN=2 Fix for Hubbard U

**Date**: 2026-02-02
**Status**: Ready for testing ✅

---

## Quick Test

```bash
# 1. Daemon already restarted ✓
verdi daemon status

# 2. Run NiO test
python examples/lego/hubbard_u/run_hubbard_u_nio.py

# 3. Monitor (replace <PK> with actual process ID)
verdi process show <PK>
```

---

## Expected Output Changes

### Console Output

Look for these new messages:
```
INFO: ISPIN=2 detected, using first 'total charge' section (spin-summed)
INFO: ISPIN=2 detected (spin-polarized calculation)
  Extracted total d-occupation: 7.8234
  Per-atom: [7.8234]
```

### D-Occupation Values

**Before fix** (WRONG):
- Ground state: ~3-5 e⁻ per Ni
- Response calculations: similar low values

**After fix** (CORRECT):
- Ground state: ~7-9 e⁻ per Ni
- Response calculations: proper spin-summed values

### Chi/Chi_0 Slopes

**Before fix**:
```
chi = 0.032 eV⁻¹
chi_0 = 0.604 eV⁻¹
Ratio = 0.053
U = 29.6 eV
```

**After fix**:
```
chi ≈ 0.064 eV⁻¹
chi_0 ≈ 0.267 eV⁻¹
Ratio ≈ 0.24
U ≈ 12 eV
```

---

## Verification Script

After workflow completes:

```python
from teros.core.lego import get_stage_results

# Replace with your actual PK
pk = 12345

# Get ground state results
gs_result = get_stage_results(pk, 'ground_state')
print("=== Ground State ===")
print(f"D-occupation (total): {gs_result['vasp']['misc']['total_d_occupation']:.3f} e⁻")
print(f"Per-atom: {gs_result['vasp']['misc']['per_atom_d_occupation']}")
print(f"ISPIN: {gs_result['vasp']['misc']['ispin']}")

# Get U calculation results
u_result = get_stage_results(pk, 'analysis')
print("\n=== Hubbard U Results ===")
print(f"U = {u_result['summary']['hubbard_u_eV']:.3f} eV")
print(f"chi = {u_result['linear_fit']['chi_scf']['slope']:.4f} eV⁻¹")
print(f"chi_0 = {u_result['linear_fit']['chi_0_nscf']['slope']:.4f} eV⁻¹")

# Calculate ratio
chi = u_result['linear_fit']['chi_scf']['slope']
chi_0 = u_result['linear_fit']['chi_0_nscf']['slope']
ratio = chi / chi_0
print(f"Ratio chi/chi_0 = {ratio:.3f}")

# Validation checks
print("\n=== Validation ===")
checks = {
    'ISPIN=2 detected': gs_result['vasp']['misc']['ispin'] == 2,
    'D-occ in range (7-9 e⁻)': 7.0 <= gs_result['vasp']['misc']['per_atom_d_occupation'][0] <= 9.0,
    'Chi positive': chi > 0,
    'Chi_0 positive': chi_0 > 0,
    'Chi < chi_0': chi < chi_0,
    'Ratio near 0.24': 0.20 <= ratio <= 0.30,
    'U reasonable (8-15 eV)': 8.0 <= u_result['summary']['hubbard_u_eV'] <= 15.0,
}

for check, passed in checks.items():
    status = '✓' if passed else '✗'
    print(f"{status} {check}")

# Overall assessment
if all(checks.values()):
    print("\n🎉 All checks passed! ISPIN=2 fix is working correctly.")
else:
    print("\n⚠️  Some checks failed. Review results above.")
```

---

## Success Criteria

✅ All these should be true:
1. ISPIN=2 detection message appears
2. D-occupation per Ni: 7-9 e⁻ (not 3-5 e⁻)
3. Chi slope positive
4. Chi_0 slope positive
5. Chi < chi_0
6. Chi/chi_0 ratio: 0.20-0.30 (not 0.053)
7. U value: 8-15 eV (not 29.6 eV)
8. No validation warnings

---

## If Tests Fail

### Symptom: Still seeing low d-occupations (~3-5 e⁻)

**Cause**: Code changes not loaded
**Fix**:
```bash
verdi daemon restart
verdi daemon status  # Verify running
```

### Symptom: U still ~29 eV

**Cause**: May be using cached/old results
**Fix**: Run fresh calculation, don't reuse previous nodes

### Symptom: Ratio still ~0.05

**Cause**: ISPIN not being detected correctly
**Fix**: Check OUTCAR manually
```python
from aiida import orm
calc = orm.load_node(36001)  # Replace with your GS calc PK
outcar = calc.outputs.retrieved.get_object_content('OUTCAR')
import re
ispin = re.search(r'ISPIN\s*=\s*(\d)', outcar)
print(f"ISPIN = {ispin.group(1) if ispin else 'NOT FOUND'}")
```

### Symptom: Workflow fails

**Cause**: Check daemon logs
```bash
verdi daemon logshow
```

---

## Fine-Tuning (If Needed)

If U is still outside 6-8 eV target but ratio is correct (~0.24):

### Option 1: Smaller Potentials

Edit `examples/lego/hubbard_u/run_hubbard_u_nio.py`:
```python
potential_values = [-0.1, -0.05, 0.05, 0.1]  # Instead of [-0.2, -0.1, 0.1, 0.2]
```

### Option 2: Better k-points

```python
incar = {
    'encut': 520,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
}
kpoints_spacing = 0.025  # Instead of 0.03 (finer grid)
```

### Option 3: Proper AFM Ordering

NiO is antiferromagnetic. Set initial magnetic moments:
```python
# For 2-atom NiO unit cell
structure.set_initial_magnetic_moments([5.0, -5.0, 0.0, 0.0])  # Ni↑, Ni↓, O, O
```

---

## Documentation Updated

✅ `HUBBARD_U_FIX_SUMMARY.md` - Updated with ISPIN=2 fix
✅ `ISPIN2_FIX_DETAILS.md` - Comprehensive technical documentation
✅ `teros/core/u_calculation/tasks.py` - Code comments and docstrings
✅ `TEST_ISPIN2_FIX.md` - This testing guide

---

## Next Steps

1. Run test calculation
2. Verify all checks pass
3. If successful, update examples with best practices
4. Consider extending to other test systems (CoO, FeO)

---

## Contact

If issues persist after following this guide, check:
- `verdi daemon logshow` for errors
- `verdi process report <PK>` for workflow details
- OUTCAR files manually for ISPIN detection
