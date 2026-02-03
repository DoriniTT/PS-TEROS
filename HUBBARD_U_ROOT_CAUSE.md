# Hubbard U Root Cause: LDAU=False in Ground State

**Date**: 2026-02-02
**Status**: ROOT CAUSE IDENTIFIED ✅

---

## Executive Summary

The ISPIN=2 fix was correct and is working. The real problem is that **the ground state has LDAU=.FALSE., but it should have LDAU=.TRUE. with LDAUU=0**.

When LDAU is OFF in the ground state, then suddenly ON in the NSCF (ICHARG=11) calculations, the charge density is incompatible and the response is incorrect.

---

## Evidence

### Diagnostic Results

**Ground state** (PK 36949):
```python
'ldau': False  # ← WRONG!
'd_occupation': 31.888 e⁻ (4 Ni × ~8 e⁻) # ← Correct
```

**NSCF responses** (ICHARG=11, LDAU=.TRUE.):
```python
V=-0.2: Δn = -0.016 e⁻ (tiny, correct sign)
V=+0.2: Δn = -0.004 e⁻ (tiny, WRONG sign!)
```

**SCF responses** (normal, LDAU=.TRUE.):
```python
V=-0.2: Δn = -0.116 e⁻ (larger, correct sign)
V=+0.2: Δn = +0.100 e⁻ (larger, correct sign)
```

### Result

- **chi (NSCF)** = 0.032 eV⁻¹ (should be χ₀, bare, LARGER)
- **chi_0 (SCF)** = 0.604 eV⁻¹ (should be χ, screened, SMALLER)
- **Magnitudes are backwards!**

---

## According to VASP Wiki

From [https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U):

### Correct Procedure

1. **Ground state**: LDAU=.TRUE., LDAUTYPE=3, LDAUU=0 (U=0)
2. **Response (NSCF)**: LDAU=.TRUE., LDAUTYPE=3, LDAUU=α, ICHARG=11
3. **Response (SCF)**: LDAU=.TRUE., LDAUTYPE=3, LDAUU=α

### Response Magnitudes

- **χ₀ (bare, NSCF)** = LARGER (e.g., 0.50 eV⁻¹ for NiO)
- **χ (screened, SCF)** = SMALLER (e.g., 0.12 eV⁻¹ for NiO)

### Formula

U = χ⁻¹ - χ₀⁻¹

Where χ > χ₀ in magnitude is WRONG! Should have χ₀ > χ.

---

## Why LDAU=False is Wrong

### What Happens

When **LDAU=.FALSE.** in ground state:
- VASP runs standard DFT (no +U correction)
- Charge density and wavefunctions reflect no on-site repulsion
- CHGCAR/WAVECAR are from a "non-interacting" Hamiltonian

When **NSCF with ICHARG=11** restarts from this:
- ICHARG=11 freezes the charge density
- But now LDAU=.TRUE. is suddenly applied
- The frozen charge density is incompatible with +U Hamiltonian
- Only wavefunctions can adjust → minimal response

Result: **NSCF response is artificially suppressed!**

### What Should Happen

When **LDAU=.TRUE., LDAUU=0** in ground state:
- VASP runs DFT+U with U=0 (equivalent to standard DFT for energies)
- But the formalism includes the +U machinery
- CHGCAR/WAVECAR are from a "+U Hamiltonian with U=0"

When **NSCF with ICHARG=11** restarts from this:
- ICHARG=11 freezes the charge density from a +U calculation
- Applying LDAUU=α is now a smooth perturbation
- Charge is frozen but wavefunctions adjust properly
- Result: Proper "bare" response (no charge relaxation)

---

## The Fix

### Change Ground State INCAR

**Before** (WRONG):
```python
{
    'ldau': False,
    'lorbit': 11,
    'lwave': True,
    'lcharg': True,
}
```

**After** (CORRECT):
```python
{
    'ldau': True,
    'ldautype': 3,
    'ldaul': [2, -1],       # d-orbitals for Ni, none for O
    'ldauj': [0.0, 0.0],
    'ldauu': [0.0, 0.0],    # U=0 for ground state
    'lorbit': 11,
    'lwave': True,
    'lcharg': True,
}
```

### Expected Results

**After fix**:
- **χ₀ (NSCF)** ≈ 0.50 eV⁻¹ (large, bare response)
- **χ (SCF)** ≈ 0.12 eV⁻¹ (small, screened response)
- **Ratio**: χ/χ₀ ≈ 0.24 ✓
- **U**: 1/0.12 - 1/0.50 ≈ 8.3 - 2.0 ≈ 6.3 eV ✓

---

## Code Changes Needed

### 1. Update Example Script

File: `examples/lego/hubbard_u/run_hubbard_u_nio.py`

Change ground state stage (lines 88-103):
```python
{
    'name': 'ground_state',
    'type': 'vasp',
    'structure_from': 'relax',
    'incar': {
        **base_incar,
        'nsw': 0,
        'ibrion': -1,
        'ldau': True,           # CHANGED: was False
        'ldautype': 3,          # ADDED
        'ldaul': [2, -1],       # ADDED
        'ldauj': [0.0, 0.0],    # ADDED
        'ldauu': [0.0, 0.0],    # ADDED
        'lorbit': 11,
        'lwave': True,
        'lcharg': True,
    },
    'restart': None,
    'retrieve': ['OUTCAR'],
},
```

### 2. Update Hubbard Response Brick

File: `teros/core/lego/bricks/hubbard_response.py`

**Line 123**: Change validation message:
```python
# OLD:
"(ground state has LDAU=False, response has LDAU=True)"

# NEW:
"(ground state has LDAU=True with LDAUU=0, response has LDAU=True with LDAUU≠0)"
```

**Add validation warning** in `validate_stage()`:
```python
# After line 92:
print(
    f"WARNING: Ground state stage '{ground_state_from}' MUST have:\n"
    f"  - ldau: True\n"
    f"  - ldautype: 3\n"
    f"  - ldaul: [2, -1]  (or appropriate for target species)\n"
    f"  - ldauj: [0.0, ...]\n"
    f"  - ldauu: [0.0, ...]  (U=0 for all species)\n"
    f"Otherwise NSCF response (ICHARG=11) will be incorrect!"
)
```

### 3. Update Documentation

**Files to update**:
- `CLAUDE.md`: Document LDAU=True requirement
- `HUBBARD_U_FIX_SUMMARY.md`: Add this as Fix #4
- `examples/lego/hubbard_u/README.md`: Explain ground state setup

---

## Testing Procedure

### 1. Run Fixed Example

```bash
verdi daemon restart
python examples/lego/hubbard_u/run_hubbard_u_nio_FIXED.py
```

### 2. Verify Ground State

```python
from aiida import orm

# Get ground state calculation
gs_calc = orm.load_node(<ground_state_PK>)
incar = gs_calc.inputs.parameters.get_dict()['incar']

# Check LDAU settings
assert incar['ldau'] == True
assert incar['ldautype'] == 3
assert incar['ldauu'] == [0.0, 0.0]
print("✓ Ground state has correct LDAU=True with LDAUU=0")
```

### 3. Check NSCF Response Magnitude

```python
# After workflow completes
u_result = get_stage_results(<PK>, 'analysis')

chi_0 = u_result['linear_fit']['chi_0_nscf']['slope']  # Should be chi_0
chi = u_result['linear_fit']['chi_scf']['slope']       # Should be chi

print(f"chi_0 (NSCF, bare): {chi_0:.4f} eV⁻¹")
print(f"chi (SCF, screened): {chi:.4f} eV⁻¹")
print(f"Ratio chi/chi_0: {chi/chi_0:.3f}")

# Expected:
# chi_0 ≈ 0.50 eV⁻¹ (large)
# chi ≈ 0.12 eV⁻¹ (small)
# Ratio ≈ 0.24
# U ≈ 6-8 eV
```

### 4. Verify U Value

```python
U = u_result['summary']['hubbard_u_eV']
print(f"U = {U:.3f} eV")

# Expected: 6-8 eV for NiO
assert 5.0 <= U <= 10.0, f"U out of range: {U:.2f} eV"
print("✓ U is in expected range!")
```

---

## Summary

### What Was Wrong

1. ✅ **ISPIN=2 parsing** - FIXED (2026-02-02)
2. ✅ **Chi/chi_0 labels** - FIXED (2026-02-02)
3. ✅ **Sign convention** - FIXED (2026-02-02)
4. ❌ **Ground state LDAU=False** - ROOT CAUSE IDENTIFIED (2026-02-02)

### Why This Matters

The NSCF calculation (ICHARG=11) with LDAU=.TRUE. cannot properly restart from a ground state with LDAU=.FALSE. The charge density is incompatible, causing artificially small (and sometimes wrong-sign) responses.

### The Solution

**Ground state MUST have LDAU=.TRUE. with LDAUU=0**, not LDAU=.FALSE.!

This ensures the charge density and wavefunctions are generated with the +U machinery active (even though U=0), making the NSCF restart compatible.

### Expected Outcome

With this fix:
- NSCF response will be LARGE (χ₀ ≈ 0.50 eV⁻¹, bare)
- SCF response will be SMALL (χ ≈ 0.12 eV⁻¹, screened)
- χ/χ₀ ratio ≈ 0.24 (matches literature!)
- U ≈ 6-8 eV (correct for NiO!)

---

## References

- [VASP Wiki: Calculate U for LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U)
- [VASP Wiki: NiO LSDA+U](https://www.vasp.at/wiki/index.php/NiO_LSDA+U)
- Cococcioni & de Gironcoli, PRB 71, 035105 (2005)
