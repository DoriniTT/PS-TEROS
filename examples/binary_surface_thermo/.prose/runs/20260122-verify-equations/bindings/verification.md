# verification

kind: let

source:
```prose
let verification = session: equation_verifier
  prompt: "Compare the theoretical equations..."
  context: { reference_theory, code_implementation }
```

---

# DETAILED EQUATION VERIFICATION REPORT

## Summary of Findings

| Component | Verdict | Issue |
|-----------|---------|-------|
| **PHI calculation** | **INCORRECT** | Missing term: $+\frac{y N_M}{x} E_O^{\text{ref}} / (2A)$ |
| **Γ_O calculation** | **CORRECT** | Matches theory |
| **γ(Δμ_O) formula** | **CORRECT** | Sign convention correct |
| **Δμ_O range** | **CORRECT** | Uses reduced stoichiometry |

## 1. PHI CALCULATION - CRITICAL ERROR

### Correct Formula (Reuter & Scheffler)

Starting from:
$$\gamma = \frac{1}{2A}\left[E_{\text{slab}} - \frac{N_M E_{\text{bulk}}}{x} + \left(\frac{y N_M}{x} - N_O\right) \mu_O\right]$$

At Δμ_O = 0, where $\mu_O = E_O^{\text{ref}}$:
$$\phi = \frac{1}{2A}\left[E_{\text{slab}} - \frac{N_M E_{\text{bulk}}}{x} + \frac{y N_M}{x} E_O^{\text{ref}} - N_O E_O^{\text{ref}}\right]$$

### Current Code (INCORRECT)

```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    - N_O_slab * E_O_ref           # <-- Only this oxygen term present
) / (2 * area)
```

### Missing Term

The code is MISSING: `+ (y / x) * N_M_slab * E_O_ref`

### Corrected Code

```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    + (y / x) * N_M_slab * E_O_ref   # <-- ADD THIS TERM
    - N_O_slab * E_O_ref
) / (2 * area)
```

Or equivalently (factoring the oxygen terms):
```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    - (N_O_slab - (y / x) * N_M_slab) * E_O_ref  # Combined oxygen term
) / (2 * area)
```

## 2. Impact of the Bug

The error introduces a slab-dependent offset:
$$\text{Error} = \frac{y N_M}{x} \cdot \frac{E_O^{\text{ref}}}{2A}$$

Since $E_O^{\text{ref}} \approx -4.9$ eV (half of O2 energy) and the missing term should be **added**, the current code:

- **Underestimates** φ for all slabs
- The error is **larger for slabs with more metal atoms**
- This **artificially favors metal-rich terminations**

This explains why T1 appears more stable than T2 - if T1 has more metal atoms, it gets a bigger (more negative) error, making it appear more stable.

## 3. Verification of Other Components

### Γ_O - CORRECT ✓
```python
stoichiometric_imbalance = expected_O - N_O_slab  # = (y/x)*N_M - N_O
Gamma_O = -stoichiometric_imbalance / (2 * area)  # = (N_O - (y/x)*N_M) / (2A)
```
This matches theory exactly.

### γ(Δμ_O) formula - CORRECT ✓
```python
gamma = phi - Gamma_O * delta_mu_O
```
Uses correct minus sign per Reuter & Scheffler convention.

### Δμ_O range - CORRECT ✓
```python
delta_mu_O_min = delta_h / y_reduced  # O-poor limit (negative)
delta_mu_O_max = 0.0                   # O-rich limit
```
Correct bounds using reduced stoichiometry.

## 4. Physical Interpretation of the Fix

After fixing phi, the surface energy will be:
$$\gamma(\Delta\mu_O) = \phi_{\text{correct}} - \Gamma_O \cdot \Delta\mu_O$$

Where:
- O-rich terminations (Γ_O > 0) have **negative slope** → become **more stable** at O-rich (Δμ_O → 0)
- O-poor terminations (Γ_O < 0) have **positive slope** → become **more stable** at O-poor (Δμ_O → negative)
- Stoichiometric terminations (Γ_O = 0) have **zero slope** → **constant** stability

## 5. Recommendation

**FIX REQUIRED in `/home/trevizam/git/PS-TEROS/teros/core/thermodynamics.py`**

Change lines 456-462 from:
```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    - N_O_slab * E_O_ref
) / (2 * area)
```

To:
```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    + (y / x) * N_M_slab * E_O_ref  # Missing term for bulk O contribution
    - N_O_slab * E_O_ref
) / (2 * area)
```
