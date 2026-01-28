# recommendations

kind: output

source:
```prose
output recommendations = session: literature_expert
  prompt: "Based on the comparison analysis..."
```

---

## Executive Summary

The PS-TEROS binary oxide implementation has a **critical error** in the chemical potential convention. The code computes **absolute** oxygen chemical potentials (mu_O) but labels them as **deviations** (Delta_mu_O). The surface energy VALUES are mathematically correct, but the x-axis (chemical potential range) is on the wrong scale.

---

## 1. Detailed Issue Analysis

### Issue 1: Absolute vs Deviation Chemical Potential

**What is wrong (Reference: Reuter & Scheffler PRB 65, 035406, Eq. 7-9):**

The standard convention defines:
```
Delta_mu_O = mu_O(T,p) - (1/2) E^total_O2
```

Where (1/2) E^total_O2 is the DFT energy of half an O2 molecule, set as the **zero reference**.

The physical bounds are:
- **O-rich limit:** Delta_mu_O = 0 (equilibrium with O2 gas)
- **O-poor limit:** Delta_mu_O = (1/y) * Delta_H_f (decomposition into metal)

For SnO2 (y=2): Delta_mu_O should range from approximately **-2.9 eV to 0 eV**.

**What the code currently does (thermodynamics.py, lines 437-445):**

```python
# Lower bound (O-poor): decomposition limit
mu_O_min = (bulk_energy.value - x * E_M_ref) / y

# Upper bound (O-rich): from formation energy
mu_O_from_formation = mu_O_min + delta_h / y
mu_O_max = min(0.0, mu_O_from_formation)

# Generate chemical potential range
delta_mu_O_range = np.linspace(mu_O_min, mu_O_max, grid_points)
```

The formula `mu_O_min = (E_bulk - x * E_M_ref) / y` computes an **absolute** chemical potential:
```
mu_O^abs = (E_MxOy - x*E_M) / y
```

This is NOT the deviation from O2 reference. The observed values (-7.4 to -8.9 eV for SnO2) are absolute mu_O values, not Delta_mu_O.

**What it should do instead:**

The code should compute deviations directly from the formation enthalpy, as the ternary implementation correctly does (lines 179-180):

```python
delta_mu_O_range = np.linspace(delta_h / z_O, 0, grid_points)
```

---

### Issue 2: Surface Energy Reference Point

**What is wrong (Reference: Eq. 20 in Reuter & Scheffler):**

The surface energy at the Gibbs-Dividing-Surface formulation is:
```
gamma = (1/2A) * [E_slab - N_M * (E_bulk/x) - N_O * mu_O]
     = (1/2A) * [E_slab - N_M * (E_bulk/x) - N_O * ((1/2)E_O2 + Delta_mu_O)]
```

**What the code currently does (lines 447-460):**

```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    + stoichiometric_imbalance * mu_O_min
) / (2 * area)

Gamma_O = stoichiometric_imbalance / (2 * area)

gamma = phi - Gamma_O * (float(mu_O) - mu_O_min)
```

The code uses `mu_O_min` as the reference point for phi, which is mathematically consistent internally but incompatible with the standard convention where **Delta_mu_O = 0** is the reference.

**What it should do instead:**

Use the O2 energy as the reference (E_O_ref = (1/2) E_O2), and compute:
```python
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    - stoichiometric_imbalance * E_O_ref
) / (2 * area)

gamma = phi - Gamma_O * delta_mu_O
```

This gives phi as the surface energy at Delta_mu_O = 0 (O-rich limit).

---

## 2. Corrected Equations

### Standard Reuter & Scheffler Convention for Binary Oxide M_x O_y

**Chemical potential definition:**
$$\Delta\mu_O = \mu_O(T,p) - \frac{1}{2}E_{O_2}^{DFT}$$

**Thermodynamic bounds:**
$$\frac{\Delta H_f}{y} \leq \Delta\mu_O \leq 0$$

Where:
- O-rich limit: $\Delta\mu_O = 0$ (equilibrium with O2 gas at T=0K reference)
- O-poor limit: $\Delta\mu_O = \Delta H_f / y$ (decomposition into metal + oxide)

**Surface energy formula:**
$$\gamma(\Delta\mu_O) = \frac{1}{2A}\left[ E_{slab} - N_M \frac{E_{bulk}}{x} - N_O \left(\frac{1}{2}E_{O_2} + \Delta\mu_O\right) \right]$$

Rearranging:
$$\gamma(\Delta\mu_O) = \phi - \Gamma_O \cdot \Delta\mu_O$$

Where:
- $\phi = \frac{1}{2A}\left[ E_{slab} - N_M \frac{E_{bulk}}{x} - N_O \cdot \frac{1}{2}E_{O_2} \right]$ is the surface energy at $\Delta\mu_O = 0$
- $\Gamma_O = \frac{N_O - (y/x)N_M}{2A}$ is the surface oxygen excess per unit area

**Note:** The stoichiometric imbalance is defined as:
$$\Delta N_O = N_O^{expected} - N_O^{actual} = \frac{y}{x}N_M - N_O$$

So $\Gamma_O = -\Delta N_O / (2A)$, and:
$$\gamma(\Delta\mu_O) = \phi + \frac{\Delta N_O}{2A} \cdot \Delta\mu_O$$

---

## 3. Specific Code Changes

### Change 1: Fix chemical potential range (lines 436-445)

**Current code:**
```python
# Oxygen chemical potential bounds
# Lower bound (O-poor): decomposition limit
mu_O_min = (bulk_energy.value - x * E_M_ref) / y

# Upper bound (O-rich): from formation energy
mu_O_from_formation = mu_O_min + delta_h / y
mu_O_max = min(0.0, mu_O_from_formation)

# Generate chemical potential range
delta_mu_O_range = np.linspace(mu_O_min, mu_O_max, grid_points)
```

**Corrected code:**
```python
# Delta_mu_O bounds (Reuter & Scheffler convention)
# O-poor limit: decomposition into metal
# Delta_mu_O_min = Delta_H_f / y (negative for stable oxide)
delta_mu_O_min = delta_h / y

# O-rich limit: equilibrium with O2 gas (reference state)
delta_mu_O_max = 0.0

# Generate chemical potential range (O-poor to O-rich)
delta_mu_O_range = np.linspace(delta_mu_O_min, delta_mu_O_max, grid_points)
```

### Change 2: Fix reference surface energy phi (lines 447-452)

**Current code:**
```python
# Calculate surface energy at O-poor limit (reference)
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    + stoichiometric_imbalance * mu_O_min
) / (2 * area)
```

**Corrected code:**
```python
# Calculate reference surface energy phi at Delta_mu_O = 0 (O-rich limit)
# phi = (1/2A) * [E_slab - N_M*(E_bulk/x) - N_O*(1/2)*E_O2]
# Note: stoichiometric_imbalance = expected_O - actual_O = (y/x)*N_M - N_O
# So: phi = (1/2A) * [E_slab - N_M*(E_bulk/x) - N_O*E_O_ref]
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    - N_O_slab * E_O_ref
) / (2 * area)
```

### Change 3: Fix gamma calculation (lines 457-461)

**Current code:**
```python
# Compute gamma(Delta_mu_O) - 1D array
gamma_array = []
for mu_O in delta_mu_O_range:
    gamma = phi - Gamma_O * (float(mu_O) - mu_O_min)
    gamma_array.append(float(gamma))
```

**Corrected code:**
```python
# Compute gamma(Delta_mu_O) - 1D array
# gamma = phi - Gamma_O * Delta_mu_O
gamma_array = []
for delta_mu_O in delta_mu_O_range:
    gamma = phi - Gamma_O * float(delta_mu_O)
    gamma_array.append(float(gamma))
```

### Change 4: Update special values (lines 463-465)

**Current code:**
```python
# Special values
gamma_O_poor = gamma_array[0]  # At mu_O_min
gamma_O_rich = gamma_array[-1]  # At mu_O_max
```

**Corrected code (no change needed, just verify labels):**
```python
# Special values
gamma_O_poor = gamma_array[0]   # At delta_mu_O_min (decomposition limit)
gamma_O_rich = gamma_array[-1]  # At delta_mu_O_max = 0 (O2 equilibrium)
```

### Change 5: Remove obsolete variables from output dict

Remove `mu_O_min` and `mu_O_max` from the output dictionary (lines 506-507) as they are no longer meaningful. Replace with:

```python
'delta_mu_O_min': float(delta_mu_O_min),  # = delta_h / y
'delta_mu_O_max': float(delta_mu_O_max),  # = 0.0
```

---

## 4. Summary of Findings

### Surface Energy VALUES (gamma): **CORRECT**

The gamma values themselves (0.06 - 0.30 J/m^2 for SnO2) are physically reasonable for oxide surfaces. The mathematical formula is internally consistent.

### Chemical Potential AXIS: **INCORRECT**

The x-axis labeled "delta_mu_O_eV" contains absolute chemical potential values (-7.4 to -8.9 eV) instead of deviations (should be -2.9 to 0 eV).

**Relationship between incorrect and correct values:**
```
Correct Delta_mu_O = Incorrect "delta_mu_O" - (1/2) E_O2
                   = (-7.4 eV) - (-4.9 eV)
                   ~ -2.5 eV  (close to expected O-poor limit)
```

---

## 5. Physical Interpretation

After correction:

| Condition | Delta_mu_O (eV) | Physical Meaning |
|-----------|-----------------|------------------|
| O-rich | 0 | Equilibrium with O2 gas at T=0K |
| O-poor | Delta_H_f/y ~ -2.9 | Metal oxide just stable vs. decomposition |

The temperature scale in plotting scripts should then correctly map:
- Delta_mu_O = 0 corresponds to very low T (O2 stable)
- Delta_mu_O = Delta_H_f/y corresponds to high T (reducing conditions)

---

## 6. Verification After Fix

After implementing the changes, verify:

1. **Delta_mu_O range:** Should be approximately [-3, 0] eV for SnO2
2. **gamma at Delta_mu_O = 0:** This is the O-rich surface energy
3. **Consistency with ternary:** Set one metal as reference in ternary; results should match binary
4. **Slope sign:** For O-excess slabs (Gamma_O > 0), gamma should decrease as Delta_mu_O increases (O-richer conditions stabilize O-excess surfaces)

---

## 7. Complete Corrected Function

For reference, here is the complete corrected `calculate_surface_energy_binary` function with all changes applied:

```python
@task.calcfunction
def calculate_surface_energy_binary(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Dict,
    sampling: orm.Int,
) -> orm.Dict:
    """
    Compute gamma(Delta_mu_O) surface energy for a single binary oxide slab.

    Following Reuter & Scheffler convention (PRB 65, 035406):
    - Delta_mu_O = mu_O - (1/2)E_O2 (deviation from O2 reference)
    - O-rich limit: Delta_mu_O = 0
    - O-poor limit: Delta_mu_O = Delta_H_f / y

    Surface energy formula:
    gamma(Delta_mu_O) = phi - Gamma_O * Delta_mu_O

    where phi is the surface energy at Delta_mu_O = 0 (O-rich limit).
    """
    grid_points = sampling.value
    if grid_points <= 0:
        raise ValueError('Sampling must be a positive integer')

    # Extract bulk composition
    bulk_ase = bulk_structure.get_ase()
    bulk_counts = Counter(bulk_ase.get_chemical_symbols())

    if 'O' not in bulk_counts:
        raise ValueError('The bulk structure contains no oxygen; expected a binary oxide.')

    metal_elements = [element for element in bulk_counts if element != 'O']
    if len(metal_elements) != 1:
        raise ValueError(f'Expected exactly one metal species; found: {metal_elements}')

    element_M = metal_elements[0]
    element_O = 'O'

    # Get stoichiometry (x and y in M_x O_y)
    x = bulk_counts[element_M]
    y = bulk_counts[element_O]

    # Get reduced stoichiometry
    common_divisor = gcd(x, y)
    x_reduced = x // common_divisor
    y_reduced = y // common_divisor

    # Extract reference energies and formation enthalpy
    ref_data = reference_energies.get_dict()
    formation_data = formation_enthalpy.get_dict()

    E_M_ref = float(ref_data['metal_energy_per_atom'])
    E_O_ref = float(ref_data['oxygen_energy_per_atom'])  # = (1/2) E_O2
    delta_h = float(formation_data['formation_enthalpy_ev'])

    # Calculate surface area
    slab_ase = slab_structure.get_ase()
    cell = slab_ase.get_cell()
    a_vec = cell[0]
    b_vec = cell[1]
    cross = np.cross(a_vec, b_vec)
    area = float(np.linalg.norm(cross))

    # Slab atom counts
    slab_counts = Counter(slab_ase.get_chemical_symbols())
    N_M_slab = slab_counts.get(element_M, 0)
    N_O_slab = slab_counts.get(element_O, 0)

    # Stoichiometric imbalance: expected_O - actual_O
    expected_O = (y / x) * N_M_slab
    stoichiometric_imbalance = expected_O - N_O_slab

    # ========== CORRECTED: Delta_mu_O bounds (Reuter & Scheffler convention) ==========
    # O-poor limit: decomposition into metal
    # Delta_mu_O_min = Delta_H_f / y (negative for stable oxide)
    delta_mu_O_min = delta_h / y_reduced

    # O-rich limit: equilibrium with O2 gas (reference state)
    delta_mu_O_max = 0.0

    # Generate chemical potential range (O-poor to O-rich)
    delta_mu_O_range = np.linspace(delta_mu_O_min, delta_mu_O_max, grid_points)

    # ========== CORRECTED: Reference surface energy phi at Delta_mu_O = 0 ==========
    # phi = (1/2A) * [E_slab - N_M*(E_bulk/x) - N_O*(1/2)*E_O2]
    phi = (
        slab_energy.value
        - N_M_slab * (bulk_energy.value / x)
        - N_O_slab * E_O_ref
    ) / (2 * area)

    # Surface oxygen excess (per unit area)
    # Gamma_O = (N_O - (y/x)*N_M) / (2A) = -stoichiometric_imbalance / (2A)
    Gamma_O = -stoichiometric_imbalance / (2 * area)

    # ========== CORRECTED: Compute gamma(Delta_mu_O) ==========
    gamma_array = []
    for delta_mu_O in delta_mu_O_range:
        gamma = phi - Gamma_O * float(delta_mu_O)
        gamma_array.append(float(gamma))

    # Special values
    gamma_O_poor = gamma_array[0]   # At delta_mu_O_min
    gamma_O_rich = gamma_array[-1]  # At delta_mu_O_max = 0

    # Calculate bulk energy per formula unit
    formula_units_in_bulk = bulk_counts[element_M] / x_reduced
    bulk_energy_per_fu = bulk_energy.value / formula_units_in_bulk

    return orm.Dict(
        dict={
            'primary': {
                'phi': float(phi),
                'Gamma_O': float(Gamma_O),
                'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
                'gamma_array': gamma_array,
                'gamma_O_poor': float(gamma_O_poor),
                'gamma_O_rich': float(gamma_O_rich),
                'gamma_at_reference': float(phi),  # gamma at Delta_mu_O = 0
                'element_M': element_M,
            },
            'oxide_type': 'binary',
            'area_A2': float(area),
            'bulk_stoichiometry': {
                f'x_{element_M}': int(x_reduced),
                'y_O': int(y_reduced),
            },
            'slab_atom_counts': {
                f'N_{element_M}': int(N_M_slab),
                'N_O': int(N_O_slab),
            },
            'reference_energies_per_atom': {
                element_M: float(E_M_ref),
                'O': float(E_O_ref),
            },
            'E_slab_eV': float(slab_energy.value),
            'E_bulk_per_fu_eV': float(bulk_energy_per_fu),
            'formation_enthalpy_eV': float(delta_h),
            'stoichiometric_imbalance': float(stoichiometric_imbalance),
            'delta_mu_O_min': float(delta_mu_O_min),
            'delta_mu_O_max': float(delta_mu_O_max),
            # Legacy keys for backward compatibility
            'phi': float(phi),
            'Gamma_O': float(Gamma_O),
            'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
            'gamma_array': gamma_array,
            'gamma_O_poor': float(gamma_O_poor),
            'gamma_O_rich': float(gamma_O_rich),
            'gamma_at_reference': float(phi),
            'element_M': element_M,
            'E_bulk_eV': float(bulk_energy.value),
            'bulk_stoichiometry_MxOy': {
                f'x_{element_M}': int(x_reduced),
                f'y_O': int(y_reduced),
            },
        }
    )
```

---

## 8. Impact on Existing Results

**For existing calculations:**
- The gamma VALUES do not need to be recalculated
- The x-axis values can be corrected post-hoc by:
  ```python
  corrected_delta_mu_O = reported_delta_mu_O - (1/2) * E_O2_DFT
  ```
  Where E_O2_DFT is the total energy of the O2 molecule used in the formation enthalpy calculation.

**For future calculations:**
- After implementing the fix, new results will be on the correct scale
- Plots will show Delta_mu_O from approximately -3 to 0 eV (for SnO2)
- Temperature mapping in plotting scripts will work correctly
