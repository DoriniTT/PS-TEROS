# binary_impl

kind: let

source:
```prose
let binary_impl = session: code_analyst
  prompt: "Analyze the BINARY oxide surface thermodynamics implementation in PS-TEROS..."
```

---

## Analysis of Binary Oxide Surface Thermodynamics Implementation

### File Location
`/home/trevizam/git/PS-TEROS/teros/core/thermodynamics.py`

### Function: `calculate_surface_energy_binary`
Lines 348-524

---

## 1. Surface Gibbs Free Energy Formula for Binary Oxides

The implementation uses the formula (lines 362-367):

```
gamma(Delta_mu_O) = phi - Gamma_O * Delta_mu_O
```

Where:
- **gamma**: Surface Gibbs free energy (eV/Angstrom^2)
- **phi**: Reference surface energy at O-poor limit
- **Gamma_O**: Surface oxygen excess relative to bulk stoichiometry
- **Delta_mu_O**: Oxygen chemical potential deviation

The actual computation (lines 457-461):
```python
# Compute gamma(Delta_mu_O) - 1D array
gamma_array = []
for mu_O in delta_mu_O_range:
    gamma = phi - Gamma_O * (float(mu_O) - mu_O_min)
    gamma_array.append(float(gamma))
```

**Note**: The implementation uses `(mu_O - mu_O_min)` as the effective Delta_mu_O, meaning the reference point is shifted to the O-poor limit (decomposition limit).

---

## 2. Delta_mu_O Range Calculation

### Lower Bound (O-poor limit): mu_O_min
Lines 437-438:
```python
# Lower bound (O-poor): decomposition limit
mu_O_min = (bulk_energy.value - x * E_M_ref) / y
```

This corresponds to the **decomposition limit** where the oxide would decompose into pure metal M and oxygen. Mathematically:
```
mu_O_min = (E_bulk - x * E_M) / y
```

Where:
- `E_bulk`: Total DFT energy of the bulk M_x O_y structure
- `x`: Number of metal atoms in bulk unit cell
- `E_M`: Energy per atom of elemental metal M
- `y`: Number of oxygen atoms in bulk unit cell

### Upper Bound (O-rich limit): mu_O_max
Lines 440-442:
```python
# Upper bound (O-rich): from formation energy
mu_O_from_formation = mu_O_min + delta_h / y
mu_O_max = min(0.0, mu_O_from_formation)
```

The O-rich limit is calculated as:
```
mu_O_max = min(0, mu_O_min + Delta_H_f / y)
```

Where:
- `Delta_H_f`: Formation enthalpy of the oxide (negative for stable compounds)
- The `min(0.0, ...)` ensures mu_O does not exceed the O2 reference (mu_O = 0 corresponds to equilibrium with O2 gas)

### Range Generation
Lines 444-445:
```python
# Generate chemical potential range
delta_mu_O_range = np.linspace(mu_O_min, mu_O_max, grid_points)
```

---

## 3. Reference Surface Energy (phi) Calculation

Lines 447-452:
```python
# Calculate surface energy at O-poor limit (reference)
phi = (
    slab_energy.value
    - N_M_slab * (bulk_energy.value / x)
    + stoichiometric_imbalance * mu_O_min
) / (2 * area)
```

The formula:
```
phi = [E_slab - N_M * (E_bulk / x) + Delta_O * mu_O_min] / (2 * A)
```

Where:
- `E_slab`: Total DFT energy of the slab
- `N_M`: Number of metal atoms in the slab
- `E_bulk / x`: Bulk energy per metal atom
- `Delta_O`: Stoichiometric imbalance (see below)
- `mu_O_min`: Oxygen chemical potential at decomposition limit
- `A`: Surface area (factor of 2 accounts for two surfaces)

---

## 4. Oxygen Surface Excess (Gamma_O) Calculation

### Stoichiometric Imbalance
Lines 431-434:
```python
# Stoichiometric imbalance: Delta_O = expected_O - actual_O
# Expected oxygen based on metal count and bulk stoichiometry
expected_O = (y / x) * N_M_slab
stoichiometric_imbalance = expected_O - N_O_slab
```

The stoichiometric imbalance:
```
Delta_O = (y/x) * N_M - N_O
```

Where:
- `y/x`: Bulk O/M ratio
- `N_M`: Actual metal atoms in slab
- `N_O`: Actual oxygen atoms in slab

### Surface Excess per Unit Area
Lines 454-455:
```python
# Surface oxygen excess (per unit area)
Gamma_O = stoichiometric_imbalance / (2 * area)
```

Formula:
```
Gamma_O = Delta_O / (2 * A)
```

---

## 5. Summary of Key Variables

| Variable | Symbol | Definition | Code Reference |
|----------|--------|------------|----------------|
| Bulk stoichiometry | x, y | M_x O_y reduced formula | Lines 402-408 |
| Metal reference energy | E_M_ref | DFT energy per atom of elemental M | Line 414 |
| Oxygen reference energy | E_O_ref | DFT energy per atom (1/2 * E_O2) | Line 415 |
| Formation enthalpy | Delta_H | E_bulk - sum(n_i * E_i) | Line 416 |
| Surface area | A | Cross product of cell vectors a,b | Lines 418-424 |
| Stoichiometric imbalance | Delta_O | Expected O - Actual O in slab | Lines 433-434 |
| O-poor limit | mu_O_min | (E_bulk - x*E_M) / y | Line 438 |
| O-rich limit | mu_O_max | min(0, mu_O_min + Delta_H/y) | Lines 441-442 |
| Reference surface energy | phi | Surface energy at O-poor limit | Lines 448-452 |
| Oxygen surface excess | Gamma_O | Delta_O / (2*A) | Line 455 |

---

## 6. Comparison with Reuter & Scheffler Convention

The PS-TEROS implementation differs slightly from the Reuter & Scheffler (PRB 65, 035406) convention:

### Reuter & Scheffler (Eq. 20):
- Delta_mu_O is referenced to (1/2) * E(O2) at T=0K
- O-poor limit: Delta_mu_O = (1/2) * Delta_G_f
- O-rich limit: Delta_mu_O = 0
- Range: (1/2)*Delta_G_f < Delta_mu_O < 0

### PS-TEROS Implementation:
- Uses absolute mu_O values (not deviations)
- mu_O_min corresponds to decomposition limit
- mu_O_max = min(0, mu_O_min + Delta_H_f/y)
- The `(mu_O - mu_O_min)` term in gamma calculation effectively shifts to the decomposition reference

### Equivalence:
Both formulations are mathematically equivalent when the reference point convention is properly accounted for. The PS-TEROS approach directly computes absolute chemical potential bounds from DFT energies, then uses the O-poor limit as the reference point for the linear gamma(mu_O) relationship.

---

## 7. Complete Code Extract

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

    For a binary oxide M_x O_y, the surface energy is calculated as:
    gamma(Delta_mu_O) = phi - Gamma_O * Delta_mu_O

    where:
    - phi: reference surface energy at O-poor limit
    - Gamma_O: surface oxygen excess relative to bulk stoichiometry
    - Delta_mu_O: oxygen chemical potential deviation (from decomposition limit to O2)
    """
    # ... [see full implementation in thermodynamics.py lines 348-524]
```

---

## 8. Output Dictionary Structure

The function returns an `orm.Dict` containing:

```python
{
    'primary': {
        'phi': float,           # Reference surface energy (eV/A^2)
        'Gamma_O': float,       # Oxygen surface excess (1/A^2)
        'delta_mu_O_range': list[float],  # mu_O values sampled
        'gamma_array': list[float],       # gamma at each mu_O
        'gamma_O_poor': float,  # gamma at mu_O_min
        'gamma_O_rich': float,  # gamma at mu_O_max
        'element_M': str,       # Metal element symbol
    },
    'oxide_type': 'binary',
    'area_A2': float,           # Surface area in Angstrom^2
    'bulk_stoichiometry': {...},
    'slab_atom_counts': {...},
    'formation_enthalpy_eV': float,
    'mu_O_min': float,          # O-poor chemical potential
    'mu_O_max': float,          # O-rich chemical potential
    # ... additional legacy keys
}
```
