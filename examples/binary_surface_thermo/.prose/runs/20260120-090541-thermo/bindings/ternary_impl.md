# ternary_impl

kind: let

source:
```prose
let ternary_impl = session: code_analyst
  prompt: "Analyze the TERNARY oxide surface thermodynamics implementation in PS-TEROS..."
```

---

## Analysis of Ternary Oxide Surface Thermodynamics Implementation

### 1. Surface Gibbs Free Energy Formula

The implementation in `/home/trevizam/git/PS-TEROS/teros/core/thermodynamics.py` uses the following formula for ternary oxides (lines 75-76):

```
gamma(Delta_mu_M, Delta_mu_O) = phi - Gamma_M * Delta_mu_M - Gamma_O * Delta_mu_O
```

Where:
- **phi**: Reference surface energy at bulk equilibrium (when Delta_mu_M = Delta_mu_O = 0)
- **Gamma_M**: Surface excess of metal M relative to reference metal N
- **Gamma_O**: Surface excess of oxygen relative to reference metal N
- **Delta_mu_M, Delta_mu_O**: Chemical potential deviations from reference state

### 2. Chemical Potential Range Calculation

The Delta_mu ranges are computed based on formation enthalpy bounds (lines 179-180):

```python
# Delta_mu ranges from decomposition limit to element-rich limit
delta_mu_M_range = np.linspace(delta_h / x_M, 0, grid_points)
delta_mu_O_range = np.linspace(delta_h / z_O, 0, grid_points)
```

Where:
- **delta_h**: Formation enthalpy of the oxide (negative for stable compounds)
- **x_M, z_O**: Stoichiometric coefficients in the reduced formula M_x N_y O_z

Physical convention (lines 175-178):
- **Element-rich limit**: Delta_mu = 0 (equilibrium with pure element)
- **Decomposition limit**: Delta_mu = Delta_H_f / stoich (stability boundary)

Since Delta_H_f < 0 for stable oxides, Delta_mu ranges from negative values to 0.

### 3. Reference Energies Definition

Reference energies are extracted from the `reference_energies` dictionary (lines 137-141):

```python
ref_energies = {
    element_M: float(ref_data['metal_energy_per_atom']),
    element_N_ref: float(ref_data['nonmetal_energy_per_atom']),
    element_O: float(ref_data['oxygen_energy_per_atom']),
}
```

From `/home/trevizam/git/PS-TEROS/teros/core/hf.py`:
- **E_M**: Energy per atom of metal M from DFT relaxation of bulk metal
- **E_N**: Energy per atom of metal N from DFT relaxation of bulk metal
- **E_O**: Energy per atom from O2 molecule: `oxygen_energy.value / oxygen_count` (line 142)

This means E_O = (1/2) * E(O2), consistent with the Reuter & Scheffler convention.

### 4. Surface Excess Terms

Surface excesses are calculated relative to the reference metal N (lines 160-161):

```python
# Surface excess relative to N (per unit area, considering both surfaces)
gamma_M = (N_M_slab - (x_M / y_N) * N_N_slab) / (2 * area)
gamma_O = (N_O_slab - (z_O / y_N) * N_N_slab) / (2 * area)
```

Where:
- **N_M_slab, N_N_slab, N_O_slab**: Atom counts in the slab
- **x_M, y_N, z_O**: Stoichiometric coefficients from reduced bulk formula
- **area**: Surface area from cross product of lattice vectors a and b
- **Factor of 2**: Accounts for two surfaces in a slab model

The surface excess Gamma_i represents the deviation from stoichiometric composition normalized by area.

### 5. Formation Enthalpy Integration

Formation enthalpy enters the reference surface energy phi calculation (lines 163-172):

```python
# Reference surface energy phi
# Formula: phi = [E_slab - N_M*E_M - N_N*E_N - N_O*E_O - (N_N/y_N)*Delta_H_f] / (2A)
# This corresponds to gamma at Delta_mu_M = Delta_mu_O = 0 (equilibrium with pure elements)
phi = (
    slab_energy.value
    - N_M_slab * ref_energies[element_M]
    - N_N_slab * ref_energies[element_N_ref]
    - N_O_slab * ref_energies[element_O]
    - (N_N_slab / y_N) * delta_h
) / (2 * area)
```

The term `(N_N/y_N) * Delta_H_f` accounts for the bulk formation energy contribution, ensuring that the reference surface energy phi corresponds to the equilibrium state where Delta_mu_i = 0 for all elements.

### 6. Formation Enthalpy Calculation

From `/home/trevizam/git/PS-TEROS/teros/core/hf.py` (lines 145-158):

```python
# Calculate formation energy
e_bulk = bulk_energy.value
formation_energy = e_bulk
formation_energy -= element_counts[metal_symbol] * metal_energy_per_atom
formation_energy -= element_counts[oxygen_symbol] * oxygen_energy_per_atom

# For ternary, subtract nonmetal contribution
if use_nonmetal:
    nonmetal_energy_per_atom_val = nonmetal_energy.value / nonmetal_count
    formation_energy -= element_counts[nonmetal_symbol] * nonmetal_energy_per_atom_val
```

Formula (from docstring lines 35-36):
```
Delta_H_f = E_bulk - (n_M * E_M/atom + n_N * E_N/atom + n_O * E_O2/atom)
```

### 7. Alternative B-Based Formulation

The implementation also provides an alternative formulation with element B (N) as the independent variable (lines 204-269):

```python
# Surface excesses for B-based formulation (per unit area)
gamma_N = (N_N_stoich - N_N_slab) / (2 * area)  # Gamma_B
gamma_O_B = (N_O_stoich_B - N_O_slab) / (2 * area)  # Gamma_O (with respect to A as reference)

# Reference surface energy for B-based formulation (at Delta_mu_B=0, Delta_mu_O=0)
phi_B = (
    slab_energy.value
    - N_M_slab * ref_energies[element_M]
    - N_N_slab * ref_energies[element_N_ref]
    - N_O_slab * ref_energies['O']
    - (N_M_slab / x_M) * delta_h
) / (2 * area)
```

This allows computing gamma(Delta_mu_B, Delta_mu_O) where B is the independent variable and A is eliminated via bulk equilibrium.

### 8. Key Differences from Binary Oxide Implementation

For binary oxides (lines 348-524), the formula simplifies to:

```
gamma(Delta_mu_O) = phi - Gamma_O * Delta_mu_O
```

The oxygen chemical potential bounds are computed differently (lines 438-442):

```python
# Lower bound (O-poor): decomposition limit
mu_O_min = (bulk_energy.value - x * E_M_ref) / y

# Upper bound (O-rich): from formation energy
mu_O_from_formation = mu_O_min + delta_h / y
mu_O_max = min(0.0, mu_O_from_formation)
```

### Summary Table

| Quantity | Symbol | Formula/Definition |
|----------|--------|-------------------|
| Surface energy | gamma | phi - Gamma_M * Delta_mu_M - Gamma_O * Delta_mu_O |
| Reference energy | phi | [E_slab - sum(N_i * E_i) - (N_ref/y_ref) * Delta_H_f] / (2A) |
| Surface excess M | Gamma_M | (N_M - (x_M/y_N) * N_N) / (2A) |
| Surface excess O | Gamma_O | (N_O - (z_O/y_N) * N_N) / (2A) |
| O reference | E_O | (1/2) * E(O2) |
| Delta_mu_O range | - | [Delta_H_f/z_O, 0] |
| Delta_mu_M range | - | [Delta_H_f/x_M, 0] |
| Formation enthalpy | Delta_H_f | E_bulk - sum(n_i * E_i/atom) |

### Files Analyzed

- `/home/trevizam/git/PS-TEROS/teros/core/thermodynamics.py` - Main surface energy calculations
- `/home/trevizam/git/PS-TEROS/teros/core/hf.py` - Formation enthalpy calculation
