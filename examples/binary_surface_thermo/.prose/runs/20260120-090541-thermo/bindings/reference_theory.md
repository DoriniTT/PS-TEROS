# reference_theory

kind: let

source:
```prose
let reference_theory = session: literature_expert
  prompt: "Read and analyze the reference article on ab initio thermodynamics of oxide surfaces..."
```

---

## Analysis of Ab Initio Thermodynamics Framework for Oxide Surfaces

**Reference:** K. Reuter and M. Scheffler, "Composition, structure, and stability of RuO2(110) as a function of oxygen pressure," Phys. Rev. B **65**, 035406 (2001).

---

### 1. General Formula for Surface Gibbs Free Energy

The surface free energy gamma(T,p) is defined as the excess Gibbs free energy per unit area required to create the surface:

**Equation (1):**
```
gamma(T,p) = (1/A) [G(T,p,N_M,N_O) - N_M * mu_M(T,p) - N_O * mu_O(T,p)]
```

Where:
- A = surface area
- G(T,p,N_M,N_O) = total Gibbs free energy of the system
- N_M, N_O = number of metal and oxygen atoms
- mu_M, mu_O = chemical potentials of metal and oxygen atoms

For a slab with two equivalent surfaces:

**Equation (2):**
```
gamma(T,p) = (1/2A) [G^slab(T,p,N_M,N_O) - N_M * mu_M(T,p) - N_O * mu_O(T,p)]
```

Using the bulk oxide equilibrium constraint (for binary oxide M_x O_y):

**Equation (3):**
```
mu_M(T,p) + 2*mu_O(T,p) = g^bulk_MO2(T,p)
```

This leads to the **working equation**:

**Equation (4):**
```
gamma(T,p) = (1/2A) [G^slab(T,p,N_M,N_O) - N_M * g^bulk_MO2(T,p) + (2*N_M - N_O) * mu_O(T,p)]
```

**Physical interpretation:** The surface free energy depends only on the oxygen chemical potential, as the metal chemical potential is constrained by bulk equilibrium.

---

### 2. Definition of Delta_mu_O (Oxygen Chemical Potential Deviation)

The oxygen chemical potential mu_O(T,p) is referenced to the total energy of an isolated O2 molecule:

**Reference state (Eq. 7 context):**
```
mu_O(0K, p) = (1/2) E^total_O2 := 0
```

Thus, **Delta_mu_O** represents the deviation from this zero reference:

```
Delta_mu_O = mu_O(T,p) - (1/2) E^total_O2
```

In practice, the paper uses the oxygen-poor limit as a "safe reference" involving only bulk quantities:

**Equation (10):**
```
(1/2)[g^bulk_MO2(0,0) - g^bulk_M(0,0)] < mu_O(T,p) < (1/2)[g^bulk_MO2(0,0) - g^bulk_M(0,0)] + (1/2)*Delta_G_f(0,0)
```

---

### 3. Bounds/Limits for Delta_mu_O

**Oxygen-poor limit (Metal-rich limit):**

**Equation (6):**
```
min[mu_O(T,p)] = (1/2) [g^bulk_MO2(0,0) - g^bulk_M(0,0)]
```

Physical meaning: Below this limit, the oxide would decompose into metal + oxygen gas.

**Oxygen-rich limit:**

**Equation (7):**
```
max[mu_O(T,p)] = (1/2) E^total_O2
```

Physical meaning: Above this limit, gas-phase O2 would start to condense on the sample.

**Range of allowed chemical potentials:**

**Equation (9):**
```
(1/2) * Delta_G_f(0,0) < mu_O(T,p) - (1/2) E^total_O2 < 0
```

For RuO2: Delta_G_f(0,0) = -3.35 eV (calculated), compared to experimental value of -3.19 eV.

---

### 4. Reference State for Oxygen

The zero reference state for mu_O(T,p) is defined as:

**From Section II.E:**
```
mu_O(0K, p) = (1/2) E^total_O2 := 0
```

This means **half the total energy of a free, isolated O2 molecule at T=0 K** serves as the reference.

**Important note:** The authors emphasize using the "oxygen-poor limit as a safe reference" because DFT total energies for extended systems (bulk) are typically more accurate than those for atoms and molecules.

---

### 5. Relationship Between mu_O, Temperature T, and Partial Pressure p_O2

For an ideal gas reservoir, the oxygen chemical potential follows:

**Equation (21):**
```
mu_O(T,p) = mu_O(T,p0) + (1/2) * k*T * ln(p/p0)
```

Where p0 = 1 atm is the standard pressure.

The temperature dependence at standard pressure is given by:

**Equation (22):**
```
mu_O(T,p0) = mu_O^{O-rich}(0K, p0) + (1/2) * Delta_G(Delta_T, p0, O2)

         = (1/2) [H(T,p0,O2) - H(0K,p0,O2)] - (1/2) * T * [S(T,p0,O2) - S(0K,p0,O2)]
```

**Table I provides mu_O(T,p0) values:**

| T (K) | mu_O(T,p0) (eV) | T (K) | mu_O(T,p0) (eV) |
|-------|-----------------|-------|-----------------|
| 100   | -0.08           | 600   | -0.61           |
| 200   | -0.17           | 700   | -0.73           |
| 300   | -0.27           | 800   | -0.85           |
| 400   | -0.38           | 900   | -0.98           |
| 500   | -0.50           | 1000  | -1.10           |

These values are obtained from JANAF thermochemical tables (Ref. 10).

---

### 6. Practical DFT Implementation

**Approximation (Eq. 20) - replacing Gibbs free energies with DFT total energies:**

```
gamma_{O-poor}(T,p) ~ (1/2A) [E^slab(V,N_M,N_O) - (N_O/2) * E^bulk_MO2(V) - (N_M - N_O/2) * E^bulk_M(V)]
```

**Justification for neglecting vibrational contributions:**
- pV term contribution: ~10^-3 meV/A^2 (negligible)
- Vibrational contribution: within +/- 10 meV/A^2 (same order as numerical uncertainty)

**Slope interpretation (Eq. 13):**
```
Slope = (1/2A) * (N_M - N_O/2) * Delta_G_f(0,0)
```

- Stoichiometric termination (N_M = N_O/2): constant gamma(T,p) vs mu_O
- O-excess termination: gamma decreases with increasing mu_O (more favorable at high O pressure)
- O-deficient termination: gamma increases with increasing mu_O (less favorable at high O pressure)

---

### 7. Generalization for Binary Oxides M_x O_y

The formalism generalizes directly to any binary oxide. The key equations become:

**Bulk equilibrium constraint:**
```
x * mu_M(T,p) + y * mu_O(T,p) = g^bulk_{M_x O_y}(T,p)
```

**Surface free energy:**
```
gamma(T,p) = (1/2A) [G^slab - N_M * g^bulk_{M_x O_y}/x + (y*N_M/x - N_O) * mu_O(T,p)]
```

**Oxygen-poor limit:**
```
min[mu_O] = (1/y) [g^bulk_{M_x O_y}(0,0) - x * g^bulk_M(0,0)]
```

---

### 8. Key Physical Insights

1. **Surface termination depends on experimental conditions:** The stable surface structure changes with (T, p). UHV conditions may stabilize different terminations than atmospheric conditions.

2. **Polar surfaces can be stable:** Traditional electrostatic arguments against polar terminations are invalid because surface atoms can deviate significantly from bulk-like states through structural relaxation and electronic rearrangement.

3. **Practical implementation:**
   - Use DFT total energies for bulk and slab
   - Use experimental thermochemical data (JANAF tables) for mu_O(T,p0)
   - Plot gamma vs mu_O with secondary axes showing T at fixed p (or p at fixed T)

4. **Numerical accuracy:** Within +/- 10 meV/A^2 from basis set and neglected vibrations. GGA vs LDA affects absolute values but not crossover points between terminations.

---

### Summary of Key Equations for Implementation

| Quantity | Equation |
|----------|----------|
| Surface energy | gamma = (1/2A)[E^slab - N_M * E^bulk_oxide + (stoich_coeff) * mu_O] |
| mu_O(T,p) | mu_O(T,p0) + (1/2)*k*T*ln(p/p0) |
| O-poor limit | (1/2)[E^bulk_oxide - E^bulk_metal] |
| O-rich limit | (1/2) E^total_O2 |
| Allowed range | (1/2)*Delta_G_f < Delta_mu_O < 0 |
