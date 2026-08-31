# reference_theory

kind: let

source:
```prose
let reference_theory = session: theory_expert
  prompt: "Read the reference article at /home/trevizam/PhysRevB.65.035406.pdf
           Extract ALL equations relevant to binary oxide surface thermodynamics..."
```

---

## Complete Extraction of Reuter & Scheffler Thermodynamic Formalism
### From: Phys. Rev. B 65, 035406 (2001)

---

## 1. The General Surface Gibbs Free Energy Formula

### Equation (1) - Fundamental Definition

$$\gamma(T,p) = \frac{1}{A}\left[G(T,p,N_{\text{Ru}},N_{\text{O}}) - N_{\text{Ru}}\mu_{\text{Ru}}(T,p) - N_{\text{O}}\mu_{\text{O}}(T,p)\right]$$

**Definition of terms:**
- $\gamma(T,p)$: Surface free energy (energy per unit area)
- $A$: Surface area
- $G(T,p,N_{\text{Ru}},N_{\text{O}})$: Total Gibbs free energy of the system
- $N_{\text{Ru}}$, $N_{\text{O}}$: Number of metal (Ru) and oxygen atoms in the sample
- $\mu_{\text{Ru}}$, $\mu_{\text{O}}$: Chemical potentials of metal and oxygen atoms

**Physical meaning:** The surface free energy is the excess Gibbs free energy of the system compared to having all atoms in their respective bulk/reservoir phases, normalized by surface area.

### Equation (2) - Slab Model with Two Equivalent Surfaces

$$\gamma(T,p) = \frac{1}{2A}\left[G^{\text{slab}}(T,p,N_{\text{Ru}},N_{\text{O}}) - N_{\text{Ru}}\mu_{\text{Ru}}(T,p) - N_{\text{O}}\mu_{\text{O}}(T,p)\right]$$

**Note:** Factor of 2 appears because slab has TWO equivalent surfaces.

---

## 2. Chemical Potential Constraint - Bulk Equilibrium

### Equation (3) - Thermodynamic Constraint

$$\mu_{\text{Ru}}(T,p) + 2\mu_{\text{O}}(T,p) = g^{\text{bulk}}_{\text{RuO}_2}(T,p)$$

**Definition of terms:**
- $g^{\text{bulk}}_{\text{RuO}_2}(T,p)$: Gibbs free energy per formula unit of bulk oxide

**Physical meaning:** When the surface is in equilibrium with the bulk oxide, the chemical potentials are constrained by the bulk formation reaction. This eliminates one independent variable.

**For a general binary oxide $M_x O_y$:** The constraint becomes:
$$x\mu_M + y\mu_O = g^{\text{bulk}}_{M_x O_y}$$

---

## 3. Surface Free Energy as Function of Oxygen Chemical Potential Only

### Equation (4) - Key Working Equation

$$\gamma(T,p) = \frac{1}{2A}\left[G^{\text{slab}}(T,p,N_{\text{Ru}},N_{\text{O}}) - N_{\text{Ru}}g^{\text{bulk}}_{\text{RuO}_2}(T,p) + (2N_{\text{Ru}} - N_{\text{O}})\mu_{\text{O}}(T,p)\right]$$

**Derivation:** Substituting Eq. (3) into Eq. (2) to eliminate $\mu_{\text{Ru}}$.

**Physical meaning:** After applying the bulk equilibrium constraint, the surface free energy depends on only ONE independent variable: the oxygen chemical potential $\mu_{\text{O}}$.

---

## 4. Bounds of the Oxygen Chemical Potential

### 4.1 Oxygen-Poor Limit (Metal-Rich Limit)

**Physical condition:** If $\mu_{\text{O}}$ becomes too low, the oxide decomposes into metal + O$_2$ gas.

### Equation (5) - Metal Precipitation Condition

$$\max[\mu_{\text{Ru}}(T,p)] = g^{\text{bulk}}_{\text{Ru}}(T,p)$$

### Equation (6) - Oxygen-Poor Limit

$$\min[\mu_{\text{O}}(T,p)] \doteq \frac{1}{2}\left[g^{\text{bulk}}_{\text{RuO}_2}(0,0) - g^{\text{bulk}}_{\text{Ru}}(0,0)\right]$$

**Physical meaning:** At the oxygen-poor limit, the metal chemical potential equals the bulk metal Gibbs free energy. Combined with Eq. (3), this sets the minimum $\mu_{\text{O}}$.

**CRITICAL SIGN CONVENTION:** This is the LOWER bound. The notation $\doteq$ indicates this is a well-defined theoretical reference point.

### 4.2 Oxygen-Rich Limit

**Physical condition:** Beyond this limit, gas-phase O$_2$ would condense on the sample.

### Equation (7) - Oxygen-Rich Limit

$$\max[\mu_{\text{O}}(T,p)] \doteq \frac{1}{2}E^{\text{total}}_{\text{O}_2}$$

**Definition:** $E^{\text{total}}_{\text{O}_2}$ is the total energy of a free, isolated O$_2$ molecule at T=0 K.

**Physical meaning:** At the oxygen-rich limit, the oxygen chemical potential equals half the O$_2$ molecule energy (i.e., the energy of one O atom in the molecule).

---

## 5. Gibbs Free Energy of Formation

### Equation (8) - Definition

$$\Delta G_f(T,p) = g^{\text{bulk}}_{\text{RuO}_2}(T,p) - g^{\text{bulk}}_{\text{Ru}}(T,p) - g^{\text{gas}}_{\text{O}_2}(T,p)$$

**Definition of terms:**
- $\Delta G_f(T,p)$: Gibbs free energy of formation of the oxide
- $g^{\text{gas}}_{\text{O}_2}(T,p)$: Gibbs free energy of an O$_2$ molecule

**Sign convention:** $\Delta G_f < 0$ for stable oxides.

**Example from paper:** $\Delta G_f(0,0) = -3.35$ eV for RuO$_2$

---

## 6. Range of Allowed Chemical Potentials

### Equation (9) - Complete Range

$$\frac{1}{2}\Delta G_f(0,0) < \mu_{\text{O}}(T,p) - \frac{1}{2}E^{\text{total}}_{\text{O}_2} < 0$$

**Defining $\Delta\mu_{\text{O}}$:** If we define:
$$\Delta\mu_{\text{O}} = \mu_{\text{O}} - \frac{1}{2}E^{\text{total}}_{\text{O}_2}$$

Then the bounds become:
- **O-rich limit:** $\Delta\mu_{\text{O}} = 0$
- **O-poor limit:** $\Delta\mu_{\text{O}} = \frac{1}{2}\Delta G_f(0,0)$

**CRITICAL:** Since $\Delta G_f < 0$ for stable oxides, the O-poor limit is at a NEGATIVE value of $\Delta\mu_{\text{O}}$.

---

## 7. Rewritten Bounds Using Only Bulk Quantities

### Equation (10) - Safe Reference Range

$$\frac{1}{2}\left[g^{\text{bulk}}_{\text{RuO}_2}(0,0) - g^{\text{bulk}}_{\text{Ru}}(0,0)\right] < \mu_{\text{O}}(T,p) < \frac{1}{2}\left[g^{\text{bulk}}_{\text{RuO}_2}(0,0) - g^{\text{bulk}}_{\text{Ru}}(0,0)\right] + \frac{1}{2}\Delta G_f(0,0)$$

**Physical meaning:** This rewrites Eq. (9) to avoid explicit dependence on the O$_2$ molecule energy.

---

## 8. Surface Free Energy at Reference Points

### Equation (11) - Oxygen-Poor Limit Reference

$$\gamma_{\text{O-poor}}(T,p) = \frac{1}{2A}\left[G^{\text{slab}}(T,p,N_{\text{Ru}},N_{\text{O}}) - N_{\text{Ru}}g^{\text{bulk}}_{\text{RuO}_2}(T,p) - \left(N_{\text{Ru}} - \frac{N_{\text{O}}}{2}\right)\left[g^{\text{bulk}}_{\text{Ru}}(0,0) - g^{\text{bulk}}_{\text{RuO}_2}(0,0)\right]\right]$$

**Physical meaning:** This is the surface free energy evaluated at the oxygen-poor limit. It involves ONLY bulk and slab quantities - no molecular energies.

### Equation (12) - Relationship to Oxygen-Rich Limit

$$\gamma_{\text{O-rich}}(T,p) = \gamma_{\text{O-poor}}(T,p) - \frac{1}{2A}\left(N_{\text{Ru}} - \frac{N_{\text{O}}}{2}\right)\Delta G_f(0,0)$$

**Physical meaning:** The difference between O-rich and O-poor surface energies depends on the stoichiometry deviation $(N_{\text{Ru}} - N_{\text{O}}/2)$ and the formation enthalpy.

---

## 9. Slope of $\gamma(\Delta\mu_{\text{O}})$ - THE SURFACE EXCESS

### Equation (13) - Slope Definition

$$\text{slope} = \frac{1}{2A}\left(N_{\text{Ru}} - \frac{N_{\text{O}}}{2}\right)\Delta G_f(0,0)$$

**CRITICAL ANALYSIS - Defining the Surface Excess $\Gamma_{\text{O}}$:**

From Eq. (4), the general form is:
$$\gamma = \text{const} + \frac{(2N_{\text{Ru}} - N_{\text{O}})}{2A}\mu_{\text{O}}$$

This can be rewritten as:
$$\gamma = \text{const} - \frac{(N_{\text{O}} - 2N_{\text{Ru}})}{2A}\mu_{\text{O}}$$

**Identifying the surface excess:**

If we write the standard thermodynamic form:
$$\gamma = \phi - \Gamma_{\text{O}} \cdot \Delta\mu_{\text{O}}$$

Then by comparison:
$$\Gamma_{\text{O}} = \frac{N_{\text{O}} - 2N_{\text{Ru}}}{2A} = \frac{N_{\text{O}} - y \cdot N_M / x}{2A}$$

where for $M_x O_y$, the stoichiometric ratio is $y/x$.

---

## 10. SIGN CONVENTION FOR SURFACE EXCESS - DETAILED ANALYSIS

### Definition from Thermodynamics

The surface excess $\Gamma_{\text{O}}$ represents the number of oxygen atoms per unit area IN EXCESS of what would be present in a stoichiometric slab.

### Sign Interpretation

| Surface Type | $N_{\text{O}} - 2N_{\text{Ru}}$ | $\Gamma_{\text{O}}$ | Meaning |
|-------------|--------------------------------|---------------------|---------|
| Stoichiometric | 0 | 0 | No excess oxygen |
| O-rich termination | $> 0$ | $> 0$ (POSITIVE) | Excess oxygen at surface |
| O-poor termination | $< 0$ | $< 0$ (NEGATIVE) | Oxygen deficiency at surface |

### Physical Behavior from Paper (Page 3, bottom)

> "A termination with an O excess (deficiency) will become more favorable (unfavorable) with increasing $\mu_{\text{O}}(T,p)$"

This confirms:
- **O-excess surface** ($\Gamma_{\text{O}} > 0$): $\gamma$ DECREASES as $\mu_{\text{O}}$ increases (negative slope in $\gamma$ vs $\mu_{\text{O}}$)
- **O-deficient surface** ($\Gamma_{\text{O}} < 0$): $\gamma$ INCREASES as $\mu_{\text{O}}$ increases (positive slope)
- **Stoichiometric surface** ($\Gamma_{\text{O}} = 0$): $\gamma$ is CONSTANT (zero slope)

---

## 11. The Linear Relationship $\gamma(\Delta\mu_{\text{O}})$

### Standard Form

$$\gamma(\Delta\mu_{\text{O}}) = \phi - \Gamma_{\text{O}} \cdot \Delta\mu_{\text{O}}$$

**Definition of $\phi$:**
$\phi$ is the surface free energy at the reference chemical potential, typically chosen as:
- $\Delta\mu_{\text{O}} = 0$ (O-rich limit), OR
- Some other convenient reference point

**CRITICAL SIGN CHECK:**

With $\gamma = \phi - \Gamma_{\text{O}} \cdot \Delta\mu_{\text{O}}$:
- For O-rich surface ($\Gamma_{\text{O}} > 0$): slope $d\gamma/d(\Delta\mu_{\text{O}}) = -\Gamma_{\text{O}} < 0$ (negative slope)
- As $\Delta\mu_{\text{O}}$ increases (toward O-rich), $\gamma$ decreases for O-rich surfaces

This is CONSISTENT with Figure 3 in the paper, where:
- RuO$_2$(110)-O$^{\text{cus}}$ (O-excess) has negative slope
- RuO$_2$(110)-O$^{\text{bridge}}$ (stoichiometric) is horizontal
- RuO$_2$(110)-Ru (O-deficient) has positive slope

---

## 12. Summary of Key Equations for Binary Oxide $M_x O_y$

### Generalized Form

**Bulk equilibrium constraint:**
$$x\mu_M + y\mu_O = g^{\text{bulk}}_{M_x O_y}$$

**Surface free energy:**
$$\gamma = \frac{1}{2A}\left[E^{\text{slab}} - \frac{N_M}{x}E^{\text{bulk}}_{M_x O_y} + \left(\frac{y \cdot N_M}{x} - N_O\right)\mu_O\right]$$

**Surface excess:**
$$\Gamma_O = \frac{N_O - (y/x) \cdot N_M}{2A}$$

**Linear form:**
$$\gamma(\Delta\mu_O) = \phi - \Gamma_O \cdot \Delta\mu_O$$

**Bounds of $\Delta\mu_O$:**
- O-rich: $\Delta\mu_O = 0$
- O-poor: $\Delta\mu_O = \frac{1}{y}\Delta H_f$ (where $\Delta H_f$ is formation enthalpy of $M_x O_y$)

---

## 13. Practical DFT Implementation

### Equation (20) - Working Equation for DFT

$$\gamma_{\text{O-poor}}(T,p) \approx \frac{1}{2A}\left[E^{\text{slab}}(V,N_{\text{Ru}},N_{\text{O}}) - \frac{N_{\text{O}}}{2}E^{\text{bulk}}_{\text{RuO}_2}(V) - \left(N_{\text{Ru}} - \frac{N_{\text{O}}}{2}\right)E^{\text{bulk}}_{\text{Ru}}(V)\right]$$

**Approximation:** Neglects vibrational contributions to Gibbs free energy (shown to be < 10 meV/A^2 in the paper).

---

## 14. Temperature and Pressure Dependence

### Equation (21) - Ideal Gas Expression

$$\mu_O(T,p) = \mu_O(T,p^\circ) + \frac{1}{2}kT\ln\left(\frac{p}{p^\circ}\right)$$

**Physical meaning:** For an ideal gas reservoir, the chemical potential depends logarithmically on pressure.

---

## 15. Verification Checklist

1. **Sign of $\Gamma_O$**: POSITIVE for O-excess surfaces, NEGATIVE for O-deficient
2. **Slope of $\gamma$ vs $\Delta\mu_O$**: Equals $-\Gamma_O$
3. **O-rich limit**: $\Delta\mu_O = 0$ (by definition)
4. **O-poor limit**: $\Delta\mu_O = \frac{1}{2}\Delta G_f < 0$ (since $\Delta G_f < 0$ for stable oxides)
5. **Stoichiometric surface**: Horizontal line (constant $\gamma$)
6. **Formation enthalpy**: Negative for thermodynamically stable oxides

---

## 16. IMPORTANT CLARIFICATIONS

### On the Minus Sign in $\gamma = \phi - \Gamma_O \cdot \Delta\mu_O$

The minus sign ensures correct physical behavior:
- When $\Gamma_O > 0$ (excess O) and $\Delta\mu_O$ increases (more O-rich conditions), $\gamma$ decreases
- This makes physical sense: surfaces with excess O become MORE stable when O is abundant

### Alternative Forms in Literature

Some sources write:
$$\gamma = \phi + \Gamma_O \cdot \Delta\mu_O$$

In this case, $\Gamma_O$ is defined with OPPOSITE sign (negative for O-excess). Always check the sign convention!

### The Paper's Convention

The paper uses $(2N_{\text{Ru}} - N_O)$ in Eq. (4), which equals $-2 \cdot \Gamma_O$ (using our definition).

This gives the term: $+\frac{(2N_{\text{Ru}} - N_O)}{2A}\mu_O = -\Gamma_O \cdot \mu_O$

Confirming the form: $\gamma = \phi - \Gamma_O \cdot \Delta\mu_O$
