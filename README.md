# PS-TEROS

**PS-TEROS** (**P**redicting **S**tability of **TER**minations of **O**xide **S**urfaces) is a Python framework for automating *ab initio* surface thermodynamics workflows in [AiiDA](https://www.aiida.net/).

### The Problem

Under realistic operating conditions (temperature and oxygen pressure, $T$ and $p_{\text{O}_2}$), metal oxides expose multiple possible surface terminations. Determining which termination is thermodynamically stable requires calculating the surface free energy $\gamma(T, p_{\text{O}_2})$. In practice, this demands coordinating dozens of interdependent DFT simulations—bulk references, gas reservoirs, and various slab terminations—that all require strictly identical numerical settings. Managing this web of calculations manually is fragile, tedious, and hard to reproduce.

### The Solution

PS-TEROS automates the path from crystal structure to thermodynamic stability:

- **Slab Builders:** Programmatically generates bulk references and symmetric/asymmetric oxide slabs (e.g., SnO₂).
- **Typed DFT Recipes:** Enforces parameter harmony across bulk, slab, and gas calculations in **Quantum ESPRESSO** and **VASP**.
- **AiiDA WorkGraphs:** Orchestrates multi-stage workflows (relaxation $\to$ static SCF) with bounded job concurrency and full provenance tracking.
- **Pure-Python Thermodynamics:** Calculates surface free energies (J/m², eV/Å²) and phase diagrams directly from converged energies.

## Start here

- **Install the library:** [installation guide](docs/source/installation.rst)
- **Build one graph without submitting work:** [first tutorial](docs/source/tutorial.rst)
- **Understand the pieces of a calculation:** [core concepts](docs/source/concepts.rst)
- **Understand the SnO2 surface-energy model:**
  [model explanation](docs/source/examples.rst)
- **Prepare a relaxation followed by a static QE calculation:**
  [QE guide](docs/source/qe-first-workflow.rst)
- **Look up public names and inputs:** [reference](docs/source/api.rst)

## A safe first interaction

The structure helpers can be used before you connect to a scheduler or submit a
calculation. This example creates an in-memory starting model for one rutile
SnO2(110) surface; it does not use an AiiDA profile or external compute time.

```python
import psteros

slab, identity = psteros.sno2_110_slab(
    termination="o",
    triple_layers=9,
    vacuum_angstrom=20.0,
)

print(identity.termination)  # o
print(slab.composition.reduced_formula)  # SnO2
```

A starting structure is not yet a converged result. The SnO2 model page explains
what a slab and a termination are, which calculations provide the needed
energies, and which assumptions must be checked before interpreting a surface
energy.

## Install

Follow the [installation guide](docs/source/installation.rst) for the source
checkout, a clean virtual environment, package verification, and the AiiDA
computer, code, and pseudopotential setup needed before the tutorial.
