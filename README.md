# PS-TEROS

**PS-TEROS** (**P**redicting **S**tability of **TER**minations of **O**xide **S**urfaces) is a Python framework tailored for automating *ab initio* surface thermodynamics of **metal oxide surfaces** in [AiiDA](https://www.aiida.net/).

### The Problem

In metal oxides, a single surface orientation rarely exposes just one atomic arrangement; it can exhibit several distinct surface terminations with varying metal-to-oxygen stoichiometries (e.g., stoichiometric, oxygen-poor, or oxygen-rich cuts). *Ab initio* atomistic thermodynamics determines the relative stability of these oxide terminations by coupling slab models to chemical reservoirs of the constituent species—most notably the oxygen reservoir (*Δμ*<sub>O</sub>).

Determining which termination is thermodynamically favored requires coordinating a complex set of interdependent DFT simulations: generating and relaxing multiple oxide slab terminations, calculating matching bulk oxide references, and evaluating gas-phase reference states under strictly identical numerical settings. Managing this multi-structure workflow manually is tedious, error-prone, and hard to reproduce.

### The Solution

PS-TEROS automates the pathway from oxide crystal structures to thermodynamic stability:

- **Oxide Slab Builders:** Programmatically generates bulk oxide references and multiple symmetric/asymmetric surface terminations (e.g., rutile SnO₂).
- **Typed DFT Recipes:** Enforces strict parameter harmony across all oxide terminations, bulk references, and gas reservoirs in **Quantum ESPRESSO** and **VASP**.
- **AiiDA WorkGraphs:** Orchestrates multi-stage workflows (relaxation → static SCF) with bounded job concurrency and full provenance tracking.
- **Pure-Python Thermodynamics:** Evaluates surface free energies (J/m², eV/Å²) and termination phase diagrams directly as a function of oxygen chemical potential.

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
