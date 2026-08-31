# PS-TEROS

**PS-TEROS** (**P**redicting **S**tability of **TER**minations of **O**xide **S**urfaces) is a Python framework for automating and managing *ab initio* surface thermodynamics workflows in [AiiDA](https://www.aiida.net/).

## The Scientific Question

Oxide surfaces are rarely static or single-state systems. Under realistic operating conditions—such as varying temperatures, oxygen partial pressures, or chemical environments—metal oxides expose a multitude of possible surface terminations, ranging from stoichiometric cuts to oxygen-deficient or oxygen-rich reconstructions.

Determining which surface termination is stable under given environmental conditions requires *ab initio* atomistic thermodynamics: calculating the grand canonical surface free energy $\gamma(T, p_{\mathrm{O}_2})$ as a function of the chemical potential of the constituent species.

## The Computational Challenge

Evaluating surface stability is never a single calculation. It requires a carefully coordinated network of interdependent simulations:

1. **Diverse Slab Geometries:** Generating symmetric and asymmetric slabs across multiple crystallographic facets and terminations, ensuring sufficient thickness to recover bulk-like interiors and adequate vacuum padding to prevent spurious periodic interactions.
2. **Strict Reference Consistency:** Computing bulk reference energies per formula unit and gas-phase chemical potential reservoirs (e.g., $\text{O}_2, \text{H}_2\text{O}$). A subtle mismatch in DFT parameters (exchange-correlation functionals, cutoff energies, k-point sampling, or pseudopotentials) between bulk and slab runs will corrupt the chemical potential differences and produce unphysical phase boundaries.
3. **Multi-Stage Orchestration & Provenance:** Submitting and tracking multi-step DFT stages (relaxation $\to$ static high-precision SCF) across HPC clusters while managing job dependencies, queue limits, and calculation provenance.
4. **Thermodynamic Assembly:** Combining all converged total energies, atom counts, surface areas, and chemical potentials into thermodynamic stability equations.

Managing this workflow manually—through disjointed shell scripts, manual geometry generation, and separate post-processing spreadsheets—is tedious, fragile, and difficult to reproduce.

## How PS-TEROS Solves It

PS-TEROS bridges the gap between atomic crystal structures, high-throughput DFT calculations, and surface thermodynamics:

- **Programmatic Slab Building:** Generates consistent bulk references and stoichiometric/non-stoichiometric slab models (e.g., rutile $\text{SnO}_2(110)$ and other oxide facets) directly in Python.
- **Typed Calculation Recipes:** Provides verified calculation recipes for **Quantum ESPRESSO** and **VASP** that enforce strict numerical and physical parameter harmony across all related bulk, slab, and gas-phase calculations.
- **Inspectable AiiDA WorkGraphs:** Orchestrates multi-stage workflows as explicit directed acyclic graphs in AiiDA, featuring built-in concurrency control (preventing HPC queue overload) and automated data provenance.
- **Decoupled Thermodynamic Analysis:** Offers pure-Python analysis helpers to directly calculate surface free energies ($\text{J/m}^2$, $\text{eV/\AA}^2$) and map surface phase diagrams from calculation outputs.

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
