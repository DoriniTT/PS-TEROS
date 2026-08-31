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

## Quickstart

Generate starting slab models for multiple SnO₂(110) surface terminations directly in Python:

```python
import psteros

# Build stoichiometric and oxygen-deficient SnO2(110) slabs
slab_stoich, id_stoich = psteros.sno2_110_slab(
    termination="o", triple_layers=9, vacuum_angstrom=20.0
)
slab_reduced, id_reduced = psteros.sno2_110_slab(
    termination="sno", triple_layers=9, vacuum_angstrom=20.0
)

print(id_stoich.termination, slab_stoich.composition)    # o    Sn18 O36
print(id_reduced.termination, slab_reduced.composition)  # sno  Sn18 O34
```

Evaluate surface free energies directly from converged DFT outputs:

```python
# Calculate surface energy for an oxide in equilibrium with bulk
point = psteros.surface_energy_oxide_equilibrium(
    slab_energy_ev=-498.0,
    n_metal=18,
    n_oxygen=36,
    bulk_formula_energy_ev=-28.0,
    oxygen_reference_energy_ev=-9.8,
    delta_mu_oxygen_ev=-1.0,
    surface_area_angstrom2=45.2,
    surfaces=2,
)

print(f"{point.gamma_j_per_m2:.3f} J/m²")  # 1.063 J/m²
```

## Installation

```bash
pip install .
```

For configuring AiiDA computers, codes, and pseudopotential families, see the [Installation Guide](docs/source/installation.rst).

## Documentation & Tutorials

- **[First Tutorial](docs/source/tutorial.rst):** Build your first unsubmitted AiiDA WorkGraph.
- **[Core Concepts](docs/source/concepts.rst):** How structures, calculation recipes, execution policies, and provenance connect.
- **[SnO₂ Surface Model](docs/source/examples.rst):** Deep dive into terminations and thermodynamic reference states.
- **[Quantum ESPRESSO Guide](docs/source/qe-first-workflow.rst):** Setting up a two-stage relaxation → static SCF workflow.
- **[API Reference](docs/source/api.rst):** Public classes, functions, and configuration schemas.

