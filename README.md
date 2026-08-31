# PS-TEROS

**PS-TEROS** (**P**redicting **S**tability of **TER**minations of **O**xide **S**urfaces) is a Python framework tailored for automating *ab initio* surface thermodynamics of **metal oxide surfaces** in [AiiDA](https://www.aiida.net/).

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.apsusc.2025.164350-blue)](https://doi.org/10.1016/j.apsusc.2025.164350)

> **Citation:** If you use PS-TEROS in your research, please cite:  
> T. T. Dorini, M. A. San-Miguel, *Accelerating rational design of oxide surfaces: The PS-TEROS workflow for automated surface stability analysis*, **Applied Surface Science** (2025). [doi:10.1016/j.apsusc.2025.164350](https://doi.org/10.1016/j.apsusc.2025.164350)

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

### 1. Prepare Oxide Slabs (from pymatgen or ASE)

PS-TEROS accepts structures directly from **pymatgen** or **ASE** for any metal oxide (e.g., TiO₂, ZnO, Fe₂O₃, SnO₂):

```python
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator

# Load any bulk oxide and generate surface terminations
bulk = Structure.from_file("bulk_oxide.cif")  # Or from ASE: Structure.from_ase_atoms(atoms)
generator = SlabGenerator(bulk, miller_index=(1, 1, 0), min_slab_size=12.0, min_vacuum_size=15.0)
slabs = generator.get_slabs()

# Map labelled structures for the workflow
structures = {
    "bulk": bulk,
    **{f"term_{i}": slab for i, slab in enumerate(slabs)},
}
```

### 2. Build the AiiDA WorkGraph

Define typed calculation and scheduler settings, then build the inspectable WorkGraph:

```python
import psteros

# Configure DFT recipe (Quantum ESPRESSO or VASP)
calc_config = psteros.QeCalculationConfig(
    code_label="pw-7.2@cluster",
    pseudo_family="SSSP/1.3/PBE/efficiency",
    parameters={
        "CONTROL": {"calculation": "relax"},
        "SYSTEM": {"ecutwfc": 60.0, "ecutrho": 480.0},
        "ELECTRONS": {"conv_thr": 1.0e-8},
    },
    kpoints_distance=0.25,
)

exec_policy = psteros.ExecutionPolicy(
    computer="cluster",
    queue="standard",
    max_concurrent_jobs=1,  # Safe concurrency limit
)

recipe = psteros.SurfaceWorkflowConfig(
    backend="qe",
    calculation=calc_config,
    execution=exec_policy,
    name="oxide_surface_study",
)

# Build the graph (inspectable before submission)
workgraph = psteros.build_surface_workgraph(structures, recipe, submit=False)
```

### 3. Evaluate Surface Thermodynamics

Extract geometry properties directly from the slab structure and compute surface free energies from converged DFT outputs:

```python
# Programmatically extract geometry from the slab
area = slab.surface_area                      # Exposed surface area (Å²)
n_metal = int(slab.composition["Ti"])         # Metal atom count
n_oxygen = int(slab.composition["O"])        # Oxygen atom count

# Compute surface free energy (γ) as a function of oxygen chemical potential (Δμ_O)
point = psteros.surface_energy_oxide_equilibrium(
    slab_energy_ev=slab_energy_from_dft,
    n_metal=n_metal,
    n_oxygen=n_oxygen,
    bulk_formula_energy_ev=bulk_energy_from_dft / n_bulk_formula_units,
    oxygen_reference_energy_ev=e_o2_from_dft,
    delta_mu_oxygen_ev=-1.0,
    surface_area_angstrom2=area,
    surfaces=2,
)

print(f"Surface energy: {point.gamma_j_per_m2:.3f} J/m²")
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


