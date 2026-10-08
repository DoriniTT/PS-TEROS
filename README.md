# PS-TEROS

**PS-TEROS** (**P**redicting **S**tability of **TER**minations **O**f **S**urfaces) is a Python framework for automating *ab initio* surface thermodynamics in [AiiDA](https://www.aiida.net/): oxides, semiconductors, and other compounds, including polar surfaces.

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.apsusc.2025.164350-blue)](https://doi.org/10.1016/j.apsusc.2025.164350)

> **Citation:** If you use PS-TEROS in your research, please cite:  
> T. T. Dorini, M. A. San-Miguel, *Accelerating rational design of oxide surfaces: The PS-TEROS workflow for automated surface stability analysis*, **Applied Surface Science** (2025). [doi:10.1016/j.apsusc.2025.164350](https://doi.org/10.1016/j.apsusc.2025.164350)

### The Problem

In a compound, a single surface orientation rarely exposes just one atomic arrangement; it can exhibit several distinct surface terminations with different stoichiometries (for an oxide: stoichiometric, oxygen-poor, or oxygen-rich cuts; for GaAs: Ga- or As-rich ones). *Ab initio* atomistic thermodynamics determines the relative stability of these terminations by coupling slab models to chemical reservoirs of the constituent species, such as the oxygen reservoir (*Δμ*<sub>O</sub>) of an oxide.

Determining which termination is thermodynamically favored requires coordinating a complex set of interdependent DFT simulations: generating and relaxing multiple slab terminations, calculating matching bulk references, and evaluating gas-phase reference states under strictly identical numerical settings. Managing this multi-structure workflow manually is tedious, error-prone, and hard to reproduce.

### The Solution

PS-TEROS automates the pathway from crystal structures to thermodynamic stability:

- **Slab Builders:** Programmatically generates bulk references and multiple symmetric/asymmetric surface terminations (e.g., rutile SnO₂).
- **Typed DFT Recipes:** Enforces strict parameter harmony across all terminations, bulk references, and reservoirs in **VASP** (the central engine, through aiida-vasp); the same recipes also run **Quantum ESPRESSO**.
- **AiiDA WorkGraphs:** Orchestrates multi-stage workflows (relaxation → static SCF) with bounded job concurrency and full provenance tracking.
- **Polar Surfaces:** Absolute surface energies of polar faces (zinc blende (111), wurtzite (0001)) with a pseudo-hydrogen passivated bottom, on the same scale as non-polar ones ([guide](docs/POLAR_SURFACES.md)).
- **Charge-Neutral Terminations:** Symmetric, charge-neutral slabs for semiconductors and insulators, and `ChargeNeutralSurfaceStudy` to run all of them with their references in one graph ([guide](docs/CHARGE_NEUTRAL_TERMINATIONS.md)).
- **Pure-Python Thermodynamics:** Evaluates surface free energies (J/m², eV/Å²) and termination phase diagrams directly as a function of chemical potential (Δμ<sub>O</sub> for an oxide, Δμ<sub>As</sub> for GaAs, ...).

## Quickstart

The complete, tested version of these three steps is
[`examples/vasp_surface_phase_diagram`](examples/vasp_surface_phase_diagram) (SnO₂(110) with VASP;
[`examples/qe_surface_phase_diagram`](examples/qe_surface_phase_diagram) is the same campaign with Quantum ESPRESSO).

### 1. Prepare the structures (pymatgen or ASE)

Any binary metal oxide (MO, MO₂, M₂O₃, …) works; other binary compounds (GaAs, ZnS, GaN, …) use `BinaryReferences`, and ternary compounds use `TernaryOxideReferences` or `TernaryReferences` with `ternary_surface_phase_diagram` (see the [phase-diagram guide](docs/source/phase-diagram.rst)). Cut symmetric slabs from the **relaxed** bulk, so that both
faces are the same termination; different cuts expose different amounts of oxygen:

```python
import psteros
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator

bulk = Structure.from_file("bulk_oxide_relaxed.cif")  # or Structure.from_ase_atoms(atoms)
metal = Structure.from_file("metal.cif")              # sets the O-poor limit

slabs = SlabGenerator(bulk, (1, 1, 0), min_slab_size=12.0, min_vacuum_size=15.0,
                      center_slab=True).get_slabs(symmetrize=True)
terminations = {f"term_{i}": slab for i, slab in enumerate(slabs)}  # e.g. Sn6O14, Sn8O16, Sn8O14
references = {"bulk": bulk, "metal": metal, "o2": psteros.triplet_o2_cell()}
```

For semiconductors and insulators, `psteros.ChargeNeutralSurfaceStudy` cuts
the charge-neutral terminations instead and also gives the per-calculation
overrides and the analysis ([guide](docs/CHARGE_NEUTRAL_TERMINATIONS.md#in-the-ps-teros-workflow)).

### 2. Relax and compute the energies (AiiDA WorkGraph)

One recipe per stage keeps every calculation on the same numerical settings; per-structure overrides
handle the bulk cell relaxations and the triplet O₂:

```python
INCAR = {"ENCUT": 520, "PREC": "Accurate", "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05, "LREAL": False}

def recipe(ionic, overrides):
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label="vasp@my-cluster",
            incar={**INCAR, **ionic},
            potential_family="PBE",                  # uploaded with `aiida-vasp potcar uploadfamily`
            potential_mapping={"Sn": "Sn_d"},        # other elements use the POTCAR of the same name
            kpoints_spacing=0.25,
            max_iterations=3,                        # allow restarts, e.g. after the walltime
        ),
        execution=psteros.ExecutionPolicy(
            computer="my-cluster", queue="my-queue", max_concurrent_jobs=1,
            resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
        ),
        role_overrides=overrides,
    )

cell = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
triplet = psteros.CalculationOverride(
    parameters={"INCAR": {"ISPIN": 2, "MAGMOM": [1.0, 1.0]}},
    kpoints_distance=10.0,  # Gamma only for the molecule
)
graph = psteros.build_relax_static_workgraph(
    {**references, **terminations},
    recipe({"IBRION": 2, "NSW": 200, "ISIF": 2, "EDIFFG": -0.01}, {"bulk": cell, "metal": cell, "o2": triplet}),
    recipe({"IBRION": -1, "NSW": 0}, {"o2": triplet}),
    submit=True,  # False builds the graph for inspection only
)
```

Each static calculation starts from its relaxed structure, and its energy is the one used below. Slab atoms
can be fixed with `psteros.CalculationOverride(fixed_sites=psteros.central_sites(slab))`.

### 3. Build the surface phase diagram

When the graph has finished, feed the static energies to the pure-Python analysis. It writes the phase
diagram **directly as a figure**, and **exports all the data as CSV** for plotting it your own way:

```python
from aiida import orm

energy, relaxed = psteros.read_vasp_results(orm.load_node(graph.pk), [*references, *terminations])

oxide = psteros.BinaryOxideReferences(
    bulk_energy_ev=energy["bulk"],
    bulk_composition=bulk.composition,
    oxygen_molecule_energy_ev=energy["o2"],
    metal_energy_per_atom_ev=energy["metal"] / len(metal),
)
diagram = psteros.surface_phase_diagram(
    [psteros.SlabTermination.from_structure(label, energy[label], relaxed[label]) for label in terminations],
    oxide,
)
print(diagram.transitions)                  # exact Δμ_O where the stable termination changes
diagram.plot("phase_diagram.png")           # γ(Δμ_O), stability window and stable-termination strip
diagram.to_csv("phase_diagram.csv")         # delta_mu_O_eV, gamma_<termination>_Jm2, ..., stable_termination
```

The same graph runs with Quantum ESPRESSO: `backend="qe"` with `psteros.QeCalculationConfig` recipes and
`psteros.read_qe_results` (install with `pip install '.[qe]'`, see the [QE guide](docs/source/qe-workflow.rst)).

## Installation

```bash
pip install .
verdi daemon restart --reset
```

This installs aiida-vasp with PS-TEROS; `pip install '.[qe]'` adds Quantum ESPRESSO support. Install PS-TEROS
in the same Python environment as the AiiDA daemon. For configuring AiiDA computers, the VASP code and the
POTCAR family, see the [Installation Guide](docs/source/installation.rst).

## Documentation & Tutorials

- **[First Tutorial](docs/source/tutorial.rst):** Build your first unsubmitted AiiDA WorkGraph.
- **[Core Concepts](docs/source/concepts.rst):** How structures, calculation recipes, execution policies, and provenance connect.
- **[SnO₂ Surface Model](docs/source/examples.rst):** Deep dive into terminations and thermodynamic reference states.
- **[VASP Guide](docs/source/vasp-workflow.rst):** Setting up a two-stage relaxation → static workflow, with per-structure settings and fixed atoms.
- **[Quantum ESPRESSO Guide](docs/source/qe-workflow.rst):** The same workflow with Quantum ESPRESSO.
- **[Surface Phase Diagrams](docs/source/phase-diagram.rst):** From energies to γ(Δμ<sub>O</sub>) (binary oxides) or stable-termination maps over (Δμ<sub>A</sub>, Δμ<sub>O</sub>) (ternary oxides), as a figure or a CSV table.
- **[API Reference](docs/source/api.rst):** Public classes, functions, and configuration schemas.


