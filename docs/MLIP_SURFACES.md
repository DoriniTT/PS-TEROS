# Surface Phase Diagrams with a Fine-Tuned MACE Model

## Summary

Relaxing every slab termination with DFT is what limits how many
terminations a surface study can afford. The **MLIP route** of PS-TEROS
replaces those relaxations: a few DFT single-point calculations ("statics")
teach a MACE foundation model your system, the model relaxes every slab on a
GPU, and one more DFT static at each relaxed geometry gives the energies of the
surface phase diagram.

The idea in one line: **MACE finds the relaxed geometry, DFT gives its
energy.** The model's own energies are not used.

It is the route PS-TEROS recommends when you have a GPU. Without one, relax the
slabs with DFT as before (`build_relax_static_workgraph` and the QE and CP2K
builders); nothing in that route changes.

The route works for metals, semiconductors and insulators. The difference lies
in the slabs: for a semiconductor or an insulator they must have the right
electron count (charge-neutral, or passivated), so PS-TEROS cuts them for you
(see [Which study](#which-study)).

![The MLIP route of PS-TEROS](images/mlip_route.png)

## The recipe

**You need:** a DFT code (QE, VASP or CP2K, through PS-TEROS recipes), the
MACE-MH-1 foundation model and its replay data on a GPU computer, and an AiiDA
code for the Python of a MACE environment there.

1. **Relax the references with DFT, once.** The bulk, the elemental
   references (a solid, or a molecule such as O₂ in a box) and, for a ternary
   compound, the competing phases, with the recipe of the study. These are
   small and cheap. Cut the slabs at the relaxed bulk cell and **do not relax
   them**.
2. **One DFT static of every structure**, bulk, references and slabs as
   built: energies and forces, the training labels.
3. **Check the electron count.** Slabs cut by PS-TEROS
   (`ChargeNeutralSurfaceStudy`, `PolarSurfaceStudy`) pass. Your own slabs
   need a metallic bulk (partly filled states at the Fermi level): in a metal
   any cut will do, in a semiconductor it will not.
4. **Replay fine-tune of MACE-MH-1** (GPU, about 30 min). A new head learns
   your DFT labels while the original head keeps training on ~150 small
   structures of your elements labelled by MH-1 itself, so the model does not
   forget general chemistry.
5. **MACE relaxes every structure** (GPU, minutes), with the cell fixed and
   the slab centre held in place.
6. **One active-learning round.** One DFT static at each MACE-relaxed
   geometry joins the labels; the model is fine-tuned again from MH-1 on all of
   them and relaxes the slabs again from the start. The first model has only
   seen unrelaxed slabs; these statics add labels near the relaxed geometries.
7. **One final DFT static at each relaxed geometry.** These energies, for the
   slabs and the references alike, go into the surface phase diagram.
8. *Optional:* **vibrations on the GPU.** The final model gives the harmonic
   modes of every structure at its relaxed geometry, and the diagram uses free
   energies at a temperature of your choice.

DFT per termination: 3 statics (steps 2, 6 and 7) instead of a DFT
relaxation; the vibrations add no DFT.

## Which study

| Your slabs | Study | Metal check |
|---|---|---|
| Charge-neutral symmetric slabs of an element, a binary or a ternary compound, cut by PS-TEROS (one or more orientations) | `MlipRoute(ChargeNeutralSurfaceStudy(...))` | none: the slabs are charge-neutral by construction |
| Polar faces (zinc blende (111), wurtzite (0001), ...) on a pseudo-hydrogen passivated bottom, VASP | `MlipRoute(PolarSurfaceStudy(...))` | none: the bottom is passivated by construction |
| Your own slabs of a binary compound (e.g. `rutile_110_slab`) | `MlipSurfaceStudy(bulk, references, slabs)` | the bulk must be a metal, unless `build(..., allow_non_metallic=True)` |
| Any labelled structures | `build_mlip_surface_workgraph(...)` | `metallic_bulk` must be a metal, unless `allow_non_metallic=True` |

With `allow_non_metallic=True` the graph warns that the electron count of the
slabs is your responsibility. See
[`CHARGE_NEUTRAL_TERMINATIONS.md`](CHARGE_NEUTRAL_TERMINATIONS.md) and
[`POLAR_SURFACES.md`](POLAR_SURFACES.md) for how those slabs are cut.

## Quick start: your own slabs (`MlipSurfaceStudy`)

```python
import psteros
from aiida import orm

# step 1: references relaxed with the recipe below (build_cp2k_relax_static_workgraph), slabs cut at that cell
slabs = {t: psteros.rutile_110_slab(metal="Ru", termination=t, a=a, c=c, u=u, vacuum_angstrom=22.0)[0]
         for t in psteros.RUTILE_110_TERMINATIONS}
study = psteros.MlipSurfaceStudy(bulk, {"Ru": ru_metal, "O": o2_box}, slabs)
print(study)            # 27 DFT statics and 4 GPU jobs for the six RuO2(110) terminations

static = psteros.SurfaceWorkflowConfig(backend="cp2k", name="ruo2", calculation=cp2k_recipe,
                                       execution=cpu_execution, role_overrides={"ref_O": triplet_o2})
fine_tune = psteros.MaceFineTuneConfig("/work/me/models/mace-mh-1.model",
                                       replay_file="replay-data-mh-1-omat-pbe.xyz")
job = psteros.MlipJobConfig("mace-python@gpu-cluster", psteros.ExecutionPolicy(
    queue="gpu", with_mpi=False, max_wallclock_seconds=6 * 3600,
    resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1}))

graph = study.build(static, fine_tune, job=job, submit=True)              # steps 2-7

# when it has finished:
result = study.analyse(orm.load_node(graph.pk))                           # energies of step 7 only
print(result)                                                             # stable termination per Δμ_O range
result.diagram.plot("ruo2_110_mlip.png")
result.diagram.to_csv("ruo2_110_mlip.csv")
```

`MlipSurfaceStudy` labels the calculations `bulk`, `ref_<element>` and the
slab labels, fixes the sites within 1.5 Å of each slab mid-plane during the
MACE relaxations (`psteros.central_sites`), builds the graph with
`psteros.build_mlip_surface_workgraph` and turns its final DFT energies into a
`psteros.SurfacePhaseDiagram`. A complete script is
[`examples/mlip_ruo2/ruo2_mlip_graph.py`](../examples/mlip_ruo2/ruo2_mlip_graph.py).

An alternative to the active-learning round: add rattled copies of each slab
to the first training set (`MlipSurfaceStudy(..., rattled=5)`) and skip the
round (`study.build(..., active_learning_rounds=0)`). One fine-tune then runs,
with more DFT statics up front.

## Quick start: slabs cut by PS-TEROS (`MlipRoute`)

A semiconductor, with the charge-neutral slabs of two orientations and VASP
statics:

```python
study = psteros.ChargeNeutralSurfaceStudy(sno2_bulk, [(1, 1, 0), (1, 0, 0)],
                                          references={"Sn": tin, "O": psteros.triplet_o2_cell()})
route = psteros.MlipRoute(study)        # slab centres held in place (fixed_half_width_angstrom=1.5)
print(route)

static = psteros.SurfaceWorkflowConfig(
    backend="vasp", name="sno2", execution=cpu_execution,
    calculation=psteros.VaspCalculationConfig("vasp@cpu-cluster", incar={
        "ENCUT": 520, "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05, "IBRION": -1, "NSW": 0, "LREAL": False}),
    role_overrides=study.vasp_overrides("static"))             # triplet O2, Gamma-only molecule
graph = route.build(static, fine_tune, job=job, submit=True)

result = route.analyse(orm.load_node(graph.pk))                # the study's own analysis
result.diagram.plot("sno2_mlip.png")
```

`route.read(graph)` gives `(energies, relaxed_structures)`, the same pair as
`psteros.read_vasp_results`, for any analysis of your own. A ternary study
(with `competing_phases`) gives a `TernarySurfacePhaseDiagram`; an element its
surface energies.

**Polar faces** run the same way with VASP:

```python
study = psteros.PolarSurfaceStudy(gaas, faces=[(1, 1, 1), (-1, -1, -1)], references={"Ga": ga, "As": arsenic})
route = psteros.MlipRoute(study)
static = psteros.SurfaceWorkflowConfig(
    backend="vasp", name="gaas", execution=cpu_execution,
    calculation=psteros.VaspCalculationConfig("vasp@cpu-cluster", incar=static_incar,
                                              potential_mapping=study.potential_mapping({"Ga": "Ga_d"})),
    role_overrides=study.vasp_overrides("static"))             # dipole correction of the asymmetric slabs
graph = route.build(static, fine_tune, job=job, submit=True)
result = route.analyse(orm.load_node(graph.pk))                # PolarStudyResult with the bottom checks
```

Every atom of a polar slab relaxes, as in the polar method; the bottom check
then compares the MACE-relaxed bottoms. For MACE the pseudo-hydrogens are
hydrogen atoms; the fine-tune learns them from their DFT labels.

Complete scripts: [`examples/mlip_routes/`](../examples/mlip_routes/).

## Vibrations on the GPU

```python
graph = route.build(static, fine_tune, job=job, submit=True,
                    vibrations=psteros.VibrationsConfig(supercells={"bulk": (3, 3, 3), "ref_Sn": (3, 3, 3)}))
result = route.analyse(orm.load_node(graph.pk), temperature_k=800.0)
```

One more GPU job: the deployed model gives the Γ-point harmonic modes of every
structure at its relaxed geometry, with the finite differences of
`build_vibrations_workgraph` (see `docs/source/vibrations.rst`). The sites held
in place are not displaced and their region is counted as bulk, so it must have
the bulk stoichiometry (`build` checks it); a bulk is displaced in its
`supercells` entry; the gas references are molecules. `analyse(...,
temperature_k=T)` gives solids `E + F_vib(T)` and gas references `E + ZPE`
(`psteros.free_energies_with_vibrations`), with `E` still the final DFT
statics. `MlipSurfaceStudy` takes the same options.

The model's modes have not yet been compared with DFT for a series of
systems, so check a few: `route.dft_vibrations(graph, static)` builds the DFT
finite-difference graph of `build_vibrations_workgraph` (QE, VASP, and here also
CP2K) on the MACE geometries with the same fixed sites, and
`route.analyse(graph, temperature_k=T, vibrations=psteros.read_vibrations(check, labels))`
uses its modes instead. Polar slabs have no vibrations yet.

## Other MACE models

MACE-MH-1 is the default, but any MACE model file works.
`psteros.MACE_FOUNDATION_MODELS` lists the released ones that have been run
through the route, with their heads, the DFT they were trained on, licence,
download and replay data. Choose by the level of theory of your DFT:

| Your DFT | Model (head) |
|---|---|
| PBE, with the Materials Project +U on TM oxides/fluorides, or no such compound | `mace-mh-1` (`omat_pbe`), `mace-omat-0-medium` |
| PBE without +U on TM oxides/fluorides | `mace-matpes-pbe-0` |
| r2SCAN | `mace-mh-1` (`matpes_r2scan`), `mace-matpes-r2scan-0` |
| an MIT-licensed model is required | `mace-mpa-0-medium` |

```python
model = psteros.mace_foundation_model("mace-matpes-pbe-0")   # model.url, model.replay_url
fine_tune = psteros.MaceFineTuneConfig.for_foundation(
    "mace-matpes-pbe-0", "/work/me/models/" + model.filename, replay_file="matpes-pbe-replay-data.xyz")
```

With a named model the builders warn when the DFT recipe's functional or +U
differs from the head's training data: the per-element energy shifts of the
fine-tune absorb constant offsets, not a different functional. Details, and
the head names of each file, in [`docs/source/mlip.rst`](source/mlip.rst).

## Rules the route enforces

* **The right electron count.** Hand-cut slabs need a metallic bulk (the first
  fine-tune stops otherwise, `psteros.require_metallic`) unless
  `allow_non_metallic=True`, which warns; the slabs of `MlipRoute` are
  charge-neutral or passivated by construction.
* **One energy source per diagram.** A constant offset per element cancels in
  γ only if the slabs and every reference carry it, so `analyse` takes every
  energy from the final DFT statics (or, explicitly, all from the model:
  `energy_source="mlip"`), never a mix (`psteros.check_one_energy_source`).
  The vibrations do not break this: a constant offset has no curvature.
* **Every energy with its vibrations.** With `temperature_k`, every structure
  needs its modes; a frozen slab centre that is not whole bulk formula units is
  refused before anything runs.
* **No unrelaxed geometry reaches DFT.** A MACE relaxation that does not
  reach 0.01 eV/Å stops the graph instead of sending its geometry to the next
  DFT static.
* **The models stay.** Each relaxation and vibrations job links the remote
  folder of its fine-tune, where the models are; that folder is never cleaned.

## Limits

* Converge the DFT recipe (cutoff, basis, k-points, vacuum) for your system
  first, and check a new system against a few DFT relaxations before relying
  on the route.
* Polar studies need VASP (the pseudo-hydrogen POTCARs) and have no
  vibrations yet.
* `MlipSurfaceStudy` takes binary compounds; `MlipRoute` with a
  `ChargeNeutralSurfaceStudy` takes elements and ternary compounds too, and
  `build_mlip_surface_workgraph` any labelled structures.

## Reference

* `psteros.MlipRoute`, `psteros.MlipSurfaceStudy`,
  `psteros.build_mlip_surface_workgraph`, `psteros.read_mlip_energies`,
  `psteros.read_mlip_results`, `psteros.MaceFineTuneConfig`,
  `psteros.MlipJobConfig`, `psteros.free_energies_with_vibrations`,
  `psteros.prepare_vibrations_job`, `psteros.MACE_FOUNDATION_MODELS`,
  `psteros.foundation_level_warnings`: see
  [`docs/source/mlip.rst`](source/mlip.rst) and [`docs/source/api.rst`](source/api.rst).
* The flowchart: `python docs/images/mlip_route.py`.
