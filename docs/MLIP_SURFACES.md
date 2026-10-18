# Surface Phase Diagrams with a Fine-Tuned MACE Model (metallic systems)

## Summary

Relaxing every slab termination with DFT is what limits how many
terminations a surface study can afford. The **MLIP route** of PS-TEROS
replaces those relaxations: a few DFT single-point calculations ("statics")
teach a MACE foundation model your system, the model relaxes every slab on a
GPU, and one more DFT static at each relaxed geometry gives the energies of the
surface phase diagram.

The idea in one line: **MACE finds the relaxed geometry, DFT gives its
energy.** The model's own energies are not used.

The route is for **metallic systems only** (metals, alloys, metallic oxides
such as RuO₂, IrO₂). PS-TEROS checks that the bulk has partly filled states at
the Fermi level and stops otherwise. It was validated on RuO₂(110) against full
CP2K relaxations and against Reuter and Scheffler, *Phys. Rev. B* **65**,
035406 (2002), with pass/fail criteria fixed before any result (see
[`validation/ruo2_mace_lovelace/RESULTS.md`](../validation/ruo2_mace_lovelace/RESULTS.md)).

![The MLIP route of PS-TEROS](images/mlip_route.png)

## The recipe

**You need:** a DFT code (QE, VASP or CP2K, through PS-TEROS recipes), the
MACE-MH-1 foundation model and its replay data on a GPU computer, and an AiiDA
code for the Python of a MACE environment there.

1. **Relax the references with DFT, once.** The bulk (e.g. RuO₂), the metal
   (e.g. hcp Ru) and the gas molecule (O₂ in a box), with the recipe of the
   study. These are small and cheap. Cut the slabs at the relaxed bulk cell
   and **do not relax them**.
2. **One DFT static of every structure**, bulk, references and slabs as
   built: energies and forces, the training labels.
3. **Check that the bulk is a metal.** If no state is partly filled at the
   Fermi level, the graph stops.
4. **Replay fine-tune of MACE-MH-1** (GPU, about 30 min). A new head learns
   your DFT labels while the original head keeps training on ~150 small
   structures of your elements labelled by MH-1 itself, so the model does not
   forget general chemistry.
5. **MACE relaxes every structure** (GPU, minutes), with the cell fixed and
   the slab centre held in place.
6. **One active-learning round.** One DFT static at each MACE-relaxed
   geometry joins the labels; the model is fine-tuned again from MH-1 on all of
   them and relaxes the slabs again from the start. The first model has only
   seen unrelaxed slabs, so on some terminations it stops short of the DFT
   minimum; these statics show it the region near the minimum.
7. **One final DFT static at each relaxed geometry.** These energies, for the
   slabs and the references alike, go into the surface phase diagram.

DFT per termination: 3 statics (steps 2, 6 and 7) instead of a DFT
relaxation.

## Quick start

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

## Validation on RuO₂(110)

![RuO2(110): CP2K relaxations against the MLIP route](images/mlip_ruo2_validation.png)

*(a) Surface energies of the six terminations from full CP2K relaxations
(solid) and from the MLIP route (dashed): every termination within
0.5 meV/Å², phase boundaries within 0.016 eV, the narrow O^br + ½ O^cus window
(0.09 eV wide) reproduced; the bands on top show the stable termination of
each. (b) The MLIP route (solid) against the four terminations of Reuter and
Scheffler's Fig. 3 (dotted; FP-LAPW, GGA), read from the vector graphics of
[arXiv:cond-mat/0107229](https://arxiv.org/abs/cond-mat/0107229): the same
picture, O^br stable at low and O^br + O^cus at high oxygen chemical potential,
the O-poor no-O^br termination never stable, the boundary at −0.94 eV for them
and −0.91 to −0.85 eV for the MLIP route, which also resolves a narrow
half-covered window; its surface energies lie 8–10 meV/Å² lower, within
their stated accuracy of 10 meV/Å². (c) Largest error of each termination
with the model's own energies (orange) and with a CP2K static at the MACE
geometry (blue). (d) Share of the CP2K relaxation energy that the MACE
relaxation recovers, before and after the active-learning round.*

| What was compared | Result |
| --- | --- |
| CP2K reference against Reuter and Scheffler | O^br → O^cus transition, Δγ and ΔH_f within the tolerances; the metal-terminated surface never stable |
| MLIP route against Reuter and Scheffler (their Fig. 3) | same stable terminations; boundary −0.91 to −0.85 eV against −0.94 eV; γ 8–10 meV/Å² lower (their accuracy: 10 meV/Å²) |
| MLIP route against CP2K relaxations | γ within 0.49 meV/Å²; transitions −0.908 and −0.848 eV against −0.924 and −0.835 eV; 96–99 % of the relaxation energy |
| Model's own energies against CP2K | 6–9 meV/Å² off for the three seeds: enough to close or move the 0.09 eV window |
| DFT time of the slabs | 34 % of relaxing them with CP2K (these cells relax in 7–13 DFT steps; the saving grows with longer relaxations) |

What did **not** work, and why it matters:

* **The model's own energies** (6–9 meV/Å² off): use the final DFT statics.
* **Fine-tuning without replay**: the model forgot what it knew and moved
  the slabs away from the DFT minima, worse than the foundation model alone.
* **Training five times longer instead of the active-learning round**: the
  formation enthalpy became exact, but the relaxations of two terminations
  still stopped at 72–77 % of the DFT relaxation. The fix was data near the
  relaxed geometries, not more training.

An alternative that also passed: add five rattled copies of each slab to the
first training set (`MlipSurfaceStudy(..., rattled=5)`) and skip the
active-learning round (`study.build(..., active_learning_rounds=0)`). One
fine-tune is then enough, but it costs more DFT (80 % of the DFT
relaxations, against 34 %).

## Rules the route enforces

* **Metals only.** The first fine-tune stops when the bulk static has no
  partly filled states (`psteros.require_metallic`).
* **One energy source per diagram.** A constant offset per element cancels in
  γ only if the slabs and every reference carry it, so `analyse` takes every
  energy from the final DFT statics (or, explicitly, all from the model:
  `energy_source="mlip"`), never a mix (`psteros.check_one_energy_source`).
* **No unrelaxed geometry reaches DFT.** A MACE relaxation that does not
  reach 0.01 eV/Å stops the graph instead of sending its geometry to the next
  DFT static.
* **The models stay.** Each relaxation job links the remote folder of its
  fine-tune, where the models are; that folder is never cleaned.

## Limits

* Validated on one metallic oxide surface (RuO₂(110), 1×1 and 1×2 cells,
  16–38 atoms). Converge the DFT recipe (cutoff, basis, k-points, vacuum) for
  your system first, and check a new system against a few DFT relaxations.
* Semiconductors and insulators need a separate treatment and are refused.
* Binary compounds in `MlipSurfaceStudy`; the lower-level
  `build_mlip_surface_workgraph` takes any labelled structures.

## Reference

* `psteros.MlipSurfaceStudy`, `psteros.build_mlip_surface_workgraph`,
  `psteros.read_mlip_energies`, `psteros.read_mlip_results`,
  `psteros.MaceFineTuneConfig`, `psteros.MlipJobConfig`: see
  [`docs/source/mlip.rst`](source/mlip.rst) and [`docs/source/api.rst`](source/api.rst).
* The plan and its pre-registered gates:
  [`docs/plans/MACE_SURFACE_PLAN.md`](plans/MACE_SURFACE_PLAN.md).
* The figures: `python docs/images/mlip_route.py` and
  `python validation/ruo2_mace_lovelace/c9_figure.py plot` (the Reuter and
  Scheffler lines are in that script, `REUTER`).
* K. Reuter and M. Scheffler, *Composition, structure, and stability of
  RuO₂(110) as a function of oxygen pressure*, Phys. Rev. B **65**, 035406
  (2002), [arXiv:cond-mat/0107229](https://arxiv.org/abs/cond-mat/0107229);
  erratum Phys. Rev. B **75**, 049901 (2007), on the temperature dependence of
  the formation energy (the surface energies of Fig. 3 are unchanged).
