# Plan: surface phase diagrams from a fine-tuned MACE model, validated on RuO2(110)

Branch: `features/mace-surface-phase-diagram` (worktree from `origin/main`, 26c72e9). Status: `[ ]` to do,
`[~]` in progress, `[x]` done.

**Status (2026-10-17): approved.** The user took every recommended option of D1–D6 ("Use your recommendations").
The gates below are fixed before any result. Set by the user on 2026-10-17:

- **Metallic systems only**; no Ag3PO4 application; semiconductors will be planned separately ("Semiconductors we
  need to be more careful, which I will plan for the future").
- **DFT with CP2K**, not VASP ("I want to use CP2K, not vasp"). PS-TEROS has no CP2K backend yet, so this plan
  adds one (Part A).

## Idea

Relaxing every termination with DFT is what limits how many terminations we can study. Instead:

1. CP2K **static** calculations only: bulk, elemental references, the O2 molecule and every termination as built.
2. Fine-tune MACE-MH-1 on those labels with **multi-head replay** (keeps the foundation model's general knowledge,
   which a relaxation needs once it leaves the training geometries).
3. Relax every termination with the fine-tuned model on a Lovelace A100.
4. Build the surface phase diagram from the model's energies (and, as a cheap check, from one CP2K static on each
   MACE-relaxed geometry).

The question this plan answers: **does step 4 give the same diagram as the full CP2K relax -> static route?**
RuO2(110) is the test because it is a metallic oxide with a well-known DFT diagram (Reuter and Scheffler) and because
full relaxations of all its terminations are cheap enough to serve as ground truth.

## Scope: metallic systems only

The feature is for surfaces of materials whose bulk is a metal in the DFT used for the labels: metals, alloys
and metallic oxides such as RuO2 and IrO2. Semiconductors and insulators are out of scope and get their own plan
later.

- **Enforced, not only documented.** The bulk static must be metallic: after the smeared SCF, some states are
  filled between 1% and 99% (partly filled states at the Fermi level). Otherwise `require_metallic` raises
  `ValueError` naming the bulk label. One rule for every code: VASP and QE give occupations per k-point and band;
  CP2K 2024.3 (the Lovelace module) cannot print k-point eigenvalues (`MO_KP` arrived later), so for CP2K the rule
  is applied to the occupations per energy bin of its k-point DOS file (`Cp2kDos.partially_occupied_bins`). Gas
  references (O2) are not checked. Slabs are not gated, but each one's result is reported.
- **No override flag.** The semiconductor plan can add one later as a new opt-in option.

## Why it can work, and where it can fail

- **Element offsets cancel.** With γ = [E_slab − N_Ru E_RuO2 − (N_O − 2 N_Ru) μ_O] / 2A and μ_O = ½ E_O2 + Δμ_O, adding
  a constant δ_Ru per Ru and δ_O per O to every energy changes γ by
  N_Ru δ_Ru + N_O δ_O − N_Ru(δ_Ru + 2δ_O) − (N_O − 2N_Ru) δ_O = 0. So the E0s of the fine-tuned head do not matter,
  **as long as every energy in one diagram comes from the same model**. The Ag3PO4 project found that E0s fitted
  at one composition are not transferable elemental references; mixing a model's slab energy with a DFT bulk
  energy would bring that error straight into γ. Rule: one energy source per diagram, and mixing raises `ValueError`.
- **Every composition is in the training set.** Each termination has its own static label, so the model is never
  asked to predict a stoichiometry it has not seen; only the relaxation moves it away from the labels.
- **Static-only labels.** One geometry per termination gives one force label per atom. The relaxation explores
  geometries the model has not seen, which is the job of the replay head. A rattled variant (D2) tests whether a
  few extra statics per termination are needed.
- **What we know about replay so far.** The 2026-09-21 pilot on Ag3PO4 found replay worse on forces (56 against 39
  meV/Å), but the comparison was confounded (different learning rate, EMA and scheduler), and forgetting was never
  measured. Goal 07 has no replay result yet. RuO2 gives a clean comparison: same recipe, replay on or off.
- **CP2K labels and an `omat_pbe` base.** MH-1's `omat_pbe` head was trained on VASP PBE (PAW). CP2K PBE with GTH
  pseudopotentials differs from it by large per-element constants (about −253,660 eV per cell for the Ag3PO4 slab)
  and by small chemical differences. The constants are handled by the E0s of the fine-tuned head (next point);
  the multi-head set-up keeps the two levels of theory in separate heads, which is what it is made for.
- **E0s from several compositions.** δ_Ru and δ_O are fitted by least squares to E_CP2K − E_MH1 over all training
  frames and added to the base E0s of the `Default` head. Unlike Ag3PO4 (one composition), hcp Ru and O2 pin each
  element separately. By the first point, this choice cannot change γ; it only makes the fit start close.
- **Energy consistent with forces.** With Fermi-Dirac smearing, the energy whose derivatives are the forces is the
  free energy, including the electronic entropy term. Labels use that energy. B0 checks it with a finite-difference
  force on one atom, and every label reports its entropy term.
- **Gaussian basis sets.** Unlike plane waves, a MOLOPT basis has a superposition error that differs between the
  O2 molecule, the bulk and the slabs; ½ E_O2 moves the whole O-rich side of the diagram. The basis is therefore
  fixed by a convergence test on ΔH_f and on the differences in γ (B0), not taken from defaults.

## Literature reference (Reuter and Scheffler, PRB 65, 035406, 2002)

FP-LAPW, GGA (PBE), symmetric slabs of 3 O-(RuO)-O trilayers, outermost trilayer relaxed, ≈ 11 Å vacuum,
(5×10×1) k-mesh for the (1×1) cell, bulk computed in the same (110)-oriented cell as the slabs, numerical accuracy of
γ within 10 meV/Å².

- Terminations: (1×1) O-bridge (stoichiometric), O-cus (extra O atop Ru-cus, O-rich), Ru (O-bridge removed, O-poor).
- O-cus is the most stable for **Δμ_O > −0.93 eV** (LDA: −0.85 eV); O-bridge below that; Ru is never stable.
- In the O-rich limit, γ(O-cus) is **49 meV/Å²** below γ(O-bridge).
- ΔG_f(RuO2) = **−3.16 eV** per formula unit (experiment −3.19 eV), so the O-poor limit is ½ ΔG_f ≈ −1.58 eV.
- A (2×1) half-covered O-cus crosses O-bridge at the same Δμ_O: the lateral interaction between O-cus atoms is
  negligible.

Follow-ups for later systems: O2/CO (PRB 68, 045407, 2003), humid conditions with OH/H2O (Sun, Reuter and
Scheffler, PRB 67, 205424, 2003).

## Fixed choices

| Item | Value |
|---|---|
| Base model | MACE-MH-1, head `omat_pbe`, SHA-256 `a522eb7f…be47` (same file as Ag3PO4 r4; Lovelace copy `/work/dorinitt/goal03_inputs/models/mace-mh-1.model`) |
| Software on the GPU | mace-torch 0.3.16, torch 2.5.1+cu118, ASE 3.29 (the Ag3PO4 venv, D4) |
| Recipe | Ag3PO4 r4/goal-07: float64, batch 1, universal loss, E/F weights 10/10, Huber 0.01, constant lr 1e-3 (`force_mh_ft_lr: true`), EMA 0.999, clip 1.0, the final EMA model after N passes (N = 100, validation reported, never used to pick) |
| Replay | `multiheads_finetuning: true`; every structure of `replay-data-mh-1-omat-pbe.xyz` (SHA-256 `8140eea1…`) of 2–128 atoms made only of Ru and/or O; one in five as replay validation; labels are MH-1 `omat_pbe` float64 predictions computed in the job; both heads must reproduce MH-1 (1e-8 eV/atom, 1e-6 eV/Å) before training |
| Heads | fine-tuned head `Default`, replay head `pt_head`; every evaluation and relaxation names its head; validation heads ordered with `Default` last (mace 0.3.16 picks checkpoints by the last head) |
| Seeds | 3 per model: 20260920, 20260921, 20260922. The gates G1–G4 use the deployment seed 20260920 (as for Ag3PO4 r4/r5; fixed 2026-10-17, before any result); the other two give the committee spread |
| Slabs | symmetric, 3 trilayers, ≥ 15 Å vacuum, on the CP2K-relaxed bulk lattice; the central trilayer fixed (`central_sites`) in both the CP2K and the MACE relaxation |
| CP2K method | GPW, PBE without dispersion (as the paper), GTH-PBE pseudopotentials (Ru q16, O q6), MOLOPT basis and cutoff from B0 (D1); diagonalisation with Broyden mixing, Fermi-Dirac smearing and `ADDED_MOS` (no OT: RuO2 is a metal); Monkhorst-Pack k-points, the same in-plane density for every slab and for the (110)-oriented bulk |
| CP2K references | bulk RuO2 and hcp Ru by `CELL_OPT` (at a cutoff where the stress is converged, B0), then static; bulk RuO2 also in the (110)-oriented (1×1) slab cell, used as E_bulk; O2 triplet (UKS, multiplicity 3, Γ only, ≥ 15 Å box) |
| CP2K outputs | energy in eV from Hartree; forces printed by `FORCE_EVAL/PRINT/FORCES` to a file that is retrieved and parsed (as in the Ag3PO4 `cp2k_labels.py`); eigenvalues and Fermi level of the bulk for the metallic guard |
| Relaxations | the same force threshold in CP2K and MACE: 0.01 eV/Å (`GEO_OPT` `MAX_FORCE` 1.94e-4 Ha/bohr); MACE with ASE FIRE in float64, max 2000 steps; non-converged runs are reported and not used |

### Terminations of RuO2(110)

| Label | Cell | Excess O per surface | Note |
|---|---|---:|---|
| `obr` | 1×1 | 0 | stoichiometric O-bridge |
| `ocus` | 1×1 | +1 | O-rich, stable for Δμ_O > −0.93 eV in the literature |
| `ru` | 1×1 | −1 | O-poor, never stable in the literature |
| `ocus_half` | 2×1 | +½ | every second O-cus; tests lateral interaction |
| `obr_half` | 2×1 | −½ | every second O-bridge removed |
| `ocus_no_obr` | 1×1 | 0 | O-cus present, O-bridge removed; non-bulk-like, expected to relax strongly |

The last three go beyond the paper: they make the set closer to the "many terminations" case and include ones
where relaxation matters.

## What is compared

Models (identical settings except the named difference):

- **M0**: MH-1 zero-shot (`omat_pbe`), no training. Its energies are VASP-PBE-like, so against CP2K R1 it also
  carries the PAW/GTH difference; reported, not gated.
- **M1**: ordinary fine-tune on the statics.
- **M2**: replay fine-tune on the statics. **This is the proposed method.**
- **M3**: replay fine-tune on the statics plus rattled statics (only if D2 includes them).
- **M2R** (2026-10-18, the plan's one active-learning round, triggered because M2 failed routes (a) and (b)): M2's
  frames plus the six route-(b) CP2K statics at the geometries of M2 seed 20260920, M2's recipe and three seeds;
  judged with the same gates.
- **M2L** (added 2026-10-18 at the user's request, after the first M2 model was still converging at pass 100; fixed
  before any of its results; **reported, not gated**): the M2 recipe with 500 passes and a learning rate decaying
  from 1e-3 to 1e-5 (ExponentialLR, gamma 0.990832 per pass), deployment seed only; relaxed and checked by route (b)
  like the others.

Two energy routes per model:

- **(a) model only**: E_M at the model's geometry, for slabs, bulk RuO2, hcp Ru and O2.
- **(b) model geometry, CP2K energy**: one CP2K static on each model-relaxed slab; bulk, Ru and O2 from CP2K.

References:

- **R0 unrelaxed CP2K**: the training statics themselves. The floor: what you get with no relaxation at all.
- **R1 relaxed CP2K**: relax -> static with the new CP2K builder (A3). The ground truth.
- **L**: Reuter and Scheffler 2002.

## Gates (fixed before any result)

- **G0, CP2K against the literature** (checks our DFT set-up, judged before any model): O-bridge -> O-cus transition
  at −0.93 ± 0.15 eV; γ(obr) − γ(ocus) in the O-rich limit 49 ± 10 meV/Å²; ΔH_f(RuO2) −3.16 ± 0.20 eV/f.u.; `ru`
  never stable in the allowed range. A failure stops the plan for a look at the CP2K settings (basis first), not the
  model.
- **G1, same diagram**: at every Δμ_O in [½ ΔH_f, 0], the same stable termination as R1, except within 0.05 eV of
  an R1 transition; every R1 transition reproduced within ±0.05 eV.
- **G2, surface energies**: |γ_M − γ_R1| ≤ 5 meV/Å² for every termination over the whole range (route a);
  ≤ 2 meV/Å² (route b).
- **G3, relaxation recovered**: for every termination whose CP2K relaxation lowers γ by more than 2 meV/Å², the CP2K
  energy at the model geometry recovers at least 90% of it:
  E_CP2K(x0) − E_CP2K(x_M) ≥ 0.9 [E_CP2K(x0) − E_CP2K(x_R1)].
- **G4, model thermochemistry**: |ΔH_f(M) − ΔH_f(CP2K)| ≤ 0.05 eV/f.u.
- **G5, is the test discriminating?** If R0 already passes G1 and G2 against R1, relaxation does not matter for
  RuO2 and the test cannot tell a good model from a bad one. We then say so and add the humid RuO2(110)
  terminations with OH and H2O (Sun et al. 2003) before calling the method validated.

**Outcome rules** (for M2):

- Passes G1–G4 in route (a): the method works as proposed.
- Passes only in route (b): recommend "MACE relaxation + one CP2K static" (still one static instead of a relaxation
  per termination).
- Fails both: **one** active-learning round is allowed, decided now: add the route-(b) CP2K statics to the training
  set, refit M2, relax again, judge again. If it still fails, stop and report.

**Reported, not gated**: M0, M1 against M2 (does replay help, and is zero-shot already enough?); M3 against M2
(do rattled statics help?); the committee spread; the largest displacement of a surface atom against R1 and the
surface Ru–O bond lengths; MACE steps and time; CP2K core-hours of the MACE route against R1; how far each model moves
from MH-1 on the replay validation structures (forgetting); the electronic entropy term of every label.

## Progress

**2026-10-17, code (no cluster job yet):**

- Part A (CP2K backend) and Part B (MLIP layer, job scripts, rutile slabs) are written; `python -m pytest tests/unit`:
  235 passed (207 before; no existing test edited).
- The two job scripts ran end to end on CPU in a local test venv (`~/.cache/psteros_mace_ruo2/venv-cpu`: mace-torch
  0.3.16, torch 2.6.0 CPU, ASE 3.29), on RuO2 bulk, hcp Ru, O2 and two slabs with fake labels (MH-1 + known
  per-element offsets + noise), 2 epochs, 16 replay structures: the E0 shifts came back (O +0.997, Ru −4.994 eV
  for +1 and −5), the untrained heads reproduced MH-1 to 1e-15 eV/atom, `std` read from MH-1 is r4's 0.77358,
  the exported model is the last epoch's; both models relaxed both slabs below 0.01 eV/Å with the centre fixed.
  This checks the plumbing, not the accuracy.
- Ru/O replay set from `replay-data-mh-1-omat-pbe.xyz` (SHA-256 `8140eea1…3df5`, 372,355 frames): 153 distinct
  structures of 2–128 atoms (123 train, 30 validation), 4.1 atoms on average (RuO2 54, O2 47, Ru2 19, ...).
- Lovelace (read-only checks): the CP2K module `cp2k/2024.3-gnu-openmpi-mkl` has `BASIS_MOLOPT_UZH` (Ru and O:
  DZVP, TZVP, TZV2P-MOLOPT-PBE-GTH) and `POTENTIAL_UZH` (Ru GTH-PBE-q16, O q6), and k-point DOS output.
- Slab vacuum: 22 Å, so that `kpoints_spacing=0.20` gives one k-point along the normal (11×5×1 for 1×1 cells).

**2026-10-18, Part C finished (`validation/ruo2_mace_lovelace/RESULTS.md`):** M2 (the method as proposed) failed
route (a) and route (b); the pre-registered active-learning round M2R passed every gate of route (b) (γ within
0.49 meV/Å², transitions within 0.016 eV, 96–99 % of the CP2K relaxation) and failed route (a) (6–9 meV/Å² for the three seeds). The
outcome rule "only (b) → MACE relaxation + one DFT static" applies, with one active-learning round.

**2026-10-18, Part D written:** the builder D1 follows the validated route rather than the plain M2 recipe: it adds
`active_learning_rounds` (default 1: a DFT static at each relaxed geometry joins the labels, a new fine-tune from
the foundation model, a new relaxation from the starting structures) and `training_structures` (e.g. the rattled
statics of M3), and `dft_check` defaults to on. Unit tests: 247 passed. An end-to-end run of the builder on
Phileas + Lovelace (`c8_builder.py`, 2-pass fine-tunes, two terminations) checks the wiring.

**2026-10-18, the route as a PS-TEROS approach (user: "make this one of the new approaches"):** `MlipSurfaceStudy`
(structures → graph → phase diagram, binary metallic compounds), the recipe guide `docs/MLIP_SURFACES.md` with a
flowchart (`docs/images/mlip_route.py`) and the validation figure (`validation/ruo2_mace_lovelace/c9_figure.py`),
which also sets the MLIP route against Fig. 3 of Reuter and Scheffler (user: "put also the plot by Reuter").

## Steps

### Part A: a CP2K backend for PS-TEROS (needed for the DFT, kept)

Ideas from Kiln's `kiln/tiles/cp2k.py` (defaults for metals, basis and pseudopotential files, Hartree conversion,
fixed atoms) and the legacy `psteros/core/aimd_cp2k.py`; reuse ideas, not imports.

- [x] A1. **Config.** `Cp2kCalculationConfig` (code, `parameters`, basis and pseudopotential files, `kinds`,
  k-point mesh or spacing, `smearing_temperature` in K) next to the QE and VASP configs; `BackendName` gains
  `"cp2k"` and `SurfaceWorkflowConfig` accepts it. Validation names the field and the label.
- [x] A2. **Adapter** `psteros/backends/cp2k.py`: `add_cp2k_task(...)` on aiida-cp2k's `Cp2kBaseWorkChain`, run
  types `ENERGY_FORCE` and `GEO_OPT`/`CELL_OPT`, fixed atoms from `CalculationOverride.fixed_sites`, the forces file
  added to the retrieve list; calcfunctions for the energy (eV), the forces (eV/Å) and the bulk eigenvalues.
- [x] A3. **Builder and reader.** `build_cp2k_relax_static_workgraph(...)`, a new function next to
  `build_qe_relax_static_workgraph`; outputs in the QE pattern (`{label}_relaxed_structure`,
  `{label}_relax_parameters`, `{label}_static_parameters`, `_retrieved`) plus a new `{label}_static_forces`.
  `read_cp2k_results(pk)`. The existing builders dispatch with `if vasp ... else`: they get a guard that raises
  `ValueError` for a CP2K config (pointing to the new builder) instead of sending it into the QE branch. Nothing
  changes for QE or VASP inputs; a test shows their graphs are unchanged.
- [x] A4. **Packaging and docs.** `setup.py` extra `cp2k = ["aiida-cp2k>=2.1"]` (2.1.1 is in `~/envs/aiida`);
  `docs/source/cp2k-workflow.rst` next to the QE and VASP pages; `__all__`; changelog. Tests build graphs with
  `submit=False` and `pytest.importorskip("aiida_cp2k")`.

### Part B: the pure-Python MLIP layer (needed by the pilot, kept by the feature)

- [x] B1. **`psteros/mlip.py`, AiiDA-free.** `LabelledStructure` (label, structure, `energy_ev`,
  `forces_ev_per_angstrom`, role); an extxyz writer with `REF_energy`/`REF_forces`; the least-squares fit of the
  per-element E0 shifts; `MaceFineTuneConfig` and `MlipRelaxConfig` frozen dataclasses with keyword-only fields and
  units in names, validated in `__post_init__`; `mace_training_options(config, ...)` returning the `mace_run_train`
  options of the recipe above; `relax_with_calculator(structure, calculator, *, config, fixed_sites)` returning
  `MlipRelaxation` (structure, `energy_ev`, `max_force_ev_per_angstrom`, steps, converged);
  `partially_occupied_count(occupations)` and `require_metallic(label, partially_occupied=...)`, which raises
  `ValueError` with the label; `check_one_energy_source`; `select_replay_structures`; `prepare_finetune_job` and
  `prepare_relax_job` (job folders). Unit tests with ASE's EMT calculator and hand-made occupations: no torch, no
  AiiDA.
- [x] B2. **Job scripts** in `psteros/mlip_jobs/` (`finetune.py`, `relax.py`): self-contained, need only ase and
  mace-torch, shipped to the job as files. Port the ideas of the goal-07 driver (`run_finetune_replay.py`): replay
  labelling in float64, the head-reproduction check, heads `{Default, pt_head}` with `Default` last, final EMA
  model. Reuse ideas, not imports.
- [x] B3. **`rutile_110_slab(*, metal, termination, triple_layers, vacuum_angstrom, a, c, supercell)`**, a new
  function next to `sno2_110_slab`, which stays as it is. Builds the six terminations above (names
  `o_bridge`, `o_cus`, `metal`, `o_cus_half`, `o_bridge_half`, `o_cus_no_bridge` for `obr`, `ocus`, `ru`, ...).
  Tests: composition, symmetry, and an unchanged `sno2_110_slab` output. Example: `examples/mlip_ruo2/`.

### Part C: the RuO2 pilot (`validation/ruo2_mace_lovelace/`)

- [ ] C0. **Set-up and CP2K checks.** The profile and codes of D4; SHA-256 of the base model and the replay file;
  count the Ru/O replay structures. Then, on `testes` (1 h, empty queue), a smoke test of the CP2K code: a bulk RuO2
  static with k-points that checks the MPI launch (one run, not 128 copies), the energy, the forces file and the
  eigenvalues; the finite-difference force check on one atom (2 extra statics); the metallic guard accepts RuO2 and
  hcp Ru.
- [ ] C1. **CP2K convergence, before any label.** Cutoff (400, 600, 800 Ry; `REL_CUTOFF` 60), basis (DZVP-MOLOPT-SR
  against a larger basis such as TZV2P from `BASIS_MOLOPT_UZH`, where it exists for Ru), k-mesh, electronic
  temperature, vacuum. Converged when the next level changes ΔH_f by ≤ 0.05 eV/f.u. and the differences
  γ(obr) − γ(ocus) and γ(obr) − γ(ru) of unrelaxed slabs by ≤ 2 meV/Å². The settings are then frozen for all of
  Part C and recorded in `PLAN.md` of the validation.
- [ ] C2. **CP2K references and labels.** Bulk RuO2 and hcp Ru `CELL_OPT` -> static; RuO2 in the (110) cell; O2;
  statics of the six unrelaxed terminations (R0, also the training labels); rattled statics if D2.
- [ ] C3. **CP2K ground truth.** `build_cp2k_relax_static_workgraph` on the six terminations (R1). Check G0 here.
  Runs alongside C4–C5.
- [ ] C4. **Fine-tunes on Lovelace A100s.** First a 2-pass run of one seed checks the set-up and the speed (on
  `testegpu` if it fits in 30 min). Then M1, M2 (and M3) × 3 seeds, several seeds in one job when they fit, to spare
  queue slots.
- [ ] C5. **MACE relaxations** of all slabs (and bulk, Ru, O2 for route a) with M0–M3, one GPU job per model.
- [ ] C6. **Route (b) checks.** One CP2K static per model-relaxed slab.
- [ ] C7. **Analysis.** Energies typed into `SlabTermination` / `BinaryOxideReferences` / `surface_phase_diagram`
  (existing pure-Python API), diagrams for R0, R1, M0–M3 in both routes, gates G0–G5, `RESULTS.md`.

### Part D: the PS-TEROS MLIP feature (only if Part C passes, or the user decides to build it anyway)

- [x] D1. **New builder** `build_mlip_surface_workgraph(structures, static, fine_tune, relax, *, mlip_code,
  dft_check=True, submit=False)`: DFT statics through `add_cp2k_task`, `add_vasp_task` or `add_qe_task` (backend
  parity) -> calcfunction that writes the training set -> aiida-shell job for the fine-tune -> one aiida-shell job for
  all relaxations -> optional DFT statics on the relaxed geometries. All tasks wired at build time. New outputs only:
  `{label}_mlip_relaxed_structure`, `{label}_mlip_energy`, `{label}_mlip_check_static_{...}`, `mlip_model`. Existing
  builders and outputs untouched; the remote folders of the model are never cleaned. The fine-tune task waits for
  `require_metallic` on the bulk static. *(Done 2026-10-18 as `build_mlip_surface_workgraph(structures, static,
  fine_tune, relax, *, job, metallic_bulk, fixed_sites, training_structures, active_learning_rounds=1,
  dft_check=True, submit=False)`; outputs `{label}_mlip_relaxed_structure`, `{label}_mlip_energy`,
  `{label}_mlip_dft_energy`, `mlip_round{r}_...`; the remote folder of each fine-tune is linked by its relaxation.)*
- [x] D2. **Reader** `read_mlip_results(pk)` -> `{label: MlipRelaxation}`, working while the graph runs or after it
  fails; a helper that builds terminations and references from one `energy_source` (`"mlip"` or `"mlip_dft"`), and
  raises `ValueError` naming the labels when sources are mixed.
- [x] D3. **Packaging and docs.** `setup.py` extra `mlip = ["ase", "aiida-shell"]` (mace-torch is needed only on the
  GPU host); new names in `__all__`; `docs/source/mlip.rst` with the formula and the one-source rule; example
  `examples/mlip_ruo2/`; changelog entry; a test that the existing call paths give the same results.

### Out of scope

- Semiconductors and insulators, including Ag3PO4: a separate plan by the user.
- More metallic systems after RuO2 (e.g. IrO2(110), O on a metal surface) are proposed after Part C, not now.
- VASP and QE runs: the builders support them (backend parity), but this validation uses CP2K only.

## Where to run

Everything goes on **Lovelace** if its CP2K passes the smoke test (C0): one AiiDA profile, one SSH connection
(2FA once), one daemon worker, and the MACE venvs of the Ag3PO4 project are already there.

**CP2K on Lovelace** (read-only checks, 2026-10-17): the user's own builds registered in `presto`
(`CP2K-2026.1-DEEPMD@lovelace`, `CP2K-2026.2-DEEPMD-SCCS@lovelace`) point to executables that do not exist; the
install folders under `/work/dorinitt/software/cp2k/` are empty. The site provides modules `cp2k/2023.1`,
`cp2k/2023.2`, `cp2k/2024.3-gnu-openmpi-mkl` (and a GPU variant), the `cenapad-*` builds and a 2025.2 container.
Recommended: `cp2k/2024.3-gnu-openmpi-mkl`, registered in the new profile with its `module load` in the prepend
text. GPW, k-points, smearing and `GEO_OPT`/`CELL_OPT` are all long-standing CP2K features. Whatever build passes
C0 is used for every CP2K calculation of Part C.

| Task | Queue | Limits (read-only `qstat -Qf`, 2026-10-17) | Why |
|---|---|---|---|
| CP2K smoke test and the force check (C0) | `testes` | 1 h, 2 per user; empty | small bulk cells |
| All other CP2K: convergence, references, labels, R1, route-(b) statics | `par128`, `max_concurrent_jobs=3` | 6 per user; 3 already used by other projects; 34 running, 150 queued | cells of 6–36 atoms with k-points |
| 2-pass check of the fine-tune (C4) | `testegpu` | 30 min, 1 job per user; empty | measures speed and memory |
| MACE relaxations (C5) | `testegpu` | as above | slabs of ≤ 36 atoms take minutes |
| Fine-tunes M1–M3 × 3 seeds (C4) | `testegpu` if the 2-pass check projects a full fine-tune under 25 min; otherwise `duasgpus`, two seeds at once, one job per model | `duasgpus`: 2 A100, 72 h, 1 running job per user; 2 running, 1 queued | leaves both `umagpu` slots to goal 07 |
| Not used | `umagpu` | 2 running per user, 4 running, 15 queued; goal-07 job 4757 queued | goal 07 needs it for 3 × 20–30 h |

Fallbacks, only if the user approves them:

- **CP2K on Phileas** (`CP2K-2026.2@phileas`, the build that labelled Ag3PO4, 96-core nodes, browser approval at each
  login) if the Lovelace module fails C0 or `par128` waits exceed a day. Its node-hours are budgeted for the Ag3PO4
  goals 04–06.
- **Fine-tunes on Jean Zay V100s** (MACE venv exists, but the `eli` hours belong to the Ag3PO4 project).

Before each launch, the queue limits are read again; the numbers above are a snapshot.

## Cost and time (estimates)

| Item | Jobs | Notes |
|---|---:|---|
| CP2K smoke test and force check | 3 | `testes` |
| CP2K convergence | ≈ 25 | bulk, Ru, O2 and three unrelaxed 1×1 slabs at a few settings |
| CP2K references (cell optimisations, statics, O2) | ≈ 7 | |
| CP2K statics, unrelaxed terminations (R0) | 6 | 18–36 atoms |
| CP2K relax -> static, terminations (R1) | 12 | |
| CP2K statics on model geometries (route b) | 6 per model, 24 for M0–M3 | |
| Rattled statics (D2) | 30 | 5 per termination, σ = 0.05 Å |
| **CP2K total** | **≈ 75–110** | 3 at a time on `par128`; elapsed set by the queue, several days |
| GPU: fine-tunes | 9 (M1–M3 × 3 seeds) | small data (≤ 40 labels of ≤ 36 atoms + a few hundred small replay cells): an estimated 1–4 h each, measured in C4 |
| GPU: relaxations | 4 | minutes each |
| **GPU total** | **≈ 20–40 GPU-h** | on `testegpu`/`duasgpus` |

Code before the first cluster job: Part A (CP2K backend) and Part B (MLIP layer, job scripts, slab builder).

## Decisions (taken 2026-10-17: the recommended option of each)

- **D1, CP2K settings.** (Recommended) PBE without dispersion, GTH-PBE, and the basis, cutoff, k-mesh and electronic
  temperature that C1 selects with the criteria above. Or Kiln's defaults (DZVP-MOLOPT-SR, 400/60 Ry, Fermi-Dirac
  300 K) without a convergence test. Not recommended: O2 is sensitive to the basis.
- **D2, training labels.** (Recommended) the statics of the starting geometries (your idea as stated) for M1 and M2,
  plus 5 rattled statics per termination for M3, to see whether they are worth it. Or the starting geometries only.
- **D3, GPU queues.** (Recommended) `testegpu` and `duasgpus` as in "Where to run", leaving `umagpu` to goal 07.
  Or share `umagpu` with goal 07, or wait for goal 07 to finish.
- **D4, profile and codes.** (Recommended) a new `core.sqlite_dos` profile `psteros_mace_ruo2` with the Lovelace
  computer, the site module `cp2k/2024.3-gnu-openmpi-mkl` as the CP2K code (after C0), and a shell code pointing at
  the existing Ag3PO4 MACE venv (read-only reuse, versions recorded). Or Phileas for CP2K; or first finish the
  user's own CP2K 2026.2 build on Lovelace.
- **D5, order of work.** (Recommended) the PS-TEROS CP2K backend first (Part A), so the pilot runs through PS-TEROS and
  tests the backend on a real cluster at the same time. Or run the pilot's CP2K jobs through Kiln's CP2K stage first
  (faster start; the backend comes after, and the ground truth would be redone through it).
- **D6, branch name.** Rename the branch to `features/mace-surface-phase-diagram`?
