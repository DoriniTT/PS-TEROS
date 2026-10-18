# RuO2(110) with a fine-tuned MACE model: results on Lovelace

Plan, gates and decisions: [`docs/plans/MACE_SURFACE_PLAN.md`](../../docs/plans/MACE_SURFACE_PLAN.md).
AiiDA profile `psteros_mace_ruo2` (core.sqlite_dos, created 2026-10-17); computer `lovelace` copied from
`psteros_vibrations_lovelace` (scheduler `tessera.pbspro_gpu`, `import_sys_environment=False`).

**CP2K moved to Phileas (2026-10-17, user: "the C1 convergence might take forever in lovelace ... launch 2 at a time,
max").** On Lovelace `par128`, C1 had 3 of 42 jobs done after 80 min; Phileas `standard` had idle nodes and no
pending job. Every CP2K energy of the study now comes from `CP2K-2026.2@phileas` (module `cp2k/2026.2`, account
m26257, at most 2 CP2K jobs at a time, `cluster.py`); the Lovelace CP2K runs below are kept as a record only. The
MACE jobs stay on the Lovelace GPUs.

| Code | What |
|---|---|
| `CP2K-2026.2@phileas` | module `cp2k/2026.2` (gcc 15.3, MKL, OpenMPI 4.1.9), `srun`; the prepend checks the executable SHA-256 and runs `cp2k.psmp --check` on the input; copied from profile `ag3po4_h2o_o2_mace_active_learning_stoichiometric`, work directory `/scratch/users/w42008@sorbonne-universite.fr/aiida_psteros_mace_ruo2` |
| `CP2K-2024.3@lovelace` (record only) | site module `cp2k/2024.3-gnu-openmpi-mkl`, `/opt/pub/cp2k/2024.3/gcc/9.4.0/bin/cp2k.psmp`, OpenMPI `mpirun -np`, `OMP_NUM_THREADS=1` |
| `mace-python-cu118@lovelace` | Ag3PO4 MACE venv (`/home/lovelace/proj/proj962/dorinitt/envs/ag3po4_mace_cu118`, mace-torch 0.3.16, torch 2.5.1+cu118), read-only, `core.shell` |

The user's own CP2K 2026.x registrations in `presto` point to binaries that do not exist on Lovelace (checked
2026-10-17); they are not used.

## C0 on Lovelace, CP2K 2024.3 (2026-10-17, `testes`, 4 MPI ranks; record only)

`c0_smoke.py`: bulk RuO2 (6 atoms, experimental lattice) with one O moved 0.05 A along x, and ±0.01 A more;
DZVP-MOLOPT-PBE-GTH (UZH), 500 Ry, k-spacing 0.30 (5×5×7), Fermi-Dirac 300 K. WorkGraph 24, all three jobs
exit 0, about 30 s each (`c0_graph_lovelace_500.json`, `c0_check_lovelace_500.json`).

| Check | Result |
|---|---|
| One CP2K run per job | yes: "Total number of message passing processes 4", CP2K version 2024.3 |
| Energy, force file, DOS file read by the psteros readers | yes |
| Bulk RuO2 metallic (`require_metallic`) | yes: 12 partly filled DOS bins |
| CP2K "Total energy" includes the electronic entropic energy | yes: the printed terms sum to it within 4e-13 Ha; it is the free energy (−TS = −0.37 mHa per cell) |
| −(E₊ − E₋)/(2δ) against the printed force on the moved O | −0.8738 against −0.8846 eV/A: 0.011 eV/A (1.2 %), just above the 0.01 eV/A the script asked for |

The free energy is the energy of the labels, so the 1.2 % gap was attributed to the real-space grid at 500 Ry. The
same check at **800 Ry** (WorkGraph 61, `c0_check_lovelace_800.json`): −0.8085 against −0.8104 eV/A, **0.0019 eV/A (0.23 %)**.
CP2K's energy is the one its forces derive from; the remaining gap at 500 Ry is the grid.

The force itself moved from −0.885 (500 Ry) to −0.810 eV/A (800 Ry), 9 %: forces converge more slowly with the cutoff
than the energy differences of the plan's C1 criteria. C1 therefore also reports the largest force change between
levels (not a gate of the plan).

The Lovelace C1 graph (WorkGraph 335, `c1_graph_lovelace.json`) was stopped when CP2K moved to Phileas.

## C0 on Phileas, CP2K 2026.2 (2026-10-17, `standard`, 8 MPI ranks, 800 Ry)

Same three calculations as on Lovelace, at 800 Ry (WorkGraph 388, `c0_graph_phileas_800.json`,
`c0_check_phileas_800.json`); about 40 s each, started at once.

| Check | Result |
|---|---|
| One CP2K run per job | yes: 8 message-passing processes, CP2K version 2026.2 |
| Readers | yes; the 2026.2 DOS file writes `# E(Fermi) = ... a.u.` on its own line, which the parser now reads (Fermi energy 10.37 eV) |
| Bulk RuO2 metallic | yes: 13 partly filled DOS bins |
| Energy consistent with the forces | −0.7805 (finite difference) against −0.7824 eV/A printed: 0.0019 eV/A |

The force on the moved O differs between the two CP2K builds at the same settings (−0.810 on Lovelace 2024.3,
−0.782 eV/A here), as do the total energies (the data files of the two versions differ): one build for the whole
study is required, which is why every CP2K step runs on Phileas.

## C1: convergence on Phileas (WorkGraph 673, 2026-10-17)

39 CP2K calculations, 2 at a time on one 96-core node each (37 statics of 13 s to 7 min, 2 cell optimisations).
Experimental lattices, unrelaxed slabs; `c1_check.json`. Changes against the level compared with (base: DZVP,
600 Ry, k-spacing 0.20, 300 K, 22 A; t150 is compared with k015):

| Change | ΔH_f (eV/f.u.) | Δ[γ(o_bridge) − γ(o_cus)] | Δ[γ(o_bridge) − γ(metal)] (meV/A²) | Largest force change (eV/A) | Criteria |
|---|---:|---:|---:|---:|---|
| 600 → 400 Ry | +0.013 | −0.09 | +0.25 | 0.200 | pass |
| 600 → 800 Ry | +0.007 | −0.20 | +0.19 | 0.073 | pass |
| DZVP → TZV2P | +0.001 | +1.46 | −0.39 | 0.094 | pass |
| k-spacing 0.20 → 0.15 | +0.001 | −0.02 | +0.20 | 0.040 | pass |
| 300 → 150 K (at 0.15) | +0.001 | −0.04 | +0.06 | 0.031 | pass |
| vacuum 22 → 30 A | 0.000 | −0.02 | +0.04 | — | pass |

Baseline values: ΔH_f = −3.079 eV/f.u. (Reuter and Scheffler −3.16, experiment −3.19); unrelaxed γ at Δμ_O = 0:
o_bridge 83.5, o_cus 46.4, metal 204.6 meV/A² (the literature 49 meV/A² is for relaxed slabs; G0 uses R1).

Bulk CELL_OPT (analytical stress with k-points): at 800 Ry it converged in 9 min, a = 4.536 A, c = 3.136 A; at
600 Ry it stalled after two optimisation steps until the 4 h walltime, and was stopped with the graph.

**Settings chosen (`settings.json`): DZVP-MOLOPT-PBE-GTH (UZH), 800 Ry, REL_CUTOFF 60, k-spacing 0.20, Fermi-Dirac
300 K, 22 A vacuum.** Every level passes at 600 Ry; 800 Ry because the forces still move by 0.073 eV/A from 600 to
800 Ry at no extra cost and because the 600 Ry cell optimisation stalled. TZV2P is not needed by the criteria
(+1.46 meV/A² ≤ 2) and costs twice as much; its force change (0.094 eV/A) does not enter the comparison of the MACE
route with CP2K, which share one recipe.

## C2 and C3 (from 2026-10-17)

References (WorkGraph 919, 39 min, all exit 0): bulk RuO2 and hcp Ru CELL_OPT -> static, O2 GEO_OPT -> static.

| Quantity | CP2K (this work) | Literature |
|---|---|---|
| RuO2 a, c, u | 4.536 A, 3.136 A, 0.3061 | experiment 4.492 A, 3.107 A, 0.306 |
| hcp Ru a, c | 2.729 A, 4.311 A | experiment 2.706 A, 4.282 A |
| O2 bond, ⟨S²⟩ | 1.237 A, 2.006 (triplet) | PBE about 1.23 A |
| ΔH_f(RuO2) | **−3.055 eV/f.u.** (O-poor limit −1.527 eV) | Reuter and Scheffler −3.16, experiment −3.19 |
| Largest force after the static | 0.0001 (bulk), 0.0011 eV/A (O2) | |
| Metallic (partly filled DOS bins) | bulk 13, Ru 18 | |

Labels (WorkGraph 1221: the six terminations as built on this lattice, bulk110, 30 rattled copies) and ground truth
(WorkGraph 1301: GEO_OPT -> static of the six terminations, central trilayer fixed) run together, one CP2K job each.

## C4 check: 2-pass replay fine-tune on Lovelace `testegpu` (ShellJob 1415, 2026-10-17)

One seed, 2 passes, the 7 labels finished at the time (bulk, Ru, O2, o_bridge, metal, o_bridge_half, bulk110) and
the 123 + 30 Ru/O replay structures. Exit 0, 122 s for the whole job.

| Check | Result |
|---|---|
| Foundation model | SHA-256 matches `a522eb7f…` at `/work/dorinitt/goal03_inputs/models/mace-mh-1.model` |
| E0 shifts of the `Default` head (CP2K GTH minus MH-1 PAW) | O −429.35 eV, Ru −2554.04 eV; residual after the shifts 0.069 eV per frame |
| Untrained heads reproduce MH-1 | `Default` 1.5e-13 eV/atom, 2.4e-15 eV/A; `pt_head` 8.9e-16 eV/atom, 4.8e-9 eV/A |
| Validation heads | `pt_head`, `Default` (Default last) |
| Speed | 130 examples per pass in 44 s (float64, batch 1), peak 1.5 GiB |

Projected for the main study (100 passes): M1 about 6 min per seed, M2 about 75 min, M3 about 90 min; on the two
GPUs of `duasgpus` (M2 then M1 on one, M3 on the other) about 5 h, within its 24 h walltime.

## C2/C3 done; G0 and G5 from CP2K alone (2026-10-17, `c7_analysis.py --dft-only`, `c7_dft.json`)

Labels (WorkGraph 1221, 37 statics) and ground truth (WorkGraph 1301, 6 GEO_OPT -> static) finished with exit 0.
Diagrams on Δμ_O in [½ ΔH_f, 0] = [−1.527, 0] eV, references bulk110, hcp Ru and O2 (ΔH_f from bulk110 −3.0549,
from the conventional bulk −3.0546 eV/f.u.: the k-point difference is 0.3 meV/f.u.).

**G0 passes** (o_bridge, o_cus and metal only, as in the paper):

| | CP2K, relaxed (R1) | Reuter and Scheffler | Tolerance |
|---|---:|---:|---:|
| o_bridge -> o_cus transition | −0.879 eV | −0.93 eV | ±0.15 |
| γ(o_bridge) − γ(o_cus) at Δμ_O = 0 | 43.7 meV/A² | 49 meV/A² | ±10 |
| ΔH_f(RuO2) | −3.055 eV/f.u. | −3.16 eV/f.u. | ±0.20 |
| metal stable somewhere in the window | no | no | — |

γ at Δμ_O = 0 (meV/A²), relaxed (R1) and as built (R0), and how much the CP2K relaxation lowers it:

| Termination | R1 | R0 | Lowered by relaxation |
|---|---:|---:|---:|
| o_bridge | 71.9 | 85.1 | 13.2 |
| o_cus | 28.2 | 46.7 | 18.5 |
| metal | 191.0 | 202.7 | 11.8 |
| o_cus_half | 48.9 | 66.1 | 17.2 |
| o_bridge_half | 134.4 | 144.4 | 10.0 |
| o_cus_no_bridge | 132.2 | 149.1 | 16.9 |

R1 with all six terminations: o_bridge -> o_cus_half at −0.924 eV and o_cus_half -> o_cus at −0.835 eV, a window of
0.09 eV in which half coverage is the most stable (a slight repulsion between O-cus atoms; Reuter and Scheffler found
the two O-cus lines crossing o_bridge at the same point). G1 requires both transitions within ±0.05 eV. R0 has a single
transition, o_bridge -> o_cus at −0.772 eV.

**G5: the test is discriminating.** Against R1, R0 fails G1 and G2 (γ off by 10.0 to 18.5 meV/A²): relaxation matters
for every termination, so G3 is gated for all six.

## C4/C5 main study (submitted 2026-10-17 to Bohr `gpu_a100`, user: "for the mace job, I want you to launch at the
Bohr GPU")

ShellJob 1591: M2, M3, M1 (3 seeds each) one after the other on the one A100, then 10 models (M0 and the 9 fine-tunes)
relax 9 structures. Bohr computer copied from `presto` (core.pbspro, `#PBS -l place=excl` as in the earlier Ag3PO4
MACE jobs there); code `mace-python-cu118@bohr` (`ag3po4_g001_mace_cu118`: mace-torch 0.3.16, torch 2.5.1+cu118, ASE
3.29); MH-1 copied to `/home/des01/zacarias/dorinitt/models/mace-mh-1/mace-mh-1.model` (SHA-256 checked). The queue
has one A100 node; at submission another user's 168 h job had run 8.6 h, and one job was queued ahead.

**Moved again (2026-10-18).** After 13 h the Bohr job (PBS 355444) had not started: the one A100 node ran another
user's 168 h job (21.8 h in) with a second 168 h job queued ahead. Cancelled; the user chose "Split into testegpu
chunks".

## C4/C5 in chunks on Lovelace `testegpu` (from 2026-10-18, `c4_chunks.py`)

`testegpu` runs one A100 job of at most 30 min per user (two queued). Every fine-tune (M1, M2, M3 × 3 seeds, 100
passes, unchanged recipe) runs as a chain of jobs, each one chunk of one run of about 20 min of training, sized from
the measured time per pass; the training state and the models stay in `/work/dorinitt/psteros_mace_ruo2`. Then one
relaxation job per model (M0 and the nine fine-tunes).

Chunking is exact. The state carried between jobs holds the raw weights, the EMA, the optimizer and scheduler, the
shuffling generator, the random-number states and the accumulated gradient of a trainable parameter outside the
optimizer: **MACE 0.3.16 clips the gradient norm over every model parameter, while `optimizer.zero_grad` clears only
the optimised ones, so the gradient of the radial Bessel weights (`radial_embedding.bessel_fn.bessel_weights`, not in
any optimizer group) accumulates over the whole training and enters the clipping of every step** (norm 0.49 after 7
steps in a CPU test). Without it a resumed run diverged (8e-4 in the weights after one chunk); with it, three chunks
of a 4-pass replay fine-tune reproduce the uninterrupted run on CPU (energies identical, forces within 4e-19 eV/A,
both heads). The recipe is kept as it is (the Ag3PO4 r4 models were trained the same way); each chunk records that
gradient's norm, since a norm above the clip value of 1 scales every update down.

### First finished model: M2, seed 20260920 (2026-10-18, chunks 0-23, 23-80, 80-100; jobs 1599, 1648, 1667)

Model `/work/dorinitt/psteros_mace_ruo2/models/ruo2_m2_seed20260920.model` (SHA-256 `22f9a861…`). The resumed chunks
loaded their state (start at passes 23 and 80, checkpoint calls 23-79 and 80-99); the exported model is pass 99's
EMA model. Accumulated gradient norm of the Bessel weights: 0.15 after 23 passes, 0.51 after 80 (below the clip
value of 1).

Fit to its own training frames (CP2K labels; reported, not a gate):

| Frame | Energy error (meV/atom) | Force RMSE (meV/A) |
|---|---:|---:|
| bulk | −4.2 | 9.8 |
| ru | +12.1 | 6.7 |
| o2 | −2.7 | 139.9 |
| o_bridge | −13.3 | 21.4 |
| o_cus | −3.5 | 71.2 |
| metal | −20.7 | 12.0 |
| o_cus_half | −9.7 | 51.1 |
| o_bridge_half | −18.4 | 25.5 |
| o_cus_no_bridge | +8.3 | 31.0 |
| bulk110 | −4.1 | 8.1 |

Training curve of the `Default` head on these frames (EMA): energy RMSE 93 meV/atom after pass 0, 23 at pass 40,
14-21 between passes 70 and 95, 11 at pass 99; force RMSE 168 -> 30-43 meV/A. Still decreasing and noisy at pass 100
(batch 1, constant lr 1e-3): 100 passes over 10 frames are 1000 updates of the target head. The replay head drifts:
on the replay validation structures its force RMSE against MH-1 is 376 meV/A after pass 0 and 414 after pass 99
(median over structures 122 meV/A at the end). The r4 recipe was made for 71 frames of 892 atoms; here the target
head sees 10 small frames per pass against 123 replay structures. The plan's recipe stays as fixed; the gates and
the outcome rules decide.

**Extra variant M2L** (user, 2026-10-18: "yes" to a labelled extra variant; fixed before any of its results;
reported, not gated): the M2 recipe with 500 passes and a learning rate decaying from 1e-3 to 1e-5 (gamma 0.990832 per
pass), seed 20260920, after the main study's fine-tunes and next to their relaxations (`c4_chunks.py add-variant`).

## Main study finished; route (a) gates (2026-10-18, `c7_analysis.py`, `c7_analysis.json`)

All 9 fine-tunes (M1, M2, M3 × 3 seeds) and the 10 relaxation jobs finished on `testegpu` without a failure; every
one of the 90 relaxations reached 0.01 eV/A (30-60 FIRE steps). Route (a): every energy of the diagram from the
model itself, at its own geometries, against R1 (CP2K relax -> static); window [−1.527, 0] eV. The gated model is
M2 seed 20260920; the other rows are reported.

| Model | ΔH_f (eV/f.u.) | G4 (≤ 0.05) | G1 | G2: largest error (meV/A², ≤ 5) | o_bridge | o_cus | metal | o_cus_half | o_bridge_half | o_cus_no_bridge |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|
| m0 | -3.185 | fail (-0.130) | fail | 8.1 | 2.8 | 8.1 | 1.5 | 3.9 | 2.4 | 1.1 |
| m1_seed20260920 | -3.070 | pass (-0.015) | fail | 23.3 | 15.4 | 22.1 | 23.3 | 16.6 | 13.3 | 17.3 |
| m1_seed20260921 | -3.077 | pass (-0.022) | fail | 24.5 | 15.7 | 24.5 | 23.7 | 18.7 | 13.8 | 22.0 |
| m1_seed20260922 | -3.058 | pass (-0.004) | fail | 28.6 | 16.3 | 22.5 | 28.6 | 17.1 | 16.1 | 21.4 |
| **m2_seed20260920** | -3.073 | pass (-0.019) | fail | 17.8 | 6.8 | 3.2 | 17.8 | 4.6 | 9.0 | 1.2 |
| m2_seed20260921 | -3.145 | fail (-0.091) | fail | 10.1 | 2.4 | 0.1 | 6.8 | 2.3 | 1.4 | 10.1 |
| m2_seed20260922 | -3.001 | fail (+0.054) | fail | 11.2 | 1.3 | 1.7 | 11.2 | 1.2 | 4.4 | 10.1 |
| m3_seed20260920 | -3.133 | fail (-0.078) | fail | 5.3 | 4.7 | 3.7 | 4.8 | 5.3 | 2.6 | 1.5 |
| m3_seed20260921 | -3.053 | pass (+0.001) | fail | 4.0 | 1.4 | 0.6 | 4.0 | 2.2 | 2.3 | 1.3 |
| m3_seed20260922 | -2.999 | fail (+0.056) | pass | 4.2 | 1.4 | 2.8 | 4.2 | 2.4 | 1.3 | 0.1 |

CP2K: ΔH_f −3.055 eV/f.u.; R1 transitions o_bridge -> o_cus_half at -0.924 eV and o_cus_half -> o_cus at -0.835 eV.
Committee spread of γ at Δμ_O = 0 (standard deviation over the three seeds, meV/A²): M1 2.4 at most; M2 4.5 at most; M3 1.6 at most.

Reading: **the method as proposed (M2) fails route (a)** (G1, and G2 with 17.8 meV/A² on `metal`; G4 passes).
Statics of the starting geometries alone are not enough: the ordinary fine-tune (M1) is the worst model (13 to 29
meV/A² on every termination), worse than MH-1 zero-shot (M0: 1 to 8 meV/A², but ΔH_f off by 0.13 eV); replay (M2)
is much better than M1 but uneven between seeds; adding five rattled statics per termination (M3) gives the best and
most stable models (largest error 4.0 to 5.3 meV/A², spread about 1 meV/A²), with G2 passing for two seeds of three.
G1 is demanding here because of the 0.09 eV o_cus_half window of R1: M3 seed 20260922 reproduces it (−0.885 and
−0.817 eV against −0.924 and −0.835), the other models either shift it by more than 0.05 eV or miss it.


## Route (b): one CP2K static at each model geometry (2026-10-18, WorkGraph 2137, 24 statics, all exit 0)

Deployment seeds; CP2K references (bulk110, hcp Ru, O2); against R1.

| Model | G1 | transitions (eV; R1 −0.924, −0.835) | G2: largest error (meV/A², ≤ 2) | G3 (≥ 90 % of the CP2K relaxation) | recovered fraction per termination |
|---|---|---|---:|---|---|
| m0 | pass | -0.903, -0.859 | 1.52 | fail | o_bridge 0.97, o_cus 0.98, metal 0.87, o_cus_half 0.94, o_bridge_half 0.85, o_cus_no_bridge 0.94 |
| m1_seed20260920 | fail | -1.040, -0.773 | 31.21 | fail | o_bridge -0.69, o_cus -0.13, metal -0.62, o_cus_half -0.13, o_bridge_half -1.18, o_cus_no_bridge -0.84 |
| **m2_seed20260920** | pass | -0.899, -0.823 | 3.36 | fail | o_bridge 0.94, o_cus 0.90, metal 0.71, o_cus_half 0.91, o_bridge_half 0.72, o_cus_no_bridge 0.89 |
| m3_seed20260920 | pass | -0.910, -0.848 | 0.50 | pass | o_bridge 1.00, o_cus 1.00, metal 0.99, o_cus_half 0.98, o_bridge_half 0.95, o_cus_no_bridge 1.00 |

**Outcome for M2 (the plan's rules): it fails route (a) and route (b) (G2 3.36 meV/A² on `metal`, G3 72 % on
`metal` and `o_bridge_half`).** The rule fixed before the results applies: one active-learning round. M2R = M2's
frames plus the six route-(b) statics of M2 seed 20260920 (16 frames), M2's recipe and seeds, relaxed and judged
again; if it fails again the plan stops and reports. Started 2026-10-18 (`c4_chunks.py add-variant m2r`).

Reported: M3 (five rattled statics per termination) passes every gate of route (b) (γ within 0.50 meV/A²,
95-100 % of the relaxation energy) and is closest in route (a); MH-1 zero-shot (M0) with one CP2K static passes
G1 and G2 of route (b) (1.52 meV/A²) but recovers only 85-87 % of the relaxation of `metal` and `o_bridge_half`;
the ordinary fine-tune (M1) moves the slabs away from the CP2K minima (negative recovered fractions).


## Active-learning round M2R (2026-10-18, `c4_chunks.py`, route (b) WorkGraph 2396, 6 statics, all exit 0)

M2R = M2's 10 frames plus the six CP2K statics at the M2 seed-20260920 geometries (WorkGraph 2137), 16 frames, M2's
recipe (100 passes, constant lr 1e-3) and seeds; fine-tunes and relaxations on `testegpu` as before, every relaxation
converged. Route (b) with the deployment seed: one CP2K static at each M2R geometry (`c6_route_b_m2r.json`).

Route (a), against R1:

| Model | ΔH_f (eV/f.u.) | G4 (≤ 0.05) | G1 | model transitions (eV) | G2: largest error (meV/A², ≤ 5) | o_bridge | o_cus | metal | o_cus_half | o_bridge_half | o_cus_no_bridge |
|---|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| **m2r_seed20260920** | -3.036 | pass (+0.018) | fail | o_bridge -> o_cus at -0.967 (no o_cus_half window) | 7.1 | 2.1 | 2.2 | 6.5 | 1.2 | 2.0 | 7.1 |
| m2r_seed20260921 | -3.100 | pass (-0.046) | fail | -1.033, -1.019 | 6.4 | 4.2 | 3.1 | 2.1 | 1.5 | 3.3 | 6.4 |
| m2r_seed20260922 | -2.986 | fail (+0.069) | fail | o_bridge -> o_cus at -1.041 (no o_cus_half window) | 9.0 | 1.0 | 9.0 | 5.6 | 3.4 | 2.5 | 3.9 |

Committee spread of γ at Δμ_O = 0 over the three M2R seeds: 3.0 meV/A² at most (M2: 4.5).

Route (b), against R1:

| Model | G1 | transitions (eV; R1 −0.924, −0.835) | G2: largest error (meV/A², ≤ 2) | G3 (≥ 90 % of the CP2K relaxation) | recovered fraction per termination |
|---|---|---|---:|---|---|
| **m2r_seed20260920** | pass | -0.908, -0.848 | 0.49 | pass | o_bridge 0.99, o_cus 0.99, metal 0.99, o_cus_half 0.97, o_bridge_half 0.96, o_cus_no_bridge 0.99 |

**Outcome (the plan's rules, after the one active-learning round): M2R fails route (a) and passes route (b).** All
gates of route (b) pass (G1, G2 0.49 meV/A², G3 96-99 %, G4 +0.018 eV), with G0 and G5 from CP2K alone. The
validated workflow for RuO2(110) is therefore: replay fine-tune on the statics of the unrelaxed terminations,
one active-learning round with a CP2K static at each model geometry, relaxation with the retrained model on the
GPU, and one CP2K static at each final geometry for the energies of the diagram. Energies taken from the model
itself (route (a)) are not good enough: 6-9 meV/A² on single terminations over the three seeds, enough to close the 0.09 eV o_cus_half
window or move it by 0.1 eV.

DFT cost of the slabs (CP2K time from the output timing report, one Phileas node with 96 ranks; the references
bulk, hcp Ru and O2 are the same for every route and left out):

| Route | CP2K calculations on the six terminations | CP2K time (s) | against R1 | gates |
|---|---|---:|---:|---|
| R1: CP2K relax (7-13 BFGS steps) + static | 6 GEO_OPT + 6 statics | 4859 | 100 % | reference |
| M0 route (b): MH-1 zero-shot + one static | 6 statics | 609 | 13 % | G3 fails (85-87 % on metal, o_bridge_half) |
| M2 route (b) | 6 start statics + 6 final statics | 1137 | 23 % | G2, G3 fail |
| **M2R route (b)** | 6 start + 6 active-learning + 6 final statics | 1629 | 34 % | all pass |
| M3 route (b) | 6 start + 30 rattled + 6 final statics | 3874 | 80 % | all pass |

GPU time on `testegpu` per model: fine-tune 100 passes about 1900 s (M2 1977 s, M2R 1870 s, M3 2458 s in chunks),
relaxation of the six terminations 5-14 min including the queue. M2R needs two fine-tunes (M2, then M2R). The
saving in DFT is modest here because the CP2K relaxations of these small cells (16-38 atoms) converge in 7-13 BFGS
steps from good starting geometries, so one relaxation costs about 7 statics; the method needs 3 statics per
termination. The saving grows with the number of relaxation steps (larger cells, softer surfaces, adsorbates),
since the DFT cost of the method does not depend on it.

## Extra variant M2L (reported, not gated; 2026-10-18, route (b) WorkGraph 2821, 6 statics, all exit 0)

M2's frames and seed 20260920 with 500 passes and the learning rate decaying from 1e-3 to 1e-5 (gamma 0.990832 per
pass), trained in 8 chunks on `testegpu`; every relaxation converged.

| Model | ΔH_f (eV/f.u.) | G4 | route (a) G1 | route (a) transitions (eV) | route (a) G2 (meV/A²) | route (b) G1 | route (b) transitions | route (b) G2 (meV/A²) | route (b) G3: recovered fraction |
|---|---:|---|---|---|---:|---|---|---:|---|
| m2_seed20260920 | -3.073 | pass (-0.019) | fail | -0.833, -0.777 | 17.8 | pass | -0.899, -0.823 | 3.36 | fail: metal 0.71, o_bridge_half 0.72 |
| m2l_seed20260920 | -3.055 | pass (+0.000) | fail | -0.916, -0.911 | 10.6 (metal) | pass | -0.913, -0.833 | 2.78 | fail: metal 0.77, o_bridge_half 0.72 (others 0.93-0.95) |
| **m2r_seed20260920** | -3.036 | pass (+0.018) | fail | -0.967 (one transition) | 7.1 | pass | -0.908, -0.848 | 0.49 | pass: 0.96-0.99 |

Reading: five times more training with a decaying learning rate fixes the formation enthalpy (ΔH_f equal to CP2K's)
and halves the route-(a) error of `metal`, but the relaxations of `metal` and `o_bridge_half` still stop short of the
CP2K minima (72-77 % of the relaxation energy), as with M2. The active-learning statics (M2R), which add labels near
the relaxed geometries, are what make the relaxations right: the missing ingredient was data, not training.

## Part D: the builder end to end (2026-10-18, `c8_builder.py`, WorkGraph 2572, `c8_check.json`)

`psteros.build_mlip_surface_workgraph` on bulk RuO2, hcp Ru, O2, o_bridge and o_cus (CP2K recipe of the labels on
Phileas, at most 2 jobs; 2-pass replay fine-tunes with one seed on `testegpu`; one active-learning round; DFT
check). A wiring check, not a result. The graph finished (exit 0, no task failed; the first fine-tune waited for a
`testegpu` slot and was played again after qsub refused it at the per-user limit of 2). The builder's first statics
give the energies of the C2 graphs for the same structures (difference 0, 3e-10 eV for the bulk); every relaxation of
both rounds converged; the second fine-tune trained on the 5 first statics and the 5 active-learning statics;
`read_mlip_results` and both sources of `read_mlip_energies` read the finished graph. The fine-tune job found the
same 123 + 30 Ru-O replay structures as the study.
