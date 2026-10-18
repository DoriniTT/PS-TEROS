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
