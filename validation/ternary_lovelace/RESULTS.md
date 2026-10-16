# Results: SrTiO3(001) ternary smoke test of the vibrational contributions on Lovelace

Status: **finished** (run 12 to 16 October 2026; design and checks T1 to T8 in PLAN.md).

**Summary.** The vibrational module works with a ternary oxide: the vibrations graph, the readers, `solid_free_energy_ev` with the frozen centre of the slab counted as one SrTiO3
unit, and `TernaryOxideReferences` / `ternary_surface_phase_diagram` on free energies all ran on real data, the checks of the mechanics pass (T2 to T7), and the static ternary
diagram is byte-identical to the one of the code before the feature. What this smoke test cannot give are converged numbers: the Gamma-only supercells are tiny, and the cubic
SrTiO3 and the TiO2-terminated slab have imaginary modes whose handling decides the diagram (section 'Results'). T8 fails on its numeric threshold.

Same profile, computer, code, queue and settings as `validation/vibrations_lovelace/RESULTS.md` (section 1), including aiida-workgraph 0.9.0
through the overlay and the daemon started from this worktree. `kpoints_spacing=0.04`. POTCARs `Sr_sv`, `Ti_pv`, `O`.

## Process keys

| Graph | PK | State |
|---|---:|---|
| refs (5 solids, 10 jobs) | 6819 | **Finished [0]**, 10 of 10 calculations exit 0 (47 min, mostly queue) |
| slabs (2 terminations, 4 jobs; lattice of 6819, a = 3.908 A) | 6965 | **Finished [0]**, 4 of 4 exit 0 (about 5 h of queue and run time on 12 Oct) |
| vibrations, attempt 1 (282 jobs, one at a time) | 9561 | submitted 12 Oct 13:56; 11 calculations finished by 13 Oct 07:23, with queue waits of 22, 70, 114, 173, 184, 23, 4, 4, 3, 61 and 17 min (mean 61 min; `par128` had 118 waiting and 32 running); the host restarted on 13 Oct about 07:25 and the graph then stalled. Cancelled by me on 13 Oct at 10:00 (the kill of a graph that is still loading its checkpoint has to be forced, and it left one work chain, 10253, and its calculation, 10258, which I cancelled too) |
| vibrations, attempt 2 (282 jobs, `--max-concurrent 8`) | 12807 | submitted 13 Oct 10:14, ended 16 Oct 20:18 `[302]` (three mode tasks failed or were skipped, see T-2); **283 VASP calculations, all exit 0** (one task, `tio2_bulk_vib_s1_z_minus_vasp`, ran twice, work chains 16192 and 16195, because two daemon workers acted on the graph, T-4); I chose 8 at a time after attempt 1 showed about 12 days for one at a time. Stalled from 13 Oct 22:01 to 15 Oct 13:16 (see 'Problems'), then 4 to 18 jobs per hour; running |
| srtio3_bulk modes, rerun of the failed task (`rerun_modes.py`) | 15849 | **Finished [0]** (calcfunction 15850), from the 60 retrieved folders |
| ti_metal modes, rerun | 16119 | **Finished [0]**, from the 24 folders |
| tio2_bulk modes, rerun (task SKIPPED in 12807) | 16505 | **Finished [0]**, from the 36 folders |

The modes of `sr_metal`, `sro_bulk`, `slab_tio2` and `slab_sro` come from graph 12807 itself; `analyse_ternary.py --rerun-pk 15849 16119 16505` reads the other three.

O2 is taken from the SnO2 test: energy of graph 977, vibrations of graph 3400.

## References and slabs (relaxed, static energies)

| Label | E (eV) | cell |
|---|---:|---|
| srtio3_bulk (Sr Ti O3) | -40.0271 | a = 3.908 A (PBE about 3.94: 0.8 % small at ENCUT 400) |
| sr_metal (Sr4, fcc) | -6.5489 | a = 6.017 A |
| ti_metal (Ti2, hcp) | -15.5992 | a = 2.942 A, c = 4.635 A |
| sro_bulk (Sr4 O4) | -48.2419 | a = 5.151 A |
| tio2_bulk (Ti2 O4, rutile) | -53.5035 | a = 4.596 A, c = 2.955 A |
| slab_tio2 (Sr2 Ti3 O8) | -105.0667 | 3.908 x 3.908 x 17.816 A |
| slab_sro (Sr3 Ti2 O7) | -91.2045 | 3.908 x 3.908 x 17.816 A |

k-meshes written: SrTiO3 7x7x7, Sr 5x5x5, SrO 5x5x5, Ti 10x10x6, TiO2 6x6x9, slabs 7x7x2. Slab relaxations use selective dynamics; the
statics do not. Delta H_f(SrTiO3) from these energies and E(O2) = -9.8568 eV: -15.8 eV (raw PBE; experiment about -17.3 eV).

## Checks

| # | Check | Result | Evidence |
|---|---|---|---|
| T1 | Graphs finish | **PASS with deviations** | `refs` and `slabs` exit 0; in `vibrations` all 283 calculations exit 0, but three of the seven mode tasks failed or were skipped by the engine (T-2) and were run again on the same inputs (graphs 15849, 16119, 16505); `graph_outputs` was skipped and the graph ended `[302]`; one displacement ran twice (T-4) |
| T2 | Displaced calculations | **PASS** | all 282 `displace_site` outputs move one coordinate of one site by 0.01 A; no selective dynamics; `NSW = 0`; stored parameters equal the recipe; `import_sys_environment = False` |
| T3 | Mode counts | **PASS** | SrTiO3 27, TiO2 15, SrO 21, Sr 9, Ti 9, slab_tio2 24, slab_sro 21; zero modes removed 3 for every periodic solid, 0 for the slabs |
| T4 | Frozen centre | **PASS** | both slabs freeze Sr1 Ti1 O3 and `solid_free_energy_ev(..., bulk=SrTiO3)` adds exactly 1.000 bulk cell; a frozen region Sr1 Ti2 O3 is refused with a message naming the stoichiometry |
| T5 | Imaginary modes | **PASS for the handling, not small** | `read_vibrations` raises by default exactly for the two structures that have imaginary modes (SrTiO3: 3, at -61.9, -60.2, -60.2 cm-1; slab_tio2: 4, at -206.5, -206.5, -167.7, -167.7) and works with `imaginary_modes="drop"`; the modes are not small (SnO2 reference: < 50i) |
| T6 | Independent recomputation | **PASS** | slab_tio2 recomputed from the retrieved `vasprun.xml` files with pymatgen: largest difference 4e-12 cm-1 |
| T7 | Backward compatibility | **PASS** | the static ternary diagram CSV is byte-identical to the one of the commit before the vibrations (`8b74134`) run on the same numbers |
| T8 | Physical sanity | **FAIL on the threshold** | Delta H_f(SrTiO3) = -15.805 eV static (raw PBE, in range) but it moves by up to 0.357 eV at 900 K (limit 0.3; for SnO2 the shift was 0.11 eV); a stability polygon exists at every T. Analysis below |

## Results

Modes (graph 12807 and the reruns):

| Label | modes | range (cm-1) | imaginary | ZPE (eV) |
|---|---:|---|---:|---:|
| srtio3_bulk (1x1x2, 10 atoms) | 27 | -61.9 to 755.0 | 3 | 0.2246 per supercell |
| sr_metal (fcc, 4 atoms) | 9 | 77.8 to 103.4 | 0 | 0.0482 |
| ti_metal (hcp, 1x1x2, 4 atoms) | 9 | 91.2 to 211.3 | 0 | 0.0379 |
| sro_bulk (8 atoms) | 21 | 124.2 to 312.9 | 0 | 0.2688 |
| tio2_bulk (rutile, 6 atoms) | 15 | 79.1 to 806.9 | 0 | 0.3308 |
| slab_tio2 (Sr2Ti3O8, 8 free atoms) | 24 | -206.5 to 775.7 | 4 | 0.3496 |
| slab_sro (Sr3Ti2O7, 7 free atoms) | 21 | 55.6 to 583.8 | 0 | 0.3183 |

Ternary phase diagram (`results/ternary_*.png`, `ternary_table.csv`); fraction of the stability region where each termination is the most stable:

| Case | Delta G_f(SrTiO3) (eV) | TiO2-terminated | SrO-terminated |
|---|---:|---:|---:|
| static energies | -15.805 | 15.0 % | 85.0 % |
| 0 K (ZPE), imaginary modes dropped | -15.757 | 17.9 % | 82.1 % |
| 300 K, dropped | -15.778 | 14.5 % | 85.5 % |
| 600 K, dropped | -15.925 | 7.5 % | 92.5 % |
| 900 K, dropped | -16.162 | **0.0 %** | 100 % |
| 0 K, imaginary modes kept as |omega| | -15.751 | 15.9 % | 84.1 % |
| 300 K, kept | -15.826 | 16.0 % | 84.0 % |
| 600 K, kept | -16.074 | 17.9 % | 82.1 % |
| 900 K, kept | -16.433 | 20.3 % | 79.7 % |

The SrO-terminated slab is the more stable in most of the region at every temperature, as the static energies already say. With the imaginary modes dropped (the only way to
obtain a free energy from them) the TiO2-terminated region shrinks and vanishes at 900 K; with the same modes counted as real (`ternary_sensitivity_imaginary_as_real.csv`,
a sensitivity and not a recommendation) it stays at 16 to 20 %. So the disappearance is a consequence of dropping four large imaginary modes of `slab_tio2`
(-207i and -168i cm-1: the ideal 1x1 TiO2-terminated surface is unstable to a distortion that the symmetric relaxation does not follow) and of the three of the bulk, not
a result of the vibrations. `figure`: `results/ternary_comparison.png` (static, 900 K dropped, 900 K kept).

Why T8 fails and what it means: Delta G_f(T) is more sensitive here than for SnO2 because the references are 4-atom Gamma-only cells of metals (Sr: 9 modes between 78 and 103 cm-1,
no acoustic region) and because the bulk has imaginary modes whose free energy is left out when they are dropped. Both are limits of this smoke test; the
magnitude 0.36 eV is not a converged result. The threshold of 0.3 eV was my guess in the plan.

Time: refs 47 min (queue), slabs 5 h (queue), vibrations 3.5 days of which about 2 days lost to the stalls of the problems below; at the end the displaced calculations ran 4 to 18 per hour.

## Problems

### T-1 The SSH 2FA session dropped (13 Oct 22:01) and the daemon did not recover by itself

All transport tasks failed five times and paused; nothing was in PBS. After you logged in again (14 Oct 09:00) `verdi process play --all`
was not enough: the worker stayed hung for another 28 h (my status watcher printed empty lines and I did not notice) and only a daemon
restart (15 Oct 13:16, with `verdi process repair` for 7 processes that had lost their tasks) made the graph move again. The resume script
now restarts the daemon, an SSH-loss watcher alerts immediately, and the graph watcher alerts when its own status is empty.

### T-2 A mode task was started with 58 of its 60 retrieved folders (16 Oct 07:21)

`srtio3_bulk_vibrations` excepted with `KeyError: 's9_z_plus'`: its inputs lacked the folders of `s9_z_plus` and `s9_z_minus`, although all 60
work chains had finished with exit 0 (the last two at 07:17:48 and 07:20:23, recorded by one daemon worker at 07:17:50 and 07:20:58 while
the other worker, 17467, which had restored the graph from its checkpoint at 07:18:33, launched the task at 07:21:28). With
`max_concurrent_jobs=1` and one worker (the SnO2 test) this did not happen. The graph went on, but `graph_outputs` was skipped.
`rerun_modes.py` (graph 15849) ran the same calcfunction on the 60 folders: 27 modes, 3 imaginary (-61.9, -60.2, -60.2 cm-1; the cubic
perovskite is not a minimum in PBE), the others from 105 to 755 cm-1.

The same happened to `ti_metal_vibrations` (21 of 24 folders; rerun graph 16119) and `tio2_bulk_vibrations`; `sr_metal`, `sro_bulk` and `slab_tio2`
finished in the graph. Cause, found afterwards in T-4: two daemon workers on one SQLite database.

psteros change (commit `756bfc5`): `harmonic_modes` now names the missing displacements in a `ValueError` instead of a bare `KeyError`;
new test that fails without the change. The race itself is in the engine, not in psteros.

### T-3 `qsub` refused my jobs: per-user limit of the `par128` queue

`qstat -Qf par128` shows `max_queued = [u:PBS_GENERIC=6]` and `max_run = [u:PBS_GENERIC=6]` (and 12 for the group). The limit counts all of
the user's jobs in the queue, here shared with other AiiDA profiles (3 to 4 jobs of other projects in `par128` at the time). With
`max_concurrent_jobs=8` several submissions failed with `qsub: would exceed queue generic's per-user limit`; AiiDA retried 5 times and
paused each calculation (`Pausing after failed transport task: submit_calculation failed 5 times consecutively`). A graph meant to run
8 jobs at a time on this queue can run at most 6 minus the user's other jobs; **`max_concurrent_jobs=8` was a bad choice**, mine.

### T-4 Two daemon workers on a SQLite profile

I raised the number of workers to 2 (`verdi daemon incr 1`) to speed up the start of the big graph. The profile uses `core.sqlite_dos`, and the worker log
shows `sqlite3.OperationalError: database is locked`, and both workers running the same WorkGraph at the same time (every `update_task_state`,
`continue_workgraph` and `on_wait` line twice, `Restored missing live awaitable ...`). That is how `srtio3_bulk_vibrations`, `ti_metal_vibrations` and
`tio2_bulk_vibrations` were started with some of their retrieved folders missing (T-2). Back to one worker on 16 Oct 14:56.
After the restart, the play requests did not reach the worker for 4 h while it loaded the graph (the command timed out), so the last 5 calculations
(four of them already finished on the cluster, in PBS since 13 Oct) stayed unretrieved until I played them again once the worker was idle.

## What this says about the branch

* **Ternary compounds work** with the same blocks as binary ones; nothing in `HarmonicVibrations`, `solid_free_energy_ev`, `build_vibrations_workgraph` or `read_vibrations` assumes two elements.
  The user has to substitute free energies for the bulk, the per-atom energies of the two elements, O2 (plus its zero-point energy) and every competing phase, as `analyse_ternary.py` shows.
* **The imaginary-mode policy decides the result** for compounds whose ideal structure is not a minimum (SrTiO3, SrTiO3(001) surfaces). The default (raise) is the right one; `imaginary_modes="drop"`
  changes the phase diagram here from a 15 to 20 % region to none, so the documentation should say that dropping is a choice with consequences and recommend checking the sensitivity.
* **Operational lessons for long vibrations graphs**, valid for the SnO2 test as well: one daemon worker only (SQLite), `max_concurrent_jobs` at most the per-user limit of the queue minus the user's other jobs,
  a connection watcher and a recovery procedure (restart the daemon, `verdi process repair`, play the processes), and a check of the task states of the WorkGraph, not of the processes only.
