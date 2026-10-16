# Results: vibrational contributions on Lovelace (par128)

Status: **finished** (run 8 to 12 October 2026). Verdict in section 8: **not ready to merge as it stands**; the feature
works end to end on Lovelace, but `main` has moved and the physical convergence of the vibrational term of gamma is not shown.

> **Units.** The calculations below ran on the branch before it was merged with the k-point fix of `main`; there `kpoints_spacing` reached aiida-vasp unconverted, so the recorded `0.04` is 0.04 x 2 pi = 0.25 A^-1 in the units of `main` (psteros values are A^-1 with the 2 pi of VASP's `KSPACING`). The scripts in this folder use 0.25.

## 1. Environment and settings actually used

| Item | Value |
|---|---|
| Branch | `features/vibrational-contributions` (base `47b1eef`) |
| Python | 3.13.13 (`~/envs/aiida`) |
| aiida-core / aiida-vasp / pymatgen / ase | 2.8.0 / 5.0.0 / 2026.3.23 / 3.26.0 |
| aiida-workgraph | 0.8.1 installed in the shared env; **0.9.0 used for every run** through a `PYTHONPATH` overlay (see 8.1) |
| psteros | `__version__` 2.0.0, taken from this worktree via `PYTHONPATH`; the shared env has it editable from another checkout and was not modified |
| AiiDA profile | `psteros_vibrations_lovelace`, **new**, created for this test (sqlite_dos, RabbitMQ); default profile unchanged |
| Computer | `lovelace`, a copy of the registration in profile `presto` (`core.ssh` through the cenapad ControlMaster, scheduler `tessera.pbspro_gpu`, work dir `/work/dorinitt/.aiida`), **with `mpirun -np {tot_num_mpiprocs}` instead of the OpenMPI 5.0.6 path of `presto`** (see 8.2) |
| Code | `VASP-6.5.1@lovelace` (`vasp_std`, copy of the `presto` registration; intel/2023.2.1 module, OpenMPI 5.0.6 mpirun) |
| POTCARs | family `PBE` (329 POTCARs copied from `presto`), mapping `Sn -> Sn_d`, `O -> O` |
| Daemon | started with `PYTHONPATH=<overlay>:<worktree>`; one worker |
| Queue / resources | `queue="par128"`, `{"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128}`, `custom_scheduler_commands="#PBS -j oe"`, `max_concurrent_jobs=1` |
| `import_sys_environment` | `False`, through `CalculationOverride(metadata=...)` for every label |
| Walltime | 7200 s for `refs` and `slabs` (relax and static share one `ExecutionPolicy`), 3600 s for `vibrations` |
| INCAR | `ENCUT=400, PREC=Accurate, EDIFF=1e-6, ISMEAR=0, SIGMA=0.05, LREAL=False, LWAVE=LCHARG=False, NCORE=16`; relax `IBRION=2, NSW=100, EDIFFG=-0.01`, `ISIF=3` for bulk and alpha-Sn, `ISIF=2` otherwise; static `IBRION=-1, NSW=0`; displaced `EDIFF=1e-7`; O2 `ISPIN=2, MAGMOM=[1,1]`, Gamma only |
| k-points | `kpoints_spacing=0.04` (aiida-vasp units, x 2 pi; the plan's 0.3 is Gamma-only, section 8.3): alpha-Sn 4x4x4, rutile 6x6x8, slabs 8x4x2 (unrelaxed cells); O2 `kpoints_distance=10` (Gamma). Graph 862 used 0.3 |
| Restarts | `max_iterations=3` |
| Vibrations | `VibrationsConfig(displacement_angstrom=0.01, supercells={"sno2_bulk": (1, 1, 2)}, molecules=("o2",))` |

Job script produced for the first job (`_aiidasubmit.sh`, PK 805): `#PBS -q par128`,
`#PBS -l walltime=02:00:00`, `#PBS -l select=1:mpiprocs=128:ncpus=128`, `#PBS -j oe`, and no
`#PBS -V`, so `import_sys_environment=False` took effect.

## 2. Process keys

| Graph | PK | State |
|---|---:|---|
| refs, attempt 1 (workgraph 0.8.1, psteros as on the branch) | 705 | Finished [302] after 26 s: INCAR case error (8.1) |
| refs, attempt 2 (workgraph 0.9.0, psteros as on the branch) | 751 | Finished [302] after 25 s: the same error (8.1) |
| refs, attempt 3 (with `PsterosVaspWorkChain`) | 797 | Killed by me after its first job (calc 805, alpha-Sn relax) failed with exit 1002: the job ran on **1 MPI rank** (8.2); calc 815 (queued second job) cancelled with 810/815 |
| refs, attempt 4 (computer `mpirun` fixed) | 862 | **Finished [0]**, all 6 children exit 0 (PBS jobs 1028184, 1028197, 1028211, 1028247, 1028256, 1028257; calcs 870, 882, 893, 905, 916, 928) |
| refs, attempt 5 (`kpoints_spacing=0.04`) | 977 | **Finished [0]** in 18 min (queue nearly empty), all 6 children exit 0 (calcs 985, 997, 1008, 1020, 1031, 1043). This is the `refs` graph used by the rest of the test; 862 stays as the 0.3 record (8.3) |
| slabs (on the lattice of 977: a = 4.7652, c = 3.2219 A) | 1077 | **Finished [0]** in 32 min, 4 children exit 0 (calcs 1086, 1098, 1110, 1122; PBS 1028288, 1028289, 1028291, 1028296) |
| vibrations (refs 977, slabs 1077) | 3400 | submitted 8 Oct 18:48, ended 12 Oct 03:42: **252 of 252 VASP calculations exit 0**; the graph ended `[302]` because the engine marked one task FAILED after a host restart and skipped `alpha_sn_vibrations` (8.5), which graph 4962 supplies |
| alpha-Sn modes, rerun of the skipped task (`rerun_alpha_sn_modes.py`) | 4962 | **Finished [0]** (calcfunction 4963) |
| ibrion5, attempt 1 | 6697 | Finished [302]: VASP stops with "VASP internal routines have requested a change of the k-point set ... NPAR" (8.6) |
| ibrion5, attempt 2 (`ISYM=0`) | 6729 | **Finished [0]** (calc 6737) |

## 3. Timings of the `refs` graph (PBS `qtime`/`stime`/`resources_used.walltime`)

| Calc | Job | Queue wait | VASP wall time |
|---|---|---:|---:|
| 870 | alpha-Sn relax (4 ionic steps) | 25.1 min | 4 min 23 s |
| 882 | alpha-Sn static | 42.5 min | 6 s |
| 893 | O2 relax | 117 min | 51 s |
| 905 | O2 static | 16.4 min | 1 min 50 s |
| 916 | rutile SnO2 relax | 0.0 min | 50 s |
| 928 | rutile SnO2 static | 34.8 min | 1 min 21 s |

Queue wait: mean 39 min, median 30 min, range 0 to 117 min (par128 had about 110 jobs queued and 31 running when the
test started). VASP: 93 s on average, 9 min in total for 6 jobs. Graph wall time 13:15 to 17:41 (4 h 26 min); besides
queue and VASP time there is about 3.5 min per job of AiiDA overhead (upload, submit, polling, retrieve).

## 4. Results of the `refs` graph

| Quantity | Result | Plan |
|---|---|---|
| rutile a, c | 4.7674 A, 3.1686 A (V = 72.02 A^3) | a = 4.8, c = 3.2 (PBE) |
| O2 magnetisation | 1.9999987 muB | 2 muB |
| O-O bond | 1.234 A | 1.22 to 1.24 A |
| alpha-Sn a | **7.108 A** (start 6.489) | PBE about 6.6 A |
| E(SnO2), E(alpha-Sn, 8 atoms), E(O2) | -36.3567, -25.4398, -9.8568 eV | |

The alpha-Sn lattice is 7% above the PBE value and the formation energy is too negative (Delta H_f about -5.1 eV per SnO2,
PBE with Sn_d about -4.6 eV): see 8.3.

### Graph 977 (`kpoints_spacing=0.04`), the references used from here on

| Quantity | Result | Plan / reference |
|---|---|---|
| rutile a, c | 4.7652 A, 3.2219 A (V = 73.16 A^3) | a = 4.8, c = 3.2 (PBE): pass |
| alpha-Sn a | 6.647 A | PBE about 6.65 A |
| O2 magnetisation, O-O | 1.9999987 muB, 1.2343 A | 2 muB: pass |
| E(SnO2), E(alpha-Sn, 8 atoms), E(O2) | -37.3318, -30.7417, -9.8568 eV | |
| Delta H_f per SnO2 | -4.97 eV | PBE with Sn_d about -4.6 to -5.0 eV (experiment -6.0) |
| k-meshes written | alpha-Sn 4x4x4, O2 1x1x1, rutile 6x6x8 | |


## 5. Checks (PLAN.md)

| # | Check | Result | Evidence |
|---|---|---|---|
| 1 | Unit tests on the local stack | **PASS** | 148 passed on the branch as found, 154 with the 6 new tests; the same on aiida-workgraph 0.8.1 and 0.9.0 |
| 2 | Merge with `main` | **FAIL** | against `origin/main` of 8 Oct (`ff0cd3e`, one commit) only `.gitignore` conflicted and the merged tree passed 148 tests; against `origin/main` of 12 Oct (`f7bdfd7`, 22 commits beyond the common base `6d8e669`; main was rebased in between, `ff0cd3e` became `e247edd`) **9 files conflict (22 hunks)**: `.gitignore`, `AGENTS.md`, `CHANGE.md`, `docs/source/{api,index,phase-diagram}.rst`, `psteros/__init__.py`, `psteros/backends/vasp.py`, `psteros/config.py`. Tests on the merged tree were not run (conflicts unresolved). See 8 |
| 3 | All graphs finish | **PASS with a deviation** | `refs` 977 and `slabs` 1077 exit 0 for every child. `vibrations` 3400: all 252 VASP calculations and work chains exit 0, but the graph ended `[302]` (8.5): one task was marked FAILED by the engine after a host restart and its dependent mode task skipped; the same calcfunction was run on the same inputs in graph 4962 (exit 0) |
| 4 | Displaced calculations | **PASS** | all 252 `displace_site` outputs differ from their input in exactly one coordinate of one site by 0.01 A (tolerance 1e-9); every VASP work chain takes a `displace_site` output; none has selective dynamics (parameters or POSCAR); `NSW = 0` in all 252 written INCARs; the stored parameters equal the vibrations static recipe; `import_sys_environment = False` |
| 5 | Mode counts | **PASS** | slab_o 36, slab_sn2o 24, bulk 33, alpha-Sn 21, O2 1 (`zero_modes` 5; bulk and alpha-Sn 3; slabs 0) |
| 6 | O2 stretch | **PASS** | 1563.415 cm-1 (PBE 1550 to 1600); VASP `IBRION=5`, `NFREE=2`, `POTIM=0.01`: 1563.386 cm-1, difference 0.03 cm-1 (tolerance 10) |
| 7 | Solids | **PASS** | bulk modes 91.7 to 747.6 cm-1 (within 0 to 800); no imaginary mode in any of the 5 labels. Because there is none, `read_vibrations` raising by default and working with `imaginary_modes="drop"` could only be shown on real data for the case "no imaginary mode" (default read ok, same number of modes with `drop`); the raising path is covered by the unit tests only |
| 8 | Independent recomputation | **PASS** | slab_sn2o (24 modes, 48 vasprun.xml files) and alpha-Sn (21 modes, 48 files) recomputed with pymatgen's parser and masses and `harmonic_vibrations_from_forces`: largest difference 4e-12 and 1e-5 cm-1 (tolerance 0.1) |
| 9 | Backward compatibility | **PASS** | the 0 K diagram from total energies is byte-identical (a) to the one of `examples/vasp_surface_phase_diagram/phase_diagram.py` without `--vib-pk`, and (b) to the one of the commit before the vibrations (`8b74134`) run in a separate process on the same energies |
| 10 | Physical sanity | **FAIL on the first criterion** | gamma(slab_o) changes by -7.1 meV/A2 at 300 K and **-25.6 meV/A2 at 900 K** (plan: a few at most); the transition Delta mu_O shifts by 0.12 eV at most (plan: <= 0.2, **pass**); the O-poor limit moves with Delta G_f(T) (**pass**). Analysis in 6 |

## 6. Vibrational results

Static energies (eV, graphs 977 and 1077): SnO2 bulk cell -37.3318, alpha-Sn (8 atoms) -30.7417, O2 -9.8568, slab_o (Sn6O12) -109.1527, slab_sn2o (Sn6O8) -82.0010.

| Label | modes | zero modes removed | lowest, highest (cm-1) | imaginary | ZPE (eV) |
|---|---:|---:|---|---:|---:|
| sno2_bulk (1x1x2 supercell, 12 atoms) | 33 | 3 | 91.7, 747.6 | 0 | 0.4044 per supercell |
| alpha_sn (8 atoms) | 21 | 3 | 38.1, 173.9 | 0 | 0.1611 |
| o2 | 1 | 5 | 1563.4 | 0 | 0.0969 |
| slab_o (12 free sites) | 36 | 0 | 39.6, 717.5 | 0 | 0.7238 |
| slab_sn2o (8 free sites) | 24 | 0 | 34.6, 614.1 | 0 | 0.3656 |

| Case | Delta G_f(SnO2) (eV) | O-poor limit (eV) | gamma slab_o at 0 (J/m2) | gamma slab_sn2o at 0 (J/m2) | transition Delta mu_O (eV) |
|---|---:|---:|---:|---:|---:|
| static E | -4.966 | -2.483 | 1.049 | 3.793 | -1.860 |
| 0 K (ZPE) | -4.881 | -2.441 | 1.017 | 3.701 | -1.818 |
| 300 K | -4.875 | -2.438 | 0.936 | 3.637 | -1.831 |
| 600 K | -4.943 | -2.471 | 0.796 | 3.581 | -1.888 |
| 900 K | -5.075 | -2.537 | 0.638 | 3.560 | -1.980 |

The stoichiometric `o` termination is stable above the transition and the reduced `sn2o` below it at every temperature; the
transition stays inside the stability window. The ZPE alone changes the static results by -3 % (gamma slab_o) and +0.085 eV (Delta G_f).

**Why check 10 fails, and what it does and does not show** (`results/gamma_vib_decomposition.csv`). At Delta mu_O = 0 the
stoichiometric slab has gamma_vib = [F_vib(12 free atoms) - 12 x F_vib(bulk per atom)] / 2A. At 900 K: F_vib of the free atoms
is -217.6 meV/atom, that of the bulk -124.8 meV/atom. The slab has 6 modes below 100 cm-1 (lowest 39.6), the bulk 1 (lowest 91.7):
the bulk, a Gamma-only 1x1x2 supercell of the 6-atom cell, samples only Gamma and Z and has no acoustic region, which is the
documented limit of the method (docs `vibrations.rst`, Limits) at the smallest supercell. Raising the slab modes to the lowest
bulk mode lowers the effect only from -25.6 to -20.1 meV/A2 (-13.1 at a floor of 150 cm-1), so softer surface atoms are part of it.
Both parts are in the sign expected for a surface; the magnitude is not established. A 2x2x2 bulk supercell (the default of
`examples/vasp_surface_phase_diagram/vibrations.py`, 48 atoms, 288 jobs) would show how much of it is the sampling: **not run**.
It does not affect checks 4 to 9 (the numbers are consistent between the graph, the independent recomputation and the formulas).

### Figures and tables (`results/`)

`phase_diagram_total_energy.{png,csv}`, `phase_diagram_T0K`, `T300K`, `T600K`, `T900K` (`.png`, `.csv`), `gamma_table.csv`,
`gamma_vs_delta_mu.png`, `vibrational_change.png`, `frequencies_<label>.csv` (every frequency, 5 labels), `mode_summary.csv`,
`gamma_vib_decomposition.csv`, `checks.json`.

## 7. Timings

| Graph | Jobs | Wall time | Notes |
|---|---:|---|---|
| refs 977 | 6 | 18 min | queue nearly empty |
| slabs 1077 | 4 | 32 min | relaxations 6.4 and 8.8 min, statics 55 s and 1.8 min |
| vibrations 3400 | 252 | 81 h (8 Oct 18:42 to 12 Oct 03:42) | about 20 h of it lost to two host restarts; sum of create-to-end times of the calculations 48 h |
| ibrion5 6729 | 1 | 5 min | |

Displaced calculations (252): median 3.9 min from creation to retrieval, mean 11.4, 90th percentile 7.9, maximum 309 min;
235 took under 30 min (mean 4.3), 17 over 30 min and 6 over 2 h, the latter in the queue (`par128` was busy at night:
118 queued, 31 running). The engine adds about 2 min between jobs with `max_concurrent_jobs=1`. The 4.4 h of graph 862 (6 jobs) and
the 10 min of graph 977 show that the same graph takes anywhere from 0.1 to 8 days (252 jobs at 3 to 45 min each) depending on the queue.

## 8. Problems found

### 8.1 psteros VASP graphs cannot start with aiida-vasp 5.0.0 (INCAR keys in upper case): bug in psteros, fixed

Symptom: every graph failed after 25 s, before any job was submitted (graphs 705 and 751; the same on aiida-workgraph 0.8.1 and
0.9.0): `InputValidationError: Case inconsistency found in the parameters dictionary please use lower case keys: Key 'ENCUT'
converted to 'encut'`. Cause: psteros recipes, the examples and the docs write the INCAR in upper case and
`vasp_parameters` keeps it so; `aiida_vasp.common.parameters_validator` accepts upper-case keys only on an unstored `Dict`
(it lower-cases it with a warning) and raises on a stored one, and a WorkGraph stores the inputs of its tasks. This is not
specific to the vibrations: it blocks `build_surface_workgraph` and `build_relax_static_workgraph` as well, with the
`aiida-vasp>=5,<6` that `setup.py` asks for. The unit tests only build graphs (`submit=False`) and never validate the inputs
of the work chain, so they could not see it.

Fix (commit `e9be44d`): `psteros.backends.vasp_workchain.PsterosVaspWorkChain`, `VaspWorkChain` with a `parameters` validator
that applies aiida-vasp's own validator to a lower-cased copy; `add_vasp_task` runs it. The recipes, the stored
`parameters` (still upper case), the task names, links and graph outputs are unchanged, and no existing test was edited
(AGENTS.md). Six new tests in `tests/unit/test_vasp_workchain.py`: aiida-vasp refuses a stored upper-case `Dict`, the new
class accepts it without rewriting the node, bad content is still rejected, the ports, outputs and exit codes equal those of
`VaspWorkChain`, and a built graph runs the new class with unchanged parameters. Changelog entry added. Consequence: psteros
must now be installed in the environment of the AiiDA daemon for VASP graphs as well (until now only QE relaxations and the
vibration calcfunctions needed it).

Alternative that was not taken because it changes what existing tests assert: lower-casing the keys in `vasp_parameters`.

### 8.2 Wrong `mpirun` in the copied computer registration (not a psteros bug)

I copied the `lovelace` computer of profile `presto`, whose `mpirun` is `/opt/pub/openmpi/5.0.6/gcc/12.2.0/bin/mpirun`.
`VASP-6.5.1@lovelace` is an Intel build (`LinuxIFC`, run after `module load intel/2023.2.1`), so OpenMPI's `mpirun -np 128`
started 128 independent single-rank copies in the same folder: VASP printed `running 1 mpi-ranks` (and a warning that
`NCORE=16` was overwritten by `NCORE=1`), the run took 551 s of wall time for 64 s CPU, and the retrieved `vasprun.xml`/`OUTCAR`
could not be parsed (exit 1002, handler `misc` not found, work chain exit 500). Jobs of the same code in profile
`tessera_photocatalysis_industry_metallicity`, whose computer uses the plain `mpirun -np`, finished with exit 0 and
`running 128 mpi-ranks`. Fix: `mpirun -np {tot_num_mpiprocs}` on my computer. Cost: one 9 min job on a par128 node.

### 8.3 `kpoints_spacing=0.3` is a Gamma-only mesh (plan and example setting)

`VaspCalculationConfig.kpoints_spacing` is "in 1/A, with 2 pi, as in aiida-vasp", and `VaspWorkChain` calls
`set_kpoints_mesh_from_density(spacing * 2 pi)`, so `0.3` means 1.88 1/A between k-points. The meshes actually written are
alpha-Sn (6.489 A cube) 1x1x1, O2 box 1x1x1 (intended), rutile (4.737, 4.737, 3.186) 1x1x2 and, for the slabs (3.19 x 6.70 x 20.1 A),
2x1x1 (ceil of 2 pi/(L * 0.3 * 2 pi)). The default of the dataclass (0.20) is also Gamma-only for these cells. The alpha-Sn result above is the
consequence. `examples/vasp_surface_phase_diagram/campaign.py`, which sets 0.3 and calls it "coarse k-points", has the
same problem, and so has the plan. A value of 0.03 to 0.05 gives meshes of roughly 4 to 7 points per direction.

### 8.4 Operations

* **aiida-workgraph 0.8.1 versus 0.9.0.** Both build all graphs of the test, including the 510-task vibrations graph (so
  `harmonic_modes` accepts the variable `retrieved` inputs on 0.8.1, plan step 0.2); the unit tests pass on both (154).
  The test was run on 0.9.0 (the version the branch was developed with) through a `PYTHONPATH` overlay, because the shared
  environment was not to be modified.
* **Host restarts.** The computer running the daemon was restarted during the vibrations graph (8 Oct about 21:50, recovered 9 Oct 18:12; and
  10 Oct 12:57, recovered 14:22), which killed the daemon. After the first restart the continuation tasks of the two active processes had
  to be recreated with `verdi process repair` (daemon stopped, then restarted); after the second, RabbitMQ had kept them. Nothing
  was lost: no job was in flight either time, and the graph continued from its checkpoint. After a restart the worker needs
  about 30 min of CPU time to load the checkpoint of the 510-task graph before the next job starts.
* **`max_concurrent_jobs=1` was exceeded once**: right after the first repair two displaced calculations (3917 and 3922) ran at
  the same time; all others ran alone.
* **Submitting the 510-task graph takes about 5 min** (all input nodes are stored at submission).
* **Queue.** After the first day the `par128` queue was quiet (nearly no wait) for the 6 `refs` and 4 `slabs` jobs, but a
  displaced job waited 5.2, 2.3, 2.2, 1.7 and 1.2 h in the queue overnight; the typical displaced job takes 2.8 min from creation to
  retrieval. Throughput of the vibrations graph: 5 to 15 min per job depending on the queue.

### 8.5 One mode task skipped after a host restart (WorkGraph engine state, not VASP)

After the first host restart the engine marked the task `alpha_sn_vib_s0_x_minus_vasp` FAILED although its work chain (PK 3499)
and VASP calculation (3917) finished with exit status 0, and therefore skipped `alpha_sn_vibrations`. The work chain had been
created the evening before and never started when the daemon died; I revived it with `verdi process repair`, and the
engine also launched the next displacement (3820, calc 3922) which ran at the same time. The 48 alpha-Sn calculations are all
fine, only the task that turns their forces into modes did not run. `verdi process list` and a watch on the VASP calculations
do not show this: only the task states of the WorkGraph do (`WorkGraph.load(pk).tasks`), so the watcher now checks them too.

Workaround: `rerun_alpha_sn_modes.py` builds a one-task graph (PK 4962) that runs the same `harmonic_modes` calcfunction with
the same structure (PK 989), the same settings (PK 2208) and the retrieved folders of the 48 calculations; `analyse.py --rerun-pk 4962`
reads that result for alpha-Sn. Result: 21 modes (3 x 8 - 3), 38.1 to 173.9 cm-1, none imaginary, ZPE 0.1611 eV per 8-atom cell;
recomputed from the vasprun.xml files with pymatgen they agree to 0.00001 cm-1.

Whether this is a psteros or an aiida-workgraph matter is open: psteros wires the tasks at build time as AGENTS.md asks and has
no state of its own. For the merge it means: a restart of the daemon host during a long vibrations graph needs the recovery above.

### 8.6 The O2 `IBRION=5` cross-check job failed once (VASP, not psteros)

Graph 6697: with `NSW=1`, `IBRION=5`, `NCORE=16` VASP stops with "VASP internal routines have requested a change of the k-point
set. Unfortunately, this is only possible if NPAR=number of nodes. Please remove the tag NPAR". The displaced atoms lower the
symmetry and VASP wants a new k-set, which it cannot do with the NPAR that `NCORE` implies. The work chain reported exit 700
and restarted three times as an unfinished relaxation before giving up (exit 401). Fix in the recipe of `run_test.py ibrion5`:
`ISYM=0`; graph 6729 finished with exit 0.
The 252 displaced static calculations (`IBRION=-1`) ran without it.

### 8.7 `origin/main` moved during the test, and its `AGENTS.md` gained an authorship rule

* Between 8 Oct (first merge check against `ff0cd3e`) and 12 Oct main was rebased and received a parallel line of work on the same base
  `6d8e669`, 22 commits up to `f7bdfd7`: reference thermochemistry (`psteros/thermochemistry.py` with `HarmonicSolid` and
  `IdealGasMolecule`, `references.py`, and `blocks.py` with `Relax`, `Static`, `Vibrations`), a rewrite of `psteros/backends/vasp.py`
  (INCAR keys lower-cased, `kpoints_spacing / (2 pi)` passed to aiida-vasp, `ExecutionPolicy.extra_options`, e.g.
  `{"import_sys_environment": False}`), a campaign `examples/vasp_sno2_phase_diagram_vibrations/` with its own log and report, and
  an `AGENTS.md` with a section on fixing bugs and one on authorship.
* Two of the findings above were fixed there independently: the upper-case INCAR (8.1) and the k-point unit (8.3, commit
  `044fa29`: psteros values mean A^-1 with the 2 pi of VASP's `KSPACING`; with that convention the plan's 0.3 is a sensible mesh).
  `import_sys_environment` has its own field there (`extra_options`), which is the follow-up the plan suggested. The branch under
  test has none of these.
* The public names do not collide (`Vibrations`, `HarmonicSolid` on main; `VibrationsConfig`, `HarmonicVibrations`,
  `build_vibrations_workgraph` here), but the two lines implement overlapping functions (harmonic solids, vibrations of references).
* Authorship rule of `main`'s `AGENTS.md` (no agent as author, committer or co-author): the commits made on the feature branch during this test carried a trailer that the
  rule forbids. The merge into `main` is a single squash commit by the maintainer, so none of them reaches `main` (section 9b).

## 9. Merge recommendation

**Not ready to merge into `main` as it stands.** The vibrational feature itself works on a real cluster; what stops the merge
is the state of `main` and one open physical question.

What the test established:

* The whole chain runs on Lovelace (`par128`, 128 ranks, one job at a time): `build_relax_static_workgraph` for the references and
  the two terminations, `build_vibrations_workgraph` for 252 displaced calculations, `read_vibrations`, `solid_free_energy_ev`,
  `molecule_reference_energy_ev` and the phase-diagram functions. Checks 1, 3 to 9 pass; the displaced inputs are right in all 252 cases;
  the frequencies equal an independent recomputation to 1e-5 cm-1 or better; the O2 stretch agrees with VASP's own `IBRION=5` to 0.03 cm-1;
  the 0 K results without the option are byte-identical to the code before the feature.
* The numbers are physically sensible where they can be compared: rutile 4.765/3.222 A, alpha-Sn 6.647 A, O2 2.0 muB and 1.234 A,
  Delta H_f -4.97 eV, O2 stretch 1563 cm-1, no imaginary mode in five structures.

What has to happen before the merge:

1. **Integrate with the current `main` (check 2).** 9 files conflict (22 hunks), among them both rewrites of `psteros/backends/vasp.py`
   and `psteros/config.py`, so this is a design decision, not a mechanical merge: main's adapter lower-cases the INCAR keys, converts
   the k-point unit, has `extra_options`, and `ExecutionPolicy` there still has machine defaults that this branch removed; this branch's
   adapter has `fixed_sites`, `kinds_structure` (needed by the slab vibrations) and the work chain of 8.1. Also decide how `HarmonicVibrations`/
   `build_vibrations_workgraph` relate to `HarmonicSolid`/`Vibrations` of main (two implementations of harmonic solids). After the
   integration the unit tests and this campaign (or a shortened one) must be run again.
2. **Take main's k-point fix (8.3).** On this branch `kpoints_spacing` reaches aiida-vasp unconverted, so every recipe value of the
   documentation and of `examples/vasp_surface_phase_diagram/campaign.py` (0.3) is a Gamma-only mesh. Main's conversion fixes it; the `results change` entry of
   its changelog should say that earlier VASP calculations with the old meshes have to be repeated.
   Once main's lower-casing is adopted, `PsterosVaspWorkChain` (8.1) is not needed; its tests (a stored upper-case `Dict` is accepted) can stay as regression tests of the combined adapter.
3. **Check 10 is open.** The vibrational change of gamma(slab_o) at 900 K is -25.6 meV/A2, not "a few", with a bulk reference that has only Gamma and Z
   (section 6). Either run the bulk with a larger supercell (2x2x2: 288 jobs, about a day on an idle queue) and show that gamma_vib converges, or state
   in `docs/source/vibrations.rst` that a bulk supercell smaller than that gives an unreliable surface term and show the size of the effect.
4. **Authorship.** Resolved by merging as one squash commit (section 9b).

Smaller points, not blockers: the documentation should say that a long vibrations graph needs the AiiDA daemon host to stay up and how to recover
(`verdi daemon stop`, `verdi process repair`, `verdi daemon start`; check `WorkGraph.load(pk).tasks` states afterwards, 8.4, 8.5), that psteros
must be installed in the environment of the daemon for VASP graphs too, and that `max_concurrent_jobs=1` makes a 252-job graph take from 1 to
8 days depending on the queue (7).

## 10. How to reproduce

```bash
export PYTHONPATH=<aiida-workgraph 0.9.0 overlay>:<worktree>      # daemon and client
python run_test.py refs --submit                                    # 6 jobs, then
python run_test.py slabs --refs-pk <REFS> --submit                  # 4 jobs, then
python run_test.py vibrations --refs-pk <REFS> --slabs-pk <SLABS> --submit    # 252 jobs
python run_test.py ibrion5 --refs-pk <REFS> --submit                # optional, 1 job
python rerun_alpha_sn_modes.py --vib-pk <VIB> --submit              # only if a mode task was skipped (8.5)
python analyse.py --refs-pk <REFS> --slabs-pk <SLABS> --vib-pk <VIB> [--rerun-pk <RERUN>] [--ibrion5-pk <I>] --legacy-src <checkout of 8b74134>
```

PKs of this run: refs 977, slabs 1077, vibrations 3400, alpha-Sn rerun 4962, ibrion5 6729 (profile `psteros_vibrations_lovelace`).

## 9b. Merged into `main`

After this report the feature was merged into `main` as one squash commit (the feature branch keeps its history). The integration, done against `origin/main` of 13 Oct (`a1e8f9c`, which by then already
contained the charge-neutral-terminations work this branch was built on, so 5 files conflicted instead of 9):

* `AGENTS.md`: `main`'s version; `docs/source/index.rst`, `psteros/__init__.py`, `docs/source/phase-diagram.rst`: both sides kept (no name collides: `Vibrations`, `HarmonicSolid` on main;
  `VibrationsConfig`, `HarmonicVibrations`, `build_vibrations_workgraph` here).
* `psteros/backends/vasp.py`: **`main`'s adapter, unchanged.** It already lower-cases the INCAR tags and converts the k-point unit (items 1 and 2 of section 9), so `PsterosVaspWorkChain`
  and its tests were left out; the change to `psteros/` is purely additive (4 files).
* Two tests of this feature asserted upper-case INCAR keys in the stored parameters, which pinned the bug that `main` had fixed; they now compare the lower-case keys
  (`test_vasp_vibrations_workgraph`, `test_example_vibrations_graph`).
* The five 10,000-line ternary grid CSVs (`ternary_*.csv`) are not in `main`; `analyse_ternary.py` regenerates them. The validation scripts use `kpoints_spacing = 0.25` (the value used here, 0.04 in aiida-vasp units).
* Still open from section 9: item 3 (convergence of the vibrational change of gamma with the bulk supercell) and the documentation notes.
