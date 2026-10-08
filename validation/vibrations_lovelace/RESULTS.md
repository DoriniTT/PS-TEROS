# Results: vibrational contributions on Lovelace (par128)

Status: **in progress** (see the log at the end of each section).

## 1. Environment and settings actually used

| Item | Value |
|---|---|
| Branch | `features/vibrational-contributions` (base `47b1eef`) |
| Python | 3.13.13 (`~/envs/aiida`) |
| aiida-core / aiida-vasp / pymatgen / ase | 2.8.0 / 5.0.0 / 2026.3.23 / 3.26.0 |
| aiida-workgraph | 0.8.1 installed in the shared env; **0.9.0 used for every run** through a `PYTHONPATH` overlay (see 5.1) |
| psteros | `__version__` 2.0.0, taken from this worktree via `PYTHONPATH`; the shared env has it editable from another checkout and was not modified |
| AiiDA profile | `psteros_vibrations_lovelace`, **new**, created for this test (sqlite_dos, RabbitMQ); default profile unchanged |
| Computer | `lovelace`, a copy of the registration in profile `presto` (`core.ssh` through the cenapad ControlMaster, scheduler `tessera.pbspro_gpu`, work dir `/work/dorinitt/.aiida`), **with `mpirun -np {tot_num_mpiprocs}` instead of the OpenMPI 5.0.6 path of `presto`** (see 5.2) |
| Code | `VASP-6.5.1@lovelace` (`vasp_std`, copy of the `presto` registration; intel/2023.2.1 module, OpenMPI 5.0.6 mpirun) |
| POTCARs | family `PBE` (329 POTCARs copied from `presto`), mapping `Sn -> Sn_d`, `O -> O` |
| Daemon | started with `PYTHONPATH=<overlay>:<worktree>`; one worker |
| Queue / resources | `queue="par128"`, `{"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128}`, `custom_scheduler_commands="#PBS -j oe"`, `max_concurrent_jobs=1` |
| `import_sys_environment` | `False`, through `CalculationOverride(metadata=...)` for every label |
| Walltime | 7200 s for `refs` and `slabs` (relax and static share one `ExecutionPolicy`), 3600 s for `vibrations` |
| INCAR | `ENCUT=400, PREC=Accurate, EDIFF=1e-6, ISMEAR=0, SIGMA=0.05, LREAL=False, LWAVE=LCHARG=False, NCORE=16`; relax `IBRION=2, NSW=100, EDIFFG=-0.01`, `ISIF=3` for bulk and alpha-Sn, `ISIF=2` otherwise; static `IBRION=-1, NSW=0`; displaced `EDIFF=1e-7`; O2 `ISPIN=2, MAGMOM=[1,1]`, Gamma only |
| k-points | `kpoints_spacing=0.04` (aiida-vasp units, x 2 pi; the plan's 0.3 is Gamma-only, section 5.3): alpha-Sn 4x4x4, rutile 6x6x8, slabs 8x4x2 (unrelaxed cells); O2 `kpoints_distance=10` (Gamma). Graph 862 used 0.3 |
| Restarts | `max_iterations=3` |
| Vibrations | `VibrationsConfig(displacement_angstrom=0.01, supercells={"sno2_bulk": (1, 1, 2)}, molecules=("o2",))` |

Job script produced for the first job (`_aiidasubmit.sh`, PK 805): `#PBS -q par128`,
`#PBS -l walltime=02:00:00`, `#PBS -l select=1:mpiprocs=128:ncpus=128`, `#PBS -j oe`, and no
`#PBS -V`, so `import_sys_environment=False` took effect.

## 2. Process keys

| Graph | PK | State |
|---|---:|---|
| refs, attempt 1 (workgraph 0.8.1, psteros as on the branch) | 705 | Finished [302] after 26 s: INCAR case error (5.1) |
| refs, attempt 2 (workgraph 0.9.0, psteros as on the branch) | 751 | Finished [302] after 25 s: the same error (5.1) |
| refs, attempt 3 (with `PsterosVaspWorkChain`) | 797 | Killed by me after its first job (calc 805, alpha-Sn relax) failed with exit 1002: the job ran on **1 MPI rank** (5.2); calc 815 (queued second job) cancelled with 810/815 |
| refs, attempt 4 (computer `mpirun` fixed) | 862 | **Finished [0]**, all 6 children exit 0 (PBS jobs 1028184, 1028197, 1028211, 1028247, 1028256, 1028257; calcs 870, 882, 893, 905, 916, 928) |
| refs, attempt 5 (`kpoints_spacing=0.04`) | 977 | running; this is the `refs` graph used by the rest of the test (862 stays as the 0.3 record, section 5.3) |
| slabs | | |
| vibrations | | |
| ibrion5 | | |

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
PBE with Sn_d about -4.6 eV): see 5.3.

## 5. Checks

(to be filled in)

## 6. Problems found

### 5.1 psteros VASP graphs cannot start with aiida-vasp 5.0.0 (INCAR keys in upper case)

(to be filled in)

### 5.2 Wrong `mpirun` in the copied computer registration (not a psteros bug)

I copied the `lovelace` computer of profile `presto`, whose `mpirun` is `/opt/pub/openmpi/5.0.6/gcc/12.2.0/bin/mpirun`.
`VASP-6.5.1@lovelace` is an Intel build (`LinuxIFC`, run after `module load intel/2023.2.1`), so OpenMPI's `mpirun -np 128`
started 128 independent single-rank copies in the same folder: VASP printed `running 1 mpi-ranks` (and a warning that
`NCORE=16` was overwritten by `NCORE=1`), the run took 551 s of wall time for 64 s CPU, and the retrieved `vasprun.xml`/`OUTCAR`
could not be parsed (exit 1002, handler `misc` not found, work chain exit 500). Jobs of the same code in profile
`tessera_photocatalysis_industry_metallicity`, whose computer uses the plain `mpirun -np`, finished with exit 0 and
`running 128 mpi-ranks`. Fix: `mpirun -np {tot_num_mpiprocs}` on my computer. Cost: one 9 min job on a par128 node.

### 5.3 `kpoints_spacing=0.3` is a Gamma-only mesh (plan and example setting)

`VaspCalculationConfig.kpoints_spacing` is "in 1/A, with 2 pi, as in aiida-vasp", and `VaspWorkChain` calls
`set_kpoints_mesh_from_density(spacing * 2 pi)`, so `0.3` means 1.88 1/A between k-points. The meshes actually written are
alpha-Sn (6.489 A cube) 1x1x1, O2 box 1x1x1 (intended), rutile (4.737, 4.737, 3.186) 1x1x2 and, for the slabs (3.19 x 6.70 x 20.1 A),
2x1x1 (ceil of 2 pi/(L * 0.3 * 2 pi)). The default of the dataclass (0.20) is also Gamma-only for these cells. The alpha-Sn result above is the
consequence. `examples/vasp_surface_phase_diagram/campaign.py`, which sets 0.3 and calls it "coarse k-points", has the
same problem, and so has the plan. A value of 0.03 to 0.05 gives meshes of roughly 4 to 7 points per direction.
