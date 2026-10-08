# Results: vibrational contributions on Lovelace (par128)

Status: **in progress** (see the log at the end of each section).

## 1. Environment and settings actually used

| Item | Value |
|---|---|
| Branch | `features/vibrational-contributions` (base `47b1eef`) |
| Python | 3.13.13 (`~/envs/aiida`) |
| aiida-core / aiida-vasp / pymatgen / ase | 2.8.0 / 5.0.0 / 2026.3.23 / 3.26.0 |
| aiida-workgraph | 0.8.1 installed in the shared env; **0.9.0 used for every run** through a `PYTHONPATH` overlay (see 4.1) |
| psteros | `__version__` 2.0.0, taken from this worktree via `PYTHONPATH`; the shared env has it editable from another checkout and was not modified |
| AiiDA profile | `psteros_vibrations_lovelace`, **new**, created for this test (sqlite_dos, RabbitMQ); default profile unchanged |
| Computer | `lovelace`, a copy of the registration in profile `presto` (`core.ssh` through the cenapad ControlMaster, scheduler `tessera.pbspro_gpu`, work dir `/work/dorinitt/.aiida`) |
| Code | `VASP-6.5.1@lovelace` (`vasp_std`, copy of the `presto` registration; intel/2023.2.1 module, OpenMPI 5.0.6 mpirun) |
| POTCARs | family `PBE` (329 POTCARs copied from `presto`), mapping `Sn -> Sn_d`, `O -> O` |
| Daemon | started with `PYTHONPATH=<overlay>:<worktree>`; one worker |
| Queue / resources | `queue="par128"`, `{"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128}`, `custom_scheduler_commands="#PBS -j oe"`, `max_concurrent_jobs=1` |
| `import_sys_environment` | `False`, through `CalculationOverride(metadata=...)` for every label |
| Walltime | 7200 s for `refs` and `slabs` (relax and static share one `ExecutionPolicy`), 3600 s for `vibrations` |
| INCAR | `ENCUT=400, PREC=Accurate, EDIFF=1e-6, ISMEAR=0, SIGMA=0.05, LREAL=False, LWAVE=LCHARG=False, NCORE=16`; relax `IBRION=2, NSW=100, EDIFFG=-0.01`, `ISIF=3` for bulk and alpha-Sn, `ISIF=2` otherwise; static `IBRION=-1, NSW=0`; displaced `EDIFF=1e-7`; O2 `ISPIN=2, MAGMOM=[1,1]`, Gamma only |
| k-points | `kpoints_spacing=0.3` (1/A, with 2 pi); O2 `kpoints_distance=10` |
| Restarts | `max_iterations=3` |
| Vibrations | `VibrationsConfig(displacement_angstrom=0.01, supercells={"sno2_bulk": (1, 1, 2)}, molecules=("o2",))` |

Job script produced for the first job (`_aiidasubmit.sh`, PK 805): `#PBS -q par128`,
`#PBS -l walltime=02:00:00`, `#PBS -l select=1:mpiprocs=128:ncpus=128`, `#PBS -j oe`, and no
`#PBS -V`, so `import_sys_environment=False` took effect.

## 2. Process keys

| Graph | PK | State |
|---|---:|---|
| refs, attempt 1 (workgraph 0.8.1, psteros as on the branch) | 705 | Finished [302] after 26 s: INCAR case error (4.1) |
| refs, attempt 2 (workgraph 0.9.0, psteros as on the branch) | 751 | Finished [302] after 25 s: the same error (4.1) |
| refs, attempt 3 (with `PsterosVaspWorkChain`) | 797 | running |
| slabs | | |
| vibrations | | |
| ibrion5 | | |

## 3. Checks

(to be filled in)

## 4. Problems found

### 4.1 psteros VASP graphs cannot start with aiida-vasp 5.0.0 (INCAR keys in upper case)

(to be filled in)
