# Log: SnO2(110) with vibrational references (VASP, Lovelace par128)

Plan: `PLAN.md`. Branch: `features/references-vibrational-contributions`.

## Environment

| Item | Value |
|---|---|
| Python / aiida-core / aiida-vasp / aiida-workgraph | 3.13.13 / 2.8.0 / 5.0.0 / 0.8.1 (`~/envs/aiida`, shared) |
| psteros | this worktree, through `PYTHONPATH` for my shell and my daemon; the editable install of the shared env points to another checkout and was not modified |
| AiiDA profile | `psteros_sno2_vibrations` (new, sqlite_dos + RabbitMQ). The default profile was not changed. `psteros_vibrations_lovelace` belongs to another session working on `features/vibrational-contributions` and was not touched |
| Computer / code | `lovelace` and `VASP-6.5.1@lovelace`, copied from profile `presto` (`verdi computer export`, `verdi code export`) |
| POTCARs | family `PBE` (329 POTCARs) exported from `presto` and imported with `verdi archive`; `Sn -> Sn_d`, `O -> O` resolve |
| Queue | `par128`: walltime max 168 h, max 6 queued / 6 running per user. The user already had other jobs there |
| Daemon | `verdi -p psteros_sno2_vibrations daemon start 1` with `PYTHONPATH=<worktree>` |
| Tests before any change | `python -m pytest`: 328 passed, 7 skipped |

Every `verdi` call uses `-p psteros_sno2_vibrations`: the default profile of this machine is somebody else's.

## Timeline

- Scripts: `campaign.py` (phases `o2`, `refs`, `slabs`), built and printed with `submit=False`; INCAR,
  k-points, options, settings, supercells (2x2x3 for SnO2 = 72 atoms, 2x2x2 for alpha-Sn = 64 atoms),
  Gamma-only O2 (`kpoints_spacing=5.0`), `IBRION` 5/6, `NSW=1`, `import_sys_environment=False` checked
  against the tables of `PLAN.md`.
- O2 smoke test, attempt 1: graph PK 353, failed in 1 minute before anything reached the cluster.
  `VaspCalculation` 362 excepted in the pre-submit step (`NotExistent: No PotcarFileData nodes found`).
  Cause: my new profile had the 329 `PotcarData` nodes of family `PBE` but not their `PotcarFileData`
  nodes (aiida-vasp keeps the POTCAR text in a second node, found by sha512). Setup problem of the
  profile, not psteros: imported the 329 `PotcarFileData` nodes from `presto` and checked that
  `PotcarData.find_file_node()` resolves for `Sn_d` and `O`. No psteros change.
- O2 smoke test, attempt 2: graph PK 714 (job `aiida-723`), killed by me after 44 min of a relaxation that should take
  seconds. The VASP output printed `running 1 mpi-ranks` 128 times: 128 independent serial copies of VASP wrote the
  same files. Cause: `vasp_std` 6.5.1 is linked with Intel MPI (`libmpi.so.12`, `mpiifort`), but the computer
  copied from `presto` has `mpirun_command = /opt/pub/openmpi/5.0.6/.../mpirun`, an OpenMPI launcher, so no rank
  talked to another. Cluster/AiiDA setup problem, not psteros. Fix in my profile only: computer `lovelace`
  `mpirun_command = mpirun -np {tot_num_mpiprocs}`; after `module load intel/2023.2.1` (code prepend text)
  `mpirun` is Intel MPI's. The job was cancelled (gone from `qstat`). Note for the user: profile
  `psteros_vibrations_lovelace` (other session) was set up with the same OpenMPI `mpirun` and may have the same problem.
- O2 smoke test, attempt 3: graph PK 747 (psteros as on the branch, MPI fixed). Relax (751 / calc 756) and
  static (766 / calc 771) finished OK; `reference_results(747)` worked while the graph was running:
  O-O bond 1.233 A, E(O2) = -9.8847 eV (sigma -> 0). Vibrations calc 785 ended with exit 700
  ("Calculation did not reach the end of execution") and the work chain restarted it (791) as if it were an
  unfinished relaxation; I killed 747 and 791 (the PBS job was removed).
  Root cause, in the retrieved OUTCAR of 785: VASP stopped itself with "VASP internal routines have requested a
  change of the k-point set. Unfortunately, this is only possible if NPAR=number of nodes. Please remove the tag
  NPAR". Displaced atoms lower the symmetry, so VASP changes its k-point set, and it cannot do that with band
  parallelisation (`NCORE = 16` of the shared INCAR). This is a psteros bug (the block did not account for it),
  and every IBRION 5/6 run with `NCORE > 1` would hit it.
  Fix: `Vibrations.defaults` now returns `{"isif": 2, "ncore": 1}`; tests in `tests/unit/test_blocks_references.py`
  and `tests/test_reference_workgraph.py`; `CHANGE.md`, `docs/source/api.rst` and
  `docs/source/reference-thermochemistry.rst` updated. 331 passed, 7 skipped.
- O2 smoke test, attempt 4: graph PK 815, after commit `a709819` (NCORE=1 for the vibrations block). AiiDA caching
  was switched on in my profile (`verdi -p psteros_sno2_vibrations config set caching.default_enabled true`) so the
  relax and static blocks, whose inputs are identical to attempt 3, were taken from the cache (calc 824, 839, marked
  with the cache check-mark); only the vibrations job ran (calc 853, 38 min on 128 ranks).
  It failed: OUTCAR "Error EDDDAV: Call to ZHEGV failed. Returncode = 137 2 256" (exit 700 again, then the work
  chain aborted, graph exit 302). The k-point-set error is gone. Cause: with NCORE = 1 on 128 ranks VASP uses
  NBANDS = 128 for the 12 electrons of O2 (6 occupied bands); the huge empty subspace is ill-conditioned.
  A second psteros bug in my first fix. Fix: for a gas the Vibrations block sets `ISYM = 0` (k-point set cannot
  change, recipe NCORE kept), for a solid `NCORE = 1` (supercell with hundreds of bands, symmetry-reduced
  displacements kept). Tests, CHANGE.md and docs updated (the tests of `a709819` were mine and for unreleased
  behaviour; no earlier test was edited). 331 passed, 7 skipped.
- O2 smoke test, attempt 5: **graph PK 880 finished OK** (exit 0), on psteros `169e291`. Relax 884/889 and static
  899/904 came from the cache (identical inputs to attempt 3, which ran for real as PK 747: calcs 756 and 771);
  vibrations ran for real (work chain 913, calc 918) and `vasp_frequencies` (925) parsed the OUTCAR.
  Checks (all within the ranges of PLAN.md step 2):
  - `reference_results(747)` worked while the graph was running (relax finished, static waiting).
  - 6 modes parsed: 1567.41, 0.0019, 0.0013, -0.00004, -21.23, -21.23 cm^-1; 1 vibration kept, 5 dropped
    (translations, rotations; the two -21 cm^-1 are the box rotations, not imaginary vibrations).
  - O-O bond 1.233 A; E_relax = -9.88472377 eV, E_static = -9.88472388 eV (sigma -> 0).
  - `reference_thermochemistry(880)` returns an `IdealGasMolecule`; ZPE = 0.0972 eV, S(298.15 K, 1 bar) =
    205.4 J/(mol K); `delta_mu_oxygen_ev(o2, 298.15, 1.0, include_zero_point=False)` = -0.272 eV
    (with ZPE -0.224 eV); p(O2) at 600 K for Delta mu_O = -1 eV is 4.6e-8 bar.
- References, attempt 1: graph PK 995 (all three references, psteros `169e291`). The O2 blocks came from the cache.
  alpha-Sn relax (calc 1050): a = 6.6516 A, E = -30.7648 eV for 8 atoms, static E from calc 1071 (energy parsed).
  alpha-Sn vibrations (calc 1081, 64 atoms, IBRION 6, NCORE 1, 128 ranks, 4 k-points, 640 bands): VASP finds 1
  degree of freedom (cubic diamond), so 2 displacements; about 40 min each (3-4 min per SCF step).
  Planning the SnO2 vibrations from that timing: 72 atoms, ~12 displacements of lower symmetry on a 4x4x4 mesh would
  need roughly a day or more, beyond the 12 h of the graph. Used the escape hatch of PLAN.md: a 2x2x2 mesh
  (`kpoints_distance=0.06`) and 24 h of walltime for the SnO2 vibrations block only (SnO2 is a wide-gap insulator and
  the cell is 9.7 A wide); relax and static keep the fine mesh. `campaign.py` updated. Graph 995 is stopped after
  the alpha-Sn vibrations and resubmitted; finished calculations come from the cache.
