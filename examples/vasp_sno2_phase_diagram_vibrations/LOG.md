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
- References, attempt 2: graph PK 1168 (psteros `3f71ef2`). O2 and alpha-Sn came from the cache (alpha-Sn
  vibrations had finished in graph 995: calc 1081, frequencies 1088). SnO2 relax (calc 1271) finished:
  a = 4.8301 A, c = 3.2434 A (PLAN: about 4.83 / 3.24), E = -37.30513 eV, max force 2 meV/A. SnO2 static
  (calc 1288) OK. SnO2 vibrations (calc 1302, 8 k-points, 384 bands, 30 s per SCF step) printed `DOF = 216`,
  `Found 1 space group operations` and `Total: 1/432` displacements: about 7 min each, 50 h, more than the 24 h.
  I killed graph 1168 (calc 1302 was cancelled, nothing left on the queue) after 36 min.
  Root cause: the POSCAR of the supercell had 24 element blocks (`Sn 2, O 4, Sn 2, O 4, ...`) because `ase`
  repeats cell by cell and `make_supercell` kept that order; VASP treats each block as a different ion type, so it
  found no symmetry at all (spglib finds P4_2/mnm with 192 operations in the same structure, at every tolerance
  from 1e-8 to 1e-2, so it was not numerical noise). alpha-Sn (one element) was not affected (DOF = 1).
  psteros bug in `make_supercell`. Fix: group the atoms by element in order of first appearance; regression
  test `test_supercell_groups_atoms_by_element_so_that_vasp_keeps_the_symmetry` (also checks P4_2/mnm with
  16 x 12 operations); `CHANGE.md`. 332 passed, 7 skipped. The daemon was restarted (it runs the calcfunction).
- References, attempt 3: **graph PK 1372 finished OK** (psteros `f0c0eb1`, 2x2x2 k-mesh and 24 h walltime for the SnO2
  vibrations). O2 (relax 1376, static 1391, vibrations 1405) and alpha-Sn (relax 1422, static 1437, vibrations 1453) came
  from the cache of earlier graphs (747/815/880 and 995); SnO2 ran for real: relax 1470, static 1485 and
  vibrations 1501 (frequencies 1513). VASP now found 16 space-group operations (D_2h), `DOF = 4`, 8 displacements,
  6 k-points. Sanity (PLAN step 3):
  - relaxed SnO2: a = 4.8301 A, c = 3.2434 A (plan: about 4.83 / 3.24); alpha-Sn: a = 6.6516 A (about 6.65).
  - Delta E_f(SnO2) from the static energies = -4.925 eV per formula unit (about -5; experiment -6.0).
    E(SnO2 cell) = -37.31192 eV, E(Sn cell of 8) = -30.77105 eV, E(O2) = -9.88472 eV.
  - Imaginary modes: none. SnO2: 216 modes, the 3 lowest are -0.14, -0.14, -0.01 cm^-1 (acoustic, dropped); the
    first optical one is 77.95 cm^-1. alpha-Sn: 192 modes, the 3 lowest are 0 cm^-1; the first is 34.8 cm^-1.
  - Highest modes: SnO2 706.5 cm^-1 (plan 700-800), alpha-Sn 173.3 cm^-1 (plan: about 200; PBE underestimates
    the experimental 200-ish, within the expected 10-15 %).
  - ZPE: SnO2 0.1985 eV per formula unit (0.15-0.25), alpha-Sn 0.0211 eV per atom (about 0.02), O2 0.0972 eV.
    ZPE per atom of the SnO2 cell 0.0662 eV, so the supercell modes scale to the static cell correctly.
- Slabs: **graph PK 1539 finished OK** (each relaxation converged in one VASP calculation): `slab_o` (work chain 1543,
  Sn6O12, E = -109.08139 eV), `slab_sn2o` (1556, Sn6O8, -82.00946 eV), `slab_sno` (1569, Sn6O10, -94.71028 eV), built on
  the relaxed lattice a = 4.8301 A, c = 3.2434 A, 3 triple layers, 15 A vacuum, cell 3.24 x 6.83 x 27.32 A.
- Analysis (`analysis.py --refs-pk 1372 --slabs-pk 1539 --out results/sno2_110`): the files of `results/`. gamma(slab_o)
  = 1.032 J/m^2 (QE example 1.033); window -2.462 <= Delta mu_O <= 0 eV (QE -2.43); slab_sn2o -> slab_o at
  Delta mu_O = -1.826 eV, i.e. p(O2) = 6.3e-22 bar at 600 K and 1.6e-8 bar at 1000 K; O-poor limit -2.462 eV (0 K, DFT),
  -2.414 eV (600 K) and -2.521 eV (1000 K) with free energies. Same stable terminations as the QE example at the
  five points of the comparison (O-rich: slab_o; O-poor end: slab_sn2o; slab_sno never stable in either).
- `main` had moved on (3 commits, incl. 80d7998 "Fix the VASP k-point spacing unit passed to aiida-vasp") and conflicted
  with the branch in `AGENTS.md`, `docs/source/api.rst` and `psteros/backends/vasp.py`. Merged `origin/main` into
  the branch (merge commit 43f69f6, no force push): main's AGENTS.md, `kpoints_spacing` in A^-1 with the 2*pi
  converted for aiida-vasp in BOTH builders (the reference blocks did not convert before), the campaign passes
  `2*pi*0.03`, `2*pi*0.06` and `5*2*pi` so that aiida-vasp receives exactly 0.03, 0.06 and 5.0 as in every
  validated run; docs, example and changelog entries (also in `docs/source/changelog.rst`) updated. Four assertions
  of `tests/test_reference_workgraph.py` (tests of this branch's own, unreleased feature) were changed to the
  converted values (value / 2*pi) and a test that both builders agree was added.
  Check after the merge: graphs 1647 (refs) and 1679 (slabs), same recipe, resubmitted: all 12 VaspCalculations came
  from the cache (inputs identical), nothing ran on the cluster, and `analysis.py` on them gives CSV files
  byte-identical to the committed ones. 337 passed, 7 skipped; flake8 clean on the changed files.

## Final report

**What ran** (profile `psteros_sno2_vibrations`, Lovelace `par128`, one job at a time, VASP 6.5.1, 128 ranks):

| Graph | PK | Result |
|---|---:|---|
| O2 alone, attempt 1 / 2 / 3 / 4 | 353 / 714 / 747 / 815 | failed: POTCAR file nodes missing / wrong MPI launcher / VASP k-point-set error / ZHEGV failure |
| O2 alone, attempt 5 | 880 | finished OK |
| references (O2, SnO2, alpha-Sn), attempt 1 / 2 | 995 / 1168 | stopped by me after alpha-Sn / after a wrong 432-displacement plan |
| references, attempt 3 | **1372** | finished OK (O2 and alpha-Sn from the cache, SnO2 ran) |
| slabs (3 terminations) | **1539** | finished OK |
| references / slabs on the merged code | 1647 / 1679 | finished OK, all from the cache |

The real VASP jobs: O2 relax 756 and static 771 (graph 747), O2 vibrations 918 (graph 880), alpha-Sn relax 1050, static
1071, vibrations 1081 (graph 995), SnO2 relax 1271, static 1288 (graph 1168), SnO2 vibrations 1506 (graph 1372), slabs
1548, 1561, 1574 (graph 1539).

**Bugs found and fixed in psteros** (each with a regression test, `CHANGE.md` and docs):

1. `a709819` + `169e291` `Vibrations` block: displaced atoms lower the symmetry, VASP must change its k-point set and
   refuses with band parallelisation ("requested a change of the k-point set ... remove NPAR", exit 700, restarted by
   aiida-vasp as an unfinished relaxation). First fix `NCORE = 1`, which failed for the molecule (128 bands for 6 occupied
   ones, "ZHEGV failed"); final: `ISYM = 0` for a gas, `NCORE = 1` for a solid.
2. `f0c0eb1` `make_supercell`: `ase` repeats cell by cell, the POSCAR had 24 element blocks, VASP saw 24 ion types and no
   symmetry (216 degrees of freedom, 432 displacements, about 50 h). Atoms are now grouped by element: 16 operations,
   4 degrees of freedom, 8 displacements.
3. `43f69f6` (merge) the reference blocks did not convert `kpoints_spacing` like the surface builder of `main`.

**Not psteros bugs, worked around in my profile only** (not in the repository): the POTCAR family needs the
`PotcarFileData` nodes as well as the `PotcarData` ones; the computer copied from `presto` used an OpenMPI `mpirun`, but
this VASP binary is linked with Intel MPI, so 128 serial copies ran (fix: `mpirun -np {tot_num_mpiprocs}` after
`module load intel/2023.2.1`). The profile `psteros_vibrations_lovelace` of another session was set up from the same
`presto` registration and may have the same launcher problem.

**Open issues**

- The SnO2 vibrations use a 2x2x2 k-mesh (the PLAN's escape hatch, for walltime); not compared with a finer mesh. alpha-Sn
  highest mode 173 cm^-1 is below the PLAN's "about 200" (expected for PBE).
- The slabs are three triple layers thick: a workflow test, not converged surface energies. The diagram uses DFT energies
  for slabs and bulk, as the docs warn; the free energies enter only the T, p reading and the stability window.
- AiiDA caching is on in my profile and the daemon of `psteros_sno2_vibrations` is still running (stop it with
  `verdi -p psteros_sno2_vibrations daemon stop`). The shared environment's editable install of psteros was not touched;
  the daemon takes psteros from this worktree through `PYTHONPATH`.
- Pushing needed `gh auth git-credential` as a one-off credential helper (no stored git credentials in this shell).

**Ready to merge into `main`?** Yes, from the code side: the branch (merge commit 43f69f6 and the commits after it)
merges cleanly into `origin/main` as of the last fetch, 337 tests pass, flake8 is clean, and the end-to-end run reproduces
itself through the cache after the merge. Not merged and no pull request opened, as instructed. Before merging, glance at
the four edited assertions in `tests/test_reference_workgraph.py` and at the unit change in the reference blocks.
