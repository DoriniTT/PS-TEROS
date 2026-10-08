# End-to-end test of the vibrational contributions on Lovelace (par128)

Branch: `features/vibrational-contributions`.
Goal: run the new feature from start to finish on a real cluster, for one small
system, and decide whether the branch is ready to merge into `main`.

This is a **functional smoke test**, not a converged study: small cutoff, coarse
k-points, a thin slab and a small Γ-only supercell. It checks that every piece
works on a real machine and gives physically sensible numbers. It does not
produce publishable values.

## Fixed choices

| Item | Value |
|---|---|
| Cluster / queue | Lovelace, `par128` (1 node, 128 cores, 128 MPI ranks) |
| Concurrency | `max_concurrent_jobs=1` in every graph; graphs run one after the other |
| Engine | VASP (CPU build), through the psteros public API only |
| System | SnO₂(110): bulk rutile, α-Sn, triplet O₂, and **two** terminations: `o` (Sn₆O₁₂, stoichiometric) and `sn2o` (Sn₆O₈, reduced) |
| Slab | `psteros.sno2_110_slab(termination=..., triple_layers=3, vacuum_angstrom=10.0, a=..., c=...)` on the relaxed bulk lattice; central triple layer fixed with `psteros.central_sites(slab, half_width=1.5)` (one Sn₂O₄ bulk cell) |
| Electronic settings | as `examples/vasp_surface_phase_diagram/campaign.py` (`ENCUT=400`, `PREC=Accurate`, `ISMEAR=0`, `SIGMA=0.05`, `LREAL=False`, `kpoints_spacing=0.3`, POTCARs `Sn_d`, `O`), plus `NCORE` suited to 128 ranks (start with `NCORE=16`) |
| Relaxation | `IBRION=2`, `NSW=100`, `EDIFF=1e-6`, `EDIFFG=-0.01` (tighter than the example, to limit imaginary modes); `ISIF=3` for the bulk and α-Sn |
| Static and displaced calculations | `IBRION=-1`, `NSW=0`, `EDIFF=1e-7` for the displaced calculations |
| Vibrations | `VibrationsConfig(displacement_angstrom=0.01, supercells={"sno2_bulk": (1, 1, 2)}, molecules=("o2",))`; α-Sn uses its 8-atom cell; slabs displace only their free sites |
| Temperatures analysed | 0 K (ZPE only), 300, 600, 900 K |

### Expected number of VASP jobs (all serial)

| Graph | Jobs |
|---|---:|
| refs (bulk, α-Sn, O₂: relax + static) | 6 |
| slabs (2 terminations: relax + static) | 4 |
| vibrations: `slab_o` 12 free sites × 6 | 72 |
| vibrations: `slab_sn2o` 8 free sites × 6 | 48 |
| vibrations: bulk 1×1×2 supercell, 12 sites × 6 | 72 |
| vibrations: α-Sn, 8 sites × 6 | 48 |
| vibrations: O₂, 2 sites × 6 | 12 |
| optional cross-check: O₂ with VASP `IBRION=5` | 1 |
| **Total** | **≈ 263** |

With one job at a time, the queue wait per job dominates. Measure it on the first
graph and report an estimate before submitting the vibrations graph.

## Lovelace settings to verify first

From the kiln cluster notes (`kiln/clusters.yaml`, `docs/clusters.md`):

- **Code:** the CPU VASP code is registered only in some profiles. Check with
  `verdi code list` and pick `VASP-6.5.1@lovelace` or `VASP-6.6.0-HDF5-VTST214@lovelace`,
  whichever exists. Do not use the GPU or `_gamma` builds.
- **Resources:** `{"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128}`.
- **Queue:** `queue="par128"` (the scheduler plugin writes `#PBS -q`) and
  `custom_scheduler_commands="#PBS -j oe"`. Do not put the queue in both places.
- **`import_sys_environment=False` is required** (Lmod functions break the
  `qstat` parser). `ExecutionPolicy` has no field for it. Pass it per label with
  `CalculationOverride(metadata={"import_sys_environment": False}, ...)` in every
  recipe's `role_overrides`, which the backends merge into the job options. Record
  this in `RESULTS.md` as a finding. An additive `ExecutionPolicy` field could be a
  follow-up, but don't add it during this test unless the run cannot work without it.
- **POTCAR family:** find the uploaded family name in the profile; it must contain `Sn_d` and `O`.
- **Walltime:** set `max_wallclock_seconds` explicitly, e.g. 2 h for relaxations and
  1 h for static and displaced calculations.

## Steps

### 0. Environment (no cluster jobs)

1. Check out `features/vibrational-contributions`. Install it in the **daemon's**
   environment (`pip install -e .`) and run `verdi daemon restart --reset`. The
   vibration calcfunctions run on the daemon.
2. Record the versions of aiida-core, aiida-workgraph, aiida-vasp and psteros. The
   branch was developed with aiida-workgraph 0.9. With 0.8, check that the
   `harmonic_modes` task accepts its variable `retrieved` inputs.
3. Run `python -m pytest tests/unit -q`. Everything must pass, including
   `tests/unit/test_vibrations.py`. Stop and report if it doesn't.
4. **Merge check:** in a scratch worktree, run `git merge --no-commit --no-ff origin/main`
   and check for conflicts. Then run `pytest tests/unit` on the merged tree. Note that
   this branch is built on `feature/charge-neutral-terminations`, so merging it brings
   all of v2 into `main`.

### 1. Driver scripts in this folder

Write the scripts here (`validation/vibrations_lovelace/`) using only the public
`psteros` API, following `docs/source/vibrations.rst`:

- `run_test.py refs|slabs|vibrations|ibrion5 [--submit]`: builds each graph with the
  Lovelace settings above. Without `--submit`, it builds the graph and prints the
  task names and the VASP job count, matching the table above.
  - `refs` and `slabs` use `build_relax_static_workgraph`.
  - `vibrations` uses `build_vibrations_workgraph` on the relaxed structures read with
    `read_vasp_results`, with the same fixed sites that were used in the relaxation.
    Rebuild them from the unrelaxed slabs on the same lattice, as
    `examples/vasp_surface_phase_diagram/vibrations.py` does.
- `analyse.py`: reads the three graphs and writes the following to `results/`:
  - the 0 K diagram from total energies, which must equal what the pipeline gave
    before this feature;
  - the diagrams at 0 K (ZPE only), 300, 600 and 900 K, using
    `solid_free_energy_ev(..., bulk=...)` and `molecule_reference_energy_ev`;
  - a table of γ per termination at Δμ_O = 0 and of the transition Δμ_O at each T;
  - all frequencies per label as CSV.

Keep each PK in `RESULTS.md` as soon as it exists.

### 2. Run

1. Submit `refs` and wait until it finishes with exit status 0 for every child.
   Check the relaxed lattice (rutile a ≈ 4.8, c ≈ 3.2 Å with PBE) and the O₂ magnetisation (≈ 2 μB).
2. Submit `slabs` (it reads the lattice from `refs`) and wait.
3. Build `vibrations` without `--submit`, check the job count, then submit and wait.
   After the first few displaced jobs, confirm in the job files that:
   - each runs `NSW=0`;
   - the INCAR has no selective dynamics;
   - O₂ keeps `ISPIN=2`;
   - `import_sys_environment` took effect.
4. Optional: submit `ibrion5`, a single VASP `IBRION=5`, `NFREE=2`, `POTIM=0.01` job on
   the relaxed O₂, as an independent frequency reference.

Do not resubmit a failed job blindly. Find the cause (VASP error, scheduler,
parser) and record it. If a psteros bug is found, fix it on this branch following
`AGENTS.md`: additive changes only, plus a test, then commit and push.

### 3. Checks (pass/fail, recorded in `RESULTS.md`)

| # | Check | Pass when |
|---|---|---|
| 1 | Unit tests on the local stack | all pass |
| 2 | Merge with `main` | no conflicts; unit tests pass on the merged tree |
| 3 | All graphs finish | every child process finished with exit status 0 |
| 4 | Displaced calculations | each displaced structure differs from the relaxed one in exactly one coordinate of one site, by ±0.01 Å; no selective dynamics; INCAR equals the vibrations static recipe |
| 5 | Mode counts | `slab_o` 36, `slab_sn2o` 24 (3 × free sites); bulk 33 (3 × 12 − 3); α-Sn 21 (3 × 8 − 3); O₂ 1 (`zero_modes` 5) |
| 6 | O₂ stretch | 1500–1650 cm⁻¹ (PBE ≈ 1550–1600); within ~10 cm⁻¹ of the `IBRION=5` job if it was run |
| 7 | Solids | bulk SnO₂ modes within 0–800 cm⁻¹; imaginary modes absent or small (< 50i cm⁻¹). Record any, and show that `read_vibrations` raises by default and works with `imaginary_modes="drop"` |
| 8 | Independent recomputation | for one slab, recompute the frequencies from the retrieved `vasprun.xml` files with an independent parser (ASE or pymatgen) and `psteros.harmonic_vibrations_from_forces`; they must equal the graph output to < 0.1 cm⁻¹ |
| 9 | Backward compatibility | the 0 K diagram from total energies equals the one from the unmodified functions (and from `examples/vasp_surface_phase_diagram/phase_diagram.py` logic without `--vib-pk`) |
| 10 | Physical sanity | vibrational change of γ for `slab_o` of a few meV/Å² at most up to 900 K; transition Δμ_O shifts by ≲ 0.2 eV; the O-poor limit moves with ΔG_f(T) |

## Deliverables in this folder

- `run_test.py`, `analyse.py`
- `results/`: figures (`.png`) and CSV files for each temperature, and the frequency tables
- `RESULTS.md`: versions, PKs, the code, queue and option settings actually used,
  wall and queue times, the pass/fail table above, problems found with their causes
  and fixes, and a final **merge recommendation**: ready, ready after the listed
  fixes, or not ready, with reasons

Commit and push to `features/vibrational-contributions` only. Do not push to `main`
and do not open a pull request.
