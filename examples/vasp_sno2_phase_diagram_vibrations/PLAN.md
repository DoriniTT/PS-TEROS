# End-to-end test: SnO2(110) with vibrational reference free energies (VASP)

Branch: `features/references-vibrational-contributions` (PS-TEROS).
Folder: `examples/vasp_sno2_phase_diagram_vibrations/` (this folder; all scripts,
results and notes of this test live here).

## Goal

Run the new feature from start to finish on a real cluster, fix every bug found,
and leave the branch in a state that can be merged into `main`:

1. reference calculations (bulk SnO2, bulk alpha-Sn, triplet O2) with
   `psteros.build_vasp_reference_workgraph`: relax -> static -> vibrations;
2. three SnO2(110) slab terminations with `psteros.build_surface_workgraph`
   (VASP path, fixed on this branch);
3. the analysis: reference free energies, the surface phase diagram, its
   Delta mu_O axis read as T and p(O2), and the stability window at T.

Everything uses the **public psteros API** only, like
`examples/qe_surface_phase_diagram/` (the QE version of the same system; use it
as the model for the scripts and README).

## Rules

- Read `AGENTS.md` first and follow it: fixes must not break existing behaviour;
  a fix that changes results gets its own commit and a `CHANGE.md` entry.
- Work only on branch `features/references-vibrational-contributions`. Commit
  small, descriptive commits and push them. Do **not** merge into `main`, and
  do not open a pull request unless asked.
- After every code change: `python -m pytest` (all green) and
  `flake8 --max-line-length=120 --ignore=E501,W503,E402,F401` on changed files.
  A bug found on the cluster gets a regression test (built without submitting).
- psteros must be installed (editable) in the Python environment of the AiiDA
  daemon, and the daemon restarted after every code change
  (`verdi daemon restart --reset`): the daemon imports
  `psteros.backends.vasp_tasks` (energy, frequencies, supercell calcfunctions).
- Keep a short log in `LOG.md` here: graph PKs, what ran, what failed, why,
  what was fixed (commit hash). Tick the boxes of this plan as you go.
- Never submit more than one job at a time (`max_concurrent_jobs=1`,
  enforced by `ExecutionPolicy`) and never submit a graph you have not first
  built with `submit=False` and inspected.

## Cluster and scheduler (Lovelace, `par128`)

From `kiln/clusters.yaml` (check against `verdi code list` / `verdi computer show lovelace`):

```python
EXECUTION = psteros.ExecutionPolicy(
    computer="lovelace",
    queue="par128",
    resources={"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128},
    max_wallclock_seconds=12 * 3600,   # check the par128 limit (qstat -Qf par128) and adjust
    max_concurrent_jobs=1,
    extra_options={"import_sys_environment": False},  # Lmod functions break the qstat parser
)
CODE = "VASP-6.5.1@lovelace"
```

`ExecutionPolicy` already writes `queue_name="par128"` (`#PBS -q par128`) and
`#PBS -j oe`. Find the AiiDA profile (`verdi profile list`) and the PBE POTCAR
family name registered in it (aiida-vasp's `verdi data vasp-potcar` commands).

## System and settings

SnO2 is the maintained psteros system, so starting structures come from psteros
(`rutile_sno2_bulk`, `alpha_sn_bulk`, `triplet_o2_cell`, `sno2_110_slab`).

Shared INCAR (one recipe for references and slabs; same ENCUT, PREC, smearing
and POTCARs everywhere):

| tag | value | note |
|---|---|---|
| ENCUT / PREC | 520 / Accurate | |
| EDIFF | 1e-7 (vibrations block: 1e-8) | |
| ISMEAR / SIGMA | 0 / 0.05 | |
| LREAL / LASPH | False / True | |
| IBRION / NSW / EDIFFG | 2 / 100 / -0.005 (refs), -0.02 (slabs) | |
| NCORE | 16 | tune if VASP complains with 128 ranks |
| LWAVE / LCHARG | False / False | |

POTCARs: `{"Sn": "Sn_d", "O": "O"}`.

**k-points:** aiida-vasp's `kpoints_spacing` is in units of 2*pi/A. Use
`kpoints_spacing=0.03` (about 0.19 1/A). The O2 molecule uses
`kpoints_distance=5.0` in its override (Gamma only).

References (`ReferenceSystem`):

| label | structure | phase | overrides | vibrations |
|---|---|---|---|---|
| `o2` | `triplet_o2_cell(cell_length=12.0)` | gas, `symmetry_number=2`, `spin=1.0` | all blocks: `ISPIN=2`, `NUPDOWN=2`, Gamma only | IBRION 5 in the box |
| `sno2` | `rutile_sno2_bulk()` (6 atoms) | solid | relax: `ISIF=3` | IBRION 6, `supercell=(2, 2, 3)` (72 atoms); if too slow, coarser k-points for the vibrations block only (`block_overrides={"vibrations": CalculationOverride(kpoints_distance=0.06)}`) |
| `sn` | `alpha_sn_bulk()` (8 atoms) | solid | relax: `ISIF=3` | IBRION 6, `supercell=(2, 2, 2)` (64 atoms) |

Blocks: `Relax()`, `Static()`, `Vibrations(incar={"ediff": 1e-8})`.

Slabs: `sno2_110_slab(termination=t, triple_layers=3, vacuum_angstrom=15.0, a=a, c=c)`
for `t` in `o`, `sno`, `sn2o`, on the **relaxed** SnO2 lattice (a, c) read from
the reference graph (see `relaxed_bulk_lattice` in
`examples/qe_surface_phase_diagram/campaign.py`). Labels `slab_o`, `slab_sno`,
`slab_sn2o`. One relaxation per slab (`ISIF=2`, all atoms free; the slabs are
symmetric). The slab energy is `misc["total_energies"]["energy_extrapolated"]`,
the same quantity `vasp_energy` uses for the references.

## Steps

### 0. Set up
- [x] `git checkout features/references-vibrational-contributions && git pull`.
- [x] `pip install -e .` in the daemon's environment, `verdi daemon restart --reset`.
- [x] `python -m pytest` passes before anything else (record the count in LOG.md).
- [x] Profile, computer, code `VASP-6.5.1@lovelace`, POTCAR family found and recorded.

### 1. Scripts (no submission yet)
- [x] `campaign.py` with subcommands `o2`, `refs` and `slabs` (argparse, like the QE
      example; `--submit` to submit, otherwise only build and print).
- [ ] `analysis.py` (see step 5).
- [x] Build every graph with `submit=False` and print, per VASP task: INCAR,
      kpoints spacing, options, settings, structure formula and atom count.
      Check against the tables above (IBRION 5/6, NSW=1, POTIM, NFREE,
      supercell sizes, Gamma-only O2, `import_sys_environment=False`).

### 2. Smoke test: O2 alone
O2 is cheap and runs every code path (relax, static, vibrations, frequency
parsing, readers).
- [x] Submit `campaign.py o2` (references `{"o2": ...}` only).
- [x] While it runs: `psteros.reference_results(pk)` works and shows states.
- [x] After it ends: all three blocks finished; `reference_thermochemistry(pk)`
      returns an `IdealGasMolecule`; 6 modes parsed, 5 near zero dropped.
- [x] Sanity: O2 bond about 1.23 A (PBE), frequency about 1550-1600 cm^-1,
      S(298.15 K, 1 bar) about 205 J/(mol K),
      `delta_mu_oxygen_ev(o2, 298.15, 1.0, include_zero_point=False)` about -0.27 eV.
- [x] Fix any bug found (see "When something fails"), then go on.

### 3. All references
- [ ] Submit `campaign.py refs` (o2, sno2, sn; the O2 results may be reused or rerun).
- [ ] Sanity, recorded in LOG.md:
  - relaxed SnO2 a about 4.83, c about 3.24 A; alpha-Sn a about 6.65 A (PBE);
  - formation energy of SnO2 from the static energies about -5 eV per formula
    unit (PBE underbinds it against the experimental -6.0 eV);
  - no imaginary bulk mode above about 20 cm^-1 after the 3 acoustic modes
    are dropped (if there are, check the relaxation and EDIFF before raising
    `imaginary_tolerance_cm1`);
  - highest SnO2 mode about 700-800 cm^-1, highest alpha-Sn mode about 200 cm^-1;
    ZPE about 0.15-0.25 eV per SnO2 formula unit and about 0.02 eV per Sn atom.

### 4. Slabs
- [ ] `campaign.py slabs --refs-pk <PK>`: slabs built on the relaxed lattice,
      built and inspected first, then submitted.
- [ ] All three relaxations finished; energies and relaxed structures readable.

### 5. Analysis (`analysis.py --refs-pk ... --slabs-pk ... --out results/sno2_110`)
- [ ] `results/reference_free_energies.csv`: `FreeEnergy.as_dict()` of each
      reference at 0, 298.15, 600, 1000 K (p = 1 bar).
- [ ] `results/sno2_110_phase_diagram.png/.csv`: `surface_phase_diagram` with
      **DFT energies** for slabs and bulk (see the warning in
      `docs/source/reference-thermochemistry.rst`: no bulk free energy inside gamma
      while slabs have none).
- [ ] Transitions printed as Delta mu_O and, through `oxygen_pressure_bar`, as
      p(O2) at 600 and 1000 K; Delta mu_O of O2 at (T, p) for a few conditions.
- [ ] O-poor limit at 0 K (DFT) and at 600 / 1000 K (free energies of SnO2 and Sn,
      bare E(O2) on the axis), compared.
- [ ] Compare the stable terminations with the QE example
      (`examples/qe_surface_phase_diagram/results/`): same ordering expected
      (O-terminated at O-rich, reduced terminations toward O-poor); numbers differ.

### 6. Finish (merge readiness)
- [ ] `README.md` here (like the QE example): what each script does, how to run
      it, the settings, the results and what they show.
- [ ] Results (`results/*.csv`, `*.png`) and `LOG.md` committed.
- [ ] Every bug fixed with a test; `CHANGE.md`, `docs/source/api.rst` and
      `docs/source/reference-thermochemistry.rst` match the final behaviour.
- [ ] Full `python -m pytest` green; flake8 clean on changed files.
- [ ] Final report in LOG.md: what was run, PKs, bugs and fixes, open issues,
      and whether the branch is ready to merge.

## When something fails

1. Find the root cause (`verdi process report <pk>`, `verdi calcjob outputcat`,
   the retrieved OUTCAR/vasprun.xml). A crash is not a flake.
2. Decide: psteros bug (fix it here), input/setting problem (change the
   campaign script), or cluster/AiiDA problem (record it, work around it in
   the script, do not patch psteros for it).
3. psteros fix: smallest change, regression test, `pytest`, daemon restart,
   commit, push, record the hash in LOG.md, rerun only what is needed.

Points to check in particular (not yet exercised on a real cluster):

- aiida-vasp with `IBRION=5/6, NSW=1`: the ionic-convergence check is disabled
  in the vibrations block; make sure no other parser or work-chain check fails
  the run, and that OUTCAR is in the retrieved folder.
- `vasp_frequencies` parses the OUTCAR printed by VASP 6.5.1 (format of the
  `f  =` / `f/i=` lines; 3N modes for the displaced cell).
- `reference_results` finds the children by link label (= task name) in this
  aiida-workgraph version; `reference_thermochemistry` scales the supercell
  modes to the static cell correctly (compare ZPE per atom of cell and supercell).
- `WorkGraph.max_number_jobs` really limits the graph to one running job.
- The INCAR fix of `build_surface_workgraph` (tags under `incar`) and the
  `Slab` conversion work on the cluster.
- Slab outputs: aiida-workgraph attaches graph outputs only when every task
  finished; read the work chain of each task otherwise (as the QE example does).
