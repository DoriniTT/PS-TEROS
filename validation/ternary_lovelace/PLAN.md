# Smoke test of the vibrational contributions for a ternary oxide: SrTiO3(001) on Lovelace

Follows `validation/vibrations_lovelace/PLAN.md` and `RESULTS.md`: same profile (`psteros_vibrations_lovelace`), code
(`VASP-6.5.1@lovelace`), queue (`par128`, 128 ranks, `max_concurrent_jobs=1`), electronic settings and
`import_sys_environment=False`. A **functional** smoke test: it asks whether the pieces compose for a ternary compound and
gives sensible numbers, not converged values.

## Question

Do `build_vibrations_workgraph`, `read_vibrations`, `solid_free_energy_ev`, `molecule_reference_energy_ev` and
`TernaryOxideReferences` / `ternary_surface_phase_diagram` work together on a real ternary oxide, with the frozen centre of a slab
counted as bulk SrTiO3?

## System

| Item | Value |
|---|---|
| Bulk | cubic perovskite SrTiO3, 5-atom cell, `ISIF=3` |
| Elements | fcc Sr (4-atom cell), hcp Ti (2-atom cell), both `ISIF=3`; O2 is **reused** from the SnO2 test (graphs 977 and 3400, same settings) |
| Competing phases | rocksalt SrO (8-atom cell), rutile TiO2 (6-atom cell), `ISIF=3` |
| Slabs | SrTiO3(001), 1x1, 5 layers, 10 A vacuum, on the relaxed bulk lattice: `slab_tio2` = TiO2/SrO/TiO2/SrO/TiO2 (Sr2Ti3O8, 13 atoms) and `slab_sro` = SrO/TiO2/SrO/TiO2/SrO (Sr3Ti2O7, 12 atoms). The frozen centre is layers 2 and 3 (one SrO + one TiO2 layer = one SrTiO3 formula unit); the other atoms are free |
| POTCARs | `Sr_sv`, `Ti_pv`, `O` (family `PBE`) |
| Settings | as in `validation/vibrations_lovelace/run_test.py` (`ENCUT=400`, `EDIFF=1e-6`, relaxation `EDIFFG=-0.01`, `NCORE=16`, `kpoints_spacing=0.04`, displaced calculations `EDIFF=1e-7`, displacement 0.01 A) |
| Gamma-only supercells | SrTiO3 1x1x2 (10 atoms), TiO2 1x1x1 (6), SrO 1x1x1 (8), Sr 1x1x1 (4), Ti 1x1x2 (4): small on purpose, so the bulk sampling is poor (see RESULTS.md of the SnO2 test, check 10) |
| Temperatures analysed | static energies, 0 K (ZPE only), 300, 600, 900 K |

## Jobs (all serial)

| Graph | Jobs |
|---|---:|
| refs: 5 solids, relax + static | 10 |
| slabs: 2 terminations, relax + static | 4 |
| vibrations: SrTiO3 10 atoms x 6 | 60 |
| vibrations: TiO2 6 x 6 | 36 |
| vibrations: SrO 8 x 6 | 48 |
| vibrations: Sr 4 x 6 | 24 |
| vibrations: Ti 4 x 6 | 24 |
| vibrations: slab_tio2, 8 free atoms x 6 | 48 |
| vibrations: slab_sro, 7 free atoms x 6 | 42 |
| **Total** | **296** |

On an idle queue about 4 min per job (about 20 h); with the busy queue of the SnO2 test about 16 min per job (about 3 days).

## Steps

1. `run_ternary.py refs [--submit]`; check the relaxed lattices (SrTiO3 a about 3.94 A, PBE), all children exit 0.
2. `run_ternary.py slabs --refs-pk R [--submit]`.
3. `run_ternary.py vibrations --refs-pk R --slabs-pk S [--submit]`; build first without `--submit` and compare the job counts with the table.
4. `analyse_ternary.py`: free-energy ternary diagrams, tables, frequencies and the checks below.

## Checks

| # | Check | Passes when |
|---|---|---|
| T1 | Graphs finish | every VASP calculation exits 0 and no WorkGraph task is FAILED or SKIPPED |
| T2 | Displaced calculations | one coordinate of one site moved by +-0.01 A, no selective dynamics, `NSW=0`, INCAR = recipe |
| T3 | Mode counts | 3 x displaced sites, minus 3 for a periodic cell with every site displaced: SrTiO3 27, TiO2 15, SrO 21, Sr 9, Ti 9, slab_tio2 24, slab_sro 21 |
| T4 | Frozen centre | `solid_free_energy_ev(..., bulk=SrTiO3)` counts 1 SrTiO3 unit for both slabs and refuses a non-stoichiometric frozen region |
| T5 | Imaginary modes | none, or all small and `read_vibrations` raises by default and works with `imaginary_modes="drop"` |
| T6 | Independent recomputation | one slab recomputed from the retrieved `vasprun.xml` files with pymatgen agrees to < 0.1 cm-1 |
| T7 | Backward compatibility | the static ternary diagram equals the one of the commit before the vibrations (CSV bytes) |
| T8 | Physical sanity | SrTiO3 Delta H_f about -17 eV at the static level and moving by less than 0.3 eV up to 900 K; a stability polygon exists at every T; the stable termination at each point is the same set as the static one or the change is explained |

Deliverables: `run_ternary.py`, `analyse_ternary.py`, `RESULTS.md`, `results/` (figures, CSV, frequencies). Commit and push to
`features/vibrational-contributions` only, without agent attribution in commits (see `AGENTS.md` of `main`).
