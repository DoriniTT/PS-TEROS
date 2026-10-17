# RuO2(110) with a fine-tuned MACE model: results on Lovelace

Plan, gates and decisions: [`docs/plans/MACE_SURFACE_PLAN.md`](../../docs/plans/MACE_SURFACE_PLAN.md).
AiiDA profile `psteros_mace_ruo2` (core.sqlite_dos, created 2026-10-17); computer `lovelace` copied from
`psteros_vibrations_lovelace` (scheduler `tessera.pbspro_gpu`, `import_sys_environment=False`).

| Code | What |
|---|---|
| `CP2K-2024.3@lovelace` | site module `cp2k/2024.3-gnu-openmpi-mkl`, `/opt/pub/cp2k/2024.3/gcc/9.4.0/bin/cp2k.psmp`, OpenMPI `mpirun -np`, `OMP_NUM_THREADS=1` |
| `mace-python-cu118@lovelace` | Ag3PO4 MACE venv (`/home/lovelace/proj/proj962/dorinitt/envs/ag3po4_mace_cu118`, mace-torch 0.3.16, torch 2.5.1+cu118), read-only, `core.shell` |

The user's own CP2K 2026.x registrations in `presto` point to binaries that do not exist on Lovelace (checked
2026-10-17); they are not used.

## C0: CP2K smoke test (2026-10-17, `testes`, 4 MPI ranks)

`c0_smoke.py`: bulk RuO2 (6 atoms, experimental lattice) with one O moved 0.05 A along x, and ±0.01 A more;
DZVP-MOLOPT-PBE-GTH (UZH), 500 Ry, k-spacing 0.30 (5×5×7), Fermi-Dirac 300 K. WorkGraph 24, all three jobs
exit 0, about 30 s each.

| Check | Result |
|---|---|
| One CP2K run per job | yes: "Total number of message passing processes 4", CP2K version 2024.3 |
| Energy, force file, DOS file read by the psteros readers | yes |
| Bulk RuO2 metallic (`require_metallic`) | yes: 12 partly filled DOS bins |
| CP2K "Total energy" includes the electronic entropic energy | yes: the printed terms sum to it within 4e-13 Ha; it is the free energy (−TS = −0.37 mHa per cell) |
| −(E₊ − E₋)/(2δ) against the printed force on the moved O | −0.8738 against −0.8846 eV/A: 0.011 eV/A (1.2 %), just above the 0.01 eV/A the script asked for |

The free energy is the energy of the labels, so the 1.2 % gap was attributed to the real-space grid at 500 Ry. The
same check at **800 Ry** (WorkGraph 61, `c0_check_800.json`): −0.8085 against −0.8104 eV/A, **0.0019 eV/A (0.23 %)**.
CP2K's energy is the one its forces derive from; the remaining gap at 500 Ry is the grid.

The force itself moved from −0.885 (500 Ry) to −0.810 eV/A (800 Ry), 9 %: forces converge more slowly with the cutoff
than the energy differences of the plan's C1 criteria. C1 therefore also reports the largest force change between
levels (not a gate of the plan).
