# RuO2(110) with a fine-tuned MACE model: results on Lovelace

Plan, gates and decisions: [`docs/plans/MACE_SURFACE_PLAN.md`](../../docs/plans/MACE_SURFACE_PLAN.md).
AiiDA profile `psteros_mace_ruo2` (core.sqlite_dos, created 2026-10-17); computer `lovelace` copied from
`psteros_vibrations_lovelace` (scheduler `tessera.pbspro_gpu`, `import_sys_environment=False`).

**CP2K moved to Phileas (2026-10-17, user: "the C1 convergence might take forever in lovelace ... launch 2 at a time,
max").** On Lovelace `par128`, C1 had 3 of 42 jobs done after 80 min; Phileas `standard` had idle nodes and no
pending job. Every CP2K energy of the study now comes from `CP2K-2026.2@phileas` (module `cp2k/2026.2`, account
m26257, at most 2 CP2K jobs at a time, `cluster.py`); the Lovelace CP2K runs below are kept as a record only. The
MACE jobs stay on the Lovelace GPUs.

| Code | What |
|---|---|
| `CP2K-2026.2@phileas` | module `cp2k/2026.2` (gcc 15.3, MKL, OpenMPI 4.1.9), `srun`; the prepend checks the executable SHA-256 and runs `cp2k.psmp --check` on the input; copied from profile `ag3po4_h2o_o2_mace_active_learning_stoichiometric`, work directory `/scratch/users/w42008@sorbonne-universite.fr/aiida_psteros_mace_ruo2` |
| `CP2K-2024.3@lovelace` (record only) | site module `cp2k/2024.3-gnu-openmpi-mkl`, `/opt/pub/cp2k/2024.3/gcc/9.4.0/bin/cp2k.psmp`, OpenMPI `mpirun -np`, `OMP_NUM_THREADS=1` |
| `mace-python-cu118@lovelace` | Ag3PO4 MACE venv (`/home/lovelace/proj/proj962/dorinitt/envs/ag3po4_mace_cu118`, mace-torch 0.3.16, torch 2.5.1+cu118), read-only, `core.shell` |

The user's own CP2K 2026.x registrations in `presto` point to binaries that do not exist on Lovelace (checked
2026-10-17); they are not used.

## C0 on Lovelace, CP2K 2024.3 (2026-10-17, `testes`, 4 MPI ranks; record only)

`c0_smoke.py`: bulk RuO2 (6 atoms, experimental lattice) with one O moved 0.05 A along x, and ±0.01 A more;
DZVP-MOLOPT-PBE-GTH (UZH), 500 Ry, k-spacing 0.30 (5×5×7), Fermi-Dirac 300 K. WorkGraph 24, all three jobs
exit 0, about 30 s each (`c0_graph_lovelace_500.json`, `c0_check_lovelace_500.json`).

| Check | Result |
|---|---|
| One CP2K run per job | yes: "Total number of message passing processes 4", CP2K version 2024.3 |
| Energy, force file, DOS file read by the psteros readers | yes |
| Bulk RuO2 metallic (`require_metallic`) | yes: 12 partly filled DOS bins |
| CP2K "Total energy" includes the electronic entropic energy | yes: the printed terms sum to it within 4e-13 Ha; it is the free energy (−TS = −0.37 mHa per cell) |
| −(E₊ − E₋)/(2δ) against the printed force on the moved O | −0.8738 against −0.8846 eV/A: 0.011 eV/A (1.2 %), just above the 0.01 eV/A the script asked for |

The free energy is the energy of the labels, so the 1.2 % gap was attributed to the real-space grid at 500 Ry. The
same check at **800 Ry** (WorkGraph 61, `c0_check_lovelace_800.json`): −0.8085 against −0.8104 eV/A, **0.0019 eV/A (0.23 %)**.
CP2K's energy is the one its forces derive from; the remaining gap at 500 Ry is the grid.

The force itself moved from −0.885 (500 Ry) to −0.810 eV/A (800 Ry), 9 %: forces converge more slowly with the cutoff
than the energy differences of the plan's C1 criteria. C1 therefore also reports the largest force change between
levels (not a gate of the plan).

The Lovelace C1 graph (WorkGraph 335, `c1_graph_lovelace.json`) was stopped when CP2K moved to Phileas.

## C0 on Phileas, CP2K 2026.2 (2026-10-17, `standard`, 8 MPI ranks, 800 Ry)

Same three calculations as on Lovelace, at 800 Ry (WorkGraph 388, `c0_graph_phileas_800.json`,
`c0_check_phileas_800.json`); about 40 s each, started at once.

| Check | Result |
|---|---|
| One CP2K run per job | yes: 8 message-passing processes, CP2K version 2026.2 |
| Readers | yes; the 2026.2 DOS file writes `# E(Fermi) = ... a.u.` on its own line, which the parser now reads (Fermi energy 10.37 eV) |
| Bulk RuO2 metallic | yes: 13 partly filled DOS bins |
| Energy consistent with the forces | −0.7805 (finite difference) against −0.7824 eV/A printed: 0.0019 eV/A |

The force on the moved O differs between the two CP2K builds at the same settings (−0.810 on Lovelace 2024.3,
−0.782 eV/A here), as do the total energies (the data files of the two versions differ): one build for the whole
study is required, which is why every CP2K step runs on Phileas.

## C1: convergence on Phileas (WorkGraph 673, 2026-10-17)

39 CP2K calculations, 2 at a time on one 96-core node each (37 statics of 13 s to 7 min, 2 cell optimisations).
Experimental lattices, unrelaxed slabs; `c1_check.json`. Changes against the level compared with (base: DZVP,
600 Ry, k-spacing 0.20, 300 K, 22 A; t150 is compared with k015):

| Change | ΔH_f (eV/f.u.) | Δ[γ(o_bridge) − γ(o_cus)] | Δ[γ(o_bridge) − γ(metal)] (meV/A²) | Largest force change (eV/A) | Criteria |
|---|---:|---:|---:|---:|---|
| 600 → 400 Ry | +0.013 | −0.09 | +0.25 | 0.200 | pass |
| 600 → 800 Ry | +0.007 | −0.20 | +0.19 | 0.073 | pass |
| DZVP → TZV2P | +0.001 | +1.46 | −0.39 | 0.094 | pass |
| k-spacing 0.20 → 0.15 | +0.001 | −0.02 | +0.20 | 0.040 | pass |
| 300 → 150 K (at 0.15) | +0.001 | −0.04 | +0.06 | 0.031 | pass |
| vacuum 22 → 30 A | 0.000 | −0.02 | +0.04 | — | pass |

Baseline values: ΔH_f = −3.079 eV/f.u. (Reuter and Scheffler −3.16, experiment −3.19); unrelaxed γ at Δμ_O = 0:
o_bridge 83.5, o_cus 46.4, metal 204.6 meV/A² (the literature 49 meV/A² is for relaxed slabs; G0 uses R1).

Bulk CELL_OPT (analytical stress with k-points): at 800 Ry it converged in 9 min, a = 4.536 A, c = 3.136 A; at
600 Ry it stalled after two optimisation steps until the 4 h walltime, and was stopped with the graph.

**Settings chosen (`settings.json`): DZVP-MOLOPT-PBE-GTH (UZH), 800 Ry, REL_CUTOFF 60, k-spacing 0.20, Fermi-Dirac
300 K, 22 A vacuum.** Every level passes at 600 Ry; 800 Ry because the forces still move by 0.073 eV/A from 600 to
800 Ry at no extra cost and because the 600 Ry cell optimisation stalled. TZV2P is not needed by the criteria
(+1.46 meV/A² ≤ 2) and costs twice as much; its force change (0.094 eV/A) does not enter the comparison of the MACE
route with CP2K, which share one recipe.

## C2 and C3 (from 2026-10-17)

References (WorkGraph 919, 39 min, all exit 0): bulk RuO2 and hcp Ru CELL_OPT -> static, O2 GEO_OPT -> static.

| Quantity | CP2K (this work) | Literature |
|---|---|---|
| RuO2 a, c, u | 4.536 A, 3.136 A, 0.3061 | experiment 4.492 A, 3.107 A, 0.306 |
| hcp Ru a, c | 2.729 A, 4.311 A | experiment 2.706 A, 4.282 A |
| O2 bond, ⟨S²⟩ | 1.237 A, 2.006 (triplet) | PBE about 1.23 A |
| ΔH_f(RuO2) | **−3.055 eV/f.u.** (O-poor limit −1.527 eV) | Reuter and Scheffler −3.16, experiment −3.19 |
| Largest force after the static | 0.0001 (bulk), 0.0011 eV/A (O2) | |
| Metallic (partly filled DOS bins) | bulk 13, Ru 18 | |

Labels (WorkGraph 1221: the six terminations as built on this lattice, bulk110, 30 rattled copies) and ground truth
(WorkGraph 1301: GEO_OPT -> static of the six terminations, central trilayer fixed) run together, one CP2K job each.
