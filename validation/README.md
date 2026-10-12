# Cluster validation of `feature/charge-neutral-terminations`

End-to-end test of the polar-surface workflow (GaAs, VASP, Bohr `par40`):
relaxed bulk and references, pseudo-hydrogen, polar slabs, consistency checks,
phase diagram. Nothing here is part of the package.

## Environment (reproducible setup)

| Item | Value |
|---|---|
| AiiDA profile | `psteros_test` (own SQLite storage, created with `verdi presto`) |
| Computer | `bohr` (PBS Pro, ssh via `gate.ifi.unicamp.br`), work dir `~/.aiida_psteros_test` |
| Code | `vasp-std-cpu-h5old@bohr` = `~/bin_HDF5_FFTW/vasp_std` (Intel MPI 2021.7.1, MKL 2022, ifort HDF5 1.10.8 at `/home/sw/vasp/hdf5-1.10.8/install_old/lib`) |
| POTCARs | family `PBE` uploaded from a local `potpaw_PBE` library (has `Ga_d`, `As`, `H.75`, `H1.25`) |
| Queue | `par40`: 1 node, 40 cores, `ncpus=40` must be explicit (`num_cores_per_machine`) |
| psteros | this branch via `PYTHONPATH=$PWD` (the editable install points at another checkout) |

Daemon: `PYTHONPATH=$PWD verdi -p psteros_test daemon start 1` (restart after any code change).

## Problems found so far

1. **aiida-vasp 5.0.0 rejects upper-case INCAR tags** in a stored `Dict`
   (`ENCUT` -> "please use lower case keys"). psteros passed them upper-case, so no
   live VASP run could start. Fixed in `psteros/backends/vasp.py::vasp_parameters`.
2. `par40` rejects `select=1:mpiprocs=40` ("violates queue limits"); needs
   `ncpus=40`. Set in the driver via `resources["num_cores_per_machine"]`.
3. `~/bin_HDF5_FFTW/vasp_std` needs the ifort-built HDF5 (`install_old`), not `install`.
4. `verdi process kill` fails in this env (`plumpy.process_comms.MESSAGE_KEY`),
   a plumpy / aiida-core version mismatch unrelated to psteros.
5. **`kpoints_spacing` was not in the unit the docs say.** aiida-vasp multiplies it
   by 2*pi, so it takes 1/A *without* the 2*pi of VASP's `KSPACING`. The documented
   values (0.2 in `docs/source/vasp-workflow.rst`, 0.3 in the SnO2 campaign, 0.2/0.25 in
   the GaAs/ZnO examples) gave a **single Gamma k-point** for cells of 5 A or more:
   GaAs bulk (5.75 A) came out as `1 1 1`, the cell relaxed to 6.27 A (PBE: ~5.75 A),
   and every energy was wrong without any error. **Fixed** in psteros
   (`aiida_vasp_kpoints_spacing`, commit `80d7998` on `main`, `6ac9c6f` here): the recipe
   keeps the documented unit and the adapter converts. Jobs submitted before the fix
   (the four polar-slab graphs and the checks graph) were given the converted value
   by the driver, so their meshes are the intended ones.
6. **Relaxed structures lost their kind names** (`KeyError: 'H'` at every static step with
   pseudo-hydrogen): aiida-vasp returns the relaxed structure with element names only.
   Found only after the relaxations had run, because the unit tests built graphs but never
   ran one. **Fixed** by a `restore_vasp_kinds` calcfunction (commit `55e9366`). The finished
   relaxations were reused through AiiDA caching (`caching.enabled_for` set for
   `aiida.calculations:vasp.vasp` in the `psteros_test` profile only).
7. `examples/polar_surfaces/gaas_polar_vasp.py` still uses the single-stage
   `build_surface_workgraph` and `vasp_overrides()` without a stage, unlike the
   other examples that now use relax -> static.

## Run

```bash
export PATH=$HOME/envs/aiida/bin:$PATH PYTHONPATH=$PWD
python validation/gaas_polar.py refs --submit      # bulk, Ga, As: 3 jobs
python validation/gaas_polar.py checks --submit    # pseudo-H molecules, Eq. 7, (110) pair
python validation/gaas_polar.py slabs --submit     # (111) and (-1-1-1) slabs, 2 graphs
python validation/gaas_polar.py analyse            # report.txt, csv, png in validation/gaas/
```

## Results (GaAs, PBE, VASP 6.4.1, 400 eV, 0.15 1/A k-spacing, 9 bilayers, 2x2, all atoms relaxed)

Run on Bohr `par40`, AiiDA profile `psteros_test`; graph PKs in `validation/gaas/runs.json`.
Surface energies at the anion-rich limit (Delta mu_As = 0), pseudo-molecule method:

| Termination | gamma (meV/A^2) | Paper (cluster method) |
|---|---|---|
| (111) ideal Ga-terminated | 88.5 | - |
| **(111) Ga vacancy, 2x2** | **42.9** | **39.2** (+3.7, +9 %) |
| (-1-1-1) ideal As-terminated | 70.5 | - |
| **(-1-1-1) As vacancy, 2x2** | **59.2** | **51.2** (+8.0, +16 %) |

Other checks: Delta H_f(GaAs) = -0.698 eV (experiment -0.74 eV); GaAs a = 5.745 A; bottom check
passed for both faces (RMSD <= 0.001 A); Eq. 7 check -10.7 meV/A^2 and non-polar (110) check
+6.2 meV/A^2, against <= 6 meV/A^2 quoted for pseudo-molecules (<= 1 for clusters).
The (111) Ga vacancy is the stable termination over the whole stability window.

Not yet done: convergence (ENCUT, k-points, slab thickness, vacuum), the cluster method,
the wurtzite builder, and reproducing the numbers on a second machine.

## Results: ZnO (wurtzite builder, ideal 1x1 slabs, 9 bilayers, PBE, VASP 6.4.1, 500 eV, 0.15 1/A k-spacing)

Pseudo-molecule method, O-rich limit (Delta mu_O = 0), graph PKs in `validation/zno/runs.json`:

| Surface | gamma (meV/A^2) | Paper (cluster method) |
|---|---|---|
| ZnO(0001), Zn-terminated | 147.7 (147.692) | 147.7 (+0.0) |
| ZnO(000-1), O-terminated | 67.5 | 63.1 (+4.4, +7 %) |

The Zn-terminated face is favoured below Delta mu_O = -1.49 eV and the O-terminated face above it.
No benchmark value enters any calculation. The (0001) agreement is probably partly fortuitous:
the pseudo-molecule method carries an error of a few meV/A^2 (Eq. 7 check below).

Other checks: a = 3.286 A, c = 5.304 A (PBE literature 3.28/5.30); Delta H_f(ZnO) = -2.894 eV
(experiment -3.6 eV, the usual PBE O2-overbinding result); O2 bond 1.233 A; bottom check passed for
both faces (RMSD <= 0.001 A); Eq. 7 check -5.99 meV/A^2 and non-polar (10-10) check -0.60 meV/A^2,
inside the <= 6 meV/A^2 quoted for pseudo-molecules.

Problems found on the way (ZnO): the static step of the ideal polar slab needs `NELM=500`
(the SCF energy drifts slowly as charge reorganises across the slab; 60-170 iterations were not
enough). `AMIX=0.2, BMIX=1e-4` made the SCF diverge and is not recommended.

## Results: GaAs with the cluster method (tetrahedral clusters n = 2, 3, 4, 5; same slabs and settings)

Eq. 9 has four unknowns, so four sizes determine it exactly (as in the paper, which used n = 2, 3, 8, 9;
there is no residual to check the fit). n = 6 (175 atoms) was dropped: about 3-4 h per ionic step.

| | Pseudo-molecules | **Clusters** | Paper (clusters) |
|---|---|---|---|
| Eq. 7 check | -10.71 meV/A^2 | **+0.40** | <= 1 |
| non-polar (110) check | +6.23 meV/A^2 | **-0.58** | <= 1 |
| (111) V_Ga 2x2, gamma | 42.9 | **39.0** | 39.2 (-0.2, -1 %) |
| (-1-1-1) V_As 2x2, gamma | 59.2 | **52.0** | 51.2 (+0.8, +2 %) |
| (111) ideal Ga-terminated | 88.5 | 84.6 | - |
| (-1-1-1) ideal As-terminated | 70.5 | 63.3 | - |

The roughly 10 meV/A^2 error of the molecule method disappears with clusters, so it came from the
molecule shortcut and not from the settings. The (111) Ga vacancy is stable over the whole window.

Operational notes: caching did not reuse `cluster_Ga_n5` (its first run needed an automatic restart, so the
final calculation has a different input hash); the duplicate was cancelled by hand and the energies of
a killed graph are read from its finished tasks (`read_from_tasks` in the drivers). A graph publishes
its outputs only when it finishes normally, so `psteros.read_vasp_results` cannot read a killed or failed
graph even when most of its tasks finished: a task-level fallback in the library would help.

## ZnO convergence (cutoff and k-mesh)

The whole ZnO workflow (references, checks, slabs) repeated at 600 eV and k-spacing 0.10 1/A
(main run: 500 eV, 0.15 1/A); output in `validation/zno_conv/`, graph PKs in its `runs.json`.

| | main (500 eV, 0.15) | converged (600 eV, 0.10) | difference |
|---|---|---|---|
| gamma ZnO(0001), O-rich, meV/A^2 | 147.692 | 147.105 | -0.59 |
| gamma ZnO(000-1), O-rich, meV/A^2 | 67.512 | 67.471 | -0.04 |
| Eq. 7 check, meV/A^2 | -5.99 | -6.10 | -0.11 |
| non-polar (10-10) check, meV/A^2 | -0.60 | -0.59 | +0.01 |
| bulk a, c (A) | 3.2857, 5.3043 | 3.2872, 5.3054 | +0.04 % |

The surface energies are converged to better than 1 meV/A^2 with respect to cutoff and k-mesh,
so the difference to the paper (+4.4 meV/A^2 for (000-1), pseudo-molecules) is not a convergence
effect. Not repeated: GaAs, and the slab thickness / vacuum.

Problems met in this run: the ZnO bulk cell relaxation failed with `ZBRENT: fatal error in bracketing`
(conjugate gradient, ISIF=3; fixed with IBRION=1 for the solids, `ZNO_QN=1`), and the cold-start SCF of the
(000-1) slab static step diverged at 600 eV (fixed with ALGO=All for that step, `--algo-all-static`).
Both are VASP numerical-robustness issues, not psteros bugs; psteros has no way to retry a failed
calculation with different solver settings automatically.

## Results: ZnO with the cluster method (tetrahedral clusters n = 2, 3, 4, 5, zinc-blende analogue)

| | Pseudo-molecules | **Clusters** | Paper (clusters) |
|---|---|---|---|
| gamma ZnO(0001), O-rich, meV/A^2 | 147.69 | **147.92** | 147.7 (+0.22, +0.2 %) |
| gamma ZnO(000-1), O-rich, meV/A^2 | 67.51 | **62.22** | 63.1 (-0.88, -1.4 %) |
| Eq. 7 check, meV/A^2 | -5.99 | **-0.93** | <= 1 |
| non-polar (10-10) check, meV/A^2 | -0.60 | **-3.31** (-6.2 %) | <= 1 |

Both surface energies agree with the paper to within 1 meV/A^2 and the Eq. 7 check closes. The non-polar
check does not: -3.31 meV/A^2 with clusters, better (-0.60) with pseudo-molecules, so the molecule value
there is not evidence that the molecule method is better. It is a different test (it compares a one-sided
passivated (10-10) slab with a symmetric one and depends on the thickness of those two slabs, 16.7 A and
the same without the passivated side), and I did not investigate it further.
Fit: four sizes for the four unknowns of Eq. 9, so there is no residual.
