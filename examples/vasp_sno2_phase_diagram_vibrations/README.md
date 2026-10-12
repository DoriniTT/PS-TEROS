# SnO₂(110) with vibrational reference free energies (VASP)

An end-to-end psteros example on a real cluster (Lovelace, queue `par128`): VASP
reference calculations with relaxation, static energy and vibrations
(`psteros.build_vasp_reference_workgraph`), three SnO₂(110) terminations
(`psteros.build_surface_workgraph`), and the analysis that turns them into
reference free energies, γ(Δμ_O), and Δμ_O read as T and p(O₂). It is the VASP
counterpart of `examples/qe_surface_phase_diagram/`; only the public psteros API
is used.

| File | What it does |
|---|---|
| `campaign.py` | Phases `o2`, `refs`, `slabs`: builds, prints and (with `--submit`) submits a WorkGraph. |
| `analysis.py` | Reads the finished graphs and writes everything in `results/`. |
| `PLAN.md`, `LOG.md` | The test plan and the log of the run (PKs, failures, fixes, final report). |
| `results/` | The numbers and the figure of the run described below. |

## Before you start

* An AiiDA profile with a `vasp.vasp` code, a POTCAR family (`PBE`, with `Sn_d`
  and `O`) and a running daemon. The code's MPI launcher must match the one
  VASP was built with (on Lovelace VASP 6.5.1 is linked against Intel MPI).
* psteros importable by the **daemon** (editable install, or `PYTHONPATH`), and
  the daemon restarted after every psteros change: it runs the calcfunctions
  `vasp_energy`, `vasp_frequencies` and `make_supercell`.
* `campaign.py` always builds the graph and prints, for every VASP task, its
  INCAR, k-point spacing, scheduler options, settings and structure. Read that
  before adding `--submit`.

## Run it

```bash
python campaign.py o2 --submit                 # smoke test: triplet O2 alone
python campaign.py refs --submit               # O2, SnO2, alpha-Sn: relax -> static -> vibrations
python campaign.py slabs --refs-pk <PK of refs> --submit
python analysis.py --refs-pk <PK of refs> --slabs-pk <PK of slabs> --out results/sno2_110
```

`slabs` builds the slabs on the relaxed SnO₂ lattice of the refs graph. Every
graph runs one job at a time (`max_concurrent_jobs=1`).

## Settings

One INCAR for the references and the slabs: `ENCUT=520`, `PREC=Accurate`,
`EDIFF=1e-7` (`1e-8` for the vibrations block), `ISMEAR=0`, `SIGMA=0.05`,
`LREAL=False`, `LASPH=True`, `NCORE=16`, `LWAVE=LCHARG=False`; POTCARs
`Sn_d` and `O`; `kpoints_spacing = 2π·0.03 = 0.1885 Å⁻¹` (with the 2π, as VASP's
`KSPACING`; meshes: SnO₂ 8×8×11 for the relaxation and 7×7×11 for the static energy
on the relaxed cell, α-Sn 6×6×6). Relaxations: `IBRION=2`, `NSW=100`, `EDIFFG=-0.005` for the
references (tight, the frequencies are computed there) and `-0.02` for the
slabs, `ISIF=3` for the two bulks, `ISIF=2` for the slabs.

| Reference | Cell | Vibrations |
|---|---|---|
| triplet O₂ | 12 Å box, `ISPIN=2`, `NUPDOWN=2`, Gamma only | `IBRION=5` (and `ISYM=0`, see below) |
| rutile SnO₂ | 6 atoms | `IBRION=6`, 2×2×3 supercell (72 atoms), 2×2×2 k-mesh (0.377 Å⁻¹) |
| α-Sn | 8 atoms | `IBRION=6`, 2×2×2 supercell (64 atoms), 3×3×3 k-mesh |

Slabs: three triple layers, 15 Å vacuum, terminations `o` (Sn₆O₁₂), `sno`
(Sn₆O₁₀) and `sn2o` (Sn₆O₈), symmetric, all atoms free, one relaxation each.

The SnO₂ vibrations use a coarser k-mesh (`kpoints_distance=0.06`) than the
rest and 24 h of walltime: SnO₂ is a wide-gap insulator and the supercell is
9.7 Å wide. Relaxation and static energy of SnO₂ keep the fine mesh. This was
not tested against a finer mesh.

## Results

(`results/`, graphs: references PK 1372, slabs PK 1539, profile `psteros_sno2_vibrations`)

| Quantity | Value | Expectation |
|---|---|---|
| SnO₂ lattice | a = 4.830 Å, c = 3.243 Å | about 4.83 / 3.24 Å (PBE) |
| α-Sn lattice | a = 6.652 Å | about 6.65 Å (PBE) |
| O₂ bond, ω | 1.233 Å, 1567 cm⁻¹ | 1.23 Å, 1550–1600 cm⁻¹ |
| S(O₂, 298.15 K, 1 bar) | 205.4 J/(mol K) | about 205 |
| Δμ_O(298.15 K, 1 bar), no ZPE | −0.272 eV | about −0.27 eV |
| ΔE_f(SnO₂), static energies | −4.925 eV per formula unit | about −5 eV (experiment −6.0) |
| Imaginary modes | none (3 acoustic modes at 0) | none above about 20 cm⁻¹ |
| Highest mode, SnO₂ / α-Sn | 706 / 173 cm⁻¹ | 700–800 / about 200 cm⁻¹ |
| ZPE, SnO₂ / α-Sn | 0.199 eV per formula unit / 0.021 eV per atom | 0.15–0.25 / about 0.02 eV |

Surface energies (J/m², DFT energies): `slab_o` 1.03 at every Δμ_O; `slab_sno`
2.65 and `slab_sn2o` 3.67 at Δμ_O = 0 (O-rich), 0.87 and 0.11 at the O-poor limit
(Δμ_O = −2.46 eV). `slab_o` is stable from the O-rich limit down to
Δμ_O = −1.83 eV, `slab_sn2o` below it. This is the same ordering and nearly the
same numbers as the QE example (`slab_o` 1.033 J/m²; window −2.43 eV).

Reading the Δμ_O axis as temperature and pressure (`oxygen_conditions.csv`,
`sno2_110_transitions.csv`), with zero-point energy:

| O₂ at | Δμ_O (eV) |
|---|---|
| 298.15 K, 1 bar | −0.224 |
| 600 K, 1 bar | −0.563 |
| 1000 K, 1 bar | −1.052 |
| 600 K, 10⁻¹⁰ bar | −1.159 |
| 1000 K, 10⁻⁶ bar | −1.648 |

The `slab_sn2o` → `slab_o` transition at Δμ_O = −1.83 eV is at p(O₂) = 6×10⁻²² bar
at 600 K and 1.6×10⁻⁸ bar at 1000 K. The O-poor limit of the stability window is
−2.462 eV with DFT energies (0 K), −2.414 eV at 600 K and −2.521 eV at 1000 K
with the free energies of SnO₂ and Sn (bare E(O₂) on the axis).

The diagram uses DFT energies for the slabs *and* the bulk: the slabs have no
vibrational free energy, and a bulk free energy inside γ without one for the
slabs would be inconsistent (see `docs/source/reference-thermochemistry.rst`).
The free energies are used for the T, p reading of the axis and for the
stability window. The slabs are thin (three triple layers) and the numbers are
a test of the workflow, not converged surface energies.

## What the run found in psteros

The run found three bugs in the new blocks, fixed with tests (see `LOG.md` and
`CHANGE.md`):

1. A displaced cell changes VASP's k-point set, which VASP refuses under band
   parallelisation (`NCORE > 1`): the vibrations block now defaults to
   `NCORE = 1` for a solid.
2. `NCORE = 1` on 128 ranks gives a 6-band molecule 128 bands and VASP fails in
   the diagonalisation: a gas gets `ISYM = 0` instead and keeps the recipe's `NCORE`.
3. `make_supercell` kept the atoms cell by cell (`Sn Sn O O O O Sn Sn …`), VASP read
   24 ion types, found no symmetry and planned 432 displacements for SnO₂; the
   supercell is now grouped by element (4 degrees of freedom, 8 displacements).
