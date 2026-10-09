# SnO₂(110) campaign in one WorkGraph (VASP)

The VASP campaign of `examples/vasp_sno2_phase_diagram_vibrations/`, built as one graph with
`psteros.build_vasp_campaign_workgraph`: the references O₂, rutile SnO₂ and α-Sn (relax → static →
vibrations) and the three SnO₂(110) terminations (relax → static) run in a single WorkGraph with one PK.
`psteros.campaign_results` and `psteros.campaign_terminations` read it back by PK, and `analysis.py`
writes the phase diagram. Only the public psteros API is used. The INCAR, the overrides, the k-point
choices and the execution policy are those of the earlier example.

| File | What it does |
|---|---|
| `campaign.py` | Builds, prints and (with `--submit`) submits the graph. |
| `analysis.py` | Reads the finished graph by `--pk` and writes the phase diagram (PNG, CSV) and the transitions. |

## Before you start

* An AiiDA profile with a `vasp.vasp` code, a POTCAR family (`PBE`, with `Sn_d` and `O`) and a running daemon.
  The code's MPI launcher must match the one VASP was built with (on Lovelace, VASP 6.5.1 is linked against Intel MPI).
* psteros importable by the **daemon** (editable install, or `PYTHONPATH`), and the daemon restarted after every
  psteros change: it runs the calcfunctions `vasp_energy`, `vasp_frequencies` and `make_supercell`.
* `campaign.py` always builds the graph and prints, for every VASP task, its INCAR, k-point spacing, scheduler
  options, settings and structure. Read that before adding `--submit`.

## Run it

```bash
python campaign.py                              # build and print only
python campaign.py --smoke --submit             # triplet O2 alone: a cheap test of every code path
python campaign.py --submit                     # the whole campaign; note the PK it prints
python analysis.py --pk <PK> --out results/sno2_110
```

The graph runs one job at a time (`max_concurrent_jobs=1`). `analysis.py` needs the full campaign, not `--smoke`.

## Settings

One INCAR for every block: `ENCUT=520`, `PREC=Accurate`, `EDIFF=1e-7`, `ISMEAR=0`, `SIGMA=0.05`, `LREAL=False`,
`LASPH=True`, `NCORE=16`, `LWAVE=LCHARG=False`, `IBRION=2`, `NSW=100`; POTCARs `Sn_d` and `O`;
`kpoints_spacing = 2π·0.03 = 0.1885 Å⁻¹` (with the 2π, as VASP's `KSPACING`).

| Group | Overrides |
|---|---|
| references | `EDIFFG = -0.005` (tight: the frequencies are computed at this minimum); `Vibrations` with `EDIFF = 1e-8` |
| O₂ (gas) | 12 Å box, `ISPIN = 2`, `NUPDOWN = 2`, Gamma only |
| SnO₂ (solid, 2×2×3 supercell) | `ISIF = 3` for the relaxation; `Vibrations` at `kpoints_distance` 0.06·2π Å⁻¹ with 24 h walltime |
| α-Sn (solid, 2×2×2 supercell) | `ISIF = 3` for the relaxation |
| slabs | `ISIF = 2` and `EDIFFG = -0.02`, all atoms free; three triple layers, 15 Å vacuum; terminations `o` (Sn₆O₁₂), `sno` (Sn₆O₁₀), `sn2o` (Sn₆O₈) |

## The slabs are an input

The graph does not cut the slabs from the bulk it relaxes. `campaign.py` cuts them from the SnO₂ lattice given by
`--a` and `--c`, by default the relaxed lattice of the earlier run (a = 4.8301 Å, c = 3.2434 Å). If your settings
relax the bulk to another lattice, pass that one, so that slab and bulk energies refer to the same cell.

## Differences from the two-graph example

* One graph and one PK, instead of a references graph and a slabs graph.
* The slab energy is the static energy of the relaxed slab, as for the references. The earlier example used the
  σ → 0 energy of the relaxation.
* `--smoke` replaces the `o2` phase of the earlier example.

Results of this graph are not recorded here yet. The numbers of the two-graph run are in
`../vasp_sno2_phase_diagram_vibrations/`.
