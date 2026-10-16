# SnO₂(110) surface phase diagram with VASP

An end-to-end psteros example: from starting structures to γ(Δμ_O) of three
SnO₂(110) terminations, written both as a figure and as a CSV table.

| File | What it does |
|---|---|
| `campaign.py` | Builds (and optionally submits) the two AiiDA WorkGraphs. |
| `phase_diagram.py` | Reads the finished graphs and writes `<out>.png` and `<out>.csv`. |
| `vibrations.py` | Optional: Γ-point vibrations of the relaxed structures, for free energies at a temperature. |

## Before you start

* An AiiDA profile with a `vasp.vasp` code, an uploaded POTCAR family
  (`aiida-vasp potcar uploadfamily --path=... --name=PBE`) containing `Sn_d`
  and `O`, and a running daemon.
* psteros installed in the same Python environment as the daemon
  (`pip install .` from the repository root, then `verdi daemon restart`).

The numerical settings are deliberately small (3 triple-layer slabs, 10 Å
vacuum, 400 eV, `kpoints_spacing=0.3`) so that every job fits a short debug
queue. Converge the cutoff, k-points, slab thickness and vacuum before
interpreting the numbers.

## 1. Reference calculations

```bash
python campaign.py refs --profile MY_PROFILE --code vasp@my-cluster \
    --potential-family PBE --computer cluster --queue debug --submit
```

One graph, six jobs run one after the other: bulk rutile SnO₂ and α-Sn
(cell relaxation with `ISIF=3`, then a static calculation) and a triplet O₂
molecule at Γ (relaxation, then static). Omit `--submit` to build and inspect
the graph without running anything.

## 2. Slab terminations

```bash
python campaign.py slabs --profile MY_PROFILE --code vasp@my-cluster \
    --potential-family PBE --refs-pk <REFS_PK> --submit
```

The three terminations — `o` (Sn₆O₁₂, bridging O), `sno` (Sn₆O₁₀) and `sn2o`
(Sn₆O₈) — are cut from the **relaxed** bulk lattice of step 1, so bulk and
slab energies refer to the same cell. The central triple layer is fixed with
selective dynamics (`CalculationOverride(fixed_sites=psteros.central_sites(slab))`).

## 3. The phase diagram

```bash
python phase_diagram.py --profile MY_PROFILE --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> \
    --out results/sno2_110_phase_diagram
```

This prints ΔH_f, the stability window and the transitions, and writes the
figure drawn by `SurfacePhaseDiagram.plot` (use `--format pdf` or `svg` for
vector output) and the CSV table with one row per Δμ_O.

## 4. Optional: vibrational free energies

```bash
python vibrations.py --profile MY_PROFILE --code vasp@my-cluster --potential-family PBE \
    --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> --parallel-jobs 4 --submit
python phase_diagram.py --profile MY_PROFILE --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> \
    --vib-pk <VIB_PK> --temperature 800 --out results/sno2_110_800K
```

`vibrations.py` displaces the sites that were free in the slab relaxations
(the fixed central triple layer is one Sn₂O₄ bulk cell and is counted as bulk),
every site of a Γ-only rutile supercell (`--bulk-supercell`, 2×2×2 by default)
and of α-Sn, and the O₂ molecule, with `EDIFF = 1e-7`. Every displaced site
costs 6 static calculations; `--parallel-jobs` runs several at once.
`phase_diagram.py --vib-pk` then uses E + F_vib(T) for the slabs and solids
and E + ZPE for O₂. Without `--vib-pk` it gives the 0 K diagram as before. See
[the guide](../../docs/source/vibrations.rst).

The same campaign with Quantum ESPRESSO is
[`examples/qe_surface_phase_diagram`](../qe_surface_phase_diagram).
