# SnO₂(110) surface phase diagram with Quantum ESPRESSO

An end-to-end psteros example: from starting structures to γ(Δμ_O) of three
SnO₂(110) terminations, written both as a figure and as a CSV table.

| File | What it does |
|---|---|
| `campaign.py` | Builds (and optionally submits) the two AiiDA WorkGraphs. |
| `phase_diagram.py` | Reads the finished graphs and writes `<out>.png` and `<out>.csv`. |
| `results/` | Output of the smoke test on Bohr's `teste` queue (see below). |

## Before you start

* An AiiDA profile with a `quantumespresso.pw` code, an `aiida-pseudo` family
  (the test used `SSSP/1.3/PBE/efficiency`) and a running daemon.
* psteros installed in the **same Python environment as the daemon**
  (`pip install .` from the repository root, then `verdi daemon restart`):
  the relaxation stage runs a psteros work chain on the daemon.

The numerical settings are deliberately small (3 triple-layer slabs, 10 Å
vacuum, 40/320 Ry, `kpoints_distance=0.5`) so that every job fits a 20-minute,
4-core debug queue. Converge cutoffs, k-points, slab thickness and vacuum
before interpreting the numbers.

## 1. Reference calculations

```bash
python campaign.py refs --profile MY_PROFILE --code QE-7.3.1@bohr \
    --pseudo-family SSSP/1.3/PBE/efficiency --submit
```

One graph, six jobs run one after the other: bulk rutile SnO₂ and α-Sn
(`vc-relax` → static SCF) and a triplet O₂ molecule (`relax` → static SCF).
Omit `--submit` to build and inspect the graph without running anything.
Use `--computer` and `--queue` for a machine other than Bohr's `teste` queue.

## 2. Slab terminations

```bash
python campaign.py slabs --profile MY_PROFILE --code QE-7.3.1@bohr \
    --pseudo-family SSSP/1.3/PBE/efficiency --refs-pk <REFS_PK> --submit
```

The three terminations — `o` (Sn₆O₁₂, bridging O), `sno` (Sn₆O₁₀) and `sn2o`
(Sn₆O₈) — are cut from the **relaxed** bulk lattice of step 1, so bulk and
slab energies refer to the same cell. The central triple layer is frozen with
`psteros.qe_fixed_coordinate_flags`. Each relaxation may need a walltime
restart, which `PwBaseWorkChain` performs automatically (`max_iterations=4`).

## 3. The phase diagram

```bash
python phase_diagram.py --profile MY_PROFILE --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> \
    --out results/sno2_110_phase_diagram
```

This prints ΔH_f, the stability window and the transitions, and writes:

* `results/sno2_110_phase_diagram.png` — the figure drawn by
  `SurfacePhaseDiagram.plot` (use `--format pdf` or `svg` for vector output);
* `results/sno2_110_phase_diagram.csv` — one row per Δμ_O value with columns
  `delta_mu_O_eV`, `gamma_<termination>_Jm2`, `stable_termination` and
  `in_stability_window`, for plotting in any other tool.

The analysis itself is four psteros calls, usable with energies from any source:

```python
references = psteros.BinaryOxideReferences(
    bulk_energy_ev=e_bulk, bulk_composition=bulk.composition,
    oxygen_molecule_energy_ev=e_o2, metal_energy_per_atom_ev=e_sn / len(sn_cell),
)
terminations = [psteros.SlabTermination.from_structure(label, energy, structure) for ...]
diagram = psteros.surface_phase_diagram(terminations, references)
diagram.plot("sno2_110.png")
diagram.to_csv("sno2_110.csv")
```

## Smoke-test results (Bohr `teste`, QE 7.3.1, 4 MPI ranks)

| Quantity | Value |
|---|---|
| Relaxed rutile a, c | 4.784 Å, 3.231 Å |
| Relaxed α-Sn a | 6.713 Å |
| O₂ bond, magnetisation | 1.222 Å, 2.00 μ_B |
| ΔH_f(SnO₂) | −4.87 eV per formula unit |
| Stability window | −2.43 ≤ Δμ_O ≤ 0 eV |
| γ, O-bridge termination | 1.03 J/m² |
| Transition | O-bridge → Sn₆O₈ at Δμ_O = −1.79 eV (Sn₆O₁₀ is never the most stable) |

About 2 h 40 min of wall time for 15 `pw.x` jobs, including queue waits.
These are smoke-test numbers from thin, unconverged slabs, useful to check
that the workflow runs end to end, not as physical results.
