# SrTiO₃(001) validation of the ternary phase diagram

A physical validation of `psteros.ternary_surface_phase_diagram`: real Quantum
ESPRESSO PBE calculations of SrTiO₃, its elements and competing oxides, and the
SrO- and TiO₂-terminated (001) surfaces, compared with experiment and with the
literature on the same slab models.

| File | What it does |
|---|---|
| `campaign.py` | Builds (and optionally submits) the three graphs below. |
| `analysis.py` | Builds the phase diagram and writes `validation.md`, `validation.json`, `phase_diagram.png`, `phase_diagram.csv`. |
| `compare_literature.py` | Puts the results and published LDA/B3PW calculations on one convention; writes `literature_comparison.md` and `.png`. |
| `results/` | Output of the validation run (one 44-core node per job). |

## Calculations

All with SSSP 1.3 PBE efficiency pseudopotentials, 60/480 Ry, Marzari–Vanderbilt
smearing of 0.01 Ry and `kpoints_distance = 0.2` Å⁻¹ (0.15 for the metals, Γ for
O₂). Every bulk phase is `vc-relax`ed and then recomputed in a static SCF.

1. `refs` — SrTiO₃ (cubic), α-Sr (fcc), α-Ti (hcp), SrO (rock salt), TiO₂ rutile
   and anatase, and a triplet O₂ molecule: 14 jobs.
2. `slabs` — symmetric 7-layer 1×1 slabs, SrO-terminated Sr₄Ti₃O₁₀ and
   TiO₂-terminated Sr₃Ti₄O₁₁, cut from the relaxed lattice with 22 Å of
   vacuum, central plane fixed; relax → static: 4 jobs.
3. `unrelaxed` — static SCF of the same slabs before relaxation, for the
   cleavage energy: 2 jobs.

```bash
python campaign.py refs      --profile P --code pw@my-cluster --queue my-queue --ranks 32 --submit
python campaign.py slabs     --profile P --code pw@my-cluster --refs-pk <REFS_PK> --submit
python campaign.py unrelaxed --profile P --code pw@my-cluster --refs-pk <REFS_PK> --submit
python analysis.py --profile P --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> --unrelaxed-pk <UNRELAXED_PK> --out results
```

Set `--ranks`, `--walltime` and `--queue` for your computer (the number of
ranks must be divisible by the k-point pools, 4 and 2). A QE wrapper that
launches MPI itself takes `--no-mpi` and its settings through `--prepend`, for
example `--no-mpi --prepend "export QE_MPI_RANKS=44"`. The recorded results
used one 44-core node per job with 44 ranks, k-point pools through `-nk` and
one active job at a time. Run the three graphs one after the other.

## What is compared, and with what

| Quantity | Reference | Source |
|---|---|---|
| ΔH_f(SrTiO₃) | −1668.99 ± 9.2 kJ/mol | bomb calorimetry, *J. Chem. Thermodyn.* study of strontium titanates (Sr₂TiO₄, Sr₃Ti₂O₇, Sr₄Ti₃O₁₀, SrTiO₃) |
| ΔH_f(SrO) | −592.04 kJ/mol | NIST-JANAF tables (Chase, 1998), via the NIST Chemistry WebBook |
| ΔH_f(TiO₂, rutile) | −944.0 ± 0.8 kJ/mol | CODATA (Cox, Wagman et al., 1984), via the NIST Chemistry WebBook |
| ΔH(SrO + TiO₂ → SrTiO₃) | −1.38 ± 0.10 eV | from the three values above |
| a(SrTiO₃) | 3.89 Å | experiment extrapolated to 0 K, quoted by Eglitis & Vanderbilt |
| Cleavage, relaxation and surface energies of the same 7-layer slabs | 1.39; −0.24 / −0.16; 1.15 / 1.23 eV per surface cell | R. I. Eglitis and D. Vanderbilt, *Phys. Rev. B* **77**, 195408 (2008), Table VII (B3PW hybrid) |
| Surface rumpling *s*, Δd₁₂, Δd₂₃ | B3PW, LEED and RHEED values | same paper, Table III |
| O₂ bond length | 1.23 Å | plane-wave PBE, V. Alexandrov et al., arXiv:1005.4833 |

What each comparison tests:

* **Formation enthalpies.** PBE overbinds O₂, so PBE formation enthalpies of
  oxides are less negative than experiment; a deviation of a few tenths of an
  eV per oxygen atom is expected, not a failure. The reaction energy between
  oxides, SrO + TiO₂ → SrTiO₃, largely cancels that error and is the sharper test.
* **Stability region.** Bounded by SrO and TiO₂, the SrTiO₃ strip in
  Δμ_SrO = Δμ_Sr + Δμ_O is exactly as wide as −ΔH(SrO + TiO₂ → SrTiO₃): this
  checks the psteros polygon against the reaction energy.
* **Slab energetics.** The SrO- and TiO₂-terminated slabs together contain seven
  bulk units, so the average of their surface energies does not depend on the
  chemical potentials. Using Eglitis and Vanderbilt's definitions
  (E_cleav = ¼[E_unrel(SrO) + E_unrel(TiO₂) − 7 E_bulk], E_rel = ½[E_relaxed − E_unrel],
  E_surf = E_cleav + E_rel), the numbers are directly comparable with their
  hybrid-functional results; the psteros γ planes must reproduce the same average
  at every point of the diagram.
* **Relaxation geometry.** Rumpling and interlayer changes, in % of the lattice
  constant and from the metal-ion positions, as in that paper.

The functional differs from Eglitis and Vanderbilt (PBE here, B3PW there),
and the phase diagram includes only the elements, SrO and the two TiO₂
polymorphs as competing phases; Ruddlesden–Popper phases such as Sr₂TiO₄ would
narrow the SrO-rich side further.
