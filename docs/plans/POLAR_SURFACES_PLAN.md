# Plan: any-material thermodynamics and absolute polar surface energies

Branch: `feature/charge-neutral-terminations`. One step at a time; each step
ends with tests passing and a commit. Status: `[ ]` to do, `[~]` in progress,
`[x]` done.

## Decisions (agreed)

- Scope of the polar method: III-V and II-VI compounds first (zinc blende
  (111)/(-1-1-1), wurtzite (0001)/(000-1)).
- Bottom: pseudo-hydrogen passivation (charge 2 - Z/4 per broken bond).
- Relaxation: relax every atom, as in the papers, then check that the bottom
  region is the same in every slab after relaxation.
- Pseudo chemical potential of pseudo-H: pseudo-molecules by default, tetrahedral
  clusters optional.
- Energies: absolute, so polar and non-polar orientations can be compared.
- DFT code: VASP only for now.
- Same-bottom requirement: all terminations of one polar face are built on one
  identical bottom (same atoms, pseudo-H, cell); enforced by a fingerprint and
  by the post-relaxation check.
- Thermodynamics work goes into the public API (`psteros/phase_diagram.py`,
  `psteros/thermodynamics.py`, `psteros/workflow.py`); the legacy
  `psteros/core/` modules stay as they are.
- Papers in `references/` (git-ignored, never committed):
  Zhang et al., Sci. Rep. 6, 20055 (2016) and Zhang et al., arXiv:1510.08961.

## Method (from the papers)

- sigma_top = [E_slab - sum_i n_i mu_i - sum_k n_Hk muhat_Hk] / A
  (Sci. Rep. Eqs. 4-5); muhat absorbs the passivated bottom energy.
- Pseudo-molecule: muhat_HX = [E(X H_X,4) - mu_X] / 4 (Eq. 8); molecule of the
  bonded atom X with four pseudo-H; 8 valence electrons, Gamma only.
- Clusters: E(n) = n(n+1)(n+2)/6 mu_A + (n-1)n(n+1)/6 (E_AB - mu_A)
  + 2(n-2)(n-3) muhat_face + 12(n-2) muhat_edge + 12 muhat_corner (Eq. 9);
  solve or least-squares fit over n = 2..9. Wurtzite uses zinc-blende clusters.
- muhat depends on Delta mu: d(muhat_HX)/d(mu_X) = -1/4.
- Checks: both-faces-passivated slab gives muhat_HA + muhat_HB (Eq. 7);
  a non-polar surface computed both ways must agree.
- Settings used: VASP PBE, 400-500 eV, >= 15 A vacuum, 9-10 bilayers, 1x1,
  ~10x10x1 to 15x15x1 k-points, molecules/clusters at Gamma, forces 0.005 eV/A,
  all atoms relaxed.
- Benchmarks (GGA, anion-rich, cluster method, meV/A^2): ZnO(0001) 147.7,
  ZnO(000-1) 63.1, GaN(0001) 168.3, GaN(000-1) 198.2; GaAs V_Ga (111)-2x2 39.2.

## Steps

- [x] 1. **Name and wording.** "Predicting Stability of TERminations Of
  Surfaces" in README, docs and package docstrings; describe oxides,
  semiconductors and other compounds. Keep the PS-TEROS name and the citation.
- [x] 2. **Any binary compound.** `BinaryReferences` for A_xB_y with elemental
  references for A and B (solid or half a molecule); the axis is Delta mu of a
  chosen element (default: the anion). `surface_phase_diagram` accepts it;
  `BinaryOxideReferences` stays as the oxide special case with identical
  results, CSV columns and figure. Labels become element-generic.
- [x] 3. **Any ternary compound.** Same generalisation for
  `TernaryOxideReferences` (any third element), oxide results unchanged.
- [x] 4. **Pseudo-hydrogen model.** Charge 2 - Z/4, formal charge of a pseudo-H
  (minus the oxidation state of its partner over 4), VASP POTCAR names
  (H.5, H.75, H1.25, H1.5, ...), AiiDA kind names via the `kind_name` site
  property.
- [x] 5. **Polar slab builder.** `find_polar_terminations` for zinc blende
  (111)/(-1-1-1) and wurtzite (0001)/(000-1): N bilayers, one shared bottom
  passivated with pseudo-H, the ideal top plus electron-counting variants
  (e.g. 2x2 vacancy), all in one common cell; bottom fingerprint; summary
  table and plots.
- [x] 6. **Pseudo chemical potentials.** Pseudo-molecule builder; optional
  tetrahedral-cluster builder and solver; `PseudoHydrogenReferences` holding
  muhat(Delta mu); both-faces-passivated check slab (Eq. 7).
- [x] 7. **Polar thermodynamics.** One-face terminations with pseudo-H counts
  in `surface_phase_diagram`, absolute gamma on the same grid as symmetric
  slabs; refuses terminations of one face with different bottoms.
- [ ] 8. **Post-relaxation bottom check.** Compare the relaxed bottom region of
  every slab with the shared reference (RMSD after removing a rigid shift);
  failing slabs are flagged and left out of comparisons.
- [ ] 9. **VASP workflow.** Structures, potential mapping and INCAR overrides
  (dipole correction for asymmetric slabs, Gamma-only molecules and clusters)
  for `build_surface_workgraph`; a collector from finished energies to the
  phase diagram.
- [ ] 10. **Validation helpers.** Eq. 7 self-consistency and the non-polar
  comparison, reported in meV/A^2.
- [ ] 11. **Docs and example.** Guide section, example script and a benchmark
  recipe against the paper values.
- [ ] 12. **Final check.** Full test suite, lint, docs build references,
  changelog, push.
