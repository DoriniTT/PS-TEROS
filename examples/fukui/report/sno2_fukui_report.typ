// PS-TEROS Report: SnO2 (110) Surface Thermodynamics and Fukui Analysis
// Author: Generated with PS-TEROS
// Date: January 2026

#set document(
  title: "SnO2 (110) Surface Thermodynamics and Fukui Function Analysis",
  author: "PS-TEROS",
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set heading(numbering: "1.")

#set par(justify: true)

// Title
#align(center)[
  #text(size: 18pt, weight: "bold")[
    SnO2 (110) Surface Thermodynamics and Fukui Function Analysis
  ]

  #v(0.5cm)

  #text(size: 12pt)[
    PS-TEROS Validation Report
  ]

  #v(0.3cm)

  #text(size: 10pt, style: "italic")[
    January 2026
  ]
]

#v(1cm)

= Introduction

This report presents results from PS-TEROS (Predicting Stability of TERminations of Oxide Surfaces) for the SnO#sub[2] (110) surface. The study validates two main components: surface thermodynamics calculations using ab initio atomistic thermodynamics, and Fukui function analysis using the charge interpolation method. These results are compared against published literature and provide input for the FukuiGrid paper @barrera2025.

= Surface Thermodynamics

== Methodology

The surface Gibbs free energy was calculated following the Reuter & Scheffler formalism @reuter2001. For a binary oxide M#sub[x]O#sub[y], the surface energy depends on the oxygen chemical potential deviation:

$ gamma(Delta mu_"O") = phi - Gamma_"O" dot Delta mu_"O" $

where:
- $phi$ is the reference surface energy at the O-rich limit ($Delta mu_"O" = 0$)
- $Gamma_"O"$ is the surface oxygen excess per unit area
- $Delta mu_"O"$ ranges from 0 (O-rich, corresponding to equilibrium with O#sub[2] gas) to $Delta H_f slash y$ (O-poor, metal-rich limit)

Three distinct terminations of the SnO#sub[2] (110) surface were considered, shown in @fig-structures. T1 is oxygen-deficient with exposed Sn atoms, T2 is the stoichiometric termination with bridging oxygens, and T3 is oxygen-rich with additional surface oxygen.

#figure(
  box(width: 65%)[
    #grid(
      columns: 3,
      gutter: 0.5em,
      image("T1.png", width: 100%),
      image("T2.png", width: 100%),
      image("T3.png", width: 100%),
    )
  ],
  caption: [
    The three SnO#sub[2] (110) surface terminations studied: T1 (oxygen-deficient), T2 (stoichiometric), and T3 (oxygen-rich). Grey spheres represent Sn atoms and red spheres represent O atoms.
  ],
) <fig-structures>

== Results

@fig-surface-energy shows the calculated surface energies compared with the reference from literature @sno2_abinit. The stoichiometric termination T2 is the most stable throughout nearly the entire thermodynamically accessible range, from O-rich to O-poor conditions. This result is consistent with the ab initio thermodynamics study by Batzill and Diebold, which found the stoichiometric (110) termination to be thermodynamically favored under typical conditions. The surface energy values are in the range of 1.0--1.5 J/m#super[2], typical for stable oxide surfaces.

#figure(
  grid(
    columns: 2,
    gutter: 1em,
    image("sno2_110_gamma_vs_muO.pdf", width: 100%),
    image("ref_surf_diagram.png", width: 100%),
  ),
  caption: [
    Surface energy $gamma$ as a function of oxygen chemical potential deviation $Delta mu_"O"$ for SnO#sub[2] (110). Left: PS-TEROS calculation showing T2 (stoichiometric) as the most stable termination. Right: Reference diagram from literature @sno2_abinit showing the same stability ordering.
  ],
) <fig-surface-energy>

= Fukui Function Analysis

== Methodology

The Fukui function for nucleophilic attack (f#super[+]) was calculated using the charge interpolation method implemented in PS-TEROS. The Fukui function is defined as the derivative of electron density with respect to the number of electrons:

$ f^+(arrow(r)) = (partial rho(arrow(r))) / (partial N)|_(v(arrow(r))) $

where:
- $rho(arrow(r))$ is the electron density at position $arrow(r)$
- $N$ is the total number of electrons
- $v(arrow(r))$ is the external potential (held constant)

The calculation procedure involves:
- Static DFT calculations at multiple fractional charge states ($Delta N = 0.00, 0.05, 0.10, 0.15$ electrons removed)
- Collection of charge density (CHGCAR) from each calculation
- Interpolation using FukuiGrid to compute the Fukui function via numerical differentiation

== Results

The calculated Fukui function identifies sites susceptible to nucleophilic attack on the SnO#sub[2] (110) surface. @fig-fukui-3d shows the 3D isosurface visualization compared with the reference from the FukuiGrid paper. The yellow isosurfaces indicate regions where electron density decreases upon removal of electrons, corresponding to electrophilic sites. The reactive regions are localized on the bridging oxygen atoms and the five-fold coordinated Sn atoms, in agreement with the published results.

#figure(
  stack(
    dir: ttb,
    spacing: 1em,
    grid(
      columns: 2,
      gutter: 1em,
      image("top_fukui_plus_01_interpolation_result.png", width: 100%),
      image("side_fukui_plus_01_interpolation_result.png", width: 100%),
    ),
    align(center)[
      #image("article_fukui.png", width: 85%)
    ],
  ),
  caption: [
    3D visualization of the Fukui function f#super[+] for SnO#sub[2] (110). Top row: PS-TEROS calculation (left: top view, right: side view). Bottom: Reference from the FukuiGrid paper @barrera2025. The isosurfaces show excellent agreement in identifying reactive sites on bridging oxygens and Sn atoms.
  ],
) <fig-fukui-3d>

The planar average of the Fukui function along the surface normal direction (@fig-fukui-planar) provides additional insight into the spatial distribution of reactivity. The peaks indicate the z-positions where the surface is most reactive toward nucleophilic species, with the profile consistent with the 3D visualization showing reactivity concentrated at the surface layer.

#figure(
  image("fukui_planar_average.png", width: 95%),
  caption: [
    Planar averages along the surface normal (z-direction) for SnO#sub[2] (110). Left: Fukui function f#super[+], indicating regions of highest reactivity toward nucleophilic species. Right: Fukui potential f#super[+], representing the electrostatic potential associated with the Fukui function.
  ],
) <fig-fukui-planar>

@fig-fukui-potential-paper shows the reference Fukui potential from the FukuiGrid paper for comparison. The calculated planar averages show excellent agreement with the published results.

#figure(
  image("fukui_potential_paper.png", width: 70%),
  caption: [
    Reference Fukui potential from the FukuiGrid paper @barrera2025 for comparison with the calculated results in @fig-fukui-planar.
  ],
) <fig-fukui-potential-paper>

The Fukui function for electrophilic attack (f#super[-]) was also calculated, identifying sites susceptible to electron donation. @fig-fukui-minus shows the 3D isosurface visualization of f#super[-], where the isosurfaces indicate regions where electron density increases upon addition of electrons.

#figure(
  image("fukui_minus.png", width: 75%),
  caption: [
    3D visualization of the Fukui function f#super[-] for SnO#sub[2] (110), identifying sites susceptible to electrophilic attack. The isosurfaces indicate regions where electron density increases upon addition of electrons.
  ],
) <fig-fukui-minus>

= Conclusions

This report demonstrates the successful validation of PS-TEROS for surface science applications. The surface thermodynamics module correctly identifies the stoichiometric T2 termination as the most stable SnO#sub[2] (110) surface, consistent with published literature. The Fukui function module accurately predicts reactive sites on the surface using the charge interpolation method, with results matching the FukuiGrid paper. Both modules are validated and provide reliable predictions for oxide surface studies.

#bibliography("references.bib", style: "american-physics-society")
