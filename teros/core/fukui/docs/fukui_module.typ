// Fukui Module Documentation
// PS-TEROS - Predicting Stability of TERminations of Oxide Surfaces

#set document(
  title: "Fukui Module Documentation",
  author: "PS-TEROS Development Team",
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
  numbering: "1",
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set heading(numbering: "1.1")

#show link: underline

// Title
#align(center)[
  #text(size: 24pt, weight: "bold")[Fukui Module]

  #text(size: 14pt)[PS-TEROS Documentation]

  #v(0.5cm)

  #text(size: 11pt, style: "italic")[
    Calculating Fukui Functions via Finite Difference Interpolation
  ]

  #v(1cm)
]

// Table of Contents
#outline(
  title: "Contents",
  indent: 1.5em,
)

#pagebreak()

= Introduction

== What are Fukui Functions?

Fukui functions are conceptual DFT descriptors that measure the local reactivity of a molecular or surface system. They quantify how the electron density at a given point changes when electrons are added or removed from the system.

The Fukui functions are defined as:

$ f^+(arrow(r)) = (partial rho(arrow(r))) / (partial N) |_(v(arrow(r))) approx rho_(N+delta N)(arrow(r)) - rho_N(arrow(r)) $

$ f^-(arrow(r)) = (partial rho(arrow(r))) / (partial N) |_(v(arrow(r))) approx rho_N(arrow(r)) - rho_(N-delta N)(arrow(r)) $

Where:
- $f^+(arrow(r))$ is the *nucleophilic Fukui function* — indicates where electrons will be added (electrophilic attack sites)
- $f^-(arrow(r))$ is the *electrophilic Fukui function* — indicates where electrons will be removed (nucleophilic attack sites)
- $rho(arrow(r))$ is the electron density
- $N$ is the number of electrons
- $v(arrow(r))$ is the external potential (kept constant)

== Why Use Fukui Functions?

Fukui functions are powerful tools for:

- *Predicting reactive sites* on surfaces and molecules
- *Understanding adsorption* — where adsorbates will preferentially bind
- *Catalyst design* — identifying active sites for reactions
- *Surface chemistry* — predicting oxidation/reduction behavior

== The Interpolation Method

Instead of using integer electron differences (which can cause convergence issues), we use *fractional charges* with interpolation:

1. Calculate charge densities at multiple fractional charge states:
   - $delta N = 0.00$ (neutral)
   - $delta N = 0.05$
   - $delta N = 0.10$
   - $delta N = 0.15$

2. Fit a polynomial to the charge density vs. $delta N$ data

3. Extract the derivative at $delta N = 0$ to obtain the Fukui function

This approach provides:
- Better numerical stability
- Smoother Fukui functions
- More accurate results for extended systems (surfaces, slabs)

#pagebreak()

= Module Architecture

== Overview

The Fukui module in PS-TEROS automates the calculation of Fukui functions using AiiDA-WorkGraph for workflow management and VASP for DFT calculations.

#figure(
  ```
  teros/core/fukui/
  ├── __init__.py      # Public API exports
  ├── workgraph.py     # Main WorkGraph builder
  ├── tasks.py         # AiiDA calcfunctions
  ├── utils.py         # Utility functions
  └── docs/            # Documentation
  ```,
  caption: "Module structure"
)

== Components

=== `utils.py` — Utility Functions

#table(
  columns: (auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Function*], [*Description*],
  [`make_delta_label(delta_n)`], [Converts float (0.05) to valid label ("delta_0_05")],
  [`validate_fukui_inputs(...)`], [Validates input parameters before workflow execution],
  [`DEFAULT_DELTA_N_VALUES`], [Default: \[0.0, 0.05, 0.10, 0.15\]],
)

=== `tasks.py` — AiiDA Calcfunctions

#table(
  columns: (auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Function*], [*Description*],
  [`extract_total_energy(misc)`], [Extracts energy from VASP output],
  [`collect_chgcar_files(...)`], [Collects CHGCAR files from multiple calculations into single FolderData],
  [`generate_fukui_summary(...)`], [Creates summary Dict with calculation metadata],
)

=== `workgraph.py` — WorkGraph Builder

#table(
  columns: (auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Function*], [*Description*],
  [`FukuiCalculationScatter`], [`@task.graph` for parallel VASP calculations],
  [`build_fukui_workgraph(...)`], [Main entry point — builds complete workflow],
  [`get_fukui_results(wg)`], [Extracts results from completed WorkGraph],
  [`print_fukui_summary(wg)`], [Prints formatted summary of results],
)

#pagebreak()

= Usage Guide

== Basic Usage

```python
from aiida import orm, load_profile
from teros.core.fukui import build_fukui_workgraph, get_fukui_results

load_profile('your_profile')

# Load or create structure
structure = orm.load_node(12345)  # or from file

# Build the workflow
wg = build_fukui_workgraph(
    structure=structure,
    nelect_neutral=192,           # REQUIRED: verify from VASP OUTCAR
    delta_n_values=[0.0, 0.05, 0.10, 0.15],
    code_label='VASP-6.5.1@cluster',
    builder_inputs={
        'parameters': {
            'incar': {
                'encut': 500,
                'ediff': 1e-6,
                'algo': 'All',
                'ispin': 2,
                # ... other settings
            }
        },
        'options': {...},
        'kpoints_spacing': 0.25,
        'potential_family': 'PBE.54',
        'potential_mapping': {'Sn': 'Sn_d', 'O': 'O'},
    },
    fukui_type='plus',            # 'plus' or 'minus'
    max_concurrent_jobs=4,
)

# Submit
wg.submit(wait=False)
print(f"Submitted: PK={wg.pk}")
```

== Parameters

#table(
  columns: (auto, auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Parameter*], [*Required*], [*Description*],
  [`structure`], [Yes], [Input structure (StructureData or PK)],
  [`nelect_neutral`], [Yes], [Number of electrons in neutral system],
  [`delta_n_values`], [No], [List of δN values (default: \[0.0, 0.05, 0.10, 0.15\])],
  [`code_label`], [Yes], [VASP code label],
  [`builder_inputs`], [Yes], [VASP calculation parameters],
  [`relax_first`], [No], [Relax structure before charge calculations (default: False)],
  [`relax_builder_inputs`], [No], [Separate inputs for relaxation step],
  [`fukui_type`], [No], [`'plus'` (nucleophilic) or `'minus'` (electrophilic)],
  [`max_concurrent_jobs`], [No], [Limit parallel calculations],
)

== Determining NELECT

The `nelect_neutral` parameter is *critical* and must be determined from your pseudopotentials:

```bash
# Run a test VASP calculation and check OUTCAR:
grep "NELECT" OUTCAR
```

For common pseudopotentials:
- Standard PBE: Sum of valence electrons per element
- Example: SnO₂ with 12 Sn + 24 O
  - Sn (standard): 4 valence electrons
  - Sn_d: 14 valence electrons
  - O: 6 valence electrons
  - NELECT = 12×4 + 24×6 = 192 (standard) or 12×14 + 24×6 = 312 (Sn_d)

== VASP Settings for Fukui Calculations

The module automatically enforces these settings for charge-state calculations:

```python
# Automatically set by the module:
'nsw': 0,        # Static calculation
'ibrion': -1,    # No ionic relaxation
'lcharg': True,  # Write CHGCAR
'lwave': False,  # Don't write WAVECAR
'nelect': <calculated>,  # Fractional electron count
```

*Recommended user settings* for robust SCF convergence:

```python
'incar': {
    'prec': 'Accurate',
    'encut': 500,           # Well-converged cutoff
    'algo': 'All',          # Most robust algorithm
    'nelm': 300,            # Allow many SCF steps
    'ediff': 1e-6,          # Tight convergence
    'icharg': 2,            # Start from atomic charges
    'amix': 0.2,            # Conservative mixing
    'bmix': 0.0001,         # Kerker mixing
    'lmaxmix': 4,           # For d-electrons
    'ismear': 0,            # Gaussian smearing
    'sigma': 0.05,
    'ispin': 2,             # Spin-polarized
}
```

#pagebreak()

== Extracting Results

After the workflow completes:

```python
from teros.core.fukui import get_fukui_results, print_fukui_summary

# Get results
results = get_fukui_results(wg.pk)

# Access CHGCAR files
chgcar_folder = results['chgcar_folder']
print(chgcar_folder.list_object_names())
# ['CHGCAR_0.00', 'CHGCAR_0.05', 'CHGCAR_0.10', 'CHGCAR_0.15']

# Access summary
summary = results['summary']
print(summary['fukui_type'])
print(summary['calculations'])

# Print formatted summary
print_fukui_summary(wg.pk)
```

== Exporting CHGCAR Files

To export CHGCAR files for processing with FukuiGrid.py:

```python
import os
from pathlib import Path

output_dir = Path('./fukui_chgcars')
output_dir.mkdir(exist_ok=True)

chgcar_folder = results['chgcar_folder']
for fname in chgcar_folder.list_object_names():
    content = chgcar_folder.get_object_content(fname)
    output_path = output_dir / fname

    if isinstance(content, bytes):
        output_path.write_bytes(content)
    else:
        output_path.write_text(content)

print(f"Exported to: {output_dir}")
```

#pagebreak()

= Workflow Architecture

== Data Flow Diagram

#figure(
  ```
  Input: structure, nelect_neutral, builder_inputs
           │
           ▼
  ┌─────────────────────────────────┐
  │ (Optional) VaspTask relaxation  │
  │ if relax_first=True             │
  └────────────┬────────────────────┘
               │ relaxed_structure
               ▼
  ┌─────────────────────────────────┐
  │ FukuiCalculationScatter         │
  │ @task.graph (parallel)          │
  │                                 │
  │  ┌─────────┐ ┌─────────┐       │
  │  │VaspTask │ │VaspTask │ ...   │
  │  │NELECT=  │ │NELECT=  │       │
  │  │N        │ │N-0.05   │       │
  │  └────┬────┘ └────┬────┘       │
  │       │           │            │
  │  retrieved   retrieved         │
  └───────┴───────────┴────────────┘
          │           │
          ▼           ▼
  ┌─────────────────────────────────┐
  │ collect_chgcar_files            │
  │ @task.calcfunction              │
  │                                 │
  │ CHGCAR_0.00, CHGCAR_0.05, ...   │
  └────────────┬────────────────────┘
               │
               ▼
  ┌─────────────────────────────────┐
  │ generate_fukui_summary          │
  │ @task.calcfunction              │
  └────────────┬────────────────────┘
               │
               ▼
        WorkGraph Outputs:
        - chgcar_files (FolderData)
        - summary (Dict)
        - relaxed_structure (optional)
  ```,
  caption: "Fukui workflow data flow"
)

== NELECT Calculation

For *Fukui+* (nucleophilic attack — where electrons are added):
$ "NELECT" = N_"neutral" - delta N $

For *Fukui-* (electrophilic attack — where electrons are removed):
$ "NELECT" = N_"neutral" + delta N $

#figure(
  table(
    columns: (auto, auto, auto),
    stroke: 0.5pt,
    inset: 8pt,
    [*δN*], [*Fukui+ NELECT*], [*Fukui- NELECT*],
    [0.00], [192.00], [192.00],
    [0.05], [191.95], [192.05],
    [0.10], [191.90], [192.10],
    [0.15], [191.85], [192.15],
  ),
  caption: [Example NELECT values for $N_"neutral" = 192$]
)

#pagebreak()

= Future Development

== Planned Features

The Fukui module will be extended with the following capabilities:

=== Phase 2: Fukui Function Computation

Integration with FukuiGrid.py for computing the actual Fukui function from CHGCAR files:

- *Polynomial interpolation* of charge densities
- *Derivative extraction* at δN = 0
- Output: `CHGCAR_FUKUI.vasp` file

=== Phase 3: Fukui Potential

Calculation of the Fukui potential $v_f(arrow(r))$ using different methods:

#table(
  columns: (auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Method*], [*Description*],
  [Interpolation], [Direct computation from charge densities],
  [Electrodes], [Fourier transform-based electrostatic corrections],
  [SCPC], [Self-consistent potential correction method],
)

The Fukui potential is useful for:
- Predicting interaction energies with charged adsorbates
- Understanding electrostatic contributions to reactivity

=== Phase 4: Perturbative Expansion

Implementation of the perturbative expansion model:

$ Delta U(arrow(r)) = q Phi(arrow(r)) - q dot Delta N dot v_f^plus.minus(arrow(r)) $

Where:
- $q$ is the charge of the active site
- $Phi(arrow(r))$ is the electrostatic potential
- $Delta N$ is the electron transfer
- $v_f^plus.minus$ is the Fukui potential

This allows prediction of:
- Adsorption site preferences
- Interaction energy landscapes
- Reactivity maps

=== Phase 5: Workflow Integration

- Integration with main PS-TEROS `build_core_workgraph()`
- Workflow preset: `workflow_preset='fukui_analysis'`
- Automatic Fukui analysis after slab relaxation

#pagebreak()

= Troubleshooting

== Common Issues

#table(
  columns: (auto, 1fr),
  stroke: 0.5pt,
  inset: 8pt,
  [*Issue*], [*Solution*],
  [SCF not converging], [Use `algo: 'All'`, reduce `amix`, increase `nelm`],
  [CHGCAR not retrieved], [Check `ADDITIONAL_RETRIEVE_LIST` includes `'CHGCAR'`],
  [Wrong NELECT], [Verify from VASP OUTCAR of a test calculation],
  [Memory issues], [Reduce `ncore`, check node memory limits],
  [Workflow stuck], [Check `verdi process report <PK>` for sub-process status],
)

== Monitoring Workflows

```bash
# Check workflow status
verdi process show <PK>

# Detailed task hierarchy
verdi process report <PK>

# Check daemon logs
verdi daemon logshow

# List recent processes
verdi process list -a -p 1
```

== Validating Results

1. *Check convergence*: All calculations should show `electronic_converged: True`
2. *Energy consistency*: Energies should vary smoothly with δN
3. *CHGCAR integrity*: Files should have correct grid dimensions
4. *Physical reasonableness*: Fukui functions should be localized at reactive sites

= References

1. Parr, R. G., & Yang, W. (1984). Density functional approach to the frontier-electron theory of chemical reactivity. _Journal of the American Chemical Society_, 106(14), 4049-4050.

2. Yang, W., & Mortier, W. J. (1986). The use of global and local molecular parameters for the analysis of the gas-phase basicity of amines. _Journal of the American Chemical Society_, 108(19), 5708-5711.

3. FukuiGrid documentation and examples: https://github.com/rubenfiszel/FukuiGrid

4. PS-TEROS documentation: See `CLAUDE.md` in the repository root.
