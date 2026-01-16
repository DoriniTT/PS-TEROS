# Fukui Module TODO

## Completed

- [x] Phase 1: Interpolation method
  - [x] Parallel VASP calculations at different charge states
  - [x] Collect CHGCAR files into FolderData
  - [x] FukuiGrid integration (`compute_fukui=True` parameter)
  - [x] aiida-shell based execution for FukuiGrid interpolation
  - [x] Output: CHGCAR_FUKUI.vasp as SinglefileData

- [x] Phase 2: Fukui potential via electrodes method
  - [x] DFPT calculation for dielectric constant (LEPSILON=.TRUE.)
  - [x] Extract epsilon from OUTCAR (arithmetic mean of tensor diagonal)
  - [x] FukuiGrid electrodes integration (`compute_fukui_potential=True`)
  - [x] Output: LOCPOT_FUKUI.vasp as SinglefileData
  - [x] Output: dielectric_constant as Float

- [x] Phase 4: Perturbative expansion model
  - [x] Retrieve LOCPOT from neutral (delta_n=0.0) calculation
  - [x] `extract_locpot_from_retrieved` calcfunction
  - [x] `run_perturbative_expansion_calcfunc` using FukuiGrid.Perturbative_point()
  - [x] Integration via `compute_perturbative_expansion=True` parameter
  - [x] Output: LOCPOT (electrostatic potential) as SinglefileData
  - [x] Output: MODELPOT_LOCPOT.vasp (interaction energy map) as SinglefileData

## Pending

- [ ] Phase 3: Fukui potential via SCPC (requires VASP with SCPC patch)
- [ ] Phase 5: Integration with main workflow (`workflow_preset='fukui_analysis'`)

---

# Implementation Notes

## Phase 1 Implementation (Completed)

FukuiGrid interpolation was integrated using **aiida-shell**:

- **Wrapper script**: `teros/core/fukui/scripts/fukui_interpolation_wrapper.py`
- **FukuiGrid location**: `teros/external/FukuiGrid/` (cloned, gitignored)
- **Integration**: `build_fukui_workgraph(..., compute_fukui=True)`

### Data Flow

```
structure → VaspTask×4 (parallel) → collect_chgcar_files → FolderData
                                                               │
                                          compute_fukui=True   │
                                                               ▼
                                    ┌──────────────────────────────────┐
                                    │ FukuiInterpolationTask           │
                                    │   - aiida-shell ShellJob         │
                                    │   - Calls FukuiGrid.py           │
                                    │   - Returns CHGCAR_FUKUI.vasp    │
                                    └──────────────────────────────────┘
                                                  │
                                                  ▼
                                         SinglefileData
                                      (CHGCAR_FUKUI.vasp)
```

### Usage

```python
from teros.core.fukui import build_fukui_workgraph, run_fukui_interpolation

# Option 1: Integrated workflow
wg = build_fukui_workgraph(
    structure=my_structure,
    nelect_neutral=312,
    code_label='VASP-6.5.1@cluster',
    builder_inputs={...},
    fukui_type='plus',
    compute_fukui=True,  # <-- Enable FukuiGrid interpolation
)
wg.submit()

# Option 2: Standalone (post-hoc analysis)
from aiida import orm
chgcar_folder = orm.load_node(31972)  # FolderData from previous run
result = run_fukui_interpolation(
    chgcar_folder,
    delta_n_values=[0.0, 0.05, 0.10, 0.15],
    fukui_type='plus',
)
```

### Visualization

Open `CHGCAR_FUKUI.vasp` in VESTA with isovalue ~0.0015 a₀⁻³

## Phase 2 Implementation (Completed)

Fukui potential via electrodes method was integrated:

- **Wrapper script**: `teros/core/fukui/scripts/fukui_electrodes_wrapper.py`
- **Integration**: `build_fukui_workgraph(..., compute_fukui_potential=True, bulk_structure=...)`

### Data Flow

```
                     ┌─────────────────────────────────┐
                     │      bulk_structure             │
                     └────────────────┬────────────────┘
                                      │
                                      ▼
                     ┌─────────────────────────────────┐
                     │   DFPT Calculation              │
                     │   (LEPSILON=.TRUE.)             │
                     └────────────────┬────────────────┘
                                      │
                                      ▼
                     ┌─────────────────────────────────┐
                     │   extract_dielectric_constant   │
                     │   ε = (εxx + εyy + εzz) / 3     │
                     └────────────────┬────────────────┘
                                      │
 ┌─────────────────────────────────┐  │
 │ Phase 1: FukuiCalculationScatter│  │
 │ → CHGCAR_0.00, CHGCAR_FUKUI.vasp│  │
 └────────────────┬────────────────┘  │
                  │                   │
                  ▼                   ▼
         ┌─────────────────────────────────────┐
         │   run_fukui_electrodes_calcfunc     │
         │   FukuiGrid.fukui_electrodes()      │
         └────────────────┬────────────────────┘
                          │
                          ▼
         ┌─────────────────────────────────────┐
         │   Output: LOCPOT_FUKUI.vasp         │
         │   (Fukui potential file)            │
         └─────────────────────────────────────┘
```

### Usage

```python
from teros.core.fukui import build_fukui_workgraph

# Phase 1 + Phase 2 integrated workflow
wg = build_fukui_workgraph(
    structure=slab_structure,
    nelect_neutral=312,
    code_label='VASP-6.5.1@cluster',
    builder_inputs={
        'parameters': {'incar': {'encut': 520, 'ediff': 1e-6}},
        'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4}},
        'potential_family': 'PBE',
        'potential_mapping': {'Sn': 'Sn_d', 'O': 'O'},
    },
    fukui_type='plus',
    compute_fukui=True,               # Phase 1: interpolation
    compute_fukui_potential=True,     # Phase 2: electrodes
    bulk_structure=bulk_structure,    # Required for DFPT
)
wg.submit()

# After completion, extract results
from teros.core.fukui import get_fukui_results, print_fukui_summary
results = get_fukui_results(wg.pk)
print(f"Dielectric constant: {results['dielectric_constant']:.4f}")
print_fukui_summary(wg.pk)
```

### Visualization

Open `LOCPOT_FUKUI.vasp` in VESTA to visualize the Fukui potential

## Phase 4 Implementation (Completed)

Perturbative expansion model for interaction energy prediction:

- **Wrapper script**: `teros/core/fukui/scripts/perturbative_expansion_wrapper.py`
- **Integration**: `build_fukui_workgraph(..., compute_perturbative_expansion=True, probe_charge=..., electron_transfer=...)`

### Theory

The c-DFT perturbative expansion predicts the interaction energy between a
charged adsorbate and the surface:

```
ΔU(r) = q·Φ(r) - q·ΔN·vf±(r)
```

Where:
- `Φ(r)` = electrostatic potential (LOCPOT from neutral slab)
- `vf±(r)` = Fukui potential (LOCPOT_FUKUI.vasp from Phase 2)
- `q` = charge of the probe (point charge model)
- `ΔN` = electron transfer (positive = oxidation, negative = reduction)

### Data Flow

```
                     ┌─────────────────────────────────┐
                     │      Phase 1 + Phase 2          │
                     │  → LOCPOT (from delta_n=0.0)    │
                     │  → LOCPOT_FUKUI.vasp            │
                     └────────────────┬────────────────┘
                                      │
                                      ▼
             ┌─────────────────────────────────────────────────┐
             │        run_perturbative_expansion_calcfunc      │
             │        FukuiGrid.Perturbative_point()           │
             │        Inputs:                                  │
             │        - locpot_neutral (SinglefileData)        │
             │        - fukui_potential (SinglefileData)       │
             │        - probe_charge: q (Float)                │
             │        - electron_transfer: ΔN (Float)          │
             └────────────────────┬────────────────────────────┘
                                  │
                                  ▼
             ┌─────────────────────────────────────────────────┐
             │        Output: MODELPOT_LOCPOT.vasp             │
             │        (Interaction energy map)                 │
             └─────────────────────────────────────────────────┘
```

### Usage

```python
from teros.core.fukui import build_fukui_workgraph

# Full workflow with Phase 4
wg = build_fukui_workgraph(
    structure=slab_structure,
    nelect_neutral=312,
    code_label='VASP-6.5.1@cluster',
    builder_inputs={
        'parameters': {'incar': {'encut': 520, 'ediff': 1e-6}},
        'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4}},
        'potential_family': 'PBE',
        'potential_mapping': {'Sn': 'Sn_d', 'O': 'O'},
    },
    fukui_type='plus',
    compute_fukui=True,                    # Phase 1: interpolation
    compute_fukui_potential=True,          # Phase 2: electrodes
    bulk_structure=bulk_structure,         # Required for DFPT
    compute_perturbative_expansion=True,   # Phase 4: perturbative expansion
    probe_charge=0.3,                      # Charge of adsorbate in |e|
    electron_transfer=-0.3,                # Electron donation to surface
)
wg.submit()

# Post-hoc analysis with different q/ΔN values
from teros.core.fukui import run_perturbative_expansion_calcfunc
from aiida import orm

locpot = orm.load_node(12345)     # LOCPOT from neutral
fukui_pot = orm.load_node(12346)  # LOCPOT_FUKUI.vasp

result = run_perturbative_expansion_calcfunc(
    locpot_neutral=locpot,
    fukui_potential=fukui_pot,
    probe_charge=orm.Float(0.5),
    electron_transfer=orm.Float(-0.5),
)
```

### Visualization

Open `MODELPOT_LOCPOT.vasp` in VESTA:
- Negative values (blue) = favorable adsorption sites
- Positive values (red) = unfavorable adsorption sites
