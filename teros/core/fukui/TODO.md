# Fukui Module TODO

## Completed

- [x] Phase 1: Interpolation method
  - [x] Parallel VASP calculations at different charge states
  - [x] Collect CHGCAR files into FolderData
  - [x] FukuiGrid integration (`compute_fukui=True` parameter)
  - [x] aiida-shell based execution for FukuiGrid interpolation
  - [x] Output: CHGCAR_FUKUI.vasp as SinglefileData

## Pending

- [ ] Phase 2: Fukui potential calculation (FukuiGrid electrodes method)
- [ ] Phase 3: Fukui potential via SCPC
- [ ] Phase 4: Perturbative expansion model
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
