# Fukui Module TODO

## Completed

- [x] Phase 1: Interpolation method - parallel VASP calculations at different charge states, collect CHGCAR files

## Pending

- [ ] Phase 2: Fukui function computation (integrate with FukuiGrid.py)
- [ ] Phase 3: Fukui potential calculation
- [ ] Phase 4: Perturbative expansion model
- [ ] Phase 5: Integration with main workflow (`workflow_preset='fukui_analysis'`)

---

# Phase 2 Implementation Plan

## Overview

Integrate FukuiGrid.py's interpolation function to compute Fukui functions from the CHGCAR files collected in Phase 1.

**Input**: `chgcar_files` FolderData (CHGCAR_0.00, CHGCAR_0.05, CHGCAR_0.10, CHGCAR_0.15)
**Output**: `CHGCAR_FUKUI.vasp` (Fukui function on charge density grid)

---

## User Decisions

- **FukuiGrid location**: Clone to `teros/external/FukuiGrid/`
- **Execution**: Run locally (CHGCAR files already in AiiDA FolderData)
- **Approach**: Wrapper script + aiida-shell (with calcfunction fallback)
- **Integration**: Add to existing `build_fukui_workgraph()` with `compute_fukui=True`

---

## Implementation Steps

### Step 1: Clone FukuiGrid Repository

```bash
mkdir -p teros/external
cd teros/external
git clone https://github.com/cacarden/FukuiGrid.git
```

Add to `.gitignore`: `teros/external/FukuiGrid/`

### Step 2: Manual Test (Phase 2a)

Test FukuiGrid directly with existing CHGCAR files (PK 18110):

```python
from aiida import orm, load_profile
from pathlib import Path
import tempfile
import os
import sys
import numpy as np

load_profile()

# Extract files from completed workflow
folder = orm.load_node(18110)  # chgcar_files FolderData

with tempfile.TemporaryDirectory() as tmpdir:
    for name in folder.list_object_names():
        content = folder.get_object_content(name, mode='rb')
        (Path(tmpdir) / name).write_bytes(content)

    # Run FukuiGrid
    os.chdir(tmpdir)
    sys.path.insert(0, '/path/to/teros/external/FukuiGrid')
    from FukuiGrid import Fukui_interpolation

    # For f+ (nucleophilic): files in order of decreasing δN
    Fukui_interpolation(
        'CHGCAR_0.15', 'CHGCAR_0.10', 'CHGCAR_0.05', 'CHGCAR_0.00',
        np.array([0.15, 0.10, 0.05, 0.0])
    )

    # Check result
    result_file = Path(tmpdir) / 'CHGCAR_FUKUI.vasp'
    print(f"Output exists: {result_file.exists()}")
    print(f"Size: {result_file.stat().st_size / 1e6:.1f} MB")
```

### Step 3: Create Calcfunction (Phase 2b - Recommended)

**File**: `teros/core/fukui/tasks.py` (add to existing file)

```python
@task.calcfunction
def compute_fukui_function(
    chgcar_files: orm.FolderData,
    delta_n_values: orm.List,
    fukui_type: orm.Str,
) -> orm.SinglefileData:
    """
    Compute Fukui function using FukuiGrid interpolation.

    Args:
        chgcar_files: FolderData with CHGCAR_0.00, CHGCAR_0.05, etc.
        delta_n_values: List of δN values [0.0, 0.05, 0.10, 0.15]
        fukui_type: 'plus' (nucleophilic f+) or 'minus' (electrophilic f-)

    Returns:
        SinglefileData containing CHGCAR_FUKUI.vasp
    """
    import tempfile
    import numpy as np
    import sys
    import os
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract CHGCAR files
        for name in chgcar_files.list_object_names():
            content = chgcar_files.get_object_content(name, mode='rb')
            (Path(tmpdir) / name).write_bytes(content)

        # Sort files by delta_n (descending for interpolation)
        dn_list = delta_n_values.get_list()
        sorted_dn = sorted(dn_list, reverse=True)
        files = [f'CHGCAR_{dn:.2f}' for dn in sorted_dn]

        # Call FukuiGrid
        fukui_path = Path(__file__).parent.parent.parent / 'external' / 'FukuiGrid'
        sys.path.insert(0, str(fukui_path))
        from FukuiGrid import Fukui_interpolation

        cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            dn = np.array(sorted_dn)
            if fukui_type.value == 'minus':
                dn = -dn  # Negative for electrophilic
            Fukui_interpolation(files[0], files[1], files[2], files[3], dn)
        finally:
            os.chdir(cwd)

        return orm.SinglefileData(file=Path(tmpdir) / 'CHGCAR_FUKUI.vasp')
```

### Step 4: Update workgraph.py

**Add parameter to `build_fukui_workgraph()`:**

```python
def build_fukui_workgraph(
    ...
    compute_fukui: bool = False,  # NEW: Enable Phase 2 computation
    ...
)
```

**Add task after collect_chgcar_files:**

```python
if compute_fukui:
    fukui_compute = wg.add_task(
        compute_fukui_function,
        name='fukui_compute',
        chgcar_files=collect_task.outputs.result,
        delta_n_values=orm.List(list=delta_n_values),
        fukui_type=orm.Str(fukui_type),
    )
    wg.outputs.fukui_chgcar = fukui_compute.outputs.result
```

### Step 5: Update Example

```python
wg = build_fukui_workgraph(
    ...
    compute_fukui=True,  # Enable Fukui function computation
)
```

---

## Data Flow (Updated)

```
Phase 1 (existing):
    structure → VaspTask×4 (parallel) → collect_chgcar_files → FolderData
                                                                   │
Phase 2 (new):                                                     │
                                                                   ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │ compute_fukui_function (@task.calcfunction)                     │
    │   - Extract CHGCAR files to temp directory                     │
    │   - Call FukuiGrid.Fukui_interpolation()                       │
    │   - Return CHGCAR_FUKUI.vasp as SinglefileData                 │
    └─────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
                            SinglefileData
                         (CHGCAR_FUKUI.vasp)
```

---

## Files to Create/Modify

| File | Action |
|------|--------|
| `teros/external/FukuiGrid/` | Clone repository |
| `.gitignore` | Add `teros/external/FukuiGrid/` |
| `teros/core/fukui/tasks.py` | Add `compute_fukui_function()` |
| `teros/core/fukui/workgraph.py` | Add `compute_fukui` parameter |
| `teros/core/fukui/__init__.py` | Export `compute_fukui_function` |
| `examples/fukui/01_Interpolation/run_fukui_plus.py` | Add `compute_fukui=True` |

---

## Verification

1. **Clone FukuiGrid**: `ls teros/external/FukuiGrid/FukuiGrid.py`

2. **Manual test** (Step 2 above) produces `CHGCAR_FUKUI.vasp`

3. **Integrated workflow**:
   ```python
   wg = build_fukui_workgraph(..., compute_fukui=True)
   wg.submit(wait=True)
   verdi process show <PK>  # Should show fukui_chgcar output
   ```

4. **Visualize result**: Open `CHGCAR_FUKUI.vasp` in VESTA

---

## Alternative: aiida-shell Approach

If calcfunction has issues (e.g., daemon import problems), use aiida-shell:

1. Create `fukui_wrapper.py` CLI script
2. Configure `python@localhost` code
3. Use `ShellJob` task in WorkGraph

See full aiida-shell plan in `.claude/plans/` if needed.
