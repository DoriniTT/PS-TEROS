# design

kind: let

source:
```prose
let design = session: module-designer
  prompt: "Based on the analysis, design enhancements..."
  context: analysis
```

---

## Thickness Convergence Visualization Enhancement Design

### Overview

This design specification details the addition of three visualization functions to make the thickness convergence module fully standalone and consistent with the ENCUT/k-points convergence visualization patterns already established in `teros/core/convergence/visualization.py`.

### Data Structure Reference

The `get_thickness_convergence_results()` function returns:

```python
{
    'bulk_energy': float,           # Bulk total energy (eV)
    'convergence_results': {
        'miller_indices': [1, 1, 1],
        'results': {
            'layers_3': {'n_layers': 3, 'gamma_J_m2': 1.234, ...},
            'layers_5': {'n_layers': 5, 'gamma_J_m2': 1.156, ...},
            ...
        },
        'summary': {
            'thicknesses': [3, 5, 7, 9, 11],
            'surface_energies_J_m2': [1.234, 1.156, 1.089, 1.045, 1.032],
            'converged': True,
            'recommended_layers': 9,
            'max_tested_layers': 11,
            'convergence_threshold': 0.01,
        }
    },
    'recommended_layers': int,      # Recommended slab thickness
    'converged': bool,              # Whether convergence was achieved
    'surface_energies': {3: 1.234, 5: 1.156, 7: 1.089, ...}  # Simple mapping
}
```

### Function Specifications

#### 1. `print_thickness_convergence_summary()`

**Purpose:** Formatted console output of thickness convergence results with ASCII table.

**Signature:**
```python
def print_thickness_convergence_summary(
    workgraph: Union[int, str, WorkGraph]
) -> None:
```

**Pattern to follow from visualization.py:**
- Use `_load_workgraph()` helper for PK/string/WorkGraph input handling
- Create a structure formula extraction helper specific to bulk structure
- Print bordered header with title
- Display table with columns: Layers, Surface Energy (J/m^2), Delta from reference (mJ/m^2), Converged
- Show convergence status and recommendation at bottom

**Output format:**
```
======================================================================
              SLAB THICKNESS CONVERGENCE RESULTS
======================================================================
Structure: Au (bulk)
Miller indices: (1, 1, 1)
Threshold: 10.0 mJ/m^2
----------------------------------------------------------------------

Thickness Convergence:
+---------+------------------+------------------+-----------+
| Layers  | Surface Energy   | Delta E          | Converged |
|         | (J/m^2)          | (mJ/m^2)         |           |
+---------+------------------+------------------+-----------+
|       3 |           1.234  |          202.00  | x         |
|       5 |           1.156  |          124.00  | x         |
|       7 |           1.089  |           57.00  | x         |
|       9 |           1.045  |           13.00  | x         |
|      11 |           1.032  |            0.00  | check (ref)|
+---------+------------------+------------------+-----------+
check Converged at 9 layers (Delta = 13.0 mJ/m^2 < 10.0 mJ/m^2)
-> Recommended: 9 layers

======================================================================
SUMMARY: Recommended thickness = 9 layers for (1, 1, 1) surface
======================================================================
```

#### 2. `plot_thickness_convergence()`

**Purpose:** Matplotlib visualization of surface energy vs. slab thickness.

**Signature:**
```python
def plot_thickness_convergence(
    workgraph: Union[int, str, WorkGraph],
    save_path: Optional[str] = None,
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 150,
) -> 'matplotlib.figure.Figure':
```

**Pattern to follow from visualization.py:**
- Import matplotlib with try/except for optional dependency
- Use `_load_workgraph()` helper
- Single panel plot (unlike ENCUT/k-points which needs two panels)
- Plot surface energy vs. number of layers
- Add horizontal line at threshold level (relative to final value)
- Mark converged thickness with vertical line
- Include legend, grid, axis labels

**Visual elements:**
- X-axis: Number of atomic layers
- Y-axis: Surface energy (J/m^2)
- Data points: Circle markers with solid line
- Threshold band: Horizontal dashed red lines at reference +/- threshold
- Converged point: Vertical green dashed line at recommended layers
- Title: "Thickness Convergence: {formula} ({h}{k}{l})"

#### 3. `export_thickness_convergence_data()`

**Purpose:** Export thickness convergence data to CSV and JSON files.

**Signature:**
```python
def export_thickness_convergence_data(
    workgraph: Union[int, str, WorkGraph],
    output_dir: str,
    prefix: str = 'thickness_convergence',
) -> dict[str, str]:
```

**Pattern to follow from visualization.py:**
- Use `_load_workgraph()` helper
- Create output directory if not exists
- Return dict mapping file types to created file paths

**Output files:**

1. `{prefix}_surface_energies.csv`:
```csv
# Thickness Convergence Data
# Structure: Au
# Miller indices: (1, 1, 1)
# Threshold: 10.0 mJ/m^2
Layers,Surface_Energy_J_m2,Surface_Energy_eV_A2,Area_A2,E_slab_eV,N_slab,Delta_mJ_m2,Converged
3,1.234,0.0771,24.56,-123.45,12,202.00,False
5,1.156,0.0722,24.56,-205.78,20,124.00,False
...
```

2. `{prefix}_summary.json`:
```json
{
  "structure": {
    "formula": "Au",
    "bulk_atoms": 4
  },
  "miller_indices": [1, 1, 1],
  "bulk_energy_eV": -15.678,
  "convergence_threshold_J_m2": 0.01,
  "convergence_threshold_mJ_m2": 10.0,
  "convergence": {
    "achieved": true,
    "recommended_layers": 9,
    "max_tested_layers": 11,
    "final_surface_energy_J_m2": 1.032
  },
  "surface_energies": {
    "3": 1.234,
    "5": 1.156,
    "7": 1.089,
    "9": 1.045,
    "11": 1.032
  }
}
```

### Implementation Details

#### Helper Function: `_get_thickness_structure_info()`

Add a new helper function specific to thickness convergence:

```python
def _get_thickness_structure_info(workgraph: WorkGraph) -> dict:
    """Extract structure and miller info from thickness convergence WorkGraph."""
    info = {
        'formula': 'Unknown',
        'n_bulk_atoms': 0,
        'miller_indices': None,
    }

    try:
        # Get miller indices from generate_slabs task
        if 'generate_slabs' in workgraph.tasks:
            task = workgraph.tasks['generate_slabs']
            if hasattr(task.inputs, 'miller_indices'):
                miller = task.inputs.miller_indices.value
                if miller:
                    info['miller_indices'] = miller.get_list()

        # Get structure from bulk_relax output or input
        if 'bulk_relax' in workgraph.tasks:
            task = workgraph.tasks['bulk_relax']
            # Try output first (relaxed structure)
            if hasattr(task.outputs, 'structure'):
                structure = task.outputs.structure.value
                if structure:
                    info['formula'] = structure.get_formula()
                    info['n_bulk_atoms'] = len(structure.sites)
            # Fall back to input
            elif hasattr(task.inputs, 'structure'):
                structure = task.inputs.structure.value
                if structure:
                    info['formula'] = structure.get_formula()
                    info['n_bulk_atoms'] = len(structure.sites)
    except Exception:
        pass

    return info
```

### Module Updates Required

#### 1. Update `visualization.py`

Add the three new functions at the end of the file, after `export_convergence_data()`:

```python
# =============================================================================
# THICKNESS CONVERGENCE VISUALIZATION
# =============================================================================

def _get_thickness_structure_info(workgraph: WorkGraph) -> dict:
    """..."""

def print_thickness_convergence_summary(
    workgraph: Union[int, str, WorkGraph]
) -> None:
    """..."""

def plot_thickness_convergence(
    workgraph: Union[int, str, WorkGraph],
    save_path: Optional[str] = None,
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 150,
):
    """..."""

def export_thickness_convergence_data(
    workgraph: Union[int, str, WorkGraph],
    output_dir: str,
    prefix: str = 'thickness_convergence',
) -> dict[str, str]:
    """..."""
```

#### 2. Update `__init__.py`

Add exports for the new functions:

```python
from .visualization import (
    print_convergence_summary,
    plot_convergence,
    export_convergence_data,
    # Thickness convergence visualization
    print_thickness_convergence_summary,
    plot_thickness_convergence,
    export_thickness_convergence_data,
)

__all__ = [
    # ... existing exports ...
    # Thickness convergence visualization
    'print_thickness_convergence_summary',
    'plot_thickness_convergence',
    'export_thickness_convergence_data',
]
```

### Usage Examples

```python
from teros.core.convergence import (
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
    print_thickness_convergence_summary,
    plot_thickness_convergence,
    export_thickness_convergence_data,
)

# Build and run workflow
wg = build_thickness_convergence_workgraph(
    bulk_structure_path='/path/to/Au.cif',
    miller_indices=[1, 1, 1],
    layer_counts=[3, 5, 7, 9, 11],
    ...
)
wg.submit()

# After completion, visualize results
print_thickness_convergence_summary(wg)           # Console output
print_thickness_convergence_summary(12345)        # Using PK

fig = plot_thickness_convergence(wg)              # Display plot
fig = plot_thickness_convergence(wg, save_path='thickness_conv.png')  # Save

files = export_thickness_convergence_data(wg, '/path/to/output')
# Returns: {'csv': '/path/to/output/thickness_convergence_surface_energies.csv',
#           'json': '/path/to/output/thickness_convergence_summary.json'}
```

### Consistency Checklist

- [x] Use `_load_workgraph()` helper for input handling
- [x] Support PK (int), PK as string, and WorkGraph object inputs
- [x] Follow same docstring format with Args, Returns, Example sections
- [x] Use same logging patterns
- [x] Use same matplotlib style (colors, markers, grid)
- [x] Use same CSV header format with comments
- [x] Use same JSON structure for summary export
- [x] Return figure object from plot function
- [x] Return dict of created files from export function

### Testing Recommendations

1. Unit test `_get_thickness_structure_info()` with mock WorkGraph
2. Test `print_thickness_convergence_summary()` output format
3. Test `plot_thickness_convergence()` returns Figure object
4. Test `export_thickness_convergence_data()` creates expected files
5. Integration test with completed thickness convergence WorkGraph
