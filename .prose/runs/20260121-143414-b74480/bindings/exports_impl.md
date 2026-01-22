# exports_impl

kind: let

source:
```prose
let exports_impl = session: implementer
  prompt: "Update teros/core/convergence/__init__.py..."
  context: visualization_impl
```

---

## Complete Updated __init__.py Content

```python
"""Convergence testing module for PS-TEROS.

This module provides:
1. Automated ENCUT and k-points convergence testing using aiida-vasp
2. Slab thickness convergence testing for surface energy calculations
3. Visualization and export utilities for convergence analysis

Usage:
    # ENCUT/k-points convergence
    from teros.core.convergence import build_convergence_workgraph
    wg = build_convergence_workgraph(structure=..., code_label=..., ...)

    # After completion, visualize results
    from teros.core.convergence import print_convergence_summary, plot_convergence
    print_convergence_summary(wg)  # Formatted console output
    plot_convergence(wg, save_path='convergence.png')  # Plot curves

    # Thickness convergence
    from teros.core.convergence import build_thickness_convergence_workgraph
    wg = build_thickness_convergence_workgraph(
        bulk_structure_path='/path/to/bulk.cif',
        miller_indices=[1, 1, 1],
        layer_counts=[3, 5, 7, 9, 11],
        ...
    )

    # After thickness convergence completion, visualize results
    from teros.core.convergence import (
        print_thickness_convergence_summary,
        plot_thickness_convergence,
        export_thickness_convergence_data,
    )
    print_thickness_convergence_summary(wg)  # Formatted console output
    plot_thickness_convergence(wg, save_path='thickness_conv.png')  # Plot curve
    export_thickness_convergence_data(wg, save_path='thickness_data.json')  # Export data
"""

from .workgraph import (
    # ENCUT/k-points convergence
    build_convergence_workgraph,
    get_convergence_results,
    # Thickness convergence
    build_thickness_convergence_workgraph,
    get_thickness_convergence_results,
    extract_total_energy,
    calculate_surface_energy,
    analyze_thickness_convergence,
    relax_thickness_series,
    compute_surface_energies,
    gather_surface_energies,
)
from .slabs import (
    generate_thickness_series,
    extract_recommended_layers,
)
from .visualization import (
    # ENCUT/k-points visualization
    print_convergence_summary,
    plot_convergence,
    export_convergence_data,
    # Thickness convergence visualization
    print_thickness_convergence_summary,
    plot_thickness_convergence,
    export_thickness_convergence_data,
)

__all__ = [
    # ENCUT/k-points convergence
    'build_convergence_workgraph',
    'get_convergence_results',
    # ENCUT/k-points visualization and export
    'print_convergence_summary',
    'plot_convergence',
    'export_convergence_data',
    # Thickness convergence
    'build_thickness_convergence_workgraph',
    'get_thickness_convergence_results',
    'generate_thickness_series',
    'extract_recommended_layers',
    'extract_total_energy',
    'calculate_surface_energy',
    'analyze_thickness_convergence',
    'relax_thickness_series',
    'compute_surface_energies',
    'gather_surface_energies',
    # Thickness convergence visualization and export
    'print_thickness_convergence_summary',
    'plot_thickness_convergence',
    'export_thickness_convergence_data',
]
```

## Summary of Changes

1. **Added imports** from `visualization` module:
   - `print_thickness_convergence_summary`
   - `plot_thickness_convergence`
   - `export_thickness_convergence_data`

2. **Updated import section** with comments to separate:
   - ENCUT/k-points visualization functions
   - Thickness convergence visualization functions

3. **Updated `__all__` list** with:
   - New section for "Thickness convergence visualization and export"
   - Added all three new visualization functions
   - Reorganized comments for clarity

4. **Enhanced module docstring** with:
   - Usage example for thickness convergence visualization
   - Clear separation between ENCUT/k-points and thickness convergence workflows

The updated file maintains the existing structure while properly exposing the new thickness convergence visualization functions, making them available for import via `from teros.core.convergence import ...`.
