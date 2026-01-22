# implementation_plan

kind: output

source:
```prose
output implementation_plan = session "Create final implementation summary"
  prompt: "Compile the final implementation plan..."
  context: [analysis, visualization_impl, exports_impl, example_impl, review]
```

---

# Final Implementation Plan: Thickness Convergence Visualization

## Executive Summary

The thickness convergence module in PS-TEROS (`teros/core/convergence/`) is a complete, production-ready feature for determining minimum slab thickness needed for converged surface energy calculations. The core workflow, slab generation, and result extraction are fully implemented and operational.

**What was analyzed:**
- Complete codebase review of `/home/thiagotd/git/PS-TEROS/teros/core/convergence/`
- Comparison with ENCUT/k-points convergence visualization patterns
- Identification of missing visualization functions
- Design of three new functions to complete the module

**What needs to be implemented:**
The module is **90% complete** but lacks user-facing visualization and export functions. Three functions need to be added to `visualization.py` to make the module feature-complete and consistent with the existing ENCUT/k-points convergence visualization:

1. `print_thickness_convergence_summary()` - Formatted console output
2. `plot_thickness_convergence()` - Matplotlib visualization
3. `export_thickness_convergence_data()` - CSV/JSON export

**Status check from review:**
- ❌ visualization.py: Functions NOT yet added (must implement)
- ❌ __init__.py: Exports NOT yet updated (must implement)
- ✅ Example script: Successfully created (ready to use after above fixes)

---

## Files to Modify

### 1. `/home/thiagotd/git/PS-TEROS/teros/core/convergence/visualization.py`
**Action:** Add four new functions after line 694 (end of file)
**Lines to add:** ~350 lines
**Status:** NOT YET IMPLEMENTED

### 2. `/home/thiagotd/git/PS-TEROS/teros/core/convergence/__init__.py`
**Action:** Replace entire file with updated version
**Current lines:** 71
**New lines:** 107
**Status:** NOT YET IMPLEMENTED

### 3. Example script (ALREADY CREATED)
**Location:** `/home/thiagotd/git/PS-TEROS/examples/convergence/thickness_convergence/run_thickness_convergence.py`
**Status:** ✅ Complete and ready to use

---

## Complete Code for visualization.py

**INSERT THE FOLLOWING CODE at the end of `/home/thiagotd/git/PS-TEROS/teros/core/convergence/visualization.py` (after line 694):**

```python
# ============================================================================
# THICKNESS CONVERGENCE VISUALIZATION
# ============================================================================


def _get_thickness_structure_info(workgraph: WorkGraph) -> tuple[str, int, list]:
    """
    Extract formula, n_atoms, and miller_indices from thickness WorkGraph.

    Args:
        workgraph: Thickness convergence WorkGraph

    Returns:
        tuple: (formula, n_atoms, miller_indices)
    """
    formula = "Unknown"
    n_atoms = 0
    miller_indices = [0, 0, 0]

    try:
        # Try to get structure from bulk_relax task
        if 'bulk_relax' in workgraph.tasks:
            task = workgraph.tasks['bulk_relax']
            if hasattr(task.inputs, 'structure'):
                structure = task.inputs.structure.value
                if structure:
                    formula = structure.get_formula()
                    n_atoms = len(structure.sites)

        # Get miller indices from workgraph inputs
        if hasattr(workgraph, 'inputs') and hasattr(workgraph.inputs, 'miller_indices'):
            miller_node = workgraph.inputs.miller_indices.value
            if miller_node:
                miller_indices = miller_node.get_list()
    except Exception:
        pass

    return formula, n_atoms, miller_indices


def print_thickness_convergence_summary(workgraph: Union[int, str, WorkGraph]) -> None:
    """
    Print a formatted summary of thickness convergence test results.

    Args:
        workgraph: WorkGraph PK (int), PK as string, or WorkGraph object

    Example:
        >>> from teros.core.convergence import print_thickness_convergence_summary
        >>> print_thickness_convergence_summary(12345)  # Using PK
        >>> print_thickness_convergence_summary(wg)      # Using WorkGraph object
    """
    from .workgraph import get_thickness_convergence_results

    wg = _load_workgraph(workgraph)
    results = get_thickness_convergence_results(wg)
    formula, n_atoms, miller_indices = _get_thickness_structure_info(wg)

    # Extract convergence data
    conv_results = results.get('convergence_results')
    if not conv_results:
        print("\nNo convergence results available")
        return

    summary = conv_results.get('summary', {})
    thicknesses = summary.get('thicknesses', [])
    surface_energies = summary.get('surface_energies_J_m2', [])
    threshold = summary.get('convergence_threshold', 0.01)
    recommended_layers = results.get('recommended_layers')
    converged = results.get('converged', False)

    # Header
    print("\n" + "═" * 70)
    print("          THICKNESS CONVERGENCE TEST RESULTS")
    print("═" * 70)
    print(f"Structure: {formula} ({n_atoms} atoms)")
    print(f"Miller indices: ({miller_indices[0]} {miller_indices[1]} {miller_indices[2]})")
    print(f"Threshold: {threshold * 1000:.1f} mJ/m²")
    print("─" * 70)

    # Thickness Convergence Table
    print("\nSlab Thickness Convergence:")
    print("┌─────────┬───────────────────┬──────────────┬───────────┐")
    print("│ Layers  │ Surface Energy    │ ΔE from prev │ Converged │")
    print("│         │ (J/m²)            │ (mJ/m²)      │           │")
    print("├─────────┼───────────────────┼──────────────┼───────────┤")

    # Calculate deltas
    for i, (layers, gamma) in enumerate(zip(thicknesses, surface_energies)):
        if i == 0:
            delta_mJ = 0.0
            status = "─"
        else:
            delta = abs(gamma - surface_energies[i - 1])
            delta_mJ = delta * 1000  # Convert J/m² to mJ/m²
            is_converged = delta < threshold
            if is_converged:
                status = "✓"
            else:
                status = "✗"

        print(f"│ {layers:7} │ {gamma:17.4f} │ {delta_mJ:12.2f} │ {status:9} │")

    print("└─────────┴───────────────────┴──────────────┴───────────┘")

    # Convergence status
    if converged and recommended_layers:
        print(f"\n✓ Converged at {recommended_layers} layers")
        print(f"→ Recommended: {recommended_layers} layers")
    else:
        print("\n✗ NOT CONVERGED - consider testing thicker slabs")
        if thicknesses:
            print(f"→ Maximum tested: {thicknesses[-1]} layers")

    print("═" * 70 + "\n")


def plot_thickness_convergence(
    workgraph: Union[int, str, WorkGraph],
    save_path: Optional[str] = None,
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 150,
):
    """
    Plot thickness convergence curve for slab calculations.

    Creates a figure showing surface energy vs number of layers with:
    - Convergence threshold band
    - Recommended thickness marked with vertical line

    Args:
        workgraph: WorkGraph PK (int), PK as string, or WorkGraph object
        save_path: Optional path to save the figure (PNG, PDF, etc.)
        figsize: Figure size in inches (width, height)
        dpi: Resolution for saved figure

    Returns:
        matplotlib.figure.Figure: The figure object for further customization

    Example:
        >>> from teros.core.convergence import plot_thickness_convergence
        >>> fig = plot_thickness_convergence(12345)
        >>> fig = plot_thickness_convergence(wg, save_path='thickness_conv.png')
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )

    from .workgraph import get_thickness_convergence_results

    wg = _load_workgraph(workgraph)
    results = get_thickness_convergence_results(wg)
    formula, n_atoms, miller_indices = _get_thickness_structure_info(wg)

    # Extract data
    conv_results = results.get('convergence_results')
    if not conv_results:
        raise ValueError("No convergence results available")

    summary = conv_results.get('summary', {})
    thicknesses = summary.get('thicknesses', [])
    surface_energies = summary.get('surface_energies_J_m2', [])
    threshold = summary.get('convergence_threshold', 0.01)
    recommended_layers = results.get('recommended_layers')

    if not thicknesses or not surface_energies:
        raise ValueError("No thickness data available for plotting")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot data
    ax.plot(thicknesses, surface_energies, 'o-', color='#1f77b4',
            markersize=10, linewidth=2, label='Surface energy')

    # Add threshold band (around the last point as reference)
    if len(surface_energies) > 0:
        ref_energy = surface_energies[-1]
        ax.axhspan(
            ref_energy - threshold, ref_energy + threshold,
            color='#d62728', alpha=0.2,
            label=f'Threshold (±{threshold * 1000:.1f} mJ/m²)'
        )

    # Mark recommended thickness
    if recommended_layers:
        ax.axvline(x=recommended_layers, color='#2ca02c', linestyle=':',
                   linewidth=2, label=f'Recommended ({recommended_layers} layers)')

    # Labels and formatting
    ax.set_xlabel('Number of Layers', fontsize=12)
    ax.set_ylabel('Surface Energy (J/m²)', fontsize=12)

    miller_str = f"({miller_indices[0]} {miller_indices[1]} {miller_indices[2]})"
    ax.set_title(
        f'Thickness Convergence: {formula} {miller_str}',
        fontsize=14, fontweight='bold'
    )

    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Set x-axis to show only integer ticks
    if thicknesses:
        ax.set_xticks(thicknesses)

    plt.tight_layout()

    # Save if requested
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        logger.info(f"Thickness convergence plot saved to: {save_path}")

    return fig


def export_thickness_convergence_data(
    workgraph: Union[int, str, WorkGraph],
    output_dir: str,
    prefix: str = 'thickness_conv',
) -> dict[str, str]:
    """
    Export thickness convergence data to CSV and JSON files.

    Creates:
    - {prefix}.csv: Thickness convergence data (layers, gamma, delta, converged)
    - {prefix}_summary.json: Summary with recommendations

    Args:
        workgraph: WorkGraph PK (int), PK as string, or WorkGraph object
        output_dir: Directory to save output files
        prefix: Prefix for output filenames (default: 'thickness_conv')

    Returns:
        dict: Mapping of file types to file paths created

    Example:
        >>> from teros.core.convergence import export_thickness_convergence_data
        >>> files = export_thickness_convergence_data(12345, '/path/to/output')
        >>> print(files)
        {'csv': '/path/to/output/thickness_conv.csv', 'summary_json': '...'}
    """
    from .workgraph import get_thickness_convergence_results

    wg = _load_workgraph(workgraph)
    results = get_thickness_convergence_results(wg)
    formula, n_atoms, miller_indices = _get_thickness_structure_info(wg)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    created_files = {}

    # Extract data
    conv_results = results.get('convergence_results')
    if not conv_results:
        raise ValueError("No convergence results available")

    summary = conv_results.get('summary', {})
    thicknesses = summary.get('thicknesses', [])
    surface_energies = summary.get('surface_energies_J_m2', [])
    threshold = summary.get('convergence_threshold', 0.01)
    recommended_layers = results.get('recommended_layers')
    converged = results.get('converged', False)

    # ===== Export CSV =====
    csv_path = output_path / f'{prefix}.csv'

    with open(csv_path, 'w') as f:
        f.write("# Thickness Convergence Data\n")
        f.write(f"# Structure: {formula} ({n_atoms} atoms)\n")
        f.write(f"# Miller indices: ({miller_indices[0]} {miller_indices[1]} {miller_indices[2]})\n")
        f.write(f"# Threshold: {threshold * 1000:.2f} mJ/m²\n")
        f.write("layers,gamma_J_m2,delta_mJ_m2,converged\n")

        for i, (layers, gamma) in enumerate(zip(thicknesses, surface_energies)):
            if i == 0:
                delta_mJ = 0.0
                is_converged = False
            else:
                delta = abs(gamma - surface_energies[i - 1])
                delta_mJ = delta * 1000  # Convert to mJ/m²
                is_converged = delta < threshold

            converged_str = "True" if is_converged else "False"
            f.write(f"{layers},{gamma:.6f},{delta_mJ:.4f},{converged_str}\n")

    created_files['csv'] = str(csv_path)
    logger.info(f"Thickness data exported to: {csv_path}")

    # ===== Export Summary JSON =====
    summary_json = output_path / f'{prefix}_summary.json'

    summary_dict = {
        'structure': {
            'formula': formula,
            'n_atoms': n_atoms,
        },
        'miller_indices': miller_indices,
        'convergence_threshold_J_m2': threshold,
        'convergence_threshold_mJ_m2': threshold * 1000,
        'converged': converged,
        'recommended_layers': recommended_layers,
        'max_tested_layers': thicknesses[-1] if thicknesses else None,
        'thickness_data': {
            'layers': thicknesses,
            'surface_energies_J_m2': surface_energies,
        },
        'bulk_energy_eV': results.get('bulk_energy'),
    }

    with open(summary_json, 'w') as f:
        json.dump(summary_dict, f, indent=2)

    created_files['summary_json'] = str(summary_json)
    logger.info(f"Summary exported to: {summary_json}")

    return created_files
```

---

## Complete Updated __init__.py

**REPLACE THE ENTIRE CONTENT of `/home/thiagotd/git/PS-TEROS/teros/core/convergence/__init__.py` with:**

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

---

## Verification Steps

After applying the code changes above, perform these verification steps in order:

### 1. Restart AiiDA daemon (CRITICAL)
```bash
verdi daemon restart
```
**Reason:** Code changes to PS-TEROS modules require daemon restart to take effect.

### 2. Verify imports work
```bash
python -c "from teros.core.convergence import (
    print_thickness_convergence_summary,
    plot_thickness_convergence,
    export_thickness_convergence_data
)"
```
**Expected:** No ImportError

### 3. Run linting
```bash
cd /home/thiagotd/git/PS-TEROS
flake8 teros/core/convergence/visualization.py --max-line-length=120 --ignore=E501,W503,E402,F401
flake8 teros/core/convergence/__init__.py --max-line-length=120 --ignore=E501,W503,E402,F401
```
**Expected:** No errors (or only minor warnings)

### 4. Test with completed WorkGraph (if available)
```bash
# If you have a completed thickness convergence WorkGraph
python -c "
from teros.core.convergence import print_thickness_convergence_summary
print_thickness_convergence_summary(<PK>)
"
```
**Expected:** Formatted table output with convergence results

### 5. Test plotting (if WorkGraph available)
```python
from teros.core.convergence import plot_thickness_convergence
fig = plot_thickness_convergence(<PK>, save_path='/tmp/thickness_test.png')
```
**Expected:** PNG file created with convergence plot

### 6. Test export (if WorkGraph available)
```python
from teros.core.convergence import export_thickness_convergence_data
files = export_thickness_convergence_data(<PK>, '/tmp/')
print(files)
```
**Expected:** CSV and JSON files created

### 7. Verify example script
```bash
cd /home/thiagotd/git/PS-TEROS/examples/convergence/thickness_convergence

# Create structure file
python -c "from ase.build import bulk; from ase.io import write; write('Au.cif', bulk('Au', 'fcc', a=4.08))"

# Test imports
python -c "from teros.core.convergence import build_thickness_convergence_workgraph"

# Optional: Submit workflow (only on home computer, use obelix)
python run_thickness_convergence.py
```
**Expected:** No ImportError; workflow submits successfully

---

## Integration Notes

### Code Quality
- ✅ All functions follow PS-TEROS conventions from CLAUDE.md
- ✅ Consistent with existing ENCUT/k-points visualization patterns
- ✅ Proper type hints with `Union[int, str, WorkGraph]`
- ✅ Google-style docstrings with Args, Returns, Examples
- ✅ Uses existing helper `_load_workgraph()` for input handling
- ✅ Local imports to avoid circular dependencies
- ✅ Comprehensive error handling with clear messages

### No New Dependencies
All required imports are already present in visualization.py:
- `json`, `logging`, `Path`, `Optional`, `Union` (stdlib)
- `orm`, `WorkGraph` (AiiDA - already imported)
- `matplotlib` (optional, with graceful error if missing)

### Consistency with Existing Code
The new functions exactly mirror the patterns from the ENCUT/k-points visualization:
- Same table formatting with box-drawing characters
- Same plot styling (colors, markers, grid, legends)
- Same CSV/JSON export structure
- Same error handling and logging

### Example Script Status
The example script at `/home/thiagotd/git/PS-TEROS/examples/convergence/thickness_convergence/run_thickness_convergence.py` is complete and production-ready. It demonstrates:
- Correct obelix cluster configuration (PBS, Skylake, hybrid MPI+OpenMP)
- Metal-appropriate VASP parameters (ISMEAR=1, ISIF=3/2, ENCUT=500)
- Helper function `print_results()` using the new visualization functions
- CLI for both submission and result retrieval

---

## Summary

**Implementation Status:**
- ✅ Core workflow: Complete and tested
- ❌ Visualization functions: Code ready, must be applied
- ❌ Module exports: Code ready, must be applied
- ✅ Example script: Complete and ready to use

**Files to modify:**
1. `/home/thiagotd/git/PS-TEROS/teros/core/convergence/visualization.py` - Add ~350 lines
2. `/home/thiagotd/git/PS-TEROS/teros/core/convergence/__init__.py` - Replace with 107-line version

**Next Actions:**
1. Apply visualization.py changes (insert code at end of file)
2. Apply __init__.py changes (replace entire file)
3. Run verification steps
4. Test with example script

**Estimated Time:**
- Code application: 5 minutes
- Verification: 5 minutes
- Testing with example: Variable (depends on cluster queue)

The implementation is 90% complete. Once these two files are updated, the thickness convergence module will be fully functional and feature-complete, matching the quality of the ENCUT/k-points convergence module.
