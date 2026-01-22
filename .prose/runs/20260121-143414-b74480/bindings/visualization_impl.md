# visualization_impl

kind: let

source:
```prose
let visualization_impl = session: implementer
  prompt: "Implement the thickness convergence visualization functions..."
  context: design
```

---

## Complete Python code to ADD to visualization.py

Add the following functions after the existing `export_convergence_data()` function (after line 437):

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

## Integration Notes

1. **Import placement**: These functions should be added after the existing `export_convergence_data()` function (line 437)

2. **No new imports needed**: All required imports are already present at the top of the file:
   - `json`, `logging`, `Path`, `Optional`, `Union` (stdlib)
   - `orm`, `WorkGraph` (AiiDA)
   - The function imports `get_thickness_convergence_results` locally to avoid circular imports

3. **Consistent patterns**:
   - Uses `_load_workgraph()` helper (already exists)
   - Follows exact table formatting style from `print_convergence_summary()`
   - Matches plot styling from `plot_convergence()`
   - Uses same CSV/JSON export patterns from `export_convergence_data()`

4. **Error handling**: Functions raise clear errors if data is missing, matching existing code style

5. **Docstrings**: Google-style docstrings with Args, Returns, and Example sections

6. **Type hints**: All parameters and returns properly typed

## Testing Checklist

After adding this code:

```bash
# 1. Restart daemon (CRITICAL)
verdi daemon restart

# 2. Check imports
python -c "from teros.core.convergence import print_thickness_convergence_summary, plot_thickness_convergence, export_thickness_convergence_data"

# 3. Linting
flake8 teros/core/convergence/visualization.py --max-line-length=120 --ignore=E501,W503,E402,F401

# 4. Test with a completed WorkGraph
python -c "
from teros.core.convergence import print_thickness_convergence_summary
print_thickness_convergence_summary(<PK>)
"
```
