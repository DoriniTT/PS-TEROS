#!/usr/bin/env python
"""
Fukui Potential Analysis Script

Analyzes Phase 1 (interpolation) + Phase 2 (electrodes) results
from the completed Fukui potential workflow.

This script:
1. Extracts all results from the completed WorkGraph
2. Exports CHGCAR_FUKUI.vasp and LOCPOT_FUKUI.vasp for VESTA visualization
3. Creates visualization if planar average data is available
4. Provides a summary with interpretation

Usage:
    python analyze_fukui_potential.py [WORKGRAPH_PK]

If no PK is provided, reads from ../last_pk.txt
"""

import os
import sys
from pathlib import Path

# Default WorkGraph PK
DEFAULT_PK = 37028


def load_aiida():
    """Load AiiDA profile."""
    from aiida import load_profile
    load_profile('presto')


def get_workgraph_pk():
    """Get WorkGraph PK from command line or last_pk.txt file."""
    if len(sys.argv) > 1:
        return int(sys.argv[1])

    # Try to read from last_pk.txt
    last_pk_file = Path(__file__).parent.parent / 'last_pk.txt'
    if last_pk_file.exists():
        with open(last_pk_file) as f:
            content = f.read().strip()
            if content.isdigit():
                return int(content)

    return DEFAULT_PK


def extract_results(pk):
    """Extract results from completed Fukui WorkGraph."""
    from teros.core.fukui import get_fukui_results, print_fukui_summary

    print(f"\n{'='*60}")
    print(f"ANALYZING FUKUI WORKGRAPH PK={pk}")
    print('='*60)

    # Print formatted summary
    print_fukui_summary(pk)

    # Get results dict for further processing
    results = get_fukui_results(pk)

    return results


def export_fukui_chgcar(results, output_dir):
    """Export Fukui function CHGCAR file."""
    fukui_file = results.get('fukui_chgcar')

    if fukui_file is None:
        print("\nFukui CHGCAR not available (compute_fukui was not enabled)")
        return None

    output_path = output_dir / 'CHGCAR_FUKUI.vasp'
    content = fukui_file.get_content()

    mode = 'wb' if isinstance(content, bytes) else 'w'
    with open(output_path, mode) as f:
        f.write(content)

    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"\nExported: {output_path.name} ({size_mb:.1f} MB)")

    return output_path


def export_fukui_potential(results, output_dir):
    """Export Fukui potential LOCPOT file."""
    potential_file = results.get('fukui_potential')

    if potential_file is None:
        print("\nFukui potential not available (compute_fukui_potential was not enabled)")
        return None

    output_path = output_dir / 'LOCPOT_FUKUI.vasp'
    content = potential_file.get_content()

    mode = 'wb' if isinstance(content, bytes) else 'w'
    with open(output_path, mode) as f:
        f.write(content)

    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"Exported: {output_path.name} ({size_mb:.1f} MB)")

    return output_path


def compute_planar_average(volumetric_data, axis=2):
    """
    Compute planar average of volumetric data along specified axis.

    Args:
        volumetric_data: pymatgen VolumetricData object (Chgcar or Locpot)
        axis: Axis to average along (0=a, 1=b, 2=c/z)

    Returns:
        tuple: (z_coordinates in Angstrom, planar_average values)
    """
    import numpy as np

    # Get the data grid
    data = volumetric_data.data['total']

    # Get lattice vector length along the axis
    lattice = volumetric_data.structure.lattice
    axis_length = lattice.abc[axis]

    # Number of grid points along the axis
    n_points = data.shape[axis]

    # Compute planar average by averaging over the other two axes
    axes_to_average = tuple(i for i in range(3) if i != axis)
    planar_avg = np.mean(data, axis=axes_to_average)

    # Create z-coordinates array
    z_coords = np.linspace(0, axis_length, n_points, endpoint=False)

    return z_coords, planar_avg


def create_planar_average_plot(output_dir):
    """
    Create visualization of planar averages for both Fukui function and potential.

    Reads the exported CHGCAR_FUKUI.vasp and LOCPOT_FUKUI.vasp files and
    computes planar averages along the z-axis (perpendicular to surface).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    chgcar_path = output_dir / 'CHGCAR_FUKUI.vasp'
    locpot_path = output_dir / 'LOCPOT_FUKUI.vasp'

    has_chgcar = chgcar_path.exists()
    has_locpot = locpot_path.exists()

    if not has_chgcar and not has_locpot:
        print("\nNo files available for planar average plot")
        return None

    # Create figure with subplots
    n_plots = sum([has_chgcar, has_locpot])
    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5))

    if n_plots == 1:
        axes = [axes]

    plot_idx = 0

    # Plot Fukui function (CHGCAR)
    if has_chgcar:
        print("\nComputing planar average for CHGCAR_FUKUI.vasp...")
        from pymatgen.io.vasp import Chgcar

        chgcar = Chgcar.from_file(str(chgcar_path))
        z_coords, z_avg = compute_planar_average(chgcar, axis=2)

        ax = axes[plot_idx]
        ax.plot(z_coords, z_avg, 'b-', linewidth=1.5)
        ax.set_xlabel(r'z ($\AA$)', fontsize=12)
        ax.set_ylabel(r'Fukui function $f^+(r)$ (e/$\AA^3$)', fontsize=12)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.set_title('Fukui Function Planar Average')
        ax.grid(True, alpha=0.3)

        # Add statistics
        ax.text(0.02, 0.98, f'Max: {z_avg.max():.4f}\nMin: {z_avg.min():.4f}',
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plot_idx += 1

    # Plot Fukui potential (LOCPOT)
    if has_locpot:
        print("Computing planar average for LOCPOT_FUKUI.vasp...")
        from pymatgen.io.vasp import Locpot

        locpot = Locpot.from_file(str(locpot_path))
        z_coords, z_avg = compute_planar_average(locpot, axis=2)

        ax = axes[plot_idx]
        ax.plot(z_coords, z_avg, 'r-', linewidth=1.5)
        ax.set_xlabel(r'z ($\AA$)', fontsize=12)
        ax.set_ylabel(r'Fukui potential $v_f^+(r)$ (eV)', fontsize=12)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.set_title('Fukui Potential Planar Average')
        ax.grid(True, alpha=0.3)

        # Add statistics
        ax.text(0.02, 0.98, f'Max: {z_avg.max():.4f}\nMin: {z_avg.min():.4f}',
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_path = output_dir / 'fukui_planar_average.png'
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"Saved: {output_path.name}")

    return output_path


def print_interpretation(results):
    """Print interpretation of results."""
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)

    # Dielectric constant interpretation
    epsilon = results.get('dielectric_constant')
    if epsilon is not None:
        print(f"\n1. Dielectric Constant: {epsilon:.4f}")
        print("   - Measures the material's ability to screen electric fields")
        print("   - Used to correct the Fukui potential for electrode effects")
        if 7.0 <= epsilon <= 10.0:
            print("   - Value is reasonable for SnO2 (typical range: 7-10)")
        elif epsilon < 7.0:
            print("   - Value seems low for SnO2 - verify calculation converged")
        else:
            print("   - Value seems high for SnO2 - check DFPT settings")

    # Fukui function interpretation
    if results.get('fukui_chgcar') is not None:
        print("\n2. Fukui Function (CHGCAR_FUKUI.vasp):")
        print("   - Shows electron density change upon reduction: f+(r) = drho/dN")
        print("   - High positive values: sites susceptible to nucleophilic attack")
        print("   - Visualize in VESTA with isosurface value ~0.001-0.01 e/A^3")

    # Fukui potential interpretation
    if results.get('fukui_potential') is not None:
        print("\n3. Fukui Potential (LOCPOT_FUKUI.vasp):")
        print("   - Electrostatic potential weighted by Fukui function")
        print("   - Accounts for electrode screening effects")
        print("   - Negative values: favorable for cation adsorption")
        print("   - Positive values: favorable for anion adsorption")

    # Next steps
    print("\n" + "="*60)
    print("NEXT STEPS")
    print("="*60)
    print("""
1. Open CHGCAR_FUKUI.vasp in VESTA:
   - Boundary > Show atoms outside boundary
   - Properties > Isosurfaces > set level to 0.005

2. Open LOCPOT_FUKUI.vasp in VESTA:
   - Use same settings, adjust isosurface for best visualization

3. Consider Phase 4 (perturbative expansion):
   - Predicts interaction energy for specific adsorbates
   - Set compute_perturbative_expansion=True with probe_charge and electron_transfer

4. Compare reactive sites with experimental data if available
""")


def main():
    """Main analysis workflow."""
    # Load AiiDA
    load_aiida()

    # Get WorkGraph PK
    pk = get_workgraph_pk()

    # Set output directory
    output_dir = Path(__file__).parent
    os.chdir(output_dir)

    # Extract results
    results = extract_results(pk)

    # Export files
    print("\n" + "-"*60)
    print("EXPORTING FILES")
    print("-"*60)

    export_fukui_chgcar(results, output_dir)
    export_fukui_potential(results, output_dir)

    # Create planar average plot from exported files
    print("\n" + "-"*60)
    print("CREATING PLANAR AVERAGE PLOT")
    print("-"*60)
    create_planar_average_plot(output_dir)

    # Print interpretation
    print_interpretation(results)

    # List exported files
    print("\n" + "="*60)
    print("EXPORTED FILES")
    print("="*60)
    for f in output_dir.glob('*'):
        if f.suffix in ['.vasp', '.png']:
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"  {f.name}: {size_mb:.2f} MB")


if __name__ == '__main__':
    main()
