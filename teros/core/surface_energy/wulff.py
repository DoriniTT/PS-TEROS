"""
Wulff Shape Generation for Metal and Intermetallic Surface Energies.

This module generates Wulff shapes from computed surface energies using
pymatgen's WulffShape class. It produces:
- A PDF visualization of the Wulff shape
- A text report with all numerical data
- A standalone Python script for customizing the plot

The Wulff construction determines the equilibrium crystal shape by minimizing
total surface energy for a given volume. Facets with lower surface energy
occupy larger fractions of the crystal surface.

References:
- G. Wulff, Z. Kristallogr. 34, 449 (1901)
- https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.wulff.WulffShape
"""

from __future__ import annotations

import io
import tempfile
from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
from aiida import orm
from aiida_workgraph import task

if TYPE_CHECKING:
    from pymatgen.analysis.wulff import WulffShape


def hkl_key_to_tuple(hkl_key: str) -> tuple[int, int, int]:
    """
    Convert hkl key string to Miller indices tuple.

    Args:
        hkl_key: String like 'hkl_111' or 'hkl_110'

    Returns:
        Tuple like (1, 1, 1) or (1, 1, 0)
    """
    # Extract digits after 'hkl_'
    digits = hkl_key.replace('hkl_', '')
    return tuple(int(d) for d in digits)


def select_lowest_energy_terminations(surface_energies: dict) -> dict:
    """
    Select the lowest energy termination for each Miller index.

    Args:
        surface_energies: Dictionary with structure:
            {
                'hkl_111': {
                    'term_0': {'gamma_J_m2': 0.79, ...},
                    'term_1': {'gamma_J_m2': 0.85, ...}
                },
                ...
            }

    Returns:
        Dictionary with structure:
            {
                (1, 1, 1): {
                    'gamma_J_m2': 0.79,
                    'termination': 'term_0',
                    'full_data': {...}
                },
                ...
            }
    """
    selected = {}

    for hkl_key, terminations in surface_energies.items():
        miller = hkl_key_to_tuple(hkl_key)

        # Find termination with lowest surface energy
        lowest_energy = float('inf')
        best_term = None
        best_data = None

        for term_key, term_data in terminations.items():
            gamma = term_data.get('gamma_J_m2', float('inf'))
            if gamma < lowest_energy:
                lowest_energy = gamma
                best_term = term_key
                best_data = term_data

        if best_data is not None:
            selected[miller] = {
                'gamma_J_m2': lowest_energy,
                'termination': best_term,
                'full_data': best_data,
            }

    return selected


def create_wulff_shape(
    lattice,
    miller_list: list[tuple],
    e_surf_list: list[float],
) -> 'WulffShape':
    """
    Create a WulffShape object from surface energies.

    Args:
        lattice: Pymatgen Lattice object
        miller_list: List of Miller indices as tuples, e.g., [(1,1,1), (1,0,0)]
        e_surf_list: List of surface energies in J/m² (same order as miller_list)

    Returns:
        WulffShape object
    """
    from pymatgen.analysis.wulff import WulffShape

    return WulffShape(lattice, miller_list, e_surf_list)


def generate_wulff_plot_pdf(
    wulff: 'WulffShape',
    formula: str,
    output_path: str,
    figsize: tuple = (10, 8),
    dpi: int = 150,
) -> None:
    """
    Generate a PDF plot of the Wulff shape.

    Args:
        wulff: WulffShape object
        formula: Chemical formula for title
        output_path: Path to save the PDF
        figsize: Figure size in inches
        dpi: Resolution
    """
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt

    # Get the plot from pymatgen
    fig = wulff.get_plot(
        color_set='RdBu',
        grid_off=True,
        axis_off=True,
        show_area=True,
        alpha=0.9,
        off_color='white',
        direction=(1, 1, 1),
    )

    # Add title
    fig.suptitle(
        f'Wulff Shape: {formula}',
        fontsize=14,
        fontweight='bold',
        y=0.95,
    )

    # Add weighted surface energy annotation
    weighted_se = wulff.weighted_surface_energy
    fig.text(
        0.5, 0.02,
        f'Weighted Surface Energy: {weighted_se:.3f} J/m²',
        ha='center',
        fontsize=11,
    )

    # Save
    fig.savefig(output_path, format='pdf', dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def generate_wulff_report_txt(
    selected_terminations: dict,
    wulff: 'WulffShape',
    formula: str,
    lattice_param: float,
    output_path: str,
) -> None:
    """
    Generate a text report with all Wulff shape data.

    Args:
        selected_terminations: Dict from select_lowest_energy_terminations()
        wulff: WulffShape object
        formula: Chemical formula
        lattice_param: Lattice parameter in Angstroms
        output_path: Path to save the text file
    """
    lines = []
    sep = '=' * 80
    lines.append(sep)
    lines.append('                         WULFF SHAPE ANALYSIS REPORT')
    lines.append(sep)
    lines.append(f'Material: {formula}')
    lines.append(f'Lattice Parameter: {lattice_param:.4f} Å')
    lines.append(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    lines.append('')

    # Surface energies table
    lines.append('SURFACE ENERGIES (selected lowest-energy terminations):')
    lines.append('-' * 80)
    lines.append(f'  {"Miller Index":<15} {"γ (J/m²)":<12} {"Area Fraction (%)":<20} {"Termination Used"}')
    lines.append('-' * 80)

    # Sort by surface energy
    sorted_millers = sorted(
        selected_terminations.keys(),
        key=lambda m: selected_terminations[m]['gamma_J_m2']
    )

    # Get area fractions from Wulff shape
    area_fractions = wulff.area_fraction_dict

    for miller in sorted_millers:
        data = selected_terminations[miller]
        gamma = data['gamma_J_m2']
        term = data['termination']

        # Find matching area fraction (miller indices might be in different tuple format)
        area_frac = 0.0
        for wulff_miller, frac in area_fractions.items():
            if tuple(wulff_miller) == miller:
                area_frac = frac * 100  # Convert to percentage
                break

        miller_str = str(miller)
        lines.append(f'  {miller_str:<15} {gamma:<12.4f} {area_frac:<20.2f} {term}')

    lines.append('-' * 80)
    lines.append('')

    # Wulff shape properties
    lines.append('WULFF SHAPE PROPERTIES:')
    lines.append(f'  Weighted Surface Energy:  {wulff.weighted_surface_energy:.4f} J/m²')
    lines.append(f'  Shape Anisotropy:         {wulff.anisotropy:.4f}')
    lines.append(f'  Shape Factor:             {wulff.shape_factor:.4f} (1.0 = sphere)')
    lines.append(f'  Total Surface Area:       {wulff.surface_area:.4f} (normalized units)')
    lines.append(f'  Volume:                   {wulff.volume:.6f} (normalized units)')
    lines.append(f'  Effective Radius:         {wulff.effective_radius:.4f} (normalized units)')
    lines.append('')

    # Facet details
    lines.append('FACET DETAILS:')
    lines.append('-' * 80)
    lines.append(f'  {"Miller Index":<15} {"Normal Distance":<18} {"Area":<15} {"# Facets"}')
    lines.append('-' * 80)

    miller_area = wulff.miller_area_dict
    miller_energy = wulff.miller_energy_dict

    for miller in sorted_millers:
        area = miller_area.get(miller, 0.0)
        # Normal distance is proportional to surface energy (Wulff theorem)
        normal_dist = miller_energy.get(miller, 0.0)
        # Count facets for this orientation (including symmetry equivalents)
        n_facets = len([f for f in wulff.facets if tuple(f.normal) == miller or
                       any(np.allclose(f.normal, m) for m in [miller, tuple(-np.array(miller))])])

        miller_str = str(miller)
        lines.append(f'  {miller_str:<15} {normal_dist:<18.4f} {area:<15.4f} {n_facets}')

    lines.append('-' * 80)
    lines.append('')
    lines.append(sep)

    # Write to file
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def generate_standalone_script(
    selected_terminations: dict,
    formula: str,
    lattice_matrix: list,
    output_path: str,
) -> None:
    """
    Generate a standalone Python script for recreating/customizing the Wulff plot.

    Args:
        selected_terminations: Dict from select_lowest_energy_terminations()
        formula: Chemical formula
        lattice_matrix: 3x3 lattice matrix as nested list
        output_path: Path to save the Python script
    """
    # Prepare data for the script
    miller_list = list(selected_terminations.keys())
    e_surf_list = [selected_terminations[m]['gamma_J_m2'] for m in miller_list]

    script_content = f'''#!/usr/bin/env python
"""
Standalone Wulff Shape Visualization Script

This script recreates the Wulff shape plot for {formula} using the surface
energies computed by PS-TEROS. You can customize the visualization parameters
below without needing to rerun the AiiDA workflow.

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

Usage:
    python {output_path.split('/')[-1]}
"""

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.wulff import WulffShape

# ============================================================================
# DATA FROM PS-TEROS CALCULATION
# Do not modify unless you want to use different surface energies
# ============================================================================

# Chemical formula
FORMULA = "{formula}"

# Lattice matrix (Angstroms)
# Row vectors: a, b, c
LATTICE_MATRIX = {lattice_matrix}

# Miller indices and corresponding surface energies (J/m²)
# These are the lowest-energy terminations selected from the calculation
MILLER_INDICES = {miller_list}
SURFACE_ENERGIES = {e_surf_list}  # J/m²

# ============================================================================
# CUSTOMIZATION PARAMETERS
# Modify these to change the appearance of the plot
# ============================================================================

# Figure settings
FIGSIZE = (10, 8)      # Figure size in inches (width, height)
DPI = 150              # Resolution for saved figure

# Color scheme options: 'RdBu', 'YlGn', 'PuOr', 'BrBG', 'RdGy', 'PiYG'
# Or provide a list of colors matching the number of unique facets
COLOR_SET = 'RdBu'

# Transparency (0.0 = fully transparent, 1.0 = fully opaque)
ALPHA = 0.9

# Grid and axis display
SHOW_GRID = False
SHOW_AXIS = False

# Show area fractions on facets
SHOW_AREA = True

# Background color for facets not shown (due to 0 area)
OFF_COLOR = 'white'

# Viewing direction (controls the initial orientation)
# Try (1,1,1), (1,0,0), (1,1,0), etc.
VIEW_DIRECTION = (1, 1, 1)

# Title settings
TITLE_FONTSIZE = 14
TITLE_FONTWEIGHT = 'bold'

# Annotation settings
ANNOTATION_FONTSIZE = 11

# Output filename (change extension to .png for PNG output)
OUTPUT_FILENAME = "{formula.lower()}_wulff_shape.pdf"

# ============================================================================
# MAIN SCRIPT - No need to modify below this line
# ============================================================================

def main():
    """Generate the Wulff shape plot."""

    print(f"Generating Wulff shape for {{FORMULA}}...")

    # Create lattice
    lattice = Lattice(LATTICE_MATRIX)
    print(f"  Lattice parameter a = {{lattice.a:.4f}} Å")

    # Create Wulff shape
    wulff = WulffShape(lattice, MILLER_INDICES, SURFACE_ENERGIES)

    # Print summary
    print(f"\\n  Surface Energies:")
    for miller, gamma in zip(MILLER_INDICES, SURFACE_ENERGIES):
        area_frac = wulff.area_fraction_dict.get(miller, 0.0) * 100
        print(f"    {{str(miller):<12}} γ = {{gamma:.4f}} J/m²  ({{area_frac:.1f}}% of surface)")

    print(f"\\n  Wulff Shape Properties:")
    print(f"    Weighted Surface Energy: {{wulff.weighted_surface_energy:.4f}} J/m²")
    print(f"    Shape Anisotropy:        {{wulff.anisotropy:.4f}}")
    print(f"    Shape Factor:            {{wulff.shape_factor:.4f}} (1.0 = sphere)")

    # Generate plot
    fig = wulff.get_plot(
        color_set=COLOR_SET,
        grid_off=not SHOW_GRID,
        axis_off=not SHOW_AXIS,
        show_area=SHOW_AREA,
        alpha=ALPHA,
        off_color=OFF_COLOR,
        direction=VIEW_DIRECTION,
    )

    # Add title
    fig.suptitle(
        f'Wulff Shape: {{FORMULA}}',
        fontsize=TITLE_FONTSIZE,
        fontweight=TITLE_FONTWEIGHT,
        y=0.95,
    )

    # Add weighted surface energy annotation
    fig.text(
        0.5, 0.02,
        f'Weighted Surface Energy: {{wulff.weighted_surface_energy:.3f}} J/m²',
        ha='center',
        fontsize=ANNOTATION_FONTSIZE,
    )

    # Save figure
    output_format = OUTPUT_FILENAME.split('.')[-1]
    fig.savefig(OUTPUT_FILENAME, format=output_format, dpi=DPI, bbox_inches='tight')
    print(f"\\n  Plot saved to: {{OUTPUT_FILENAME}}")

    # Optionally show the plot interactively
    # Uncomment the following line to display the plot
    # plt.show()

    plt.close(fig)
    print("\\nDone!")


if __name__ == '__main__':
    main()
'''

    with open(output_path, 'w') as f:
        f.write(script_content)


@task.calcfunction
def generate_wulff_shape_data(
    surface_energies: orm.Dict,
    bulk_structure: orm.StructureData,
) -> dict:
    """
    Generate Wulff shape data, plot, and report from surface energies.

    This calcfunction:
    1. Selects the lowest-energy termination for each Miller index
    2. Constructs the Wulff shape using pymatgen
    3. Generates a PDF visualization
    4. Generates a text report with all numerical data
    5. Generates a standalone Python script for plot customization

    Args:
        surface_energies: Dict containing surface energies for all orientations
            with structure: {'hkl_111': {'term_0': {...}, 'term_1': {...}}, ...}
        bulk_structure: Relaxed bulk structure (for lattice parameters)

    Returns:
        Dictionary with outputs:
        - wulff_data: Dict with numerical Wulff shape properties
        - wulff_plot_pdf: SinglefileData with PDF visualization
        - wulff_report_txt: SinglefileData with text report
        - wulff_script_py: SinglefileData with standalone Python script
    """
    from pymatgen.core.lattice import Lattice
    from pymatgen.io.ase import AseAtomsAdaptor

    # Get surface energies dict
    surf_dict = surface_energies.get_dict()

    # Select lowest energy terminations
    selected = select_lowest_energy_terminations(surf_dict)

    if not selected:
        raise ValueError("No valid surface energies found in input data")

    # Get bulk structure info
    bulk_ase = bulk_structure.get_ase()

    # Convert to pymatgen structure for lattice
    adaptor = AseAtomsAdaptor()
    pmg_structure = adaptor.get_structure(bulk_ase)
    lattice = pmg_structure.lattice

    # Get chemical formula
    formula = pmg_structure.composition.reduced_formula

    # Prepare data for WulffShape
    miller_list = list(selected.keys())
    e_surf_list = [selected[m]['gamma_J_m2'] for m in miller_list]

    # Create Wulff shape
    wulff = create_wulff_shape(lattice, miller_list, e_surf_list)

    # Create temporary files for outputs
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        # Generate PDF plot
        pdf_path = os.path.join(tmpdir, f'{formula.lower()}_wulff_shape.pdf')
        generate_wulff_plot_pdf(wulff, formula, pdf_path)

        # Generate text report
        txt_path = os.path.join(tmpdir, f'{formula.lower()}_wulff_report.txt')
        generate_wulff_report_txt(selected, wulff, formula, lattice.a, txt_path)

        # Generate standalone script
        py_path = os.path.join(tmpdir, f'{formula.lower()}_wulff_plot.py')
        lattice_matrix = lattice.matrix.tolist()
        generate_standalone_script(selected, formula, lattice_matrix, py_path)

        # Create SinglefileData nodes
        wulff_plot_pdf = orm.SinglefileData(file=pdf_path)
        wulff_report_txt = orm.SinglefileData(file=txt_path)
        wulff_script_py = orm.SinglefileData(file=py_path)

    # Prepare numerical data dictionary
    # Convert miller tuples to strings for JSON serialization
    miller_energy_dict = {str(m): e for m, e in zip(miller_list, e_surf_list)}

    area_fraction_dict = {}
    for miller, frac in wulff.area_fraction_dict.items():
        area_fraction_dict[str(tuple(miller))] = float(frac)

    miller_area_dict = {}
    for miller, area in wulff.miller_area_dict.items():
        miller_area_dict[str(tuple(miller))] = float(area)

    # Selected terminations info
    termination_info = {}
    for miller, data in selected.items():
        termination_info[str(miller)] = {
            'termination': data['termination'],
            'gamma_J_m2': data['gamma_J_m2'],
        }

    wulff_data = orm.Dict(dict={
        'formula': formula,
        'lattice_parameter_a_A': float(lattice.a),
        'lattice_parameter_b_A': float(lattice.b),
        'lattice_parameter_c_A': float(lattice.c),
        'miller_energy_dict_J_m2': miller_energy_dict,
        'area_fraction_dict': area_fraction_dict,
        'miller_area_dict': miller_area_dict,
        'weighted_surface_energy_J_m2': float(wulff.weighted_surface_energy),
        'anisotropy': float(wulff.anisotropy),
        'shape_factor': float(wulff.shape_factor),
        'surface_area': float(wulff.surface_area),
        'volume': float(wulff.volume),
        'effective_radius': float(wulff.effective_radius),
        'selected_terminations': termination_info,
        'n_facets': len(miller_list),
    })

    return {
        'wulff_data': wulff_data,
        'wulff_plot_pdf': wulff_plot_pdf,
        'wulff_report_txt': wulff_report_txt,
        'wulff_script_py': wulff_script_py,
    }
