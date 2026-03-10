"""
Wulff Shape Construction for Metal and Intermetallic Surfaces.

This module provides functionality to:
1. Expand calculated Miller indices to symmetry-equivalent orientations
2. Construct Wulff shapes from surface energies using Pymatgen

The Wulff construction requires surface energies for multiple Miller indices.
By leveraging crystal symmetry, we can expand a small set of calculated
orientations (e.g., 3 unique families) to the full set needed for accurate
Wulff shape construction.

Example for FCC metals:
- Calculate: (111), (110), (001)
- Expand via symmetry to get 26+ equivalent orientations
- Construct Wulff shape showing equilibrium crystal morphology
"""

from __future__ import annotations

import typing as t
import numpy as np

from aiida import orm
from aiida_workgraph import task
from .wulff_geometry import extract_wulff_geometry


def get_symmetrically_equivalent_miller_indices(
    structure: orm.StructureData,
    miller_index: t.Union[list, tuple],
) -> list[tuple[int, int, int]]:
    """
    Get all symmetrically equivalent Miller indices for a given orientation.

    Uses the crystal's point group symmetry operations to generate all
    equivalent Miller indices. For cubic crystals, this can expand a single
    index to up to 48 equivalents (for general hkl).

    Args:
        structure: Crystal structure (for determining symmetry)
        miller_index: Miller index as list or tuple, e.g., [1, 1, 1]

    Returns:
        List of unique equivalent Miller indices as tuples

    Example:
        For FCC Au with (001):
        >>> equivalents = get_symmetrically_equivalent_miller_indices(au_structure, (0, 0, 1))
        >>> equivalents
        [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.ase import AseAtomsAdaptor

    # Convert AiiDA structure to pymatgen
    ase_atoms = structure.get_ase()
    pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)

    # Get symmetry operations
    sga = SpacegroupAnalyzer(pmg_structure)
    symm_ops = sga.get_symmetry_operations()

    # Apply all symmetry operations to the Miller index
    equivalent = set()
    miller = np.array(miller_index, dtype=float)

    for op in symm_ops:
        # Apply rotation matrix (transpose for reciprocal space)
        new_miller = np.dot(op.rotation_matrix.T, miller)
        # Round to integers
        new_miller_int = tuple(int(round(m)) for m in new_miller)
        equivalent.add(new_miller_int)

    return sorted(list(equivalent))


def _parse_miller_key(hkl_key: str) -> tuple[int, int, int]:
    """
    Parse a Miller index key string to a tuple.

    Handles formats like 'hkl_111', 'hkl_-111', 'hkl_11-1'.

    Note:
        This parser assumes single-digit Miller indices (0-9), which covers
        the vast majority of practical surface calculations. For high-index
        surfaces with indices >= 10, use a different key format.

    Args:
        hkl_key: Key string like 'hkl_111' or 'hkl_-111'

    Returns:
        Tuple of (h, k, l) integers
    """
    hkl_str = hkl_key.replace('hkl_', '')
    miller = []
    i = 0
    while i < len(hkl_str):
        if hkl_str[i] == '-':
            miller.append(-int(hkl_str[i + 1]))
            i += 2
        else:
            miller.append(int(hkl_str[i]))
            i += 1
    return tuple(miller)


def expand_surface_energies_with_symmetry(
    structure: orm.StructureData,
    surface_energies: dict[str, dict],
) -> dict[tuple[int, int, int], float]:
    """
    Expand calculated surface energies to all symmetry-equivalent orientations.

    Takes surface energies calculated for a few unique Miller indices and
    expands them to all equivalent orientations using crystal symmetry.
    For surface energy calculations, equivalent orientations have identical
    surface energies.

    Args:
        structure: Bulk crystal structure
        surface_energies: Dictionary from gather_surface_energies output
            Format: {'hkl_111': {'term_0': {'gamma_J_m2': 0.79, ...}}, ...}

    Returns:
        Dictionary mapping Miller index tuples to surface energies in J/m²
        Format: {(1, 1, 1): 0.79, (-1, -1, -1): 0.79, ...}

    Note:
        - For orientations with multiple terminations, uses the lowest energy
          termination (most stable surface).
        - Assumes single-digit Miller indices (0-9). High-index surfaces
          with indices >= 10 are not supported.
        - Negative surface energies are skipped with a warning (indicate
          unstable surfaces or calculation errors).
    """
    import warnings

    expanded = {}

    for hkl_key, terminations in surface_energies.items():
        # Parse hkl_key using helper function
        miller = _parse_miller_key(hkl_key)

        # Get the minimum surface energy across terminations (most stable)
        min_gamma = float('inf')
        for term_data in terminations.values():
            if isinstance(term_data, dict) and 'gamma_J_m2' in term_data:
                gamma = term_data['gamma_J_m2']
                if gamma < min_gamma:
                    min_gamma = gamma

        if min_gamma == float('inf'):
            continue

        # Validate: negative surface energies indicate unstable surfaces
        if min_gamma < 0:
            warnings.warn(
                f"Negative surface energy for {miller}: {min_gamma:.4f} J/m². "
                "This indicates an unstable surface or calculation error. Skipping."
            )
            continue

        # Expand to all equivalent orientations
        equivalents = get_symmetrically_equivalent_miller_indices(structure, miller)
        for eq_miller in equivalents:
            # Only add if not already present (avoid overwriting)
            if eq_miller not in expanded:
                expanded[eq_miller] = min_gamma

    return expanded


@task.calcfunction
def build_wulff_shape(
    bulk_structure: orm.StructureData,
    surface_energies: orm.Dict,
    symprec: orm.Float = None,
    only_stoichiometric_symmetric: orm.Bool = None,
) -> orm.Dict:
    """
    Build Wulff shape from calculated surface energies.

    Uses Pymatgen's WulffShape class to construct the equilibrium crystal
    shape from surface energies. Automatically expands calculated orientations
    to symmetry equivalents.

    Args:
        bulk_structure: Bulk crystal structure (for lattice and symmetry)
        surface_energies: Output from gather_surface_energies
            Format: {'hkl_111': {'term_0': {'gamma_J_m2': 0.79, ...}}, ...}
        symprec: Symmetry precision for WulffShape (default: 1e-5)
        only_stoichiometric_symmetric: If True, only include terminations that are
            marked as stoichiometric (is_stoichiometric=True in the surface energy data).
            This filters out non-stoichiometric surfaces from the Wulff shape.

    Returns:
        Dictionary containing:
        - miller_energy_dict: Dict mapping Miller indices to surface energies
        - expanded_miller_energy_dict: Same with symmetry equivalents
        - n_calculated_orientations: Number of unique orientations calculated
        - n_expanded_orientations: Number after symmetry expansion
        - shape_factor: Wulff shape factor (1 = sphere, >1 = anisotropic)
        - anisotropy: Degree of anisotropy
        - weighted_surface_energy: Area-weighted average surface energy (J/m²)
        - total_surface_area: Total surface area of Wulff shape (arbitrary units)
        - volume: Volume of Wulff shape (arbitrary units)
        - facet_fractions: Dict of Miller index to fraction of total surface area
        - dominant_facet: Miller index with largest surface area fraction
        - wulff_shape_valid: Whether Wulff shape was successfully constructed
        - only_stoichiometric: Whether stoichiometric filter was applied
        - n_filtered_out: Number of terminations filtered out (if filter applied)
        - wulff_geometry: Serializable polyhedron (vertices, faces, normals, bounding radius)

    Example:
        For FCC Au with (111), (110), (001):
        - dominant_facet: (1, 1, 1) with ~60% surface area
        - shape_factor: ~1.02 (nearly spherical due to low anisotropy)
    """
    from pymatgen.analysis.wulff import WulffShape
    from pymatgen.io.ase import AseAtomsAdaptor

    # Set default symprec
    if symprec is None:
        symprec_val = 1e-5
    else:
        symprec_val = symprec.value

    # Parse stoichiometric filter flag
    filter_stoichiometric = False
    if only_stoichiometric_symmetric is not None:
        filter_stoichiometric = only_stoichiometric_symmetric.value

    # Get surface energies dict
    surf_energies_dict = surface_energies.get_dict()

    # Apply stoichiometric filter if requested
    n_filtered_out = 0
    if filter_stoichiometric:
        filtered_dict = {}
        for hkl_key, terminations in surf_energies_dict.items():
            filtered_terms = {}
            for term_key, term_data in terminations.items():
                if isinstance(term_data, dict):
                    is_stoich = term_data.get('is_stoichiometric', True)
                    if is_stoich:
                        filtered_terms[term_key] = term_data
                    else:
                        n_filtered_out += 1
                else:
                    # Can't check, keep it
                    filtered_terms[term_key] = term_data

            if filtered_terms:
                filtered_dict[hkl_key] = filtered_terms

        surf_energies_dict = filtered_dict

    # Expand to symmetry equivalents
    expanded_energies = expand_surface_energies_with_symmetry(
        bulk_structure, surf_energies_dict
    )

    if len(expanded_energies) == 0:
        return orm.Dict(dict={
            'wulff_shape_valid': False,
            'error': 'No valid surface energies found',
            'n_calculated_orientations': 0,
            'n_expanded_orientations': 0,
        })

    # Convert to pymatgen structure and get lattice
    ase_atoms = bulk_structure.get_ase()
    pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    lattice = pmg_structure.lattice

    # Prepare inputs for WulffShape
    miller_list = list(expanded_energies.keys())
    e_surf_list = [expanded_energies[m] for m in miller_list]

    # Build original (non-expanded) miller_energy_dict for reference
    original_miller_energy = {}
    for hkl_key, terminations in surf_energies_dict.items():
        miller_tuple = _parse_miller_key(hkl_key)

        min_gamma = float('inf')
        for term_data in terminations.values():
            if isinstance(term_data, dict) and 'gamma_J_m2' in term_data:
                gamma = term_data['gamma_J_m2']
                if gamma < min_gamma:
                    min_gamma = gamma

        if min_gamma != float('inf') and min_gamma > 0:
            original_miller_energy[str(miller_tuple)] = min_gamma

    try:
        # Construct Wulff shape
        wulff = WulffShape(lattice, miller_list, e_surf_list, symprec=symprec_val)

        # Extract facet fractions (area of each facet / total area)
        facet_fractions = {}
        total_area = sum(wulff.color_area)
        for i, miller in enumerate(wulff.miller_list):
            miller_str = str(tuple(miller))
            if wulff.on_wulff[i] and total_area > 0:
                facet_fractions[miller_str] = wulff.color_area[i] / total_area

        # Find dominant facet
        dominant_facet = None
        max_fraction = 0
        for miller_str, fraction in facet_fractions.items():
            if fraction > max_fraction:
                max_fraction = fraction
                dominant_facet = miller_str

        # Build result dictionary
        geometry = extract_wulff_geometry(wulff)

        result = {
            'wulff_shape_valid': True,
            'miller_energy_dict': original_miller_energy,
            'expanded_miller_energy_dict': {str(k): v for k, v in expanded_energies.items()},
            'n_calculated_orientations': len(original_miller_energy),
            'n_expanded_orientations': len(expanded_energies),
            'shape_factor': float(wulff.shape_factor),
            'anisotropy': float(wulff.anisotropy),
            'weighted_surface_energy': float(wulff.weighted_surface_energy),
            'total_surface_area': float(wulff.surface_area),
            'volume': float(wulff.volume),
            'facet_fractions': facet_fractions,
            'dominant_facet': dominant_facet,
            'dominant_facet_fraction': max_fraction if dominant_facet else 0.0,
            'only_stoichiometric': filter_stoichiometric,
            'n_filtered_out': n_filtered_out,
            'wulff_geometry': geometry,
        }

        return orm.Dict(dict=result)

    except Exception as e:
        return orm.Dict(dict={
            'wulff_shape_valid': False,
            'error': str(e),
            'miller_energy_dict': original_miller_energy,
            'expanded_miller_energy_dict': {str(k): v for k, v in expanded_energies.items()},
            'n_calculated_orientations': len(original_miller_energy),
            'n_expanded_orientations': len(expanded_energies),
        })


def visualize_wulff_shape(
    bulk_structure: orm.StructureData,
    surface_energies: t.Union[orm.Dict, dict],
    color_set: str = 'PuBu',
    direction: tuple = (1, 1, 1),
    show_area: bool = False,
    alpha: float = 1.0,
    save_path: str = None,
):
    """
    Visualize the Wulff shape from surface energies.

    This is a helper function for interactive visualization, not a calcfunction.
    It creates a 3D plot of the Wulff shape using matplotlib.

    Args:
        bulk_structure: Bulk crystal structure
        surface_energies: Output from gather_surface_energies (Dict or dict)
        color_set: Matplotlib colormap name (default: 'PuBu')
        direction: Viewing direction as (x, y, z) tuple
        show_area: Whether to show area labels on facets
        alpha: Transparency (0-1)
        save_path: If provided, save figure to this path

    Returns:
        matplotlib Axes3D object
    """
    from pymatgen.analysis.wulff import WulffShape
    from pymatgen.io.ase import AseAtomsAdaptor

    # Handle Dict node or plain dict
    if isinstance(surface_energies, orm.Dict):
        surf_energies_dict = surface_energies.get_dict()
    else:
        surf_energies_dict = surface_energies

    # Expand to symmetry equivalents
    expanded_energies = expand_surface_energies_with_symmetry(
        bulk_structure, surf_energies_dict
    )

    if len(expanded_energies) == 0:
        raise ValueError("No valid surface energies found")

    # Convert to pymatgen structure and get lattice
    ase_atoms = bulk_structure.get_ase()
    pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    lattice = pmg_structure.lattice

    # Prepare inputs
    miller_list = list(expanded_energies.keys())
    e_surf_list = [expanded_energies[m] for m in miller_list]

    # Create Wulff shape
    wulff = WulffShape(lattice, miller_list, e_surf_list)

    # Get plot
    ax = wulff.get_plot(
        color_set=color_set,
        direction=direction,
        show_area=show_area,
        alpha=alpha,
    )

    if save_path:
        import matplotlib.pyplot as plt
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    return ax


def get_wulff_shape_summary(wulff_result: t.Union[orm.Dict, dict]) -> str:
    """
    Generate a human-readable summary of Wulff shape results.

    Args:
        wulff_result: Output from build_wulff_shape

    Returns:
        Formatted string summary
    """
    if isinstance(wulff_result, orm.Dict):
        data = wulff_result.get_dict()
    else:
        data = wulff_result

    lines = []
    lines.append("=" * 60)
    lines.append("WULFF SHAPE ANALYSIS")
    lines.append("=" * 60)

    if not data.get('wulff_shape_valid', False):
        lines.append(f"Error: {data.get('error', 'Unknown error')}")
        return '\n'.join(lines)

    lines.append(f"\nOrientations:")
    lines.append(f"  Calculated: {data['n_calculated_orientations']}")
    lines.append(f"  After symmetry expansion: {data['n_expanded_orientations']}")

    lines.append(f"\nShape Properties:")
    lines.append(f"  Shape factor: {data['shape_factor']:.4f} (1.0 = sphere)")
    lines.append(f"  Anisotropy: {data['anisotropy']:.4f}")
    lines.append(f"  Weighted surface energy: {data['weighted_surface_energy']:.4f} J/m²")

    lines.append(f"\nDominant Facet:")
    lines.append(f"  {data['dominant_facet']} ({data['dominant_facet_fraction']*100:.1f}% of surface)")

    lines.append(f"\nFacet Distribution:")
    facet_fractions = data.get('facet_fractions', {})
    # Sort by fraction descending
    sorted_facets = sorted(facet_fractions.items(), key=lambda x: x[1], reverse=True)
    for miller_str, fraction in sorted_facets[:10]:  # Top 10
        if fraction > 0.001:  # Only show significant facets
            lines.append(f"  {miller_str}: {fraction*100:.1f}%")

    lines.append(f"\nSurface Energies (J/m²):")
    miller_energy = data.get('miller_energy_dict', {})
    for miller_str, energy in sorted(miller_energy.items()):
        lines.append(f"  {miller_str}: {energy:.4f}")

    lines.append("=" * 60)

    return '\n'.join(lines)
