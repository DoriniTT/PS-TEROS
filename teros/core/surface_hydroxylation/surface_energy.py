"""
Surface energy calculator for hydroxylated surfaces.

Implements equations 4-10 from Section S3 for calculating surface Gibbs free energies.
"""

from collections import Counter


def _get_formula_dict(structure):
    """
    Extract chemical formula as dictionary from structure.

    Args:
        structure: AiiDA StructureData

    Returns:
        dict: Element symbols to counts, e.g. {'Ag': 12, 'P': 4, 'O': 16}
    """
    # Get chemical symbols
    symbols = [site.kind_name for site in structure.sites]

    # Count occurrences
    formula_dict = dict(Counter(symbols))

    return formula_dict


def analyze_composition(slab_structure, bulk_structure):
    """
    Determine stoichiometry parameters (n, x, y) from slab composition.

    Algorithm:
        1. Parse bulk formula (e.g., Ag3PO4 → {'Ag': 3, 'P': 1, 'O': 4})
        2. Parse slab formula (e.g., Ag12P4O18H2)
        3. Determine n from element ratios
        4. Calculate x from hydrogen count: x = n_h / 2
        5. Calculate y from oxygen deficit: y = (n_o_deficit + x) / 2

    Args:
        slab_structure: AiiDA StructureData for modified surface
        bulk_structure: AiiDA StructureData for bulk unit cell

    Returns:
        dict: {
            'n': int,              # Bulk formula units in slab
            'x': int,              # OH groups added
            'y': int,              # Net O atoms removed
            'n_h': int,            # Total H atoms
            'n_o_deficit': int,    # O deficit vs n×bulk
            'formulas': {
                'bulk': str,
                'slab': str
            }
        }

    Raises:
        ValueError: If composition analysis fails or gives non-integer results
    """
    # Parse formulas
    bulk_formula = _get_formula_dict(bulk_structure)
    slab_formula = _get_formula_dict(slab_structure)

    # Determine n (number of bulk units in slab)
    # Use elements present in bulk (exclude H which isn't in bulk)
    ratios = []
    for element, bulk_count in bulk_formula.items():
        if element in slab_formula:
            ratio = slab_formula[element] / bulk_count
            ratios.append((element, ratio))

    # All ratios should be equal for stoichiometric slab
    # Use most common ratio as n
    n_values = [r[1] for r in ratios]
    n = round(sum(n_values) / len(n_values))  # Average and round

    # Validate n is positive integer
    if n <= 0:
        raise ValueError(f"Invalid n={n}, must be positive")

    # Calculate hydrogen addition
    n_h = slab_formula.get('H', 0)

    # Calculate oxygen deficit
    expected_o = n * bulk_formula['O']
    actual_o = slab_formula['O']
    n_o_deficit = expected_o - actual_o

    # Derive x and y
    # From modified formula: Ag_{3n}P_nO_{4n+x-2y}H_x
    # n_h = x → x = n_h
    # n_o_deficit = -(x - 2y) → y = (n_o_deficit + x) / 2
    x = n_h
    y = (n_o_deficit + x) // 2

    # Validate results
    if (n_o_deficit + x) % 2 != 0:
        raise ValueError(f"O balance gives non-integer y: n_o_deficit={n_o_deficit}, x={x}")

    if (n_o_deficit + x) != 2 * y:
        raise ValueError(f"O balance inconsistent: n_o_deficit={n_o_deficit}, x={x}, y={y}")

    # Format output
    bulk_str = ''.join(f"{elem}{count}" for elem, count in sorted(bulk_formula.items()))
    slab_str = ''.join(f"{elem}{count}" for elem, count in sorted(slab_formula.items()))

    return {
        'n': n,
        'x': x,
        'y': y,
        'n_h': n_h,
        'n_o_deficit': n_o_deficit,
        'formulas': {
            'bulk': bulk_str,
            'slab': slab_str,
        }
    }
