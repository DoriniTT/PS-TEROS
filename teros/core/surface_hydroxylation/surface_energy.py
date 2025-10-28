"""
Surface energy calculator for hydroxylated surfaces.

Implements equations 4-10 from Section S3 for calculating surface Gibbs free energies.
"""

from collections import Counter
from aiida.engine import calcfunction
from aiida.orm import Float, Int


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

    # Validate n_h is even (required for H_{2x} formula)
    if n_h % 2 != 0:
        raise ValueError(f"Number of H atoms must be even for H_{{2x}} formula, got n_h={n_h}")

    # Calculate oxygen deficit
    expected_o = n * bulk_formula['O']
    actual_o = slab_formula['O']
    n_o_deficit = expected_o - actual_o

    # Derive x and y
    # From modified formula: Ag_{3n}P_nO_{4n+x-2y}H_{2x}
    # n_h = 2x → x = n_h / 2
    # n_o_deficit = -(x - 2y) → y = (n_o_deficit + x) / 2
    # TODO: Verify H_{2x} formula interpretation with experimental data
    x = n_h // 2
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


@calcfunction
def calc_delta_g_reaction1(E_slab, E_bulk, n, x, y, delta_mu_h2o, delta_mu_o2):
    """
    Calculate formation energy using Reaction 1: H2O/O2 reservoirs.

    Reaction: n Ag3PO4 + x H2O ↔ surface + y O2
    Use case: Water-rich oxidizing environment

    Equation 4: ΔG = E + 2y·μ_O - n·E_bulk - x·μ_H2O

    Args:
        E_slab: Float - Total energy of modified slab (eV)
        E_bulk: Float - Bulk energy per formula unit (eV/f.u.)
        n: Int - Number of bulk formula units in slab
        x: Int - Number of OH groups added
        y: Int - Net oxygen atoms removed (can be negative)
        delta_mu_h2o: Float - Δμ(H2O) from JANAF (eV/molecule)
        delta_mu_o2: Float - Δμ(O2) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV

    Raises:
        ValueError: If validation fails (n <= 0 or E_bulk >= 0)
    """
    # Validate inputs
    if n.value <= 0:
        raise ValueError(f"n must be positive, got {n.value}")

    if E_bulk.value >= 0:
        raise ValueError(
            f"E_bulk should be negative (binding energy), got {E_bulk.value} eV. "
            f"Check that bulk energy is per formula unit, not total energy."
        )

    # Calculate atomic oxygen chemical potential from O2
    # μ(O) = (1/2)·μ(O2) = (1/2)·Δμ(O2)
    mu_o = 0.5 * delta_mu_o2.value

    # Equation 4: ΔG = E + 2y·μ_O - n·E_bulk - x·μ_H2O
    delta_g = (
        E_slab.value +
        2 * y.value * mu_o -
        n.value * E_bulk.value -
        x.value * delta_mu_h2o.value
    )

    return Float(delta_g)
