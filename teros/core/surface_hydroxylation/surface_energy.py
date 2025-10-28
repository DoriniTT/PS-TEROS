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
