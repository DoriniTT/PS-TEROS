"""Utility functions for Fukui calculations.

Pure Python utilities for label sanitization and input validation.
"""

from typing import List


def make_delta_label(delta_n: float) -> str:
    """
    Create a sanitized label for a delta_n value.

    Converts float values to valid Python identifiers for AiiDA link labels.

    Args:
        delta_n: Float value (e.g., 0.05)

    Returns:
        Sanitized label (e.g., 'delta_0_05')

    Examples:
        >>> make_delta_label(0.0)
        'delta_0_00'
        >>> make_delta_label(0.05)
        'delta_0_05'
        >>> make_delta_label(0.15)
        'delta_0_15'
    """
    # Format as string with 2 decimal places
    formatted = f"{delta_n:.2f}"
    # Replace dot with underscore for AiiDA link label compatibility
    sanitized = formatted.replace('.', '_')
    return f"delta_{sanitized}"


def validate_fukui_inputs(
    nelect_neutral: int,
    delta_n_values: List[float],
    fukui_type: str,
) -> None:
    """
    Validate Fukui workflow inputs.

    Args:
        nelect_neutral: Number of electrons in neutral system
        delta_n_values: List of delta_n values for fractional charges.
                       **Maximum 4 values supported** due to internal WorkGraph
                       limitations.
        fukui_type: Type of Fukui calculation ('plus' or 'minus')

    Raises:
        ValueError: If inputs are invalid

    Warns:
        UserWarning: If delta_n=0.0 is not included (required for proper interpolation)
        UserWarning: If very small delta_n values may cause SCF convergence issues
        UserWarning: If large delta_n values may yield unphysical results
    """
    import warnings

    # Check nelect_neutral is positive
    if nelect_neutral <= 0:
        raise ValueError(f"nelect_neutral must be positive, got {nelect_neutral}")

    # Check delta_n_values are non-negative
    for delta_n in delta_n_values:
        if delta_n < 0:
            raise ValueError(f"delta_n values must be non-negative, got {delta_n}")

    # Check fukui_type
    if fukui_type not in ('plus', 'minus'):
        raise ValueError(f"fukui_type must be 'plus' or 'minus', got '{fukui_type}'")

    # Check that delta_n values won't result in invalid NELECT
    max_delta = max(delta_n_values)
    if fukui_type == 'plus':
        # Removing electrons: NELECT = nelect_neutral - delta_n
        min_nelect = nelect_neutral - max_delta
    else:
        # Adding electrons: NELECT = nelect_neutral + delta_n
        min_nelect = nelect_neutral  # Adding electrons won't reduce NELECT

    if min_nelect < 1:
        raise ValueError(
            f"delta_n values too large: would result in NELECT={min_nelect:.2f} < 1. "
            f"Reduce delta_n values or increase nelect_neutral."
        )

    # Check maximum number of delta_n values (internal WorkGraph limitation)
    if len(delta_n_values) > 4:
        raise ValueError(
            f"Maximum of 4 delta_n values supported, got {len(delta_n_values)}. "
            f"This limitation is due to the internal @task.graph implementation. "
            f"For higher-order polynomial fitting, consider running multiple workflows."
        )

    # Warning: delta_n=0.0 should be included for proper interpolation
    if 0.0 not in delta_n_values:
        warnings.warn(
            "delta_n=0.0 (neutral reference) is not in delta_n_values. "
            "Including the neutral state is recommended for proper interpolation.",
            UserWarning
        )

    # Warning: very small delta_n values may cause SCF convergence issues
    small_dn = [dn for dn in delta_n_values if 0 < dn < 0.01]
    if small_dn:
        warnings.warn(
            f"Very small delta_n values ({small_dn}) may cause SCF convergence issues "
            f"with fractional electron occupation. Consider using delta_n >= 0.01.",
            UserWarning
        )

    # Warning: large delta_n values may yield unphysical results
    large_dn = [dn for dn in delta_n_values if dn > 0.5]
    if large_dn:
        warnings.warn(
            f"Large delta_n values ({large_dn}) may yield unphysical results. "
            f"The linear response regime is typically valid for delta_n <= 0.2.",
            UserWarning
        )


# Default delta_n values for Fukui calculations (from FukuiGrid documentation)
DEFAULT_DELTA_N_VALUES = [0.0, 0.05, 0.10, 0.15]


def calculate_nelect(
    structure,
    potential_family: str,
    potential_mapping: dict,
) -> int:
    """
    Calculate NELECT (number of valence electrons) from structure and POTCARs.

    Automatically looks up ZVAL from POTCAR files stored in AiiDA and
    calculates the total number of valence electrons.

    Args:
        structure: AiiDA StructureData or pymatgen Structure
        potential_family: POTCAR family name (e.g., 'PBE')
        potential_mapping: Element to POTCAR symbol mapping
                          (e.g., {'Sn': 'Sn_d', 'O': 'O'})

    Returns:
        Total number of valence electrons (NELECT)

    Example:
        >>> from aiida import orm
        >>> structure = orm.load_node(12345)
        >>> nelect = calculate_nelect(
        ...     structure,
        ...     potential_family='PBE',
        ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'}
        ... )
        >>> print(f"NELECT = {nelect}")
        NELECT = 312
    """
    import re
    from aiida import orm
    from aiida_vasp.data.potcar import PotcarData

    # Handle both AiiDA StructureData and pymatgen Structure
    if hasattr(structure, 'get_pymatgen_structure'):
        pmg_structure = structure.get_pymatgen_structure()
    else:
        pmg_structure = structure

    # Get element counts from structure
    composition = pmg_structure.composition.as_dict()

    # Get elements list
    elements = list(composition.keys())

    # Get POTCAR nodes
    potcars_dict = PotcarData.get_potcars_dict(
        elements=elements,
        family_name=potential_family,
        mapping=potential_mapping,
    )

    # Extract ZVAL from each POTCAR and calculate total NELECT
    total_nelect = 0.0
    zval_info = {}

    for element, potcar in potcars_dict.items():
        # Get POTCAR content
        content = potcar.get_content()
        if isinstance(content, bytes):
            content = content.decode('utf-8')

        # Parse ZVAL from POTCAR content
        # Format: "POMASS =  118.710; ZVAL   =   14.000    mass and valenz"
        zval_match = re.search(r'ZVAL\s*=\s*([\d.]+)', content)
        if zval_match:
            zval = float(zval_match.group(1))
        else:
            raise ValueError(f"Could not find ZVAL in POTCAR for {element}")

        count = composition[element]
        contribution = zval * count
        total_nelect += contribution
        zval_info[element] = {'zval': zval, 'count': count, 'contribution': contribution}

    return int(round(total_nelect))


def print_nelect_breakdown(
    structure,
    potential_family: str,
    potential_mapping: dict,
) -> None:
    """
    Print detailed NELECT breakdown showing contribution from each element.

    Args:
        structure: AiiDA StructureData or pymatgen Structure
        potential_family: POTCAR family name
        potential_mapping: Element to POTCAR symbol mapping
    """
    import re
    from aiida_vasp.data.potcar import PotcarData

    # Handle both AiiDA StructureData and pymatgen Structure
    if hasattr(structure, 'get_pymatgen_structure'):
        pmg_structure = structure.get_pymatgen_structure()
    else:
        pmg_structure = structure

    composition = pmg_structure.composition.as_dict()
    elements = list(composition.keys())

    potcars_dict = PotcarData.get_potcars_dict(
        elements=elements,
        family_name=potential_family,
        mapping=potential_mapping,
    )

    print("\nNELECT Breakdown:")
    print("-" * 60)
    print(f"{'Element':<10} {'POTCAR':<15} {'ZVAL':<10} {'Count':<10} {'Subtotal':<10}")
    print("-" * 60)

    total_nelect = 0.0
    for element in sorted(composition.keys()):
        potcar = potcars_dict[element]
        content = potcar.get_content()
        if isinstance(content, bytes):
            content = content.decode('utf-8')

        zval_match = re.search(r'ZVAL\s*=\s*([\d.]+)', content)
        zval = float(zval_match.group(1)) if zval_match else 0

        count = composition[element]
        subtotal = zval * count
        total_nelect += subtotal

        print(f"{element:<10} {potential_mapping.get(element, element):<15} {zval:<10.1f} {count:<10.0f} {subtotal:<10.0f}")

    print("-" * 60)
    print(f"{'TOTAL NELECT:':<47} {int(round(total_nelect))}")
    print("-" * 60)
