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
        delta_n_values: List of delta_n values for fractional charges
        fukui_type: Type of Fukui calculation ('plus' or 'minus')

    Raises:
        ValueError: If inputs are invalid
    """
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


# Default delta_n values for Fukui calculations (from FukuiGrid documentation)
DEFAULT_DELTA_N_VALUES = [0.0, 0.05, 0.10, 0.15]
