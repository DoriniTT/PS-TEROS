"""Utility functions for AIMD module."""


def validate_stage_sequence(stages: list[dict]) -> None:
    """
    Validate AIMD stage sequence format.

    Args:
        stages: List of stage dicts, each must contain 'temperature' and 'steps'

    Raises:
        ValueError: If stages empty or missing required keys
    """
    if not stages:
        raise ValueError("aimd_stages must contain at least one stage")

    for idx, stage in enumerate(stages):
        if 'temperature' not in stage:
            raise ValueError(
                f"Stage {idx} missing required key 'temperature'. "
                f"Each stage must contain {{'temperature': K, 'steps': N}}"
            )
        if 'steps' not in stage:
            raise ValueError(
                f"Stage {idx} missing required key 'steps'. "
                f"Each stage must contain {{'temperature': K, 'steps': N}}"
            )


def validate_supercell_spec(spec: list[int]) -> None:
    """
    Validate supercell specification.

    Args:
        spec: [nx, ny, nz] supercell dimensions

    Raises:
        ValueError: If spec not valid 3D integer list with positive values
    """
    if not isinstance(spec, list):
        raise ValueError(f"Supercell spec must be a list, got {type(spec).__name__}")

    if len(spec) != 3:
        raise ValueError(
            f"Supercell spec must be a 3-element list [nx, ny, nz], got {len(spec)} elements"
        )

    for idx, val in enumerate(spec):
        if not isinstance(val, int):
            raise ValueError(
                f"Supercell spec elements must be positive integers, "
                f"element {idx} is {type(val).__name__}"
            )
        if val <= 0:
            raise ValueError(
                f"Supercell spec elements must be positive integers, "
                f"element {idx} is {val}"
            )
