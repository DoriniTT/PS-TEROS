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
