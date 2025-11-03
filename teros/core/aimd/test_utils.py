"""Tests for AIMD utils."""
import pytest
from teros.core.aimd.utils import validate_stage_sequence


def test_validate_stage_sequence_valid():
    """Valid stage sequence passes."""
    stages = [
        {'temperature': 300, 'steps': 100},
        {'temperature': 500, 'steps': 200},
    ]
    validate_stage_sequence(stages)  # Should not raise


def test_validate_stage_sequence_missing_temperature():
    """Stage missing temperature raises ValueError."""
    stages = [{'steps': 100}]
    with pytest.raises(ValueError, match="missing required key 'temperature'"):
        validate_stage_sequence(stages)


def test_validate_stage_sequence_missing_steps():
    """Stage missing steps raises ValueError."""
    stages = [{'temperature': 300}]
    with pytest.raises(ValueError, match="missing required key 'steps'"):
        validate_stage_sequence(stages)


def test_validate_stage_sequence_empty():
    """Empty stage list raises ValueError."""
    with pytest.raises(ValueError, match="at least one stage"):
        validate_stage_sequence([])
