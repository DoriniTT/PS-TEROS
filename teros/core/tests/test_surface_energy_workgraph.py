"""Tests for WorkGraph integration of surface energy calculations."""
import pytest
from unittest.mock import Mock
from aiida import orm


def test_create_surface_energy_task():
    """Test that surface energy task is created with correct structure."""
    from teros.core.surface_hydroxylation.surface_energy_workgraph import (
        create_surface_energy_task
    )

    # Mock upstream tasks
    aggregate_task = Mock()
    aggregate_task.outputs = Mock()
    aggregate_task.outputs.structures = Mock()
    aggregate_task.outputs.energies = Mock()
    bulk_task = Mock()
    bulk_task.outputs = Mock()
    bulk_task.outputs.structure = Mock()
    bulk_task.outputs.energy = Mock()
    pristine_task = Mock()
    pristine_task.outputs = Mock()

    # Create task (returns TaskSocketNamespace representing outputs)
    task_outputs = create_surface_energy_task(
        aggregate_results_task=aggregate_task,
        bulk_relax_task=bulk_task,
        pristine_relax_task=pristine_task,
    )

    # Verify task outputs structure (WorkGraph tasks return outputs namespace directly)
    assert task_outputs is not None
    assert hasattr(task_outputs, 'reaction1_results')
    assert hasattr(task_outputs, 'reaction2_results')
    assert hasattr(task_outputs, 'reaction3_results')
