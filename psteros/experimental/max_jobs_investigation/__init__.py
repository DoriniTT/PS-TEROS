"""max_number_jobs investigation module."""

from .workgraph_functions import (
    process_structures_approach_1,
    process_structures_approach_3,
)
from .mock_tasks import mock_vasp_calculation

__all__ = [
    'process_structures_approach_1',
    'process_structures_approach_3',
    'mock_vasp_calculation',
]
