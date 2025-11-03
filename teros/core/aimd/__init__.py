"""AIMD standalone module for PS-TEROS.

This module provides standalone AIMD functionality with full control over:
- Sequential multi-stage AIMD on multiple structures
- Optional supercell transformations
- Restart chaining between stages
- Concurrency control with max_concurrent_jobs

NOTE: This package (aimd/) coexists with teros.core.aimd module (aimd.py).
      aimd_single_stage_scatter must be imported directly from parent where needed.
"""

# Import from the standalone module
from .workgraph import build_aimd_workgraph
from .utils import organize_aimd_results

__all__ = [
    'build_aimd_workgraph',
    'organize_aimd_results',
]
