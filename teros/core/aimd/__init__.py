"""AIMD standalone module for PS-TEROS.

This module provides standalone AIMD functionality with full control over:
- Sequential multi-stage AIMD on multiple structures
- Optional supercell transformations
- Restart chaining between stages
- Concurrency control with max_concurrent_jobs

NOTE: This package (aimd/) coexists with the teros.core.aimd module (aimd.py).
Due to Python's package precedence rules, we must explicitly load functions
from the aimd.py file to avoid circular imports.
"""

# Import from the standalone module
from .workgraph import build_aimd_workgraph
from .utils import organize_aimd_results

# Import aimd_single_stage_scatter from parent aimd.py module
# We need this because workgraph.py uses it to build the AIMD stages
# Due to package/module naming conflict, we use explicit file loading
from pathlib import Path
import importlib.util

def _load_aimd_module():
    """Load the parent aimd.py module to avoid naming conflict."""
    aimd_file = Path(__file__).parent.parent / 'aimd.py'
    spec = importlib.util.spec_from_file_location('teros.core._aimd_file', aimd_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

_aimd_file_module = _load_aimd_module()
aimd_single_stage_scatter = _aimd_file_module.aimd_single_stage_scatter

__all__ = [
    'build_aimd_workgraph',
    'organize_aimd_results',
    'aimd_single_stage_scatter',
]
