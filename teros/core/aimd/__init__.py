"""AIMD standalone module for PS-TEROS."""

# Import from the standalone module
from .workgraph import build_aimd_workgraph
from .utils import organize_aimd_results

# Import from the original aimd.py file (sibling module) and re-export
# This resolves the naming conflict between aimd.py file and aimd/ directory
import importlib.util
import os
_aimd_file_path = os.path.join(os.path.dirname(__file__), '..', 'aimd.py')
_spec = importlib.util.spec_from_file_location("_teros_core_aimd_file", _aimd_file_path)
_aimd_file = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_aimd_file)
aimd_single_stage_scatter = _aimd_file.aimd_single_stage_scatter

__all__ = [
    'build_aimd_workgraph',
    'organize_aimd_results',
    'aimd_single_stage_scatter',
]
