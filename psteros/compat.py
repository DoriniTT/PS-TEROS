"""Narrow compatibility access to the pre-1.0 VASP workflow builder.

New code should use :func:`psteros.build_surface_workgraph`.  This module is
kept separate so the historic wide VASP builder is never the default API.
"""

from psteros.core.workgraph import build_core_workgraph

__all__ = ["build_core_workgraph"]

