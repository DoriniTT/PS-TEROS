"""
Core module for the Teros package.

Contains the main workgraph implementation for surface energy calculations.
"""

try:
    from .workgraph import create_teros_workgraph
except ImportError:
    # Fallback for development/testing
    pass

__all__ = ['create_teros_workgraph']