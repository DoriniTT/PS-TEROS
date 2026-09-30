"""Legacy implementation modules retained behind the psteros 1.0 API.

New projects import the compact public API from :mod:`psteros`.  The modules in
this namespace remain importable for audited VASP calculations, but the former
generic composition API has been removed.
"""

from .exceptions import (
    ConfigurationError,
    ConvergenceError,
    StructureError,
    ValidationError,
    WorkflowError,
)

__all__ = [
    "ConfigurationError",
    "ConvergenceError",
    "StructureError",
    "ValidationError",
    "WorkflowError",
]
