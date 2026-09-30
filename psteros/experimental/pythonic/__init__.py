"""
Pythonic Scatter-Gather Pattern for AiiDA-WorkGraph
====================================================

This module demonstrates the scatter-gather pattern for parallel workflow execution
using AiiDA-WorkGraph's pythonic approach (@task.graph decorator).

Applied Use Case: Slab Surface Generation and Analysis
-------------------------------------------------------
- Generate multiple slab terminations from bulk structure
- Relax all slabs in parallel with VASP
- Compute ab initio atomistic thermodynamics (AIAT) for each surface

Key Features:
- Dynamic parallelization (unknown number of tasks at definition time)
- Full AiiDA provenance tracking
- Plain Python control flow (no DSL required)
- Conditional workflow steps
- Nested graph tasks

Status: ✅ Production Ready
---------------------------
Successfully tested with real VASP calculations (4 Ag₃PO₄ terminations)
See THERMODYNAMICS_COMPLETE.md for test results.

Documentation:
--------------
- README.md                    : Complete implementation guide (architecture, patterns)
- WORKGRAPH_CHEATSHEET.md      : Quick reference for AiiDA-WorkGraph
- THERMODYNAMICS_COMPLETE.md   : AIAT implementation details and examples
- AIAT_IMPLEMENTATION.md       : Thermodynamics integration guide

Files:
------
- slabs_relax.py    : CLI driver (supports --mock and --with-thermodynamics)
- workgraph.py      : Main workflow using scatter-gather pattern
- aiat_ternary.py   : Ab initio atomistic thermodynamics module

Usage:
------
# Quick test with mock data
python slabs_relax.py --mock --with-thermodynamics

# Full workflow with VASP
source ~/envs/psteros/bin/activate
python slabs_relax.py --with-thermodynamics --sampling 100

For detailed usage, see individual documentation files.
"""

from .workgraph import (
    slab_relaxation_scatter_gather,
    build_pythonic_workgraph,
)

from .aiat_ternary import (
    calculate_surface_energy_ternary,
    compute_surface_energies_scatter,
    create_mock_bulk_energy,
    create_mock_reference_energies,
    create_mock_formation_enthalpy,
)

__all__ = [
    'slab_relaxation_scatter_gather',
    'build_pythonic_workgraph',
    'calculate_surface_energy_ternary',
    'compute_surface_energies_scatter',
    'create_mock_bulk_energy',
    'create_mock_reference_energies',
    'create_mock_formation_enthalpy',
]

__version__ = '1.0.0'
__author__ = 'PS-TEROS Development Team'
