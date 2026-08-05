"""psteros: reproducible surface thermodynamics with AiiDA.

Quantum ESPRESSO through ``aiida-quantumespresso`` is the primary backend.
The VASP adapter remains available for established VASP workflows.
"""

from .config import (
    CalculationOverride,
    ExecutionPolicy,
    QeCalculationConfig,
    SurfaceWorkflowConfig,
    VaspCalculationConfig,
    qe_fixed_coordinate_flags,
)
from .structures import (
    SlabIdentity,
    alpha_sn_bulk,
    litharge_sno_bulk,
    rutile_sno2_bulk,
    sno2_110_slab,
    triplet_o2_cell,
)
from .thermodynamics import (
    EV_PER_ANGSTROM2_TO_J_PER_M2,
    SurfaceEnergyPoint,
    stable_termination,
    surface_energy_elemental,
    surface_energy_oxide_equilibrium,
)
from .workflow import build_qe_relax_static_workgraph, build_surface_workgraph

__version__ = "1.0.0"

__all__ = [
    "CalculationOverride",
    "ExecutionPolicy",
    "QeCalculationConfig",
    "SurfaceWorkflowConfig",
    "VaspCalculationConfig",
    "qe_fixed_coordinate_flags",
    "SlabIdentity",
    "alpha_sn_bulk",
    "litharge_sno_bulk",
    "rutile_sno2_bulk",
    "sno2_110_slab",
    "triplet_o2_cell",
    "EV_PER_ANGSTROM2_TO_J_PER_M2",
    "SurfaceEnergyPoint",
    "stable_termination",
    "surface_energy_elemental",
    "surface_energy_oxide_equilibrium",
    "build_qe_relax_static_workgraph",
    "build_surface_workgraph",
]
