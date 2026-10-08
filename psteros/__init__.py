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
from .phase_diagram import (
    BinaryOxideReferences,
    BinaryReferences,
    SlabTermination,
    SurfacePhaseDiagram,
    surface_phase_diagram,
)
from .phase_diagram_ternary import (
    CompetingPhase,
    TernaryOxideReferences,
    TernaryReferences,
    TernarySurfacePhaseDiagram,
    ternary_surface_phase_diagram,
)
from .polar import (
    PolarTermination,
    PolarTerminationSet,
    PseudoHydrogen,
    PseudoHydrogenReferences,
    find_polar_terminations,
    fit_cluster_pseudo_chemical_potentials,
    pseudo_hydrogen_charge,
    pseudo_hydrogens,
    pseudo_molecule,
    tetrahedral_cluster,
)
from .thermodynamics import (
    EV_PER_ANGSTROM2_TO_J_PER_M2,
    SurfaceEnergyPoint,
    stable_termination,
    surface_energy_binary_equilibrium,
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
    "BinaryOxideReferences",
    "BinaryReferences",
    "SlabTermination",
    "SurfacePhaseDiagram",
    "surface_phase_diagram",
    "CompetingPhase",
    "TernaryOxideReferences",
    "TernaryReferences",
    "TernarySurfacePhaseDiagram",
    "ternary_surface_phase_diagram",
    "PolarTermination",
    "PolarTerminationSet",
    "PseudoHydrogen",
    "PseudoHydrogenReferences",
    "find_polar_terminations",
    "fit_cluster_pseudo_chemical_potentials",
    "pseudo_hydrogen_charge",
    "pseudo_hydrogens",
    "pseudo_molecule",
    "tetrahedral_cluster",
    "EV_PER_ANGSTROM2_TO_J_PER_M2",
    "SurfaceEnergyPoint",
    "stable_termination",
    "surface_energy_binary_equilibrium",
    "surface_energy_elemental",
    "surface_energy_oxide_equilibrium",
    "build_qe_relax_static_workgraph",
    "build_surface_workgraph",
]
