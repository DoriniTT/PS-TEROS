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
    SlabTermination,
    SurfacePhaseDiagram,
    surface_phase_diagram,
)
from .phase_diagram_ternary import (
    CompetingPhase,
    TernaryOxideReferences,
    TernarySurfacePhaseDiagram,
    ternary_surface_phase_diagram,
)
from .thermodynamics import (
    EV_PER_ANGSTROM2_TO_J_PER_M2,
    SurfaceEnergyPoint,
    stable_termination,
    surface_energy_elemental,
    surface_energy_oxide_equilibrium,
)
from .blocks import Relax, Static, Vibrations
from .references import (
    ReferenceSystem,
    build_vasp_reference_workgraph,
    reference_results,
    reference_thermochemistry,
)
from .thermochemistry import (
    FreeEnergy,
    HarmonicSolid,
    IdealGasMolecule,
    delta_mu_oxygen_ev,
    free_energies,
    oxygen_pressure_bar,
    parse_vasp_frequencies_cm1,
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
    "SlabTermination",
    "SurfacePhaseDiagram",
    "surface_phase_diagram",
    "CompetingPhase",
    "TernaryOxideReferences",
    "TernarySurfacePhaseDiagram",
    "ternary_surface_phase_diagram",
    "EV_PER_ANGSTROM2_TO_J_PER_M2",
    "SurfaceEnergyPoint",
    "stable_termination",
    "surface_energy_elemental",
    "surface_energy_oxide_equilibrium",
    "FreeEnergy",
    "HarmonicSolid",
    "IdealGasMolecule",
    "delta_mu_oxygen_ev",
    "free_energies",
    "oxygen_pressure_bar",
    "parse_vasp_frequencies_cm1",
    "build_qe_relax_static_workgraph",
    "build_surface_workgraph",
    "Relax",
    "Static",
    "Vibrations",
    "ReferenceSystem",
    "build_vasp_reference_workgraph",
    "reference_results",
    "reference_thermochemistry",
]
