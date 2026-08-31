"""Typed configuration for reproducible psteros calculations.

The public objects in this module deliberately contain no AiiDA nodes.  That
makes a calculation recipe inspectable, serialisable, and testable before it
is submitted to a computer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping


BackendName = Literal["qe", "vasp"]


def _require_positive(name: str, value: float | int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value!r}")


def qe_fixed_coordinate_flags(
    number_of_sites: int, fixed_site_indices: set[int] | list[int] | tuple[int, ...]
) -> list[list[bool]]:
    """Return ``aiida-quantumespresso`` ``FIXED_COORDS`` flags.

    In that plugin, ``True`` fixes a Cartesian coordinate and is rendered as
    QE's ``0`` positional flag; ``False`` leaves it free and is rendered as
    ``1``.  This helper makes the otherwise easy-to-invert convention explicit
    for constrained slab relaxations.
    """

    if number_of_sites <= 0:
        raise ValueError("number_of_sites must be positive")
    fixed = set(fixed_site_indices)
    if any(not isinstance(index, int) for index in fixed):
        raise TypeError("fixed_site_indices must contain integers")
    invalid = sorted(index for index in fixed if index < 0 or index >= number_of_sites)
    if invalid:
        raise ValueError(
            f"fixed_site_indices contains values outside [0, {number_of_sites}): {invalid}"
        )
    return [
        [True, True, True] if index in fixed else [False, False, False]
        for index in range(number_of_sites)
    ]


@dataclass(frozen=True)
class ExecutionPolicy:
    """Scheduler contract for one psteros WorkGraph.

    The defaults encode the Bohr A100 policy used by the maintained SnO2
    campaign: one machine, one MPI rank, the ``gpu_a100`` queue, and no more
    than one active calculation in a graph.
    """

    computer: str = "bohr"
    queue: str = "gpu_a100"
    max_concurrent_jobs: int = 1
    resources: Mapping[str, int] = field(
        default_factory=lambda: {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
    )
    max_wallclock_seconds: int = 86_400
    with_mpi: bool = True

    def __post_init__(self) -> None:
        if self.max_concurrent_jobs != 1:
            raise ValueError(
                "psteros surface workflows currently require "
                "max_concurrent_jobs=1"
            )
        _require_positive("max_wallclock_seconds", self.max_wallclock_seconds)
        if not self.computer:
            raise ValueError("computer must not be empty")
        if not self.queue:
            raise ValueError("queue must not be empty")
        for key, value in self.resources.items():
            _require_positive(f"resources[{key!r}]", value)

    def scheduler_options(self) -> dict[str, Any]:
        """Return AiiDA metadata options without mutating the source recipe."""

        return {
            "resources": dict(self.resources),
            "max_wallclock_seconds": self.max_wallclock_seconds,
            "withmpi": self.with_mpi,
            # Bohr's ``gpu_a100`` queue is a fixed one-A100 resource. Queue
            # selection is consequently the GPU request itself.
            "custom_scheduler_commands": f"#PBS -q {self.queue}\n#PBS -j oe",
        }


@dataclass(frozen=True)
class QeCalculationConfig:
    """Inputs shared by Quantum ESPRESSO ``PwBaseWorkChain`` calculations."""

    code_label: str
    pseudo_family: str
    parameters: Mapping[str, Mapping[str, Any]]
    kpoints_distance: float = 0.20
    max_iterations: int = 1
    clean_workdir: bool = False

    def __post_init__(self) -> None:
        if not self.code_label:
            raise ValueError("QE code_label must not be empty")
        if not self.pseudo_family:
            raise ValueError("QE pseudo_family must not be empty")
        _require_positive("kpoints_distance", self.kpoints_distance)
        _require_positive("max_iterations", self.max_iterations)
        if not self.parameters:
            raise ValueError("QE parameters must include CONTROL, SYSTEM, and ELECTRONS")
        required = {"CONTROL", "SYSTEM", "ELECTRONS"}
        missing = required.difference(self.parameters)
        if missing:
            raise ValueError(f"QE parameters missing namelists: {sorted(missing)}")
        # QE reads these relaxation convergence controls from ``&CONTROL``.
        # Rejecting a common misplaced spelling prevents a remote job that
        # fails immediately in ``read_namelists`` without useful physics.
        misplaced = {
            str(key).lower()
            for key in self.parameters.get("IONS", {})
        }.intersection({"forc_conv_thr", "etot_conv_thr", "nstep"})
        if misplaced:
            raise ValueError(
                "QE relaxation controls "
                f"{sorted(misplaced)} must be in CONTROL, not IONS"
            )


@dataclass(frozen=True)
class VaspCalculationConfig:
    """Inputs shared by an aiida-vasp ``VaspWorkChain`` calculation."""

    code_label: str
    incar: Mapping[str, Any]
    potential_family: str = "PBE"
    potential_mapping: Mapping[str, str] = field(default_factory=dict)
    kpoints_spacing: float = 0.20
    clean_workdir: bool = False

    def __post_init__(self) -> None:
        if not self.code_label:
            raise ValueError("VASP code_label must not be empty")
        if not self.incar:
            raise ValueError("VASP INCAR must not be empty")
        if not self.potential_family:
            raise ValueError("VASP potential_family must not be empty")
        _require_positive("kpoints_spacing", self.kpoints_spacing)


@dataclass(frozen=True)
class CalculationOverride:
    """Per-structure changes to a shared calculation recipe."""

    parameters: Mapping[str, Mapping[str, Any]] | None = None
    kpoints_distance: float | None = None
    settings: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.kpoints_distance is not None:
            _require_positive("kpoints_distance", self.kpoints_distance)


@dataclass(frozen=True)
class SurfaceWorkflowConfig:
    """Backend-neutral recipe accepted by :func:`psteros.build_surface_workgraph`.

    ``role_overrides`` is keyed by the structure label passed to the builder.
    It is especially useful for references such as spin-polarised O2 while
    retaining one auditable base recipe for bulk and slab calculations.
    """

    backend: BackendName
    calculation: QeCalculationConfig | VaspCalculationConfig
    execution: ExecutionPolicy = field(default_factory=ExecutionPolicy)
    name: str = "psteros_surface"
    role_overrides: Mapping[str, CalculationOverride] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.backend == "qe" and not isinstance(self.calculation, QeCalculationConfig):
            raise TypeError("backend='qe' requires QeCalculationConfig")
        if self.backend == "vasp" and not isinstance(self.calculation, VaspCalculationConfig):
            raise TypeError("backend='vasp' requires VaspCalculationConfig")
        if not self.name or not self.name.replace("_", "").replace("-", "").isalnum():
            raise ValueError("name must contain letters, numbers, '_' or '-'")
