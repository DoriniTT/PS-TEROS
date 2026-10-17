"""Typed configuration for reproducible psteros calculations.

The public objects in this module deliberately contain no AiiDA nodes.  That
makes a calculation recipe inspectable, serialisable, and testable before it
is submitted to a computer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping


BackendName = Literal["qe", "vasp", "cp2k"]


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
    """Scheduler settings for the calculations of one psteros WorkGraph.

    Every field describes *your* computer, so nothing is assumed: options
    that are not given are left out of the job and the AiiDA computer or the
    scheduler uses its own default.

    ``computer``
        Name of the AiiDA computer, for the record; the code label of the
        calculation recipe selects the computer that runs the job.
    ``queue``
        Queue or partition (AiiDA ``queue_name``), written by the scheduler
        plugin in its own syntax (``#PBS -q``, ``#SBATCH --partition``, ...).
    ``max_concurrent_jobs``
        Largest number of calculations of the graph running at once
        (``None``: no limit). The default of 1 runs them one after the other.
    ``resources``
        AiiDA resources, e.g. ``{"num_machines": 1, "num_mpiprocs_per_machine": 32}``.
    ``max_wallclock_seconds``, ``account``, ``custom_scheduler_commands``
        Passed to AiiDA when given.
    ``prepend_text``
        Shell lines placed in the job script before the executable, such as
        module loads.
    ``extra_options``
        Adds or replaces AiiDA scheduler options, for example
        ``{"import_sys_environment": False}``; it is applied last.
    """

    computer: str | None = None
    queue: str | None = None
    max_concurrent_jobs: int | None = 1
    resources: Mapping[str, int] = field(default_factory=lambda: {"num_machines": 1})
    max_wallclock_seconds: int | None = None
    with_mpi: bool = True
    prepend_text: str = ""
    extra_options: Mapping[str, Any] = field(default_factory=dict)
    account: str | None = None
    custom_scheduler_commands: str = ""

    def __post_init__(self) -> None:
        if self.max_concurrent_jobs is not None:
            _require_positive("max_concurrent_jobs", self.max_concurrent_jobs)
        if self.max_wallclock_seconds is not None:
            _require_positive("max_wallclock_seconds", self.max_wallclock_seconds)
        for name in ("computer", "queue", "account"):
            value = getattr(self, name)
            if value is not None and not str(value).strip():
                raise ValueError(f"{name} must not be empty; leave it out instead")
        if not self.resources:
            raise ValueError("resources must not be empty")
        for key, value in self.resources.items():
            _require_positive(f"resources[{key!r}]", value)
        if any(not isinstance(key, str) or not key for key in self.extra_options):
            raise ValueError("extra_options keys must be AiiDA option names")

    def scheduler_options(self) -> dict[str, Any]:
        """Return AiiDA metadata options without mutating the source recipe."""

        options: dict[str, Any] = {"resources": dict(self.resources), "withmpi": self.with_mpi}
        if self.max_wallclock_seconds is not None:
            options["max_wallclock_seconds"] = self.max_wallclock_seconds
        if self.queue is not None:
            options["queue_name"] = self.queue
        if self.account is not None:
            options["account"] = self.account
        if self.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.custom_scheduler_commands
        if self.prepend_text:
            options["prepend_text"] = self.prepend_text
        options.update(dict(self.extra_options))
        return options


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


class _DefaultKpointsSpacing(float):
    """The default of ``kpoints_spacing``: a float that tells "left out" from "written by the user" (it equals 0.20)."""


_DEFAULT_KPOINTS_SPACING = _DefaultKpointsSpacing(0.20)


@dataclass(frozen=True)
class VaspCalculationConfig:
    """Inputs shared by aiida-vasp ``VaspWorkChain`` calculations.

    ``incar`` is the INCAR, either flat (``{"ENCUT": 520, ...}``) or already
    under an ``"incar"`` namespace beside other aiida-vasp namespaces
    (``{"incar": {...}, "dynamics": {...}}``).
    ``potential_mapping`` maps an element or kind name to its POTCAR
    (``{"Sn": "Sn_d"}``); elements that are not listed use the POTCAR of the
    same name. ``kpoints_spacing`` is in A^-1 with the 2*pi of VASP's ``KSPACING`` (a
    Gamma-centred mesh has ``ceil(|b_i| / kpoints_spacing)`` points along each
    reciprocal vector); the VASP adapter converts it for aiida-vasp.
    ``max_iterations`` limits the restarts of the work chain (aiida-vasp's
    default when ``None``).

    ``precision`` is ``"fast"``, ``"moderate"`` or ``"precise"`` (see
    :mod:`psteros.precision`): the level fills in what ``incar`` and
    ``kpoints_spacing`` leave out (cutoff, k-point spacing, convergence, the
    POTCAR of the elements that ``potential_mapping`` does not name), and what
    you write explicitly wins. Without ``precision`` nothing is filled in, as
    before. Leaving ``incar`` out selects ``precision="moderate"``, so
    ``VaspCalculationConfig("vasp@cluster")`` is a complete recipe; to relax,
    give the ionic tags (``incar={"IBRION": 2, "NSW": 100, "ISIF": 2}``) and the
    level adds the rest, including ``EDIFFG``.
    """

    code_label: str
    incar: Mapping[str, Any] | None = None
    potential_family: str = "PBE"
    potential_mapping: Mapping[str, str] = field(default_factory=dict)
    kpoints_spacing: float = _DEFAULT_KPOINTS_SPACING
    clean_workdir: bool = False
    max_iterations: int | None = None
    precision: str | None = None

    @property
    def effective_kpoints_spacing(self) -> float:
        """``kpoints_spacing``, or the spacing of the precision level when ``precision`` is set and the spacing was left out."""

        if self.precision is not None and isinstance(self.kpoints_spacing, _DefaultKpointsSpacing):
            from psteros.precision import precision_level

            return precision_level(self.precision).kpoints_spacing
        return float(self.kpoints_spacing)

    def __post_init__(self) -> None:
        if not self.code_label:
            raise ValueError("VASP code_label must not be empty")
        if self.incar is None:  # nothing given: the default precision level is the recipe
            object.__setattr__(self, "incar", {})
            if self.precision is None:
                from psteros.precision import DEFAULT_PRECISION

                object.__setattr__(self, "precision", DEFAULT_PRECISION)
        if self.precision is not None:
            from psteros.precision import precision_level

            precision_level(self.precision)
        if not self.incar and self.precision is None:
            raise ValueError("VASP INCAR must not be empty")
        if not self.potential_family:
            raise ValueError("VASP potential_family must not be empty")
        _require_positive("kpoints_spacing", self.kpoints_spacing)
        if self.max_iterations is not None:
            _require_positive("max_iterations", self.max_iterations)
        if not any(str(key).lower() == "incar" for key in self.incar):
            nested = sorted(str(key) for key, value in self.incar.items() if isinstance(value, Mapping))
            if nested:
                raise ValueError(
                    f"VASP incar must be a flat INCAR mapping or hold it under an 'incar' key; "
                    f"{nested} hold nested mappings"
                )


@dataclass(frozen=True)
class Cp2kCalculationConfig:
    """Inputs shared by aiida-cp2k ``Cp2kBaseWorkChain`` calculations (Quickstep, GPW).

    ``parameters`` is the CP2K input as the nested mapping of aiida-cp2k, with
    upper-case section names; it must hold ``FORCE_EVAL/DFT`` (method, XC,
    ``MGRID``, ``SCF`` with smearing for a metal, ...). The run type is
    ``GLOBAL/RUN_TYPE`` (``ENERGY_FORCE`` when left out); a relaxation needs
    ``GEO_OPT`` or ``CELL_OPT`` and its ``MOTION`` section. psteros adds what
    depends on the structure: the cell and coordinates (through aiida-cp2k),
    one ``KIND`` per kind with the names of ``basis_sets`` and ``potentials``
    (element -> name, e.g. ``{"Ru": "DZVP-MOLOPT-PBE-GTH-q16"}``), the k-point
    mesh, the fixed atoms of a relaxation, and, for ``ENERGY_FORCE``, the
    force and density-of-states files that the readers use.

    ``basis_set_file`` and ``potential_file`` are the files CP2K reads, from
    its data directory or from ``input_files`` (file name -> local path,
    copied into every job, which keeps them in the provenance).
    ``kpoints_spacing`` is in A^-1 with the 2*pi, as for VASP: a
    Monkhorst-Pack mesh with ``ceil(|b_i| / kpoints_spacing)`` points along
    each reciprocal vector; ``None`` or a 1x1x1 mesh runs at Gamma only, with
    no ``KPOINTS`` section. ``max_iterations`` limits the restarts of the
    work chain.
    """

    code_label: str
    parameters: Mapping[str, Any]
    basis_sets: Mapping[str, str]
    potentials: Mapping[str, str]
    basis_set_file: str = "BASIS_MOLOPT"
    potential_file: str = "GTH_POTENTIALS"
    kpoints_spacing: float | None = 0.20
    input_files: Mapping[str, str] = field(default_factory=dict)
    max_iterations: int = 3
    clean_workdir: bool = False

    def __post_init__(self) -> None:
        if not self.code_label:
            raise ValueError("CP2K code_label must not be empty")
        force_eval = self.parameters.get("FORCE_EVAL") if self.parameters else None
        if not isinstance(force_eval, Mapping) or not isinstance(force_eval.get("DFT"), Mapping):
            raise ValueError("CP2K parameters must hold FORCE_EVAL/DFT (upper-case section names)")
        if "KPOINTS" in force_eval["DFT"]:
            raise ValueError("give the CP2K k-points with kpoints_spacing, not FORCE_EVAL/DFT/KPOINTS")
        placed = sorted({"CELL", "COORD", "TOPOLOGY"}.intersection(force_eval.get("SUBSYS", {})))
        if placed:
            raise ValueError(f"CP2K FORCE_EVAL/SUBSYS must not set {placed}: psteros writes the structure")
        if not self.basis_sets or not self.potentials:
            raise ValueError("CP2K basis_sets and potentials must name a basis set and a potential per element")
        missing = sorted(set(self.basis_sets).symmetric_difference(self.potentials))
        if missing:
            raise ValueError(f"CP2K basis_sets and potentials must cover the same elements; differ for {missing}")
        if not self.basis_set_file or not self.potential_file:
            raise ValueError("CP2K basis_set_file and potential_file must not be empty")
        if self.kpoints_spacing is not None:
            _require_positive("kpoints_spacing", self.kpoints_spacing)
        _require_positive("max_iterations", self.max_iterations)
        bad = sorted(name for name in self.input_files if not name or "/" in name)
        if bad:
            raise ValueError(f"CP2K input_files keys must be plain file names, got {bad}")


@dataclass(frozen=True)
class CalculationOverride:
    """Per-structure changes to a shared calculation recipe.

    ``parameters`` updates the recipe: ``{"INCAR": {...}}`` for VASP, namelists
    such as ``{"SYSTEM": {...}}`` for QE, input sections such as
    ``{"FORCE_EVAL": {"DFT": {"UKS": True}}}`` for CP2K (merged section by
    section). ``fixed_sites`` holds the indices of the sites kept fixed during
    a relaxation (VASP selective dynamics, QE ``FIXED_COORDS``, CP2K
    ``FIXED_ATOMS``), e.g. the central layers of a slab.
    """

    parameters: Mapping[str, Mapping[str, Any]] | None = None
    kpoints_distance: float | None = None
    settings: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)
    fixed_sites: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if self.kpoints_distance is not None:
            _require_positive("kpoints_distance", self.kpoints_distance)
        fixed = tuple(self.fixed_sites)
        if any(not isinstance(index, int) or isinstance(index, bool) or index < 0 for index in fixed):
            raise ValueError("fixed_sites must contain non-negative integer site indices")
        object.__setattr__(self, "fixed_sites", tuple(sorted(set(fixed))))


@dataclass(frozen=True)
class SurfaceWorkflowConfig:
    """Backend-neutral recipe accepted by :func:`psteros.build_surface_workgraph`.

    ``execution`` describes the computer and scheduler (see
    :class:`ExecutionPolicy`). ``role_overrides`` is keyed by the structure
    label passed to the builder.
    It is especially useful for references such as spin-polarised O2 while
    retaining one auditable base recipe for bulk and slab calculations.
    """

    backend: BackendName
    calculation: QeCalculationConfig | VaspCalculationConfig | Cp2kCalculationConfig
    execution: ExecutionPolicy
    name: str = "psteros_surface"
    role_overrides: Mapping[str, CalculationOverride] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.backend not in ("qe", "vasp", "cp2k"):
            raise ValueError(f"backend must be 'qe' or 'vasp' (or 'cp2k' with Cp2kCalculationConfig), got {self.backend!r}")
        if self.backend == "qe" and not isinstance(self.calculation, QeCalculationConfig):
            raise TypeError("backend='qe' requires QeCalculationConfig")
        if self.backend == "vasp" and not isinstance(self.calculation, VaspCalculationConfig):
            raise TypeError("backend='vasp' requires VaspCalculationConfig")
        if self.backend == "cp2k" and not isinstance(self.calculation, Cp2kCalculationConfig):
            raise TypeError("backend='cp2k' requires Cp2kCalculationConfig")
        if not self.name or not self.name.replace("_", "").replace("-", "").isalnum():
            raise ValueError("name must contain letters, numbers, '_' or '-'")
