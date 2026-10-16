"""Gamma-point finite-difference vibrations as an AiiDA WorkGraph.

:func:`build_vibrations_workgraph` takes relaxed structures and a static
(SCF) recipe, displaces every free site by ``+d`` and ``-d`` along x, y and
z, runs one static calculation per displacement (VASP or QE, with the same
recipe as the energies) and turns the forces into harmonic modes.  Each
label gives a ``<label>_vibrations`` output that :func:`read_vibrations`
reads as :class:`psteros.HarmonicVibrations`.

The sites kept fixed in the relaxation (``CalculationOverride.fixed_sites``
of the recipe) are not displaced: a slab with a frozen centre gets the
partial Hessian of its free sites.  A bulk is calculated in a Gamma-only
supercell (``VibrationsConfig.supercells``).  A future implementation may
add the bulk vibrations from phonopy on a q-point mesh.

The builders that already exist are not changed; this is a separate graph,
usually run after the relaxation graph on its relaxed structures.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride, SurfaceWorkflowConfig

#: Ry/bohr in eV/A, for the forces printed by pw.x.
EV_PER_ANGSTROM_PER_RY_PER_BOHR = 13.605_693_122_994 / 0.529_177_210_903

AXES = "xyz"


@dataclass(frozen=True)
class VibrationsConfig:
    """How :func:`build_vibrations_workgraph` calculates the harmonic modes.

    ``displacement_angstrom``
        Finite-difference step ``d`` (central differences, ``+d`` and ``-d``).
    ``displaced_sites``
        Per label, the indices of the sites to displace. By default every
        site except the ``fixed_sites`` of the recipe's override for that
        label (the frozen centre of a slab).
    ``supercells``
        Per label, a Gamma-only supercell ``(na, nb, nc)`` for a bulk; every
        site of the supercell is displaced and the energies stay per input
        cell. The electronic k-points follow the supercell through the
        recipe's k-point spacing.
    ``molecules``
        Labels of gas molecules (O2, H2O, ...): their 5 (linear) or 6
        rotational and translational zero modes are removed. A periodic cell
        with every site displaced loses its 3 translations; a slab with
        frozen sites has no zero mode.
    """

    displacement_angstrom: float = 0.01
    displaced_sites: Mapping[str, tuple[int, ...]] = field(default_factory=dict)
    supercells: Mapping[str, tuple[int, int, int]] = field(default_factory=dict)
    molecules: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.displacement_angstrom > 0:
            raise ValueError("displacement_angstrom must be positive")
        if self.displacement_angstrom > 0.1:
            raise ValueError(
                f"displacement_angstrom={self.displacement_angstrom} A is too large for the harmonic "
                "approximation; use 0.005-0.03 A"
            )
        sites = {}
        for label, indices in dict(self.displaced_sites).items():
            indices = tuple(indices)
            if not indices or any(
                not isinstance(index, int) or isinstance(index, bool) or index < 0 for index in indices
            ):
                raise ValueError(f"displaced_sites[{label!r}] must be non-negative site indices")
            if len(set(indices)) != len(indices):
                raise ValueError(f"displaced_sites[{label!r}] repeats a site")
            sites[label] = tuple(sorted(indices))
        object.__setattr__(self, "displaced_sites", sites)
        supercells = {}
        for label, size in dict(self.supercells).items():
            size = tuple(size)
            if len(size) != 3 or any(not isinstance(n, int) or isinstance(n, bool) or n <= 0 for n in size):
                raise ValueError(f"supercells[{label!r}] must be three positive integers, got {size!r}")
            supercells[label] = size
        object.__setattr__(self, "supercells", supercells)
        object.__setattr__(self, "molecules", tuple(self.molecules))
        for label in self.molecules:
            if label in supercells:
                raise ValueError(f"{label!r} cannot be both a molecule and a supercell")
        for label in supercells:
            if label in sites:
                raise ValueError(f"{label!r}: a supercell displaces every site; do not give displaced_sites")

    def sites_to_displace(
        self, label: str, number_of_sites: int, override: CalculationOverride | None = None
    ) -> tuple[int, ...]:
        """Indices of the sites of ``label`` that are displaced (before any supercell)."""

        if label in self.displaced_sites:
            sites = self.displaced_sites[label]
            outside = [index for index in sites if index >= number_of_sites]
            if outside:
                raise ValueError(f"displaced_sites[{label!r}] has indices outside [0, {number_of_sites}): {outside}")
            return sites
        fixed = set(override.fixed_sites) if override else set()
        if label in self.supercells and fixed:
            raise ValueError(f"{label!r}: a supercell displaces every site, but the recipe fixes sites {sorted(fixed)}")
        sites = tuple(index for index in range(number_of_sites) if index not in fixed)
        if not sites:
            raise ValueError(f"{label!r}: every site is fixed, nothing to displace")
        return sites


def build_vibrations_workgraph(
    structures: Mapping[str, Any],
    static: SurfaceWorkflowConfig,
    vibrations: VibrationsConfig | None = None,
    *,
    submit: bool = False,
) -> Any:
    """Build or submit a WorkGraph of Gamma-point finite-difference vibrations.

    Parameters
    ----------
    structures
        Mapping from labels to *relaxed* structures (pymatgen ``Structure``,
        AiiDA ``StructureData`` or PK), e.g. the ``<label>_relaxed_structure``
        outputs of :func:`psteros.build_relax_static_workgraph`.
    static
        Static recipe of the displaced calculations, normally the one of the
        energies with a tight electronic convergence (VASP ``EDIFF <= 1e-7``,
        QE ``conv_thr <= 1e-10``). VASP INCARs must not relax; QE runs with
        ``tprnfor`` switched on. ``role_overrides`` apply per label as usual,
        and their ``fixed_sites`` select the sites that are *not* displaced.
    vibrations
        :class:`VibrationsConfig` (displacement, supercells, molecules).
    submit
        Submit the graph when true; otherwise return it unsubmitted.

    The graph runs ``6 x (displaced sites)`` static calculations per label,
    at most ``static.execution.max_concurrent_jobs`` at once (one by
    default; raise it to run them in parallel). Graph outputs:
    ``<label>_vibrations`` (a ``Dict``, see :func:`read_vibrations`).
    """

    vibrations = vibrations or VibrationsConfig()
    if not structures:
        raise ValueError("at least one labelled structure is required")
    for label in structures:
        if not label or not label.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"invalid calculation label: {label!r}")
    unknown = sorted(
        (set(vibrations.displaced_sites) | set(vibrations.supercells) | set(vibrations.molecules)).difference(structures)
    )
    if unknown:
        raise ValueError(f"VibrationsConfig names labels that are not in structures: {unknown}")
    if static.backend == "vasp":
        from psteros.backends.vasp import vasp_incar
        from psteros.workflow import _relaxation_check

        for label in structures:
            _relaxation_check(vasp_incar(static.calculation, static.role_overrides.get(label)), label, "static")

    from aiida import orm
    from aiida_workgraph import WorkGraph

    from psteros.backends.qe import as_aiida_structure
    from psteros.backends.vibration_tasks import displace_site, harmonic_modes, repeat_structure

    workgraph = WorkGraph(name=f"{static.name}_vibrations")
    if static.execution.max_concurrent_jobs is not None:
        workgraph.max_number_jobs = static.execution.max_concurrent_jobs
    for label, source in structures.items():
        structure = as_aiida_structure(source)
        override = static.role_overrides.get(label)
        sites = vibrations.sites_to_displace(label, len(structure.sites), override)
        size = vibrations.supercells.get(label, (1, 1, 1))
        supercell_size = size[0] * size[1] * size[2]
        displaced_structure = structure
        if supercell_size > 1:
            displaced_structure = workgraph.add_task(
                repeat_structure, name=f"{label}_supercell", structure=structure, size=orm.List(list(size))
            ).outputs.result
            sites = tuple(range(len(structure.sites) * supercell_size))
        calculation_override = _displacement_override(override, static.backend)

        retrieved = {}
        for site in sites:
            for axis in range(3):
                for sign, name in ((1.0, "plus"), (-1.0, "minus")):
                    key = displacement_key(site, axis, name)
                    displaced = workgraph.add_task(
                        displace_site,
                        name=f"{label}_displace_{key}",
                        structure=displaced_structure,
                        site=orm.Int(site),
                        vector=orm.List([sign * vibrations.displacement_angstrom if a == axis else 0.0 for a in range(3)]),
                    )
                    retrieved[key] = _add_static_task(
                        workgraph, static, f"{label}_vib_{key}", displaced.outputs.result, structure, calculation_override
                    ).outputs.retrieved
        modes = workgraph.add_task(
            harmonic_modes,
            name=f"{label}_vibrations",
            structure=displaced_structure,
            settings=orm.Dict({
                "backend": static.backend,
                "displaced_sites": list(sites),
                "displacement_angstrom": vibrations.displacement_angstrom,
                "supercell_size": supercell_size,
                "molecule": label in vibrations.molecules,
            }),
            retrieved=retrieved,
        )
        workgraph.outputs.__setattr__(f"{label}_vibrations", modes.outputs.result)
    if submit:
        workgraph.submit()
    return workgraph


def displacement_key(site: int, axis: int, sign: str) -> str:
    """Name of one displacement, e.g. ``s12_x_plus``."""

    return f"s{site}_{AXES[axis]}_{sign}"


def _displacement_override(override: CalculationOverride | None, backend: str) -> CalculationOverride:
    """The recipe override of a displaced calculation: no fixed sites, and forces printed by QE."""

    override = replace(override, fixed_sites=()) if override else CalculationOverride()
    if backend == "qe":
        parameters = {namelist: dict(values) for namelist, values in dict(override.parameters or {}).items()}
        parameters.setdefault("CONTROL", {})["tprnfor"] = True
        override = replace(override, parameters=parameters)
    return override


def _add_static_task(
    workgraph: Any,
    static: SurfaceWorkflowConfig,
    label: str,
    structure: Any,
    reference: Any,
    override: CalculationOverride,
) -> Any:
    if static.backend == "vasp":
        from psteros.backends import add_vasp_task

        return add_vasp_task(
            workgraph, label=label, structure=structure, kinds_structure=reference,
            config=static.calculation, execution=static.execution, override=override,
        )
    from psteros.backends import add_qe_task

    return add_qe_task(
        workgraph, label=label, structure=structure, pseudo_structure=reference,
        config=static.calculation, execution=static.execution, override=override,
    )


def read_vibrations(
    graph: Any,
    labels: Sequence[str],
    *,
    imaginary_modes: str = "raise",
    low_frequency_cutoff_cm1: float | None = None,
) -> dict[str, Any]:
    """Harmonic modes of a finished :func:`build_vibrations_workgraph` graph.

    ``graph`` is the WorkGraph process node (``orm.load_node(pk)``); the
    result maps each label to a :class:`psteros.HarmonicVibrations`.
    ``imaginary_modes`` and ``low_frequency_cutoff_cm1`` are passed to it.
    """

    from psteros.surface_study import _graph_output
    from psteros.vibrations import HarmonicVibrations

    result = {}
    for label in labels:
        data = _graph_output(graph.outputs, (f"{label}_vibrations",))
        if data is None:
            raise ValueError(f"no vibrations output for {label}")
        try:
            result[label] = HarmonicVibrations.from_dict(
                data, imaginary_modes=imaginary_modes, low_frequency_cutoff_cm1=low_frequency_cutoff_cm1
            )
        except ValueError as error:
            raise ValueError(f"{label}: {error}") from error
    return result


def vasprun_forces(text: str) -> tuple[list[str], list[list[float]]]:
    """Element symbols and forces (eV/A) of the last ionic step in a ``vasprun.xml``."""

    import xml.etree.ElementTree as ElementTree

    root = ElementTree.fromstring(text)
    symbols = [
        row.find("c").text.strip()
        for row in root.findall("./atominfo/array[@name='atoms']/set/rc")
    ]
    calculations = root.findall("calculation")
    for calculation in reversed(calculations):
        varray = calculation.find("varray[@name='forces']")
        if varray is not None:
            forces = [[float(value) for value in row.text.split()] for row in varray.findall("v")]
            break
    else:
        raise ValueError("vasprun.xml has no forces")
    if len(forces) != len(symbols):
        raise ValueError(f"vasprun.xml has {len(forces)} forces for {len(symbols)} atoms")
    return symbols, forces


_QE_FORCE = re.compile(r"atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*(\S+)\s+(\S+)\s+(\S+)")


def qe_output_forces(text: str, number_of_sites: int) -> list[list[float]]:
    """Total forces (eV/A) printed by pw.x (``tprnfor``), in site order."""

    start = text.rfind("Forces acting on atoms")
    if start < 0:
        raise ValueError("the pw.x output has no forces; run the SCF with tprnfor = .true.")
    forces: dict[int, list[float]] = {}
    for match in _QE_FORCE.finditer(text, start):
        index = int(match.group(1)) - 1
        if index in forces:
            break  # the decomposition into contributions that follows the total forces
        forces[index] = [float(match.group(k)) * EV_PER_ANGSTROM_PER_RY_PER_BOHR for k in (2, 3, 4)]
        if len(forces) == number_of_sites:
            break
    if sorted(forces) != list(range(number_of_sites)):
        raise ValueError(f"the pw.x output has forces for {len(forces)} of {number_of_sites} atoms")
    return [forces[index] for index in range(number_of_sites)]
