"""CP2K adapter built on aiida-cp2k (Quickstep, GPW).

The input assembly and the parsers of the force and density-of-states files
are plain Python, so they can be checked without an AiiDA profile; AiiDA is
imported only while a graph is built or read.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from typing import Any

import numpy as np

from psteros.config import CalculationOverride, Cp2kCalculationConfig, ExecutionPolicy

HARTREE_TO_EV = 27.211386245988  # CODATA 2018
BOHR_TO_ANGSTROM = 0.529177210903  # CODATA 2018
FORCES_FILE = "psteros_forces.dat"
DOS_FILE = "psteros.dos"
RELAX_RUN_TYPES = ("GEO_OPT", "CELL_OPT")
STATIC_RUN_TYPES = ("ENERGY_FORCE", "ENERGY")
_KIND_NAME = re.compile(r"^([A-Z][a-z]?)(\d*)$")


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return deepcopy(value)


def _merge(base: dict[str, Any], update: Mapping[str, Any]) -> dict[str, Any]:
    for key, value in update.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = _plain(value)
    return base


def _recipe(config: Cp2kCalculationConfig, override: CalculationOverride | None) -> dict[str, Any]:
    return _merge(_plain(config.parameters), (override.parameters or {}) if override else {})


def cp2k_run_type(config: Cp2kCalculationConfig, override: CalculationOverride | None = None) -> str:
    """``GLOBAL/RUN_TYPE`` of the recipe plus override, upper case (``ENERGY_FORCE`` when not set)."""

    parameters = _recipe(config, override)
    return str(parameters.get("GLOBAL", {}).get("RUN_TYPE", "ENERGY_FORCE")).upper()


def cp2k_parameters(
    config: Cp2kCalculationConfig,
    override: CalculationOverride | None,
    *,
    kinds: Mapping[str, str],
    label: str = "",
) -> dict[str, Any]:
    """The CP2K input of one calculation, before aiida-cp2k adds the cell, coordinates and k-points.

    ``kinds`` maps each kind name of the structure to its element. The recipe
    is updated section by section by ``override.parameters``; psteros then
    fills in, unless they are already given: ``GLOBAL/RUN_TYPE``
    (``ENERGY_FORCE``) and ``PRINT_LEVEL`` (``MEDIUM``, which the aiida-cp2k
    parser reads), the basis and potential files, one ``KIND`` per kind, and
    for ``ENERGY_FORCE`` the force file (``FORCE_EVAL/PRINT/FORCES``) and the
    density-of-states file (``FORCE_EVAL/DFT/PRINT/DOS``) of the readers. The
    ``fixed_sites`` of the override become ``MOTION/CONSTRAINT/FIXED_ATOMS``
    for ``GEO_OPT`` and ``CELL_OPT``.
    """

    where = f"{label}: " if label else ""
    parameters = _recipe(config, override)
    settings = parameters.setdefault("GLOBAL", {})
    settings.setdefault("RUN_TYPE", "ENERGY_FORCE")
    settings.setdefault("PRINT_LEVEL", "MEDIUM")
    run_type = str(settings["RUN_TYPE"]).upper()
    force_eval = parameters["FORCE_EVAL"]
    dft = force_eval["DFT"]
    dft.setdefault("BASIS_SET_FILE_NAME", config.basis_set_file)
    dft.setdefault("POTENTIAL_FILE_NAME", config.potential_file)

    elements = {}
    for name, symbol in kinds.items():
        match = _KIND_NAME.match(name)
        if match is None or match.group(1) != symbol:
            raise ValueError(f"{where}CP2K kind {name!r} must be its element ({symbol}) followed by digits, if any")
        elements[name] = symbol
    missing = sorted({symbol for symbol in elements.values() if symbol not in config.basis_sets})
    if missing:
        raise ValueError(f"{where}no CP2K basis set or potential for {missing}")
    subsys = force_eval.setdefault("SUBSYS", {})
    if "KIND" not in subsys:
        subsys["KIND"] = [
            {"_": name, "ELEMENT": symbol, "BASIS_SET": config.basis_sets[symbol], "POTENTIAL": config.potentials[symbol]}
            for name, symbol in elements.items()
        ]

    if run_type in STATIC_RUN_TYPES:
        force_eval.setdefault("PRINT", {}).setdefault(
            "FORCES", {"_": "ON", "FILENAME": f"={FORCES_FILE}", "NDIGITS": 12})
        dft.setdefault("PRINT", {}).setdefault(
            "DOS", {"_": "ON", "FILENAME": f"={DOS_FILE}", "DELTA_E": 0.0005, "NDIGITS": 8})
    fixed = tuple(override.fixed_sites) if override else ()
    if fixed and run_type in RELAX_RUN_TYPES:
        constraint = parameters.setdefault("MOTION", {}).setdefault("CONSTRAINT", {})
        if "FIXED_ATOMS" in constraint:
            raise ValueError(f"{where}give fixed atoms either as fixed_sites or in MOTION/CONSTRAINT, not both")
        constraint["FIXED_ATOMS"] = {"LIST": " ".join(str(index + 1) for index in fixed), "COMPONENTS_TO_FIX": "XYZ"}
    return parameters


def cp2k_kpoints_mesh(cell: Sequence[Sequence[float]], spacing: float | None) -> tuple[int, int, int]:
    """Monkhorst-Pack mesh ``ceil(|b_i| / spacing)`` (A^-1 with the 2*pi); ``None`` gives Gamma only."""

    if spacing is None:
        return (1, 1, 1)
    if spacing <= 0:
        raise ValueError(f"kpoints spacing must be positive, got {spacing!r}")
    reciprocal = 2.0 * math.pi * np.linalg.inv(np.asarray(cell, dtype=float)).T
    lengths = np.linalg.norm(reciprocal, axis=1)
    return tuple(max(1, math.ceil(length / spacing - 1e-9)) for length in lengths)


def parse_cp2k_forces(text: str) -> np.ndarray:
    """Forces (eV/A, one row per atom) from a CP2K ``FORCE_EVAL/PRINT/FORCES`` file.

    The unit is read from the header (``ATOMIC FORCES in [a.u.]`` is
    Hartree/bohr), not taken from the input, which CP2K does not always
    follow.
    """

    headers = re.findall(r"ATOMIC FORCES in \[([^\]]+)\]", text)
    if len(headers) != 1:
        raise ValueError(f"expected one ATOMIC FORCES block in the CP2K force file, found {len(headers)}")
    unit = headers[0].strip().lower()
    if unit == "a.u.":
        scale = HARTREE_TO_EV / BOHR_TO_ANGSTROM
    elif unit in ("ev/angstrom", "ev/a"):
        scale = 1.0
    else:
        raise ValueError(f"unknown unit of the CP2K forces: [{headers[0]}]")
    rows = [line.split() for line in text.splitlines() if re.match(r"^\s*\d+\s+\d+\s+[A-Za-z]+\s+\S", line)]
    if not rows or [int(row[0]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("the CP2K force file does not list atoms 1, 2, ... in order")
    if any(len(row) != 6 for row in rows):
        raise ValueError("a row of the CP2K force file does not have three force components")
    forces = np.array([[float(value) for value in row[3:]] for row in rows]) * scale
    if not np.isfinite(forces).all():
        raise ValueError("the CP2K force file holds non-finite forces")
    return forces


@dataclass(frozen=True)
class Cp2kDos:
    """A CP2K density-of-states file (``FORCE_EVAL/DFT/PRINT/DOS``).

    ``states[s][i]`` is the (k-weighted) fraction of the orbitals of spin
    ``s`` in energy bin ``i`` and ``occupations[s][i]`` their summed
    occupation, as CP2K writes them. ``fermi_energy_ev`` is in the header of
    Gamma-point runs only.
    """

    energies_ev: tuple[float, ...]
    states: tuple[tuple[float, ...], ...]
    occupations: tuple[tuple[float, ...], ...]
    fermi_energy_ev: float | None = None

    def filling(self) -> list[np.ndarray]:
        """Per spin, the filling of every bin that holds states: 0 (empty) to 1 (full)."""

        ratios = []
        for states, occupations in zip(self.states, self.occupations):
            states, occupations = np.asarray(states), np.asarray(occupations)
            held = states > 0
            ratios.append(occupations[held] / states[held])
        full = max((float(ratio.max()) for ratio in ratios if ratio.size), default=0.0)
        if full <= 0:
            raise ValueError("the CP2K density of states holds no occupied states")
        return [ratio / full for ratio in ratios]

    def partially_occupied_bins(self, tolerance: float = 0.01) -> int:
        """Number of bins whose filling is between ``tolerance`` and ``1 - tolerance``."""

        if not 0 < tolerance < 0.5:
            raise ValueError("tolerance must be between 0 and 0.5")
        return int(sum(((ratio > tolerance) & (ratio < 1 - tolerance)).sum() for ratio in self.filling()))


def parse_cp2k_dos(text: str) -> Cp2kDos:
    """Read the last block of a CP2K ``.dos`` file (energies converted from Hartree to eV)."""

    starts = [match.start() for match in re.finditer(r"#\s*DOS at iteration step", text)]
    if not starts:
        raise ValueError("no 'DOS at iteration step' block in the CP2K density-of-states file")
    block = text[starts[-1]:]
    header = block.splitlines()[0]
    fermi = re.search(r"E_Fermi\[a\.u\.\]\s*=\s*(\S+)", header)
    rows = []
    for line in block.splitlines()[1:]:
        if line.strip().startswith("#") or not line.strip():
            continue
        rows.append([float(value) for value in line.split()])
    if not rows or len({len(row) for row in rows}) != 1 or len(rows[0]) not in (3, 5):
        raise ValueError("the CP2K density-of-states rows must have 3 (one spin) or 5 (two spins) columns")
    data = np.array(rows)
    spins = (data.shape[1] - 1) // 2
    return Cp2kDos(
        energies_ev=tuple(data[:, 0] * HARTREE_TO_EV),
        states=tuple(tuple(data[:, 1 + 2 * spin]) for spin in range(spins)),
        occupations=tuple(tuple(data[:, 2 + 2 * spin]) for spin in range(spins)),
        fermi_energy_ev=float(fermi.group(1)) * HARTREE_TO_EV if fermi else None,
    )


def add_cp2k_task(
    workgraph: Any,
    *,
    label: str,
    structure: Any,
    config: Cp2kCalculationConfig,
    execution: ExecutionPolicy,
    override: CalculationOverride | None = None,
    kinds_structure: Any | None = None,
) -> Any:
    """Add one aiida-cp2k ``Cp2kBaseWorkChain`` task and return it.

    The input is :func:`cp2k_parameters`; the k-point mesh follows
    ``config.kpoints_spacing`` or ``override.kpoints_distance`` and the cell of
    ``kinds_structure`` (given when ``structure`` is the output of an earlier
    task). ``override.settings`` is merged into the aiida-cp2k ``settings``,
    ``override.metadata`` into the job options. The output is parsed with
    ``cp2k_advanced_parser``; the force and density-of-states files are kept
    in ``retrieved``. The restart handler of aiida-cp2k is switched on.
    """

    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    from psteros.backends.qe import as_aiida_structure

    structure = as_aiida_structure(structure)
    reference = as_aiida_structure(kinds_structure if kinds_structure is not None else structure)
    if not isinstance(reference, orm.StructureData):
        raise TypeError("a deferred CP2K structure requires kinds_structure to be an AiiDA StructureData node")
    kinds = {kind.name: kind.symbol for kind in reference.kinds}
    parameters = cp2k_parameters(config, override, kinds=kinds, label=label)
    metadata = {"parser_name": "cp2k_advanced_parser", **execution.scheduler_options()}
    settings: dict[str, Any] = {"additional_retrieve_list": [FORCES_FILE, DOS_FILE]}
    spacing = config.kpoints_spacing
    if override:
        metadata.update(dict(override.metadata))
        settings.update(dict(override.settings))
        if override.kpoints_distance is not None:
            spacing = override.kpoints_distance
    cp2k_inputs: dict[str, Any] = {
        "code": orm.load_code(config.code_label),
        "structure": structure,
        "parameters": orm.Dict(dict=parameters),
        "settings": orm.Dict(dict=settings),
        "metadata": {"options": metadata},
    }
    mesh = cp2k_kpoints_mesh(reference.cell, spacing)
    if mesh != (1, 1, 1):
        kpoints = orm.KpointsData()
        kpoints.set_kpoints_mesh(list(mesh))
        cp2k_inputs["kpoints"] = kpoints
    if config.input_files:
        cp2k_inputs["file"] = {
            re.sub(r"\W", "_", name): orm.SinglefileData(file=str(path), filename=name)
            for name, path in config.input_files.items()
        }
    inputs = {
        "name": f"{label}_cp2k",
        "cp2k": cp2k_inputs,
        "max_iterations": orm.Int(config.max_iterations),
        "clean_workdir": orm.Bool(config.clean_workdir),
        "handler_overrides": orm.Dict(dict={"restart_incomplete_calculation": {"enabled": True}}),
    }
    return workgraph.add_task(task(WorkflowFactory("cp2k.base")), **inputs)


def _output(outputs: Any, names: Sequence[str]) -> Any:
    for name in names:
        try:
            return outputs[name]
        except (KeyError, AttributeError):
            continue
    return None


def _retrieved_file(graph: Any, label: str, filename: str) -> str:
    retrieved = _output(graph.outputs, (f"{label}_static_retrieved", f"{label}_retrieved"))
    if retrieved is None:
        raise ValueError(f"no retrieved folder for {label}")
    try:
        return retrieved.base.repository.get_object_content(filename)
    except (FileNotFoundError, OSError) as error:
        raise ValueError(f"{label}: {filename} was not retrieved (was it an ENERGY_FORCE run?)") from error


def read_cp2k_results(graph: Any, labels: Sequence[str]) -> tuple[dict[str, float], dict[str, Any]]:
    """Energies (eV) and structures of a finished CP2K psteros graph.

    Works for :func:`psteros.build_cp2k_relax_static_workgraph` (energy of the
    static calculation, structure of the relaxation) and for
    :func:`psteros.build_cp2k_surface_workgraph` (one calculation per label).
    The energy is CP2K's total ``FORCE_EVAL`` energy; with smearing it is the
    free energy, whose derivatives are the forces. ``graph`` is the WorkGraph
    process node (``orm.load_node(pk)``).
    """

    energies, structures = {}, {}
    for label in labels:
        parameters = _output(graph.outputs, (f"{label}_static_parameters", f"{label}_parameters"))
        if parameters is None:
            raise ValueError(f"no output parameters for {label}")
        values = parameters.get_dict() if hasattr(parameters, "get_dict") else dict(parameters)
        if "energy" not in values:
            raise ValueError(f"no total energy in the CP2K output parameters of {label}")
        if values.get("energy_units", "a.u.") != "a.u.":
            raise ValueError(f"{label}: unexpected CP2K energy unit {values['energy_units']!r}")
        energies[label] = float(values["energy"]) * HARTREE_TO_EV
        structure = _output(graph.outputs, (f"{label}_relaxed_structure", f"{label}_structure"))
        if structure is not None:
            structures[label] = structure
    return energies, structures


def read_cp2k_forces(graph: Any, labels: Sequence[str]) -> dict[str, np.ndarray]:
    """Forces (eV/A) of the ``ENERGY_FORCE`` calculations of a CP2K psteros graph, from their force files."""

    return {label: parse_cp2k_forces(_retrieved_file(graph, label, FORCES_FILE)) for label in labels}


def read_cp2k_dos(graph: Any, labels: Sequence[str]) -> dict[str, Cp2kDos]:
    """Densities of states of the ``ENERGY_FORCE`` calculations of a CP2K psteros graph."""

    return {label: parse_cp2k_dos(_retrieved_file(graph, label, DOS_FILE)) for label in labels}
