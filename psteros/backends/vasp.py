"""VASP adapter built on aiida-vasp, the central psteros backend."""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride, ExecutionPolicy, VaspCalculationConfig


def lower_keys(mapping: Mapping[str, Any]) -> dict[str, Any]:
    """Return ``mapping`` with lower-case keys (INCAR tags are case-insensitive)."""

    return {str(key).lower(): value for key, value in mapping.items()}


def split_recipe_incar(incar: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split a recipe INCAR into (INCAR tags, other aiida-vasp parameter namespaces).

    A recipe may give the INCAR flat or, as aiida-vasp 5 expects it, under
    ``"incar"`` beside other namespaces such as ``"dynamics"``.
    """

    lowered = lower_keys(incar)
    if "incar" in lowered:
        inner = lowered.pop("incar")
        return lower_keys(inner), lowered
    return lowered, {}


def aiida_vasp_kpoints_spacing(spacing_per_angstrom: float) -> float:
    """Convert a psteros k-point spacing to the value aiida-vasp expects.

    psteros states ``kpoints_spacing`` in A^-1 *including* the 2*pi, as VASP's
    ``KSPACING`` does: a Gamma-centred mesh has ``N_i = ceil(|b_i| / spacing)``
    points along reciprocal vector ``b_i`` (``|b_i| = 2*pi / a`` for a cubic
    cell of edge ``a``). aiida-vasp multiplies its ``kpoints_spacing`` input by
    2*pi before building the mesh, so it must receive ``spacing / (2*pi)``;
    passing the psteros value unchanged gives ``N = ceil(1 / (a * spacing))``,
    which is a single k-point for ``a = 5.75 A`` and ``spacing = 0.2 A^-1``.
    """

    return float(spacing_per_angstrom) / (2.0 * math.pi)


def has_custom_kinds(structure: Any) -> bool:
    """Whether an AiiDA ``StructureData`` has kinds whose name is not their element (``H0p75``, ``Fe1``)."""

    return any(kind.name != kind.symbol for kind in structure.kinds)


def assign_kinds(
    template_symbols: Sequence[str],
    template_fractional: Any,
    relaxed_symbols: Sequence[str],
    relaxed_fractional: Any,
    cell: Any,
    *,
    max_distance: float = 2.0,
) -> list[int]:
    """For every relaxed site, the index of the template site it came from.

    aiida-vasp rebuilds the relaxed structure from element names, so the kind
    names of the starting structure (``H0p75``, ``H1p25``) are lost, and VASP
    may also reorder the atoms. Atoms are matched one to one within each
    element by their minimum-image distance (a relaxation moves them little);
    ``cell`` (rows are lattice vectors) converts fractional differences to A.
    """

    import numpy as np
    from scipy.optimize import linear_sum_assignment

    if sorted(template_symbols) != sorted(relaxed_symbols):
        raise ValueError(
            f"the relaxed structure has different atoms ({sorted(relaxed_symbols)}) "
            f"than the starting structure ({sorted(template_symbols)})"
        )
    cell = np.asarray(cell, dtype=float)
    template_fractional = np.asarray(template_fractional, dtype=float)
    relaxed_fractional = np.asarray(relaxed_fractional, dtype=float)
    assignment: list[int] = [-1] * len(relaxed_symbols)
    for symbol in sorted(set(relaxed_symbols)):
        relaxed_index = [i for i, value in enumerate(relaxed_symbols) if value == symbol]
        template_index = [i for i, value in enumerate(template_symbols) if value == symbol]
        delta = relaxed_fractional[relaxed_index][:, None, :] - template_fractional[template_index][None, :, :]
        delta -= np.round(delta)
        distance = np.linalg.norm(delta @ cell, axis=2)
        for row, column in zip(*linear_sum_assignment(distance)):
            if distance[row, column] > max_distance:
                raise ValueError(
                    f"{symbol} site {relaxed_index[row]} of the relaxed structure is {distance[row, column]:.2f} A "
                    f"from the nearest unmatched {symbol} of the starting structure (limit {max_distance} A)"
                )
            assignment[relaxed_index[row]] = template_index[column]
    return assignment


def recipe_incar(incar: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """(INCAR tags, other aiida-vasp namespaces) of a recipe INCAR, which may be flat or namespaced.

    Like :func:`split_recipe_incar`, but the spelling of the tags is kept.
    """

    for key in incar:
        if str(key).lower() == "incar":
            return dict(incar[key]), {k: v for k, v in incar.items() if k != key}
    return dict(incar), {}


def vasp_incar(
    config: VaspCalculationConfig,
    override: CalculationOverride | None = None,
    site_symbols: Sequence[str] | None = None,
) -> dict[str, Any]:
    """The INCAR of one structure: the shared recipe updated by ``override.parameters["INCAR"]``.

    With ``config.precision`` the tags of the precision level that the recipe and the override leave out come first
    (see :mod:`psteros.precision`; ``site_symbols`` lets the level decide ``ISPIN``); without it nothing is added.
    """

    incar, _ = recipe_incar(config.incar)
    if override and override.parameters:
        for namespace, values in override.parameters.items():
            if str(namespace).upper() != "INCAR":
                raise ValueError(f"VASP overrides take an 'INCAR' namespace, got {namespace!r}")
            incar.update(dict(values))
    if config.precision is not None:
        from psteros.precision import level_incar

        incar = {**level_incar(config.precision, incar, site_symbols=site_symbols), **incar}
    return incar


def vasp_parameters(
    config: VaspCalculationConfig,
    override: CalculationOverride | None = None,
    number_of_sites: int | None = None,
    site_symbols: Sequence[str] | None = None,
) -> dict[str, Any]:
    """The ``parameters`` input of an aiida-vasp ``VaspWorkChain``.

    aiida-vasp expects the INCAR in an ``incar`` namespace with lower-case tags
    (it refuses a stored ``Dict`` with ``ENCUT``); fixed sites
    (``override.fixed_sites``) become selective dynamics in ``dynamics``.
    """

    _, namespaces = recipe_incar(config.incar)
    parameters: dict[str, Any] = {
        **lower_keys(namespaces), "incar": lower_keys(vasp_incar(config, override, site_symbols)),
    }
    if override and override.fixed_sites:
        if number_of_sites is None:
            raise ValueError("fixed_sites needs the number of sites of the structure")
        parameters["dynamics"] = {"positions_dof": vasp_positions_dof(number_of_sites, override.fixed_sites)}
    return parameters


def vasp_positions_dof(number_of_sites: int, fixed_site_indices: Sequence[int]) -> list[list[bool]]:
    """Selective-dynamics flags: ``True`` lets a coordinate relax, ``False`` fixes it (VASP's T/F)."""

    from psteros.config import qe_fixed_coordinate_flags

    return [[not flag for flag in row] for row in qe_fixed_coordinate_flags(number_of_sites, list(fixed_site_indices))]


def vasp_potential_mapping(config: VaspCalculationConfig, kinds: Mapping[str, str]) -> dict[str, str]:
    """POTCAR of every kind: the recipe's mapping, else the element symbol of the kind.

    ``kinds`` maps kind name to element symbol. Kinds whose name is not an
    element (e.g. pseudo-hydrogen ``H0p75``) must be in the recipe's mapping.
    """

    mapping = dict(config.potential_mapping)
    missing = []
    for kind, symbol in kinds.items():
        if kind in mapping:
            continue
        if kind == symbol:
            potcar = mapping.get(symbol)
            if potcar is None and config.precision is not None:
                from psteros.precision import protocol_potential

                potcar = protocol_potential(symbol)  # the aiida-vasp protocol's POTCAR for the element
            mapping[kind] = potcar or symbol
        else:
            missing.append(kind)
    if missing:
        raise ValueError(f"potential_mapping has no POTCAR for the kinds {sorted(missing)}")
    return mapping


def restart_potential_aliases(mapping: Mapping[str, str], kinds: Mapping[str, str]) -> dict[str, str]:
    """``mapping`` plus the element name of each custom kind whose POTCAR is unambiguous.

    aiida-vasp restarts a work chain (after a walltime, or a failed line search)
    from the relaxed structure it returned, whose kinds are element names only,
    and looks the POTCAR of ``H`` up in the mapping: ``KeyError: 'H'`` when the
    structure had a pseudo-hydrogen ``H0p75``. An element with exactly one
    custom kind and no plain kind of that name gets that kind's POTCAR under its
    element name too. With two custom kinds of one element (``H0p75`` and
    ``H1p25`` in one slab) the element name stays out and a restart fails loudly.
    """

    custom: dict[str, list[str]] = {}
    for kind, symbol in kinds.items():
        if kind != symbol:
            custom.setdefault(symbol, []).append(kind)
    result = dict(mapping)
    for symbol, names in custom.items():
        if len(names) == 1 and symbol not in kinds and symbol not in result:
            result[symbol] = mapping[names[0]]
    return result


def as_vasp_structure(structure: Any) -> Any:
    """StructureData from a pymatgen structure or PK; WorkGraph sockets pass through."""

    from aiida import orm

    if isinstance(structure, int):
        return orm.load_node(structure)
    if isinstance(structure, orm.StructureData):
        return structure
    if structure.__class__.__module__.startswith("aiida_workgraph"):
        return structure
    from psteros.backends.qe import as_aiida_structure

    return as_aiida_structure(structure)


def add_vasp_task(
    workgraph: Any,
    *,
    label: str,
    structure: Any,
    config: VaspCalculationConfig,
    execution: ExecutionPolicy,
    override: CalculationOverride | None = None,
    kinds_structure: Any | None = None,
) -> Any:
    """Add one aiida-vasp ``VaspWorkChain`` task and return it.

    The INCAR reaches aiida-vasp under its ``incar`` namespace, whether the
    recipe gives it flat or already namespaced (other namespaces such as
    ``dynamics`` are kept); ``override.parameters["INCAR"]`` updates it for
    this structure. ``override.kpoints_distance`` replaces the k-point spacing
    (a large value gives a Gamma-only mesh for molecules),
    ``override.fixed_sites`` fixes atoms with selective dynamics and
    ``override.settings`` is merged into the aiida-vasp ``settings``. OUTCAR,
    vasprun.xml, CONTCAR and OSZICAR are kept in the ``retrieved`` output.
    ``kinds_structure`` gives the kinds and site count when ``structure`` is
    the output of an earlier task.
    """

    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    from psteros.backends.vasp_tasks import RETRIEVE

    structure = as_vasp_structure(structure)
    reference = as_vasp_structure(kinds_structure if kinds_structure is not None else structure)
    if not isinstance(reference, orm.StructureData):
        raise TypeError("a deferred VASP structure requires kinds_structure to be an AiiDA StructureData node")
    metadata = execution.scheduler_options()
    spacing = config.effective_kpoints_spacing
    settings: dict[str, Any] = {"ADDITIONAL_RETRIEVE_LIST": list(RETRIEVE)}
    if override:
        metadata.update(dict(override.metadata))
        if override.kpoints_distance is not None:
            spacing = override.kpoints_distance
        settings.update(dict(override.settings))
    kinds = {kind.name: kind.symbol for kind in reference.kinds}
    inputs = {
        "name": f"{label}_vasp",
        "structure": structure,
        "code": orm.load_code(config.code_label),
        "parameters": orm.Dict(
            dict=vasp_parameters(config, override, len(reference.sites), [kinds[site.kind_name] for site in reference.sites])
        ),
        # ``spacing`` is the recipe's value or the override's, both in psteros units (A^-1 with 2*pi).
        "kpoints_spacing": orm.Float(aiida_vasp_kpoints_spacing(spacing)),
        "potential_family": orm.Str(config.potential_family),
        "potential_mapping": orm.Dict(dict=restart_potential_aliases(vasp_potential_mapping(config, kinds), kinds)),
        "options": orm.Dict(dict=metadata),
        "settings": orm.Dict(dict=settings),
        "clean_workdir": orm.Bool(config.clean_workdir),
    }
    if config.max_iterations is not None:
        inputs["max_iterations"] = orm.Int(config.max_iterations)
    return workgraph.add_task(task(WorkflowFactory("vasp.v2.vasp")), **inputs)
