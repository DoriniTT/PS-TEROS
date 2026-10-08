"""VASP adapter built on aiida-vasp, the central psteros backend."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride, ExecutionPolicy, VaspCalculationConfig


def vasp_incar(config: VaspCalculationConfig, override: CalculationOverride | None = None) -> dict[str, Any]:
    """The INCAR of one structure: the shared recipe updated by ``override.parameters["INCAR"]``."""

    incar = dict(config.incar)
    if override and override.parameters:
        for namespace, values in override.parameters.items():
            if str(namespace).upper() != "INCAR":
                raise ValueError(f"VASP overrides take an 'INCAR' namespace, got {namespace!r}")
            incar.update(dict(values))
    return incar


def vasp_parameters(
    config: VaspCalculationConfig,
    override: CalculationOverride | None = None,
    number_of_sites: int | None = None,
) -> dict[str, Any]:
    """The ``parameters`` input of an aiida-vasp ``VaspWorkChain``.

    aiida-vasp expects the INCAR in an ``incar`` namespace; fixed sites
    (``override.fixed_sites``) become selective dynamics in ``dynamics``.
    """

    parameters: dict[str, Any] = {"incar": vasp_incar(config, override)}
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
            mapping[kind] = mapping.get(symbol, symbol)
        else:
            missing.append(kind)
    if missing:
        raise ValueError(f"potential_mapping has no POTCAR for the kinds {sorted(missing)}")
    return mapping


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

    ``override.parameters["INCAR"]`` updates the INCAR of this structure,
    ``override.kpoints_distance`` replaces the k-point spacing (a large value
    gives a Gamma-only mesh for molecules) and ``override.fixed_sites``
    fixes atoms with selective dynamics. ``kinds_structure`` gives the kinds
    and site count when ``structure`` is the output of an earlier task.
    """

    from aiida import orm
    from aiida_workgraph import task

    from psteros.backends.vasp_workchain import PsterosVaspWorkChain

    structure = as_vasp_structure(structure)
    reference = as_vasp_structure(kinds_structure if kinds_structure is not None else structure)
    if not isinstance(reference, orm.StructureData):
        raise TypeError("a deferred VASP structure requires kinds_structure to be an AiiDA StructureData node")
    metadata = execution.scheduler_options()
    if override:
        metadata.update(dict(override.metadata))
    kinds = {kind.name: kind.symbol for kind in reference.kinds}
    inputs = {
        "name": f"{label}_vasp",
        "structure": structure,
        "code": orm.load_code(config.code_label),
        "parameters": orm.Dict(dict=vasp_parameters(config, override, len(reference.sites))),
        "kpoints_spacing": orm.Float(
            override.kpoints_distance if override and override.kpoints_distance else config.kpoints_spacing
        ),
        "potential_family": orm.Str(config.potential_family),
        "potential_mapping": orm.Dict(dict=vasp_potential_mapping(config, kinds)),
        "options": orm.Dict(dict=metadata),
        "clean_workdir": orm.Bool(config.clean_workdir),
    }
    if config.max_iterations is not None:
        inputs["max_iterations"] = orm.Int(config.max_iterations)
    if override and override.settings:
        inputs["settings"] = orm.Dict(dict=dict(override.settings))
    return workgraph.add_task(task(PsterosVaspWorkChain), **inputs)
