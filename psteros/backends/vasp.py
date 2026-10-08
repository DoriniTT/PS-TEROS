"""VASP adapter retained as a tested secondary psteros backend."""

from __future__ import annotations

from typing import Any, Mapping

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


def add_vasp_task(workgraph: Any, *, label: str, structure: Any,
                  config: VaspCalculationConfig, execution: ExecutionPolicy,
                  override: CalculationOverride | None = None) -> Any:
    """Add one standard aiida-vasp workchain task and return it.

    The INCAR reaches aiida-vasp under its ``incar`` namespace, whether the
    recipe gives it flat or already namespaced; an override's ``"INCAR"``
    tags are merged into it.  ``override.kpoints_distance`` replaces the
    recipe's ``kpoints_spacing`` and ``override.settings`` is merged into the
    aiida-vasp ``settings``.  OUTCAR, vasprun.xml, CONTCAR and OSZICAR are
    kept in the ``retrieved`` output.

    Changed after 1.0.0: version 1.0.0 passed a flat INCAR, which aiida-vasp 5
    rejects, put override tags beside it, ignored ``kpoints_distance`` and
    ``settings``, and did not keep OUTCAR or vasprun.xml (see CHANGE.md).
    """

    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    from psteros.backends.vasp_tasks import RETRIEVE

    if isinstance(structure, int):
        structure = orm.load_node(structure)
    elif not isinstance(structure, orm.StructureData):
        structure = orm.StructureData(pymatgen=structure)
    metadata = execution.scheduler_options()
    incar, namespaces = split_recipe_incar(config.incar)
    spacing = config.kpoints_spacing
    settings: dict[str, Any] = {"ADDITIONAL_RETRIEVE_LIST": list(RETRIEVE)}
    if override:
        metadata.update(dict(override.metadata))
        if override.parameters:
            # Overrides name the INCAR as "INCAR", as QE overrides name namelists.
            incar.update(lower_keys(override.parameters.get("INCAR", {})))
        if override.kpoints_distance is not None:
            spacing = override.kpoints_distance
        settings.update(dict(override.settings))
    vasp = task(WorkflowFactory("vasp.v2.vasp"))
    return workgraph.add_task(
        vasp,
        name=f"{label}_vasp",
        structure=structure,
        code=orm.load_code(config.code_label),
        parameters=orm.Dict(dict={**namespaces, "incar": incar}),
        kpoints_spacing=orm.Float(spacing),
        potential_family=orm.Str(config.potential_family),
        potential_mapping=orm.Dict(dict=dict(config.potential_mapping)),
        options=orm.Dict(dict=metadata),
        settings=orm.Dict(dict=settings),
        clean_workdir=orm.Bool(config.clean_workdir),
    )
