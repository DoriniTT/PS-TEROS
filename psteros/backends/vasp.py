"""VASP adapter retained as a tested secondary psteros backend."""

from __future__ import annotations

import math
from typing import Any

from psteros.config import CalculationOverride, ExecutionPolicy, VaspCalculationConfig


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


def add_vasp_task(workgraph: Any, *, label: str, structure: Any,
                  config: VaspCalculationConfig, execution: ExecutionPolicy,
                  override: CalculationOverride | None = None) -> Any:
    """Add one standard aiida-vasp workchain task and return it."""

    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    if isinstance(structure, int):
        structure = orm.load_node(structure)
    elif not isinstance(structure, orm.StructureData):
        structure = orm.StructureData(pymatgen=structure)
    metadata = execution.scheduler_options()
    if override:
        metadata.update(dict(override.metadata))
    parameters = dict(config.incar)
    if override and override.parameters:
        # VASP uses a flat INCAR mapping. Accept one "INCAR" namespace for
        # symmetry with the QE override representation.
        parameters.update(dict(override.parameters.get("INCAR", {})))
    vasp = task(WorkflowFactory("vasp.v2.vasp"))
    return workgraph.add_task(
        vasp,
        name=f"{label}_vasp",
        structure=structure,
        code=orm.load_code(config.code_label),
        parameters=orm.Dict(dict=parameters),
        kpoints_spacing=orm.Float(aiida_vasp_kpoints_spacing(config.kpoints_spacing)),
        potential_family=orm.Str(config.potential_family),
        potential_mapping=orm.Dict(dict=dict(config.potential_mapping)),
        options=orm.Dict(dict=metadata),
        clean_workdir=orm.Bool(config.clean_workdir),
    )
