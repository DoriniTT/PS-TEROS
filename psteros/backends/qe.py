"""Quantum ESPRESSO adapter built on aiida-quantumespresso."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from psteros.config import CalculationOverride, ExecutionPolicy, QeCalculationConfig


def as_aiida_structure(structure: Any) -> Any:
    """Convert a pymatgen structure or PK to StructureData when needed.

    An AiiDA WorkGraph output socket is already a valid deferred input. It
    must not be materialised as a pymatgen structure during graph construction.
    """

    from aiida import orm

    if isinstance(structure, int):
        return orm.load_node(structure)
    if isinstance(structure, orm.StructureData):
        return structure
    if structure.__class__.__module__.startswith("aiida_workgraph.sockets"):
        return structure
    # pymatgen's Slab is a Structure subclass with a different converter
    # dispatch name. Convert it explicitly to a plain Structure without
    # changing lattice, species, or fractional coordinates before AiiDA sees it.
    if structure.__class__.__name__ == "Slab":
        from pymatgen.core import Structure

        structure = Structure(
            lattice=structure.lattice,
            species=structure.species_and_occu,
            coords=structure.frac_coords,
            coords_are_cartesian=False,
        )
    return orm.StructureData(pymatgen=structure)


def _merged_parameters(
    base: dict[str, dict[str, Any]], override: CalculationOverride | None
) -> dict[str, dict[str, Any]]:
    result = deepcopy(base)
    if override and override.parameters:
        for namelist, values in override.parameters.items():
            result.setdefault(namelist, {}).update(values)
    return result


def add_qe_task(
    workgraph: Any,
    *,
    label: str,
    structure: Any,
    config: QeCalculationConfig,
    execution: ExecutionPolicy,
    override: CalculationOverride | None = None,
    pseudo_structure: Any | None = None,
) -> Any:
    """Add one ``PwBaseWorkChain`` task and return it.

    Imports are local so pure configuration/analysis use does not require an
    AiiDA profile.  Pseudopotentials are resolved by aiida-pseudo rather than
    hard-coded into psteros.
    """

    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    structure = as_aiida_structure(structure)
    pseudo_structure = as_aiida_structure(pseudo_structure or structure)
    if not isinstance(pseudo_structure, orm.StructureData):
        raise TypeError(
            "a deferred QE structure requires pseudo_structure to be an "
            "AiiDA StructureData node"
        )
    group = orm.load_group(config.pseudo_family)
    if not hasattr(group, "get_pseudos"):
        raise TypeError(f"{config.pseudo_family!r} is not an AiiDA pseudo family")
    code = orm.load_code(config.code_label)
    distance = override.kpoints_distance if override and override.kpoints_distance else config.kpoints_distance
    metadata = execution.scheduler_options()
    if override:
        metadata.update(dict(override.metadata))
    pw = task(WorkflowFactory("quantumespresso.pw.base"))
    pw_inputs = {
        "code": code,
        "structure": structure,
        "parameters": orm.Dict(dict=_merged_parameters(dict(config.parameters), override)),
        "pseudos": group.get_pseudos(structure=pseudo_structure),
        "metadata": {"options": metadata},
    }
    if override and override.settings:
        # aiida-quantumespresso's ``FIXED_COORDS`` setting is needed for
        # symmetric slabs with frozen central trilayers.  Keep settings
        # per-structure so a bulk/reference calculation cannot inherit them.
        pw_inputs["settings"] = orm.Dict(dict=dict(override.settings))
    inputs = {
        "name": f"{label}_qe",
        "pw": pw_inputs,
        "kpoints_distance": orm.Float(distance),
        "max_iterations": orm.Int(config.max_iterations),
        "clean_workdir": orm.Bool(config.clean_workdir),
    }
    return workgraph.add_task(pw, **inputs)
