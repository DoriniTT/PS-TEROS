"""One WorkGraph for a whole surface campaign, with its outputs grouped by structure.

:func:`build_vasp_campaign_workgraph` runs the reference calculations (bulks
and molecules, by default relax -> static -> vibrations) and the slab
terminations (by default relax -> static) of a phase diagram in a single
graph, with one recipe for all of them.  Its outputs are nested by group,
structure and block::

    node.outputs.references.o2.static.energy        # eV
    node.outputs.references.o2.vibrations.frequencies
    node.outputs.slabs.slab_o.relax.structure
    node.outputs.slabs.slab_o.static.energy

:func:`campaign_results` reads the same values back as plain Python, also
while the graph runs or after it fails, and :func:`campaign_terminations`
turns the slabs into :class:`~psteros.phase_diagram.SlabTermination` objects.
:func:`~psteros.references.reference_results` and
:func:`~psteros.references.reference_thermochemistry` work on a campaign
graph as on a reference graph.

The slabs are given as finished structures (cut from a bulk relaxed
beforehand); the graph does not rebuild them from the bulk it relaxes.
VASP (aiida-vasp ``vasp.v2.vasp``) is the only backend.
"""

from __future__ import annotations

import keyword
import re
from dataclasses import dataclass, field
from typing import Any, ClassVar, Mapping, Sequence

from psteros.backends.vasp import split_recipe_incar
from psteros.blocks import DEFAULT_BLOCKS, Block, Relax, Static, Vibrations, check_blocks
from psteros.config import CalculationOverride, SurfaceWorkflowConfig, VaspCalculationConfig
from psteros.references import EXTRA as REFERENCES_EXTRA
from psteros.references import ReferenceSystem, _block_inputs, _description, _state

EXTRA = "psteros_campaign"

GROUPS = ("references", "slabs")

DEFAULT_SLAB_BLOCKS: tuple[Block, ...] = (Relax(), Static())

# Labels and block names become output namespaces (``slabs.<label>.<block>``)
# and AiiDA link labels, where ``__`` separates namespaces.
_NAMESPACE = re.compile(r"[A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)*", re.ASCII)


def _is_namespace(name: Any) -> bool:
    return isinstance(name, str) and bool(_NAMESPACE.fullmatch(name)) and not keyword.iskeyword(name)


@dataclass(frozen=True)
class SlabSystem:
    """One slab termination of a campaign.

    ``override`` changes the recipe for every block of this slab and
    ``block_overrides`` for one block, by block name (for example a looser
    ``EDIFFG`` for the relaxation).  VASP overrides put INCAR tags under
    ``"INCAR"``.  A slab is computed as a solid, in its own cell.
    """

    structure: Any
    override: CalculationOverride | None = None
    block_overrides: Mapping[str, CalculationOverride] = field(default_factory=dict)

    phase: ClassVar[str] = "solid"

    def __post_init__(self) -> None:
        named = {"override": self.override}
        named.update({f"block_overrides[{name!r}]": value for name, value in self.block_overrides.items()})
        for where, override in named.items():
            if override is None:
                continue
            if not isinstance(override, CalculationOverride):
                raise TypeError(f"{where}: overrides must be CalculationOverride objects, got {type(override).__name__}")
            other = set(override.parameters or {}).difference({"INCAR"})
            if other:
                raise ValueError(f"{where}: VASP overrides take INCAR tags under 'INCAR', got namespaces {sorted(other)}")


def _check_labels(references: Mapping[str, Any], slabs: Mapping[str, Any]) -> None:
    invalid = [label for label in (*references, *slabs) if not _is_namespace(label)]
    if invalid:
        raise ValueError(
            f"invalid labels {invalid}: use letters, digits and single '_' between them, starting with a letter "
            "(Python keywords are not allowed)"
        )
    shared = sorted(set(references).intersection(slabs))
    if shared:
        raise ValueError(f"labels {shared} are used both as reference and as slab; labels must be unique")


def _check_block_names(group: str, blocks: Sequence[Block]) -> None:
    invalid = [block.name for block in blocks if not _is_namespace(block.name)]
    if invalid:
        raise ValueError(
            f"{group} block names {invalid} cannot name an output namespace: use letters, digits and "
            "single '_' between them, starting with a letter"
        )


def build_vasp_campaign_workgraph(
    references: Mapping[str, ReferenceSystem],
    slabs: Mapping[str, SlabSystem | Any],
    config: SurfaceWorkflowConfig,
    *,
    reference_blocks: Sequence[Block] = DEFAULT_BLOCKS,
    slab_blocks: Sequence[Block] = DEFAULT_SLAB_BLOCKS,
    submit: bool = False,
) -> Any:
    """Build or submit one VASP WorkGraph for the references and the slabs of a campaign.

    Every reference runs ``reference_blocks`` and every slab ``slab_blocks``
    with the same recipe (``config.calculation``), changed only by the
    overrides of each :class:`~psteros.references.ReferenceSystem` or
    :class:`SlabSystem`, so that all energies share one set of numerical
    settings.  A bare structure in ``slabs`` stands for ``SlabSystem(structure)``.
    ``config.role_overrides`` is not used; giving it is an error.

    For a label ``<label>`` and a block ``<block>`` the graph has the same
    tasks as :func:`~psteros.references.build_vasp_reference_workgraph`
    (``<label>_<block>_vasp``, ``<label>_<block>_energy`` or
    ``<label>_<block>_frequencies``, ``<label>_<block>_supercell``) and the
    nested outputs ``<group>.<label>.<block>.<port>``, where ``<group>`` is
    ``references`` or ``slabs`` and ``<port>`` is ``energy`` (eV, sigma -> 0)
    or ``frequencies`` (cm^-1, imaginary negative), ``structure``
    (relaxations), ``misc``, ``remote`` and ``retrieved``.

    With ``submit=True`` the graph is submitted and its description is stored
    on its node, so :func:`campaign_results`, :func:`campaign_terminations`,
    :func:`~psteros.references.reference_results` and
    :func:`~psteros.references.reference_thermochemistry` can read it back by PK.
    """

    if config.backend != "vasp" or not isinstance(config.calculation, VaspCalculationConfig):
        raise TypeError("build_vasp_campaign_workgraph requires a VASP recipe (backend='vasp')")
    if config.role_overrides:
        raise ValueError(
            "build_vasp_campaign_workgraph does not use config.role_overrides "
            f"(given for {sorted(config.role_overrides)}); put overrides on the ReferenceSystem or SlabSystem"
        )
    if not references and not slabs:
        raise ValueError("at least one reference or slab is required")
    for label, reference in references.items():
        if not isinstance(reference, ReferenceSystem):
            raise TypeError(f"{label}: references must be ReferenceSystem objects")
    for label, slab in slabs.items():
        if isinstance(slab, (ReferenceSystem, CalculationOverride)):
            raise TypeError(f"{label}: slabs must be SlabSystem objects or structures, got {type(slab).__name__}")
    slabs = {label: slab if isinstance(slab, SlabSystem) else SlabSystem(slab) for label, slab in slabs.items()}
    _check_labels(references, slabs)
    reference_blocks = check_blocks(reference_blocks)
    slab_blocks = check_blocks(slab_blocks)
    vibrations = [block.name for block in slab_blocks if isinstance(block, Vibrations)]
    if vibrations:
        raise ValueError(f"slab blocks {vibrations} are vibrations, which are not supported for slabs")

    groups = {"references": (references, reference_blocks), "slabs": (slabs, slab_blocks)}
    calculation = config.calculation
    base_incar, namespaces = split_recipe_incar(calculation.incar)
    planned = {}
    for group, (systems, blocks) in groups.items():
        _check_block_names(group, blocks)
        names = {block.name for block in blocks}
        for label, system in systems.items():
            unknown = set(system.block_overrides).difference(names)
            if unknown:
                raise ValueError(
                    f"{label}: block_overrides name unknown blocks {sorted(unknown)}; {group} blocks are {sorted(names)}"
                )
            for block in blocks:
                planned[(label, block.name)] = _block_inputs(system, block, base_incar, calculation)
                incar = planned[(label, block.name)][0]
                if isinstance(block, Relax) and int(incar.get("nsw", 0) or 0) <= 0:
                    raise ValueError(f"{label}: relaxation block {block.name!r} needs NSW > 0 in its INCAR")

    from aiida import orm
    from aiida_workgraph import WorkGraph

    from psteros.backends.qe import as_aiida_structure
    from psteros.backends.vasp_tasks import (
        RETRIEVE,
        add_vasp_block_task,
        make_supercell,
        vasp_energy,
        vasp_frequencies,
    )

    workgraph = WorkGraph(name=f"{config.name}_campaign")
    workgraph.max_number_jobs = config.execution.max_concurrent_jobs
    code = orm.load_code(calculation.code_label)
    for group, (systems, blocks) in groups.items():
        for label, system in systems.items():
            produced: dict[str, Any] = {}
            current = as_aiida_structure(system.structure)
            for block in blocks:
                source = produced[block.structure_from] if block.structure_from else current
                incar, spacing, metadata, extra_settings = planned[(label, block.name)]
                options = config.execution.scheduler_options()
                options.update(metadata)
                settings = {"ADDITIONAL_RETRIEVE_LIST": list(RETRIEVE)}
                structure = source
                prefix = f"{label}_{block.name}"
                if isinstance(block, Vibrations):
                    # A displacement run is no relaxation: aiida-vasp would report
                    # it as an unconverged one (more ionic steps than NSW).
                    settings["CHECK_IONIC_CONVERGENCE"] = False
                    settings["parser_settings"] = {"check_ionic_convergence": False}
                    if system.supercell != (1, 1, 1):
                        structure = workgraph.add_task(
                            make_supercell,
                            name=f"{prefix}_supercell",
                            structure=source,
                            size=orm.List(list(system.supercell)),
                        ).outputs.result
                settings.update(extra_settings)
                vasp = add_vasp_block_task(
                    workgraph,
                    name=f"{prefix}_vasp",
                    structure=structure,
                    code=code,
                    incar=incar,
                    namespaces=namespaces,
                    kpoints_spacing=spacing,
                    potential_family=calculation.potential_family,
                    potential_mapping=calculation.potential_mapping,
                    options=options,
                    settings=settings,
                    clean_workdir=calculation.clean_workdir,
                )
                outputs = {
                    "misc": vasp.outputs.misc,
                    "remote": vasp.outputs.remote_folder,
                    "retrieved": vasp.outputs.retrieved,
                }
                if isinstance(block, Vibrations):
                    frequencies = workgraph.add_task(
                        vasp_frequencies, name=f"{prefix}_frequencies", retrieved=vasp.outputs.retrieved
                    )
                    outputs["frequencies"] = frequencies.outputs.result
                else:
                    energy = workgraph.add_task(
                        vasp_energy, name=f"{prefix}_energy", misc=vasp.outputs.misc, retrieved=vasp.outputs.retrieved
                    )
                    outputs["energy"] = energy.outputs.result
                if block.moves_ions:
                    outputs["structure"] = vasp.outputs.structure
                    current = vasp.outputs.structure
                else:
                    current = source
                produced[block.name] = current
                for port, socket in outputs.items():
                    setattr(workgraph.outputs, f"{group}.{label}.{block.name}.{port}", socket)
    if submit:
        workgraph.submit()
        extras = {EXTRA: _campaign_description(references, reference_blocks, slabs, slab_blocks)}
        if references:
            extras[REFERENCES_EXTRA] = _description(references, reference_blocks)
        workgraph.process.base.extras.set_many(extras)
    return workgraph


def _campaign_description(
    references: Mapping[str, ReferenceSystem],
    reference_blocks: Sequence[Block],
    slabs: Mapping[str, SlabSystem],
    slab_blocks: Sequence[Block],
) -> dict[str, Any]:
    def blocks(sequence: Sequence[Block]) -> list[dict[str, str]]:
        return [{"name": block.name, "kind": block.kind} for block in sequence]

    return {
        "version": 1,
        "references": {"blocks": blocks(reference_blocks), "labels": list(references)},
        "slabs": {"blocks": blocks(slab_blocks), "labels": list(slabs)},
    }


def _load_campaign(pk: int) -> tuple[Any, dict[str, Any]]:
    from aiida import orm

    node = orm.load_node(pk)
    description = node.base.extras.get(EXTRA, None)
    if description is None:
        raise ValueError(
            f"node {node.pk} has no {EXTRA!r} extra; it was not submitted by "
            "build_vasp_campaign_workgraph(..., submit=True)"
        )
    return node, description


def _block_result(children: Mapping[str, Any], prefix: str, kind: str) -> dict[str, Any]:
    vasp = children.get(f"{prefix}_vasp")
    result: dict[str, Any] = {
        "kind": kind,
        "state": _state(vasp),
        "pk": vasp.pk if vasp is not None else None,
        "structure": None,
        "misc": None,
        "remote": None,
        "retrieved": None,
    }
    if vasp is not None:
        if kind != "relax":
            result["structure"] = vasp.inputs.structure
        elif "structure" in vasp.outputs:
            result["structure"] = vasp.outputs.structure
        if "misc" in vasp.outputs:
            result["misc"] = vasp.outputs.misc.get_dict()
        if "remote_folder" in vasp.outputs:
            result["remote"] = vasp.outputs.remote_folder
        if "retrieved" in vasp.outputs:
            result["retrieved"] = vasp.outputs.retrieved
    key = "frequencies" if kind == "vibrations" else "energy"
    child = children.get(f"{prefix}_{key}")
    value = child.outputs.result if child is not None and "result" in child.outputs else None
    if key == "frequencies":
        result[key] = tuple(value.get_list()) if value is not None else None
    else:
        result[key] = value.value if value is not None else None
    return result


def campaign_results(pk: int) -> dict[str, dict[str, dict[str, dict[str, Any]]]]:
    """``{group: {label: {block: result}}}`` of a campaign WorkGraph, while it runs or after it ends.

    ``group`` is ``"references"`` or ``"slabs"``.  Each result has ``kind``,
    ``state``, ``pk`` (of the VASP work chain), ``energy`` (eV) or
    ``frequencies`` (cm^-1), ``structure`` (the block's input structure, or
    the relaxed one for a relaxation), ``misc`` (a dict), ``remote`` and
    ``retrieved``.  Values are ``None`` until the block produces them.
    """

    from aiida.common.links import LinkType

    node, description = _load_campaign(pk)
    links = node.base.links.get_outgoing(link_type=(LinkType.CALL_WORK, LinkType.CALL_CALC)).all()
    # Sorting by PK keeps the newest child when a task was re-run.
    children = {link.link_label: link.node for link in sorted(links, key=lambda link: link.node.pk)}
    return {
        group: {
            label: {
                block["name"]: _block_result(children, f"{label}_{block['name']}", block["kind"])
                for block in description[group]["blocks"]
            }
            for label in description[group]["labels"]
        }
        for group in GROUPS
    }


def campaign_terminations(pk: int, *, energy_block: str | None = None, surfaces: int = 2) -> list[Any]:
    """One :class:`~psteros.phase_diagram.SlabTermination` per slab of a finished campaign graph.

    The energy (eV) comes from ``energy_block`` (by default the last static
    slab block, else the last relaxation) and the composition and area from
    the structure that block computed.  ``surfaces`` is passed to
    :meth:`~psteros.phase_diagram.SlabTermination.from_structure`.
    """

    from psteros.phase_diagram import SlabTermination

    _node, description = _load_campaign(pk)
    blocks = description["slabs"]["blocks"]
    kinds = {block["name"]: block["kind"] for block in blocks}

    def last(kind: str) -> str | None:
        names = [block["name"] for block in blocks if block["kind"] == kind]
        return names[-1] if names else None

    energy_block = energy_block or last("static") or last("relax")
    if energy_block not in kinds or kinds[energy_block] == "vibrations":
        raise ValueError(f"energy_block must name a relax or static slab block, got {energy_block!r}")
    if not description["slabs"]["labels"]:
        raise ValueError(f"campaign graph {pk} has no slabs")
    terminations = []
    for label, results in campaign_results(pk)["slabs"].items():
        result = results[energy_block]
        if result["energy"] is None or result["structure"] is None:
            raise ValueError(f"{label}: block {energy_block!r} has not finished ({result['state']})")
        terminations.append(SlabTermination.from_structure(label, result["energy"], result["structure"], surfaces=surfaces))
    return terminations
