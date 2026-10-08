"""Reference calculations (bulks and molecules) as blocks, with their thermochemistry.

:func:`build_vasp_reference_workgraph` runs the same blocks (by default
relax -> static -> vibrations, see :mod:`psteros.blocks`) for every labelled
:class:`ReferenceSystem`.  :func:`reference_results` reads the energies,
structures and frequencies back, also while the graph runs, and
:func:`reference_thermochemistry` turns them into
:class:`~psteros.thermochemistry.IdealGasMolecule` and
:class:`~psteros.thermochemistry.HarmonicSolid` objects whose free energies
can be evaluated at any temperature and pressure.

VASP (aiida-vasp ``vasp.v2.vasp``) is the only backend of these blocks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from psteros.blocks import DEFAULT_BLOCKS, Block, Phase, Relax, Static, Vibrations, check_blocks
from psteros.config import CalculationOverride, SurfaceWorkflowConfig, VaspCalculationConfig

EXTRA = "psteros_references"


@dataclass(frozen=True)
class ReferenceSystem:
    """One reference of a phase diagram: a bulk (``"solid"``) or a molecule (``"gas"``).

    ``override`` changes the recipe for every block of this reference (for
    example ``ISPIN = 2`` and a Gamma-only mesh for triplet O2) and
    ``block_overrides`` for one block, by block name (for example
    ``{"relax": CalculationOverride(parameters={"INCAR": {"isif": 3}})}``
    for a bulk cell relaxation).  VASP overrides put INCAR tags under
    ``"INCAR"``.

    ``supercell`` is the repetition of a solid in which its vibrations are
    computed; use a supercell of at least about 10 A per side, because the
    Gamma point of a small cell misses most of the phonons.  A gas needs its
    rotational ``symmetry_number`` (2 for O2, H2, H2O) and total ``spin``
    (1 for triplet O2, 0 for a closed shell).
    """

    structure: Any
    phase: Phase
    override: CalculationOverride | None = None
    block_overrides: Mapping[str, CalculationOverride] = field(default_factory=dict)
    supercell: tuple[int, int, int] = (1, 1, 1)
    symmetry_number: int | None = None
    spin: float | None = None

    def __post_init__(self) -> None:
        if self.phase not in ("gas", "solid"):
            raise ValueError(f"phase must be 'gas' or 'solid', got {self.phase!r}")
        size = tuple(self.supercell)
        if len(size) != 3 or not all(isinstance(n, int) and not isinstance(n, bool) and n > 0 for n in size):
            raise ValueError(f"supercell must be three positive integers, got {self.supercell!r}")
        object.__setattr__(self, "supercell", size)
        if self.phase == "gas":
            if size != (1, 1, 1):
                raise ValueError("a gas reference is computed in its own box; supercell must be (1, 1, 1)")
            if self.symmetry_number is None or self.spin is None:
                raise ValueError(
                    "a gas reference needs symmetry_number and spin (O2: 2 and 1.0, H2O: 2 and 0.0)"
                )
            if isinstance(self.symmetry_number, bool) or not isinstance(self.symmetry_number, int) or self.symmetry_number <= 0:
                raise ValueError(f"symmetry_number must be a positive integer, got {self.symmetry_number!r}")
            if self.spin < 0 or not float(2 * self.spin).is_integer():
                raise ValueError(f"spin must be a non-negative multiple of 1/2, got {self.spin!r}")
        elif self.symmetry_number is not None or self.spin is not None:
            raise ValueError("symmetry_number and spin describe a gas; leave them unset for a solid")
        for override in (self.override, *self.block_overrides.values()):
            if override is None:
                continue
            if not isinstance(override, CalculationOverride):
                raise TypeError("overrides must be CalculationOverride objects")
            other = set(override.parameters or {}).difference({"INCAR"})
            if other:
                raise ValueError(f"VASP overrides take INCAR tags under 'INCAR', got namespaces {sorted(other)}")


def _lower(mapping: Mapping[str, Any]) -> dict[str, Any]:
    return {str(key).lower(): value for key, value in mapping.items()}


def _recipe_incar(incar: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split a recipe INCAR into (INCAR tags, other aiida-vasp parameter namespaces).

    A recipe may give the INCAR flat or, as aiida-vasp expects it, under
    ``"incar"`` beside other namespaces.
    """

    lowered = _lower(incar)
    if "incar" in lowered:
        inner = lowered.pop("incar")
        return _lower(inner), lowered
    return lowered, {}


def _block_inputs(
    reference: ReferenceSystem, block: Block, base_incar: Mapping[str, Any], config: VaspCalculationConfig
) -> tuple[dict[str, Any], float, dict[str, Any], dict[str, Any]]:
    """INCAR, k-point spacing, scheduler-option and settings overrides of one block."""

    incar = dict(base_incar)
    incar.update(block.defaults(reference.phase))
    incar.update(block.incar)
    spacing = config.kpoints_spacing
    metadata: dict[str, Any] = {}
    settings: dict[str, Any] = {}
    for override in (reference.override, reference.block_overrides.get(block.name)):
        if override is None:
            continue
        incar.update(_lower((override.parameters or {}).get("INCAR", {})))
        if override.kpoints_distance is not None:
            spacing = override.kpoints_distance
        metadata.update(dict(override.metadata))
        settings.update(dict(override.settings))
    incar.update(block.required(reference.phase))
    return incar, spacing, metadata, settings


def build_vasp_reference_workgraph(
    references: Mapping[str, ReferenceSystem],
    config: SurfaceWorkflowConfig,
    *,
    blocks: Sequence[Block] = DEFAULT_BLOCKS,
    submit: bool = False,
) -> Any:
    """Build or submit a VASP WorkGraph running ``blocks`` for every reference.

    Every reference runs the same blocks with the same recipe
    (``config.calculation``, a :class:`~psteros.config.VaspCalculationConfig`),
    changed only by its own overrides, so that all reference energies share
    one set of numerical settings.  ``config.role_overrides`` is not used
    here; per-reference changes live in each :class:`ReferenceSystem`.

    For a label ``<label>`` and a block ``<block>`` the graph has the tasks
    ``<label>_<block>_vasp`` and ``<label>_<block>_energy`` (or
    ``<label>_<block>_frequencies`` for :class:`~psteros.blocks.Vibrations`)
    and the outputs ``<label>_<block>_energy`` (eV, sigma -> 0),
    ``<label>_<block>_misc``, ``<label>_<block>_retrieved``,
    ``<label>_<block>_structure`` (relaxations) and
    ``<label>_<block>_frequencies`` (cm^-1, imaginary negative).

    With ``submit=True`` the graph is submitted and the reference description
    is stored on its node, so :func:`reference_results` and
    :func:`reference_thermochemistry` can read it back by PK.
    """

    if config.backend != "vasp" or not isinstance(config.calculation, VaspCalculationConfig):
        raise TypeError("build_vasp_reference_workgraph requires a VASP recipe (backend='vasp')")
    if not references:
        raise ValueError("at least one labelled reference is required")
    invalid = [label for label in references if not label or not label.replace("_", "").isalnum()]
    if invalid:
        raise ValueError(f"invalid reference labels (letters, digits and '_'): {invalid}")
    for label, reference in references.items():
        if not isinstance(reference, ReferenceSystem):
            raise TypeError(f"{label}: references must be ReferenceSystem objects")
    blocks = check_blocks(blocks)
    names = {block.name for block in blocks}
    for label, reference in references.items():
        unknown = set(reference.block_overrides).difference(names)
        if unknown:
            raise ValueError(f"{label}: block_overrides name unknown blocks {sorted(unknown)}; blocks are {sorted(names)}")
    calculation = config.calculation
    base_incar, namespaces = _recipe_incar(calculation.incar)
    planned = {
        (label, block.name): _block_inputs(reference, block, base_incar, calculation)
        for label, reference in references.items()
        for block in blocks
    }
    for (label, name), (incar, *_rest) in planned.items():
        block = next(block for block in blocks if block.name == name)
        if isinstance(block, Relax) and int(incar.get("nsw", 0) or 0) <= 0:
            raise ValueError(f"{label}: relaxation block {name!r} needs NSW > 0 in its INCAR")

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

    workgraph = WorkGraph(name=f"{config.name}_references")
    workgraph.max_number_jobs = config.execution.max_concurrent_jobs
    code = orm.load_code(calculation.code_label)
    for label, reference in references.items():
        produced: dict[str, Any] = {}
        current = as_aiida_structure(reference.structure)
        for block in blocks:
            source = produced[block.structure_from] if block.structure_from else current
            incar, spacing, metadata, extra_settings = planned[(label, block.name)]
            options = config.execution.scheduler_options()
            options.update(metadata)
            settings = {"ADDITIONAL_RETRIEVE_LIST": list(RETRIEVE)}
            structure = source
            if isinstance(block, Vibrations):
                # A displacement run is no relaxation: aiida-vasp would report
                # it as an unconverged one (more ionic steps than NSW).
                settings["CHECK_IONIC_CONVERGENCE"] = False
                settings["parser_settings"] = {"check_ionic_convergence": False}
                if reference.supercell != (1, 1, 1):
                    structure = workgraph.add_task(
                        make_supercell,
                        name=f"{label}_{block.name}_supercell",
                        structure=source,
                        size=orm.List(list(reference.supercell)),
                    ).outputs.result
            settings.update(extra_settings)
            prefix = f"{label}_{block.name}"
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
            setattr(workgraph.outputs, f"{prefix}_misc", vasp.outputs.misc)
            setattr(workgraph.outputs, f"{prefix}_retrieved", vasp.outputs.retrieved)
            if isinstance(block, Vibrations):
                frequencies = workgraph.add_task(
                    vasp_frequencies, name=f"{prefix}_frequencies", retrieved=vasp.outputs.retrieved
                )
                setattr(workgraph.outputs, f"{prefix}_frequencies", frequencies.outputs.result)
            else:
                energy = workgraph.add_task(
                    vasp_energy, name=f"{prefix}_energy", misc=vasp.outputs.misc, retrieved=vasp.outputs.retrieved
                )
                setattr(workgraph.outputs, f"{prefix}_energy", energy.outputs.result)
            if block.moves_ions:
                setattr(workgraph.outputs, f"{prefix}_structure", vasp.outputs.structure)
                current = vasp.outputs.structure
            else:
                current = source
            produced[block.name] = current
    if submit:
        workgraph.submit()
        workgraph.process.base.extras.set(EXTRA, _description(references, blocks))
    return workgraph


def _description(references: Mapping[str, ReferenceSystem], blocks: Sequence[Block]) -> dict[str, Any]:
    return {
        "blocks": [{"name": block.name, "kind": block.kind} for block in blocks],
        "references": {
            label: {
                "phase": reference.phase,
                "supercell": list(reference.supercell),
                "symmetry_number": reference.symmetry_number,
                "spin": reference.spin,
            }
            for label, reference in references.items()
        },
    }


def _stored_description(node: Any) -> dict[str, Any]:
    description = node.base.extras.get(EXTRA, None)
    if description is None:
        raise ValueError(
            f"node {node.pk} has no {EXTRA!r} extra; it was not submitted by "
            "build_vasp_reference_workgraph(..., submit=True)"
        )
    return description


def _state(process: Any) -> str:
    if process is None:
        return "not started"
    state = process.process_state.value if process.process_state is not None else "unknown"
    if state == "finished" and process.exit_status != 0:
        return f"failed [{process.exit_status}]"
    return state


def reference_results(pk: int) -> dict[str, dict[str, dict[str, Any]]]:
    """``{label: {block: result}}`` of a reference WorkGraph, while it runs or after it ends.

    Each result has ``state``, ``pk`` (of the VASP work chain), ``energy``
    (eV) or ``frequencies`` (cm^-1), ``structure`` (the block's input
    structure, or the relaxed one for a relaxation) and ``misc``.  Values are
    ``None`` until the block produces them.
    """

    from aiida import orm
    from aiida.common.links import LinkType

    node = orm.load_node(pk)
    description = _stored_description(node)
    links = node.base.links.get_outgoing(link_type=(LinkType.CALL_WORK, LinkType.CALL_CALC)).all()
    # Sorting by PK keeps the newest child when a task was re-run.
    children = {link.link_label: link.node for link in sorted(links, key=lambda link: link.node.pk)}
    results: dict[str, dict[str, dict[str, Any]]] = {}
    for label in description["references"]:
        results[label] = {}
        for block in description["blocks"]:
            prefix = f"{label}_{block['name']}"
            vasp = children.get(f"{prefix}_vasp")
            result: dict[str, Any] = {
                "kind": block["kind"],
                "state": _state(vasp),
                "pk": vasp.pk if vasp is not None else None,
                "structure": None,
                "misc": None,
            }
            if vasp is not None:
                if block["kind"] != "relax":
                    result["structure"] = vasp.inputs.structure
                elif "structure" in vasp.outputs:
                    result["structure"] = vasp.outputs.structure
                if "misc" in vasp.outputs:
                    result["misc"] = vasp.outputs.misc.get_dict()
            key = "frequencies" if block["kind"] == "vibrations" else "energy"
            child = children.get(f"{prefix}_{key}")
            value = child.outputs.result if child is not None and "result" in child.outputs else None
            if key == "frequencies":
                result[key] = tuple(value.get_list()) if value is not None else None
            else:
                result[key] = value.value if value is not None else None
            results[label][block["name"]] = result
    return results


def reference_thermochemistry(
    pk: int,
    *,
    energy_block: str | None = None,
    vibrations_block: str | None = None,
    imaginary_tolerance_cm1: float = 0.0,
    corrections_ev: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """``{label: IdealGasMolecule | HarmonicSolid}`` from a finished reference WorkGraph.

    The electronic energy comes from ``energy_block`` (by default the last
    static block, else the last relaxation) and the frequencies from
    ``vibrations_block`` (by default the last vibrations block).  A solid's
    vibrational terms are scaled from the displaced supercell to the cell of
    the energy block.  ``corrections_ev`` adds an explicit correction to the
    named references (for example an O2 binding correction); it is reported
    as its own term.
    """

    from psteros.thermochemistry import HarmonicSolid, IdealGasMolecule

    from aiida import orm

    description = _stored_description(orm.load_node(pk))
    blocks = description["blocks"]

    def last(kinds: tuple[str, ...]) -> str | None:
        names = [block["name"] for block in blocks if block["kind"] in kinds]
        return names[-1] if names else None

    energy_block = energy_block or last(("static",)) or last(("relax",))
    vibrations_block = vibrations_block or last(("vibrations",))
    kinds = {block["name"]: block["kind"] for block in blocks}
    if energy_block not in kinds or kinds[energy_block] == "vibrations":
        raise ValueError(f"energy_block must name a relax or static block, got {energy_block!r}")
    if vibrations_block not in kinds or kinds[vibrations_block] != "vibrations":
        raise ValueError(f"vibrations_block must name a vibrations block, got {vibrations_block!r}")
    corrections = dict(corrections_ev or {})
    unknown = set(corrections).difference(description["references"])
    if unknown:
        raise ValueError(f"corrections_ev names unknown references {sorted(unknown)}")

    results = reference_results(pk)
    systems: dict[str, Any] = {}
    for label, spec in description["references"].items():
        energy = results[label][energy_block]
        vibrations = results[label][vibrations_block]
        if energy["energy"] is None or vibrations["frequencies"] is None:
            raise ValueError(
                f"{label}: blocks {energy_block!r} ({energy['state']}) and "
                f"{vibrations_block!r} ({vibrations['state']}) have not both finished"
            )
        energy_structure = orm.load_node(energy["pk"]).inputs.structure
        vibrated = orm.load_node(vibrations["pk"]).inputs.structure
        if spec["phase"] == "gas":
            systems[label] = IdealGasMolecule.from_structure(
                vibrated,
                electronic_energy_ev=energy["energy"],
                frequencies_cm1=vibrations["frequencies"],
                symmetry_number=spec["symmetry_number"],
                spin=spec["spin"],
                correction_ev=corrections.get(label, 0.0),
                imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            )
        else:
            systems[label] = HarmonicSolid(
                electronic_energy_ev=energy["energy"],
                atoms_in_cell=len(energy_structure.sites),
                frequencies_cm1=vibrations["frequencies"],
                atoms_in_supercell=len(vibrated.sites),
                correction_ev=corrections.get(label, 0.0),
                imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            )
    return systems
