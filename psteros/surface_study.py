"""Calculations for the charge-neutral symmetric slabs of any compound.

:class:`ChargeNeutralSurfaceStudy` cuts the charge-neutral, symmetric
terminations of one or more orientations with
:func:`psteros.find_charge_neutral_terminations`, collects them with the bulk,
the elemental references and (for a ternary compound) the competing phases,
gives the QE or VASP settings of each calculation, and turns the finished
energies into a surface phase diagram:

    study = ChargeNeutralSurfaceStudy(bulk, [(1, 0, 0), (1, 1, 0)],
                                      references={"Ag": ag, "P": p, "O": o2},
                                      competing_phases={"Ag2O": ag2o, "P2O5": p2o5},
                                      unit_bonds={("P", "O"): 1.9})
    relax = SurfaceWorkflowConfig(backend="qe", ..., role_overrides=study.qe_overrides("relax"))
    static = SurfaceWorkflowConfig(backend="qe", ..., role_overrides=study.qe_overrides("static"))
    graph = build_qe_relax_static_workgraph(study.structures, relax, static, submit=True)
    ...
    result = study.analyse(*read_qe_results(graph, study.structures))
    result.diagram.plot("ag3po4_surfaces.png")

The slabs are cut from ``bulk`` at its cell, so ``bulk`` should be relaxed
already; the bulk and the slabs are then relaxed at that cell.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride

#: Elements whose reference is a molecule in a box (half its energy per atom).
GAS_REFERENCES = ("H", "N", "O", "F", "Cl")

#: k-point spacing that gives a Gamma-only mesh for molecules and clusters.
GAMMA_ONLY_SPACING = 10.0

#: QE settings of the triplet O2 reference.
QE_TRIPLET_O2 = {"SYSTEM": {"nspin": 2, "tot_magnetization": 2, "starting_magnetization": {"O": 0.5}}}

#: QE cell relaxation of a solid reference.
QE_VC_RELAX = {"CONTROL": {"calculation": "vc-relax"}, "CELL": {"press_conv_thr": 0.5}}


def hkl_text(miller: Sequence[int]) -> str:
    """Miller index as a calculation label: ``(1, -1, 0)`` -> ``"1m10"``."""

    return "".join(str(v) if v >= 0 else f"m{-v}" for v in miller)


def merge_overrides(
    overrides: Mapping[str, CalculationOverride],
    extra: Mapping[str, CalculationOverride] | None,
) -> dict[str, CalculationOverride]:
    """Apply ``extra`` overrides by label on top of ``overrides``, namelist by namelist."""

    merged = dict(overrides)
    for label, override in dict(extra or {}).items():
        if label not in merged:
            raise ValueError(f"unknown calculation label {label!r}")
        base = merged[label]
        parameters = {name: dict(values) for name, values in dict(base.parameters or {}).items()}
        for name, values in dict(override.parameters or {}).items():
            parameters.setdefault(name, {}).update(values)
        merged[label] = CalculationOverride(
            parameters=parameters,
            kpoints_distance=override.kpoints_distance or base.kpoints_distance,
            settings=override.settings or base.settings,
            metadata=override.metadata or base.metadata,
        )
    return merged


def reference_energies_per_atom(
    elements: Sequence[str], references: Mapping[str, Any], energies_ev: Mapping[str, float]
) -> dict[str, float]:
    """Energy per atom of each ``ref_<element>`` calculation."""

    result = {}
    for element in elements:
        label = f"ref_{element}"
        structure = references[element]
        count = sum(1 for site in structure if site.specie.symbol == element)
        if count != len(structure):
            raise ValueError(f"reference {label} must contain only {element}")
        result[element] = energies_ev[label] / count
    return result


def reservoir_labels(elements: Sequence[str]) -> dict[str, str]:
    """Figure labels of the references: ``1/2 O2`` for gases, ``Ag bulk`` for solids."""

    return {e: ("$\\frac{1}{2}$" + e + "$_2$" if e in GAS_REFERENCES else f"{e} bulk") for e in elements}


@dataclass
class ChargeNeutralSurfaceStudy:
    """Structures, settings and analysis of the charge-neutral slabs of a compound.

    Args:
        bulk: Relaxed bulk (pymatgen ``Structure``); the slabs are cut from it
            and its energy is taken at this cell.
        miller_indices: Orientations, e.g. ``[(1, 0, 0), (1, 1, 0)]``. Polar
            orientations without a charge-neutral symmetric slab are rejected
            (see :func:`psteros.find_polar_terminations` for those).
        references: Structure of the reference phase of each element: the
            elemental solid or a molecule in a box (O2, N2). Needed for every
            element of a compound; an element is its own reference (pass
            ``{}``).
        competing_phases: For a ternary compound, ``{label: structure}`` of the
            bulk phases that bound its stability region (e.g. Ag2O and P2O5
            for Ag3PO4).
        min_slab_thickness: Minimum slab thickness (Å).
        vacuum: Vacuum thickness (Å).
        stoichiometric_only: Keep only the stoichiometric terminations.
        oxidation_states, unit_bonds, supercell: Passed to
            :func:`psteros.find_charge_neutral_terminations`.
        termination_options: Further options of that function.
    """

    bulk: Any
    miller_indices: Sequence[Sequence[int]]
    references: Mapping[str, Any]
    competing_phases: Mapping[str, Any] = field(default_factory=dict)
    min_slab_thickness: float = 12.0
    vacuum: float = 15.0
    stoichiometric_only: bool = False
    oxidation_states: Mapping[str, float] | None = None
    unit_bonds: Any = None
    supercell: Sequence[int] | None = None
    termination_options: Mapping[str, Any] = field(default_factory=dict)
    structures: dict = field(init=False)
    roles: dict = field(init=False)
    termination_sets: dict = field(init=False)

    def __post_init__(self) -> None:
        from psteros.core.terminations import find_charge_neutral_terminations

        self.elements = sorted({site.specie.symbol for site in self.bulk})
        if len(self.elements) > 3:
            raise ValueError(f"{len(self.elements)}-component compounds are not supported (at most ternary)")
        missing = sorted(set(self.elements).difference(self.references)) if len(self.elements) > 1 else []
        if missing:
            raise ValueError(f"references lacks a reference structure for {missing}")
        if self.competing_phases and len(self.elements) != 3:
            raise ValueError("competing_phases are used for ternary compounds only")
        if not self.miller_indices:
            raise ValueError("at least one Miller index is required")
        self.formula = self.bulk.composition.reduced_formula
        self.structures, self.roles, self.termination_sets = {}, {}, {}

        self._add("bulk", self.bulk, "bulk")
        if len(self.elements) > 1:
            for element in self.elements:
                self._add(f"ref_{element}", self.references[element], "reference")
        for label, structure in self.competing_phases.items():
            self._add(f"phase_{label}", structure, "competing_phase")
        for miller in self.miller_indices:
            miller = tuple(int(v) for v in miller)
            terminations = find_charge_neutral_terminations(
                self.bulk, miller, self.min_slab_thickness, self.vacuum,
                oxidation_states=self.oxidation_states, unit_bonds=self.unit_bonds,
                supercell=None if self.supercell is None else tuple(self.supercell),
                **dict(self.termination_options),
            )
            prefix = f"{self.formula}_{hkl_text(miller)}"
            self.termination_sets[prefix] = terminations
            for termination in terminations:
                if self.stoichiometric_only and not termination.is_stoichiometric:
                    continue
                self._add(f"{prefix}_{termination.label}", termination.structure, "slab")
        if not any(role == "slab" for role in self.roles.values()):
            raise ValueError("no slab to calculate (all terminations are non-stoichiometric)")

    def _add(self, label: str, structure: Any, role: str) -> None:
        if not label.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"invalid calculation label {label!r}")
        if label in self.structures:
            raise ValueError(f"duplicate calculation label {label!r}")
        self.structures[label] = structure
        self.roles[label] = role

    @property
    def slabs(self) -> list[str]:
        """Labels of the slab calculations."""

        return [label for label, role in self.roles.items() if role == "slab"]

    def summary(self) -> str:
        lines = [f"{self.formula}: {len(self.structures)} calculations "
                 f"({len(self.slabs)} slabs, {len(self.elements)} references"
                 + (f", {len(self.competing_phases)} competing phases" if self.competing_phases else "") + ")"]
        for prefix, terminations in self.termination_sets.items():
            lines += ["", f"{prefix}:", terminations.summary()]
        return "\n".join(lines)

    __str__ = summary

    # ------------------------------------------------------------ settings

    def _gas(self, label: str) -> str | None:
        element = label[len("ref_"):]
        return element if self.roles[label] == "reference" and element in GAS_REFERENCES else None

    def qe_overrides(
        self,
        stage: str = "relax",
        *,
        gamma_only_spacing: float = GAMMA_ONLY_SPACING,
        extra: Mapping[str, CalculationOverride] | None = None,
    ) -> dict[str, CalculationOverride]:
        """Per-structure QE overrides for one stage of :func:`psteros.build_qe_relax_static_workgraph`.

        ``stage="relax"``: solid references and competing phases ``vc-relax``;
        the bulk and the slabs keep the base (fixed-cell) relaxation. Gas
        references are Gamma only and O2 is a triplet, in both stages.
        For :func:`psteros.build_surface_workgraph` (one calculation per
        structure) use the stage that matches its recipe. ``extra``
        overrides by label are applied on top.
        """

        if stage not in ("relax", "static"):
            raise ValueError("stage must be 'relax' or 'static'")
        overrides: dict[str, CalculationOverride] = {}
        for label, role in self.roles.items():
            parameters: dict[str, dict] = {}
            kpoints = None
            gas = self._gas(label)
            if gas is not None:
                kpoints = gamma_only_spacing
                if gas == "O":
                    parameters = {name: dict(values) for name, values in QE_TRIPLET_O2.items()}
            elif role in ("reference", "competing_phase") and stage == "relax":
                parameters = {name: dict(values) for name, values in QE_VC_RELAX.items()}
            overrides[label] = CalculationOverride(parameters=parameters, kpoints_distance=kpoints)
        return merge_overrides(overrides, extra)

    def vasp_overrides(
        self,
        *,
        gamma_only_spacing: float = GAMMA_ONLY_SPACING,
        extra: Mapping[str, CalculationOverride] | None = None,
    ) -> dict[str, CalculationOverride]:
        """Per-structure INCAR and k-point overrides for :func:`psteros.build_surface_workgraph`.

        The bulk and the slabs relax at fixed cell (``ISIF=2``); solid
        references and competing phases relax their cell (``ISIF=3``); gas
        references are Gamma only, and O2 is spin-polarised (triplet).
        ``extra`` overrides by label are applied on top.
        """

        overrides: dict[str, CalculationOverride] = {}
        for label, role in self.roles.items():
            incar: dict[str, Any] = {"ISIF": 2}
            kpoints = None
            gas = self._gas(label)
            if gas is not None:
                kpoints = gamma_only_spacing
                if gas == "O":
                    incar.update({"ISPIN": 2, "MAGMOM": [1.0] * len(self.structures[label])})
            elif role in ("reference", "competing_phase"):
                incar["ISIF"] = 3
            overrides[label] = CalculationOverride(parameters={"INCAR": incar}, kpoints_distance=kpoints)
        return merge_overrides(overrides, extra)

    def potential_mapping(self, base: Mapping[str, str] | None = None) -> dict[str, str]:
        """POTCAR of every element: ``base`` where given, the element symbol otherwise."""

        mapping = dict(base or {})
        for structure in self.structures.values():
            for site in structure:
                mapping.setdefault(site.specie.symbol, site.specie.symbol)
        return mapping

    # ------------------------------------------------------------ analysis

    def chemical_references(self, energies_ev: Mapping[str, float]):
        """:class:`psteros.BinaryReferences` or :class:`psteros.TernaryReferences` from the energies.

        For an element, its bulk energy per atom (``{element: eV}``).
        """

        if len(self.elements) == 1:
            return {self.elements[0]: energies_ev["bulk"] / len(self.bulk)}
        per_atom = reference_energies_per_atom(self.elements, self.references, energies_ev)
        if len(self.elements) == 2:
            from psteros.phase_diagram import BinaryReferences

            return BinaryReferences(
                bulk_energy_ev=energies_ev["bulk"], bulk_composition=self.bulk.composition,
                reference_energies_per_atom_ev=per_atom, reservoir_labels=reservoir_labels(self.elements),
            )
        from psteros.phase_diagram_ternary import CompetingPhase, TernaryReferences

        phases = tuple(
            CompetingPhase(label, energies_ev[f"phase_{label}"], structure.composition)
            for label, structure in self.competing_phases.items()
        )
        return TernaryReferences(
            bulk_energy_ev=energies_ev["bulk"], bulk_composition=self.bulk.composition,
            reference_energies_per_atom_ev=per_atom, competing_phases=phases,
            reservoir_labels=reservoir_labels(self.elements),
        )

    def analyse(
        self,
        energies_ev: Mapping[str, float],
        relaxed_structures: Mapping[str, Any] | None = None,
        *,
        points: int | None = None,
        **diagram_options,
    ) -> "ChargeNeutralStudyResult":
        """Surface phase diagram of every slab.

        Composition and area are read from the relaxed slabs when given (a
        fixed-cell relaxation keeps both), otherwise from the built slabs.
        A binary compound gives a :class:`psteros.SurfacePhaseDiagram`, a
        ternary one a :class:`psteros.TernarySurfacePhaseDiagram`; for an
        element the surface energies are in ``result.surface_energies_j_per_m2``.
        """

        from psteros.phase_diagram import SlabTermination

        missing = sorted(set(self.structures).difference(energies_ev))
        if missing:
            raise ValueError(f"no energy for {missing}")
        relaxed_structures = relaxed_structures or {}
        references = self.chemical_references(energies_ev)
        terminations = [
            SlabTermination.from_structure(label, energies_ev[label], relaxed_structures.get(label, self.structures[label]))
            for label in self.slabs
        ]
        if points is not None:
            diagram_options["points"] = points
        if len(self.elements) == 1:
            from psteros.thermodynamics import EV_PER_ANGSTROM2_TO_J_PER_M2, surface_energy_elemental

            energies = {
                t.label: surface_energy_elemental(
                    slab_energy_ev=t.slab_energy_ev, stoichiometry=t.composition,
                    chemical_potentials_ev=references, surface_area_angstrom2=t.surface_area_angstrom2,
                ) * EV_PER_ANGSTROM2_TO_J_PER_M2
                for t in terminations
            }
            return ChargeNeutralStudyResult(self.formula, references, None, tuple(terminations), energies)
        if len(self.elements) == 2:
            from psteros.phase_diagram import surface_phase_diagram

            diagram = surface_phase_diagram(terminations, references, **diagram_options)
        else:
            from psteros.phase_diagram_ternary import ternary_surface_phase_diagram

            diagram = ternary_surface_phase_diagram(terminations, references, **diagram_options)
        return ChargeNeutralStudyResult(self.formula, references, diagram, tuple(terminations))


@dataclass(frozen=True)
class ChargeNeutralStudyResult:
    """References, the phase diagram and the terminations that entered it."""

    formula: str
    references: Any
    diagram: Any
    terminations: tuple
    surface_energies_j_per_m2: Mapping[str, float] = field(default_factory=dict)

    def summary(self) -> str:
        if self.diagram is None:
            lines = [f"{self.formula}: surface energies (J/m²)"]
            lines += [f"  {label}: {gamma:.3f}" for label, gamma in self.surface_energies_j_per_m2.items()]
            return "\n".join(lines)
        if hasattr(self.diagram, "transitions"):
            references = self.diagram.references
            grid = self.diagram.delta_mu_ev
            lines = [f"{self.formula}: {len(self.terminations)} terminations, "
                     f"Delta mu_{references.variable} from {grid[0]:.3f} to {grid[-1]:.3f} eV"]
            edges = [grid[0], *(t[0] for t in self.diagram.transitions), grid[-1]]
            stable = [self.diagram.stable[0], *(t[2] for t in self.diagram.transitions)]
            for low, high, label in zip(edges, edges[1:], stable):
                lines.append(f"  {low:7.3f} to {high:7.3f} eV: {label}")
            return "\n".join(lines)
        stable = [label for label, region in self.diagram.regions.items() if region]
        return (f"{self.formula}: {len(self.terminations)} terminations; stable somewhere in the "
                f"stability region: {', '.join(stable) or 'none'}")

    __str__ = summary


def _graph_output(outputs: Any, names: Sequence[str]) -> Any:
    for name in names:
        try:
            return outputs[name]
        except (KeyError, AttributeError):
            continue
    return None


def read_qe_results(graph: Any, labels: Sequence[str]) -> tuple[dict[str, float], dict[str, Any]]:
    """Energies (eV) and relaxed structures of a finished QE psteros graph.

    Works for :func:`psteros.build_qe_relax_static_workgraph` (energy of the
    static SCF, structure of the relaxation) and for
    :func:`psteros.build_surface_workgraph` (one calculation per label).
    ``graph`` is the WorkGraph process node (``orm.load_node(pk)``).
    """

    energies, structures = {}, {}
    outputs = graph.outputs
    for label in labels:
        parameters = _graph_output(outputs, (f"{label}_static_parameters", f"{label}_parameters"))
        if parameters is None:
            raise ValueError(f"no output parameters for {label}")
        values = parameters.get_dict() if hasattr(parameters, "get_dict") else dict(parameters)
        if "energy" not in values:
            raise ValueError(f"no total energy in the output parameters of {label}")
        energies[label] = float(values["energy"])
        structure = _graph_output(outputs, (f"{label}_relaxed_structure", f"{label}_structure"))
        if structure is not None:
            structures[label] = structure
    return energies, structures
