"""VASP calculations for absolute surface energies of polar faces.

:class:`PolarSurfaceStudy` collects every structure the method needs (bulk,
elemental references, the slabs of each face on their shared passivated
bottom, the pseudo-molecules or tetrahedral clusters for the pseudo-hydrogen,
and optional consistency checks), the VASP settings that go with them, and
turns the finished energies into a surface phase diagram:

    study = PolarSurfaceStudy(bulk, faces=[(1, 1, 1), (-1, -1, -1)],
                              references={"Ga": ga_bulk, "As": as_bulk})
    config = SurfaceWorkflowConfig(backend="vasp", calculation=VaspCalculationConfig(
        ..., potential_mapping=study.potential_mapping({"Ga": "Ga_d", "As": "As"})),
        role_overrides=study.vasp_overrides())
    graph = build_surface_workgraph(study.structures, config, submit=True)
    ...
    result = study.analyse(*read_vasp_results(graph, study.structures))
    result.diagram.plot("gaas_111.png")

Every atom is relaxed at fixed cell, as in Zhang et al. (Sci. Rep. 6, 20055
(2016)); the bottom of each face is then checked to be the same in every
slab.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from psteros.config import CalculationOverride
from psteros.surface_study import (
    GAMMA_ONLY_SPACING,
    GAS_REFERENCES,
    hkl_text as _hkl_text,
    merge_overrides,
    reference_energies_per_atom,
    reservoir_labels,
)

#: INCAR of every slab with two different faces: dipole correction along c,
#: centred on the slab (which psteros places at the middle of the cell).
ASYMMETRIC_SLAB_INCAR = {"ISIF": 2, "LDIPOL": True, "IDIPOL": 3, "DIPOL": [0.5, 0.5, 0.5]}


@dataclass
class PolarSurfaceStudy:
    """Structures, settings and analysis of a polar-surface calculation set.

    Args:
        bulk: Relaxed bulk (pymatgen ``Structure``); the slabs are cut from it
            and its energy is taken at this cell.
        faces: Top faces to study, e.g. ``[(1, 1, 1), (-1, -1, -1)]``. Each face
            gets its own shared bottom.
        references: Structure of the reference phase of each element: the
            elemental solid (Ga, As, Zn) or a molecule in a box (O2, N2).
        bilayers: Slab thickness (9 in the Sci. Rep. paper).
        pseudo_hydrogen_method: ``"molecules"`` (default) or ``"clusters"``.
        cluster_sizes: Tetrahedral cluster sizes for the cluster method.
        eq7_check: Add the slab with both faces passivated (Sci. Rep. Eq. 7).
        nonpolar_check: A non-polar orientation, e.g. ``(1, 1, 0)``, computed
            both as a symmetric slab and with a passivated bottom; the two
            energies must agree.
        electron_counting: Also build the tops that satisfy electron counting.
        oxidation_states, vacuum: Passed to the slab builders.
    """

    bulk: Any
    faces: Sequence[Sequence[int]]
    references: Mapping[str, Any]
    bilayers: int = 9
    pseudo_hydrogen_method: str = "molecules"
    cluster_sizes: Sequence[int] = (2, 3, 8, 9)
    eq7_check: bool = True
    nonpolar_check: Sequence[int] | None = None
    electron_counting: bool = True
    oxidation_states: Mapping[str, float] | None = None
    vacuum: float = 15.0
    structures: dict = field(init=False)
    roles: dict = field(init=False)
    face_sets: dict = field(init=False)

    def __post_init__(self) -> None:
        from psteros.polar import (
            doubly_passivated_slab,
            find_polar_terminations,
            pseudo_hydrogens,
            pseudo_molecule,
            tetrahedral_cluster,
        )

        if self.pseudo_hydrogen_method not in ("molecules", "clusters"):
            raise ValueError("pseudo_hydrogen_method must be 'molecules' or 'clusters'")
        elements = sorted({site.specie.symbol for site in self.bulk})
        missing = sorted(set(elements).difference(self.references))
        if missing:
            raise ValueError(f"references lacks a reference structure for {missing}")
        if not self.faces:
            raise ValueError("at least one face is required")
        self.formula = self.bulk.composition.reduced_formula
        self.hydrogens = pseudo_hydrogens(self.bulk, self.oxidation_states)
        self.structures, self.roles, self.face_sets = {}, {}, {}

        self._add("bulk", self.bulk, "bulk")
        for element in elements:
            self._add(f"ref_{element}", self.references[element], "reference")
        passivated = set()
        for miller in self.faces:
            miller = tuple(int(v) for v in miller)
            slabs = find_polar_terminations(
                self.bulk, miller, bilayers=self.bilayers, vacuum=self.vacuum,
                oxidation_states=self.oxidation_states, electron_counting=self.electron_counting,
            )
            prefix = f"{self.formula}_{_hkl_text(miller)}"
            self.face_sets[prefix] = slabs
            for termination in slabs:
                self._add(f"{prefix}_{termination.label}", termination.structure, "slab")
                passivated.update(termination.pseudo_hydrogen_counts)
        if self.eq7_check:
            miller = tuple(int(v) for v in self.faces[0])
            both = doubly_passivated_slab(self.bulk, miller, bilayers=self.bilayers, vacuum=self.vacuum,
                                          oxidation_states=self.oxidation_states)
            self._add(f"{self.formula}_{_hkl_text(miller)}_both_passivated", both, "eq7_check")
            passivated.update(p for p in both.site_properties["pseudo_hydrogen"] if p)
        if self.nonpolar_check is not None:
            from psteros.core.terminations import find_charge_neutral_terminations

            miller = tuple(int(v) for v in self.nonpolar_check)
            prefix = f"{self.formula}_{_hkl_text(miller)}"
            one_sided = find_polar_terminations(self.bulk, miller, vacuum=self.vacuum,
                                                oxidation_states=self.oxidation_states, electron_counting=False)
            if one_sided.polar:
                raise ValueError(f"nonpolar_check {miller} is a polar direction")
            reference = one_sided[0]
            # The symmetric slab as thick as the passivated one (without its pseudo-H).
            atoms = [site.coords for site, p in zip(reference.structure, reference.structure.site_properties[
                "pseudo_hydrogen"]) if p is None]
            normal = reference.structure.lattice.matrix[2] / reference.structure.lattice.c
            heights = [float(x @ normal) for x in atoms]
            symmetric = [t for t in find_charge_neutral_terminations(
                self.bulk, miller, max(heights) - min(heights) - 0.1, self.vacuum,
                oxidation_states=self.oxidation_states) if t.is_stoichiometric]
            if not symmetric:
                raise ValueError(f"no stoichiometric symmetric slab for {miller}; choose another nonpolar_check")
            self.face_sets[f"{prefix}_check"] = one_sided
            self._add(f"{prefix}_symmetric", symmetric[0].structure, "nonpolar_symmetric")
            self._add(f"{prefix}_passivated", reference.structure, "nonpolar_passivated")
            passivated.update(reference.pseudo_hydrogen_counts)
        for element in sorted(passivated):
            if self.pseudo_hydrogen_method == "molecules":
                self._add(f"pseudo_molecule_{element}", pseudo_molecule(self.bulk, element), "pseudo_molecule")
            else:
                for size in self.cluster_sizes:
                    self._add(f"cluster_{element}_n{size}", tetrahedral_cluster(self.bulk, element, int(size)),
                              "cluster")

    def _add(self, label: str, structure: Any, role: str) -> None:
        if label in self.structures:
            raise ValueError(f"duplicate calculation label {label!r}")
        self.structures[label] = structure
        self.roles[label] = role

    # ------------------------------------------------------------------ VASP

    def potential_mapping(self, base: Mapping[str, str] | None = None) -> dict[str, str]:
        """POTCAR of every kind: ``base`` for the elements, ``H.75``, ``H1.25``... for pseudo-H."""

        from psteros.polar import _potcar_for_kind

        mapping = dict(base or {})
        for structure in self.structures.values():
            kinds = structure.site_properties.get("kind_name") or [site.specie.symbol for site in structure]
            passivates = structure.site_properties.get("pseudo_hydrogen") or [None] * len(structure)
            for kind, site, passivated in zip(kinds, structure, passivates):
                if passivated is not None:
                    mapping.setdefault(kind, _potcar_for_kind(kind))
                else:
                    mapping.setdefault(kind, site.specie.symbol)
        return mapping

    def vasp_overrides(
        self,
        stage: str = "relax",
        *,
        gamma_only_spacing: float = GAMMA_ONLY_SPACING,
        extra: Mapping[str, CalculationOverride] | None = None,
    ) -> dict[str, CalculationOverride]:
        """Per-structure INCAR and k-point overrides for :func:`psteros.build_surface_workgraph`.

        * slabs with two different faces: fixed cell and dipole correction
          (``ISIF=2``, ``LDIPOL``, ``IDIPOL=3``, ``DIPOL`` at the slab centre);
        * the bulk and the symmetric check slab: fixed cell;
        * solid references: cell relaxation (``ISIF=3``);
        * molecules, clusters and gas references: Gamma only, fixed cell;
          O2 is spin-polarised (triplet).

        ``extra`` overrides by label are applied on top.

        ``stage="static"`` gives the overrides of the static calculations of
        :func:`psteros.build_relax_static_workgraph`: the same without ``ISIF``.
        """

        if stage not in ("relax", "static"):
            raise ValueError("stage must be 'relax' or 'static'")
        overrides: dict[str, CalculationOverride] = {}
        for label, role in self.roles.items():
            incar: dict[str, Any] = {}
            kpoints = None
            if role in ("slab", "eq7_check", "nonpolar_passivated"):
                incar.update(ASYMMETRIC_SLAB_INCAR)
            elif role in ("bulk", "nonpolar_symmetric"):
                incar["ISIF"] = 2
            elif role == "reference":
                element = label[len("ref_"):]
                if element in GAS_REFERENCES:
                    incar["ISIF"] = 2
                    kpoints = gamma_only_spacing
                    if element == "O":
                        count = len(self.structures[label])
                        incar.update({"ISPIN": 2, "MAGMOM": [1.0] * count})
                else:
                    incar["ISIF"] = 3
            elif role in ("pseudo_molecule", "cluster"):
                incar["ISIF"] = 2
                kpoints = gamma_only_spacing
            if stage == "static":
                incar.pop("ISIF", None)
            overrides[label] = CalculationOverride(parameters={"INCAR": incar}, kpoints_distance=kpoints)
        return merge_overrides(overrides, extra)

    # -------------------------------------------------------------- analysis

    def binary_references(self, energies_ev: Mapping[str, float]):
        """:class:`psteros.BinaryReferences` from the bulk and reference energies."""

        from psteros.phase_diagram import BinaryReferences

        elements = sorted({site.specie.symbol for site in self.bulk})
        references = reference_energies_per_atom(elements, self.references, energies_ev)
        labels = reservoir_labels(elements)
        return BinaryReferences(
            bulk_energy_ev=energies_ev["bulk"], bulk_composition=self.bulk.composition,
            reference_energies_per_atom_ev=references, reservoir_labels=labels,
        )

    def pseudo_hydrogen_references(self, energies_ev: Mapping[str, float], references=None):
        """:class:`psteros.PseudoHydrogenReferences` from the molecule or cluster energies."""

        from psteros.polar import PseudoHydrogenReferences, fit_cluster_pseudo_chemical_potentials

        if self.pseudo_hydrogen_method == "molecules":
            return PseudoHydrogenReferences.from_pseudo_molecules({
                label[len("pseudo_molecule_"):]: energies_ev[label]
                for label, role in self.roles.items() if role == "pseudo_molecule"
            })
        references = references or self.binary_references(energies_ev)
        mu = references.chemical_potentials_ev(0.0)
        result = PseudoHydrogenReferences()
        elements = sorted({label.split("_")[1] for label, role in self.roles.items() if role == "cluster"})
        for element in elements:
            energies = {int(label.rsplit("_n", 1)[1]): energies_ev[label]
                        for label, role in self.roles.items() if role == "cluster" and label.split("_")[1] == element}
            result[element] = fit_cluster_pseudo_chemical_potentials(element, energies, mu[element]).reference
        return result

    def analyse(
        self,
        energies_ev: Mapping[str, float],
        relaxed_structures: Mapping[str, Any],
        *,
        delta_mu_range: tuple[float, float] | None = None,
        points: int = 201,
        **check_options,
    ) -> "PolarStudyResult":
        """Phase diagram of all faces, with the bottom check of every face.

        Slabs whose bottom changed during relaxation are left out (see
        ``result.bottom_checks``).
        """

        from psteros.phase_diagram import surface_phase_diagram
        from psteros.polar import check_bottoms

        missing = sorted(set(self.structures).difference(energies_ev))
        if missing:
            raise ValueError(f"no energy for {missing}")
        references = self.binary_references(energies_ev)
        hydrogen = self.pseudo_hydrogen_references(energies_ev, references)
        terminations, reports = [], {}
        from dataclasses import replace

        from psteros.phase_diagram import SlabTermination

        for prefix, slabs in self.face_sets.items():
            if prefix.endswith("_check"):
                continue
            relaxed = {t.label: relaxed_structures[f"{prefix}_{t.label}"] for t in slabs}
            report = check_bottoms(slabs, relaxed, **check_options)
            reports[prefix] = report
            for termination in slabs:
                if termination.label in report.passed:
                    terminations.append(SlabTermination.from_polar(
                        termination, energies_ev[f"{prefix}_{termination.label}"],
                        label=f"{prefix}_{termination.label}",
                    ))
        diagram = surface_phase_diagram(terminations, references, pseudo_hydrogen=hydrogen,
                                        delta_mu_range=delta_mu_range, points=points)
        return PolarStudyResult(references, hydrogen, diagram, reports,
                                self._consistency(energies_ev, references, hydrogen))

    def _consistency(self, energies_ev, references, hydrogen) -> tuple:
        """Eq. 7 and non-polar checks for the calculations that are in the set."""

        from psteros.phase_diagram import SlabTermination
        from psteros.polar import eq7_check, nonpolar_check

        checks = []
        for label, role in self.roles.items():
            if role == "eq7_check":
                checks.append(eq7_check(self.structures[label], energies_ev[label], references, hydrogen))
        symmetric = [label for label, role in self.roles.items() if role == "nonpolar_symmetric"]
        if symmetric:
            label = symmetric[0]
            prefix = label[: -len("_symmetric")]
            structure = self.structures[label]
            one_sided = self.face_sets[f"{prefix}_check"][0]
            checks.append(nonpolar_check(
                SlabTermination.from_structure(label, energies_ev[label], structure),
                SlabTermination.from_polar(one_sided, energies_ev[f"{prefix}_passivated"]),
                references, hydrogen,
            ))
        return tuple(checks)


@dataclass(frozen=True)
class PolarStudyResult:
    """References, pseudo chemical potentials, the phase diagram and the bottom checks."""

    references: Any
    pseudo_hydrogen: Any
    diagram: Any
    bottom_checks: Mapping[str, Any]
    consistency: tuple = ()

    def summary(self) -> str:
        lines = [f"{self.references.formula}: Delta mu_{self.references.variable} from "
                 f"{self.references.poor_limit_ev:.3f} to 0 eV"]
        for element, reference in self.pseudo_hydrogen.items():
            lines.append(f"  muhat(H on {element}) = {reference.constant_ev:.4f} eV - mu_{element}/4 ({reference.method})")
        for prefix, report in self.bottom_checks.items():
            lines.append(f"  {prefix}: bottom check {'passed' if report.all_passed else 'FAILED for ' + ', '.join(report.failed)}")
        for check in self.consistency:
            lines.append(f"  {check.summary()}")
        return "\n".join(lines)

    __str__ = summary


def read_vasp_results(graph: Any, labels: Sequence[str]) -> tuple[dict[str, float], dict[str, Any]]:
    """Energies (eV) and relaxed structures of a finished VASP psteros graph.

    Works for :func:`psteros.build_relax_static_workgraph` (energy of the
    static calculation, structure of the relaxation) and for
    :func:`psteros.build_surface_workgraph` (one calculation per label).
    ``graph`` is the WorkGraph process node (``orm.load_node(pk)``); the
    energy is ``energy_extrapolated`` of each calculation's ``misc`` output.
    """

    from psteros.surface_study import _graph_output

    energies, structures = {}, {}
    outputs = graph.outputs
    for label in labels:
        misc = _graph_output(outputs, (f"{label}_static_misc", f"{label}_misc"))
        if misc is None:
            raise ValueError(f"no misc output for {label}")
        values = misc.get_dict() if hasattr(misc, "get_dict") else dict(misc)
        values = values.get("total_energies", values)
        for key in ("energy_extrapolated", "energy_no_entropy", "energy"):
            if key in values:
                energies[label] = float(values[key])
                break
        else:
            raise ValueError(f"no total energy in the misc output of {label}")
        structure = _graph_output(outputs, (f"{label}_relaxed_structure", f"{label}_structure"))
        if structure is not None:
            structures[label] = structure
    return energies, structures
