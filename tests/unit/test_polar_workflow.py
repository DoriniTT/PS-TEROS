"""Tier-1/2 tests for the VASP calculation set of polar surfaces."""

from __future__ import annotations

import pytest

pytest.importorskip("pymatgen")

from pymatgen.core import Lattice, Structure  # noqa: E402

import psteros  # noqa: E402
from psteros import polar  # noqa: E402
from psteros.polar_workflow import PolarSurfaceStudy, read_vasp_results  # noqa: E402


def _sg(group, lattice, species, coords):
    return Structure.from_spacegroup(group, lattice, species, coords)


def gaas():
    return _sg("F-43m", Lattice.cubic(5.653), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def ga_metal():
    return _sg("Cmce", Lattice.orthorhombic(4.52, 7.66, 4.53), ["Ga"], [[0, 0.1549, 0.081]])


def as_solid():
    return _sg("R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]])


def zno():
    return _sg("P6_3mc", Lattice.hexagonal(3.25, 5.21), ["Zn", "O"], [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])


def o2_molecule():
    return Structure(Lattice.cubic(12.0), ["O", "O"], [[6, 6, 5.4], [6, 6, 6.6]], coords_are_cartesian=True)


def zn_metal():
    return _sg("P6_3/mmc", Lattice.hexagonal(2.66, 4.95), ["Zn"], [[1 / 3, 2 / 3, 0.25]])


def gaas_study(**options):
    values = dict(faces=[(1, 1, 1), (-1, -1, -1)], references={"Ga": ga_metal(), "As": as_solid()}, bilayers=3)
    values.update(options)
    return PolarSurfaceStudy(gaas(), **values)


def fake_energies(study):
    """Energies from a per-atom model (any consistent numbers will do)."""
    per_atom = {"Ga": -3.0, "As": -4.7, "Zn": -1.3, "O": -4.9}
    energies = {}
    for label, structure in study.structures.items():
        energy = 0.0
        for site, passivates in zip(structure, structure.site_properties.get("pseudo_hydrogen") or [None] * len(structure)):
            energy += -1.1 if passivates else per_atom[site.specie.symbol]
        if label == "bulk":
            energy -= 0.4 * len(structure)  # makes the compound stable
        energies[label] = energy
    return energies


def test_study_collects_every_calculation():
    study = gaas_study(nonpolar_check=(1, 1, 0))
    roles = study.roles
    assert roles["bulk"] == "bulk" and roles["ref_Ga"] == roles["ref_As"] == "reference"
    assert [label for label, role in roles.items() if role == "slab"] == [
        "GaAs_111_term_0", "GaAs_111_term_1", "GaAs_m1m1m1_term_0", "GaAs_m1m1m1_term_1",
    ]
    assert roles["GaAs_111_both_passivated"] == "eq7_check"
    assert roles["GaAs_110_symmetric"] == "nonpolar_symmetric"
    assert roles["GaAs_110_passivated"] == "nonpolar_passivated"
    assert {label for label, role in roles.items() if role == "pseudo_molecule"} == {
        "pseudo_molecule_As", "pseudo_molecule_Ga",
    }
    symmetric = study.structures["GaAs_110_symmetric"].composition
    passivated = study.structures["GaAs_110_passivated"].composition
    assert symmetric["Ga"] == passivated["Ga"] and symmetric["As"] == passivated["As"]
    for label in study.structures:
        assert label.replace("_", "").replace("-", "").isalnum()


def test_potential_mapping_adds_the_pseudo_hydrogen_potcars():
    mapping = gaas_study().potential_mapping({"Ga": "Ga_d"})
    assert mapping == {"Ga": "Ga_d", "As": "As", "H0p75": "H.75", "H1p25": "H1.25"}


def test_vasp_overrides_by_role():
    study = gaas_study()
    overrides = study.vasp_overrides()
    slab = overrides["GaAs_111_term_0"]
    assert dict(slab.parameters["INCAR"]) == {"ISIF": 2, "LDIPOL": True, "IDIPOL": 3, "DIPOL": [0.5, 0.5, 0.5]}
    assert slab.kpoints_distance is None
    assert dict(overrides["bulk"].parameters["INCAR"]) == {"ISIF": 2}
    assert dict(overrides["ref_Ga"].parameters["INCAR"]) == {"ISIF": 3}
    molecule = overrides["pseudo_molecule_As"]
    assert dict(molecule.parameters["INCAR"]) == {"ISIF": 2} and molecule.kpoints_distance == 10.0
    tuned = study.vasp_overrides(extra={"bulk": psteros.CalculationOverride(parameters={"INCAR": {"NSW": 0}})})
    assert dict(tuned["bulk"].parameters["INCAR"]) == {"ISIF": 2, "NSW": 0}
    with pytest.raises(ValueError, match="unknown calculation label"):
        study.vasp_overrides(extra={"nothing": psteros.CalculationOverride()})


def test_oxide_with_an_o2_reference_and_clusters():
    study = PolarSurfaceStudy(
        zno(), faces=[(0, 0, 1)], references={"Zn": zn_metal(), "O": o2_molecule()}, bilayers=2,
        pseudo_hydrogen_method="clusters", cluster_sizes=(2, 3, 4, 5), eq7_check=False,
    )
    overrides = study.vasp_overrides()
    o2 = overrides["ref_O"]
    assert dict(o2.parameters["INCAR"]) == {"ISIF": 2, "ISPIN": 2, "MAGMOM": [1.0, 1.0]}
    assert o2.kpoints_distance == 10.0
    assert sorted(label for label, role in study.roles.items() if role == "cluster") == [
        "cluster_O_n2", "cluster_O_n3", "cluster_O_n4", "cluster_O_n5",
    ]
    references = study.binary_references(fake_energies(study))
    assert references.variable == "O" and references.reservoir_labels["O"].startswith("$\\frac{1}{2}$")


def test_study_rejects_missing_references_and_bad_options():
    with pytest.raises(ValueError, match="lacks a reference"):
        PolarSurfaceStudy(gaas(), faces=[(1, 1, 1)], references={"Ga": ga_metal()})
    with pytest.raises(ValueError, match="pseudo_hydrogen_method"):
        gaas_study(pseudo_hydrogen_method="wedges")
    with pytest.raises(ValueError, match="polar direction"):
        gaas_study(nonpolar_check=(1, 1, 1))


def test_analysis_builds_the_phase_diagram_and_checks_bottoms():
    study = gaas_study()
    energies = fake_energies(study)
    relaxed = dict(study.structures)
    result = study.analyse(energies, relaxed, points=5)
    assert set(result.diagram.curves) == {
        "GaAs_111_term_0", "GaAs_111_term_1", "GaAs_m1m1m1_term_0", "GaAs_m1m1m1_term_1",
    }
    assert result.pseudo_hydrogen["As"].constant_ev == pytest.approx(energies["pseudo_molecule_As"] / 4)
    assert all(report.all_passed for report in result.bottom_checks.values())
    assert "bottom check passed" in result.summary()

    # A slab whose bottom moved is left out.
    moved = relaxed["GaAs_111_term_1"].copy()
    index = study.face_sets["GaAs_111"][1].bottom_indices[0]
    moved.translate_sites([index], [0.3, 0.0, 0.0], frac_coords=False, to_unit_cell=False)
    relaxed["GaAs_111_term_1"] = moved
    result = study.analyse(energies, relaxed, points=5)
    assert "GaAs_111_term_1" not in result.diagram.curves
    assert result.bottom_checks["GaAs_111"].failed == ["term_1"]


def test_cluster_analysis_uses_the_eq_9_fit():
    study = PolarSurfaceStudy(
        gaas(), faces=[(1, 1, 1)], references={"Ga": ga_metal(), "As": as_solid()}, bilayers=2,
        pseudo_hydrogen_method="clusters", cluster_sizes=(2, 3, 4, 5), eq7_check=False,
    )
    energies = fake_energies(study)
    references = study.binary_references(energies)
    mu = references.chemical_potentials_ev(0.0)
    face, edge, corner = -1.3, -1.2, -1.0
    for size in (2, 3, 4, 5):
        c = polar.cluster_counts(size)
        energies[f"cluster_As_n{size}"] = (c["outer"] * mu["As"] + c["inner"] * (references.bulk_energy_per_formula_unit_ev - mu["As"])
                                           + c["face"] * face + c["edge"] * edge + c["corner"] * corner)
    hydrogen = study.pseudo_hydrogen_references(energies, references)
    assert hydrogen["As"].method == "cluster"
    assert hydrogen["As"].value_ev(mu["As"]) == pytest.approx(face)


def test_read_vasp_results_from_graph_outputs():
    class Graph:
        outputs = {
            "a_misc": {"total_energies": {"energy_extrapolated": -12.5}},
            "a_structure": "relaxed-a",
            "b_misc": {"energy_no_entropy": -3.0},
            "b_structure": "relaxed-b",
        }

    energies, structures = read_vasp_results(Graph(), ["a", "b"])
    assert energies == {"a": -12.5, "b": -3.0} and structures == {"a": "relaxed-a", "b": "relaxed-b"}


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vasp_workgraph_carries_pseudo_hydrogen_kinds_and_settings(tmp_path):
    pytest.importorskip("aiida_vasp")
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    study = gaas_study(faces=[(1, 1, 1)], bilayers=2, eq7_check=False)
    config = psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=f"vasp@{computer.label}", incar={"ENCUT": 400, "IBRION": 2, "NSW": 100, "EDIFFG": -0.005},
            potential_family="PBE", potential_mapping=study.potential_mapping({"Ga": "Ga_d"}),
            kpoints_spacing=0.2,
        ),
        execution=psteros.ExecutionPolicy(computer=computer.label, queue="debug", max_concurrent_jobs=1,
                                          resources={"num_machines": 1, "num_mpiprocs_per_machine": 1}),
        role_overrides=study.vasp_overrides(),
    )
    graph = psteros.build_surface_workgraph(study.structures, config, submit=False)
    slab = graph.tasks["GaAs_111_term_0_vasp"].inputs
    assert sorted(slab.structure.value.get_kind_names()) == ["As", "Ga", "H0p75"]
    assert slab.parameters.value.get_dict()["LDIPOL"] is True
    assert slab.parameters.value.get_dict()["ENCUT"] == 400
    assert slab.potential_mapping.value.get_dict()["H0p75"] == "H.75"
    molecule = graph.tasks["pseudo_molecule_As_vasp"].inputs
    assert molecule.kpoints_spacing.value.value == pytest.approx(10.0)
    assert graph.tasks["bulk_vasp"].inputs.kpoints_spacing.value.value == pytest.approx(0.2)
