"""Tier-1 tests for surface phase diagrams of any binary compound."""

from __future__ import annotations

import csv

import pytest

import psteros

A = 20.0  # surface area of one face, A^2


def oxide_pair():
    """The same SnO2 references as an oxide and as a general binary compound."""

    oxide = psteros.BinaryOxideReferences(
        bulk_energy_ev=-20.0, bulk_composition={"Sn": 2, "O": 4},
        oxygen_molecule_energy_ev=-5.0, metal_energy_per_atom_ev=-1.0,
    )
    general = psteros.BinaryReferences(
        bulk_energy_ev=-20.0, bulk_composition={"Sn": 2, "O": 4},
        reference_energies_per_atom_ev={"O": -2.5, "Sn": -1.0},
    )
    return oxide, general


def sno2_terminations():
    return [
        psteros.SlabTermination("o", -122.0, {"Sn": 6, "O": 12}, A),
        psteros.SlabTermination("sno", -112.0, {"Sn": 6, "O": 10}, A),
        psteros.SlabTermination("sn2o", -103.0, {"Sn": 6, "O": 8}, A),
    ]


def gaas_references(**overrides) -> psteros.BinaryReferences:
    # E(GaAs) = -8.5 per formula unit, E(Ga) = -3.0, E(As) = -4.7 -> Delta H_f = -0.8 eV
    values = dict(
        bulk_energy_ev=-17.0,
        bulk_composition={"Ga": 2, "As": 2},
        reference_energies_per_atom_ev={"Ga": -3.0, "As": -4.7},
        reservoir_labels={"Ga": "Ga bulk", "As": "As bulk"},
    )
    values.update(overrides)
    return psteros.BinaryReferences(**values)


def test_general_references_reproduce_the_oxide_diagram(tmp_path) -> None:
    oxide, general = oxide_pair()
    assert general.variable == "O" and general.other == "Sn"
    assert general.formula_unit == oxide.formula_unit == (1, 2)
    assert general.formation_enthalpy_ev == pytest.approx(oxide.formation_enthalpy_ev)
    assert general.poor_limit_ev == pytest.approx(oxide.oxygen_poor_limit_ev)

    first = psteros.surface_phase_diagram(sno2_terminations(), oxide, points=11)
    second = psteros.surface_phase_diagram(sno2_terminations(), general, points=11)
    assert first.delta_mu_ev == pytest.approx(second.delta_mu_ev)
    for label in first.curves:
        assert [p.gamma_ev_per_angstrom2 for p in first.curves[label]] == pytest.approx(
            [p.gamma_ev_per_angstrom2 for p in second.curves[label]]
        )
    assert first.stable == second.stable
    assert [t[0] for t in first.transitions] == pytest.approx([t[0] for t in second.transitions])
    header_oxide = first.to_csv(tmp_path / "oxide.csv").read_text().splitlines()[0]
    header_general = second.to_csv(tmp_path / "general.csv").read_text().splitlines()[0]
    assert header_oxide == header_general
    assert header_oxide.startswith("delta_mu_O_eV,")


def test_oxide_aliases_stay_available() -> None:
    oxide, _ = oxide_pair()
    diagram = psteros.surface_phase_diagram(sno2_terminations(), oxide, points=5)
    assert diagram.delta_mu_oxygen_ev == diagram.delta_mu_ev
    assert diagram.curves["o"][0].delta_mu_oxygen_ev == diagram.curves["o"][0].delta_mu_ev


def test_axis_element_defaults_to_the_most_electronegative() -> None:
    references = gaas_references()
    assert references.variable == "As" and references.other == "Ga"
    assert references.formula == "GaAs"
    assert references.formation_enthalpy_ev == pytest.approx(-0.8)
    assert references.poor_limit_ev == pytest.approx(-0.8)
    zno = psteros.BinaryReferences(-9.0, {"Zn": 1, "O": 1}, {"O": -4.9, "Zn": -1.3})
    assert zno.variable == "O"
    gan = psteros.BinaryReferences(-24.0, {"Ga": 2, "N": 2}, {"N": -8.3, "Ga": -3.0}, variable="Ga")
    assert gan.variable == "Ga" and gan.other == "N"


def test_chemical_potentials_satisfy_bulk_equilibrium() -> None:
    references = gaas_references()
    for delta_mu in (-0.8, -0.3, 0.0):
        mu = references.chemical_potentials_ev(delta_mu)
        assert mu["As"] == pytest.approx(-4.7 + delta_mu)
        assert mu["Ga"] + mu["As"] == pytest.approx(-8.5)
    al2s3 = psteros.BinaryReferences(-50.0, {"Al": 4, "S": 6}, {"Al": -3.7, "S": -4.1})
    mu = al2s3.chemical_potentials_ev(-0.2)
    assert 2 * mu["Al"] + 3 * mu["S"] == pytest.approx(-25.0)


def test_gaas_phase_diagram_slopes_follow_the_surface_excess(tmp_path) -> None:
    references = gaas_references()
    terminations = [
        psteros.SlabTermination("stoich", -8.5 * 6 + 2 * A * 0.05, {"Ga": 6, "As": 6}, A),
        psteros.SlabTermination("as_rich", -8.5 * 6 - 4.7 * 2 + 2 * A * 0.06, {"Ga": 6, "As": 8}, A),
        psteros.SlabTermination("ga_rich", -8.5 * 6 - 3.0 * 2 + 2 * A * 0.06, {"Ga": 8, "As": 6}, A),
    ]
    diagram = psteros.surface_phase_diagram(terminations, references, points=9)
    assert diagram.delta_mu_ev[0] == pytest.approx(-0.8)
    assert diagram.delta_mu_ev[-1] == pytest.approx(0.0)
    flat = [p.gamma_ev_per_angstrom2 for p in diagram.curves["stoich"]]
    assert max(flat) - min(flat) == pytest.approx(0.0, abs=1e-12)
    rich = diagram.curves["as_rich"]
    assert rich[-1].gamma_ev_per_angstrom2 == pytest.approx(0.06)  # As reference at Delta mu = 0
    assert rich[0].gamma_ev_per_angstrom2 > rich[-1].gamma_ev_per_angstrom2  # As-rich costs more As-poor
    assert diagram.stable[-1] == "stoich"

    path = diagram.to_csv(tmp_path / "gaas.csv", units="eV/A2")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0])[0] == "delta_mu_As_eV"
    assert {row["in_stability_window"] for row in rows} == {"True"}


def test_figure_uses_element_labels(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    references = gaas_references()
    terminations = [
        psteros.SlabTermination("stoich", -50.0, {"Ga": 6, "As": 6}, A),
        psteros.SlabTermination("as_rich", -60.0, {"Ga": 6, "As": 8}, A),
    ]
    diagram = psteros.surface_phase_diagram(terminations, references)
    figure = diagram.figure()
    texts = [text.get_text() for axis in figure.axes for text in axis.texts]
    assert any("As-poor limit" in text and "Ga bulk" in text for text in texts)
    assert any("As-rich limit" in text and "As bulk" in text for text in texts)
    assert "Delta\\mu_\\mathrm{As}" in figure.axes[1].get_xlabel()
    assert diagram.plot(tmp_path / "gaas.png").stat().st_size > 1000


def test_invalid_general_references_are_rejected() -> None:
    with pytest.raises(ValueError, match="binary"):
        psteros.BinaryReferences(-10.0, {"Ag": 3, "P": 1, "O": 4}, {"O": -4.9})
    with pytest.raises(ValueError, match="axis element"):
        psteros.BinaryReferences(-17.0, {"Ga": 2, "As": 2}, {"Ga": -3.0})
    with pytest.raises(ValueError, match="not in the compound"):
        psteros.BinaryReferences(-17.0, {"Ga": 2, "As": 2}, {"As": -4.7, "N": -8.0})
    with pytest.raises(ValueError, match="unstable"):
        gaas_references(reference_energies_per_atom_ev={"Ga": -4.0, "As": -5.0})
    with pytest.raises(ValueError, match="variable element"):
        gaas_references(variable="O")


def test_terminations_must_match_the_compound() -> None:
    references = gaas_references()
    wrong = [psteros.SlabTermination("zno", -50.0, {"Zn": 6, "O": 6}, A)]
    with pytest.raises(ValueError, match="expected only Ga and As"):
        psteros.surface_phase_diagram(wrong, references)
    no_ga = psteros.BinaryReferences(-17.0, {"Ga": 2, "As": 2}, {"As": -4.7})
    with pytest.raises(ValueError, match="no reference energy of Ga"):
        psteros.surface_phase_diagram(
            [psteros.SlabTermination("s", -50.0, {"Ga": 6, "As": 6}, A)], no_ga
        )
    diagram = psteros.surface_phase_diagram(
        [psteros.SlabTermination("s", -50.0, {"Ga": 6, "As": 6}, A)], no_ga, delta_mu_range=(-1.0, 0.0)
    )
    assert diagram.in_stability_window(-1.0)


def test_binary_equilibrium_helper_matches_the_oxide_helper() -> None:
    oxide = psteros.surface_energy_oxide_equilibrium(
        slab_energy_ev=-112.0, n_metal=6, n_oxygen=10, bulk_formula_energy_ev=-10.0,
        oxygen_reference_energy_ev=-5.0, delta_mu_oxygen_ev=-0.7, surface_area_angstrom2=A,
    )
    general = psteros.surface_energy_binary_equilibrium(
        slab_energy_ev=-112.0, n_other=6, n_variable=10, bulk_formula_energy_ev=-10.0,
        variable_reference_energy_ev=-2.5, delta_mu_ev=-0.7, surface_area_angstrom2=A,
        formula_unit=(1, 2),
    )
    assert general == oxide
    corrected = psteros.surface_energy_binary_equilibrium(
        slab_energy_ev=-112.0, n_other=6, n_variable=10, bulk_formula_energy_ev=-10.0,
        variable_reference_energy_ev=-2.5, delta_mu_ev=-0.7, surface_area_angstrom2=A,
        formula_unit=(1, 2), surfaces=1, reservoir_correction_ev=-3.0,
    )
    assert corrected.gamma_ev_per_angstrom2 == pytest.approx(2 * oxide.gamma_ev_per_angstrom2 + 3.0 / A)
