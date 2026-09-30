"""Tier-1 tests for binary-oxide surface phase diagrams (plot and CSV export)."""

from __future__ import annotations

import csv

import pytest

import psteros

A = 20.0  # surface area of one face, A^2


def sno2_references(**overrides) -> psteros.BinaryOxideReferences:
    values = dict(
        bulk_energy_ev=-20.0,  # Sn2O4 cell -> -10 eV per SnO2
        bulk_composition={"Sn": 2, "O": 4},
        oxygen_molecule_energy_ev=-10.0,
        metal_energy_per_atom_ev=-3.0,  # Delta H_f = -10 - (-3) - (-10) = +3 eV: unstable oxide
    )
    values.update(overrides)
    return psteros.BinaryOxideReferences(**values)


def stable_sno2_references() -> psteros.BinaryOxideReferences:
    # E_fu = -10, E_Sn = -1, E_O2 = -5  ->  Delta H_f = -10 + 1 + 5 = -4 eV, O-poor limit -2 eV
    return psteros.BinaryOxideReferences(
        bulk_energy_ev=-20.0,
        bulk_composition={"Sn": 2, "O": 4},
        oxygen_molecule_energy_ev=-5.0,
        metal_energy_per_atom_ev=-1.0,
    )


def sno2_terminations(gamma_sno: float = 0.10) -> list[psteros.SlabTermination]:
    # Energies chosen so gamma(Delta mu_O = 0) = 0.05, gamma_sno, 0.12 eV/A^2 for O excess 0, -2, -4.
    references = stable_sno2_references()
    mu_o = references.oxygen_molecule_energy_ev / 2
    rows = []
    for label, n_o, gamma0 in (("o", 12, 0.05), ("sno", 10, gamma_sno), ("sn2o", 8, 0.12)):
        energy = 2 * A * gamma0 + 6 * references.bulk_energy_per_formula_unit_ev + (n_o - 12) * mu_o
        rows.append(psteros.SlabTermination(label, energy, {"Sn": 6, "O": n_o}, A))
    return rows


def test_references_reduce_the_bulk_cell_to_a_formula_unit() -> None:
    references = stable_sno2_references()
    assert references.metal == "Sn"
    assert references.formula == "SnO2"
    assert references.formula_unit == (1, 2)
    assert references.bulk_energy_per_formula_unit_ev == pytest.approx(-10.0)
    assert references.formation_enthalpy_ev == pytest.approx(-4.0)
    assert references.oxygen_poor_limit_ev == pytest.approx(-2.0)
    corundum = psteros.BinaryOxideReferences(-60.0, {"Al": 4, "O": 6}, -5.0, -2.0)
    assert (corundum.formula, corundum.formula_unit) == ("Al2O3", (2, 3))
    assert corundum.bulk_energy_per_formula_unit_ev == pytest.approx(-30.0)
    # Delta H_f = -30 - 2(-2) - 3(-5)/2 = -18.5 eV; O-poor limit = Delta H_f / 3
    assert corundum.oxygen_poor_limit_ev == pytest.approx(-18.5 / 3)


def test_references_reject_unstable_or_non_binary_oxides() -> None:
    with pytest.raises(ValueError, match="not negative"):
        sno2_references()
    with pytest.raises(ValueError, match="binary oxide"):
        psteros.BinaryOxideReferences(-10.0, {"Ag": 3, "P": 1, "O": 4}, -5.0)
    with pytest.raises(ValueError, match="binary oxide"):
        psteros.BinaryOxideReferences(-10.0, {"Sn": 1}, -5.0)


def test_mo2_phase_diagram_matches_the_oxide_equilibrium_helper() -> None:
    references = stable_sno2_references()
    diagram = psteros.surface_phase_diagram(sno2_terminations(), references, points=11)
    assert diagram.delta_mu_oxygen_ev[0] == pytest.approx(-2.0)
    assert diagram.delta_mu_oxygen_ev[-1] == 0.0
    for termination in sno2_terminations():
        for point in diagram.curves[termination.label]:
            expected = psteros.surface_energy_oxide_equilibrium(
                slab_energy_ev=termination.slab_energy_ev,
                n_metal=6,
                n_oxygen=termination.composition["O"],
                bulk_formula_energy_ev=-10.0,
                oxygen_reference_energy_ev=-5.0,
                delta_mu_oxygen_ev=point.delta_mu_oxygen_ev,
                surface_area_angstrom2=A,
            )
            assert point.gamma_ev_per_angstrom2 == pytest.approx(expected.gamma_ev_per_angstrom2)


def test_transitions_are_exact_crossings_of_the_lower_envelope() -> None:
    # gamma_o = 0.05; gamma_sno = 0.10 + dmu/20; gamma_sn2o = 0.12 + dmu/10 (eV/A^2).
    # sn2o crosses o at -0.7 eV, where sno (0.065) is above both: sno is never stable.
    diagram = psteros.surface_phase_diagram(sno2_terminations(), stable_sno2_references(), points=7)
    assert len(diagram.transitions) == 1
    delta_mu, below, above = diagram.transitions[0]
    assert (below, above) == ("sn2o", "o")
    assert delta_mu == pytest.approx(-0.7)
    assert diagram.stable[0] == "sn2o" and diagram.stable[-1] == "o"
    assert "sno" not in diagram.stable

    # With gamma_sno(0) = 0.08, sno is lowest between -0.8 eV (vs sn2o) and -0.6 eV (vs o).
    diagram = psteros.surface_phase_diagram(sno2_terminations(0.08), stable_sno2_references())
    assert [(round(x, 12), b, a) for x, b, a in diagram.transitions] == [
        (-0.8, "sn2o", "sno"), (-0.6, "sno", "o"),
    ]
    assert "sno" in diagram.stable


def test_general_oxides_keep_stoichiometric_slabs_independent_of_oxygen() -> None:
    # Regression: the MO2-only formula gave a Delta mu_O dependent gamma for ZnO.
    zno = psteros.BinaryOxideReferences(-16.0, {"Zn": 2, "O": 2}, -9.0, -1.0)
    slab = psteros.SlabTermination("zn8o8", -60.0, {"Zn": 8, "O": 8}, A)
    diagram = psteros.surface_phase_diagram([slab], zno, points=5)
    gammas = {round(p.gamma_ev_per_angstrom2, 12) for p in diagram.curves["zn8o8"]}
    assert gammas == {round((-60.0 - 8 * -8.0) / (2 * A), 12)}

    # M2O3: compare with explicit elemental chemical potentials.
    alumina = psteros.BinaryOxideReferences(-60.0, {"Al": 4, "O": 6}, -5.0, -2.0)
    slab = psteros.SlabTermination("al4o7", -70.0, {"Al": 4, "O": 7}, A)
    point = psteros.surface_phase_diagram([slab], alumina, delta_mu_range=(-1.0, 0.0), points=2).curves["al4o7"][0]
    mu_o = -2.5 - 1.0
    mu_al = (-30.0 - 3 * mu_o) / 2
    expected = psteros.surface_energy_elemental(
        slab_energy_ev=-70.0,
        stoichiometry={"Al": 4, "O": 7},
        chemical_potentials_ev={"Al": mu_al, "O": mu_o},
        surface_area_angstrom2=A,
    )
    assert point.gamma_ev_per_angstrom2 == pytest.approx(expected)


def test_oxide_equilibrium_helper_accepts_a_formula_unit() -> None:
    kwargs = dict(
        slab_energy_ev=-60.0, n_metal=8, n_oxygen=8, bulk_formula_energy_ev=-8.0,
        oxygen_reference_energy_ev=-9.0, surface_area_angstrom2=A, formula_unit=(1, 1),
    )
    low = psteros.surface_energy_oxide_equilibrium(delta_mu_oxygen_ev=-2.0, **kwargs)
    high = psteros.surface_energy_oxide_equilibrium(delta_mu_oxygen_ev=0.0, **kwargs)
    assert low.gamma_ev_per_angstrom2 == pytest.approx(high.gamma_ev_per_angstrom2)
    with pytest.raises(ValueError, match="formula_unit"):
        psteros.surface_energy_oxide_equilibrium(delta_mu_oxygen_ev=0.0, **{**kwargs, "formula_unit": (0, 1)})


def test_range_defaults_to_the_window_and_limits_join_a_wider_grid() -> None:
    no_metal = psteros.BinaryOxideReferences(-20.0, {"Sn": 2, "O": 4}, -5.0)
    with pytest.raises(ValueError, match="delta_mu_range is required"):
        psteros.surface_phase_diagram(sno2_terminations(), no_metal)
    diagram = psteros.surface_phase_diagram(
        sno2_terminations(), stable_sno2_references(), delta_mu_range=(-2.55, 0.45), points=4
    )
    grid = diagram.delta_mu_oxygen_ev
    assert -2.0 in grid and 0.0 in grid and len(grid) == 6
    flags = [diagram.in_stability_window(x) for x in grid]
    assert flags == [False, True, True, True, True, False]


def test_csv_export_holds_every_curve_and_the_stable_termination(tmp_path) -> None:
    diagram = psteros.surface_phase_diagram(sno2_terminations(), stable_sno2_references(), points=5)
    path = diagram.to_csv(tmp_path / "diagram.csv")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "delta_mu_O_eV", "gamma_o_Jm2", "gamma_sno_Jm2", "gamma_sn2o_Jm2",
        "stable_termination", "in_stability_window",
    ]
    assert len(rows) == 5
    assert float(rows[0]["delta_mu_O_eV"]) == pytest.approx(-2.0)
    assert float(rows[-1]["gamma_o_Jm2"]) == pytest.approx(0.05 * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2)
    assert [row["stable_termination"] for row in rows] == list(diagram.stable)
    assert {row["in_stability_window"] for row in rows} == {"True"}

    ev_path = diagram.to_csv(tmp_path / "diagram_ev.csv", units="eV/A2")
    header = ev_path.read_text().splitlines()[0]
    assert "gamma_o_eVA2" in header
    with pytest.raises(ValueError, match="units"):
        diagram.to_csv(tmp_path / "bad.csv", units="kcal")


def test_plot_writes_a_figure_file(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    diagram = psteros.surface_phase_diagram(sno2_terminations(), stable_sno2_references())
    for suffix in ("png", "pdf"):
        path = diagram.plot(tmp_path / f"diagram.{suffix}", title="SnO2(110)")
        assert path.stat().st_size > 1000
    # More terminations than categorical colours still render (dashed second cycle).
    references = stable_sno2_references()
    many = [
        psteros.SlabTermination(f"t{i}", -125.0 + i, {"Sn": 6, "O": 12 - (i % 3)}, A) for i in range(10)
    ]
    figure = psteros.surface_phase_diagram(many, references).figure(units="eV/A2")
    assert len(figure.axes) == 2


def test_termination_from_structure_reads_composition_and_area() -> None:
    slab, _ = psteros.sno2_110_slab(termination="sno", triple_layers=3, vacuum_angstrom=10.0)
    termination = psteros.SlabTermination.from_structure("sno", -100.0, slab)
    assert termination.composition == {"Sn": 6, "O": 10}
    assert termination.surface_area_angstrom2 == pytest.approx(slab.surface_area)


def test_invalid_termination_sets_are_rejected() -> None:
    references = stable_sno2_references()
    slab = sno2_terminations()[0]
    with pytest.raises(ValueError, match="unique"):
        psteros.surface_phase_diagram([slab, slab], references)
    with pytest.raises(ValueError, match="expected only Sn and O"):
        psteros.surface_phase_diagram(
            [psteros.SlabTermination("oh", -1.0, {"Sn": 6, "O": 12, "H": 2}, A)], references
        )
    with pytest.raises(ValueError, match="integer"):
        psteros.SlabTermination("x", -1.0, {"Sn": 6.5, "O": 12}, A)
    with pytest.raises(ValueError, match="increasing"):
        psteros.surface_phase_diagram([slab], references, delta_mu_range=(0.0, -1.0))
