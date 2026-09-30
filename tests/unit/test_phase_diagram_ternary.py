"""Tier-1 tests for ternary-oxide surface phase diagrams over (Delta mu_A, Delta mu_O).

The SrTiO3 numbers mimic the textbook case: Delta H_f = -17 eV per formula
unit, SrO at -6.1 eV and TiO2 at -9.8 eV, which leave the classic narrow strip
-7.2 <= Delta mu_Sr + Delta mu_O <= -6.1 eV.
"""

from __future__ import annotations

import csv

import pytest

import psteros
from psteros.phase_diagram_ternary import _area

A = 15.25  # SrTiO3(001) 1x1 face, A^2
ELEMENTS = {"Sr": -1.0, "Ti": -2.0}
E_O2 = -8.0
SRO = psteros.CompetingPhase("SrO", -11.1, {"Sr": 1, "O": 1})  # Delta H = -6.1 eV
TIO2 = psteros.CompetingPhase("TiO2", -19.8, {"Ti": 1, "O": 2})  # Delta H = -9.8 eV


def references(**overrides) -> psteros.TernaryOxideReferences:
    values = dict(
        bulk_energy_ev=-64.0,  # two formula units: Delta H_f = -32 + 1 + 2 + 12 = -17 eV each
        bulk_composition={"Sr": 2, "Ti": 2, "O": 6},
        oxygen_molecule_energy_ev=E_O2,
        element_energies_per_atom_ev=ELEMENTS,
        competing_phases=(SRO, TIO2),
    )
    values.update(overrides)
    return psteros.TernaryOxideReferences(**values)


def terminations(refs=None) -> list[psteros.SlabTermination]:
    """SrO- and TiO2-terminated 7-layer slabs plus an O-deficient SrO termination.

    Energies are set from gamma at Delta mu = 0 (eV/A^2), so the SrO/TiO2 boundary
    lies at Delta mu_Sr + Delta mu_O = -6.6 eV.
    """
    refs = refs or references()
    mu = refs.chemical_potentials_ev(0.0, 0.0)
    rows = []
    for label, counts, gamma0 in (
        ("SrO", {"Sr": 4, "Ti": 3, "O": 10}, -0.14139),
        ("TiO2", {"Sr": 3, "Ti": 4, "O": 11}, 0.29139),
        ("SrO-VO", {"Sr": 4, "Ti": 3, "O": 8}, 0.19),
    ):
        energy = 2 * A * gamma0 + sum(n * mu[element] for element, n in counts.items())
        rows.append(psteros.SlabTermination(label, energy, counts, A))
    return rows


def same_polygon(actual, expected, tolerance=1e-9) -> bool:
    return len(actual) == len(expected) and all(
        any(abs(p[0] - q[0]) < tolerance and abs(p[1] - q[1]) < tolerance for q in actual) for p in expected
    )


def test_references_reduce_the_cell_and_default_to_the_first_cation() -> None:
    refs = references()
    assert (refs.independent, refs.eliminated) == ("Sr", "Ti")
    assert refs.formula == "SrTiO3" and refs.formula_unit == (1, 1, 3)
    assert refs.bulk_energy_per_formula_unit_ev == pytest.approx(-32.0)
    assert refs.formation_enthalpy_ev == pytest.approx(-17.0)
    assert refs.delta_mu_eliminated_ev(-3.0, -1.0) == pytest.approx(-17.0 + 3.0 + 3.0)
    swapped = references(independent="Ti")
    assert (swapped.independent, swapped.eliminated, swapped.formula) == ("Ti", "Sr", "TiSrO3")


def test_elements_alone_bound_a_triangle() -> None:
    refs = references(competing_phases=())
    assert same_polygon(refs.stability_region, [(0.0, 0.0), (-17.0, 0.0), (0.0, -17.0 / 3)])
    assert set(refs.stability_boundaries) == {"Sr", "O2", "Ti"}


def test_competing_phases_cut_the_region_to_the_srtio3_strip() -> None:
    refs = references()
    assert same_polygon(
        refs.stability_region, [(-6.1, 0.0), (-7.2, 0.0), (-2.3, -4.9), (-0.65, -5.45)]
    )
    assert set(refs.stability_boundaries) == {"O2", "SrO", "TiO2", "Ti"}
    assert refs.in_stability_region(-6.5, -0.1)
    assert not refs.in_stability_region(-3.0, -1.0)  # SrO would form
    assert not refs.in_stability_region(-6.0, -1.5)  # TiO2 would form


def test_unstable_oxide_and_invalid_inputs_are_rejected() -> None:
    too_stable_tio2 = psteros.CompetingPhase("TiO2", -21.5, {"Ti": 1, "O": 2})  # Delta H = -11.5 eV
    with pytest.raises(ValueError, match="no stability region.*TiO2"):
        references(competing_phases=(SRO, too_stable_tio2))
    with pytest.raises(ValueError, match="ternary oxide"):
        references(bulk_composition={"Sn": 1, "O": 2})
    with pytest.raises(ValueError, match="lacks"):
        references(element_energies_per_atom_ev={"Sr": -1.0})
    with pytest.raises(ValueError, match="absent from the bulk"):
        references(competing_phases=(psteros.CompetingPhase("BaO", -10.0, {"Ba": 1, "O": 1}),))
    with pytest.raises(ValueError, match="not negative"):
        references(bulk_energy_ev=-20.0)
    with pytest.raises(ValueError, match="independent"):
        references(independent="O")
    refs = references()
    with pytest.raises(ValueError, match="expected only"):
        psteros.ternary_surface_phase_diagram(
            [psteros.SlabTermination("Ba", -1.0, {"Ba": 1, "O": 1}, A)], refs
        )


def test_planes_match_the_legacy_surface_excess_formulation() -> None:
    refs = references()
    diagram = psteros.ternary_surface_phase_diagram(terminations(refs), refs, points=5)
    x, y, z = refs.formula_unit
    for termination in terminations(refs):
        n_sr, n_ti, n_o = (termination.composition[e] for e in ("Sr", "Ti", "O"))
        # psteros.core.thermodynamics.calculate_surface_energy_ternary (primary form)
        phi = (
            termination.slab_energy_ev - n_sr * ELEMENTS["Sr"] - n_ti * ELEMENTS["Ti"] - n_o * E_O2 / 2
            - (n_ti / y) * refs.formation_enthalpy_ev
        ) / (2 * A)
        gamma_sr = (n_sr - (x / y) * n_ti) / (2 * A)
        gamma_o = (n_o - (z / y) * n_ti) / (2 * A)
        assert diagram.planes[termination.label] == pytest.approx((phi, -gamma_sr, -gamma_o))


def test_gamma_does_not_depend_on_the_axis_element() -> None:
    slabs = terminations()  # same slab energies; only the axes change
    by_sr = psteros.ternary_surface_phase_diagram(slabs, references())
    by_ti = psteros.ternary_surface_phase_diagram(slabs, references(independent="Ti"))
    delta_mu_sr, delta_mu_o = -5.0, -1.5
    delta_mu_ti = references().delta_mu_eliminated_ev(delta_mu_sr, delta_mu_o)
    for label in by_sr.planes:
        assert by_sr.gamma_ev_per_angstrom2(label, delta_mu_sr, delta_mu_o) == pytest.approx(
            by_ti.gamma_ev_per_angstrom2(label, delta_mu_ti, delta_mu_o)
        )


def test_regions_tile_the_stability_region_with_exact_boundaries() -> None:
    refs = references()
    diagram = psteros.ternary_surface_phase_diagram(terminations(refs), refs)
    assert sum(_area(polygon) for polygon in diagram.regions.values()) == pytest.approx(
        _area(refs.stability_region)
    )
    shared = [p for p in diagram.regions["SrO"] if any(abs(p[0] - q[0]) + abs(p[1] - q[1]) < 1e-9 for q in diagram.regions["TiO2"])]
    boundary = A * (-0.14139 - 0.29139)  # gamma_SrO = gamma_TiO2, about -6.6 eV
    assert shared and all(p[0] + p[1] == pytest.approx(boundary) for p in shared)
    # The O-deficient termination wins only in the O-poor corner, below Delta mu_O = -5.054 eV.
    assert diagram.regions["SrO-VO"]
    assert max(p[1] for p in diagram.regions["SrO-VO"]) == pytest.approx(A * (-0.14139 - 0.19))
    for label, polygon in diagram.regions.items():
        cx = sum(p[0] for p in polygon) / len(polygon)
        cy = sum(p[1] for p in polygon) / len(polygon)
        assert diagram.stable_termination(cx, cy) == label


def test_csv_holds_every_grid_point(tmp_path) -> None:
    refs = references()
    diagram = psteros.ternary_surface_phase_diagram(terminations(refs), refs, points=6)
    with diagram.to_csv(tmp_path / "srtio3.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "delta_mu_Sr_eV", "delta_mu_O_eV", "delta_mu_Ti_eV",
        "gamma_SrO_Jm2", "gamma_TiO2_Jm2", "gamma_SrO-VO_Jm2",
        "stable_termination", "in_stability_region",
    ]
    assert len(rows) == 36
    for row in rows:
        sr, o = float(row["delta_mu_Sr_eV"]), float(row["delta_mu_O_eV"])
        assert float(row["delta_mu_Ti_eV"]) == pytest.approx(refs.delta_mu_eliminated_ev(sr, o))
        assert row["in_stability_region"] == str(refs.in_stability_region(sr, o))
        assert row["stable_termination"] == diagram.stable_termination(sr, o)
    assert {row["in_stability_region"] for row in rows} == {"True", "False"}
    header = diagram.to_csv(tmp_path / "ev.csv", units="eV/A2").read_text().splitlines()[0]
    assert "gamma_SrO_eVA2" in header


def test_plot_writes_a_figure_file(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    refs = references()
    diagram = psteros.ternary_surface_phase_diagram(terminations(refs), refs)
    for suffix in ("png", "pdf"):
        assert diagram.plot(tmp_path / f"srtio3.{suffix}", title="SrTiO3(001)").stat().st_size > 1000
