"""Tier-1 cross-check of the campaign analysis against the phase-diagram formulas.

The energies are typed by hand (0 K values of realistic size, not taken from a
calculation) and no AiiDA profile is needed.  Every expected number is written
out from the formulas of ``docs/source/theory.rst`` and
``docs/source/phase-diagram.rst``, never taken from the function under test.
Energies are in eV, areas in A^2 (one face), chemical potentials in eV per atom.
"""

from __future__ import annotations

import pytest

import psteros
from psteros import (
    BinaryOxideReferences,
    CampaignEntry,
    CompetingPhase,
    SlabTermination,
    TernaryOxideReferences,
    campaign_chemical_potentials,
    campaign_references,
    campaign_surface_energies,
    surface_phase_diagram,
    ternary_surface_phase_diagram,
)

# 1 eV = 1.602176634e-19 J and 1 A^2 = 1e-20 m^2
ONE_EV_PER_A2_IN_J_PER_M2 = 1.602_176_634e-19 / 1e-20


def reference(label, composition, energy_ev, phase="solid", state="finished"):
    return CampaignEntry(label, "references", phase, composition, energy_ev, state=state)


def slab(label, composition, energy_ev, area_angstrom2):
    return CampaignEntry(label, "slabs", "solid", composition, energy_ev, surface_area_angstrom2=area_angstrom2)


def point_at(curve, delta_mu_oxygen_ev):
    point = min(curve, key=lambda p: abs(p.delta_mu_oxygen_ev - delta_mu_oxygen_ev))
    assert point.delta_mu_oxygen_ev == pytest.approx(delta_mu_oxygen_ev, abs=1e-9)
    return point


def same_polygon(actual, expected, tolerance=1e-9) -> bool:
    return len(actual) == len(expected) and all(
        any(abs(p[0] - q[0]) < tolerance and abs(p[1] - q[1]) < tolerance for q in actual) for p in expected
    )


def polygon_area(points) -> float:
    """Shoelace formula, A^2 in the (dmu_A, dmu_O) plane."""
    return abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1]))) / 2


# --- 1a. Unary Au: two slabs against the bulk per atom ----------------------------------

AU_BULK_ENERGY_EV = -13.1200  # 4-atom fcc cell, -3.2800 eV per atom
AU111_ENERGY_EV, AU111_AREA_ANGSTROM2 = -48.9627, 28.8300  # 2x2 (111) slab, 16 atoms, one face
AU100_ENERGY_EV, AU100_AREA_ANGSTROM2 = -34.5658, 33.2928  # 2x2 (100) slab, 12 atoms, one face


def gold_entries():
    return (
        reference("au_bulk", {"Au": 4}, AU_BULK_ENERGY_EV),
        slab("au111", {"Au": 16}, AU111_ENERGY_EV, AU111_AREA_ANGSTROM2),
        slab("au100", {"Au": 12}, AU100_ENERGY_EV, AU100_AREA_ANGSTROM2),
    )


def test_unary_gamma_matches_the_hand_formula_in_eV_and_J_per_m2() -> None:
    gammas = campaign_surface_energies(gold_entries())
    # gamma = (E_slab - N E_bulk / N_bulk) / (2 A), with E_bulk / N_bulk = -13.12 / 4 eV per atom
    gamma_111 = (-48.9627 - 16 * (-13.1200 / 4)) / (2 * 28.8300)
    gamma_100 = (-34.5658 - 12 * (-13.1200 / 4)) / (2 * 33.2928)
    assert gamma_111 == pytest.approx(0.0610, abs=1e-5)  # the value the energies were typed from
    assert gamma_100 == pytest.approx(0.0720, abs=1e-5)
    assert gammas["au111"] == pytest.approx(gamma_111, rel=1e-12)
    assert gammas["au100"] == pytest.approx(gamma_100, rel=1e-12)
    assert psteros.EV_PER_ANGSTROM2_TO_J_PER_M2 == pytest.approx(ONE_EV_PER_A2_IN_J_PER_M2, rel=1e-8)
    assert gammas["au111"] * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2 == pytest.approx(
        gamma_111 * ONE_EV_PER_A2_IN_J_PER_M2, rel=1e-12
    )


def test_unary_default_limit_is_the_elemental_one_and_surfaces_count_faces() -> None:
    entries = gold_entries()
    assert campaign_chemical_potentials(entries) == {"Au": pytest.approx(-3.28, rel=1e-12)}
    default = campaign_surface_energies(entries)
    explicit = campaign_surface_energies(entries, chemical_potentials_ev={"Au": -3.28})
    assert default == pytest.approx(explicit, rel=1e-12)
    one_face = campaign_surface_energies(entries, surfaces=1)
    assert one_face["au111"] == pytest.approx(2 * default["au111"], rel=1e-12)
    assert one_face["au100"] == pytest.approx(2 * default["au100"], rel=1e-12)


def test_reservoirs_and_chemical_potentials_are_not_both_given() -> None:
    with pytest.raises(ValueError, match="not both"):
        campaign_surface_energies(gold_entries(), chemical_potentials_ev={"Au": -3.28}, reservoirs={"Au": "au_bulk"})


def test_an_unfinished_unrelated_reference_does_not_block_a_unary_surface_energy() -> None:
    # Only the elements of the slabs are needed; a running Pt reference is not one of them.
    running_pt = CampaignEntry("pt_bulk", "references", "solid", {"Pt": 4}, None, state="running")
    entries = gold_entries() + (running_pt,)
    assert campaign_chemical_potentials(entries, elements={"Au"}) == {"Au": pytest.approx(-3.28, rel=1e-12)}
    assert campaign_surface_energies(entries) == pytest.approx(campaign_surface_energies(gold_entries()), rel=1e-12)


def test_two_references_of_one_element_need_an_explicit_choice() -> None:
    entries = gold_entries() + (reference("au_alt", {"Au": 2}, -6.2000),)
    with pytest.raises(ValueError, match="several references of Au") as info:
        campaign_surface_energies(entries)
    assert "au_alt" in str(info.value) and "au_bulk" in str(info.value)
    chosen = campaign_surface_energies(entries, reservoirs={"Au": "au_alt"})
    assert chosen["au111"] == pytest.approx((-48.9627 - 16 * (-3.1000)) / (2 * 28.8300), rel=1e-12)


def test_missing_chemical_potential_names_the_slab_and_the_element() -> None:
    with pytest.raises(ValueError, match=r"au111: no chemical potential for \['Au'\]"):
        campaign_surface_energies(gold_entries(), chemical_potentials_ev={"Pd": -5.3})


def test_unfinished_slab_is_named_with_its_state() -> None:
    running = CampaignEntry(
        "au111", "slabs", "solid", {"Au": 16}, None, state="running", surface_area_angstrom2=28.83
    )
    with pytest.raises(ValueError, match="au111.*running"):
        campaign_surface_energies((gold_entries()[0], running))


# --- 1b. Binary SnO2: references by hand, three (110) terminations ----------------------

SNO2_AREA_ANGSTROM2 = 42.8  # one face
SNO2_E_O2_EV = -9.8600  # triplet O2
SNO2_E_SN_CELL_EV = -8.3000  # 2-atom Sn cell, -4.15 eV per atom
SNO2_E_BULK_EV = -38.4200  # Sn2O4 cell (2 SnO2): dH_f = -38.42/2 - (-4.15) - (-9.86) = -5.20 eV


def sno2_entries(extra=()):
    return (
        reference("sn", {"Sn": 2}, SNO2_E_SN_CELL_EV),
        reference("o2", {"O": 2}, SNO2_E_O2_EV, phase="gas"),
        reference("sno2", {"Sn": 2, "O": 4}, SNO2_E_BULK_EV),
        slab("o_term", {"Sn": 6, "O": 12}, -110.2952, SNO2_AREA_ANGSTROM2),
        slab("minus1_o", {"Sn": 6, "O": 10}, -98.1240, SNO2_AREA_ANGSTROM2),
        slab("minus2_o", {"Sn": 6, "O": 8}, -84.8400, SNO2_AREA_ANGSTROM2),
        *extra,
    )


def sno2_terminations():
    return [
        SlabTermination("o_term", -110.2952, {"Sn": 6, "O": 12}, SNO2_AREA_ANGSTROM2),
        SlabTermination("minus1_o", -98.1240, {"Sn": 6, "O": 10}, SNO2_AREA_ANGSTROM2),
        SlabTermination("minus2_o", -84.8400, {"Sn": 6, "O": 8}, SNO2_AREA_ANGSTROM2),
    ]


def hand_sno2_references() -> BinaryOxideReferences:
    return BinaryOxideReferences(
        bulk_energy_ev=-38.4200,
        bulk_composition={"Sn": 2, "O": 4},
        oxygen_molecule_energy_ev=-9.8600,
        metal_energy_per_atom_ev=-4.1500,
    )


def test_binary_references_equal_the_hand_built_object_field_by_field() -> None:
    built = campaign_references(sno2_entries(), host="sno2")
    hand = hand_sno2_references()
    assert isinstance(built, BinaryOxideReferences)
    assert built.bulk_energy_ev == pytest.approx(hand.bulk_energy_ev, rel=1e-15)
    assert dict(built.bulk_composition) == dict(hand.bulk_composition)
    assert built.oxygen_molecule_energy_ev == pytest.approx(hand.oxygen_molecule_energy_ev, rel=1e-15)
    assert built.metal_energy_per_atom_ev == pytest.approx(hand.metal_energy_per_atom_ev, rel=1e-15)
    assert (built.metal, built.formula, built.formula_unit) == (hand.metal, hand.formula, hand.formula_unit)
    assert built.formation_enthalpy_ev == pytest.approx(-5.20, rel=1e-12)
    assert built.oxygen_poor_limit_ev == pytest.approx(-2.60, rel=1e-12)
    assert built.oxygen_rich_limit_ev == 0.0
    assert built == hand


def test_binary_transitions_and_stable_terminations_agree_with_the_hand_built_object() -> None:
    window = (-2.60, 0.0)
    built = surface_phase_diagram(
        sno2_terminations(), campaign_references(sno2_entries(), host="sno2"), delta_mu_range=window, points=261
    )
    hand = surface_phase_diagram(sno2_terminations(), hand_sno2_references(), delta_mu_range=window, points=261)
    assert built.stable == hand.stable
    assert built.transitions == hand.transitions
    # minus2_o is lowest from -2.60 eV to -1.712 eV, minus1_o up to -1.1556 eV, o_term above
    assert built.stable[0] == "minus2_o" and built.stable[-1] == "o_term"
    assert [(round(d, 9), below, above) for d, below, above in built.transitions] == [
        (-1.712, "minus2_o", "minus1_o"),
        (-1.1556, "minus1_o", "o_term"),
    ]


def test_binary_gamma_at_one_point_matches_the_textbook_formula() -> None:
    diagram = surface_phase_diagram(
        sno2_terminations(), campaign_references(sno2_entries(), host="sno2"), delta_mu_range=(-2.60, 0.0), points=261
    )
    # gamma = [E_slab - n_Sn E_fu / x - (n_O - n_Sn y / x)(E(O2)/2 + dmu_O)] / (2 A), x = 1, y = 2
    e_slab, n_sn, n_o, area = -98.1240, 6, 10, 42.8
    e_formula_unit = -38.4200 / 2
    dmu_o = -1.0
    expected = (e_slab - n_sn * e_formula_unit - (n_o - n_sn * 2) * (-9.8600 / 2 + dmu_o)) / (2 * area)
    assert expected == pytest.approx(0.0616355, abs=1e-6)
    assert point_at(diagram.curves["minus1_o"], dmu_o).gamma_ev_per_angstrom2 == pytest.approx(expected, rel=1e-12)


def test_surface_energies_at_bulk_equilibrium_potentials_reproduce_the_diagram() -> None:
    # mu_O = E(O2)/2 + dmu_O and mu_Sn = (E_fu - 2 mu_O) in equilibrium with SnO2
    diagram = surface_phase_diagram(
        sno2_terminations(), campaign_references(sno2_entries(), host="sno2"), delta_mu_range=(-2.60, 0.0), points=261
    )
    for dmu_o, label in ((0.0, "o_term"), (-1.0, "minus1_o"), (-2.0, "minus2_o")):
        mu_o = -9.8600 / 2 + dmu_o
        mu_sn = (-19.2100) - 2 * mu_o
        gamma = campaign_surface_energies(sno2_entries(), chemical_potentials_ev={"Sn": mu_sn, "O": mu_o})[label]
        assert gamma == pytest.approx(point_at(diagram.curves[label], dmu_o).gamma_ev_per_angstrom2, rel=1e-12)


def test_a_compound_slab_refuses_the_default_chemical_potentials() -> None:
    # The elemental limits are not in equilibrium with SnO2, so the default is refused for every slab.
    with pytest.raises(ValueError, match="chemical_potentials_ev") as info:
        campaign_surface_energies(sno2_entries())
    assert all(label in str(info.value) for label in ("o_term", "minus1_o", "minus2_o"))


def test_elemental_limits_differ_from_the_o_rich_point_by_the_formation_enthalpy() -> None:
    # gamma(elemental limits) - gamma(O-rich) = n_Sn dH_f / (2 A), with dH_f = -5.20 eV per SnO2
    elemental = campaign_surface_energies(sno2_entries(), chemical_potentials_ev={"Sn": -4.15, "O": -4.93})
    o_rich = campaign_surface_energies(sno2_entries(), chemical_potentials_ev={"Sn": -9.35, "O": -4.93})
    for label in ("o_term", "minus1_o", "minus2_o"):
        assert elemental[label] - o_rich[label] == pytest.approx(6 * (-5.20) / (2 * 42.8), rel=1e-9)


def test_two_tin_references_need_a_choice_and_the_unchosen_one_is_excluded() -> None:
    entries = sno2_entries(extra=(reference("sn_alt", {"Sn": 1}, -3.9000),))
    with pytest.raises(ValueError, match="several references of Sn"):
        campaign_references(entries, host="sno2")
    # the unchosen tin is another solid of a binary oxide: an error until it is excluded
    with pytest.raises(ValueError, match=r"\['sn'\].*list them in exclude"):
        campaign_references(entries, host="sno2", reservoirs={"Sn": "sn_alt"})
    chosen = campaign_references(entries, host="sno2", reservoirs={"Sn": "sn_alt"}, exclude=("sn",))
    assert chosen.metal_energy_per_atom_ev == pytest.approx(-3.9000, rel=1e-12)


def test_the_oxygen_reservoir_is_the_o2_gas_in_every_reader() -> None:
    # An O atom (gas) does not compete with O2: the same reservoir in references and chemical potentials.
    entries = sno2_entries(extra=(reference("o_atom", {"O": 1}, -4.5000, phase="gas"),))
    assert campaign_references(entries, host="sno2") == hand_sno2_references()
    assert campaign_chemical_potentials(entries) == {"Sn": pytest.approx(-4.15), "O": pytest.approx(-4.93)}


def test_host_must_be_a_solid_with_a_phase_diagram_model() -> None:
    with pytest.raises(ValueError, match="must be a solid reference"):
        campaign_references(sno2_entries(), host="o2")
    with pytest.raises(ValueError, match="campaign_surface_energies"):
        campaign_references(gold_entries(), host="au_bulk")


# --- 1c. Ternary SrTiO3: SrO- and TiO2-terminated slabs ---------------------------------

SRTIO3_AREA_ANGSTROM2 = 15.2490  # one face of the 1x1 (001) slab, a^2 with a = 3.905 A


def srtio3_entries(extra=()):
    return (
        reference("srtio3", {"Sr": 1, "Ti": 1, "O": 3}, -40.7900),
        reference("sr", {"Sr": 2}, -3.2000),
        reference("ti", {"Ti": 2}, -15.8000),
        reference("o2", {"O": 2}, -9.8600, phase="gas"),
        reference("sro", {"Sr": 1, "O": 1}, -12.0300),
        reference("tio2", {"Ti": 1, "O": 2}, -27.1600),
        slab("sro_term", {"Sr": 4, "Ti": 3, "O": 10}, -133.6701, SRTIO3_AREA_ANGSTROM2),
        slab("tio2_term", {"Sr": 3, "Ti": 4, "O": 11}, -148.2001, SRTIO3_AREA_ANGSTROM2),
        *extra,
    )


def srtio3_terminations():
    return [
        SlabTermination("sro_term", -133.6701, {"Sr": 4, "Ti": 3, "O": 10}, SRTIO3_AREA_ANGSTROM2),
        SlabTermination("tio2_term", -148.2001, {"Sr": 3, "Ti": 4, "O": 11}, SRTIO3_AREA_ANGSTROM2),
    ]


def hand_srtio3_references(extra_phases=()) -> TernaryOxideReferences:
    return TernaryOxideReferences(
        bulk_energy_ev=-40.7900,
        bulk_composition={"Sr": 1, "Ti": 1, "O": 3},
        oxygen_molecule_energy_ev=-9.8600,
        element_energies_per_atom_ev={"Sr": -1.6000, "Ti": -7.9000},
        competing_phases=(
            CompetingPhase("sro", -12.0300, {"Sr": 1, "O": 1}),
            CompetingPhase("tio2", -27.1600, {"Ti": 1, "O": 2}),
            *extra_phases,
        ),
    )


def test_ternary_references_equal_the_hand_built_object() -> None:
    built = campaign_references(srtio3_entries(), host="srtio3")
    hand = hand_srtio3_references()
    assert isinstance(built, TernaryOxideReferences)
    assert built.bulk_energy_ev == pytest.approx(hand.bulk_energy_ev, rel=1e-15)
    assert dict(built.bulk_composition) == dict(hand.bulk_composition)
    assert built.oxygen_molecule_energy_ev == pytest.approx(hand.oxygen_molecule_energy_ev, rel=1e-15)
    assert dict(built.element_energies_per_atom_ev) == pytest.approx(dict(hand.element_energies_per_atom_ev))
    assert [(p.label, p.energy_ev, dict(p.composition)) for p in built.competing_phases] == [
        (p.label, p.energy_ev, dict(p.composition)) for p in hand.competing_phases
    ]
    assert (built.independent, built.eliminated, built.formula_unit) == ("Sr", "Ti", (1, 1, 3))
    # dH_f = -40.79 - (-1.6 - 7.9 + 3 * (-9.86 / 2)) = -16.5 eV per SrTiO3
    assert built.formation_enthalpy_ev == pytest.approx(-16.5, rel=1e-12)
    assert built == hand


def test_ternary_stability_polygon_matches_the_hand_clipping() -> None:
    # Triangle dmu_Sr <= 0, dmu_O <= 0, dmu_Sr + 3 dmu_O >= -16.5, cut by SrO (dmu_Sr + dmu_O <= -5.5)
    # and TiO2 (dmu_Sr + dmu_O >= -7.1): the quadrilateral below.
    region = campaign_references(srtio3_entries(), host="srtio3").stability_region
    assert same_polygon(region, [(-7.1, 0.0), (-5.5, 0.0), (0.0, -5.5), (-2.4, -4.7)])
    assert polygon_area(region) == pytest.approx(8.16, rel=1e-9)


def test_an_unchosen_lower_energy_sr_allotrope_shrinks_the_stability_region() -> None:
    # sr_fcc: 2 atoms at -3.30 eV, 0.05 eV/atom below the chosen Sr reservoir sr (2 atoms, -3.20 eV).
    # As a competing phase Sr2: 2 dmu_Sr <= E_fcc - 2 E_sr = -3.30 + 3.20 = -0.10, so dmu_Sr <= -0.05 eV.
    entries = srtio3_entries(extra=(reference("sr_fcc", {"Sr": 2}, -3.3000),))
    refs = campaign_references(entries, host="srtio3", reservoirs={"Sr": "sr"})
    assert [phase.label for phase in refs.competing_phases] == ["sro", "tio2", "sr_fcc"]
    assert refs == hand_srtio3_references((CompetingPhase("sr_fcc", -3.3000, {"Sr": 2}),))
    # Hand clipping: the corner (0, -5.5) has dmu_Sr = 0 > -0.05 and goes.  dmu_Sr = -0.05 meets the SrO line
    # at dmu_O = -5.45 and the bulk line dmu_Sr + 3 dmu_O = -16.5 at dmu_O = (-16.5 + 0.05) / 3 = -16.45 / 3.
    region = refs.stability_region
    assert same_polygon(region, [(-7.1, 0.0), (-5.5, 0.0), (-0.05, -5.45), (-0.05, -16.45 / 3), (-2.4, -4.7)])
    assert "sr_fcc" in refs.stability_boundaries
    # the sliver (0, -5.5), (-0.05, -5.45), (-0.05, -16.45/3) has area 0.5 * (1/30) * 0.05 = 1/1200
    assert polygon_area(region) == pytest.approx(8.16 - 1 / 1200, rel=1e-9)
    assert not refs.in_stability_region(-0.02, -5.49)  # inside the region without the allotrope
    assert refs.in_stability_region(-1.0, -4.8)
    # without a choice the two Sr references are an error, not a guess
    with pytest.raises(ValueError, match="several references of Sr"):
        campaign_references(entries, host="srtio3")


def test_an_unchosen_lower_energy_ti_allotrope_shrinks_the_region_through_the_eliminated_element() -> None:
    # ti_alt: 2 atoms at -16.00 eV, 0.10 eV/atom below ti.  dmu_Ti = dH_f - dmu_Sr - 3 dmu_O (x = y = 1, z = 3),
    # so Ti2 gives 2 dmu_Ti <= -0.20, i.e. dmu_Sr + 3 dmu_O >= -16.4: the bulk line moves inward.
    entries = srtio3_entries(extra=(reference("ti_alt", {"Ti": 2}, -16.0000),))
    refs = campaign_references(entries, host="srtio3", reservoirs={"Ti": "ti"})
    # hand: the bulk line -16.4 meets the SrO line at (-0.05, -5.45) and the TiO2 line at (-2.45, -4.65)
    assert same_polygon(refs.stability_region, [(-7.1, 0.0), (-5.5, 0.0), (-0.05, -5.45), (-2.45, -4.65)])
    assert polygon_area(refs.stability_region) == pytest.approx(8.08, rel=1e-9)
    assert not refs.in_stability_region(-0.025, -5.475)
    assert refs.in_stability_region(-0.3, -5.3)


def test_choice_of_independent_element_changes_the_axes_not_gamma() -> None:
    # docs/source/phase-diagram.rst: the same state has the same gamma in either representation.
    sr_axes = ternary_surface_phase_diagram(srtio3_terminations(), campaign_references(srtio3_entries(), host="srtio3"))
    ti_axes = ternary_surface_phase_diagram(
        srtio3_terminations(), campaign_references(srtio3_entries(), host="srtio3", independent="Ti")
    )
    for dmu_sr, dmu_o in ((-4.0, -2.0), (-3.9, -3.0)):
        dmu_ti = -16.5 - dmu_sr - 3 * dmu_o  # bulk equilibrium: dmu_Sr + dmu_Ti + 3 dmu_O = dH_f
        for label in ("sro_term", "tio2_term"):
            assert ti_axes.gamma_ev_per_angstrom2(label, dmu_ti, dmu_o) == pytest.approx(
                sr_axes.gamma_ev_per_angstrom2(label, dmu_sr, dmu_o), rel=1e-12
            )


def test_ternary_stable_termination_matches_the_hand_gamma_at_two_points() -> None:
    built = ternary_surface_phase_diagram(srtio3_terminations(), campaign_references(srtio3_entries(), host="srtio3"))
    hand = ternary_surface_phase_diagram(srtio3_terminations(), hand_srtio3_references())
    for dmu_sr, dmu_o, expected in ((-4.0, -2.0, "sro_term"), (-3.9, -3.0, "tio2_term")):
        # mu_Sr = E_Sr + dmu_Sr, mu_Ti = E_Ti + (dH_f - dmu_Sr - 3 dmu_O), mu_O = E(O2)/2 + dmu_O
        mu_sr = -1.6000 + dmu_sr
        mu_ti = -7.9000 + (-16.5 - dmu_sr - 3 * dmu_o)
        mu_o = -9.8600 / 2 + dmu_o
        gamma_sro = (-133.6701 - 4 * mu_sr - 3 * mu_ti - 10 * mu_o) / (2 * SRTIO3_AREA_ANGSTROM2)
        gamma_tio2 = (-148.2001 - 3 * mu_sr - 4 * mu_ti - 11 * mu_o) / (2 * SRTIO3_AREA_ANGSTROM2)
        assert built.gamma_ev_per_angstrom2("sro_term", dmu_sr, dmu_o) == pytest.approx(gamma_sro, rel=1e-12)
        assert built.gamma_ev_per_angstrom2("tio2_term", dmu_sr, dmu_o) == pytest.approx(gamma_tio2, rel=1e-12)
        winner = "sro_term" if gamma_sro < gamma_tio2 else "tio2_term"
        assert winner == expected
        assert built.stable_termination(dmu_sr, dmu_o) == expected == hand.stable_termination(dmu_sr, dmu_o)
        assert built.references.in_stability_region(dmu_sr, dmu_o)


# --- 1d. Intermetallic PdIn: no oxide model, surface energies at the limits ---------------

PDIN_AREA_ANGSTROM2 = 19.0


def pdin_entries():
    return (
        reference("pd", {"Pd": 4}, -21.2000),
        reference("in", {"In": 2}, -6.0000),
        reference("pdin", {"Pd": 1, "In": 1}, -8.9000),
        slab("pdin_pd_rich", {"Pd": 9, "In": 7}, -69.6700, PDIN_AREA_ANGSTROM2),
    )


def test_intermetallic_has_no_oxide_model_and_points_to_the_surface_energy_api() -> None:
    with pytest.raises(ValueError, match="campaign_surface_energies") as info:
        campaign_references(pdin_entries(), host="pdin")
    assert "campaign_chemical_potentials" in str(info.value)


def test_intermetallic_elemental_limits_are_the_elemental_energies_per_atom() -> None:
    assert campaign_chemical_potentials(pdin_entries()) == {"Pd": pytest.approx(-5.30), "In": pytest.approx(-3.00)}


def test_intermetallic_gamma_at_the_pd_rich_limit_matches_the_hand_formula() -> None:
    # Pd-rich: mu_Pd = E_Pd (elemental) and mu_In from PdIn equilibrium, mu_In = E_fu - mu_Pd, E_fu = -8.90 eV
    mu_pd = -21.2000 / 4
    mu_in = -8.9000 - mu_pd
    expected = (-69.6700 - 9 * mu_pd - 7 * mu_in) / (2 * 19.0)  # gamma = (E - sum N_i mu_i) / (2 A)
    assert expected == pytest.approx(0.085, abs=1e-12)
    gamma = campaign_surface_energies(pdin_entries(), chemical_potentials_ev={"Pd": mu_pd, "In": mu_in})
    assert gamma["pdin_pd_rich"] == pytest.approx(expected, rel=1e-12)


def test_intermetallic_default_is_refused_and_the_elemental_limits_give_the_hand_value() -> None:
    with pytest.raises(ValueError, match="pdin_pd_rich.*chemical_potentials_ev"):
        campaign_surface_energies(pdin_entries())
    elemental = {"Pd": -21.2000 / 4, "In": -6.0000 / 2}  # both elemental limits, eV per atom
    expected = (-69.6700 - 9 * elemental["Pd"] - 7 * elemental["In"]) / (2 * 19.0)
    gamma = campaign_surface_energies(pdin_entries(), chemical_potentials_ev=elemental)
    assert gamma["pdin_pd_rich"] == pytest.approx(expected, rel=1e-12)
