"""Tier-1 tests for the campaign analysis: entries, chemical potentials, oxide references, surface energies.

The entries are typed by hand, so nothing here needs an AiiDA profile.  Each
expected number is computed in the test, with its formula written next to it.
"""

from __future__ import annotations

import dataclasses
import inspect
import re

import pytest

import psteros


def reference(label, composition, energy_ev, phase="solid", **options):
    """A reference entry: a bulk solid unless ``phase`` says gas."""

    return psteros.CampaignEntry(label, "references", phase, composition, energy_ev, **options)


def slab(label, composition, energy_ev, area, **options):
    """A slab entry; ``area`` is the area of one face in A^2."""

    return psteros.CampaignEntry(label, "slabs", "solid", composition, energy_ev, surface_area_angstrom2=area, **options)


def as_generator(items):
    return (item for item in items)


# Binary oxide SnO2.  mu_Sn = E(Sn8) / 8 = -32.0 / 8 = -4.0 eV and mu_O = E(O2) / 2 = -9.0 / 2 = -4.5 eV.
# dH(SnO2) per formula unit = E(Sn2O4) / 2 - mu_Sn - 2 mu_O = -19.0 + 4.0 + 9.0 = -6.0 eV.
SN_BULK = reference("sn_bulk", {"Sn": 8}, -32.0)
SN_ALT = reference("sn_alt", {"Sn": 4}, -15.0)  # mu_Sn = -15.0 / 4 = -3.75 eV
O2 = reference("o2", {"O": 2}, -9.0, phase="gas")
O2_ALT = reference("o2_alt", {"O": 2}, -8.8, phase="gas")  # mu_O = -8.8 / 2 = -4.4 eV
O_ATOM = reference("o_atom", {"O": 1}, -4.0, phase="gas")
SNO2 = reference("sno2", {"Sn": 2, "O": 4}, -38.0)
SNO2_UNFINISHED = reference("sno2", {"Sn": 2, "O": 4}, None, state="not started")
SNO = reference("sno", {"Sn": 1, "O": 1}, -12.0)  # a multi-element solid: no competing phases in a binary oxide
SNO_GAS = reference("sno_gas", {"Sn": 1, "O": 1}, -6.0, phase="gas")  # a compound gas inside the host elements
H2O = reference("h2o", {"H": 2, "O": 1}, -14.0, phase="gas")  # a compound gas with a foreign element
PD_BULK = reference("pd_bulk", {"Pd": 4}, -20.0)  # a metal outside the host
SN_SLAB = slab("sn_slab", {"Sn": 4}, -10.0, area=8.0)  # unary slab: gamma = (-10.0 + 4 * 4.0) / 16.0 = 0.375
O_SLAB = slab("o_slab", {"O": 2}, -9.0, area=8.0)
O_SLAB_UNARY = slab("o_slab_u", {"O": 2}, -8.0, area=8.0)  # gamma = (-8.0 + 2 * 4.5) / 16.0 = 0.0625
SLAB_O = slab("slab_o", {"Sn": 4, "O": 8}, -40.0, area=16.0)

# Unary, intermetallic and nitride hosts.  mu_Au = -13.0 / 4 = -3.25 eV.
AU_BULK = reference("au_bulk", {"Au": 4}, -13.0)
AU_ALT = reference("au_alt", {"Au": 2}, -6.0)  # mu_Au = -3.0 eV
AU_SLAB = slab("au_slab", {"Au": 12}, -38.0, area=8.0)
PDIN = reference("pdin", {"Pd": 1, "In": 1}, -10.0)
TIN = reference("tin", {"Ti": 1, "N": 1}, -8.0)

# Ternary oxide SrTiO3.  mu_Sr = -2.0 / 2 = -1.0, mu_Ti = -4.0 / 2 = -2.0 and mu_O = -8.0 / 2 = -4.0 eV.
# dH(SrTiO3) = -32.0 - (1 * -1.0 + 1 * -2.0 + 3 * -4.0) = -17.0 eV.
# dH(SrO) = -11.1 - (-1.0 - 4.0) = -6.1 eV and dH(TiO2) = -19.8 - (-2.0 - 8.0) = -9.8 eV.
SRTIO3 = reference("srtio3", {"Sr": 1, "Ti": 1, "O": 3}, -32.0)
SR = reference("sr", {"Sr": 2}, -2.0)
TI = reference("ti", {"Ti": 2}, -4.0)
O2_TERNARY = reference("o2", {"O": 2}, -8.0, phase="gas")
SRO = reference("sro", {"Sr": 1, "O": 1}, -11.1)
TIO2 = reference("tio2", {"Ti": 1, "O": 2}, -19.8)
TERNARY = [SRTIO3, SR, TI, O2_TERNARY, SRO, TIO2]

MU_SN_O = {"Sn": -4.0, "O": -4.5}


# =============================================================================
# CampaignEntry: validation and properties
# =============================================================================


def test_defaults_are_the_static_block_a_finished_state_and_no_area() -> None:
    assert (SNO2.block, SNO2.state, SNO2.surface_area_angstrom2) == ("static", "finished", None)
    with pytest.raises(dataclasses.FrozenInstanceError):
        SNO2.energy_ev = -1.0


def test_atoms_and_energy_per_atom_follow_the_composition() -> None:
    assert SNO2.atoms == 2 + 4
    assert SNO2.energy_per_atom_ev == pytest.approx(-38.0 / 6)
    assert SNO2_UNFINISHED.energy_per_atom_ev is None


@pytest.mark.parametrize(
    ("group", "phase"),
    [("systems", "solid"), ("reference", "solid"), ("references", "liquid"), ("slabs", "liquid")],
)
def test_group_and_phase_must_be_known(group: str, phase: str) -> None:
    with pytest.raises(ValueError, match="sno2"):
        psteros.CampaignEntry("sno2", group, phase, {"Sn": 2, "O": 4}, -38.0)


def test_a_slab_is_a_solid_not_a_gas() -> None:
    with pytest.raises(ValueError, match="slab_o"):
        psteros.CampaignEntry("slab_o", "slabs", "gas", {"O": 2}, -9.0, surface_area_angstrom2=8.0)


@pytest.mark.parametrize("composition", [{}, {"Sn": 0}, {"O": -2}, {"Sn": 1.5}, {"Sn": 2, "O": -4}])
def test_composition_must_hold_positive_integer_counts(composition) -> None:
    with pytest.raises(ValueError, match="sno2"):
        psteros.CampaignEntry("sno2", "references", "solid", composition, -38.0)


def test_integer_valued_counts_are_stored_as_int() -> None:
    entry = psteros.CampaignEntry("sno2", "references", "solid", {"Sn": 2.0, "O": 4.0}, -38.0)
    assert dict(entry.composition) == {"Sn": 2, "O": 4}
    assert all(type(count) is int for count in entry.composition.values())


def test_a_zero_count_is_dropped_from_the_composition() -> None:
    # A zero count is dropped (the element is absent) rather than rejected; the spec wording allows either.
    entry = psteros.CampaignEntry("sno2", "references", "solid", {"Sn": 2, "O": 4, "H": 0}, -38.0)
    assert dict(entry.composition) == {"Sn": 2, "O": 4}


@pytest.mark.parametrize("area", [0.0, -3.0])
def test_a_slab_area_must_be_positive(area: float) -> None:
    with pytest.raises(ValueError, match="slab_o"):
        slab("slab_o", {"Sn": 4, "O": 8}, -40.0, area=area)


@pytest.mark.parametrize("area", [float("nan"), float("inf")])
def test_a_non_finite_surface_area_is_an_error(area: float) -> None:
    with pytest.raises(ValueError, match="slab_o"):
        slab("slab_o", {"Sn": 4, "O": 8}, -40.0, area=area)


@pytest.mark.parametrize("energy", [float("nan"), float("inf"), float("-inf")])
def test_a_non_finite_energy_is_an_error(energy: float) -> None:
    with pytest.raises(ValueError, match="finite") as info:
        psteros.CampaignEntry("sno2", "references", "solid", {"Sn": 2, "O": 4}, energy)
    assert "sno2" in str(info.value)


def test_a_slab_area_may_be_unknown_until_its_structure_is_known() -> None:
    entry = slab("slab_o", {"Sn": 4, "O": 8}, None, area=None, state="not started")
    assert entry.surface_area_angstrom2 is None


def test_a_reference_has_no_surface_area() -> None:
    with pytest.raises(ValueError, match="sno2"):
        reference("sno2", {"Sn": 2, "O": 4}, -38.0, surface_area_angstrom2=16.0)


# =============================================================================
# campaign_chemical_potentials
# =============================================================================


def test_a_unary_chemical_potential_is_the_energy_per_atom_of_the_bulk() -> None:
    # mu_Au = E / N = -13.0 / 4
    assert psteros.campaign_chemical_potentials([AU_BULK]) == pytest.approx({"Au": -13.0 / 4})


def test_a_binary_oxide_gives_tin_per_atom_and_oxygen_from_o2() -> None:
    # mu_Sn = E(Sn8) / 8 = -32.0 / 8 (solid); mu_O = E(O2) / 2 = -9.0 / 2 (gas)
    result = psteros.campaign_chemical_potentials([SN_BULK, O2, SNO2])
    assert result == pytest.approx({"Sn": -32.0 / 8, "O": -9.0 / 2})


def test_an_element_without_a_single_element_reference_is_absent() -> None:
    assert psteros.campaign_chemical_potentials([SN_BULK, SNO2]) == pytest.approx({"Sn": -4.0})


def test_two_references_of_one_element_need_a_reservoir_choice() -> None:
    with pytest.raises(ValueError) as info:
        psteros.campaign_chemical_potentials([SN_BULK, O2, O2_ALT, SNO2])
    assert re.search(r"\bo2\b", str(info.value)) and re.search(r"\bo2_alt\b", str(info.value))


def test_the_reservoir_choice_picks_one_candidate() -> None:
    entries = [SN_BULK, O2, O2_ALT, SNO2]
    assert psteros.campaign_chemical_potentials(entries, reservoirs={"O": "o2"}) == pytest.approx(
        {"Sn": -4.0, "O": -9.0 / 2}
    )
    assert psteros.campaign_chemical_potentials(entries, reservoirs={"O": "o2_alt"}) == pytest.approx(
        {"Sn": -4.0, "O": -8.8 / 2}
    )
    # A reservoir for another element does not settle the two oxygen candidates.
    with pytest.raises(ValueError, match="o2"):
        psteros.campaign_chemical_potentials(entries, reservoirs={"Sn": "sn_bulk"})


@pytest.mark.parametrize("label", ["sno2", "sn_bulk", "o_slab", "nope"])
def test_a_reservoir_must_be_a_single_element_reference_of_that_element(label: str) -> None:
    with pytest.raises(ValueError, match=label):
        psteros.campaign_chemical_potentials([SN_BULK, O2, SNO2, O_SLAB], reservoirs={"O": label})


def test_an_unfinished_oxygen_reference_is_named_with_its_state() -> None:
    o2_unfinished = reference("o2", {"O": 2}, None, phase="gas", state="not started")
    with pytest.raises(ValueError, match="o2") as info:
        psteros.campaign_chemical_potentials([SN_BULK, o2_unfinished])
    assert "not started" in str(info.value)


def test_an_unfinished_metal_reference_is_named_with_its_state() -> None:
    sn_running = reference("sn_bulk", {"Sn": 8}, None, state="running")
    with pytest.raises(ValueError, match="sn_bulk") as info:
        psteros.campaign_chemical_potentials([sn_running, O2])
    assert "running" in str(info.value)


def test_multi_element_references_are_not_used() -> None:
    # The unfinished SnO2 would be an error if it were a single-element reference.
    assert psteros.campaign_chemical_potentials([SN_BULK, O2, SNO2_UNFINISHED]) == pytest.approx(MU_SN_O)


def test_a_single_element_slab_is_not_a_reference() -> None:
    # A tin slab would be a second tin candidate if slabs counted as references.
    assert psteros.campaign_chemical_potentials([SN_BULK, O2, SN_SLAB]) == pytest.approx(MU_SN_O)


def test_the_oxygen_candidates_are_the_o2_gas_references_when_there_are_any() -> None:
    # An O atom does not compete with O2 for the oxygen reservoir.
    assert psteros.campaign_chemical_potentials([SN_BULK, O2, O_ATOM]) == pytest.approx(MU_SN_O)
    with pytest.raises(ValueError) as info:
        psteros.campaign_chemical_potentials([SN_BULK, O_ATOM, O2, O2_ALT])
    message = str(info.value)
    assert re.search(r"\bo2\b", message) and re.search(r"\bo2_alt\b", message)
    assert "o_atom" not in message


def test_without_o2_the_single_element_oxygen_references_are_the_candidates() -> None:
    # The O atom is then the oxygen reference: mu_O = E(O) / 1 = -4.0 eV.
    assert psteros.campaign_chemical_potentials([SN_BULK, O_ATOM]) == pytest.approx({"Sn": -4.0, "O": -4.0})
    o_atom_alt = reference("o_atom_alt", {"O": 1}, -3.9, phase="gas")
    with pytest.raises(ValueError) as info:
        psteros.campaign_chemical_potentials([SN_BULK, O_ATOM, o_atom_alt])
    assert re.search(r"\bo_atom\b", str(info.value)) and re.search(r"\bo_atom_alt\b", str(info.value))


def test_elements_limits_the_result_and_its_checks() -> None:
    # The ambiguous oxygen references and the unfinished O2 are not needed for tin.
    assert psteros.campaign_chemical_potentials([SN_BULK, O2, O2_ALT], elements=["Sn"]) == pytest.approx({"Sn": -4.0})
    unfinished_o2 = reference("o2", {"O": 2}, None, phase="gas", state="not started")
    assert psteros.campaign_chemical_potentials([SN_BULK, unfinished_o2], elements=("Sn",)) == pytest.approx(
        {"Sn": -4.0}
    )
    # Without elements the same unfinished O2 is an error.
    with pytest.raises(ValueError, match="o2"):
        psteros.campaign_chemical_potentials([SN_BULK, unfinished_o2])
    # A reservoir choice applies to its element, and the result holds that element only.
    result = psteros.campaign_chemical_potentials(
        [SN_BULK, O2, O2_ALT], elements=["O"], reservoirs={"O": "o2_alt"}
    )
    assert result == pytest.approx({"O": -4.4})


def test_an_element_asked_for_without_a_reference_is_an_error() -> None:
    with pytest.raises(ValueError, match=r"no single-element reference of \['Au'\]"):
        psteros.campaign_chemical_potentials([SN_BULK, O2], elements=["Sn", "Au"])
    # Without ``elements``, only the elements that have a reference are returned.
    assert psteros.campaign_chemical_potentials([SN_BULK, O2]) == pytest.approx({"Sn": -4.0, "O": -4.5})


# =============================================================================
# campaign_references
# =============================================================================


def test_a_binary_host_gives_binary_oxide_references_with_exact_values() -> None:
    refs = psteros.campaign_references([SN_BULK, O2, SNO2], host="sno2")
    assert isinstance(refs, psteros.BinaryOxideReferences)
    assert refs.bulk_energy_ev == -38.0
    assert dict(refs.bulk_composition) == {"Sn": 2, "O": 4}
    assert refs.oxygen_molecule_energy_ev == -9.0
    assert refs.metal_energy_per_atom_ev == -32.0 / 8
    assert refs.formation_enthalpy_ev == pytest.approx(-38.0 / 2 + 4.0 + 9.0)  # -6.0 eV per SnO2
    assert refs == psteros.BinaryOxideReferences(-38.0, {"Sn": 2, "O": 4}, -9.0, -4.0)


def test_a_binary_host_without_a_metal_reference_has_no_metal_energy() -> None:
    refs = psteros.campaign_references([O2, SNO2], host="sno2")
    assert refs.metal_energy_per_atom_ev is None
    assert refs.formation_enthalpy_ev is None
    assert refs.oxygen_poor_limit_ev is None


def test_a_ternary_host_takes_the_elements_and_the_solid_compounds() -> None:
    refs = psteros.campaign_references(TERNARY, host="srtio3")
    assert isinstance(refs, psteros.TernaryOxideReferences)
    assert refs.bulk_energy_ev == -32.0
    assert dict(refs.bulk_composition) == {"Sr": 1, "Ti": 1, "O": 3}
    assert refs.oxygen_molecule_energy_ev == -8.0
    assert dict(refs.element_energies_per_atom_ev) == {"Sr": -2.0 / 2, "Ti": -4.0 / 2}
    # SrO and TiO2 compete; the elemental references sr and ti are not competing phases.
    assert [(phase.label, phase.energy_ev, dict(phase.composition)) for phase in refs.competing_phases] == [
        ("sro", -11.1, {"Sr": 1, "O": 1}),
        ("tio2", -19.8, {"Ti": 1, "O": 2}),
    ]
    assert (refs.independent, refs.eliminated) == ("Sr", "Ti")
    assert refs.formation_enthalpy_ev == pytest.approx(-32.0 - (1 * -1.0 + 1 * -2.0 + 3 * -4.0))  # -17.0 eV
    expected = psteros.TernaryOxideReferences(
        bulk_energy_ev=-32.0,
        bulk_composition={"Sr": 1, "Ti": 1, "O": 3},
        oxygen_molecule_energy_ev=-8.0,
        element_energies_per_atom_ev={"Sr": -1.0, "Ti": -2.0},
        competing_phases=(
            psteros.CompetingPhase("sro", -11.1, {"Sr": 1, "O": 1}),
            psteros.CompetingPhase("tio2", -19.8, {"Ti": 1, "O": 2}),
        ),
    )
    assert refs == expected


def test_the_ternary_competing_phases_include_an_unchosen_elemental_phase() -> None:
    # A second Sr allotrope (sr_alt: 2 atoms at -1.8 eV) is not the chosen reservoir, so it competes.
    sr_alt = reference("sr_alt", {"Sr": 2}, -1.8)
    refs = psteros.campaign_references([*TERNARY, sr_alt], host="srtio3", reservoirs={"Sr": "sr"})
    assert sorted(phase.label for phase in refs.competing_phases) == ["sr_alt", "sro", "tio2"]
    competitor = next(phase for phase in refs.competing_phases if phase.label == "sr_alt")
    assert dict(competitor.composition) == {"Sr": 2}
    # Without a reservoir choice the two Sr references are ambiguous.
    with pytest.raises(ValueError, match="sr_alt"):
        psteros.campaign_references([*TERNARY, sr_alt], host="srtio3")


def test_the_independent_element_is_passed_through() -> None:
    refs = psteros.campaign_references(TERNARY, host="srtio3", independent="Ti")
    assert (refs.independent, refs.eliminated, refs.formula) == ("Ti", "Sr", "TiSrO3")
    with pytest.raises(ValueError, match="independent"):
        psteros.campaign_references(TERNARY, host="srtio3", independent="O")


def test_independent_is_an_error_for_a_binary_host() -> None:
    with pytest.raises(ValueError, match="independent") as info:
        psteros.campaign_references([SN_BULK, O2, SNO2], host="sno2", independent="Sn")
    assert "binary" in str(info.value)


def test_the_host_must_be_a_solid_reference() -> None:
    with pytest.raises(ValueError, match="nope"):
        psteros.campaign_references([SN_BULK, O2, SNO2], host="nope")
    with pytest.raises(ValueError, match="o2"):
        psteros.campaign_references([SN_BULK, O2, SNO2], host="o2")
    with pytest.raises(ValueError, match="slab_o"):
        psteros.campaign_references([SN_BULK, O2, SNO2, SLAB_O], host="slab_o")


def test_the_host_cannot_be_excluded() -> None:
    with pytest.raises(ValueError, match="cannot be excluded"):
        psteros.campaign_references([SN_BULK, O2, SNO2], host="sno2", exclude=("sno2",))


def test_a_reservoir_cannot_also_be_excluded() -> None:
    with pytest.raises(ValueError, match="also in exclude"):
        psteros.campaign_references(
            [SN_BULK, O2, SNO2], host="sno2", reservoirs={"Sn": "sn_bulk"}, exclude=("sn_bulk",)
        )


def test_an_oxygen_gas_reference_is_required() -> None:
    with pytest.raises(ValueError, match=r"(?i)oxygen|\bO2\b"):
        psteros.campaign_references([SN_BULK, SNO2], host="sno2")


@pytest.mark.parametrize(
    "oxygen",
    [O_ATOM, reference("o4", {"O": 4}, -18.0, phase="gas"), reference("o2_solid", {"O": 2}, -9.0)],
    ids=["atom", "o4", "solid"],
)
def test_the_oxygen_reservoir_is_a_gas_of_exactly_two_atoms(oxygen) -> None:
    with pytest.raises(ValueError, match=r"(?i)oxygen|\bO2\b"):
        psteros.campaign_references([SN_BULK, oxygen, SNO2], host="sno2")


def test_an_oxygen_atom_next_to_o2_does_not_change_the_reservoir() -> None:
    # The reservoir is the gas of exactly {"O": 2}; an O atom is not a second candidate.
    refs = psteros.campaign_references([SN_BULK, O2, O_ATOM, SNO2], host="sno2")
    assert refs.oxygen_molecule_energy_ev == -9.0


def test_an_o_atom_alone_is_no_oxygen_reservoir_of_a_binary_oxide() -> None:
    with pytest.raises(ValueError) as info:
        psteros.campaign_references([SN_BULK, O_ATOM, SNO2], host="sno2")
    assert str(info.value).startswith("'sno2': ")
    assert "O2" in str(info.value)


def test_the_oxygen_error_lists_the_gas_references() -> None:
    with pytest.raises(ValueError) as info:
        psteros.campaign_references([SN_BULK, O_ATOM, SNO2], host="sno2")
    assert "gas references: ['o_atom']" in str(info.value)
    with pytest.raises(ValueError) as info:
        psteros.campaign_references([SN_BULK, SNO2], host="sno2")
    assert "gas references: []" in str(info.value)


def test_two_o2_gas_references_need_a_reservoir_choice() -> None:
    entries = [SN_BULK, O2, O2_ALT, SNO2]
    with pytest.raises(ValueError) as info:
        psteros.campaign_references(entries, host="sno2")
    assert re.search(r"\bo2\b", str(info.value)) and re.search(r"\bo2_alt\b", str(info.value))
    refs = psteros.campaign_references(entries, host="sno2", reservoirs={"O": "o2_alt"})
    assert refs.oxygen_molecule_energy_ev == -8.8


@pytest.mark.parametrize("element", ["Sr", "Ti"])
def test_a_ternary_host_needs_both_element_references(element: str) -> None:
    entries = [entry for entry in TERNARY if entry.label != element.lower()]
    with pytest.raises(ValueError, match=element):
        psteros.campaign_references(entries, host="srtio3")


@pytest.mark.parametrize("host", [AU_BULK, PDIN, TIN], ids=["unary", "intermetallic", "nitride"])
def test_a_host_that_is_not_an_oxide_points_to_the_elemental_functions(host) -> None:
    with pytest.raises(ValueError, match=r"(?i)no phase.diagram model") as info:
        psteros.campaign_references([host, O2], host=host.label)
    assert "campaign_chemical_potentials" in str(info.value)
    assert "campaign_surface_energies" in str(info.value)


@pytest.mark.parametrize("extra", [H2O, PD_BULK, SNO_GAS, SNO], ids=["h2o", "pd_bulk", "sno_gas", "sno"])
def test_a_reference_outside_the_oxide_model_is_an_error_until_excluded(extra) -> None:
    base = [SN_BULK, O2, SNO2]
    with pytest.raises(ValueError, match=rf"\b{extra.label}\b"):
        psteros.campaign_references([*base, extra], host="sno2")
    excluded = psteros.campaign_references([*base, extra], host="sno2", exclude=(extra.label,))
    assert excluded == psteros.campaign_references(base, host="sno2")


def test_a_binary_oxide_rejects_an_unchosen_elemental_phase_until_excluded() -> None:
    entries = [SN_BULK, SN_ALT, O2, SNO2]
    with pytest.raises(ValueError, match="sn_alt") as info:
        psteros.campaign_references(entries, host="sno2", reservoirs={"Sn": "sn_bulk"})
    assert "exclude" in str(info.value)
    excluded = psteros.campaign_references(entries, host="sno2", reservoirs={"Sn": "sn_bulk"}, exclude=("sn_alt",))
    assert excluded == psteros.campaign_references([SN_BULK, O2, SNO2], host="sno2")


def test_a_binary_oxide_rejects_a_solid_oxygen_reference() -> None:
    o_solid = reference("o_solid", {"O": 2}, -9.0)
    with pytest.raises(ValueError, match="o_solid"):
        psteros.campaign_references([SN_BULK, O2, o_solid, SNO2], host="sno2")


def test_errors_of_the_model_are_prefixed_with_the_host() -> None:
    # dH(SnO2) = -20.0 / 2 + 4.0 + 9.0 = +3.0 eV: the oxide is unstable against Sn + O2.
    unstable = reference("sno2", {"Sn": 2, "O": 4}, -20.0)
    with pytest.raises(ValueError, match="not negative") as info:
        psteros.campaign_references([SN_BULK, O2, unstable], host="sno2")
    assert str(info.value).startswith("'sno2': ")
    with pytest.raises(ValueError, match="needs single-element") as info:
        psteros.campaign_references([entry for entry in TERNARY if entry.label != "ti"], host="srtio3")
    assert str(info.value).startswith("'srtio3': ")
    with pytest.raises(ValueError, match="other solids") as info:
        psteros.campaign_references([SN_BULK, SN_ALT, O2, SNO2], host="sno2", reservoirs={"Sn": "sn_bulk"})
    assert str(info.value).startswith("'sno2': ")


def test_unknown_labels_in_exclude_or_reservoirs_are_errors() -> None:
    entries = [SN_BULK, O2, SNO2]
    with pytest.raises(ValueError, match="nope"):
        psteros.campaign_references(entries, host="sno2", exclude=("nope",))
    with pytest.raises(ValueError, match="nope"):
        psteros.campaign_references(entries, host="sno2", reservoirs={"O": "nope"})


def test_an_unfinished_host_is_named_with_its_state() -> None:
    with pytest.raises(ValueError, match="sno2") as info:
        psteros.campaign_references([SN_BULK, O2, SNO2_UNFINISHED], host="sno2")
    assert "not started" in str(info.value)


# =============================================================================
# campaign_surface_energies
# =============================================================================


def test_a_unary_slab_is_exact_at_the_elemental_limit() -> None:
    # gamma = (E_slab - N mu) / (surfaces * A) = (-38.0 - 12 * (-13.0 / 4)) / (2 * 8.0) = 1.0 / 16
    assert psteros.campaign_surface_energies([AU_BULK, AU_SLAB]) == pytest.approx({"au_slab": 0.0625})


def test_given_chemical_potentials_replace_the_elemental_limits() -> None:
    # No reference is needed: gamma = (-38.0 - 12 * -3.0) / (2 * 8.0) = -2.0 / 16
    result = psteros.campaign_surface_energies([AU_SLAB], chemical_potentials_ev={"Au": -3.0})
    assert result == pytest.approx({"au_slab": -0.125})


def test_surfaces_is_the_number_of_faces() -> None:
    # gamma = (-38.0 + 39.0) / (1 * 8.0)
    assert psteros.campaign_surface_energies([AU_BULK, AU_SLAB], surfaces=1) == pytest.approx({"au_slab": 0.125})


def test_a_compound_slab_is_gamma_at_the_given_chemical_potentials() -> None:
    # gamma = (-40.0 - (4 * -4.0 + 8 * -4.5)) / (2 * 16.0) = 12.0 / 32.0
    assert psteros.campaign_surface_energies([SLAB_O], chemical_potentials_ev=MU_SN_O) == pytest.approx({"slab_o": 0.375})


def test_a_compound_slab_needs_explicit_chemical_potentials() -> None:
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_surface_energies([SN_BULK, O2, SLAB_O])
    assert "chemical_potentials_ev" in str(info.value)


def test_every_compound_slab_is_named_when_chemical_potentials_are_missing() -> None:
    other = slab("slab_x", {"Sn": 2, "O": 4}, -20.0, area=8.0)
    with pytest.raises(ValueError) as info:
        psteros.campaign_surface_energies([SN_BULK, O2, SLAB_O, other])
    message = str(info.value)
    assert re.search(r"\bslab_o\b", message) and re.search(r"\bslab_x\b", message)


def test_the_reservoir_sets_the_elemental_limit_of_a_unary_slab() -> None:
    # gamma = (E - 4 mu) / (2 * 8.0) with E = -10.0: (-10.0 + 16.0) / 16.0 = 0.375 for sn_bulk,
    # and (-10.0 + 15.0) / 16.0 = 0.3125 for sn_alt (mu = -3.75 eV).
    entries = [SN_BULK, SN_ALT, SN_SLAB]
    assert psteros.campaign_surface_energies(entries, reservoirs={"Sn": "sn_alt"}) == pytest.approx({"sn_slab": 0.3125})
    assert psteros.campaign_surface_energies(entries, reservoirs={"Sn": "sn_bulk"}) == pytest.approx({"sn_slab": 0.375})
    with pytest.raises(ValueError, match="sn_alt"):
        psteros.campaign_surface_energies(entries)


def test_the_unary_oxygen_slab_uses_o2_even_with_an_o_atom_present() -> None:
    # mu_O = E(O2) / 2 = -4.5 eV: gamma = (E - 2 mu_O) / (2 * 8.0) = (-8.0 + 9.0) / 16.0
    assert psteros.campaign_surface_energies([O2, O_ATOM, O_SLAB_UNARY]) == pytest.approx({"o_slab_u": 0.0625})


def test_references_a_slab_does_not_need_are_ignored() -> None:
    # Pd is unfinished and tin is ambiguous, but neither is an element of the unary Au slab.
    unfinished_pd = reference("pd_bulk", {"Pd": 4}, None, state="not started")
    result = psteros.campaign_surface_energies([AU_BULK, AU_SLAB, unfinished_pd, SN_BULK, SN_ALT])
    assert result == pytest.approx({"au_slab": 0.0625})


def test_an_ambiguous_reference_of_a_needed_element_is_still_an_error() -> None:
    with pytest.raises(ValueError) as info:
        psteros.campaign_surface_energies([AU_BULK, AU_ALT, AU_SLAB])
    message = str(info.value)
    assert re.search(r"\bau_bulk\b", message) and re.search(r"\bau_alt\b", message)


def test_reservoirs_and_chemical_potentials_together_are_an_error() -> None:
    with pytest.raises(ValueError, match="reservoirs"):
        psteros.campaign_surface_energies(
            [AU_BULK, AU_SLAB], chemical_potentials_ev={"Au": -3.0}, reservoirs={"Au": "au_bulk"}
        )


def test_a_missing_chemical_potential_is_named_with_slab_and_element() -> None:
    with pytest.raises(ValueError, match="au_slab") as info:
        psteros.campaign_surface_energies([AU_SLAB])  # no Au reference
    assert re.search(r"\bAu\b", str(info.value))


def test_chemical_potentials_given_by_hand_must_cover_every_element() -> None:
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_surface_energies([SLAB_O], chemical_potentials_ev={"Sn": -4.0})
    assert re.search(r"\bO\b", str(info.value))


def test_a_slab_without_energy_is_named_with_its_state() -> None:
    running = slab("slab_o", {"Sn": 4, "O": 8}, None, area=16.0, state="running")
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_surface_energies([running], chemical_potentials_ev=MU_SN_O)
    assert "running" in str(info.value)


def test_a_slab_without_area_is_named_with_its_state() -> None:
    no_area = slab("slab_o", {"Sn": 4, "O": 8}, -40.0, area=None)
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_surface_energies([no_area], chemical_potentials_ev=MU_SN_O)
    assert "finished" in str(info.value)


def test_only_slabs_get_a_surface_energy() -> None:
    assert set(psteros.campaign_surface_energies([SN_BULK, O2, SNO2, SLAB_O], chemical_potentials_ev=MU_SN_O)) == {
        "slab_o"
    }


# =============================================================================
# Sources and signatures
# =============================================================================


@pytest.mark.parametrize("wrap", [list, tuple, as_generator], ids=["list", "tuple", "generator"])
def test_every_function_takes_a_list_a_tuple_or_a_generator(wrap) -> None:
    entries = [SN_BULK, O2, SNO2, SLAB_O]
    assert psteros.campaign_chemical_potentials(wrap(entries)) == pytest.approx(MU_SN_O)
    assert psteros.campaign_references(wrap(entries), host="sno2") == psteros.BinaryOxideReferences(
        -38.0, {"Sn": 2, "O": 4}, -9.0, -4.0
    )
    assert psteros.campaign_surface_energies(wrap(entries), chemical_potentials_ev=MU_SN_O) == pytest.approx(
        {"slab_o": 0.375}
    )


def test_the_analysis_functions_take_keyword_only_options() -> None:
    chemical = inspect.signature(psteros.campaign_chemical_potentials).parameters
    assert list(chemical) == ["source", "reservoirs", "elements"]
    assert all(chemical[name].kind is inspect.Parameter.KEYWORD_ONLY for name in ("reservoirs", "elements"))
    assert chemical["reservoirs"].default is None
    assert chemical["elements"].default is None

    references = inspect.signature(psteros.campaign_references).parameters
    assert list(references) == ["source", "host", "reservoirs", "exclude", "independent"]
    assert all(references[name].kind is inspect.Parameter.KEYWORD_ONLY for name in list(references)[1:])
    assert references["host"].default is inspect.Parameter.empty
    assert references["exclude"].default == ()
    assert references["independent"].default is None
    with pytest.raises(TypeError):
        psteros.campaign_references([SN_BULK, O2, SNO2], "sno2")

    surfaces = inspect.signature(psteros.campaign_surface_energies).parameters
    assert list(surfaces) == ["source", "chemical_potentials_ev", "reservoirs", "surfaces"]
    assert all(surfaces[name].kind is inspect.Parameter.KEYWORD_ONLY for name in list(surfaces)[1:])
    assert surfaces["chemical_potentials_ev"].default is None
    assert surfaces["surfaces"].default == 2


def test_campaign_entries_takes_its_energy_blocks_by_keyword() -> None:
    parameters = inspect.signature(psteros.campaign_entries).parameters
    assert list(parameters) == ["pk", "reference_energy_block", "slab_energy_block"]
    for name in list(parameters)[1:]:
        assert parameters[name].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[name].default is None, name


def test_the_new_names_are_exported() -> None:
    names = (
        "CampaignEntry",
        "campaign_entries",
        "campaign_chemical_potentials",
        "campaign_references",
        "campaign_surface_energies",
    )
    for name in names:
        assert name in psteros.__all__, name
        assert callable(getattr(psteros, name)), name


def test_elements_must_be_a_collection_not_a_string() -> None:
    with pytest.raises(TypeError, match="not a string"):
        psteros.campaign_chemical_potentials([SN_BULK, O2], elements="Sn")
    assert psteros.campaign_chemical_potentials([SN_BULK, O2], elements=("Sn",)) == {"Sn": -32.0 / 8}


@pytest.mark.parametrize(
    "entries, message",
    [
        ([SNO2, O2, reference("h2o", {"H": 2, "O": 1}, -14.0, phase="gas")], "outside"),
        ([SNO2, O2, reference("co", {"C": 1, "O": 1}, -14.0, phase="gas")], "outside"),
        ([SNO2, O2, reference("sn_a", {"Sn": 1}, -4.0), reference("sn_b", {"Sn": 2}, -8.1)], "several references of Sn"),
    ],
)
def test_every_error_of_campaign_references_names_the_host(entries, message) -> None:
    with pytest.raises(ValueError, match=message) as info:
        psteros.campaign_references(entries, host="sno2")
    assert str(info.value).startswith("'sno2'")


def test_an_entry_without_energy_must_say_why() -> None:
    with pytest.raises(ValueError, match="sno2: an entry without energy_ev cannot be in state 'finished'"):
        psteros.CampaignEntry("sno2", "references", "solid", {"Sn": 2, "O": 4}, None)
    entry = psteros.CampaignEntry("sno2", "references", "solid", {"Sn": 2, "O": 4}, None, state="running")
    assert entry.energy_per_atom_ev is None


@pytest.mark.parametrize("field, value", [("energy_ev", "abc"), ("surface_area_angstrom2", "abc")])
def test_a_non_numeric_value_is_named_with_its_label(field, value) -> None:
    options = {"energy_ev": -10.0, "surface_area_angstrom2": 20.0, field: value}
    with pytest.raises(ValueError, match=f"slab_o: {field} must be a finite number"):
        psteros.CampaignEntry("slab_o", "slabs", "solid", {"Sn": 2, "O": 4}, **options)
