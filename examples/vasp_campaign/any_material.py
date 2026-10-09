"""The general campaign readers on energies typed by hand: unary, binary, ternary, intermetallic.

The same functions take the PK of a finished campaign graph instead of the
entries (``psteros.campaign_entries(pk)`` gives these entries from the graph).
The energies below are illustrative numbers, not results.

    python any_material.py
"""

from __future__ import annotations

import psteros
from psteros import CampaignEntry


def reference(label, phase, composition, energy_ev):
    return CampaignEntry(label=label, group="references", phase=phase, composition=composition, energy_ev=energy_ev)


def slab(label, composition, energy_ev, area_angstrom2):
    return CampaignEntry(
        label=label, group="slabs", phase="solid", composition=composition, energy_ev=energy_ev,
        surface_area_angstrom2=area_angstrom2,
    )


def unary_gold() -> None:
    entries = [
        reference("au", "solid", {"Au": 4}, -12.80),
        slab("au_111", {"Au": 12}, -36.30, 7.20),
        slab("au_100", {"Au": 12}, -35.90, 8.32),
    ]
    for label, gamma in psteros.campaign_surface_energies(entries).items():
        print(f"Au   {label}: {gamma * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2:.3f} J/m2")


def binary_tin_oxide() -> None:
    entries = [
        reference("sno2", "solid", {"Sn": 2, "O": 4}, -41.80),
        reference("sn", "solid", {"Sn": 8}, -30.40),
        reference("o2", "gas", {"O": 2}, -9.86),
        slab("slab_o", {"Sn": 6, "O": 14}, -120.10, 21.90),
        slab("slab_sno", {"Sn": 6, "O": 12}, -116.30, 21.90),
        slab("slab_sn2o", {"Sn": 6, "O": 10}, -110.20, 21.90),
    ]
    oxide = psteros.campaign_references(entries, host="sno2")
    terminations = [
        psteros.SlabTermination(entry.label, entry.energy_ev, entry.composition, entry.surface_area_angstrom2)
        for entry in entries
        if entry.group == "slabs"
    ]
    diagram = psteros.surface_phase_diagram(terminations, oxide)
    stable = list(dict.fromkeys(diagram.stable))  # in order of increasing Delta mu_O
    print(f"SnO2 stable terminations: {stable}, transitions (Delta mu_O, from, to): {diagram.transitions}")


def ternary_strontium_titanate() -> None:
    entries = [
        reference("srtio3", "solid", {"Sr": 1, "Ti": 1, "O": 3}, -40.00),
        reference("sr", "solid", {"Sr": 1}, -1.60),
        reference("ti", "solid", {"Ti": 2}, -15.50),
        reference("o2", "gas", {"O": 2}, -9.86),
        reference("sro", "solid", {"Sr": 1, "O": 1}, -12.00),  # a competing phase, found automatically
        reference("tio2", "solid", {"Ti": 2, "O": 4}, -53.50),
    ]
    oxide = psteros.campaign_references(entries, host="srtio3")
    print(f"SrTiO3 competing phases: {[phase.label for phase in oxide.competing_phases]}, "
          f"formation enthalpy {oxide.formation_enthalpy_ev:.3f} eV per formula unit")


def intermetallic_palladium_indium() -> None:
    entries = [
        reference("pdin", "solid", {"Pd": 1, "In": 1}, -9.40),
        reference("pd", "solid", {"Pd": 1}, -5.20),
        reference("indium", "solid", {"In": 2}, -5.00),  # "in" is a Python keyword, so not a label
        slab("pdin_110", {"Pd": 6, "In": 6}, -55.50, 13.40),
    ]
    # No phase-diagram model for intermetallics yet: gamma at the Pd-rich and In-rich limits
    # (equal for this stoichiometric slab, different for a slab that is not stoichiometric).
    mu = psteros.campaign_chemical_potentials(entries)  # {"Pd": E/atom of Pd, "In": E/atom of In}
    bulk = next(entry for entry in entries if entry.label == "pdin").energy_ev
    pd_rich = {"Pd": mu["Pd"], "In": bulk - mu["Pd"]}
    in_rich = {"In": mu["In"], "Pd": bulk - mu["In"]}
    for name, potentials in (("Pd-rich", pd_rich), ("In-rich", in_rich)):
        gamma = psteros.campaign_surface_energies(entries, chemical_potentials_ev=potentials)["pdin_110"]
        print(f"PdIn pdin_110 ({name}): {gamma * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2:.3f} J/m2")


if __name__ == "__main__":
    unary_gold()
    binary_tin_oxide()
    ternary_strontium_titanate()
    intermetallic_palladium_indium()
