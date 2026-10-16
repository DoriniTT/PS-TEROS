"""Write the static ternary diagram CSV with whichever psteros is on PYTHONPATH (check T7 of analyse_ternary.py).

    echo <json> | PYTHONPATH=<checkout> python legacy_ternary.py out.csv

Only functions that existed before the vibrational contributions are called.
"""

import json
import sys

import psteros

data = json.load(sys.stdin)
energies = data["energies"]
bulk = data["bulk"]["srtio3_bulk"]["composition"]


def per_atom(label, element):
    return energies[label] / data["bulk"][label]["composition"][element]


references = psteros.TernaryOxideReferences(
    bulk_energy_ev=energies["srtio3_bulk"],
    bulk_composition=bulk,
    oxygen_molecule_energy_ev=energies["o2"],
    element_energies_per_atom_ev={"Sr": per_atom("sr_metal", "Sr"), "Ti": per_atom("ti_metal", "Ti")},
    competing_phases=(
        psteros.CompetingPhase("SrO", energies["sro_bulk"], data["bulk"]["sro_bulk"]["composition"]),
        psteros.CompetingPhase("TiO2", energies["tio2_bulk"], data["bulk"]["tio2_bulk"]["composition"]),
    ),
)
terminations = [
    psteros.SlabTermination(label, energies[label], values["composition"], values["area"])
    for label, values in data["slabs"].items()
]
psteros.ternary_surface_phase_diagram(terminations, references).to_csv(sys.argv[1])
