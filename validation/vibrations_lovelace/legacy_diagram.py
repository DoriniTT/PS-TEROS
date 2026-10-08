"""Write the 0 K diagram CSV from static energies with whichever psteros is on PYTHONPATH.

Used by ``analyse.py`` (check 9) to run the code of the commit before the vibrations and the current code on the
same numbers: ``echo <json> | PYTHONPATH=<checkout> python legacy_diagram.py out.csv``. Only functions that
existed before the vibrational contributions are called.
"""

import json
import sys

import psteros

data = json.load(sys.stdin)
energies = data["energies"]
references = psteros.BinaryOxideReferences(
    bulk_energy_ev=energies["sno2_bulk"],
    bulk_composition=data["bulk_composition"],
    oxygen_molecule_energy_ev=energies["o2"],
    metal_energy_per_atom_ev=energies["alpha_sn"] / data["metal_atoms"],
)
terminations = [
    psteros.SlabTermination(label, energies[label], values["composition"], values["area"])
    for label, values in data["slabs"].items()
]
psteros.surface_phase_diagram(terminations, references).to_csv(sys.argv[1])
