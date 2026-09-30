#!/usr/bin/env python
"""
Test script to verify slab relaxations work with Zone using context manager paradigm.
"""

from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, Zone, task, spec
from aiida.plugins import WorkflowFactory
from typing import Any

load_profile()

# Get the slabs from the previous run
slabs_node = orm.load_node(4923)  # The get_slabs output from the failed run
slabs_dict = dict(slabs_node.outputs.result.slabs)

print(f"Loaded {len(slabs_dict)} slabs: {list(slabs_dict.keys())}")

# VASP parameters
code = orm.load_code('VASP-VTST-6.4.3@bohr')
slab_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 100,
    "EDIFFG": -0.02,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

slab_options = {
    "resources": {
        "num_machines": 1,
        "num_cores_per_machine": 40,
    },
    "queue_name": "par40",
}

# Build WorkGraph with Zone for parallel slab relaxations
VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
VaspTask = task(VaspWorkChain)

@task
def extract_total_energy(energies):
    """Extract total energy from VASP energies output."""
    if isinstance(energies, orm.Dict):
        energy_dict = energies.get_dict()
    else:
        energy_dict = energies

    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ['energy_extrapolated', 'energy_no_entropy', 'energy']:
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    raise ValueError(f"Could not find energy in energies dict. Available keys: {list(energy_dict.keys())}")

# Create WorkGraph with context manager and Zone
with WorkGraph(name='SlabRelaxWithZone') as wg:
    # Wrap all slab relaxations in a Zone for parallel execution
    with Zone() as slab_zone:
        vasp_results = {}
        energies = {}

        for slab_id, slab_structure in slabs_dict.items():
            # Run VASP relaxation on this slab
            vasp_out = VaspTask(
                structure=slab_structure,
                code=code,
                parameters={'incar': slab_parameters},
                options=slab_options,
                kpoints_spacing=0.3,
                potential_family='PBE',
                potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
                clean_workdir=True,
            )

            # Extract energy
            energy = extract_total_energy(energies=vasp_out.misc)

            vasp_results[slab_id] = vasp_out.structure
            energies[slab_id] = energy.result

    # Set outputs (must be outside the Zone)
    wg.outputs.relaxed_slabs = vasp_results
    wg.outputs.slab_energies = energies

# Export visualization
wg.to_html('test_zone_slabs.html')
print(f"\\n✓ WorkGraph visualization saved to: test_zone_slabs.html")

# Submit
print("\\nSubmitting WorkGraph with Zone for parallel slab relaxations...")
wg.submit(wait=False)

print(f"\\n✓ WorkGraph submitted successfully!")
print(f"WorkGraph PK: {wg.pk}")
print(f"\\nMonitor with: verdi process show {wg.pk}")
