#!/usr/bin/env python
"""Test script for Map-based slab relaxations."""

from aiida import load_profile
from teros.workgraph_map import build_core_workgraph_with_map

print("Loading AiiDA profile...")
load_profile()

# Define structures directory
structures_dir = "/home/thiagotd/git/PS-TEROS/teros/structures"

# Define calculation parameters
code_label = "VASP-VTST-6.4.3@bohr"
potential_family = "PBE"

# BULK RELAXATION PARAMETERS
bulk_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

bulk_options = {
    "resources": {
        "num_machines": 1,
        "num_cores_per_machine": 40,
    },
    "queue_name": "par40",
}

# Reference parameters (simplified - same for all)
ref_parameters = bulk_parameters.copy()
ref_options = bulk_options.copy()

# SLAB PARAMETERS
slab_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,  # Relax atoms only
    "NSW": 100,
    "EDIFFG": -0.02,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

slab_options = bulk_options.copy()

# Slab generation parameters
miller_indices = [1, 0, 0]  # (100) surface
min_slab_thickness = 10.0
min_vacuum_thickness = 15.0

print(f"\n{'='*80}")
print(f"BUILDING WORKGRAPH WITH MAP ZONE")
print(f"{'='*80}\n")

wg = build_core_workgraph_with_map(
    structures_dir=structures_dir,
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label=code_label,
    potential_family=potential_family,
    bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"P": "P"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters=bulk_parameters,
    bulk_options=bulk_options,
    metal_parameters=ref_parameters,
    metal_options=ref_options,
    nonmetal_parameters=ref_parameters,
    nonmetal_options=ref_options,
    oxygen_parameters=ref_parameters,
    oxygen_options=ref_options,
    clean_workdir=True,
    miller_indices=miller_indices,
    min_slab_thickness=min_slab_thickness,
    min_vacuum_thickness=min_vacuum_thickness,
    slab_parameters=slab_parameters,
    slab_options=slab_options,
    slab_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    slab_kpoints_spacing=0.3,
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
    in_unit_planes=False,
    max_normal_search=None,
    relax_slabs=True,
    name="Ag3PO4_SlabsRelax_Map_100",
)

# Export visualization
try:
    html_file = f"ag3po4_slabs_relax_map_100.html"
    wg.to_html(html_file)
    print(f"✓ WorkGraph visualization saved to: {html_file}\n")
except Exception as e:
    print(f"✗ Could not generate HTML visualization: {e}\n")

# Submit
print(f"{'='*80}")
print("Submitting WorkGraph...")
print(f"{'='*80}\n")
wg.submit(wait=False)

print(f"✓ WorkGraph submitted successfully!")
print(f"WorkGraph PK: {wg.pk}\n")
print(f"Monitor with: verdi process show {wg.pk}")
