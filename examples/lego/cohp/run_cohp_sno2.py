"""
COHP analysis example: SnO2 bulk (rutile).

Two-stage workflow:
  1. SCF with ISYM=-1, LWAVE=True to produce WAVECAR and other files
  2. COHP analysis using LOBSTER to compute bonding properties

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/cohp/run_cohp_sno2.py
    verdi process show <PK>

Requirements:
    - lobster binary installed at ~/.local/bin/lobster or in PATH
    - VASP calculation with ISYM=-1, LWAVE=True
"""

from pathlib import Path

from aiida import load_profile, orm
from pymatgen.core import Structure

load_profile()

# ── Load structure ──────────────────────────────────────────────────────
# Use SnO2 rutile structure
struct_path = Path(__file__).parent.parent / 'bader' / 'sno2.vasp'
if not struct_path.exists():
    # Create a simple SnO2 structure if not available
    from pymatgen.core import Lattice, Structure
    lattice = Lattice.tetragonal(4.737, 3.186)
    species = ['Sn', 'Sn', 'O', 'O', 'O', 'O']
    coords = [
        [0.0, 0.0, 0.0],      # Sn
        [0.5, 0.5, 0.5],      # Sn
        [0.306, 0.306, 0.0],  # O
        [0.694, 0.694, 0.0],  # O
        [0.194, 0.806, 0.5],  # O
        [0.806, 0.194, 0.5],  # O
    ]
    pmg_struct = Structure(lattice, species, coords)
else:
    pmg_struct = Structure.from_file(str(struct_path))

structure = orm.StructureData(pymatgen=pmg_struct)

# ── Define stages ───────────────────────────────────────────────────────
stages = [
    {
        'name': 'scf',
        'type': 'vasp',
        'incar': {
            'encut': 520,
            'ediff': 1e-6,
            'ismear': 0,
            'sigma': 0.05,
            'ibrion': -1,
            'nsw': 0,
            'isym': -1,      # Required for LOBSTER
            'lwave': True,   # Required for LOBSTER
            'prec': 'Accurate',
            'lreal': 'Auto',
            'lcharg': False,
        },
        'restart': None,
        'kpoints_spacing': 0.02,  # Dense k-mesh recommended for COHP
        'retrieve': ['WAVECAR', 'CONTCAR', 'POTCAR', 'INCAR', 'DOSCAR'],
    },
    {
        'name': 'cohp',
        'type': 'cohp',
        'cohp_from': 'scf',
        'cohp_start_energy': -15.0,
        'cohp_end_energy': 10.0,
        'basis_set': 'pbeVaspFit2015',
        'bond_cutoff': 3.0,
        'save_projection': True,
        'compute_coop': True,
    },
]

# ── Submit workflow ─────────────────────────────────────────────────────
from teros.core.lego import quick_vasp_sequential

# For testing on localwork (adjust for your cluster)
result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.5.1@localwork',
    kpoints_spacing=0.02,
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options={
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
        },
    },
    name='sno2_cohp',
)

print(f"Submitted WorkGraph PK: {result['__workgraph_pk__']}")
print(f"Stages: {result['__stage_names__']}")
print()
print("Monitor with:")
print(f"  verdi process show {result['__workgraph_pk__']}")
print(f"  verdi process report {result['__workgraph_pk__']}")
print()
print("Get results when done:")
print("  from teros.core.lego import print_sequential_results")
print(f"  print_sequential_results({result})")
