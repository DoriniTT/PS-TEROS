"""
Hubbard U calculation example: SnO2 bulk (rutile).

Calculates the Hubbard U parameter for Sn d-electrons using the linear
response method (Cococcioni & de Gironcoli).

The workflow performs:
  1. Ground state SCF (no +U, LORBIT=11)
  2. For each potential V in [-0.2, -0.1, 0.1, 0.2]:
     - Non-SCF response (ICHARG=11, LDAUTYPE=3)
     - SCF response (LDAUTYPE=3)
  3. Linear regression to extract U = 1/chi - 1/chi_0

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/hubbard_u/run_hubbard_u_sno2.py
    verdi process show <PK>

Note:
    SnO2 is not a typical system where Hubbard U is applied (Sn d-electrons
    are deep core states). This is a demonstration example. For realistic
    applications, use transition metal oxides like NiO, Fe2O3, or MnO.
"""

from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

load_profile()

# ── Load structure ──────────────────────────────────────────────────────
struct_path = Path(__file__).parent / 'sno2.vasp'
structure = orm.StructureData(ase=read(str(struct_path)))

# ── Submit Hubbard U calculation ────────────────────────────────────────
from teros.core.lego import quick_hubbard_u

result = quick_hubbard_u(
    structure=structure,
    code_label='VASP-6.5.1@localwork',
    target_species='Sn',
    incar={
        'encut': 400,
        'ediff': 1e-5,
        'ismear': 0,
        'sigma': 0.05,
        'prec': 'Accurate',
        'algo': 'Normal',
        'nelm': 100,
    },
    potential_values=[-0.2, -0.1, 0.1, 0.2],
    ldaul=2,  # d-electrons
    kpoints_spacing=0.05,
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options={
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
        },
    },
    name='sno2_hubbard_u',
)

pk = result['__workgraph_pk__']
print(f"Submitted WorkGraph PK: {pk}")
print()
print("Monitor with:")
print(f"  verdi process show {pk}")
print(f"  verdi process report {pk}")
print()
print("Get results when done:")
print("  from teros.core.lego import get_results")
print(f"  results = get_results({pk})")
