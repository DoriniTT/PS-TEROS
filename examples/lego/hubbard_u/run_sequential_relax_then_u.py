"""
Sequential workflow: Relax SnO2 bulk, then calculate Hubbard U.

Two-stage sequential workflow:
  1. Ionic relaxation (IBRION=2, NSW=50)
  2. Hubbard U calculation on the relaxed structure

This demonstrates the hubbard_u brick used as a stage in
quick_vasp_sequential, chained after a relaxation stage.

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/hubbard_u/run_sequential_relax_then_u.py
    verdi process show <PK>

Note:
    SnO2 is used here for demonstration. For realistic Hubbard U
    calculations, use transition metal oxides (NiO, Fe2O3, MnO, etc.).
"""

from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

load_profile()

# ── Load structure ──────────────────────────────────────────────────────
struct_path = Path(__file__).parent / 'sno2.vasp'
structure = orm.StructureData(ase=read(str(struct_path)))

# ── Define stages ───────────────────────────────────────────────────────
stages = [
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {
            'encut': 400,
            'ediff': 1e-5,
            'ismear': 0,
            'sigma': 0.05,
            'ibrion': 2,
            'nsw': 50,
            'isif': 3,
            'prec': 'Accurate',
            'lreal': 'Auto',
            'lwave': False,
            'lcharg': False,
        },
        'restart': None,
        'kpoints_spacing': 0.05,
        'retrieve': ['CONTCAR', 'OUTCAR'],
    },
    {
        'name': 'hubbard_u',
        'type': 'hubbard_u',
        'structure_from': 'relax',
        'target_species': 'Sn',
        'incar': {
            'encut': 400,
            'ediff': 1e-5,
            'ismear': 0,
            'sigma': 0.05,
            'prec': 'Accurate',
            'algo': 'Normal',
            'nelm': 100,
        },
        'potential_values': [-0.2, -0.1, 0.1, 0.2],
        'ldaul': 2,
        'kpoints_spacing': 0.05,
    },
]

# ── Submit sequential workflow ──────────────────────────────────────────
from teros.core.lego import quick_vasp_sequential, print_sequential_results

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.5.1@localwork',
    kpoints_spacing=0.05,
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options={
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
        },
    },
    name='sno2_relax_then_hubbard_u',
)

pk = result['__workgraph_pk__']
print(f"Submitted WorkGraph PK: {pk}")
print(f"Stages: {result['__stage_names__']}")
print()
print("Monitor with:")
print(f"  verdi process show {pk}")
print(f"  verdi process report {pk}")
print()
print("Get results when done:")
print("  from teros.core.lego import print_sequential_results, get_stage_results")
print(f"  print_sequential_results({result})")
print(f"  u_result = get_stage_results({result}, 'hubbard_u')")
print("  print(f\"U = {u_result['hubbard_u_eV']:.3f} eV\")")
