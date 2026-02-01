"""
Hubbard U calculation for NiO (rocksalt) using the Lego module.

Four-stage sequential workflow:
  1. Ionic relaxation (IBRION=2, ISIF=3)
  2. Ground state SCF (LORBIT=11, LWAVE=True, LCHARG=True, no +U)
  3. Response calculations (NSCF + SCF per perturbation potential)
  4. Analysis (linear regression to extract U = 1/chi - 1/chi_0)

NiO is a prototypical correlated oxide where Hubbard U corrections are
essential for reproducing the insulating ground state. Literature values
for U(Ni-3d) are typically 5-8 eV depending on the method and functional.

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/hubbard_u/run_hubbard_u_nio.py
    verdi process show <PK>

Viewing results:
    from teros.core.lego import print_sequential_results
    print_sequential_results(<PK>)
"""

from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

load_profile()

# ── Load structure ──────────────────────────────────────────────────────
struct_path = Path(__file__).parent / 'nio.vasp'
structure = orm.StructureData(ase=read(str(struct_path)))

# ── Common INCAR parameters ────────────────────────────────────────────
base_incar = {
    'encut': 520,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'prec': 'Accurate',
    'algo': 'Normal',
    'nelm': 200,
}

# ── Define stages ───────────────────────────────────────────────────────
stages = [
    # Stage 1: Full relaxation (ions + cell)
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {
            **base_incar,
            'ibrion': 2,
            'nsw': 100,
            'isif': 3,
            'lreal': 'Auto',
            'lwave': False,
            'lcharg': False,
        },
        'restart': None,
        'retrieve': ['CONTCAR', 'OUTCAR'],
    },
    # Stage 2: Ground state SCF (produces WAVECAR, CHGCAR, OUTCAR)
    # Prerequisites for hubbard_response: lorbit=11, lwave=True,
    # lcharg=True, and retrieve=['OUTCAR']
    {
        'name': 'ground_state',
        'type': 'vasp',
        'structure_from': 'relax',
        'incar': {
            **base_incar,
            'nsw': 0,
            'ibrion': -1,
            'ldau': False,
            'lmaxmix': 4,
            'lorbit': 11,
            'lwave': True,
            'lcharg': True,
        },
        'restart': None,
        'retrieve': ['OUTCAR'],
    },
    # Stage 3: Response calculations (NSCF + SCF per potential V)
    # Runs 4 pairs of calculations with perturbation potentials
    {
        'name': 'response',
        'type': 'hubbard_response',
        'ground_state_from': 'ground_state',
        'structure_from': 'relax',
        'target_species': 'Ni',
        'potential_values': [-0.2, -0.1, 0.1, 0.2],
        'ldaul': 2,  # d-electrons
        'incar': base_incar,
    },
    # Stage 4: Linear regression to extract U
    # No VASP calculations — pure data analysis
    {
        'name': 'analysis',
        'type': 'hubbard_analysis',
        'response_from': 'response',
        'structure_from': 'relax',
        'target_species': 'Ni',
        'ldaul': 2,
    },
]

# ── Submit ──────────────────────────────────────────────────────────────
from teros.core.lego import quick_vasp_sequential, print_sequential_results

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.5.1@localwork',
    kpoints_spacing=0.04,
    potential_family='PBE',
    potential_mapping={'Ni': 'Ni_pv', 'O': 'O'},
    options={
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8,
        },
    },
    name='nio_hubbard_u',
)

pk = result['__workgraph_pk__']
print(f"Submitted WorkGraph PK: {pk}")
print(f"Stages: {result['__stage_names__']}")
print()
print("Monitor with:")
print(f"  verdi process show {pk}")
print(f"  verdi process report {pk}")
print()
print("View results when finished:")
print("  from teros.core.lego import print_sequential_results")
print(f"  print_sequential_results({pk})")
