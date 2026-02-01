"""
Convergence parameter test: SnO2 bulk (rutile).

Finds optimal ENCUT and k-points spacing using vasp.v2.converge.

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/convergence/run_convergence_sno2.py
    verdi process show <PK>
"""
from pathlib import Path

from aiida import load_profile, orm
from pymatgen.core import Structure

load_profile()

# ── Load structure ──────────────────────────────────────────────────────
struct_path = Path(__file__).parent / 'sno2.vasp'
pmg_struct = Structure.from_file(str(struct_path))
structure = orm.StructureData(pymatgen=pmg_struct)

# ── Define stages ───────────────────────────────────────────────────────
stages = [
    {
        'name': 'conv_test',
        'type': 'convergence',
        'incar': {
            'prec': 'Accurate',
            'ismear': 0,
            'sigma': 0.05,
            'ediff': 1e-6,
            'ibrion': -1,
            'nsw': 0,
            'lreal': 'Auto',
            'lwave': False,
            'lcharg': False,
        },
        'conv_settings': {
            'cutoff_start': 300,
            'cutoff_stop': 700,
            'cutoff_step': 50,
            'kspacing_start': 0.06,
            'kspacing_stop': 0.02,
            'kspacing_step': -0.005,
            'cutoff_kconv': 520,
            'kspacing_cutconv': 0.03,
        },
        'convergence_threshold': 0.001,  # 1 meV/atom
    },
]

# ── Cluster config: bohr, par40 queue ───────────────────────────────────
from teros.core.lego import quick_vasp_sequential

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.4.3@bohr',
    kpoints_spacing=0.03,
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options={
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'custom_scheduler_commands': '#PBS -q par40\n#PBS -j oe\n#PBS -N sno2_conv',
    },
    name='sno2_convergence',
)

print(f"Submitted WorkGraph PK: {result['__workgraph_pk__']}")
print()
print("Monitor with:")
print(f"  verdi process show {result['__workgraph_pk__']}")
print(f"  verdi process report {result['__workgraph_pk__']}")
print()
print("Get results when done:")
print("  from teros.core.lego import print_sequential_results")
print(f"  print_sequential_results({result})")
