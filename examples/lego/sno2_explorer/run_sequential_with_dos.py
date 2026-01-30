"""
Sequential VASP workflow with DOS: SnO2 relaxation + DOS calculation.

This example demonstrates the new DOS stage support in quick_vasp_sequential().

Stages:
  1. relax - Relax primitive SnO2 cell
  2. dos   - Calculate DOS on relaxed structure (SCF + DOS)

Monitor:  verdi process show <PK>
Results:  print_sequential_results(result)
"""

from aiida import orm, load_profile
from ase.io import read
from pathlib import Path

load_profile('presto')

from teros.core.lego import quick_vasp_sequential, print_sequential_results

# Load structure
structure_file = Path(__file__).parent / 'sno2.vasp'
structure = orm.StructureData(ase=read(str(structure_file), format='vasp'))

# Obelix cluster configuration
code_label = 'VASP-6.5.1-idefix-4@obelix'
potential_family = 'PBE'
potential_mapping = {'Sn': 'Sn_d', 'O': 'O'}

options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,  # PROCESS_MPI=4 (hybrid MPI+OpenMP)
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N sno2_dos''',
}

# Light INCAR for testing (6-atom cell)
incar_relax = {
    'nsw': 50,
    'ibrion': 2,
    'isif': 2,
    'ediff': 1e-4,
    'encut': 400,
    'prec': 'Normal',
    'ismear': 0,
    'sigma': 0.05,
    'algo': 'Normal',
    'lwave': True,
    'lcharg': True,
}

# Stage definitions
stages = [
    # Stage 1: VASP relaxation
    {
        'name': 'relax',
        'type': 'vasp',  # Optional, default is 'vasp'
        'incar': incar_relax,
        'restart': None,
        'kpoints_spacing': 0.06,
        'retrieve': ['CONTCAR', 'OUTCAR'],
    },
    # Stage 2: DOS calculation on relaxed structure
    {
        'name': 'dos',
        'type': 'dos',
        'structure_from': 'relax',  # Use relaxed structure
        'scf_incar': {
            # SCF step: generate charge density
            # Note: BandsWorkChain handles lwave/lcharg internally
            'encut': 400,
            'ediff': 1e-5,
            'ismear': 0,
            'sigma': 0.05,
            'prec': 'Normal',
            'algo': 'Normal',
            # These are forced by the module:
            'nsw': 0,        # Static calculation
            'ibrion': -1,    # No ionic relaxation
        },
        'dos_incar': {
            # DOS step: non-SCF calculation on denser k-mesh
            # Note: BandsWorkChain handles ICHARG internally for non-SCF DOS
            'encut': 400,
            'prec': 'Normal',
            'nedos': 2000,   # Number of DOS points
            'lorbit': 11,    # Projected DOS (lm-decomposed)
            'ismear': -5,    # Tetrahedron method for accurate DOS
            # These are forced by the module:
            'nsw': 0,        # Static calculation
            'ibrion': -1,    # No ionic relaxation
        },
        'kpoints_spacing': 0.06,
        'dos_kpoints_spacing': 0.04,  # Denser k-mesh for DOS
        'retrieve': ['DOSCAR'],
    },
]

if __name__ == '__main__':
    result = quick_vasp_sequential(
        structure=structure,
        stages=stages,
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=options,
        name='sno2_sequential_with_dos',
    )

    print(f"\nWorkGraph PK: {result['__workgraph_pk__']}")
    print(f"Stage names: {result['__stage_names__']}")
    print(f"Stage types: {result['__stage_types__']}")
    print(f"\nMonitor: verdi process show {result['__workgraph_pk__']}")
    print("\nAfter completion, view results with:")
    print("  from teros.core.lego import print_sequential_results, get_stage_results")
    print(f"  result = {{'__workgraph_pk__': {result['__workgraph_pk__']}, '__stage_names__': {result['__stage_names__']}, '__stage_types__': {result['__stage_types__']}}}")
    print("  print_sequential_results(result)")
    print("  dos_result = get_stage_results(result, 'dos')")
    print("  print(dos_result)")
