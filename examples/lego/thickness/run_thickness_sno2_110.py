"""
Slab thickness convergence example: SnO2 (110) surface.

Four-stage workflow:
  1. Bulk relaxation (full cell + ions, ISIF=3)
  2. Slab generation: create slabs at 3, 5, 7, 9, 11 layers
  3. Slab relaxation: relax all slabs in parallel (ions only, ISIF=2)
  4. Thickness analysis: compute surface energies and determine convergence

Usage:
    source ~/envs/aiida/bin/activate
    verdi daemon restart
    python examples/lego/thickness/run_thickness_sno2_110.py
    verdi process show <PK>

After completion:
    from teros.core.lego import print_sequential_results
    print_sequential_results(result)

    # Or extract results programmatically:
    from teros.core.lego import get_stage_results
    thickness_results = get_stage_results(result, 'thickness_110')
    print(f"Converged: {thickness_results['converged']}")
    print(f"Recommended layers: {thickness_results['recommended_layers']}")
    for n, gamma in thickness_results['surface_energies'].items():
        print(f"  {n} layers: {gamma:.4f} J/m^2")
"""

from pathlib import Path

from aiida import load_profile, orm
from pymatgen.core import Structure

load_profile()

# ── Load bulk structure ───────────────────────────────────────────────
struct_path = Path(__file__).parent / 'sno2.vasp'
pmg_struct = Structure.from_file(str(struct_path))
structure = orm.StructureData(pymatgen=pmg_struct)

# ── VASP configuration ───────────────────────────────────────────────
CODE_LABEL = 'VASP-6.5.1@localwork'
POTENTIAL_FAMILY = 'PBE'
POTENTIAL_MAPPING = {'Sn': 'Sn_d', 'O': 'O'}
OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 8,
    },
}

# ── Shared INCAR settings ────────────────────────────────────────────
SLAB_INCAR = {
    'encut': 520,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'ibrion': 2,
    'nsw': 100,
    'isif': 2,           # Ions only (fixed cell for slab)
    'prec': 'Accurate',
    'lreal': 'Auto',
    'lwave': False,
    'lcharg': False,
}

# ── Define stages ─────────────────────────────────────────────────────
stages = [
    # Stage 1: Bulk relaxation (full cell optimization)
    {
        'name': 'bulk_relax',
        'type': 'vasp',
        'incar': {
            'encut': 520,
            'ediff': 1e-6,
            'ismear': 0,
            'sigma': 0.05,
            'ibrion': 2,
            'nsw': 100,
            'isif': 3,       # Full relaxation (ions + cell)
            'prec': 'Accurate',
            'lreal': 'Auto',
            'lwave': False,
            'lcharg': False,
        },
        'restart': None,
        'kpoints_spacing': 0.03,
    },
    # Stage 2: Generate slabs at multiple thicknesses
    {
        'name': 'gen_slabs',
        'type': 'slab_gen',
        'structure_from': 'bulk_relax',
        'miller_indices': [1, 1, 0],
        'layer_counts': [3, 5, 7, 9, 11],
        'min_vacuum_thickness': 15.0,     # Angstroms
        'termination_index': 0,           # First termination
    },
    # Stage 3: Relax all slabs in parallel
    {
        'name': 'relax_slabs',
        'type': 'batch',
        'structures_from': 'gen_slabs',
        'base_incar': SLAB_INCAR,
        'kpoints_spacing': 0.03,
        'max_concurrent_jobs': 4,
    },
    # Stage 4: Analyze thickness convergence
    {
        'name': 'thickness_110',
        'type': 'thickness',
        'relaxed_from': 'relax_slabs',
        'bulk_from': 'bulk_relax',
        'miller_indices': [1, 1, 0],
        'convergence_threshold': 0.01,    # J/m^2
    },
]

# ── Submit workflow ───────────────────────────────────────────────────
from teros.core.lego import quick_vasp_sequential

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label=CODE_LABEL,
    kpoints_spacing=0.03,
    potential_family=POTENTIAL_FAMILY,
    potential_mapping=POTENTIAL_MAPPING,
    options=OPTIONS,
    name='sno2_110_thickness',
)

print(f"Submitted WorkGraph PK: {result['__workgraph_pk__']}")
print(f"Stages: {result['__stage_names__']}")
print(f"Stage types: {result['__stage_types__']}")
print()
print("Monitor with:")
print(f"  verdi process show {result['__workgraph_pk__']}")
print(f"  verdi process report {result['__workgraph_pk__']}")
print()
print("Get results when done:")
print("  from teros.core.lego import print_sequential_results")
print(f"  print_sequential_results({result})")
