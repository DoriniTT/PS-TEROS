#!/usr/bin/env python
"""Example driver for the Map-zone slab relaxation workflow."""

from __future__ import annotations

import sys
from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

# Ensure repository root is on ``sys.path`` when executed as a script
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

from test_modules.zone_approach.workgraph import build_zone_workgraph


def load_bulk_structure(structure_path: Path) -> orm.StructureData:
    """Create a ``StructureData`` directly from a structure file."""
    atoms = read(structure_path)
    return orm.StructureData(ase=atoms)


def main() -> None:
    """Configure parameters and launch the Map-zone slab relaxation workflow."""
    load_profile()

    structures_dir = REPO_ROOT / 'structures'
    bulk_path = structures_dir / 'ag3po4.cif'

    if not bulk_path.exists():
        raise FileNotFoundError(f"Bulk structure file not found: {bulk_path}")

    bulk_structure = load_bulk_structure(bulk_path)

    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'
    potential_mapping = {'Ag': 'Ag', 'P': 'P', 'O': 'O'}

    slab_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,
        'NSW': 100,
        'EDIFFG': -0.02,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    slab_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    workflow = build_zone_workgraph(
        bulk_structure=bulk_structure,
        miller_indices=(1, 0, 0),
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=slab_parameters,
        options=slab_options,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    results = workflow.run()

    process = getattr(workflow, 'process', None)
    generated_slabs = None

    if hasattr(workflow, 'outputs') and hasattr(workflow.outputs, 'generated_slabs'):
        generated_slabs = workflow.outputs.generated_slabs.value
    elif isinstance(results, dict):
        generated_slabs = results.get('generated_slabs')

    print('\nWorkflow finished:')
    if process is not None:
        print(f'  PK: {process.pk}')
    else:
        print('  PK: unavailable')

    if isinstance(generated_slabs, dict):
        print('  Generated slabs:', sorted(generated_slabs.keys()))
    else:
        print('  Generated slabs: unavailable')


if __name__ == '__main__':
    main()
