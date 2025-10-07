#!/usr/bin/env python
"""Example driver for the pythonic scatter-gather slab relaxation workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

# Ensure repository root is on ``sys.path`` when executed as a script
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

from test_modules.pythonic.workgraph import (
    build_mock_scatter_gather_workgraph,
    build_pythonic_workgraph,
)
from test_modules.pythonic.aiat_ternary import (
    create_mock_bulk_energy,
    create_mock_formation_enthalpy,
    create_mock_reference_energies,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Return parsed command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--mock',
        action='store_true',
        help='run a lightweight scatter-gather workflow that does not require VASP.',
    )
    parser.add_argument(
        '--mock-count',
        type=int,
        default=3,
        help='number of mock items to generate when ``--mock`` is set.',
    )
    parser.add_argument(
        '--mock-delta',
        type=float,
        default=0.5,
        help='offset applied to mock values when ``--mock`` is set.',
    )
    parser.add_argument(
        '--with-thermodynamics',
        action='store_true',
        help='compute ab initio atomistic thermodynamics (surface energies).',
    )
    parser.add_argument(
        '--sampling',
        type=int,
        default=100,
        help='number of sampling points for chemical potential grid.',
    )
    return parser.parse_args(argv)


def load_bulk_structure(structure_path: Path) -> orm.StructureData:
    """Create a ``StructureData`` directly from a structure file."""
    atoms = read(structure_path)
    return orm.StructureData(ase=atoms)


def _print_thermodynamics_results(title: str, thermo_namespace) -> None:
    """Pretty-print thermodynamics results."""
    print(title)
    
    # Try direct dictionary-like access (AttributeDict after workflow.run())
    if hasattr(thermo_namespace, 'keys'):
        keys = list(thermo_namespace.keys())
        for key in sorted(keys):
            if not key.startswith('_'):
                thermo_data = thermo_namespace[key]
                if hasattr(thermo_data, 'get_dict'):
                    data = thermo_data.get_dict()
                    print(f'    {key}:')
                    phi = data.get('phi')
                    gamma_m = data.get('Gamma_M_vs_Nref')
                    gamma_o = data.get('Gamma_O_vs_Nref')
                    area = data.get('area_A2')
                    
                    if phi is not None:
                        print(f'      φ (reference surface energy): {phi:.6f} eV/Ų')
                    if gamma_m is not None:
                        print(f'      Γ_M (surface excess M): {gamma_m:.6f} atoms/Ų')
                    if gamma_o is not None:
                        print(f'      Γ_O (surface excess O): {gamma_o:.6f} atoms/Ų')
                    if area is not None:
                        print(f'      Surface area: {area:.2f} Ų')
        return
    
    # Try TaskSocketNamespace access (during building)
    if hasattr(thermo_namespace, '_sockets'):
        for key in sorted(k for k in thermo_namespace._sockets.keys() if not k.startswith('_')):
            thermo_data = getattr(thermo_namespace, key)
            if hasattr(thermo_data, 'get_dict'):
                data = thermo_data.get_dict()
                print(f'    {key}:')
                phi = data.get('phi')
                gamma_m = data.get('Gamma_M_vs_Nref')
                gamma_o = data.get('Gamma_O_vs_Nref')
                area = data.get('area_A2')
                
                if phi is not None:
                    print(f'      φ (reference surface energy): {phi:.6f} eV/Ų')
                if gamma_m is not None:
                    print(f'      Γ_M (surface excess M): {gamma_m:.6f} atoms/Ų')
                if gamma_o is not None:
                    print(f'      Γ_O (surface excess O): {gamma_o:.6f} atoms/Ų')
                if area is not None:
                    print(f'      Surface area: {area:.2f} Ų')
        return
    
    print('    (Unable to display thermodynamics results)')


def _print_scalar_namespace(title: str, namespace: orm.Dict | dict) -> None:
    """Pretty-print a namespace of scalar values."""
    print(title)
    
    # If it's a TaskSocketNamespace, iterate over its sockets
    if hasattr(namespace, '_sockets'):
        for key in sorted(namespace._sockets.keys()):
            if not key.startswith('_'):
                socket = getattr(namespace, key)
                if hasattr(socket, 'value'):
                    print(f'  {key}: {socket.value}')
                else:
                    print(f'  {key}: {socket}')
        return
    
    # Otherwise handle as dict
    if hasattr(namespace, 'value'):
        data = namespace.value
    else:
        data = namespace

    for key, value in sorted(data.items()):
        if hasattr(value, 'value'):
            numeric = value.value
        else:
            numeric = value
        print(f'  {key}: {numeric}')


    """Pretty-print a namespace of scalar values."""
    print(title)
    
    # If it's a TaskSocketNamespace, iterate over its sockets
    if hasattr(namespace, '_sockets'):
        for key in sorted(namespace._sockets.keys()):
            if not key.startswith('_'):
                socket = getattr(namespace, key)
                if hasattr(socket, 'value'):
                    print(f'  {key}: {socket.value}')
                else:
                    print(f'  {key}: {socket}')
        return
    
    # Otherwise handle as dict
    if hasattr(namespace, 'value'):
        data = namespace.value
    else:
        data = namespace

    for key, value in sorted(data.items()):
        if hasattr(value, 'value'):
            numeric = value.value
        else:
            numeric = value
        print(f'  {key}: {numeric}')


def main(argv: list[str] | None = None) -> None:
    """Configure parameters and launch the pythonic scatter-gather slab relaxation workflow."""
    args = parse_args(argv)
    load_profile()

    if args.mock:
        workflow = build_mock_scatter_gather_workgraph(
            count=args.mock_count,
            delta=args.mock_delta,
            with_thermodynamics=args.with_thermodynamics,
        )
        workflow.run()

        print('\nMock scatter-gather workflow finished:')
        _print_scalar_namespace('  Source values:', workflow.outputs.source_values)
        _print_scalar_namespace('  Shifted values:', workflow.outputs.shifted_values)
        
        # Display thermodynamics if computed
        if args.with_thermodynamics:
            # Try to get thermo_results from process outputs (after run) or workflow outputs (during build)
            thermo_ns = None
            if hasattr(workflow, 'process') and hasattr(workflow.process, 'outputs'):
                if hasattr(workflow.process.outputs, 'thermo_results'):
                    thermo_ns = workflow.process.outputs.thermo_results
            elif hasattr(workflow, 'outputs') and hasattr(workflow.outputs, 'thermo_results'):
                thermo_ns = workflow.outputs.thermo_results
            
            if thermo_ns is not None:
                print('  Thermodynamics results:')
                
                # Get keys - check _sockets FIRST (TaskSocketNamespace), then keys (AttributeDict)
                if hasattr(thermo_ns, '_sockets'):
                    keys = [k for k in thermo_ns._sockets.keys() if not k.startswith('_')]
                elif hasattr(thermo_ns, 'keys'):
                    keys = [k for k in thermo_ns.keys() if not k.startswith('_')]
                else:
                    keys = []
                
                for key in sorted(keys):
                    try:
                        if hasattr(thermo_ns, '_sockets'):
                            thermo_data = getattr(thermo_ns, key)
                        else:
                            thermo_data = thermo_ns[key]
                        
                        if hasattr(thermo_data, 'get_dict'):
                            data = thermo_data.get_dict()
                            print(f'    {key}:')
                            phi = data.get('phi')
                            if phi is not None:
                                print(f'      φ: {phi:.6f}')
                            gamma_m = data.get('Gamma_M_vs_Nref')
                            if gamma_m is not None:
                                print(f'      Γ_M: {gamma_m:.6f} atoms/Ų')
                    except Exception as e:
                        print(f'    {key}: Error - {e}')
        return

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

    # Prepare thermodynamics inputs if requested
    thermodynamics_kwargs = {}
    if args.with_thermodynamics:
        print('\nPreparing thermodynamics inputs (using mock data)...')
        # Create mock inputs as plain AiiDA data nodes (not task outputs)
        # We need to call the functions directly, not as tasks
        from aiida_workgraph import task
        
        # Get the actual Python functions (not decorated versions)
        bulk_energy = orm.Float(-100.0)  # Mock energy
        reference_energies = orm.Dict(dict={
            'ag_energy_per_atom': -2.7,
            'p_energy_per_atom': -5.0,
            'o_energy_per_atom': -4.5,
        })
        formation_enthalpy = orm.Float(-0.5)
        
        thermodynamics_kwargs = {
            'compute_thermodynamics': True,
            'bulk_energy': bulk_energy,
            'reference_energies': reference_energies,
            'formation_enthalpy': formation_enthalpy,
            'sampling': args.sampling,
        }
        print('  Mock bulk energy, reference energies, and formation enthalpy created.')

    workflow = build_pythonic_workgraph(
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
        **thermodynamics_kwargs,
    )

    results = workflow.run()

    process = getattr(workflow, 'process', None)
    generated_slabs = None

    if hasattr(workflow, 'outputs') and hasattr(workflow.outputs, 'generated_slabs'):
        # Handle TaskSocketNamespace
        if hasattr(workflow.outputs.generated_slabs, '_sockets'):
            generated_slabs = {
                key: getattr(workflow.outputs.generated_slabs, key)
                for key in workflow.outputs.generated_slabs._sockets.keys()
                if not key.startswith('_')
            }
        elif hasattr(workflow.outputs.generated_slabs, 'value'):
            generated_slabs = workflow.outputs.generated_slabs.value
    elif isinstance(results, dict):
        generated_slabs = results.get('generated_slabs')

    print('\nWorkflow finished:')
    if process is not None:
        print(f'  PK: {process.pk}')
    else:
        print('  PK: unavailable')

    if isinstance(generated_slabs, dict):
        print(f'  Generated slabs: {len(generated_slabs)} terminations')
        print(f'    Labels: {sorted(generated_slabs.keys())}')
    else:
        print('  Generated slabs: unavailable')
    
    # Print thermodynamics results if computed
    if args.with_thermodynamics:
        # Try to get thermo_results from process outputs (after run) or workflow outputs (during build)
        thermo_ns = None
        if hasattr(workflow, 'process') and hasattr(workflow.process, 'outputs'):
            if hasattr(workflow.process.outputs, 'surface_energies'):
                thermo_ns = workflow.process.outputs.surface_energies
        elif hasattr(workflow, 'outputs') and hasattr(workflow.outputs, 'surface_energies'):
            thermo_ns = workflow.outputs.surface_energies
        
        if thermo_ns is not None:
            _print_thermodynamics_results('\n  Surface energy calculations:', thermo_ns)
            print(f'\n  Use "verdi process show {process.pk if process else "PK"}" for full results')



if __name__ == '__main__':
    main(sys.argv[1:])

