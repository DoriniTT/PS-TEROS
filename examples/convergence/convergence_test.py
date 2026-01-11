#!/usr/bin/env python
"""
Example: VASP convergence testing with PS-TEROS.

This script demonstrates how to use the convergence module to determine
optimal ENCUT and k-points spacing for a given structure.
"""

from aiida import orm, load_profile
from teros.core.convergence import build_convergence_workgraph, get_convergence_results


def submit_convergence_test():
    """Submit a convergence test workflow."""
    load_profile()

    # Create simple Si structure (diamond cubic, 2 atoms)
    from ase.build import bulk
    from aiida.orm import StructureData
    ase_structure = bulk('Si', 'diamond', a=5.43)
    structure = StructureData(ase=ase_structure)

    # Lightweight VASP parameters for local testing (4 processors)
    # Note: aiida-vasp requires lowercase INCAR keys
    builder_inputs = {
        'parameters': {
            'incar': {
                'prec': 'Normal',      # Not Accurate to save time
                'ismear': 0,
                'sigma': 0.1,
                'ediff': 1e-5,         # Looser convergence
                'lreal': 'Auto',
                'lwave': False,
                'lcharg': False,
                'nelm': 60,
                'algo': 'Fast',
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
            'max_wallclock_seconds': 1800,  # 30 min max
        },
        'kpoints_spacing': 0.1,  # Coarse starting value
        'potential_family': 'PBE',
        'potential_mapping': {'Si': 'Si'},
        'clean_workdir': True,
    }

    # Small convergence scan for testing
    conv_settings = {
        'cutoff_start': 200,
        'cutoff_stop': 350,
        'cutoff_step': 50,
        'kspacing_start': 0.12,
        'kspacing_stop': 0.06,
        'kspacing_step': -0.02,
        'cutoff_kconv': 300,
        'kspacing_cutconv': 0.08,
    }

    wg = build_convergence_workgraph(
        structure=structure,
        code_label='vasp-6.5.1-std@localhost',
        builder_inputs=builder_inputs,
        conv_settings=conv_settings,
        convergence_threshold=0.002,  # 2 meV/atom (relaxed for testing)
        name='Si_convergence_test',
    )

    wg.submit(wait=False)
    print(f"Submitted convergence test: PK {wg.pk}")
    print(f"Monitor with: verdi process show {wg.pk}")

    return wg


def get_results(pk: int):
    """Get results from a completed convergence test."""
    from aiida import load_profile
    from aiida_workgraph import WorkGraph
    load_profile()

    wg = WorkGraph.from_pk(pk)
    results = get_convergence_results(wg)

    print("\n" + "=" * 60)
    print("CONVERGENCE TEST RESULTS")
    print("=" * 60)

    if results['recommended_cutoff']:
        print(f"\nRecommended ENCUT: {results['recommended_cutoff']} eV")
    else:
        print("\nENCUT: NOT CONVERGED")

    if results['recommended_kspacing']:
        print(f"Recommended k-spacing: {results['recommended_kspacing']} A^-1")
    else:
        print("K-spacing: NOT CONVERGED")

    if results['convergence_summary']:
        summary = results['convergence_summary']
        print(f"\nThreshold used: {summary['threshold_used'] * 1000:.1f} meV/atom")

    print("=" * 60)
    return results


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        pk = int(sys.argv[1])
        get_results(pk)
    else:
        submit_convergence_test()
