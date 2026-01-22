#!/usr/bin/env python
"""
Example: VASP ENCUT and k-points convergence testing for Si on bohr cluster.

This script demonstrates how to use the convergence module to determine
optimal ENCUT and k-points spacing for a silicon structure using the
bohr cluster with par40 queue.

Usage:
    # Submit the convergence test
    python convergence_si_bohr.py

    # Get results after completion (replace PK with your workflow PK)
    python convergence_si_bohr.py 12345

    # Plot results
    python convergence_si_bohr.py 12345 --plot

Cluster Configuration:
    - Computer: bohr
    - Queue: par40 (40 cores, 1 node, max 3 days)
    - VASP code: VASP-6.4.3@bohr
    - POTCAR family: PBE
"""

import sys
from aiida import orm, load_profile
from teros.core.convergence import (
    build_convergence_workgraph,
    get_convergence_results,
    print_convergence_summary,
    plot_convergence,
    export_convergence_data,
)


def create_si_structure():
    """Create Si diamond structure (2 atoms)."""
    from ase.build import bulk
    ase_structure = bulk('Si', 'diamond', a=5.431)  # Si lattice constant
    return orm.StructureData(ase=ase_structure)


def submit_convergence_test():
    """Submit a convergence test workflow for Si on bohr cluster."""
    load_profile()

    # Check daemon is running
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("ERROR: AiiDA daemon is not running!")
        print("Start it with: verdi daemon start")
        return None

    print("\n" + "=" * 70)
    print("      VASP CONVERGENCE TEST - Si on bohr (par40)")
    print("=" * 70)

    # Create Si structure
    structure = create_si_structure()
    print(f"Structure: {structure.get_formula()} ({len(structure.sites)} atoms)")

    # =====================================================
    # BOHR CLUSTER CONFIGURATION (par40 queue)
    # =====================================================
    # - par40: up to 40 cores (1 node), max walltime 3 days
    # - PBS scheduler
    # =====================================================

    builder_inputs = {
        'parameters': {
            'incar': {
                # Basic SCF settings (no relaxation for convergence test)
                'prec': 'Accurate',
                'ismear': 0,       # Gaussian smearing (good for semiconductors)
                'sigma': 0.05,
                'ediff': 1e-6,     # Tight convergence for accurate energy
                'lreal': 'Auto',
                'lwave': False,    # Don't write WAVECAR
                'lcharg': False,   # Don't write CHGCAR
                'nelm': 100,       # Max SCF iterations
                'algo': 'Normal',  # Standard SCF algorithm
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 40,  # par40 queue
            },
            'queue_name': 'par40',
            'max_wallclock_seconds': 3 * 24 * 3600,  # 3 days max
        },
        'kpoints_spacing': 0.1,  # Starting value (will be overridden in scan)
        'potential_family': 'PBE',
        'potential_mapping': {'Si': 'Si'},
        'clean_workdir': True,   # Clean remote folders after completion
    }

    # =====================================================
    # CONVERGENCE SCAN SETTINGS
    # =====================================================
    # - ENCUT: 200-600 eV in 50 eV steps (9 calculations)
    # - k-spacing: 0.08-0.02 Å⁻¹ in -0.01 steps (7 calculations)
    # - Total: ~16 single-point calculations
    # =====================================================

    conv_settings = {
        # ENCUT convergence scan
        'cutoff_start': 200,    # eV - starting ENCUT
        'cutoff_stop': 600,     # eV - ending ENCUT
        'cutoff_step': 50,      # eV - step size

        # K-points convergence scan
        'kspacing_start': 0.08,   # Å⁻¹ - coarsest grid
        'kspacing_stop': 0.02,    # Å⁻¹ - finest grid
        'kspacing_step': -0.01,   # Å⁻¹ - step (negative = getting finer)

        # Fixed values for the other scan
        'cutoff_kconv': 450,      # eV - ENCUT used during k-points scan
        'kspacing_cutconv': 0.03, # Å⁻¹ - k-spacing used during ENCUT scan
    }

    print("\nConvergence scan settings:")
    print(f"  ENCUT: {conv_settings['cutoff_start']} - {conv_settings['cutoff_stop']} eV "
          f"(step: {conv_settings['cutoff_step']} eV)")
    print(f"  k-spacing: {conv_settings['kspacing_start']} - {conv_settings['kspacing_stop']} Å⁻¹ "
          f"(step: {conv_settings['kspacing_step']} Å⁻¹)")
    print(f"  Threshold: 1.0 meV/atom")

    # Build WorkGraph
    wg = build_convergence_workgraph(
        structure=structure,
        code_label='VASP-6.4.3@bohr',
        builder_inputs=builder_inputs,
        conv_settings=conv_settings,
        convergence_threshold=0.001,  # 1 meV/atom
        name='Si_convergence_bohr',
    )

    # Submit
    wg.submit(wait=False)

    print("\n" + "-" * 70)
    print(f"Submitted! WorkGraph PK: {wg.pk}")
    print("-" * 70)
    print("\nMonitor progress with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")
    print("\nAfter completion, get results with:")
    print(f"  python {__file__} {wg.pk}")
    print(f"  python {__file__} {wg.pk} --plot")
    print("=" * 70 + "\n")

    return wg


def show_results(pk: int, do_plot: bool = False, export_dir: str = None):
    """Show results from a completed convergence test."""
    load_profile()

    print(f"\nLoading WorkGraph PK: {pk}")

    # Print formatted summary
    print_convergence_summary(pk)

    # Plot if requested
    if do_plot:
        print("Generating convergence plot...")
        # plot_convergence handles both saving and interactive display
        fig = plot_convergence(pk, save_path=f'convergence_{pk}.png', show=True)

    # Export data if requested
    if export_dir:
        print(f"\nExporting data to: {export_dir}")
        files = export_convergence_data(pk, export_dir, prefix=f'si_conv_{pk}')
        print("Created files:")
        for file_type, path in files.items():
            print(f"  {file_type}: {path}")

    # Return raw results for programmatic use
    from aiida_workgraph import WorkGraph
    wg = WorkGraph.load(pk)
    return get_convergence_results(wg)


def main():
    """Main entry point."""
    if len(sys.argv) > 1:
        # Get results for existing workflow
        pk = int(sys.argv[1])
        do_plot = '--plot' in sys.argv
        export_dir = None

        for i, arg in enumerate(sys.argv):
            if arg == '--export' and i + 1 < len(sys.argv):
                export_dir = sys.argv[i + 1]

        show_results(pk, do_plot=do_plot, export_dir=export_dir)
    else:
        # Submit new convergence test
        submit_convergence_test()


if __name__ == '__main__':
    main()
