#!/usr/bin/env python
"""
Fukui+ (Nucleophilic Attack) Calculation Example

This example demonstrates running Fukui function calculations for SnO2 (110)
surface using the interpolation method.

The workflow:
1. Runs 4 parallel static VASP calculations at different charge states
   (delta_N = 0.0, 0.05, 0.10, 0.15)
2. Collects CHGCAR files into a single FolderData output
3. The CHGCAR files can then be processed with FukuiGrid.py

Usage:
    source ~/envs/aiida/bin/activate
    python run_fukui_plus.py

Requirements:
    - AiiDA profile configured
    - VASP code registered
    - POTCAR files available
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.fukui import (
    build_fukui_workgraph,
    get_fukui_results,
    print_fukui_summary,
    calculate_nelect,
    print_nelect_breakdown,
)


def main():
    """Run Fukui+ calculation example."""

    print("\n" + "=" * 70)
    print("FUKUI+ CALCULATION - SnO2 (110) Surface")
    print("=" * 70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')  # MODIFY: use your profile name
    print("   Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   Daemon is running")

    # Load structure
    print("\n2. Loading SnO2 (110) structure...")
    structure_file = Path(__file__).parent / "sno2_110_paper.vasp"
    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        structure = orm.StructureData(ase=atoms)
        print(f"   Loaded from file: {structure_file}")
    else:
        print(f"   ERROR: Structure file not found: {structure_file}")
        return 1

    print(f"   Composition: {structure.get_composition()}")
    print(f"   Number of atoms: {len(structure.sites)}")

    # Define POTCAR family and mapping (used for both NELECT and calculations)
    potential_family = 'PBE'
    potential_mapping = {
        'Sn': 'Sn_d',  # Sn_d has 14 valence electrons
        'O': 'O',      # O has 6 valence electrons
    }

    # Calculate NELECT automatically from POTCARs
    print("\n3. Calculating NELECT from POTCAR files...")
    nelect_neutral = calculate_nelect(
        structure=structure,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
    )
    print_nelect_breakdown(structure, potential_family, potential_mapping)

    # Define builder inputs
    print("\n4. Defining VASP calculation inputs...")

    # Code label for obelix cluster
    code_label = 'VASP-6.5.1-idefix@obelix'

    builder_inputs = {
        'parameters': {
            'incar': {
                # Precision and cutoff
                'prec': 'Accurate',
                'encut': 500,       # Use well-converged cutoff for charge density

                # SCF convergence - robust settings for fractional charges
                'algo': 'All',      # Try all algorithms (most robust)
                'nelm': 300,        # Allow many SCF steps for fractional charges
                'ediff': 1e-6,      # Tight convergence for accurate charge density

                # Charge density initialization
                'icharg': 2,        # Start from superposition of atomic charges

                # Mixing parameters for better convergence
                'amix': 0.2,        # Linear mixing parameter (default 0.4)
                'bmix': 0.0001,     # Kerker mixing parameter
                'amix_mag': 0.8,    # Magnetic mixing for spin-polarized
                'bmix_mag': 0.0001,
                'lmaxmix': 4,       # For d-electrons (Sn has d-electrons)

                # Smearing - Gaussian for insulators
                'ismear': 0,
                'sigma': 0.05,

                # Real-space projection
                'lreal': 'Auto',

                # Spin settings - spin-polarized
                'ispin': 2,

                # Parallelization for obelix (4 cores per node in idefix-4)
                'ncore': 2,
                'kpar': 1,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
            'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N Fukui_SnO2_110_plus''',
        },
        'kpoints_spacing': 0.03,              # Good k-points density for charge density
        'potential_family': potential_family,
        'potential_mapping': potential_mapping,
        'clean_workdir': False,               # Keep files for debugging
    }

    # Delta N values (standard interpolation points)
    delta_n_values = [0.0, 0.05, 0.10, 0.15]

    print(f"   Code: {code_label}")
    print(f"   ENCUT: {builder_inputs['parameters']['incar']['encut']} eV")
    print(f"   k-points spacing: {builder_inputs['kpoints_spacing']} A^-1")
    print(f"   Delta N values: {delta_n_values}")

    # Build WorkGraph
    print("\n5. Building Fukui WorkGraph...")
    wg = build_fukui_workgraph(
        structure=structure,
        nelect_neutral=nelect_neutral,
        delta_n_values=delta_n_values,
        code_label=code_label,
        builder_inputs=builder_inputs,
        relax_first=False,        # Structure is already relaxed
        # Atom fixing (only used when relax_first=True):
        # fix_atoms=True,         # Enable selective dynamics
        # fix_type='center',      # Fix center atoms, relax surfaces
        # fix_thickness=5.0,      # 5 Angstroms from center
        fukui_type='plus',
        max_concurrent_jobs=4,    # obelix can handle multiple jobs
        name='Fukui_SnO2_110_plus',
    )
    print(f"   WorkGraph created: {wg.name}")

    # Print expected calculations
    print("\n6. Expected VASP calculations:")
    print("-" * 50)
    for delta_n in delta_n_values:
        nelect = nelect_neutral - delta_n  # Fukui+: remove electrons
        print(f"   delta_N = {delta_n:.2f} -> NELECT = {nelect:.2f}")

    # Submit
    print("\n7. Submitting WorkGraph...")
    wg.submit(wait=False)
    print(f"   Submitted!")
    print(f"   PK: {wg.pk}")

    # Save PK to file for reference
    pk_file = Path(__file__).parent / 'last_pk.txt'
    with open(pk_file, 'w') as f:
        f.write(f"{wg.pk}\n")
    print(f"   PK saved to: {pk_file}")

    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print(f"""
Monitor progress:
    verdi process show {wg.pk}
    verdi process report {wg.pk}

After completion, extract results:
    >>> from teros.core.fukui import get_fukui_results, print_fukui_summary
    >>> results = get_fukui_results({wg.pk})
    >>> print(results['file_names'])
    >>> print_fukui_summary({wg.pk})

The CHGCAR files will be in results['chgcar_folder']:
    >>> chgcar_folder = results['chgcar_folder']
    >>> chgcar_folder.list_object_names()
    ['CHGCAR_0.00', 'CHGCAR_0.05', 'CHGCAR_0.10', 'CHGCAR_0.15']

Export CHGCAR files for FukuiGrid.py:
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as tmpdir:
    ...     for fname in chgcar_folder.list_object_names():
    ...         content = chgcar_folder.get_object_content(fname)
    ...         with open(f"{{tmpdir}}/{{fname}}", 'wb' if isinstance(content, bytes) else 'w') as f:
    ...             f.write(content)
""")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
