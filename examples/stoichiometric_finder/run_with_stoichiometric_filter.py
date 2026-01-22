#!/usr/bin/env python
"""
Metal Surface Energy with Stoichiometric+Symmetric Requirement

This script demonstrates the surface_energy module which ONLY works with
stoichiometric AND symmetric slabs. This is the default behavior - all
slabs are validated BEFORE running expensive DFT calculations.

Material: PdIn (B2 intermetallic)
Surfaces: (110), (100) - pre-filtered for stoichiometry+symmetry

For non-stoichiometric or asymmetric surfaces, use teros.core.thermodynamics.

Usage:
    source ~/envs/aiida/bin/activate
    python run_with_stoichiometric_filter.py
"""

import sys
import os
from aiida import load_profile
from teros.core.surface_energy import (
    build_metal_surface_energy_workgraph,
    analyze_miller_feasibility,
    get_feasibility_summary,
)


def main():
    """Run metal surface energy calculation for stoichiometric+symmetric surfaces."""

    print("\n" + "=" * 70)
    print("METAL SURFACE ENERGY CALCULATION")
    print("Material: PdIn (B2 intermetallic)")
    print("=" * 70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   Profile loaded")

    # Create a simple B2 (CsCl-type) PdIn structure
    print("\n2. Creating PdIn structure...")
    from pymatgen.core import Lattice, Structure
    from aiida import orm

    lattice = Lattice.cubic(3.25)  # Approximate lattice constant
    pdin_pmg = Structure(
        lattice,
        ["Pd", "In"],
        [[0, 0, 0], [0.5, 0.5, 0.5]]
    )

    # Convert to AiiDA StructureData
    from pymatgen.io.ase import AseAtomsAdaptor
    pdin_ase = AseAtomsAdaptor.get_atoms(pdin_pmg)
    pdin_structure = orm.StructureData(ase=pdin_ase)

    print(f"   Structure: {pdin_pmg.composition.reduced_formula}")
    print(f"   Lattice: {pdin_pmg.lattice.a:.3f} A (cubic)")

    # Pre-flight feasibility check
    print("\n3. Pre-flight feasibility analysis...")
    miller_indices = [[1, 1, 0], [1, 0, 0]]

    reports = analyze_miller_feasibility(
        pdin_pmg,
        miller_indices=[(1, 1, 0), (1, 0, 0)],
        min_slab_thickness=15.0,
    )

    for miller_tuple, report in reports.items():
        status = "OK" if report.has_valid_surfaces else "FAIL"
        print(f"   {miller_tuple}: {status} ({report.n_both} valid termination(s))")

    # Check if all orientations are feasible
    all_feasible = all(r.has_valid_surfaces for r in reports.values())
    if not all_feasible:
        print("\n   WARNING: Some orientations have no stoichiometric+symmetric surfaces!")
        print("   The workflow will raise NoStoichiometricSymmetricSurfaceError.")
        print("   Use teros.core.thermodynamics for those orientations instead.")

    # Code configuration (update for your cluster)
    code_label = 'VASP-6.5.1-idefix@obelix'
    potential_family = 'PBE'

    # VASP parameters for metals
    bulk_parameters = {
        'prec': 'Accurate',
        'encut': 500,
        'ediff': 1e-6,
        'ismear': 1,      # Methfessel-Paxton for metals
        'sigma': 0.2,
        'ibrion': 2,
        'isif': 3,        # Full relaxation for bulk
        'nsw': 100,
        'ediffg': -0.01,
        'algo': 'Normal',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }

    slab_parameters = bulk_parameters.copy()
    slab_parameters['isif'] = 2  # Fix cell, relax ions

    # Scheduler options for obelix
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 4,
        },
        'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N PdIn_surf''',
    }

    print("\n4. Building workgraph...")

    # Build workflow
    # Note: require_stoichiometric_symmetric=True is the DEFAULT
    wg = build_metal_surface_energy_workgraph(
        # Structure (pass StructureData directly instead of file path)
        bulk_structure=pdin_structure,

        # Code
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping={'Pd': 'Pd', 'In': 'In'},
        kpoints_spacing=0.02,
        clean_workdir=False,

        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,

        # Slab generation
        miller_indices=miller_indices,
        min_slab_thickness=20.0,
        min_vacuum_thickness=20.0,
        lll_reduce=True,
        center_slab=True,
        primitive=True,

        # Stoichiometric finder configuration
        stoichiometric_strategies=['filter_first', 'symmetrize_check', 'thickness_scan'],
        stoichiometric_max_thickness=30.0,

        # Slab relaxation
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.03,

        # Concurrency
        max_concurrent_jobs=4,

        name='PdIn_surface_energy',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'=' * 70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nModule behavior:")
    print(f"  - Only stoichiometric+symmetric slabs are generated")
    print(f"  - Non-stoichiometric terminations filtered BEFORE DFT")
    print(f"  - Saves computational resources")
    print(f"\nFor non-stoichiometric surfaces:")
    print(f"  Use teros.core.thermodynamics instead")
    print(f"{'=' * 70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
