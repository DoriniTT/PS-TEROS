#!/home/thiagotd/envs/aiida/bin/python
"""
Intermetallic Surface Energy Calculation for PdIn (B2 structure)

This script demonstrates the surface_energy module for computing
surface energies of stoichiometric intermetallics using the simple formula:
γ = (E_slab - N·E_bulk/atom) / (2A)

Material: PdIn (B2 / CsCl-type structure)
Surfaces: (110), (100), (111)

This formula is valid for stoichiometric and symmetric surfaces where:
- The slab composition matches the bulk stoichiometry
- Both surfaces (top and bottom) are equivalent

Usage:
    source ~/envs/aiida/bin/activate
    python surface_thermo_intermetallics.py
"""

import sys
import os
from aiida import load_profile
from teros.core.surface_energy import build_metal_surface_energy_workgraph


def main():
    """Run surface energy calculation for PdIn intermetallic."""

    print("\n" + "="*70)
    print("INTERMETALLIC SURFACE ENERGY CALCULATION")
    print("Material: PdIn (B2 structure)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bulk_structure_path = os.path.join(script_dir, 'pdin.cif')

    print(f"\n2. Structure:")
    print(f"   Bulk:   {bulk_structure_path}")

    # Code configuration - using cluster02
    code_label = 'VASP-6.5.1@cluster02'
    potential_family = 'PBE'

    # VASP parameters for intermetallics
    # Note: ISMEAR=1 (Methfessel-Paxton) is recommended for metals/intermetallics
    bulk_parameters = {
        'prec': 'Accurate',
        'encut': 500,
        'ediff': 1e-6,
        'ismear': 1,      # Methfessel-Paxton for metals
        'sigma': 0.2,     # Smearing width
        'ibrion': 2,
        'isif': 3,        # Full relaxation for bulk
        'nsw': 100,
        'ediffg': -0.01,
        'algo': 'Normal',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }

    # Slab parameters - only relax ionic positions
    slab_parameters = bulk_parameters.copy()
    slab_parameters['isif'] = 2  # Fix cell, relax ions

    # Scheduler options for cluster02 (24 cores per machine)
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    }

    # Surface orientations to study - all in ONE WorkGraph!
    # For B2 structure, (110) is typically the most stable
    miller_indices = [
        [1, 1, 0],  # Most stable for B2
        [1, 0, 0],  # Second most stable
        [1, 1, 1],  # Third
    ]

    print(f"\n3. Miller indices: {miller_indices}")
    print("   All orientations in a single WorkGraph!")

    print("\n4. Building workgraph...")

    # Build ONE workflow containing all Miller indices
    wg = build_metal_surface_energy_workgraph(
        # Structure
        bulk_structure_path=bulk_structure_path,

        # Code
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping={'Pd': 'Pd', 'In': 'In'},  # Map both elements!
        kpoints_spacing=0.03,
        clean_workdir=False,

        # Bulk parameters
        bulk_parameters=bulk_parameters,
        bulk_options=common_options,

        # Slab generation - ALL Miller indices!
        miller_indices=miller_indices,
        min_slab_thickness=20.0,
        min_vacuum_thickness=20.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,  # Important for stoichiometric surfaces!
        primitive=True,

        # Slab relaxation
        slab_parameters=slab_parameters,
        slab_options=common_options,
        slab_kpoints_spacing=0.04,

        # Concurrency control - cluster02 can only run one calculation at a time
        max_concurrent_jobs=1,

        name='PdIn_surface_energy',
    )

    print("   WorkGraph built successfully")
    print(f"   Contains tasks: surface_hkl_110, surface_hkl_100, surface_hkl_111")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs per orientation:")
    print(f"  - bulk_energy, bulk_energy_per_atom, bulk_structure")
    print(f"  - slab_structures, relaxed_slabs, slab_energies")
    print(f"  - surface_energies (gamma_eV_A2, gamma_J_m2)")
    print(f"\nNew intermetallic-specific outputs:")
    print(f"  - compound_type: 'intermetallic'")
    print(f"  - formula: 'InPd'")
    print(f"  - elements: ['In', 'Pd']")
    print(f"  - composition: bulk composition dict")
    print(f"  - slab_composition: slab composition dict")
    print(f"  - is_stoichiometric: True/False")
    print(f"\nExpected surface energies for PdIn (B2):")
    print(f"  (110): ~0.5-1.0 J/m2 (typically most stable for B2)")
    print(f"  (100): ~1.0-1.5 J/m2")
    print(f"  (111): ~1.5-2.5 J/m2")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
