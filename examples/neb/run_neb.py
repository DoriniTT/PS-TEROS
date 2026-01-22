#!/usr/bin/env python
"""Example: NEB calculation for vacancy migration or surface diffusion.

This script demonstrates how to use the PS-TEROS NEB module to calculate
activation barriers for atomic transitions.

Prerequisites:
    - AiiDA profile configured with VASP code
    - VTST-compiled VASP (for IOPT optimizer support)
    - Initial and final structures prepared (e.g., from manual setup or MD)

Usage:
    1. Modify the structure loading section to use your structures
    2. Adjust builder_inputs for your cluster configuration
    3. Run: python run_neb.py
"""

from aiida import load_profile, orm

# Load AiiDA profile
load_profile()

from teros.core.neb import (
    build_neb_workgraph,
    get_neb_results,
    print_neb_summary,
    estimate_neb_resources,
    validate_neb_structures,
)

# =============================================================================
# Configuration
# =============================================================================

# Code label (must be VTST-compiled VASP)
# Use appropriate label for your cluster:
# - 'VASP-6.4.3@bohr' for bohr cluster
# - 'VASP-6.5.1-idefix@obelix' for obelix cluster
CODE_LABEL = 'VASP-6.4.3@bohr'

# Number of intermediate NEB images (typically 3-7)
# More images = more accurate path but higher cost
N_IMAGES = 5

# Builder inputs for VASP calculations
# Adjust for your cluster and system
BUILDER_INPUTS = {
    'parameters': {
        'incar': {
            # Electronic settings
            'encut': 520,        # Plane wave cutoff (eV)
            'ediff': 1e-6,       # Electronic convergence (eV)
            'ismear': 0,         # Gaussian smearing (for insulators)
            'sigma': 0.05,       # Smearing width (eV)

            # Precision and accuracy
            'prec': 'Accurate',
            'algo': 'Normal',

            # Output control
            'lwave': False,      # Don't write WAVECAR (saves disk)
            'lcharg': False,     # Don't write CHGCAR
        }
    },
    'kpoints_spacing': 0.03,    # K-point spacing (Å⁻¹)
    'potential_family': 'PBE',
    'potential_mapping': {
        # Map elements to POTCAR types
        # Adjust for your system
        'Ag': 'Ag',
        'O': 'O',
        # 'Cu': 'Cu_pv',  # Example: use _pv for transition metals
    },
    'options': {
        'resources': {
            'num_machines': 3,             # 3 nodes for parallel NEB
            'num_cores_per_machine': 40,   # 40 cores per node on bohr
        },
        'queue_name': 'par120',            # Use par120 for multi-node jobs
    },
}

# NEB-specific options
NEB_OPTIONS = {
    'relax_endpoints': True,       # Relax initial/final before NEB
    'interpolation_method': 'idpp',  # 'idpp' (recommended) or 'linear'
    'climb': True,                 # Two-stage NEB→CI-NEB for accurate TS
    'spring_constant': -5.0,       # Variable spring (negative = nudged)
    'neb_optimizer': 1,            # 1=LBFGS (robust), 3=FIRE (sometimes faster)
    'force_convergence': -0.05,    # Force convergence (eV/Å)
    'max_steps': 500,              # Max ionic steps per stage
}


# =============================================================================
# Load Structures
# =============================================================================

def load_structures():
    """
    Load or create initial and final structures for NEB.

    Modify this function for your specific use case:
    - Load from files (CIF, POSCAR, etc.)
    - Load from AiiDA database (by PK)
    - Create programmatically

    Returns:
        tuple: (initial_structure, final_structure) as StructureData
    """
    from ase.io import read

    # Option 1: Load from files
    # initial = orm.StructureData(ase=read('POSCAR_initial'))
    # final = orm.StructureData(ase=read('POSCAR_final'))

    # Option 2: Load from AiiDA database by PK
    # initial = orm.load_node(12345)
    # final = orm.load_node(12346)

    # Option 3: Create example structures (for testing)
    # This creates a simple Cu surface with a vacancy migration
    from ase.build import fcc111, add_adsorbate
    from ase.constraints import FixAtoms
    import numpy as np

    # Create Cu(111) slab
    slab_initial = fcc111('Cu', size=(3, 3, 4), vacuum=10.0)

    # Create vacancy by removing one atom (initial state)
    # Remove atom at position (1, 1) in the surface layer
    del slab_initial[8]  # Remove atom 8 (adjust based on your structure)

    # Create final state by moving a neighboring atom into vacancy
    slab_final = fcc111('Cu', size=(3, 3, 4), vacuum=10.0)
    del slab_final[9]  # Remove different atom to represent final state

    # Convert to AiiDA StructureData
    initial = orm.StructureData(ase=slab_initial)
    final = orm.StructureData(ase=slab_final)

    return initial, final


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Run NEB calculation."""
    print("=" * 60)
    print("PS-TEROS NEB Calculation")
    print("=" * 60)

    # Load structures
    print("\nLoading structures...")
    initial, final = load_structures()

    # Validate structures are compatible
    print("Validating structures...")
    try:
        validate_neb_structures(initial, final)
        print("  Structures are compatible for NEB")
    except ValueError as e:
        print(f"  ERROR: {e}")
        return

    # Print structure info
    print(f"\nInitial structure: {initial.get_formula()}")
    print(f"  Atoms: {len(initial.sites)}")
    print(f"Final structure: {final.get_formula()}")
    print(f"  Atoms: {len(final.sites)}")

    # Estimate resources
    print("\nResource estimate:")
    estimate = estimate_neb_resources(
        n_images=N_IMAGES,
        n_atoms=len(initial.sites),
    )
    print(f"  Estimated walltime: {estimate['estimated_walltime_hours']:.1f} hours")
    print(f"  Estimated nodes: {estimate['estimated_nodes']}")
    for note in estimate['notes']:
        print(f"  Note: {note}")

    # Build WorkGraph
    print(f"\nBuilding NEB WorkGraph...")
    print(f"  N images: {N_IMAGES}")
    print(f"  Code: {CODE_LABEL}")
    print(f"  Relax endpoints: {NEB_OPTIONS['relax_endpoints']}")
    print(f"  Climbing image: {NEB_OPTIONS['climb']}")

    wg = build_neb_workgraph(
        initial_structure=initial,
        final_structure=final,
        n_images=N_IMAGES,
        code_label=CODE_LABEL,
        builder_inputs=BUILDER_INPUTS,
        **NEB_OPTIONS,
        name='NEB_vacancy_migration',
    )

    # Submit workflow
    print("\nSubmitting WorkGraph...")
    wg.submit(wait=False)

    print(f"\nSubmitted NEB WorkGraph: PK={wg.pk}")
    print("\nMonitor progress with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")

    print("\nAfter completion, analyze results with:")
    print(f"  from teros.core.neb import print_neb_summary, get_neb_results")
    print(f"  print_neb_summary({wg.pk})")
    print(f"  results = get_neb_results({wg.pk})")

    return wg.pk


if __name__ == '__main__':
    pk = main()
