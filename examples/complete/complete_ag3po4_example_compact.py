#!/home/thiagotd/envs/aiida/bin/python
"""
PS-TEROS Complete Workflow Example - Ag3PO4
Testing: bulk relaxation, formation enthalpy, slab generation, relaxation energy,
         cleavage energy, and surface thermodynamics.

Usage: source ~/envs/aiida/bin/activate && python complete_ag3po4_example_compact.py
"""

import sys
import os
from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph


# =============================================================================
# CONFIGURATION - Modify parameters here
# =============================================================================

# --- File paths and structure filenames ---
STRUCTURES_DIR = 'structures'  # Relative to this script
BULK_FILE = 'ag3po4.cif'
METAL_FILE = 'Ag.cif'
NONMETAL_FILE = 'P.cif'
OXYGEN_FILE = 'O2.cif'

# Predefined slab structures (optional - leave empty to auto-generate)
SLABS_DIR = "input_structures/ag3po4"  # Relative to this script
SLAB_FILES = [
    "slab_term_0_ter.cif",
    "slab_term_1_ter.cif",
]

# --- Computational resources ---
CODE_LABEL = 'VASP-VTST-6.4.3@bohr'
POTENTIAL_FAMILY = 'PBE'
QUEUE_NAME = 'par40'
NUM_CORES = 40

# --- Slab generation parameters ---
MILLER_INDICES = [1, 0, 0]  # (100) surface
MIN_SLAB_THICKNESS = 15.0   # Angstroms
MIN_VACUUM_THICKNESS = 15.0  # Angstroms

# --- Workflow flags ---
RELAX_SLABS = True
COMPUTE_RELAXATION_ENERGY = True  # E_relaxed - E_unrelaxed
COMPUTE_CLEAVAGE = True
COMPUTE_THERMODYNAMICS = True
THERMODYNAMICS_SAMPLING = 50  # Grid points for μ sampling

# --- Other settings ---
KPOINTS_SPACING = 0.4
CLEAN_WORKDIR = False
WORKFLOW_NAME = 'Ag3PO4_Complete_Workflow'


# =============================================================================
# VASP PARAMETERS - Organized by calculation type
# =============================================================================

# Common parameters for all calculations
COMMON_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ALGO': 'Normal',
    'LWAVE': False,
    'LCHARG': False,
}

# --- Bulk relaxation (Ag3PO4) ---
BULK_PARAMS = {
    **COMMON_PARAMS,
    'ISMEAR': 0,      # Gaussian smearing (semiconductors)
    'SIGMA': 0.05,
    'IBRION': 2,      # Conjugate gradient
    'ISIF': 3,        # Relax cell + atoms
    'NSW': 100,
    'EDIFFG': -0.01,  # Force convergence
    'LREAL': 'Auto',
}

BULK_POTENTIAL_MAP = {'Ag': 'Ag', 'P': 'P', 'O': 'O'}

# --- Metal reference (Ag) ---
METAL_PARAMS = {
    **COMMON_PARAMS,
    'ISMEAR': 1,      # Methfessel-Paxton (metals)
    'SIGMA': 0.2,
    'IBRION': 2,
    'ISIF': 3,
    'NSW': 100,
    'EDIFFG': -0.01,
    'LREAL': 'Auto',
}

METAL_POTENTIAL_MAP = {'Ag': 'Ag'}

# --- Nonmetal reference (P) ---
NONMETAL_PARAMS = {
    **COMMON_PARAMS,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 3,
    'NSW': 100,
    'EDIFFG': -0.01,
    'LREAL': 'Auto',
}

NONMETAL_POTENTIAL_MAP = {'P': 'P'}

# --- Oxygen reference (O2 molecule) ---
OXYGEN_PARAMS = {
    **COMMON_PARAMS,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 2,        # Relax atoms only (molecule)
    'NSW': 100,
    'EDIFFG': -0.01,
    'LREAL': False,   # Exact for small systems
}

OXYGEN_POTENTIAL_MAP = {'O': 'O'}

# --- Slab relaxation ---
SLAB_PARAMS = {
    **COMMON_PARAMS,
    'NELM': 100,
    'ISMEAR': 0,
    'SIGMA': 0.1,
    'IBRION': 2,
    'ISIF': 2,        # Relax atoms only (fixed cell)
    'NSW': 100,
    'EDIFFG': -0.1,   # Slightly relaxed for slabs
    'LREAL': 'Auto',
    'LWAVE': True,
    'LCHARG': True,
    'LASPH': True,
    # Uncomment for dipole correction if needed:
    # 'IDIPOL': 3,
    # 'LDIPOL': True,
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_resource_options(num_cores=NUM_CORES, queue=QUEUE_NAME):
    """Create computational resource options."""
    return {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': num_cores,
        },
        'queue_name': queue,
    }


def load_slab_structures(slabs_dir, slab_files):
    """Load predefined slab structures from CIF/POSCAR files."""
    input_slabs = {}

    for idx, slab_file in enumerate(slab_files):
        try:
            slab_path = os.path.join(slabs_dir, slab_file)
            atoms = read(slab_path)
            structure = orm.StructureData(ase=atoms)
            structure.store()
            input_slabs[f"term_{idx}"] = structure
            print(f"  ✓ Loaded {slab_file}")
        except FileNotFoundError:
            print(f"  ✗ Warning: {slab_file} not found")

    if not input_slabs:
        raise FileNotFoundError(
            f"No slab structures found in {slabs_dir}\n"
            f"Expected files: {', '.join(slab_files)}"
        )

    print(f"✓ Loaded {len(input_slabs)} slab structures\n")
    return input_slabs


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    """Run the complete PS-TEROS workflow."""

    # Initialize AiiDA
    print("=" * 60)
    print("PS-TEROS Complete Workflow - Ag3PO4")
    print("=" * 60)
    load_profile(profile='psteros')

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, STRUCTURES_DIR)
    slabs_dir = os.path.join(script_dir, SLABS_DIR)

    # Load slab structures if provided
    print("\nLoading slab structures...")
    input_slabs = load_slab_structures(slabs_dir, SLAB_FILES)

    # Create resource options
    options = get_resource_options()

    # Build WorkGraph
    print("Creating WorkGraph...")
    wg = build_core_workgraph(
        # Structures
        structures_dir=structures_dir,
        bulk_name=BULK_FILE,
        metal_name=METAL_FILE,
        nonmetal_name=NONMETAL_FILE,
        oxygen_name=OXYGEN_FILE,

        # Code and potentials
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,

        # Bulk
        bulk_potential_mapping=BULK_POTENTIAL_MAP,
        bulk_parameters=BULK_PARAMS,
        bulk_options=options,

        # Metal reference
        metal_potential_mapping=METAL_POTENTIAL_MAP,
        metal_parameters=METAL_PARAMS,
        metal_options=options,

        # Nonmetal reference
        nonmetal_potential_mapping=NONMETAL_POTENTIAL_MAP,
        nonmetal_parameters=NONMETAL_PARAMS,
        nonmetal_options=options,

        # Oxygen reference
        oxygen_potential_mapping=OXYGEN_POTENTIAL_MAP,
        oxygen_parameters=OXYGEN_PARAMS,
        oxygen_options=options,

        # Slab generation
        miller_indices=MILLER_INDICES,
        min_slab_thickness=MIN_SLAB_THICKNESS,
        min_vacuum_thickness=MIN_VACUUM_THICKNESS,

        # Slab relaxation
        slab_parameters=SLAB_PARAMS,
        slab_options=options,
        input_slabs=input_slabs,
        relax_slabs=RELAX_SLABS,

        # Calculation flags
        compute_relaxation_energy=COMPUTE_RELAXATION_ENERGY,
        compute_cleavage=COMPUTE_CLEAVAGE,
        compute_thermodynamics=COMPUTE_THERMODYNAMICS,
        thermodynamics_sampling=THERMODYNAMICS_SAMPLING,

        # Other
        kpoints_spacing=KPOINTS_SPACING,
        clean_workdir=CLEAN_WORKDIR,
        name=WORKFLOW_NAME,
    )

    print(f"✓ WorkGraph created: {wg.name}")
    print(f"  Total tasks: {len(wg.tasks)}")

    # Export visualization
    try:
        html_file = f'{WORKFLOW_NAME.lower()}.html'
        wg.to_html(html_file)
        print(f"  Visualization: {html_file}")
    except Exception as e:
        print(f"  Warning: Could not create visualization: {e}")

    # Submit
    print("\nSubmitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 60}")
    print(f"✓ Workflow submitted successfully!")
    print(f"{'=' * 60}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")
    print(f"  verdi process list")

    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, metal_energy, nonmetal_energy, oxygen_energy")
    print(f"  - formation_enthalpy")
    print(f"  - slab_structures, slab_energies, unrelaxed_slab_energies")
    print(f"  - relaxation_energies (E_relaxed - E_unrelaxed)")
    print(f"  - cleavage_energies")
    print(f"  - surface_energies")
    print()

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
