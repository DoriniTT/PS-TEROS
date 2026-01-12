"""
Example: Calculate Hubbard U for NiO using linear response method.

This example demonstrates how to use the u_calculation module to determine
the Hubbard U parameter for Ni d-electrons in NiO, following the VASP wiki
methodology: https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U

The workflow performs:
1. Ground state DFT (no +U) to get baseline d-occupancy
2. Non-SCF response calculations (ICHARG=11) for multiple potentials
3. SCF response calculations for multiple potentials
4. Linear regression to extract U from responses

Expected result: U ~ 5-6 eV for Ni in NiO (literature values vary 4-8 eV)

Usage:
    python step_20_hubbard_u_calculation.py
"""

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aiida import load_profile, orm
from pymatgen.core import Structure, Lattice

# Load AiiDA profile
load_profile('presto')

from teros.core.u_calculation import build_u_calculation_workgraph, get_u_calculation_results

# ==============================================================================
# CONFIGURATION - Adjust for your setup
# ==============================================================================

VASP_CODE = 'VASP-6.5.1@localwork'
POTENTIAL_FAMILY = 'PBE'
POTENTIAL_MAPPING = {'Ni': 'Ni', 'O': 'O'}

# Light parameters for local testing (8 cores)
LIGHT_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 8,
    },
    'max_wallclock_seconds': 3600,  # 1 hour per calculation
}

# Light INCAR parameters for testing
# Using smaller ENCUT and looser convergence for speed
# Note: AiiDA-VASP requires lowercase INCAR keys
LIGHT_GROUND_STATE_PARAMS = {
    'encut': 400,
    'ediff': 1e-5,
    'ismear': 0,
    'sigma': 0.05,
    'prec': 'Normal',
    'algo': 'Fast',
    'nelm': 100,
    'ispin': 2,  # Spin-polarized for NiO (antiferromagnetic)
}

# Potential values for linear response
# Note: V=0 should NOT be included as it causes inconsistency
# (GS has LDAU=False, response has LDAU=True)
# Using symmetric non-zero values gives better linear fit
TEST_POTENTIAL_VALUES = [-0.2, -0.1, 0.1, 0.2]

# K-point spacing (larger = fewer k-points = faster)
KPOINTS_SPACING = 0.1  # Coarse for testing


def create_nio_structure() -> orm.StructureData:
    """
    Create a simple NiO rocksalt structure.

    NiO has the rocksalt (NaCl) structure:
    - Ni at (0, 0, 0) and face centers
    - O at (0.5, 0.5, 0.5) and edge centers

    Lattice parameter: ~4.17 Å (experimental)
    """
    # NiO rocksalt structure
    a = 4.17  # Lattice parameter in Angstroms

    lattice = Lattice.cubic(a)

    # Rocksalt structure: Ni at (0,0,0), O at (0.5,0.5,0.5)
    species = ['Ni', 'O']
    coords = [
        [0.0, 0.0, 0.0],  # Ni
        [0.5, 0.5, 0.5],  # O
    ]

    pmg_structure = Structure(lattice, species, coords)

    # Convert to AiiDA StructureData
    structure = orm.StructureData(pymatgen=pmg_structure)
    structure.label = 'NiO_rocksalt'
    structure.description = 'NiO rocksalt structure for Hubbard U calculation'

    return structure


def main():
    """Run the Hubbard U calculation for NiO."""

    print("=" * 70)
    print("Hubbard U Calculation for NiO")
    print("=" * 70)

    # Create NiO structure
    print("\n1. Creating NiO structure...")
    structure = create_nio_structure()
    structure.store()
    print(f"   Structure PK: {structure.pk}")
    print(f"   Formula: {structure.get_formula()}")
    print(f"   Cell volume: {structure.get_cell_volume():.2f} Å³")

    # Build the workgraph
    print("\n2. Building Hubbard U workgraph...")
    wg = build_u_calculation_workgraph(
        structure=structure,
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        potential_mapping=POTENTIAL_MAPPING,
        target_species='Ni',
        potential_values=TEST_POTENTIAL_VALUES,
        ldaul=2,  # d-electrons
        ldauj=0.0,
        ground_state_parameters=LIGHT_GROUND_STATE_PARAMS,
        response_parameters=LIGHT_GROUND_STATE_PARAMS,
        options=LIGHT_OPTIONS,
        kpoints_spacing=KPOINTS_SPACING,
        clean_workdir=False,  # Keep files for debugging
        name='NiO_HubbardU_Test',
    )

    print(f"   WorkGraph name: {wg.name}")
    print(f"   Number of tasks: {len(wg.tasks)}")
    print(f"   Potential values: {TEST_POTENTIAL_VALUES}")

    # Print task overview
    print("\n   Tasks:")
    for task in wg.tasks:
        print(f"     - {task.name}")

    # Submit the workgraph
    print("\n3. Submitting workgraph...")
    wg.submit()
    print(f"   WorkGraph PK: {wg.pk}")
    print(f"\n   Monitor with: verdi process show {wg.pk}")
    print(f"   Detailed report: verdi process report {wg.pk}")

    print("\n" + "=" * 70)
    print("WorkGraph submitted! Use the commands above to monitor progress.")
    print("=" * 70)

    return wg


if __name__ == '__main__':
    wg = main()
