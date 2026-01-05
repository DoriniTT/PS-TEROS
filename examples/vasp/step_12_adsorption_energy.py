#!/usr/bin/env python
"""
STEP 12: Adsorption Energy Calculation

This script tests the adsorption energy workflow using the simplified API.

Workflow phases:
1. Initial relaxation of complete system (optional)
2. Structure separation using connectivity analysis (pymatgen CrystalNN)
3. SCF calculations on all three components (substrate, molecule, complete)
4. Adsorption energy: E_ads = E_complete - E_substrate - E_molecule

Test case: OH radical on Ag(111) at two different adsorption sites
- Site 1: OH on hollow site (3-fold coordination)
- Site 2: OH on top site (1-fold coordination)

Usage:
    # Activate your AiiDA environment first
    python step_12_adsorption_energy.py
"""

import sys
import os

# ==============================================================================
# CONFIGURATION - Modify these values for your setup
# ==============================================================================
try:
    from config import AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY, COMMON_OPTIONS
except ImportError:
    AIIDA_PROFILE = 'presto'
    VASP_CODE = 'VASP-6.4.1@cluster'
    POTENTIAL_FAMILY = 'PBE'
    COMMON_OPTIONS = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal',
    }
# ==============================================================================

from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def create_ag_oh_structure(site_type='hollow'):
    """Create Ag(111) slab with OH adsorbate for testing.

    Args:
        site_type: Type of adsorption site ('hollow' or 'top')

    Returns:
        orm.StructureData: Structure with Ag slab + OH
    """
    from pymatgen.core import Structure, Lattice

    # Create 2x2 Ag(111) slab (simplified for testing)
    lattice = Lattice.from_parameters(
        a=5.8, b=5.8, c=20.0,
        alpha=90, beta=90, gamma=90
    )

    # Ag atoms (4 atom slab)
    ag_positions = [
        [0.0, 0.0, 10.0],
        [2.9, 0.0, 10.0],
        [0.0, 2.9, 10.0],
        [2.9, 2.9, 10.0],
    ]

    # OH adsorbate position depends on site type
    if site_type == 'hollow':
        # OH on hollow site (3-fold coordination)
        oh_positions = [
            [1.45, 1.45, 12.0],  # O atom
            [1.45, 1.45, 13.0],  # H atom
        ]
    elif site_type == 'top':
        # OH on top site (1-fold coordination)
        oh_positions = [
            [0.0, 0.0, 12.0],  # O atom
            [0.0, 0.0, 13.0],  # H atom
        ]
    else:
        raise ValueError(f"Unknown site_type: {site_type}")

    species = ['Ag'] * 4 + ['O', 'H']
    positions = ag_positions + oh_positions

    structure = Structure(
        lattice, species, positions,
        coords_are_cartesian=True
    )

    return orm.StructureData(pymatgen=structure)


def main():
    """Step 12: Test adsorption energy calculation."""

    print("\n" + "="*70)
    print("STEP 12: ADSORPTION ENERGY CALCULATION")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile: {AIIDA_PROFILE}")

    # Create test structures
    print("\n2. Creating test structures...")
    print("   System: Ag(111) + OH")
    print("   Method: Connectivity analysis (pymatgen StructureGraph)")

    structure_hollow = create_ag_oh_structure(site_type='hollow')
    structure_top = create_ag_oh_structure(site_type='top')

    adsorption_structures = {
        'oh_hollow': structure_hollow,
        'oh_top': structure_top,
    }
    adsorption_formulas = {
        'oh_hollow': 'OH',
        'oh_top': 'OH',
    }

    print(f"   Site 1: OH on hollow site ({len(structure_hollow.get_ase())} atoms)")
    print(f"   Site 2: OH on top site ({len(structure_top.get_ase())} atoms)")

    print(f"\n3. VASP configuration:")
    print(f"   Code: {VASP_CODE}")
    print(f"   Potentials: {POTENTIAL_FAMILY}")

    # Relaxation parameters
    relax_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 400,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'NSW': 100,
        'ISIF': 2,
        'EDIFFG': -0.02,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # SCF parameters
    scf_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 400,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        'NELM': 100,
    }

    adsorption_potential_mapping = {'Ag': 'Ag', 'O': 'O', 'H': 'H'}

    print("\n4. Building workgraph...")
    print("   Preset: 'adsorption_energy'")
    print("   Phases:")
    print("     1. Relaxation of complete systems")
    print("     2. Structure separation")
    print("     3. SCF calculations")
    print("     4. Adsorption energy calculation")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")

    # Build workgraph
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # VASP code
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        clean_workdir=False,

        # Adsorption structures
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        adsorption_potential_mapping=adsorption_potential_mapping,

        # Relaxation phase
        relax_before_adsorption=True,
        adsorption_relax_builder_inputs={'parameters': {'incar': relax_parameters}},

        # SCF phase
        adsorption_scf_builder_inputs={'parameters': {'incar': scf_parameters}},

        # Scheduler and k-points
        adsorption_options=COMMON_OPTIONS,
        adsorption_kpoints_spacing=0.3,

        # Concurrency control
        max_concurrent_jobs=4,

        name='Step12_AdsorptionEnergy_Ag_OH',
    )

    print("   WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 12 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - separated_structures: Substrate/molecule/complete for each site")
    print(f"  - substrate_energies: E(Ag slab)")
    print(f"  - molecule_energies: E(OH)")
    print(f"  - complete_energies: E(Ag + OH)")
    print(f"  - adsorption_energies: E_ads for each site")
    print(f"\nNegative E_ads = favorable (exothermic) adsorption")
    print(f"Expected: E_ads(hollow) < E_ads(top)")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
