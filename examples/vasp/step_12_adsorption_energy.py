#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 12: Adsorption Energy Calculation (Multi-Site Test)

This script tests the adsorption energy workflow with multiple structures:
- Structure separation (substrate + molecule identification)
- Substrate relaxation
- Molecule relaxation
- Complete system relaxation
- Adsorption energy calculation: E_ads = E_complete - E_substrate - E_molecule

The module uses connectivity analysis (pymatgen StructureGraph with CrystalNN)
to automatically identify and separate the adsorbate from the substrate.

Test case: OH radical on Ag(111) at two different adsorption sites
- Site 1: OH on hollow site (3-fold coordination)
- Site 2: OH on top site (1-fold coordination, bonded to Ag)

This tests the scatter functionality with 6 parallel VASP calculations:
  2 sites × 3 calculations each (substrate, molecule, complete)

Usage:
    source ~/envs/aiida/bin/activate
    python step_12_adsorption_energy.py
"""

import sys
import os
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

    # Create 2x2 Ag(111) slab
    lattice = Lattice.from_parameters(
        a=5.8, b=5.8, c=20.0,
        alpha=90, beta=90, gamma=90
    )

    # Ag atoms (4 atom slab - simplified for testing)
    ag_positions = [
        [0.0, 0.0, 10.0],   # Bottom layer
        [2.9, 0.0, 10.0],
        [0.0, 2.9, 10.0],
        [2.9, 2.9, 10.0],
    ]

    # OH adsorbate position depends on site type
    if site_type == 'hollow':
        # OH on hollow site (3-fold coordination)
        oh_positions = [
            [1.45, 1.45, 12.0],  # O atom on hollow site
            [1.45, 1.45, 13.0],  # H atom above O
        ]
    elif site_type == 'top':
        # OH on top site (directly above Ag atom, 1-fold coordination)
        oh_positions = [
            [0.0, 0.0, 12.0],    # O atom directly above Ag
            [0.0, 0.0, 13.0],    # H atom above O
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
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Create test structures
    print("\n2. Creating test structures...")
    print("   System: Ag(111) + OH")
    print("   Adsorbate: OH radical")
    print("   Method: Connectivity analysis (pymatgen StructureGraph)")

    # Create two structures with OH at different adsorption sites
    structure_hollow = create_ag_oh_structure(site_type='hollow')
    structure_top = create_ag_oh_structure(site_type='top')

    # Adsorption structures and formulas must be dicts with matching keys
    adsorption_structures = {
        'oh_hollow': structure_hollow,
        'oh_top': structure_top,
    }
    adsorption_formulas = {
        'oh_hollow': 'OH',
        'oh_top': 'OH',
    }

    print("   ✓ Structures created")
    print(f"   Site 1: OH on hollow site ({len(structure_hollow.get_ase())} atoms)")
    print(f"   Site 2: OH on top site ({len(structure_top.get_ase())} atoms)")
    print(f"   Composition: Ag4OH (each)")

    # Code configuration
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    print(f"\n3. VASP Configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # VASP parameters (production-quality for Ag + OH)
    adsorption_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 400,          # Adequate for Ag and light elements
        'EDIFF': 1e-5,
        'ISMEAR': 0,           # Gaussian smearing for molecules
        'SIGMA': 0.05,
        'IBRION': 2,           # Conjugate gradient
        'NSW': 100,            # Relaxation steps
        'ISIF': 2,             # Relax ions only (keep cell fixed)
        'EDIFFG': -0.02,       # Force convergence
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        'NCORE': 4,            # Parallel efficiency
    }

    # Scheduler options
    adsorption_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 5,
        },
        'queue_name': 'teste',
    }

    # Potential mapping
    adsorption_potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
        'H': 'H',
    }

    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   This will create 6 VASP calculations (3 per structure):")
    print("     For each site (hollow and top):")
    print("       1. Complete system (Ag slab + OH)")
    print("       2. Bare substrate (Ag slab only)")
    print("       3. Isolated molecule (OH in same cell)")
    print("   ")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Negative E_ads = favorable (exothermic) adsorption")
    print("   Positive E_ads = unfavorable (endothermic) adsorption")
    print("   ")
    print("   Expected: E_ads(hollow) < E_ads(top)")
    print("   (Hollow sites typically more favorable than top sites)")

    # Build workgraph using adsorption_energy preset
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,

        # Adsorption energy specific
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        adsorption_parameters=adsorption_parameters,
        adsorption_options=adsorption_options,
        adsorption_potential_mapping=adsorption_potential_mapping,
        adsorption_kpoints_spacing=0.3,

        name='Step12_AdsorptionEnergy_Ag_OH',
    )

    print("   ✓ WorkGraph built successfully")

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
    print(f"  1. Separated structures:")
    print(f"     - separated_structures (dict with 2 entries)")
    print(f"       * oh_hollow: {{substrate, molecule, complete}}")
    print(f"       * oh_top: {{substrate, molecule, complete}}")
    print(f"     Each contains: Ag slab (4 atoms), OH (2 atoms), complete (6 atoms)")
    print(f"  ")
    print(f"  2. Individual energies (dict with 2 entries each):")
    print(f"     - substrate_energies: {{oh_hollow: E(Ag), oh_top: E(Ag)}}")
    print(f"     - molecule_energies: {{oh_hollow: E(OH), oh_top: E(OH)}}")
    print(f"     - complete_energies: {{oh_hollow: E(Ag+OH), oh_top: E(Ag+OH)}}")
    print(f"  ")
    print(f"  3. Adsorption energies (dict with 2 entries):")
    print(f"     - adsorption_energies: {{oh_hollow: E_ads, oh_top: E_ads}}")
    print(f"  ")
    print(f"Expected E_ads for OH/Ag(111):")
    print(f"  Hollow site: -2.0 to -2.5 eV (DFT-PBE, more favorable)")
    print(f"  Top site:    -1.5 to -2.0 eV (DFT-PBE, less favorable)")
    print(f"  (Negative = favorable adsorption)")
    print(f"\nCell handling:")
    print(f"  All systems use the SAME simulation cell (per structure)")
    print(f"  This eliminates basis set superposition error (BSSE)")
    print(f"\nParallel execution:")
    print(f"  6 VASP calculations will run in parallel")
    print(f"  Total: 2 sites × 3 calculations each")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
