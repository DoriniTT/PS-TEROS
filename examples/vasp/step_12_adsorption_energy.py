#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 12: Adsorption Energy Calculation

This script tests the adsorption energy workflow:
- Structure separation (substrate + molecule identification)
- Substrate relaxation
- Molecule relaxation
- Complete system relaxation
- Adsorption energy calculation: E_ads = E_complete - E_substrate - E_molecule

The module uses connectivity analysis (pymatgen StructureGraph with CrystalNN)
to automatically identify and separate the adsorbate from the substrate.

Material: Ag(111) surface with OH adsorbate
Adsorbate: OH radical

Usage:
    source ~/envs/aiida/bin/activate
    python step_12_adsorption_energy.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def create_ag_oh_structure():
    """Create Ag(111) slab with OH adsorbate for testing.

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

    # OH adsorbate on top (bonded O-H cluster above surface)
    oh_positions = [
        [1.45, 1.45, 12.0],  # O atom on hollow site
        [1.45, 1.45, 13.0],  # H atom above O
    ]

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

    # Create test structure
    print("\n2. Creating test structure...")
    print("   System: Ag(111) + OH")
    print("   Adsorbate: OH radical")
    print("   Method: Connectivity analysis (pymatgen StructureGraph)")

    complete_structure = create_ag_oh_structure()

    adsorption_structures = [complete_structure]
    adsorption_formulas = ['OH']

    print("   ✓ Structure created")
    print(f"   Total atoms: {len(complete_structure.get_ase())}")
    print(f"   Composition: Ag4OH")

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
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Potential mapping
    adsorption_potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
        'H': 'H',
    }

    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   This will create 3 VASP calculations:")
    print("     1. Complete system (Ag slab + OH)")
    print("     2. Bare substrate (Ag slab only)")
    print("     3. Isolated molecule (OH in same cell)")
    print("   ")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Negative E_ads = favorable (exothermic) adsorption")
    print("   Positive E_ads = unfavorable (endothermic) adsorption")

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
    print(f"     - separated_structures (list with 1 dict)")
    print(f"       * substrate: Ag slab (4 atoms)")
    print(f"       * molecule: OH (2 atoms)")
    print(f"       * complete: Ag + OH (6 atoms)")
    print(f"  ")
    print(f"  2. Individual energies:")
    print(f"     - substrate_energies (list): E(Ag slab)")
    print(f"     - molecule_energies (list): E(OH)")
    print(f"     - complete_energies (list): E(Ag+OH)")
    print(f"  ")
    print(f"  3. Adsorption energy:")
    print(f"     - adsorption_energies (list): E_ads (eV)")
    print(f"  ")
    print(f"Expected E_ads for OH/Ag(111):")
    print(f"  Literature range: -2.0 to -2.5 eV (DFT-PBE)")
    print(f"  (Negative = favorable adsorption)")
    print(f"\nCell handling:")
    print(f"  All three systems use the SAME simulation cell")
    print(f"  This eliminates basis set superposition error (BSSE)")
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
