#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 12: Adsorption Energy Calculation (Multi-Site Test)

This script tests the adsorption energy workflow with multiple structures using
the simplified API that directly accepts INCAR parameters.

Workflow (4 phases):
1. Initial relaxation of complete system (optional, controlled by relax_before_adsorption)
2. Structure separation using connectivity analysis (pymatgen StructureGraph with CrystalNN)
3. SCF calculations on all three components (substrate, molecule, complete)
4. Adsorption energy calculation: E_ads = E_complete - E_substrate - E_molecule

Test case: OH radical on Ag(111) at two different adsorption sites
- Site 1: OH on hollow site (3-fold coordination)
- Site 2: OH on top site (1-fold coordination, bonded to Ag)

This tests the scatter functionality with parallel VASP calculations:
  If relax_before_adsorption=True:  2 relax + (2 sites × 3 SCF) = 8 calculations
  If relax_before_adsorption=False: 2 sites × 3 SCF = 6 calculations

Key API features:
- Uses only vasp.v2.vasp plugin (not vasp.v2.relax)
- Direct INCAR parameters via builder_inputs
- Relaxation controlled by NSW parameter (NSW > 0 = relax, NSW = 0 = SCF)

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

    # Relaxation parameters for initial complete system relaxation (NSW > 0)
    # These are used when relax_before_adsorption=True
    relax_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 400,          # Adequate for Ag and light elements
        'EDIFF': 1e-5,
        'ISMEAR': 0,           # Gaussian smearing for molecules
        'SIGMA': 0.05,
        'IBRION': 2,           # Conjugate gradient
        'NSW': 100,            # Relaxation steps (>0 = relaxation)
        'ISIF': 2,             # Relax ions only (keep cell fixed)
        'EDIFFG': -0.02,       # Force convergence
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        'NCORE': 4,            # Parallel efficiency
    }

    # SCF parameters for final energy calculations (NSW=0 set automatically)
    # These are used for substrate, molecule, and complete system energies
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
        'NCORE': 4,
        'NELM': 100,           # SCF iterations
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

    # ===== STRUCTURE-SPECIFIC BUILDER INPUTS (OPTIONAL - NEW FEATURE) =====
    # Override parameters for specific structures using structure indices (0, 1, 2, ...)
    # Indices match the order in adsorption_structures dict (oh_hollow=0, oh_top=1)

    # Example 1: Override relaxation parameters for structure 0 (oh_hollow)
    # structure_specific_relax = {
    #     0: {  # oh_hollow - use different algorithm
    #         'parameters': {
    #             'incar': {
    #                 'ALGO': 'All',      # Try all algorithms
    #                 'EDIFFG': -0.03,    # Looser convergence
    #             }
    #         },
    #         'kpoints_spacing': 0.4,  # Coarser k-points
    #     }
    # }

    # Example 2: Override SCF parameters for structure 1 (oh_top)
    # structure_specific_scf = {
    #     1: {  # oh_top - use tighter convergence for final energy
    #         'parameters': {
    #             'incar': {
    #                 'EDIFF': 1e-6,      # Tighter electronic convergence
    #                 'PREC': 'Accurate',
    #             }
    #         },
    #         'kpoints_spacing': 0.2,  # Denser k-points for final SCF
    #     }
    # }

    # For this example, we'll use default settings for all structures (no overrides)
    structure_specific_relax = None
    structure_specific_scf = None

    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   Workflow phases:")
    print("     1. Initial relaxation of complete systems (if enabled)")
    print("     2. Structure separation (connectivity analysis)")
    print("     3. SCF calculations on separated components")
    print("     4. Adsorption energy calculation")
    if structure_specific_relax is not None or structure_specific_scf is not None:
        print("   Structure-specific overrides: ENABLED")
        if structure_specific_relax:
            print(f"     Relaxation overrides for structures: {list(structure_specific_relax.keys())}")
        if structure_specific_scf:
            print(f"     SCF overrides for structures: {list(structure_specific_scf.keys())}")
    else:
        print("   Structure-specific overrides: DISABLED (using defaults for all)")
    print("   ")
    print("   This will create 8 VASP calculations (with relax_before_adsorption=True):")
    print("     Phase 1: 2 relaxations (one per site)")
    print("       - oh_hollow: Ag slab + OH (relaxation)")
    print("       - oh_top:    Ag slab + OH (relaxation)")
    print("     Phase 3: 6 SCF calculations (3 per site)")
    print("       For each site (hollow and top):")
    print("         1. Complete system (Ag slab + OH) - SCF")
    print("         2. Bare substrate (Ag slab only) - SCF")
    print("         3. Isolated molecule (OH in same cell) - SCF")
    print("   ")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Negative E_ads = favorable (exothermic) adsorption")
    print("   Positive E_ads = unfavorable (endothermic) adsorption")
    print("   ")
    print("   Expected: E_ads(hollow) < E_ads(top)")
    print("   (Hollow sites typically more favorable than top sites)")

    # Build workgraph using adsorption_energy preset with simplified API
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,

        # Adsorption energy structures
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        adsorption_potential_mapping=adsorption_potential_mapping,

        # Simplified API: Direct INCAR parameters
        # Phase 1: Relax complete system before separation
        relax_before_adsorption=True,
        adsorption_relax_builder_inputs={'parameters': {'incar': relax_parameters}},

        # Phase 3: SCF calculations on separated components
        adsorption_scf_builder_inputs={'parameters': {'incar': scf_parameters}},

        # Scheduler and k-points
        adsorption_options=adsorption_options,
        adsorption_kpoints_spacing=0.3,

        # NEW FEATURE: Structure-specific builder inputs (optional)
        adsorption_structure_specific_relax_builder_inputs=structure_specific_relax,
        adsorption_structure_specific_scf_builder_inputs=structure_specific_scf,
        # Concurrency control (limits simultaneous VASP calculations)
        max_concurrent_jobs=4,  # Default: 4 concurrent calculations


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
    print(f"  With relax_before_adsorption=True:")
    print(f"    Phase 1: 2 relaxations run in parallel")
    print(f"    Phase 3: 6 SCF calculations run in parallel")
    print(f"  Total: 2 relax + (2 sites × 3 SCF) = 8 VASP calculations")
    print(f"\nKey API change:")
    print(f"  Old: adsorption_parameters (single dict)")
    print(f"  New: adsorption_relax_builder_inputs + adsorption_scf_builder_inputs")
    print(f"       (separate dicts for relaxation and SCF)")
    print(f"  Plugin: Uses only vasp.v2.vasp (not vasp.v2.relax)")
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
