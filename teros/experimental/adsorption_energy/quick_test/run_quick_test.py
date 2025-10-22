#!/home/thiagotd/envs/aiida/bin/python
"""
Quick Test: Dual Adsorption Energy Calculation (OH + H2O on Ag)

This is a MINIMAL test script for the adsorption energy module.
Based on run_lamno3_oh_adsorption.py but with:
- TWO tiny structures tested in PARALLEL
- Minimal VASP parameters for FAST testing
- Same workflow: relax → separate → SCF → calculate E_ads

Systems (tested in parallel):
1. Ag(111) + OH:
   - Substrate: Ag(111) surface (2x2, single layer = 4 atoms)
   - Adsorbate: OH radical (2 atoms)
   - Total: 6 atoms

2. Ag(111) + H2O:
   - Substrate: Ag(111) surface (2x2, single layer = 4 atoms)
   - Adsorbate: H2O molecule (3 atoms)
   - Total: 7 atoms

Workflow (4 phases per structure):
1. Relax complete systems (Ag+OH and Ag+H2O in parallel)
2. Separate relaxed structures
3. Run 6 SCF calculations in parallel (3 per structure)
4. Calculate E_ads = E_complete - E_substrate - E_molecule

Total VASP jobs: 2 relax + 6 SCF = 8 calculations

Usage:
    source ~/envs/aiida/bin/activate
    python run_quick_test.py
"""

import sys
import os
from pathlib import Path
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def main():
    """Run quick Ag + OH adsorption energy test."""

    print("=" * 70)
    print("QUICK TEST: Ag + OH ADSORPTION ENERGY")
    print("=" * 70)

    # 1. Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # 2. Load structures (TWO for parallel testing)
    print("\n2. Loading structures...")
    script_dir = Path(__file__).parent

    # Structure 1: Ag + OH
    structure_file_oh = script_dir / "Ag_OH_quick_test.cif"
    if not structure_file_oh.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_file_oh}")

    # Structure 2: Ag + H2O
    structure_file_h2o = script_dir / "Ag_H2O_quick_test.cif"
    if not structure_file_h2o.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_file_h2o}")

    # Load both structures using AiiDA
    from ase.io import read

    ase_structure_oh = read(str(structure_file_oh))
    complete_structure_oh = orm.StructureData(ase=ase_structure_oh)

    ase_structure_h2o = read(str(structure_file_h2o))
    complete_structure_h2o = orm.StructureData(ase=ase_structure_h2o)

    print(f"   Structure 1 (OH):")
    print(f"     File: {structure_file_oh.name}")
    print(f"     Formula: {ase_structure_oh.get_chemical_formula()}")
    print(f"     Total atoms: {len(complete_structure_oh.get_ase())}")
    print(f"     Substrate: Ag4, Adsorbate: OH")

    print(f"   Structure 2 (H2O):")
    print(f"     File: {structure_file_h2o.name}")
    print(f"     Formula: {ase_structure_h2o.get_chemical_formula()}")
    print(f"     Total atoms: {len(complete_structure_h2o.get_ase())}")
    print(f"     Substrate: Ag4, Adsorbate: H2O")

    print("   ✓ Both structures loaded")

    # 3. VASP Configuration
    print("\n3. VASP Configuration:")
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # Adsorption structures and formulas (TWO structures for parallel testing)
    adsorption_structures = {
        'oh_site': complete_structure_oh,   # Ag + OH (6 atoms)
        'h2o_site': complete_structure_h2o,  # Ag + H2O (7 atoms)
    }

    adsorption_formulas = {
        'oh_site': 'OH',    # Adsorbate to identify and separate
        'h2o_site': 'H2O',  # Adsorbate to identify and separate
    }

    # Common settings
    common_potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
        'H': 'H',
    }

    common_options = {
        'resources': {
            'num_machines': 1,  # Minimal resources for quick test
            'num_cores_per_machine': 5,
        },
        'queue_name': 'teste',
    }

    kpoints_spacing = 1.0  # COARSE k-points for speed (vs 0.6 in production)

    # MINIMAL Relaxation parameters for QUICK testing
    # Much reduced from production parameters
    relax_builder_inputs = {
        'relax_settings': {
            'algo': 'cg',
            'force_cutoff': 1,
            'steps': 500,  # Only 500 steps (vs 500 in production)
            'positions': True,
            'shape': False,
            'volume': False,
            'convergence_on': True,
            'convergence_absolute': False,
            'convergence_max_iterations': 2,  # Only 2 iterations (vs 5)
            'convergence_positions': 0.01,
            'convergence_volume': 0.01,
            'convergence_shape_lengths': 0.1,
            'convergence_shape_angles': 0.1,
            'convergence_mode': 'last',  # REQUIRED: 'last' or 'inout'
            'perform': True,
        },
        'vasp': {
            'parameters': {'incar': {
                'PREC': 'Normal',  # Normal precision (vs Accurate)
                'ALGO': 'Fast',    # Fast algorithm
                'LREAL': 'Auto',
                'ENCUT': 250,      # Low cutoff (vs 400 eV)
                'EDIFF': 1e-4,     # Loose convergence (vs 1e-5)
                'NELM': 200,        # Fewer electronic steps (vs 40)
                'NELMIN': 3,       # Fewer min steps (vs 6)
                'ISMEAR': 0,
                'SIGMA': 0.1,      # Larger smearing for speed
                'ISPIN': 1,        # NO spin polarization for speed
                'LWAVE': False,    # Don't write WAVECAR
                'LCHARG': False,   # Don't write CHGCAR
                'NCORE': 3,
                'LASPH': False,    # No aspherical terms for speed
            }},
            'options': common_options,
            'potential_family': potential_family,
            'potential_mapping': common_potential_mapping,
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': False,
        }
    }

    # MINIMAL SCF parameters
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Normal',
            'ALGO': 'Fast',
            'ENCUT': 250,
            'EDIFF': 1e-4,
            'ISMEAR': 0,
            'NELMIN': 3,
            'SIGMA': 0.1,
            'ISPIN': 1,  # NO spin polarization
            'LWAVE': False,
            'LCHARG': False,
            'NCORE': 3,
            'LASPH': False,
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,  # REQUIRED for VaspWorkChain
    }

    # Atom fixing - test the feature even though structure is tiny
    fix_atoms = True
    fix_type = 'bottom'
    fix_thickness = 2.0  # Fix bottom 2 Å (the Ag layer)
    fix_elements = None

    print(f"   K-points spacing: {kpoints_spacing} Å⁻¹ (COARSE for speed)")
    print(f"   ENCUT: 250 eV (LOW for speed)")
    print(f"   Spin-polarized: No (for speed)")
    print(f"   Relaxation: ENABLED (relax → separate → SCF)")
    print(f"   Relax max ionic steps: {relax_builder_inputs['relax_settings']['steps']}")
    print(f"   Atom fixing: {'ENABLED' if fix_atoms else 'DISABLED'}")
    if fix_atoms:
        print(f"   Fix type: {fix_type}, thickness: {fix_thickness} Å")

    # 4. Build WorkGraph
    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   Testing: TWO structures in PARALLEL")
    print("")
    print("   Workflow phases (for each structure):")
    print("     Phase 1: Relax complete systems")
    print("       - oh_site: Ag + OH (6 atoms)")
    print("       - h2o_site: Ag + H2O (7 atoms)")
    print("     Phase 2: Separate relaxed structures")
    print("     Phase 3: SCF calculations (3 single-point per structure)")
    print("       - Substrate (Ag)")
    print("       - Molecule (OH or H2O)")
    print("       - Complete (Ag+OH or Ag+H2O)")
    print("     Phase 4: Calculate adsorption energies")
    print("")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Total VASP jobs: 2 relax + 6 SCF = 8 calculations")

    # Build workgraph
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,

        # Adsorption energy specific parameters
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        adsorption_potential_mapping=common_potential_mapping,

        # Relaxation and SCF
        relax_before_adsorption=True,
        adsorption_relax_builder_inputs=relax_builder_inputs,
        adsorption_scf_builder_inputs=scf_builder_inputs,

        # Atom fixing
        adsorption_fix_atoms=fix_atoms,
        adsorption_fix_type=fix_type,
        adsorption_fix_thickness=fix_thickness,
        adsorption_fix_elements=fix_elements,

        name='Quick_Test_Dual_Adsorption',
    )

    print("   ✓ WorkGraph built successfully")

    # 5. Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 70}")
    print("DUAL ADSORPTION TEST SUBMITTED SUCCESSFULLY")
    print(f"{'=' * 70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  - adsorption_energies['oh_site']: E_ads for OH on Ag (eV)")
    print(f"  - adsorption_energies['h2o_site']: E_ads for H2O on Ag (eV)")
    print(f"\nExpected runtime:")
    print(f"  Phase 1 (2 relax): ~5-15 min each (6-7 atoms, parallel)")
    print(f"  Phase 3 (6 SCF): ~2-5 min each (parallel)")
    print(f"  Total: ~10-20 minutes")
    print(f"\nAdvantages over LaMnO3 example:")
    print(f"  - Runtime: 10-20 min vs 12-18 hours (~40× faster)")
    print(f"  - Tests: 2 structures vs 1 (validates parallel execution)")
    print(f"  - Size: 6-7 atoms vs 150 atoms")
    print(f"{'=' * 70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
