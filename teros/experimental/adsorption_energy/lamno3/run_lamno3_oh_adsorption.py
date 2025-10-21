#!/home/thiagotd/envs/aiida/bin/python
"""
LaMnO3 + OH Adsorption Energy Calculation with Relaxation and Atom Fixing

This script demonstrates the NEW features of the adsorption energy module:
1. Full builder-based API for vasp.v2.relax (complete control)
2. Atom fixing for slab relaxations (fix bottom layers)

System:
- Substrate: LaMnO3 (100) surface (2x2x1 supercell)
- Adsorbate: OH radical
- Structure: La32Mn28O89H (complete system)

Workflow (4 phases):
1. Relax complete system (LaMnO3 + OH) using vasp.v2.relax
   - NEW: Uses full builder inputs for relax_settings control
   - NEW: Fixes bottom 7 Å of slab to simulate bulk
2. Separate relaxed structure into substrate + adsorbate
3. Run 3 SCF calculations in parallel using vasp.v2.vasp:
   - Complete system (from relaxed)
   - Bare substrate (from relaxed)
   - Isolated OH molecule (from relaxed)
4. Calculate E_ads = E_complete - E_substrate - E_molecule

NEW Features Demonstrated:
- relax_builder_inputs: Full control over vasp.v2.relax with relax_settings Dict
- adsorption_fix_atoms: Fix bottom layers of slab during relaxation
- adsorption_fix_type: 'bottom' (also supports 'top', 'center')
- adsorption_fix_thickness: 7.0 Å thickness for fixed region

Usage:
    source ~/envs/aiida/bin/activate
    python run_lamno3_oh_adsorption.py
"""

import sys
import os
from pathlib import Path
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def main():
    """Run LaMnO3 + OH adsorption energy calculation."""

    print("=" * 70)
    print("LaMnO3 + OH ADSORPTION ENERGY CALCULATION")
    print("=" * 70)

    # 1. Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # 2. Load structure
    print("\n2. Loading structure...")
    script_dir = Path(__file__).parent
    structure_file = script_dir / "LaMnO3_100_A4_surface_2x2x1_OOH.cif"

    if not structure_file.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_file}")

    # Load structure using AiiDA
    from ase.io import read
    ase_structure = read(str(structure_file))
    complete_structure = orm.StructureData(ase=ase_structure)

    print(f"   File: {structure_file.name}")
    print(f"   Formula: La32Mn28O89H")
    print(f"   Total atoms: {len(complete_structure.get_ase())}")
    print(f"   Substrate: La32Mn28O88")
    print(f"   Adsorbate: OH")
    print("   ✓ Structure loaded")

    # 3. VASP Configuration
    print("\n3. VASP Configuration:")
    code_label = 'VASPGAM-6.5.0@lovelace-parexp'
    potential_family = 'PBE'

    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # Adsorption structures and formulas (dicts with matching keys)
    adsorption_structures = {
        'oh_lamno3': complete_structure,
    }

    adsorption_formulas = {
        'oh_lamno3': 'OH',  # The adsorbate to identify and separate
    }

    # Common settings
    common_potential_mapping = {
        'La': 'La',
        'Mn': 'Mn',  # Use _pv for transition metals
        'O': 'O',
        'H': 'H',
    }

    common_options = {
        'resources': {
            'num_machines': 4,
            'num_cores_per_machine': 48,
        },
        'queue_name': 'parexp',
    }

    kpoints_spacing = 0.6  # Angstrom^-1

    # Relaxation parameters (Phase 1: relax complete system)
    # Using NEW builder-based API for full control over vasp.v2.relax workchain
    # This provides complete control over relax_settings and VASP parameters
    # NOTE: Use plain Python dicts - they will be converted to orm.Dict internally
    relax_builder_inputs = {
        'relax_settings': {  # Plain dict - converted to orm.Dict internally
            'algo': 'cg',  # Conjugate gradient
            'force_cutoff': 0.1,  # 0.1 eV/Å force convergence
            'steps': 500,  # Maximum ionic steps
            'positions': True,  # Relax atomic positions
            'shape': False,  # Don't relax cell shape
            'volume': False,  # Don't relax cell volume
            'convergence_on': True,  # Enable convergence checks
            'convergence_absolute': False,
            'convergence_max_iterations': 5,
            'convergence_positions': 0.01,  # 0.01 Å convergence
            'convergence_volume': 0.01,
            'convergence_shape_lengths': 0.1,
            'convergence_shape_angles': 0.1,
            'perform': True,
        },
        'vasp': {
            'parameters': {'incar': {
                'PREC': 'Accurate',
                'ALGO': 'All',
                'LREAL': 'Auto',
                'ENCUT': 400,           # Cutoff energy
                'EDIFF': 1e-5,          # Electronic convergence
                'NELM': 40,
                'NELMIN': 6,            # Minimum electronic steps
                'ISMEAR': 0,            # Gaussian smearing
                'SIGMA': 0.05,          # Small smearing width
                'ISPIN': 2,             # Spin-polarized
                #'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments
                'LORBIT': 11,           # DOSCAR and PROCAR
                'LWAVE': True,          # Write WAVECAR
                'LCHARG': True,         # Write CHGCAR
                'NCORE': 3,             # Parallelization
                'KPAR': 4,              # k-point parallelization
                'LASPH': True,          # Non-spherical contributions
                'IVDW': 12,             # DFT-D3 van der Waals
                #'IDIPOL': 3,           # Dipole correction for slab
                #'LDIPOL': True,        # Enable dipole correction
            }},
            'options': common_options,
            'potential_family': potential_family,
            'potential_mapping': common_potential_mapping,
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': False,
        }
    }

    # SCF inputs (Phase 3: single-point calculations)
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ALGO': 'Normal',
            'ENCUT': 400,           # Cutoff energy
            'EDIFF': 1e-5,          # Electronic convergence
            'ISMEAR': 0,            # Gaussian smearing
            'NELMIN': 6,          # Minimum electronic steps
            'SIGMA': 0.05,          # Small smearing width
            'ISPIN': 2,             # Spin-polarized
            #'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments
            'LORBIT': 11,           # DOSCAR and PROCAR
            # NSW=0 and IBRION=-1 will be set automatically by workflow
            'LWAVE': True,         # Don't write WAVECAR
            'LCHARG': True,         # Write CHGCAR
            'NCORE': 3,             # Parallelization
            'KPAR': 4,             # k-point parallelization
            'LASPH': True,          # Non-spherical contributions
            'IVDW': 12,             # DFT-D3 van der Waals
            #'IDIPOL': 3,            # Dipole correction for slab
            #'LDIPOL': True,         # Enable dipole correction
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        # Note: k-points spacing passed separately to workgraph
    }

    # Atom fixing parameters (NEW feature)
    # Fix bottom 7 Å of the slab to simulate bulk-like behavior
    fix_atoms = True
    fix_type = 'bottom'
    fix_thickness = 7.0  # Angstrom
    fix_elements = None  # None = fix all elements in the region

    print(f"   K-points spacing: {kpoints_spacing} Å⁻¹")
    print(f"   Spin-polarized: Yes")
    print(f"   Relaxation: ENABLED (relax → separate → SCF)")
    print(f"   Relax max ionic steps: {relax_builder_inputs['relax_settings']['steps']}")
    print(f"   Atom fixing: {'ENABLED' if fix_atoms else 'DISABLED'}")
    if fix_atoms:
        print(f"   Fix type: {fix_type}, thickness: {fix_thickness} Å")

    # 4. Build WorkGraph
    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   Workflow phases:")
    print("     Phase 1: Relax complete system (LaMnO3 + OH)")
    print("     Phase 2: Separate relaxed structure")
    print("     Phase 3: SCF calculations (3 single-point calculations)")
    print("       - Substrate (LaMnO3 from relaxed)")
    print("       - Molecule (OH from relaxed)")
    print("       - Complete (LaMnO3+OH from relaxed)")
    print("")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Negative E_ads = favorable (exothermic) adsorption")
    print("   Positive E_ads = unfavorable (endothermic) adsorption")

    # Build workgraph using adsorption_energy preset with relaxation
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,  # Keep files for analysis

        # Adsorption energy specific parameters
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,
        adsorption_potential_mapping=common_potential_mapping,

        # NEW: Full builder control for relaxation
        relax_before_adsorption=True,
        adsorption_relax_builder_inputs=relax_builder_inputs,  # Full builder inputs for relax
        adsorption_scf_builder_inputs=scf_builder_inputs,  # Builder inputs for SCF phase

        # NEW: Atom fixing for slabs
        adsorption_fix_atoms=fix_atoms,
        adsorption_fix_type=fix_type,
        adsorption_fix_thickness=fix_thickness,
        adsorption_fix_elements=fix_elements,

        name='LaMnO3_OH_RelaxThenAdsorption_WithFixing',
    )

    print("   ✓ WorkGraph built successfully")

    # 5. Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 70}")
    print("CALCULATION SUBMITTED SUCCESSFULLY")
    print(f"{'=' * 70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nExpected outputs:")
    print(f"  ")
    print(f"  1. Relaxed structure:")
    print(f"     - relaxed_complete_structures['oh_lamno3']: Relaxed La32Mn28O89H")
    print(f"  ")
    print(f"  2. Separated structures:")
    print(f"     - separated_structures['oh_lamno3']:")
    print(f"       * substrate: La32Mn28O88 (149 atoms)")
    print(f"       * molecule: OH (2 atoms)")
    print(f"       * complete: La32Mn28O89H (150 atoms)")
    print(f"  ")
    print(f"  3. Individual energies:")
    print(f"     - substrate_energies['oh_lamno3']: E(LaMnO3) from SCF")
    print(f"     - molecule_energies['oh_lamno3']: E(OH) from SCF")
    print(f"     - complete_energies['oh_lamno3']: E(LaMnO3+OH) from SCF")
    print(f"  ")
    print(f"  4. Adsorption energy:")
    print(f"     - adsorption_energies['oh_lamno3']: E_ads (eV)")
    print(f"  ")
    print(f"Expected E_ads for OH/LaMnO3:")
    print(f"  Literature range: -2 to -4 eV (DFT-PBE+U)")
    print(f"  (Negative = favorable adsorption)")
    print(f"\nWorkflow details:")
    print(f"  Phase 1: Relaxation of complete system")
    print(f"    - 1 VASP relaxation (NSW=100, IBRION=2)")
    print(f"  Phase 2: Structure separation (automatic)")
    print(f"  Phase 3: SCF calculations")
    print(f"    - 3 VASP SCF (NSW=0, single-point)")
    print(f"  Total VASP jobs: 4 (1 relax + 3 SCF)")
    print(f"\nCell handling:")
    print(f"  All systems use the SAME simulation cell")
    print(f"  This eliminates basis set superposition error (BSSE)")
    print(f"\nEstimated runtime:")
    print(f"  Phase 1 (relax): ~12-18 hours (150 atoms, ionic convergence)")
    print(f"  Phase 3 (SCF): ~2-4 hours each (3 parallel jobs)")
    print(f"  Total: ~12-18 hours (relax is bottleneck, SCF runs after)")
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
