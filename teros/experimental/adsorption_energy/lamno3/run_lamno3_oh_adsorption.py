#!/home/thiagotd/envs/aiida/bin/python
"""
LaMnO3 + OH Adsorption Energy Calculation with Relaxation

This script demonstrates the new relaxation + SCF workflow for adsorption energy calculations.

System:
- Substrate: LaMnO3 (100) surface (2x2x1 supercell)
- Adsorbate: OH radical
- Structure: La32Mn28O89H (complete system)

NEW Workflow (4 phases):
1. Relax complete system (LaMnO3 + OH) using vasp.v2.relax
2. Separate relaxed structure into substrate + adsorbate
3. Run 3 SCF calculations in parallel using vasp.v2.vasp:
   - Complete system (from relaxed)
   - Bare substrate (from relaxed)
   - Isolated OH molecule (from relaxed)
4. Calculate E_ads = E_complete - E_substrate - E_molecule

NOTE: This demonstrates the new builder-based API for SCF calculations.
Relaxation currently uses old-style parameters due to different input structure
of vasp.v2.relax workchain. Full builder API for relaxation is a future enhancement.

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
    structure_file = script_dir / "LaMnO3_100_A4_surface_2x2x1_OOH_relaxed.cif"

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
    code_label = 'VASP-6.5.0@lovelace-parexp'
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
        'Mn': 'Mn_pv',  # Use _pv for transition metals
        'O': 'O',
        'H': 'H',
    }

    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 48,
        },
        'queue_name': 'parexp',
    }

    kpoints_spacing = 0.3  # Angstrom^-1

    # Relaxation parameters (Phase 1: relax complete system)
    # NOTE: Using old-style dict for now - vasp.v2.relax has different input structure than vasp.v2.vasp
    # Builder-based API for relax will be added in future enhancement
    relax_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 500,           # Cutoff energy
        'EDIFF': 1e-5,          # Electronic convergence
        'ISMEAR': 0,            # Gaussian smearing
        'SIGMA': 0.05,          # Small smearing width
        'ISPIN': 2,             # Spin-polarized
        'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments (La, Mn, O, H)
        'LORBIT': 11,           # DOSCAR and PROCAR
        'IBRION': 2,            # Conjugate gradient relaxation
        'NSW': 100,             # Max ionic steps
        'ISIF': 2,              # Relax ions only (keep cell fixed)
        'LWAVE': False,         # Don't write WAVECAR
        'LCHARG': True,         # Write CHGCAR
        'NCORE': 4,             # Parallelization
        'LASPH': True,          # Non-spherical contributions
        'IVDW': 12,             # DFT-D3 van der Waals
        'IDIPOL': 3,            # Dipole correction for slab
        'LDIPOL': True,         # Enable dipole correction
    }

    # SCF inputs (Phase 3: single-point calculations)
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 500,           # Cutoff energy
            'EDIFF': 1e-5,          # Electronic convergence
            'ISMEAR': 0,            # Gaussian smearing
            'SIGMA': 0.05,          # Small smearing width
            'ISPIN': 2,             # Spin-polarized
            'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments
            'LORBIT': 11,           # DOSCAR and PROCAR
            # NSW=0 and IBRION=-1 will be set automatically by workflow
            'LWAVE': False,         # Don't write WAVECAR
            'LCHARG': True,         # Write CHGCAR
            'NCORE': 4,             # Parallelization
            'LASPH': True,          # Non-spherical contributions
            'IVDW': 12,             # DFT-D3 van der Waals
            'IDIPOL': 3,            # Dipole correction for slab
            'LDIPOL': True,         # Enable dipole correction
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        # Note: k-points spacing passed separately to workgraph
    }

    print(f"   K-points spacing: {kpoints_spacing} Å⁻¹")
    print(f"   Spin-polarized: Yes")
    print(f"   Relaxation: ENABLED (relax → separate → SCF)")
    print(f"   Relax max ionic steps: {relax_parameters['NSW']}")

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
        adsorption_parameters=relax_parameters,  # Used for relaxation (old-style, compatible with vasp.v2.relax)
        adsorption_options=common_options,
        adsorption_kpoints_spacing=kpoints_spacing,

        # NEW: Enable relaxation and use builder inputs for SCF
        relax_before_adsorption=True,
        adsorption_scf_builder_inputs=scf_builder_inputs,  # Builder inputs for SCF phase

        name='LaMnO3_OH_RelaxThenAdsorption',
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
