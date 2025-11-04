#!/home/thiagotd/envs/aiida/bin/python
"""
STEP XX: MLFF Test - Experiment with Machine Learning Force Fields

This script tests if the existing AIMD module can support MLFF workflows
without code changes, by simply adding ML-related INCAR tags.

Test Strategy:
1. Stage 0: Training phase (ML_ISTART=0) - Collect ab initio training data
2. Stage 1: Production phase (ML_ISTART=2) - Use trained ML model

We skip the refinement phase for this initial test to keep it simple.

Critical Questions This Test Answers:
- Q1: Does remote_folder restart automatically include ML_AB/ML_FFN files?
- Q2: Can ML_ISTART=2 stage successfully read and use the trained model?
- Q3: Can we reduce ENCUT in production stage?

Material: Ag (silver) - small system for fast testing
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.aimd import build_aimd_workgraph
from ase.io import read


def main():
    """Test MLFF with existing AIMD infrastructure."""

    print("\n" + "="*70)
    print("MLFF EXPERIMENTAL TEST")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Navigate to examples/vasp/structures
    structures_dir = os.path.join(script_dir, '..', '..', 'examples', 'vasp', 'structures')
    ag_cif = os.path.join(structures_dir, 'Ag.cif')

    if not os.path.exists(ag_cif):
        print(f"\n✗ Error: Structure file not found at {ag_cif}")
        print("   Please check the path or copy a suitable Ag structure.")
        sys.exit(1)

    print(f"\n2. Loading structure:")
    print(f"   File: {ag_cif}")

    # Load structure using ASE
    ag_ase = read(ag_cif)
    ag_structure = orm.StructureData(ase=ag_ase)

    print(f"   ✓ Loaded Ag structure ({len(ag_structure.sites)} atoms)")

    # Code configuration
    code_label = 'VASP6.5.0@cluster02'
    potential_family = 'PBE'

    print(f"\n3. VASP configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # MLFF stages - Testing with minimal workflow
    print("\n4. MLFF stages configuration:")

    # IMPORTANT: We're using SHORT runs for testing!
    # In production, training would need 200-500+ steps
    mlff_stages = [
        # Stage 0: Training (ab initio with ML learning)
        {
            'TEBEG': 300,
            'TEEND': 300,
            'NSW': 20,  # SHORT for testing (use 200+ in production!)

            # MLFF training parameters
            'ML_LMLFF': True,    # Enable MLFF
            'ML_ISTART': 0,      # Start training from scratch

            # ML hyperparameters (VASP defaults are usually good)
            'ML_WTOTEN': 0.1,    # Energy weight in loss function
            'ML_WTIFOR': 1.0,    # Force weight in loss function
            # 'ML_WTSIF': 0.0,   # Stress weight (default 0, not needed for slabs)
        },

        # Stage 1: Production (pure ML predictions)
        {
            'TEBEG': 300,
            'TEEND': 300,
            'NSW': 50,  # Can be much longer (1000s) since ML is fast

            # MLFF production parameters
            'ML_LMLFF': True,
            'ML_ISTART': 2,      # Use trained model, no further training

            # TEST: Can we reduce DFT parameters in ML mode?
            # If VASP still needs wavefunctions, this might fail
            # If not, this could speed up the calculation
            'ENCUT': 200,  # Reduced from training value (400)
        },
    ]

    print("   Stage 0 (Training):")
    print(f"     - Temperature: {mlff_stages[0]['TEBEG']} K")
    print(f"     - Steps: {mlff_stages[0]['NSW']} (SHORT TEST, use 200+ for real training!)")
    print(f"     - ML_ISTART: {mlff_stages[0]['ML_ISTART']} (train from scratch)")
    print(f"     - ML_WTOTEN: {mlff_stages[0]['ML_WTOTEN']}, ML_WTIFOR: {mlff_stages[0]['ML_WTIFOR']}")

    print("   Stage 1 (Production):")
    print(f"     - Temperature: {mlff_stages[1]['TEBEG']} K")
    print(f"     - Steps: {mlff_stages[1]['NSW']}")
    print(f"     - ML_ISTART: {mlff_stages[1]['ML_ISTART']} (use trained model)")
    print(f"     - ENCUT: {mlff_stages[1]['ENCUT']} (TEST: reduced from 400)")

    # Builder inputs (DFT parameters for training stage)
    builder_inputs = {
        'parameters': {
            'incar': {
                # Basic DFT settings (used in training phase)
                'PREC': 'Normal',
                'ENCUT': 400,  # Full accuracy for training
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': 'Auto',

                # AIMD settings
                'IBRION': 0,      # MD mode
                'MDALGO': 2,      # Nosé-Hoover thermostat
                'POTIM': 1.0,     # 1 fs timestep
                'SMASS': 0.0,     # Automatic Nosé mass

                # Output control
                'LWAVE': True,    # CRITICAL: Need for restart (and ML files?)
                'LCHARG': True,   # CRITICAL: Need for restart

                # Optional: ML-related output
                # 'ML_OUTBLOCK': 10,  # Output ML data every N steps
            }
        },
        'kpoints_spacing': 0.5,
        'potential_family': potential_family,
        'potential_mapping': {'Ag': 'Ag'},
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'clean_workdir': False,  # CRITICAL: Must keep workdir to preserve ML files!
    }

    print("\n5. Building MLFF workgraph...")
    print("   Testing with existing build_aimd_workgraph()...")
    print("   (No code changes - just adding ML INCAR tags)")

    # Build MLFF workgraph using EXISTING function
    # This tests if MLFF works without any code changes!
    wg = build_aimd_workgraph(
        # Single structure for testing
        structures={
            'ag_test': ag_structure,
        },

        # MLFF stages (training → production)
        aimd_stages=mlff_stages,

        # Code configuration
        code_label=code_label,

        # Builder inputs
        builder_inputs=builder_inputs,

        # Optional: Test with small supercell (more atoms = better ML training)
        # Comment out for faster testing
        # supercell_specs={
        #     'ag_test': [2, 2, 1],
        # },

        # Concurrency control
        max_concurrent_jobs=1,  # Only 1 job at a time on cluster02

        # Workgraph name
        name='MLFF_Test_Ag',
    )

    print("   ✓ WorkGraph built successfully")

    # Display expected workflow
    print("\n6. Expected workflow:")
    print("   Task 1: stage_0_aimd")
    print("           - Run DFT-MD with ML training")
    print("           - Generate ML_AB (training data)")
    print("           - Generate ML_FFN (neural network weights)")
    print("           - Save remote_folder for restart")
    print()
    print("   Task 2: stage_1_aimd")
    print("           - Restart from stage_0 remote_folder")
    print("           - TEST: Check if ML_AB/ML_FFN are automatically read")
    print("           - Run pure ML-MD (no DFT forces)")
    print("           - TEST: Check if reduced ENCUT works")

    # Submit
    print("\n7. Submitting to AiiDA daemon...")
    print(f"   ⚠️  This is an EXPERIMENTAL test!")
    print(f"   ⚠️  Monitor closely for errors!")

    # Uncomment to actually submit:
    # wg.submit(wait=False)

    # For safety, just print the command
    print(f"\n   To submit, uncomment wg.submit() in this script and run again.")
    print(f"   Or manually submit with:")
    print(f"     from aiida import orm")
    print(f"     wg = orm.load_node({wg.pk})")
    print(f"     wg.submit(wait=False)")

    print(f"\n{'='*70}")
    print("TEST SETUP COMPLETE")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")

    print(f"\nCritical things to check after running:")
    print(f"  1. Stage 0 completes successfully")
    print(f"     verdi process show <stage_0_PK>")
    print(f"     Check for ML_AB and ML_FFN in outputs")
    print()
    print(f"  2. Stage 1 reads ML files correctly")
    print(f"     Check VASP OUTCAR for 'reading ML_AB' message")
    print(f"     Check if forces come from ML model (should be much faster)")
    print()
    print(f"  3. Production stage runs successfully with reduced ENCUT")
    print(f"     Compare timing: Stage 1 should be ~10x faster than Stage 0")
    print()
    print(f"  4. Check ML model quality")
    print(f"     Look for ML_LOGFILE in calculation folders")
    print(f"     Check prediction errors (should be < 1 meV/atom)")

    print(f"\n{'='*70}")
    print("EXPECTED OUTCOMES")
    print(f"{'='*70}")
    print(f"  SUCCESS: ML_ISTART=2 runs and uses trained model")
    print(f"           → No code changes needed! Just document MLFF usage.")
    print()
    print(f"  FAILURE: Stage 1 can't find ML_AB/ML_FFN files")
    print(f"           → Need to modify restart mechanism")
    print(f"           → Explicitly copy ML files between stages")
    print()
    print(f"  FAILURE: Reduced ENCUT in Stage 1 causes crash")
    print(f"           → Document that ENCUT must stay constant")
    print(f"           → Or investigate proper ENCUT reduction procedure")

    print(f"\n{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
        print(f"✓ Test script completed successfully")
        print(f"  WorkGraph ready but NOT submitted (uncomment wg.submit() to run)")
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
