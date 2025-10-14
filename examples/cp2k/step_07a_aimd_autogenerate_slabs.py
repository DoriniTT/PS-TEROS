#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 7A: AIMD with CP2K - Auto-generate Slabs

This workflow:
1. Relaxes bulk structure (VASP)
2. Generates slabs from relaxed bulk
3. Runs AIMD on generated slabs (CP2K)

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence (equilibration + production)

Usage:
    source ~/envs/psteros/bin/activate
    python step_07a_aimd_autogenerate_slabs.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

def main():
    """Step 7A: Test AIMD with CP2K - auto-generate slabs."""

    print("\n" + "="*70)
    print("STEP 7A: AIMD WITH CP2K - AUTO-GENERATE SLABS")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, '../structures')

    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")

    # VASP parameters for bulk relaxation
    bulk_params = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # AIMD configuration
    print("\n3. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    for i, stage in enumerate(aimd_sequence):
        print(f"   Stage {i+1}: {stage['temperature']} K, {stage['steps']} steps")

    # CP2K AIMD parameters
    aimd_params = get_aimd_defaults_cp2k(
        cutoff=400,
        rel_cutoff=60,
        timestep=1.0,
        eps_scf=1e-6,
        thermostat='NOSE',
    )

    # Add KIND section for Ag and O
    if 'FORCE_EVAL' not in aimd_params:
        aimd_params['FORCE_EVAL'] = {}
    if 'SUBSYS' not in aimd_params['FORCE_EVAL']:
        aimd_params['FORCE_EVAL']['SUBSYS'] = {}

    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {
            "_": "Ag",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11",
            "POTENTIAL": "GTH-PBE-q11",
        },
        {
            "_": "O",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6",
            "POTENTIAL": "GTH-PBE-q6",
        }
    ]

    aimd_options = {
        'resources': {
            'num_machines': 2,
            'num_cores_per_machine': 128,
        },
        'queue_name': 'paralela',
    }

    print("\n4. Building workgraph...")
    print("   Workflow: Bulk relaxation → Slab generation → AIMD")
    print("   Bulk code: VASP-VTST-6.4.3@bohr")
    print("   AIMD code: CP2K-NEWCPU-2023@bohr")
    print("   Using preset: 'aimd_only'")

    # Build workgraph using preset with CP2K for AIMD
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',  # Use CP2K for AIMD

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',

        # Code for bulk relaxation (VASP)
        code_label='VASP-VTST-6.4.3@bohr',
        
        # Code for AIMD (CP2K)
        aimd_code_label='CP2K-NEWCPU-2023@bohr',
        potential_family='PBE',
        kpoints_spacing=0.4,
        clean_workdir=False,

        # Bulk parameters (REQUIRED for auto-generation)
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=bulk_params,
        bulk_options=bulk_options,

        # Slab generation (REQUIRED for auto-generation)
        miller_indices=[1, 1, 1],
        min_slab_thickness=15.0,
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,

        # AIMD (CP2K)
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        name='Step07A_AIMD_CP2K_AutoGen',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 7A SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nWorkflow stages:")
    print(f"  1. Bulk relaxation (VASP)")
    print(f"  2. Slab generation")
    print(f"  3. AIMD on all slabs (CP2K)")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures")
    print(f"  - AIMD stage outputs in tasks:")
    print(f"    wg.tasks['aimd_stage_00_300K'].outputs:")
    print(f"      - structures, remote_folders")
    print(f"      - parameters, trajectories, retrieved")
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
