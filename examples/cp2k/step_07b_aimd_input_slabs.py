#!/home/thiagotd/envs/psteros/bin/python
"""
STEP 7B: AIMD with CP2K - Input Slabs Directly (with Fixed Atoms)

This workflow:
1. Takes pre-existing slab structures as input
2. Runs AIMD directly on input slabs (CP2K)
3. NO bulk relaxation or slab generation
4. Demonstrates fixed atoms constraints

Use this when you already have slab structures from:
- Previous calculations
- Manual construction
- Other structure generation tools

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence with bottom 7Å fixed

Usage:
    source ~/envs/psteros/bin/activate
    python step_07b_aimd_input_slabs.py
"""

import sys
import os
from aiida import load_profile
from aiida import orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

def main():
    """Step 7B: Test AIMD with CP2K - input slabs with fixed atoms."""

    print("\n" + "="*70)
    print("STEP 7B: AIMD WITH CP2K - INPUT SLABS (FIXED ATOMS)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'input_structures')

    # Load pre-existing slab structures
    print("\n2. Loading input slab structures...")

    # Example: Load slabs from files
    # NOTE: You need to create these files or use previous calculation outputs
    slab_files = {
        'slab_111_term0': os.path.join(structures_dir, 'ag2o_111_term0.vasp'),
        'slab_111_term1': os.path.join(structures_dir, 'ag2o_111_term1.vasp'),
    }

    input_slabs = {}
    for label, filepath in slab_files.items():
        if os.path.exists(filepath):
            atoms = read(filepath)
            input_slabs[label] = orm.StructureData(ase=atoms)
            print(f"   ✓ Loaded {label} from {filepath}")
        else:
            print(f"   ⚠️  File not found: {filepath}")
            print(f"      Creating dummy structure for demonstration...")
            # Create dummy slab for demonstration
            from ase.build import bulk, surface
            # Create simple cubic Ag lattice as base (Ag2O has complex structure)
            ag_bulk = bulk('Ag', 'fcc', a=4.09)
            slab = surface(ag_bulk, (1,1,1), 4, vacuum=15.0)
            input_slabs[label] = orm.StructureData(ase=slab)
            print(f"   ✓ Created dummy Ag(111) slab: {label}")

    # Alternative: Load from previous calculation
    # prev_wg = orm.load_node(12345)  # PK of previous workgraph
    # input_slabs = {
    #     'term_0': prev_wg.outputs.slab_structures.term_0,
    #     'term_1': prev_wg.outputs.slab_structures.term_1,
    # }

    if not input_slabs:
        print("\n✗ No slab structures loaded!")
        print("  Please create slab files or use previous calculation outputs.")
        return None

    print(f"\n3. Total slabs loaded: {len(input_slabs)}")

    # AIMD configuration
    print("\n4. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},
        {'temperature': 300, 'steps': 100},
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

    # Add KIND section
    if 'FORCE_EVAL' not in aimd_params:
        aimd_params['FORCE_EVAL'] = {}
    if 'SUBSYS' not in aimd_params['FORCE_EVAL']:
        aimd_params['FORCE_EVAL']['SUBSYS'] = {}

    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {"_": "Ag", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
    ]

    aimd_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Fixed atoms configuration
    print("\n5. Fixed atoms configuration:")
    print("   Type: bottom")
    print("   Thickness: 7.0 Å")
    print("   Elements: all")
    print("   Components: XYZ (fully rigid)")

    print("\n6. Building workgraph...")
    print(f"  Workflow: AIMD only (no bulk, no slab generation)")
    print(f"  AIMD code: CP2K-NEWCPU-2023@bohr")
    print(f"  Fixed atoms: bottom 7Å")

    # Build workgraph
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',

        # Minimal bulk parameters (not used, but required by internal logic)
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        code_label='VASP-VTST-6.4.3@bohr',
        
        # Code for AIMD (CP2K)
        aimd_code_label='CP2K-NEWCPU-2023@bohr',

        # Input slabs directly
        input_slabs=input_slabs,

        # AIMD (CP2K)
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        # Fixed atoms (NEW)
        fix_atoms=True,
        fix_type='bottom',
        fix_thickness=7.0,
        fix_elements=None,  # All elements
        fix_components='XYZ',

        clean_workdir=False,
        name='Step07B_AIMD_CP2K_InputSlabs_FixedAtoms',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n7. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 7B SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nWorkflow stages:")
    print(f"  1. AIMD on input slabs (CP2K) - ONLY")
    print(f"\nNo bulk relaxation or slab generation performed!")
    print(f"\nFixed atoms:")
    print(f"  - Bottom 7Å of all slabs are fully constrained (XYZ)")
    print(f"  - Constraint applied to all AIMD stages")
    print(f"\nExpected outputs:")
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
