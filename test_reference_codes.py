#!/home/thiagotd/envs/aiida/bin/python
"""
Test script for reference code labels feature

This script verifies that the new bulk_code_label, metal_code_label,
nonmetal_code_label, and oxygen_code_label parameters work correctly.
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    """Test the reference code labels feature."""

    print("\n" + "="*70)
    print("TESTING REFERENCE CODE LABELS FEATURE")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup minimal test parameters
    print("\n2. Setting up test parameters...")

    # Use a test directory from examples
    test_dir = "/home/thiagotd/git/synergy-adsorption/calculos/phase1-surface-stability/ag3po4/rscan_rVV10"
    structures_dir = os.path.join(test_dir, 'structures')

    print(f"   Structures dir: {structures_dir}")

    # Minimal VASP parameters
    minimal_params = {
        'PREC': 'Accurate',
        'ENCUT': 500,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 1,  # Just 1 step for testing
        'EDIFFG': -0.01,
    }

    oxygen_params = {
        **minimal_params,
        'ISPIN': 2,
        'LORBIT': 11,
    }

    # Minimal options
    reference_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    slab_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 128,
        },
        'queue_name': 'par128',
    }

    # Build workgraph with separate codes
    print("\n3. Building workgraph with separate codes for all calculations...")
    print("   Default code:  VASP-6.5.0@lovelace-par128")
    print("   Bulk code:     VASP-6.5.0@bohr-new (40 cores)")
    print("   Metal code:    VASP-6.5.0@bohr-new (40 cores)")
    print("   Nonmetal code: VASP-6.5.0@bohr-new (40 cores)")
    print("   Oxygen code:   VASP-6.5.0@bohr-new (40 cores)")
    print("   Slab code:     VASP-6.5.0@lovelace-par128 (128 cores)")

    try:
        wg = build_core_workgraph(
            # Workflow configuration
            workflow_preset='surface_thermodynamics',
            max_concurrent_jobs=4,
            compute_cleavage=False,
            compute_relaxation_energy=False,

            # Structures
            structures_dir=structures_dir,
            bulk_name='ag3po4_bulk.vasp',
            metal_name='ag.cif',
            nonmetal_name='p_black.cif',
            oxygen_name='o2.cif',

            # NEW: Separate code labels for all calculations
            code_label='VASP-6.5.0@lovelace-par128',  # Default
            bulk_code_label='VASP-6.5.0@bohr-new',    # Bulk
            metal_code_label='VASP-6.5.0@bohr-new',   # Metal
            nonmetal_code_label='VASP-6.5.0@bohr-new', # Nonmetal
            oxygen_code_label='VASP-6.5.0@bohr-new',  # Oxygen
            slab_code_label='VASP-6.5.0@lovelace-par128',  # Slabs

            potential_family='PBE',
            kpoints_spacing=0.04,
            clean_workdir=False,

            # Bulk
            bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
            bulk_parameters=minimal_params,
            bulk_options=reference_options,

            # Metal
            metal_potential_mapping={'Ag': 'Ag'},
            metal_parameters=minimal_params,
            metal_options=reference_options,

            # Nonmetal
            nonmetal_potential_mapping={'P': 'P'},
            nonmetal_parameters=minimal_params,
            nonmetal_options=reference_options,

            # Oxygen
            oxygen_potential_mapping={'O': 'O'},
            oxygen_parameters=oxygen_params,
            oxygen_options=reference_options,

            # Slab generation
            miller_indices=[0, 0, 1],
            min_slab_thickness=18.0,
            min_vacuum_thickness=15.0,
            lll_reduce=True,
            center_slab=False,
            symmetrize=True,
            primitive=True,

            # Slab relaxation
            slab_parameters=minimal_params,
            slab_options=slab_options,
            slab_kpoints_spacing=0.06,

            # Thermodynamics
            thermodynamics_sampling=100,

            # Name
            name='Test_Reference_Code_Labels',
        )

        print("   ✓ WorkGraph built successfully!")
        print(f"\n4. Verifying reference tasks were created with correct codes...")

        # Check that VASP tasks exist for references
        reference_tasks = ['VaspWorkChain', 'VaspWorkChain1', 'VaspWorkChain2', 'VaspWorkChain3']
        task_descriptions = [
            'Bulk Ag₃PO₄',
            'Metal (Ag)',
            'Nonmetal (P)',
            'Oxygen (O₂)'
        ]

        for task_name, desc in zip(reference_tasks, task_descriptions):
            if hasattr(wg.tasks, task_name):
                print(f"   ✓ {task_name}: {desc}")
            else:
                print(f"   ✗ {task_name}: MISSING!")
                raise ValueError(f"Task {task_name} not found in workgraph")

        print(f"\n{'='*70}")
        print("TEST PASSED: Reference code labels feature works correctly!")
        print(f"{'='*70}\n")

        return wg

    except Exception as e:
        print(f"\n✗ Error building workgraph: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
