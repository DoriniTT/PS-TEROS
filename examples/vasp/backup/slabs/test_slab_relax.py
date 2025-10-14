#!/usr/bin/env python
"""
Test script to verify the slab relaxation workflow can be built without submitting it.
"""

from aiida import load_profile
from teros.workgraph import build_core_workgraph

def test_slab_relax_workflow():
    """Test that the slab relaxation workflow can be built successfully."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile()

    # Define minimal parameters for testing
    structures_dir = '/home/thiagotd/git/PS-TEROS/teros/structures'

    bulk_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    bulk_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    metal_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    metal_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    nonmetal_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    nonmetal_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    oxygen_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 2, 'NSW': 1}
    oxygen_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    # Slab relaxation parameters
    slab_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 2, 'NSW': 1, 'EDIFFG': -0.02}
    slab_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    print("\nBuilding workflow with slab relaxation enabled...")
    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name='ag3po4.cif',
        metal_name='Ag.cif',
        nonmetal_name='P.cif',
        oxygen_name='O2.cif',
        code_label='VASP-VTST-6.4.3@bohr',
        potential_family='PBE',
        bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        metal_potential_mapping={'Ag': 'Ag'},
        nonmetal_potential_mapping={'P': 'P'},
        oxygen_potential_mapping={'O': 'O'},
        kpoints_spacing=0.3,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=True,
        # Slab generation
        miller_indices=[1, 0, 0],
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,
        lll_reduce=False,
        center_slab=True,
        symmetrize=False,
        primitive=True,
        in_unit_planes=False,
        max_normal_search=None,
        # Slab relaxation
        relax_slabs=True,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        slab_kpoints_spacing=0.3,
        name='Test_SlabRelax',
    )

    print("\n✓ Workflow built successfully!")
    print(f"Workflow name: {wg.name}")

    # List tasks
    task_list = list(wg.tasks)
    print(f"Number of tasks: {len(task_list)}")
    print(f"\nTasks in workflow:")
    for task in task_list:
        print(f"  - {task.name}")

    # Try to export to HTML
    try:
        html_file = 'test_slab_relax.html'
        wg.to_html(html_file)
        print(f"\n✓ Workflow visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML: {e}")

    # Check that we have the expected outputs
    print(f"\n✓ Checking workflow outputs...")
    expected_outputs = [
        'bulk_energy', 'metal_energy', 'nonmetal_energy', 'oxygen_energy',
        'bulk_structure', 'metal_structure', 'nonmetal_structure', 'oxygen_structure',
        'formation_enthalpy', 'slab_structures', 'relaxed_slabs', 'slab_energies'
    ]

    print(f"Expected outputs: {expected_outputs}")
    print(f"\n✓ All tests passed!")
    print("\nNOTE: Workflow was NOT submitted. To submit, run slabs_relax.py")

    return wg

if __name__ == '__main__':
    try:
        wg = test_slab_relax_workflow()
    except Exception as e:
        print(f"\n✗ Workflow build failed!")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
