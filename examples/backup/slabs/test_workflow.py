#!/usr/bin/env python
"""
Test script to verify the slabs workflow can be built without submitting it.
"""

from aiida import load_profile
from teros.workgraph import build_core_workgraph

def test_workflow_build():
    """Test that the workflow can be built successfully."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile()

    # Define parameters (minimal set for testing)
    structures_dir = '/home/thiagotd/git/PS-TEROS/teros/structures'

    bulk_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    bulk_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    metal_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    metal_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    nonmetal_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'NSW': 1}
    nonmetal_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    oxygen_parameters = {'ENCUT': 520, 'IBRION': 2, 'ISIF': 2, 'NSW': 1}
    oxygen_options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}

    print("\nBuilding workflow...")
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
        # Slab parameters
        miller_indices=[1, 0, 0],
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,
        lll_reduce=False,
        center_slab=True,
        symmetrize=False,
        primitive=True,
        in_unit_planes=False,
        max_normal_search=None,
        name='Test_Workflow',
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
        html_file = 'test_workflow.html'
        wg.to_html(html_file)
        print(f"\n✓ Workflow visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML: {e}")

    print("\n✓ All workflow build tests passed!")
    print("\nNOTE: Workflow was NOT submitted. To submit, run slabs.py")

    return wg

if __name__ == '__main__':
    try:
        wg = test_workflow_build()
    except Exception as e:
        print(f"\n✗ Workflow build failed!")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
