#!/usr/bin/env python
"""
Test script demonstrating the corrected Map usage with task outputs.

Based on forum response: use map_zone.item.value instead of map_zone.item
"""

from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, task, Map
from aiida.plugins import WorkflowFactory

load_profile()

print("\n" + "="*80)
print("TESTING MAP WITH TASK OUTPUTS")
print("="*80)

# Define tasks
@task.calcfunction
def generate_structures():
    """Generate multiple structures dynamically (simulates get_slabs)."""
    from ase.build import bulk

    structures = {}
    for i in range(3):
        atoms = bulk('Cu', 'fcc', a=3.6 + i*0.1)
        structures[f'struct_{i}'] = orm.StructureData(ase=atoms)

    return orm.Dict(dict=structures)


@task
def get_num_atoms(structure):
    """Count atoms in a structure."""
    atoms = structure.get_ase()
    return orm.Int(len(atoms))


print("\n" + "-"*80)
print("Test 1: Map with Static Data (baseline)")
print("-"*80)

# Test 1: Map with static data (this should work)
static_dict = {}
for i in range(3):
    from ase.build import bulk
    atoms = bulk('Cu', 'fcc', a=3.6 + i*0.1)
    static_dict[f'struct_{i}'] = orm.StructureData(ase=atoms)

with WorkGraph('test_static_map') as wg1:
    with Map(static_dict) as map_zone:
        # Use map_zone.item.value to access each structure
        num_atoms = get_num_atoms(structure=map_zone.item.value).result

    wg1.ctx.num_atoms = num_atoms
    wg1.outputs.num_atoms = wg1.ctx.num_atoms

wg1.to_html('test_static_map.html')
print(f"✓ Static Map WorkGraph created")
print(f"  Visualization: test_static_map.html")

# Submit and check
print(f"\nSubmitting static map workflow...")
wg1.submit(wait=True)
print(f"✓ Static Map completed: PK {wg1.pk}")
print(f"  Results: {wg1.outputs.num_atoms.value}")


print("\n" + "-"*80)
print("Test 2: Map with Task Outputs (our actual use case)")
print("-"*80)

# Test 2: Map with task outputs (this is what we need)
with WorkGraph('test_dynamic_map') as wg2:
    # Generate structures dynamically
    gen_task = wg2.add_task(
        generate_structures,
        name='generate_structures'
    )

    # Map over task outputs - use map_zone.item.value
    with Map(gen_task.outputs.result) as map_zone:
        # Use map_zone.item.value to access each structure
        num_atoms = get_num_atoms(structure=map_zone.item.value).result

    wg2.ctx.num_atoms = num_atoms
    wg2.outputs.num_atoms = wg2.ctx.num_atoms

wg2.to_html('test_dynamic_map.html')
print(f"✓ Dynamic Map WorkGraph created")
print(f"  Visualization: test_dynamic_map.html")

# Submit and check
print(f"\nSubmitting dynamic map workflow...")
wg2.submit(wait=True)
print(f"✓ Dynamic Map completed: PK {wg2.pk}")
print(f"  Results: {wg2.outputs.num_atoms.value}")


print("\n" + "="*80)
print("SUCCESS: Map with task outputs works with map_zone.item.value!")
print("="*80)
