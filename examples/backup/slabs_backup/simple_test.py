#!/usr/bin/env python
"""
Simple test to validate the input_slabs parameter works correctly.
This test checks the API without building the full WorkGraph.
"""

from aiida import load_profile, orm
from ase.io import read
import inspect


def simple_test():
    """Simple API test."""
    
    print("\n" + "="*80)
    print("SIMPLE TEST: User-Provided Slabs API")
    print("="*80)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile()
    print("   ✓ Profile loaded")

    # Load slab structures
    print("\n2. Loading slab structures...")
    slabs_dir = "/home/thiagotd/git/PS-TEROS/examples/slabs/input_structures"
    
    input_slabs = {}
    for idx in range(3):
        slab_file = f"slab_term_{idx}.cif"
        slab_path = f"{slabs_dir}/{slab_file}"
        atoms = read(slab_path)
        input_slabs[f"term_{idx}"] = orm.StructureData(ase=atoms)
        print(f"   ✓ Loaded term_{idx}: {atoms.get_chemical_formula()}")
    
    print(f"\n   ✓ Total: {len(input_slabs)} slabs loaded")

    # Test imports
    print("\n3. Testing imports...")
    from teros.core.workgraph import (
        build_core_workgraph,
        build_core_workgraph_with_map,
    )
    print("   ✓ Imported workflow builders")

    # Check function signatures
    print("\n4. Checking function signatures...")
    
    sig = inspect.signature(build_core_workgraph)
    params = list(sig.parameters.keys())
    
    if 'input_slabs' in params:
        print("   ✓ 'input_slabs' parameter exists in build_core_workgraph")
        default = sig.parameters['input_slabs'].default
        print(f"   ✓ Default value: {default}")
    else:
        print("   ✗ 'input_slabs' parameter missing!")
        return False
    
    sig = inspect.signature(build_core_workgraph_with_map)
    params = list(sig.parameters.keys())
    
    if 'input_slabs' in params:
        print("   ✓ 'input_slabs' parameter exists in build_core_workgraph_with_map")
    else:
        print("   ✗ 'input_slabs' parameter missing!")
        return False

    # Test parameter binding
    print("\n5. Testing parameter binding...")
    try:
        sig = inspect.signature(build_core_workgraph_with_map)
        test_args = {
            'structures_dir': '/test',
            'bulk_name': 'test.cif',
            'metal_name': 'test.cif',
            'nonmetal_name': 'test.cif',
            'oxygen_name': 'test.cif',
            'input_slabs': input_slabs,  # Test with our loaded slabs
        }
        bound = sig.bind_partial(**test_args)
        print("   ✓ Function accepts input_slabs parameter")
        print(f"   ✓ Bound {len(bound.arguments)} arguments successfully")
    except TypeError as e:
        print(f"   ✗ Binding error: {e}")
        return False

    # Summary
    print("\n" + "="*80)
    print("API TEST RESULTS")
    print("="*80)
    print("\n✅ All API tests PASSED!")
    print("\nVerified:")
    print("  • AiiDA profile loads correctly")
    print(f"  • {len(input_slabs)} slab structures can be loaded from files")
    print("  • Workflow builders have 'input_slabs' parameter")
    print("  • Parameters can be bound correctly")
    print("\nThe API is correct. The actual workflow execution would:")
    print("  1. Skip slab generation (since input_slabs provided)")
    print("  2. Use the provided slab structures directly")
    print("  3. Relax them with VASP (if relax_slabs=True)")
    print("\nNote: Full WorkGraph building requires actual calculation submission")
    print("      or a more complex mock setup. This test validates the API layer.")
    print("\n" + "="*80 + "\n")
    
    return True


if __name__ == "__main__":
    try:
        success = simple_test()
        if success:
            print("🎉 API TEST SUCCESSFUL!\n")
            exit(0)
        else:
            print("❌ API TEST FAILED!\n")
            exit(1)
    except Exception as e:
        print(f"\n❌ TEST ERROR: {e}\n")
        import traceback
        traceback.print_exc()
        exit(1)
