#!/home/thiagotd/envs/psteros/bin/python
"""
Test script for the default_builders module.

This script tests the functionality of the default builders module without
running actual VASP calculations. It validates:
1. Default parameters can be loaded
2. Parameters can be overridden
3. Deep merging works correctly
4. All required keys are present

Usage:
    source ~/envs/psteros/bin/activate && python test_default_builders.py
"""

from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults, update_builder_params


def test_basic_defaults():
    """Test that basic defaults can be loaded."""
    print("\n" + "="*80)
    print("TEST 1: Basic Defaults Loading")
    print("="*80)
    
    defaults = get_ag3po4_defaults(
        structures_dir="/path/to/structures",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE"
    )
    
    # Check required keys exist
    required_keys = [
        'structures_dir', 'bulk_name', 'metal_name', 'nonmetal_name', 'oxygen_name',
        'code_label', 'potential_family',
        'bulk_parameters', 'bulk_options',
        'metal_parameters', 'metal_options',
        'nonmetal_parameters', 'nonmetal_options',
        'oxygen_parameters', 'oxygen_options',
        'slab_parameters', 'slab_options',
        'kpoints_spacing', 'slab_kpoints_spacing',
        'bulk_potential_mapping', 'metal_potential_mapping',
        'nonmetal_potential_mapping', 'oxygen_potential_mapping',
        'slab_potential_mapping',
        'miller_indices', 'min_slab_thickness', 'min_vacuum_thickness',  # Slab generation params
        'clean_workdir', 'relax_slabs', 'compute_thermodynamics',
        'thermodynamics_sampling'
    ]
    
    missing_keys = [key for key in required_keys if key not in defaults]
    
    if missing_keys:
        print(f"✗ FAILED: Missing keys: {missing_keys}")
        return False
    else:
        print(f"✓ PASSED: All {len(required_keys)} required keys present")
        
    # Check some default values
    assert defaults['bulk_parameters']['ENCUT'] == 520, "Default ENCUT should be 520"
    assert defaults['bulk_parameters']['PREC'] == "Accurate", "Default PREC should be Accurate"
    assert defaults['kpoints_spacing'] == 0.3, "Default kpoints_spacing should be 0.3"
    assert defaults['miller_indices'] == [1, 0, 0], "Default miller_indices should be [1, 0, 0]"
    assert defaults['min_slab_thickness'] == 10.0, "Default min_slab_thickness should be 10.0"
    assert defaults['min_vacuum_thickness'] == 15.0, "Default min_vacuum_thickness should be 15.0"
    
    print(f"✓ PASSED: Default values are correct")
    print(f"  - bulk_parameters['ENCUT'] = {defaults['bulk_parameters']['ENCUT']}")
    print(f"  - kpoints_spacing = {defaults['kpoints_spacing']}")
    
    return True


def test_override_at_creation():
    """Test that parameters can be overridden at creation."""
    print("\n" + "="*80)
    print("TEST 2: Override at Creation")
    print("="*80)
    
    defaults = get_ag3po4_defaults(
        structures_dir="/path/to/structures",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
        bulk_parameters={'ENCUT': 600},
        kpoints_spacing=0.25
    )
    
    # Check that overrides worked
    if defaults['bulk_parameters']['ENCUT'] != 600:
        print(f"✗ FAILED: ENCUT should be 600, got {defaults['bulk_parameters']['ENCUT']}")
        return False
    
    if defaults['kpoints_spacing'] != 0.25:
        print(f"✗ FAILED: kpoints_spacing should be 0.25, got {defaults['kpoints_spacing']}")
        return False
    
    # Check that other parameters were not lost
    if 'EDIFF' not in defaults['bulk_parameters']:
        print("✗ FAILED: EDIFF was lost after override")
        return False
    
    print("✓ PASSED: Overrides applied correctly")
    print(f"  - bulk_parameters['ENCUT'] = {defaults['bulk_parameters']['ENCUT']} (overridden)")
    print(f"  - bulk_parameters['EDIFF'] = {defaults['bulk_parameters']['EDIFF']} (preserved)")
    print(f"  - kpoints_spacing = {defaults['kpoints_spacing']} (overridden)")
    
    return True


def test_deep_merge():
    """Test that deep merging works correctly."""
    print("\n" + "="*80)
    print("TEST 3: Deep Merge with update_builder_params")
    print("="*80)
    
    defaults = get_ag3po4_defaults(
        structures_dir="/path/to/structures",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE"
    )
    
    # Apply overrides using update_builder_params
    overrides = {
        'bulk_parameters': {'ENCUT': 650, 'EDIFF': 1e-7},
        'slab_options': {
            'resources': {
                'num_machines': 2
            }
        }
    }
    
    updated = update_builder_params(defaults, overrides)
    
    # Check that overrides worked
    if updated['bulk_parameters']['ENCUT'] != 650:
        print(f"✗ FAILED: ENCUT should be 650, got {updated['bulk_parameters']['ENCUT']}")
        return False
    
    if updated['bulk_parameters']['EDIFF'] != 1e-7:
        print(f"✗ FAILED: EDIFF should be 1e-7, got {updated['bulk_parameters']['EDIFF']}")
        return False
    
    if updated['slab_options']['resources']['num_machines'] != 2:
        print(f"✗ FAILED: num_machines should be 2")
        return False
    
    # Check that other parameters were preserved
    if updated['bulk_parameters']['PREC'] != "Accurate":
        print(f"✗ FAILED: PREC was lost after override")
        return False
    
    if updated['slab_options']['queue_name'] != "par40":
        print(f"✗ FAILED: queue_name was lost after override")
        return False
    
    if updated['slab_options']['resources']['num_cores_per_machine'] != 40:
        print(f"✗ FAILED: num_cores_per_machine was lost after deep override")
        return False
    
    print("✓ PASSED: Deep merge works correctly")
    print(f"  - bulk_parameters['ENCUT'] = {updated['bulk_parameters']['ENCUT']} (overridden)")
    print(f"  - bulk_parameters['EDIFF'] = {updated['bulk_parameters']['EDIFF']} (overridden)")
    print(f"  - bulk_parameters['PREC'] = {updated['bulk_parameters']['PREC']} (preserved)")
    print(f"  - slab_options['resources']['num_machines'] = {updated['slab_options']['resources']['num_machines']} (overridden)")
    print(f"  - slab_options['resources']['num_cores_per_machine'] = {updated['slab_options']['resources']['num_cores_per_machine']} (preserved)")
    print(f"  - slab_options['queue_name'] = {updated['slab_options']['queue_name']} (preserved)")
    
    return True


def test_structure_info():
    """Test that structure information is correct."""
    print("\n" + "="*80)
    print("TEST 4: Structure Information")
    print("="*80)
    
    defaults = get_ag3po4_defaults(
        structures_dir="/path/to/structures",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE"
    )
    
    # Check structure file names
    expected_files = {
        'bulk_name': 'ag3po4.cif',
        'metal_name': 'Ag.cif',
        'nonmetal_name': 'P.cif',
        'oxygen_name': 'O2.cif'
    }
    
    all_correct = True
    for key, expected_value in expected_files.items():
        if defaults[key] != expected_value:
            print(f"✗ FAILED: {key} should be '{expected_value}', got '{defaults[key]}'")
            all_correct = False
        else:
            print(f"✓ {key} = {defaults[key]}")
    
    # Check potential mappings
    if defaults['bulk_potential_mapping'] != {"Ag": "Ag", "P": "P", "O": "O"}:
        print(f"✗ FAILED: bulk_potential_mapping incorrect")
        all_correct = False
    else:
        print(f"✓ bulk_potential_mapping = {defaults['bulk_potential_mapping']}")
    
    if all_correct:
        print("\n✓ PASSED: All structure information is correct")
        return True
    else:
        print("\n✗ FAILED: Some structure information is incorrect")
        return False


def test_unpacking():
    """Test that defaults can be unpacked for function calls."""
    print("\n" + "="*80)
    print("TEST 5: Dictionary Unpacking")
    print("="*80)
    
    defaults = get_ag3po4_defaults(
        structures_dir="/path/to/structures",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE"
    )
    
    # Simulate unpacking into a function
    def mock_build_workgraph(**kwargs):
        """Mock function to test unpacking."""
        required = ['structures_dir', 'code_label', 'potential_family', 'bulk_parameters']
        for key in required:
            if key not in kwargs:
                return False, f"Missing required key: {key}"
        return True, "All required keys present"
    
    success, message = mock_build_workgraph(**defaults)
    
    if success:
        print(f"✓ PASSED: {message}")
        print(f"  - Successfully unpacked {len(defaults)} parameters")
        return True
    else:
        print(f"✗ FAILED: {message}")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "="*78 + "╗")
    print("║" + " "*20 + "DEFAULT BUILDERS MODULE TESTS" + " "*28 + "║")
    print("╚" + "="*78 + "╝")
    
    tests = [
        test_basic_defaults,
        test_override_at_creation,
        test_deep_merge,
        test_structure_info,
        test_unpacking,
    ]
    
    results = []
    for test in tests:
        try:
            results.append(test())
        except Exception as e:
            print(f"\n✗ EXCEPTION in {test.__name__}: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)
    
    # Print summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    passed = sum(results)
    total = len(results)
    
    for i, (test, result) in enumerate(zip(tests, results), 1):
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"  {i}. {test.__name__}: {status}")
    
    print(f"\n  Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n  🎉 ALL TESTS PASSED! 🎉")
        print("\n" + "="*80)
        print("The default_builders module is working correctly.")
        print("You can now use it in your PS-TEROS workflows!")
        print("="*80 + "\n")
        return 0
    else:
        print("\n  ❌ SOME TESTS FAILED")
        print("\n" + "="*80)
        print("Please review the failures above and fix the issues.")
        print("="*80 + "\n")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
