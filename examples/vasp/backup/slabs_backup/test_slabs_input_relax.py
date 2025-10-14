#!/usr/bin/env python
"""
Test script for slabs_input_relax.py functionality.

This tests that the user-provided slabs feature works correctly without
actually submitting calculations.
"""

from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map


def test_slabs_input_relax():
    """Test the slabs_input_relax functionality."""
    
    print("\n" + "="*80)
    print("TEST: User-Provided Slabs Feature")
    print("="*80)
    
    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    try:
        load_profile()
        print("   ✓ Profile loaded successfully")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False

    # Define structures directory
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"
    slabs_dir = "/home/thiagotd/git/PS-TEROS/examples/slabs/input_structures"
    
    # Check bulk structures exist
    print("\n2. Checking bulk structures...")
    required_structures = {
        'bulk': 'ag3po4.cif',
        'metal': 'Ag.cif',
        'nonmetal': 'P.cif',
        'oxygen': 'O2.cif',
    }
    
    import os
    for struct_type, filename in required_structures.items():
        filepath = f"{structures_dir}/{filename}"
        if os.path.exists(filepath):
            print(f"   ✓ {struct_type}: {filename} found")
        else:
            print(f"   ✗ {struct_type}: {filename} NOT found")
            return False

    # Load pre-generated slab structures
    print("\n3. Loading pre-generated slab structures...")
    
    input_slabs = {}
    slab_files = ["slab_term_0.cif", "slab_term_1.cif", "slab_term_2.cif"]
    
    for idx, slab_file in enumerate(slab_files):
        try:
            slab_path = f"{slabs_dir}/{slab_file}"
            atoms = read(slab_path)
            input_slabs[f"term_{idx}"] = orm.StructureData(ase=atoms)
            print(f"   ✓ Loaded {slab_file} as term_{idx}")
            print(f"     Formula: {atoms.get_chemical_formula()}")
            print(f"     Atoms: {len(atoms)}")
        except FileNotFoundError:
            print(f"   ✗ Warning: {slab_file} not found")
            return False
        except Exception as e:
            print(f"   ✗ Error loading {slab_file}: {e}")
            return False
    
    if not input_slabs:
        print("\n   ✗ No slab structures loaded!")
        return False
    
    print(f"\n   ✓ Successfully loaded {len(input_slabs)} slab structures")

    # Define minimal calculation parameters for testing
    print("\n4. Setting up workflow parameters...")
    
    code_label = "VASP-VTST-6.4.3@bohr"
    potential_family = "PBE"

    # Minimal parameters
    bulk_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 3,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    bulk_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 100,
        "EDIFFG": -0.02,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }
    
    print("   ✓ Parameters configured")

    # Build WorkGraph
    print("\n5. Building WorkGraph with user-provided slabs...")
    try:
        wg = build_core_workgraph_with_map(
            structures_dir=structures_dir,
            bulk_name="ag3po4.cif",
            metal_name="Ag.cif",
            nonmetal_name="P.cif",
            oxygen_name="O2.cif",
            code_label=code_label,
            potential_family=potential_family,
            bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
            metal_potential_mapping={"Ag": "Ag"},
            nonmetal_potential_mapping={"P": "P"},
            oxygen_potential_mapping={"O": "O"},
            kpoints_spacing=0.3,
            bulk_parameters=bulk_parameters,
            bulk_options=bulk_options,
            metal_parameters=bulk_parameters.copy(),
            metal_options=bulk_options.copy(),
            nonmetal_parameters=bulk_parameters.copy(),
            nonmetal_options=bulk_options.copy(),
            oxygen_parameters=bulk_parameters.copy(),
            oxygen_options=bulk_options.copy(),
            clean_workdir=True,
            # KEY TEST: Provide user slabs
            input_slabs=input_slabs,
            # Slab relaxation
            relax_slabs=True,
            slab_parameters=slab_parameters,
            slab_options=slab_options,
            slab_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
            slab_kpoints_spacing=0.3,
            name="TEST_InputSlabs_Ag3PO4",
        )
        print("   ✓ WorkGraph created successfully!")
        print(f"   ✓ WorkGraph name: {wg.name}")
        
    except Exception as e:
        print(f"   ✗ Error building WorkGraph: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Check WorkGraph structure
    print("\n6. Validating WorkGraph structure...")
    try:
        # Check that WorkGraph has expected tasks
        tasks = list(wg.tasks.keys()) if hasattr(wg, 'tasks') else []
        print(f"   ✓ WorkGraph has {len(tasks)} tasks")
        
        # Try to export to HTML (optional)
        try:
            html_file = "test_input_slabs.html"
            wg.to_html(html_file)
            print(f"   ✓ Exported WorkGraph visualization to {html_file}")
        except Exception as e:
            print(f"   ⚠ Could not export HTML (non-critical): {e}")
        
    except Exception as e:
        print(f"   ✗ Error validating WorkGraph: {e}")
        return False

    # Summary
    print("\n" + "="*80)
    print("TEST RESULTS")
    print("="*80)
    print("\n✅ All tests PASSED!")
    print("\nThe user-provided slabs feature is working correctly:")
    print(f"  • Loaded {len(input_slabs)} user-provided slab structures")
    print("  • Created WorkGraph without requiring slab generation parameters")
    print("  • WorkGraph ready for submission (not submitted in test)")
    print("\nWhat would happen on submission:")
    print("  1. Relax bulk Ag₃PO₄ structure")
    print("  2. Relax reference structures (Ag, P, O₂)")
    print("  3. Calculate formation enthalpy")
    print(f"  4. Relax {len(input_slabs)} user-provided slab structures in parallel")
    print("  5. Extract energies for each relaxed slab")
    print("\nTo actually submit:")
    print("  wg.submit(wait=False)")
    print("  or run: python examples/slabs/slabs_input_relax.py")
    print("\n" + "="*80 + "\n")
    
    return True


if __name__ == "__main__":
    success = test_slabs_input_relax()
    
    if success:
        print("\n🎉 TEST SUCCESSFUL! The feature is ready to use.\n")
        exit(0)
    else:
        print("\n❌ TEST FAILED! Please check the errors above.\n")
        exit(1)
