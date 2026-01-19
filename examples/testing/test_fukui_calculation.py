#!/usr/bin/env python
"""
Testing Workflow Demonstration: Fukui Calculation

This script demonstrates how to use the testing framework before
submitting a real Fukui calculation. It shows the 3-level testing approach:

Level 1: Pure Python validation (no AiiDA needed)
Level 2: WorkGraph structure check (needs AiiDA)
Level 3: Dry-run input generation (optional, needs AiiDA + code)

Usage:
    # Level 1 only (instant, no AiiDA):
    python test_fukui_calculation.py --level 1

    # Levels 1 + 2 (needs AiiDA profile):
    python test_fukui_calculation.py --level 2

    # All levels (needs AiiDA + VASP code):
    python test_fukui_calculation.py --level 3 --code-label VASP-6.5.1@localwork
"""

import sys
import argparse
from pathlib import Path


def run_level1_validation(builder_inputs: dict) -> bool:
    """
    Level 1: Pure Python Validation (No AiiDA Required)

    This level validates:
    - INCAR parameter consistency
    - Physical plausibility of values
    - Required keys present
    - K-points spacing appropriate

    Returns True if validation passes (no errors).
    """
    print("\n" + "=" * 70)
    print("LEVEL 1: Pure Python Validation")
    print("=" * 70)

    from teros.core.testing import validate_builder_inputs, validate_incar

    # Validate full builder_inputs
    result = validate_builder_inputs(builder_inputs)

    # Print results
    if result.errors:
        print("\nERRORS (must fix before running):")
        for err in result.errors:
            print(f"  [X] {err}")

    if result.warnings:
        print("\nWARNINGS (review before running):")
        for warn in result.warnings:
            print(f"  [!] {warn}")

    if result.info:
        print("\nINFO:")
        for key, value in result.info.items():
            print(f"  {key}: {value}")

    # Summary
    if result.is_valid:
        print("\n[OK] Level 1 validation PASSED")
        if result.warnings:
            print(f"     ({len(result.warnings)} warnings to review)")
    else:
        print("\n[FAIL] Level 1 validation FAILED")
        print("     Fix errors before proceeding")

    return result.is_valid


def run_level2_workgraph_check(structure, builder_inputs: dict, nelect: int, code_label: str = None) -> bool:
    """
    Level 2: WorkGraph Structure Check (Needs AiiDA)

    This level validates:
    - WorkGraph builds without errors
    - All tasks are properly connected
    - No orphan tasks
    - Execution order is valid

    Returns True if WorkGraph structure is valid.
    """
    print("\n" + "=" * 70)
    print("LEVEL 2: WorkGraph Structure Check")
    print("=" * 70)

    from aiida import orm
    from aiida.orm import QueryBuilder
    from teros.core.fukui import build_fukui_workgraph
    from teros.core.testing import (
        check_workgraph_wiring,
        print_workgraph_structure,
        estimate_workgraph_complexity,
        get_workgraph_execution_order,
    )

    # Convert structure if needed
    if not isinstance(structure, orm.StructureData):
        structure_node = orm.StructureData(pymatgen=structure)
    else:
        structure_node = structure

    # Find a valid code if not specified
    if not code_label:
        # Try to find any VASP code
        qb = QueryBuilder()
        qb.append(orm.Code, filters={'label': {'ilike': '%vasp%'}})
        codes = qb.all(flat=True)
        if codes:
            code_label = codes[0].full_label
            print(f"\nUsing available code: {code_label}")
        else:
            print("\n[SKIP] No VASP code found for WorkGraph building")
            print("       Level 2 requires a configured VASP code")
            return True  # Not a failure, just skip

    # Build WorkGraph (but don't submit)
    print("\nBuilding Fukui WorkGraph...")

    try:
        wg = build_fukui_workgraph(
            structure=structure_node,
            nelect_neutral=nelect,
            delta_n_values=[0.0, 0.05, 0.10, 0.15],
            code_label=code_label,
            builder_inputs=builder_inputs,
            fukui_type='plus',
            compute_fukui=True,
            name='test_fukui_workgraph',
        )
        print(f"  WorkGraph created: {wg.name}")
    except Exception as e:
        print(f"\n[FAIL] Failed to build WorkGraph: {e}")
        return False

    # Print structure
    print_workgraph_structure(wg, show_connections=True)

    # Check wiring
    warnings = check_workgraph_wiring(wg)

    if warnings:
        print("\nWiring warnings:")
        for w in warnings:
            print(f"  [!] {w}")

    # Execution order
    order = get_workgraph_execution_order(wg)
    if order:
        print("\nExecution order (tasks in same level run in parallel):")
        for i, level in enumerate(order):
            print(f"  Level {i}: {level}")

    # Complexity estimate
    complexity = estimate_workgraph_complexity(wg)
    print("\nComplexity estimate:")
    print(f"  Total tasks: {complexity['total_tasks']}")
    print(f"  VASP jobs: {complexity['vasp_tasks']}")
    print(f"  Max parallel: {complexity['max_parallel']}")
    print(f"  Execution levels: {complexity['execution_levels']}")

    # Summary
    if not warnings:
        print("\n[OK] Level 2 WorkGraph check PASSED")
        return True
    else:
        print(f"\n[WARN] Level 2 completed with {len(warnings)} warnings")
        return True  # Warnings are not fatal


def run_level3_dry_run(structure, builder_inputs: dict, code_label: str) -> bool:
    """
    Level 3: Dry-Run Input Generation (Needs AiiDA + Code)

    This level generates:
    - INCAR file
    - POSCAR file
    - KPOINTS file
    - (POTCAR if available)

    Allows inspection of actual VASP inputs before submission.

    Returns True if dry-run completes successfully.
    """
    print("\n" + "=" * 70)
    print("LEVEL 3: Dry-Run Input Generation")
    print("=" * 70)

    from teros.core.testing import dry_run_vasp

    output_dir = Path('./test_dry_run_output')

    print(f"\nGenerating VASP inputs with code: {code_label}")
    print(f"Output directory: {output_dir}")
    print("This may take a few seconds...")

    try:
        result = dry_run_vasp(
            structure=structure,
            builder_inputs=builder_inputs,
            code_label=code_label,
            output_dir=output_dir,
            potential_family=builder_inputs.get('potential_family', 'PBE'),
            potential_mapping=builder_inputs.get('potential_mapping'),
        )

        if result.success:
            print(f"\n[OK] Dry-run completed!")
            result.list_files()

            # Show INCAR
            if result.incar_path.exists():
                print("\nGenerated INCAR:")
                print("-" * 40)
                content = result.incar_path.read_text()
                # Show first 30 lines
                lines = content.split('\n')[:30]
                for line in lines:
                    print(f"  {line}")
                if len(content.split('\n')) > 30:
                    print(f"  ... ({len(content.split(chr(10))) - 30} more lines)")
                print("-" * 40)

            return True
        else:
            print(f"\n[FAIL] Dry-run failed: {result.error}")
            return False

    except Exception as e:
        print(f"\n[FAIL] Dry-run error: {e}")
        print("\nNote: Dry-run requires a valid AiiDA code configuration.")
        print("If the code doesn't exist, Level 2 checks are sufficient.")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Test Fukui calculation inputs before submission',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--level', type=int, default=1, choices=[1, 2, 3],
        help='Testing level: 1=validation, 2=+workgraph, 3=+dry-run (default: 1)',
    )
    parser.add_argument(
        '--code-label', type=str, default=None,
        help='AiiDA code label for Level 3 (e.g., VASP-6.5.1@localwork)',
    )
    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("FUKUI CALCULATION TESTING WORKFLOW")
    print("=" * 70)
    print(f"\nRunning tests up to Level {args.level}")

    # =========================================================
    # Define the calculation parameters to test
    # (These are extracted from run_fukui_plus.py)
    # =========================================================

    builder_inputs = {
        'parameters': {
            'incar': {
                # Precision and cutoff
                'prec': 'Accurate',
                'encut': 500,

                # SCF convergence
                'algo': 'All',
                'nelm': 300,
                'ediff': 1e-6,

                # Charge density
                'icharg': 2,

                # Mixing parameters
                'amix': 0.2,
                'bmix': 0.0001,
                'amix_mag': 0.8,
                'bmix_mag': 0.0001,
                'lmaxmix': 4,

                # Smearing
                'ismear': 0,
                'sigma': 0.05,

                # Real-space projection
                'lreal': 'Auto',

                # Spin settings
                'ispin': 2,

                # Parallelization
                'ncore': 2,
                'kpar': 1,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
        },
        'kpoints_spacing': 0.03,
        'potential_family': 'PBE',
        'potential_mapping': {'Sn': 'Sn_d', 'O': 'O'},
    }

    # Example structure (simple SnO2 for testing)
    from pymatgen.core import Structure, Lattice

    # Create a simple SnO2 rutile structure for testing
    # (Real calculation would load from file)
    a, c = 4.737, 3.186
    lattice = Lattice.tetragonal(a, c)
    structure = Structure(
        lattice,
        ['Sn', 'Sn', 'O', 'O', 'O', 'O'],
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.306, 0.306, 0.0],
            [0.694, 0.694, 0.0],
            [0.806, 0.194, 0.5],
            [0.194, 0.806, 0.5],
        ]
    )

    print(f"\nTest structure: {structure.composition.reduced_formula}")
    print(f"  Atoms: {len(structure)}")
    print(f"  Volume: {structure.volume:.2f} A^3")

    # Estimated NELECT (2 Sn * 14 + 4 O * 6 = 52)
    nelect = 52
    print(f"  NELECT (estimated): {nelect}")

    # =========================================================
    # Run tests
    # =========================================================

    all_passed = True

    # Level 1: Pure Python validation
    if not run_level1_validation(builder_inputs):
        all_passed = False
        if args.level == 1:
            print("\n[STOP] Fix Level 1 errors before proceeding")
            return 1

    # Level 2: WorkGraph structure check
    if args.level >= 2:
        try:
            from aiida import load_profile
            load_profile()
            if not run_level2_workgraph_check(structure, builder_inputs, nelect, args.code_label):
                all_passed = False
        except Exception as e:
            print(f"\n[SKIP] Level 2 requires AiiDA: {e}")
            all_passed = False

    # Level 3: Dry-run
    if args.level >= 3:
        if not args.code_label:
            print("\n[SKIP] Level 3 requires --code-label")
        else:
            if not run_level3_dry_run(structure, builder_inputs, args.code_label):
                all_passed = False

    # =========================================================
    # Final summary
    # =========================================================
    print("\n" + "=" * 70)
    print("TESTING SUMMARY")
    print("=" * 70)

    if all_passed:
        print("\n[OK] All tests passed!")
        print("\nYou can now submit the calculation with confidence:")
        print("  python run_fukui_plus.py")
    else:
        print("\n[!] Some tests had warnings or failures")
        print("    Review the output above before submitting")

    print("\n" + "=" * 70)
    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
