#!/usr/bin/env python
"""
Example: Test module inputs before cluster submission.

This script demonstrates the 3-level testing approach for PS-TEROS modules:

1. Pure Python validation (no AiiDA needed) - instant feedback
2. Dry-run input generation (needs AiiDA profile) - inspect VASP files
3. WorkGraph structure check (needs AiiDA profile) - verify task wiring

Usage:
    # Level 1 only (fast, no AiiDA):
    python dry_run_example.py --validate-only

    # Level 1 + 2 (generate VASP inputs):
    python dry_run_example.py --dry-run

    # Level 1 + 2 + 3 (full check):
    python dry_run_example.py --full-check

    # Show all options:
    python dry_run_example.py --help

Example output:

    ============================================================
    LEVEL 1: Pure Python Validation
    ============================================================

    WARNINGS:
      ! ISMEAR=0 with SIGMA=0.2: consider SIGMA=0.05 for insulators

    Validation passed!

    ============================================================
    LEVEL 2: Dry-Run Input Generation
    ============================================================

    Generated files in: ./dry_run_output

    Files in ./dry_run_output:
      INCAR                      245 B
      KPOINTS                     48 B
      POSCAR                     302 B
      POTCAR                    1.2 MB

    INCAR:
    ENCUT = 520
    ISMEAR = 0
    ...
"""
import argparse
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Test PS-TEROS module inputs before cluster submission',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only run pure Python validation (Level 1)',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Generate VASP inputs without submission (Levels 1+2)',
    )
    parser.add_argument(
        '--full-check',
        action='store_true',
        help='Run all checks including WorkGraph wiring (Levels 1+2+3)',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default='./dry_run_output',
        help='Directory for generated files (default: ./dry_run_output)',
    )
    parser.add_argument(
        '--code-label',
        type=str,
        default=None,
        help='AiiDA code label (e.g., VASP-6.5.1@localwork)',
    )
    parser.add_argument(
        '--show-incar',
        action='store_true',
        help='Print generated INCAR after dry-run',
    )
    args = parser.parse_args()

    # Default to validate-only if no option specified
    if not any([args.validate_only, args.dry_run, args.full_check]):
        args.validate_only = True

    # =========================================================
    # Example builder_inputs - modify for your use case
    # =========================================================
    builder_inputs = {
        'parameters': {
            'incar': {
                'encut': 520,
                'ediff': 1e-6,
                'ismear': 0,
                'sigma': 0.05,
                'ibrion': 2,
                'nsw': 100,
                'isif': 2,
                'prec': 'Accurate',
                'lreal': 'Auto',
                'lwave': False,
                'lcharg': True,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
        },
        'kpoints_spacing': 0.03,
    }

    # =========================================================
    # LEVEL 1: Pure Python validation (no AiiDA)
    # =========================================================
    print("\n" + "=" * 60)
    print("LEVEL 1: Pure Python Validation")
    print("=" * 60)

    from teros.core.testing import validate_builder_inputs

    result = validate_builder_inputs(builder_inputs)

    if result.errors:
        print("\nERRORS:")
        for err in result.errors:
            print(f"  x {err}")

    if result.warnings:
        print("\nWARNINGS:")
        for warn in result.warnings:
            print(f"  ! {warn}")

    if result.info:
        print("\nINFO:")
        for key, value in result.info.items():
            print(f"  {key}: {value}")

    if result.is_valid:
        print("\nValidation passed!")
    else:
        print("\nValidation failed - fix errors before continuing")
        return 1

    if args.validate_only:
        return 0

    # =========================================================
    # LEVEL 2: Dry-run input generation (needs AiiDA)
    # =========================================================
    print("\n" + "=" * 60)
    print("LEVEL 2: Dry-Run Input Generation")
    print("=" * 60)

    # Check for code label
    if not args.code_label:
        print("\nNo --code-label specified.")
        print("To run dry-run, provide a code label:")
        print("  python dry_run_example.py --dry-run --code-label VASP-6.5.1@localwork")
        print("\nSkipping Level 2...")

        if not args.full_check:
            return 0
    else:
        try:
            from pymatgen.core import Structure, Lattice
            from teros.core.testing import dry_run_vasp

            # Example structure - Si diamond
            structure = Structure(
                Lattice.cubic(5.43),
                ['Si', 'Si'],
                [[0, 0, 0], [0.25, 0.25, 0.25]]
            )

            print(f"\nTest structure: {structure.composition.reduced_formula}")
            print(f"  Atoms: {len(structure)}")
            print(f"  Volume: {structure.volume:.2f} A^3")

            print(f"\nRunning dry-run with code: {args.code_label}")
            print("This may take a few seconds...")

            dry_result = dry_run_vasp(
                structure=structure,
                builder_inputs=builder_inputs,
                code_label=args.code_label,
                output_dir=args.output_dir,
                potential_family='PBE',
            )

            if dry_result.success:
                print(f"\nGenerated files in: {dry_result.output_dir}")
                dry_result.list_files()

                if args.show_incar:
                    print("\nINCARQ:")
                    dry_result.print_incar()
            else:
                print(f"\nDry-run failed: {dry_result.error}")
                print("Note: Some codes may not support dry-run mode.")

        except Exception as e:
            print(f"\nError during dry-run: {e}")
            print("Make sure AiiDA is configured and the code exists.")

    if not args.full_check:
        return 0

    # =========================================================
    # LEVEL 3: WorkGraph structure check
    # =========================================================
    print("\n" + "=" * 60)
    print("LEVEL 3: WorkGraph Structure Check")
    print("=" * 60)

    if not args.code_label:
        print("\nSkipping Level 3 (no code label provided)")
        return 0

    try:
        from aiida import orm
        from teros.core.testing import (
            check_workgraph_wiring,
            print_workgraph_structure,
            estimate_workgraph_complexity,
        )

        # Try to build a simple WorkGraph for testing
        # This example uses the Fukui module if available
        try:
            from teros.core.fukui import build_fukui_workgraph

            structure_node = orm.StructureData(pymatgen=structure)
            code = orm.load_code(args.code_label)

            print("\nBuilding Fukui WorkGraph for structure check...")

            wg = build_fukui_workgraph(
                structure=structure_node,
                code=code,
                builder_inputs=builder_inputs,
                nelect_neutral=8,  # 2 Si atoms x 4 valence electrons
                delta_n_values=[0.0, 0.05],
                fukui_type='plus',
                potential_family='PBE',
            )

            print_workgraph_structure(wg, show_connections=True)

            wiring_warnings = check_workgraph_wiring(wg)
            if wiring_warnings:
                print("\nWiring warnings:")
                for w in wiring_warnings:
                    print(f"  ! {w}")
            else:
                print("\nWorkGraph wiring looks correct!")

            complexity = estimate_workgraph_complexity(wg)
            print("\nComplexity estimate:")
            print(f"  Total tasks: {complexity['total_tasks']}")
            print(f"  VASP jobs: {complexity['vasp_tasks']}")
            print(f"  Max parallel: {complexity['max_parallel']}")
            print(f"  Execution levels: {complexity['execution_levels']}")

        except ImportError:
            print("\nFukui module not available.")
            print("To test WorkGraph checks, build your own WorkGraph:")
            print("  wg = build_my_workgraph(...)")
            print("  print_workgraph_structure(wg)")

    except Exception as e:
        print(f"\nError during WorkGraph check: {e}")

    print("\n" + "=" * 60)
    print("Testing complete!")
    print("=" * 60)

    return 0


if __name__ == '__main__':
    sys.exit(main())
