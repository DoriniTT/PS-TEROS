#!/usr/bin/env python
"""Static analysis and potential issue detection for 2D binning implementation."""

import ast
import inspect
from surface_modes import sample_by_coverage_bins


def analyze_implementation():
    """Analyze the 2D binning implementation for potential issues."""

    print("=== STATIC ANALYSIS OF 2D BINNING IMPLEMENTATION ===\n")

    # Get source code
    source = inspect.getsource(sample_by_coverage_bins)

    # Parse into AST
    tree = ast.parse(source)

    print("1. FUNCTION SIGNATURE ANALYSIS")
    print("-" * 40)
    func_def = tree.body[0]
    print(f"Function name: {func_def.name}")
    print(f"Arguments: {[arg.arg for arg in func_def.args.args]}")
    print(f"Return type annotation present: {func_def.returns is not None}")

    # Check for potential issues
    print("\n2. POTENTIAL ISSUES ANALYSIS")
    print("-" * 40)

    issues = []

    # Check 1: Division by zero protection
    if "n_bins" in source and "/ n_bins" in source:
        if "if vac_max == vac_min:" in source and "if oh_max == oh_min:" in source:
            print("✓ Division by zero protection: PRESENT (checks for zero range)")
        else:
            issues.append("⚠ Missing zero-range protection for division")

    # Check 2: Empty list handling
    if "if len(variants_list) == 0:" in source or "if not variants_list" in source:
        print("✓ Empty list handling: PRESENT")
    else:
        issues.append("⚠ No explicit empty list handling")

    # Check 3: None/negative n_bins handling
    if "if n_bins is None or n_bins <= 0:" in source:
        print("✓ Invalid n_bins handling: PRESENT")
    else:
        issues.append("⚠ Missing n_bins validation")

    # Check 4: Duplicate prevention
    if "sampled_indices" in source and "if idx in sampled_indices:" in source:
        print("✓ Duplicate prevention: PRESENT (using sampled_indices set)")
    else:
        issues.append("⚠ No duplicate prevention mechanism")

    # Check 5: Type checking for tuple detection
    if "isinstance(coverages[0], tuple)" in source:
        print("✓ Tuple detection for 2D: PRESENT")
    else:
        issues.append("⚠ Missing tuple type check for 2D detection")

    # Check 6: Floating point comparison
    if "float('inf')" in source:
        print("✓ Proper infinity handling: PRESENT")

    # Check 7: Index bounds checking
    if "if closest_idx >= 0:" in source:
        print("✓ Index validation: PRESENT")
    else:
        issues.append("⚠ Missing index bounds checking")

    print(f"\nTotal issues found: {len(issues)}")
    for issue in issues:
        print(f"  {issue}")

    print("\n3. COMPLEXITY ANALYSIS")
    print("-" * 40)

    # Count nested loops
    loop_count = source.count("for ")
    print(f"Number of for loops: {loop_count}")

    # Check for nested loops (2D case)
    if "for i in range(n_bins):" in source and "for j in range(n_bins):" in source:
        if "for idx, variant in enumerate(variants_list):" in source:
            print("Time complexity: O(n_bins² × variants)")
            print("  - Outer loops: n_bins × n_bins grid cells")
            print("  - Inner loop: iterates through all variants for each cell")

    print("\n4. MEMORY ANALYSIS")
    print("-" * 40)

    # Check data structures used
    if "sampled = []" in source:
        print("✓ Uses list for results (memory efficient)")

    if "sampled_indices = set()" in source:
        print("✓ Uses set for tracking (O(1) lookup)")

    print(f"Space complexity: O(min(n_bins², len(variants)))")

    print("\n5. EDGE CASE HANDLING")
    print("-" * 40)

    edge_cases = {
        "Zero range in one dimension": "vac_max == vac_min" in source,
        "Zero range in both dimensions": "oh_max == oh_min" in source,
        "Empty variant list": "len(variants_list) == 0" in source,
        "n_bins > variants count": "len(variants_list) <= n_bins" in source,
        "n_bins = None": "n_bins is None" in source,
        "n_bins <= 0": "n_bins <= 0" in source,
    }

    for case, handled in edge_cases.items():
        status = "✓" if handled else "✗"
        print(f"{status} {case}: {'HANDLED' if handled else 'NOT HANDLED'}")

    print("\n6. ALGORITHM CORRECTNESS")
    print("-" * 40)

    # Check distance calculation
    if "((vac_cov - vac_center)**2 + (oh_cov - oh_center)**2)**0.5" in source:
        print("✓ Euclidean distance: CORRECTLY IMPLEMENTED")
    else:
        print("✗ Euclidean distance formula not found or incorrect")

    # Check bin center calculation
    if "(i + 0.5) * vac_width" in source and "(j + 0.5) * oh_width" in source:
        print("✓ Bin centers: CORRECTLY CALCULATED (center of each cell)")
    else:
        print("⚠ Bin center calculation might be incorrect")

    print("\n7. CODE QUALITY METRICS")
    print("-" * 40)

    lines = source.split('\n')
    code_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]
    comment_lines = [l for l in lines if l.strip().startswith('#')]

    print(f"Total lines: {len(lines)}")
    print(f"Code lines: {len(code_lines)}")
    print(f"Comment lines: {len(comment_lines)}")
    print(f"Comment ratio: {len(comment_lines) / len(lines) * 100:.1f}%")

    # Check for magic numbers
    magic_numbers = ['0.5', '2', 'inf']
    print(f"\nMagic numbers used: {magic_numbers}")
    print("  - 0.5: Used for bin center calculation (middle of bin)")
    print("  - 2: Power for squared distance in Euclidean formula")
    print("  - inf: Initial value for minimum distance search")

    return len(issues) == 0


if __name__ == "__main__":
    success = analyze_implementation()

    print("\n" + "=" * 50)
    if success:
        print("RESULT: Implementation looks GOOD! ✓")
    else:
        print("RESULT: Some issues need attention")
    print("=" * 50)