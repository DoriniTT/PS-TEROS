#!/usr/bin/env python
"""Apply all required patches for the vasp_parallelization module.

This script applies patches to:
1. aiida-workgraph: Adds AcceptError exception and _accept_task_as_finished method
2. aiida-vasp: Fixes method reference typo in _calculation_sanity_checks

Usage:
    # From command line
    python -m teros.core.vasp_parallelization.patches.apply_patches

    # From Python
    from teros.core.vasp_parallelization.patches import apply_all_patches
    apply_all_patches()

    # Verify patches
    from teros.core.vasp_parallelization.patches import verify_patches
    verify_patches()

After applying patches, restart the AiiDA daemon:
    verdi daemon restart
"""

import sys


def apply_all_patches(force: bool = False) -> dict:
    """Apply all required patches.

    Args:
        force: If True, apply patches even if verification succeeds.

    Returns:
        Dict with patch names and whether they were applied (True) or
        already present (False).
    """
    from .patch_aiida_workgraph import apply_patch as apply_workgraph_patch
    from .patch_aiida_vasp import apply_patch as apply_vasp_patch

    results = {}

    print("=" * 60)
    print("Applying patches for vasp_parallelization module")
    print("=" * 60)

    print("\n[1/2] Patching aiida-workgraph...")
    try:
        results["aiida-workgraph"] = apply_workgraph_patch(force=force)
    except Exception as e:
        print(f"ERROR: Failed to patch aiida-workgraph: {e}")
        results["aiida-workgraph"] = None

    print("\n[2/2] Patching aiida-vasp...")
    try:
        results["aiida-vasp"] = apply_vasp_patch(force=force)
    except Exception as e:
        print(f"ERROR: Failed to patch aiida-vasp: {e}")
        results["aiida-vasp"] = None

    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    for name, applied in results.items():
        if applied is True:
            status = "APPLIED"
        elif applied is False:
            status = "already applied"
        else:
            status = "FAILED"
        print(f"  {name}: {status}")

    if any(v is None for v in results.values()):
        print("\nWARNING: Some patches failed to apply!")
        print("The vasp_parallelization module may not work correctly.")
    else:
        print("\nAll patches applied successfully!")
        print("\nIMPORTANT: Restart the AiiDA daemon to pick up changes:")
        print("  verdi daemon restart")

    return results


def verify_patches() -> bool:
    """Verify all patches are applied.

    Returns:
        True if all patches are applied.

    Raises:
        RuntimeError: If any patch is missing.
    """
    from .patch_aiida_workgraph import verify_patch as verify_workgraph
    from .patch_aiida_vasp import verify_patch as verify_vasp

    errors = []

    try:
        verify_workgraph()
        print("aiida-workgraph patch: OK")
    except RuntimeError as e:
        errors.append(str(e))

    try:
        verify_vasp()
        print("aiida-vasp patch: OK")
    except RuntimeError as e:
        errors.append(str(e))

    if errors:
        raise RuntimeError("Missing patches:\n" + "\n".join(f"  - {e}" for e in errors))

    return True


def main():
    """Main entry point for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Apply patches required for vasp_parallelization module"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force apply patches even if already applied",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Only verify patches, don't apply",
    )
    args = parser.parse_args()

    if args.verify:
        try:
            verify_patches()
            print("\nAll patches verified successfully!")
            sys.exit(0)
        except RuntimeError as e:
            print(f"\nPatch verification failed:\n{e}")
            sys.exit(1)
    else:
        results = apply_all_patches(force=args.force)
        if any(v is None for v in results.values()):
            sys.exit(1)
        sys.exit(0)


if __name__ == "__main__":
    main()
