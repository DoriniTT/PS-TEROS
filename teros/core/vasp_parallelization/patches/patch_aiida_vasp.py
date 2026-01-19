"""Patch for aiida-vasp vasp.py to fix method reference typo.

This patch fixes a bug in VaspWorkChain._calculation_sanity_checks where
private method names (self._check_*) are used instead of the actual public
method names (self.check_*).

Usage:
    from teros.core.vasp_parallelization.patches.patch_aiida_vasp import (
        apply_patch,
        verify_patch,
    )
    apply_patch()
    verify_patch()  # Raises if patch not applied
"""

import importlib.util
import re
from pathlib import Path


def get_vasp_workchain_path() -> Path:
    """Get the path to the aiida_vasp/workchains/v2/vasp.py file."""
    spec = importlib.util.find_spec("aiida_vasp.workchains.v2.vasp")
    if spec is None or spec.origin is None:
        raise ImportError("Could not find aiida_vasp.workchains.v2.vasp")
    return Path(spec.origin)


def verify_patch() -> bool:
    """Verify the aiida-vasp patch is applied.

    Returns:
        True if patch is applied.

    Raises:
        RuntimeError: If patch is not applied.
    """
    file_path = get_vasp_workchain_path()

    with open(file_path, "r") as f:
        content = f.read()

    # Check for the buggy pattern (should NOT exist after patching)
    if "self._check_misc_output" in content:
        raise RuntimeError(
            "aiida-vasp patch not applied: self._check_misc_output found. "
            "Run: python -m teros.core.vasp_parallelization.patches.apply_patches"
        )

    # Check for the correct pattern (should exist after patching)
    if "self.check_misc_output" not in content:
        raise RuntimeError(
            "aiida-vasp patch incomplete: self.check_misc_output not found"
        )

    return True


def apply_patch(force: bool = False) -> bool:
    """Apply the patch to aiida-vasp vasp.py.

    Args:
        force: If True, apply patch even if verification succeeds.

    Returns:
        True if patch was applied, False if already patched.
    """
    # Check if already patched
    if not force:
        try:
            verify_patch()
            print("aiida-vasp patch already applied.")
            return False
        except RuntimeError:
            pass  # Needs patching

    # Get file path
    file_path = get_vasp_workchain_path()

    # Backup original file
    backup_path = file_path.with_suffix(".py.backup")
    if not backup_path.exists():
        import shutil

        shutil.copy(file_path, backup_path)
        print(f"Backed up original to: {backup_path}")

    # Read current content
    with open(file_path, "r") as f:
        content = f.read()

    # Apply the fix: replace self._check_* with self.check_*
    # Only in the _calculation_sanity_checks method context
    # The pattern is in a list like:
    # checks = [
    #     self._check_misc_output,
    #     self._check_calc_is_finished,
    #     self._check_electronic_converged,
    #     self._check_ionic_converged,
    # ]

    replacements = [
        ("self._check_misc_output", "self.check_misc_output"),
        ("self._check_calc_is_finished", "self.check_calc_is_finished"),
        ("self._check_electronic_converged", "self.check_electronic_converged"),
        ("self._check_ionic_converged", "self.check_ionic_converged"),
    ]

    patched_content = content
    changes_made = 0

    for old, new in replacements:
        if old in patched_content:
            patched_content = patched_content.replace(old, new)
            changes_made += 1
            print(f"  Replaced: {old} -> {new}")

    if changes_made == 0:
        print("No changes needed - patch may already be applied")
        return False

    # Write patched content
    with open(file_path, "w") as f:
        f.write(patched_content)

    print(f"Patched: {file_path} ({changes_made} replacements)")

    # Verify
    verify_patch()
    print("aiida-vasp patch verification: OK")

    return True


if __name__ == "__main__":
    apply_patch()
