"""Patches for aiida-workgraph and aiida-vasp required by the vasp_parallelization module.

This module provides functions to apply necessary patches to aiida-workgraph and
aiida-vasp packages. These patches are required for the benchmark workflow to
function correctly with non-converging calculations.

Usage:
    # Apply all patches
    from teros.core.vasp_parallelization.patches import apply_all_patches
    apply_all_patches()

    # Or run from command line
    python -m teros.core.vasp_parallelization.patches.apply_patches

    # Verify patches are applied
    from teros.core.vasp_parallelization.patches import verify_patches
    verify_patches()
"""

from .apply_patches import apply_all_patches, verify_patches

__all__ = ["apply_all_patches", "verify_patches"]
