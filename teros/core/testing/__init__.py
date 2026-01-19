"""Testing utilities for PS-TEROS module development.

This module provides three levels of testing for new modules:

**Level 1: Pure Python Validation (No AiiDA Required)**

Instant validation of INCAR parameters and builder_inputs without
needing an AiiDA profile or database connection.

    >>> from teros.core.testing import validate_builder_inputs
    >>> result = validate_builder_inputs({
    ...     'parameters': {'incar': {'encut': 520, 'ibrion': 2, 'nsw': 0}}
    ... })
    >>> print(result.warnings)
    ['IBRION >= 0 requires NSW > 0 for ionic relaxation']

**Level 2: Dry-Run Input Generation (Needs AiiDA Profile)**

Generate all VASP input files (INCAR, POSCAR, KPOINTS, POTCAR) without
submitting to a cluster. Inspect files before committing to a calculation.

    >>> from teros.core.testing import dry_run_vasp
    >>> result = dry_run_vasp(structure, builder_inputs, 'vasp@localhost')
    >>> result.print_incar()  # Inspect generated INCAR
    >>> result.list_files()   # See all generated files

**Level 3: WorkGraph Structure Checks (Needs AiiDA Profile)**

Validate WorkGraph task wiring and connections before submission.
Catch orphan tasks and missing connections early.

    >>> from teros.core.testing import check_workgraph_wiring, print_workgraph_structure
    >>> wg = build_my_workgraph(...)  # Don't submit yet
    >>> print_workgraph_structure(wg)  # Visualize structure
    >>> warnings = check_workgraph_wiring(wg)  # Check for issues

**Typical Workflow:**

    # 1. Quick validation (no AiiDA needed)
    result = validate_builder_inputs(my_builder_inputs)
    if not result.is_valid:
        print(result)
        exit(1)

    # 2. Generate VASP files to inspect
    dry_result = dry_run_vasp(structure, my_builder_inputs, 'vasp@localhost')
    dry_result.print_incar()

    # 3. Check WorkGraph structure
    wg = build_my_workgraph(...)
    print_workgraph_structure(wg)
    warnings = check_workgraph_wiring(wg)

    # 4. If all looks good, submit
    wg.submit()
"""

# Level 1: Pure Python validation (always available)
from .validation import (
    ValidationResult,
    validate_incar,
    validate_builder_inputs,
    validate_structure_for_vasp,
    estimate_kpoints_mesh,
    INCAR_RULES,
)

# Level 1 extras: Pure Python INCAR/KPOINTS generation
from .dry_run import (
    generate_incar_from_dict,
    generate_kpoints_from_mesh,
)


# Level 2 & 3: AiiDA-dependent functions (lazy import)
# This allows importing validation functions without AiiDA installed

def dry_run_vasp(*args, **kwargs):
    """Generate VASP input files without running calculation.

    See :func:`teros.core.testing.dry_run.dry_run_vasp` for full documentation.
    """
    from .dry_run import dry_run_vasp as _dry_run
    return _dry_run(*args, **kwargs)


def check_workgraph_wiring(*args, **kwargs):
    """Check WorkGraph for common wiring issues.

    See :func:`teros.core.testing.workgraph_checks.check_workgraph_wiring` for full documentation.
    """
    from .workgraph_checks import check_workgraph_wiring as _check
    return _check(*args, **kwargs)


def print_workgraph_structure(*args, **kwargs):
    """Print WorkGraph structure for inspection.

    See :func:`teros.core.testing.workgraph_checks.print_workgraph_structure` for full documentation.
    """
    from .workgraph_checks import print_workgraph_structure as _print
    return _print(*args, **kwargs)


def analyze_workgraph(*args, **kwargs):
    """Perform comprehensive analysis of WorkGraph structure.

    See :func:`teros.core.testing.workgraph_checks.analyze_workgraph` for full documentation.
    """
    from .workgraph_checks import analyze_workgraph as _analyze
    return _analyze(*args, **kwargs)


def get_workgraph_execution_order(*args, **kwargs):
    """Determine task execution order (topological sort with levels).

    See :func:`teros.core.testing.workgraph_checks.get_workgraph_execution_order` for full documentation.
    """
    from .workgraph_checks import get_workgraph_execution_order as _order
    return _order(*args, **kwargs)


def estimate_workgraph_complexity(*args, **kwargs):
    """Estimate computational complexity of WorkGraph.

    See :func:`teros.core.testing.workgraph_checks.estimate_workgraph_complexity` for full documentation.
    """
    from .workgraph_checks import estimate_workgraph_complexity as _estimate
    return _estimate(*args, **kwargs)


# Classes that may be useful
from .dry_run import DryRunResult
from .workgraph_checks import WorkGraphAnalysis


__all__ = [
    # Level 1: Pure Python validation
    'ValidationResult',
    'validate_incar',
    'validate_builder_inputs',
    'validate_structure_for_vasp',
    'estimate_kpoints_mesh',
    'INCAR_RULES',
    # Level 1: Pure Python file generation
    'generate_incar_from_dict',
    'generate_kpoints_from_mesh',
    # Level 2: Dry-run (requires AiiDA)
    'dry_run_vasp',
    'DryRunResult',
    # Level 3: WorkGraph checks (requires AiiDA)
    'check_workgraph_wiring',
    'print_workgraph_structure',
    'analyze_workgraph',
    'get_workgraph_execution_order',
    'estimate_workgraph_complexity',
    'WorkGraphAnalysis',
]
