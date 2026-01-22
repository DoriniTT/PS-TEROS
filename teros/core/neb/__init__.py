"""NEB (Nudged Elastic Band) calculation module for PS-TEROS.

This module provides tools for calculating minimum energy paths (MEP) and
activation barriers using the Nudged Elastic Band method with VASP.

The workflow supports:
- Optional endpoint relaxation before NEB
- IDPP or linear interpolation for intermediate images
- Two-stage NEB workflow (NEB → CI-NEB) for accurate barriers
- Single-stage NEB for quick estimates

Key features:
- Uses aiida-vasp NEB plugin (vasp.neb)
- Follows PS-TEROS architectural patterns
- Automatic barrier calculation from energy profile
- Supports VTST optimizers (LBFGS, FIRE)

Public API:
    build_neb_workgraph: Main workflow builder
    get_neb_results: Extract results from completed workflow
    print_neb_summary: Print formatted summary
    interpolate_structures: Generate NEB images via IDPP/linear interpolation
    calculate_barrier: Compute forward/reverse barriers from energies
    get_neb_incar_parameters: Get NEB-specific INCAR defaults

Example (basic NEB):
    >>> from teros.core.neb import build_neb_workgraph, print_neb_summary
    >>> wg = build_neb_workgraph(
    ...     initial_structure=initial,
    ...     final_structure=final,
    ...     n_images=5,
    ...     code_label='VASP-6.4.3@bohr',
    ...     builder_inputs={
    ...         'parameters': {'incar': {'encut': 520, 'ismear': 0}},
    ...         'options': {
    ...             'resources': {'num_machines': 3, 'num_cores_per_machine': 40},
    ...             'queue_name': 'par120',
    ...         },
    ...         'kpoints_spacing': 0.03,
    ...         'potential_family': 'PBE',
    ...         'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    ...     },
    ... )
    >>> wg.submit(wait=False)
    >>> # After completion:
    >>> print_neb_summary(wg.pk)

Example (quick estimate without endpoint relaxation):
    >>> wg = build_neb_workgraph(
    ...     initial_structure=relaxed_initial,
    ...     final_structure=relaxed_final,
    ...     n_images=3,
    ...     code_label='VASP-6.4.3@bohr',
    ...     builder_inputs={...},
    ...     relax_endpoints=False,
    ...     climb=False,  # Single NEB stage
    ... )

Example (access results):
    >>> from teros.core.neb import get_neb_results
    >>> results = get_neb_results(wg.pk)
    >>> print(f"Forward barrier: {results['barrier']['forward_barrier']:.3f} eV")
    >>> print(f"Reaction energy: {results['barrier']['reaction_energy']:.3f} eV")

Notes:
    - VASP must be compiled with VTST extensions for IOPT optimizer support
    - CI-NEB (climb=True) provides more accurate saddle point location
    - IDPP interpolation generally provides better initial paths than linear
    - For large barriers (>1 eV), consider using more images (7-9)
"""

from .workgraph import (
    build_neb_workgraph,
    get_neb_results,
    print_neb_summary,
)

from .tasks import (
    interpolate_structures,
    create_single_neb_image,
    extract_neb_energies,
    calculate_barrier,
    extract_neb_trajectory,
    extract_total_energy,
    create_neb_summary,
)

from .utils import (
    get_neb_incar_parameters,
    validate_neb_structures,
    compute_reaction_coordinate,
    format_neb_results,
    estimate_neb_resources,
    get_endpoint_relax_incar,
    get_stage2_parameters,
    DEFAULT_NEB_INCAR,
)

__all__ = [
    # Main workflow
    'build_neb_workgraph',
    'get_neb_results',
    'print_neb_summary',

    # Calcfunctions
    'interpolate_structures',
    'create_single_neb_image',
    'extract_neb_energies',
    'calculate_barrier',
    'extract_neb_trajectory',
    'extract_total_energy',
    'create_neb_summary',

    # Utilities
    'get_neb_incar_parameters',
    'validate_neb_structures',
    'compute_reaction_coordinate',
    'format_neb_results',
    'estimate_neb_resources',
    'get_endpoint_relax_incar',
    'get_stage2_parameters',

    # Constants
    'DEFAULT_NEB_INCAR',
]
