"""Input validation without running calculations.

This module provides pure Python validation utilities for VASP inputs.
No AiiDA imports are used, allowing instant validation without a database.

Example:
    >>> from teros.core.testing import validate_builder_inputs
    >>> result = validate_builder_inputs({
    ...     'parameters': {'incar': {'encut': 520, 'ibrion': 2, 'nsw': 0}}
    ... })
    >>> print(result.warnings)
    ['IBRION >= 0 requires NSW > 0 for ionic relaxation']
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Callable, Tuple


@dataclass
class ValidationResult:
    """Result of input validation.

    Attributes:
        errors: Critical issues that will cause calculation failure
        warnings: Non-critical issues that may cause problems
        info: Additional information about the validation
    """

    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    info: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_valid(self) -> bool:
        """Check if validation passed (no errors)."""
        return len(self.errors) == 0

    def __str__(self) -> str:
        """Human-readable summary."""
        lines = []
        if self.errors:
            lines.append("ERRORS:")
            for e in self.errors:
                lines.append(f"  - {e}")
        if self.warnings:
            lines.append("WARNINGS:")
            for w in self.warnings:
                lines.append(f"  - {w}")
        if not lines:
            lines.append("Validation passed!")
        return "\n".join(lines)


# Type alias for INCAR validation rules
IncarRule = Tuple[Callable[[Dict[str, Any]], bool], str]

# Common INCAR consistency rules
# Each rule is (condition_fn, warning_message)
# condition_fn returns True if there's a problem
INCAR_RULES: List[IncarRule] = [
    # Relaxation settings
    (
        lambda i: i.get('ibrion', -1) >= 0 and i.get('nsw', 0) == 0,
        "IBRION >= 0 requires NSW > 0 for ionic relaxation",
    ),
    (
        lambda i: i.get('isif', 2) == 3 and i.get('ibrion', -1) < 0,
        "ISIF=3 (cell relaxation) requires IBRION >= 0",
    ),
    (
        lambda i: i.get('isif', 2) in (3, 4, 5, 6, 7) and i.get('nsw', 0) == 0,
        "ISIF > 2 (cell shape/volume relaxation) requires NSW > 0",
    ),
    # Wavefunction/charge handling
    (
        lambda i: (
            i.get('lwave', True) is True
            and i.get('istart', 0) == 0
            and i.get('icharg', 0) == 0
            and i.get('nsw', 0) == 0
        ),
        "LWAVE=True with ISTART=0 and single SCF: WAVECAR generated but won't be used",
    ),
    (
        lambda i: i.get('istart', 0) > 0 and i.get('lwave', True) is False,
        "ISTART > 0 expects WAVECAR but LWAVE=False (no WAVECAR saved)",
    ),
    (
        lambda i: i.get('icharg', 0) > 0 and i.get('lcharg', True) is False,
        "ICHARG > 0 expects CHGCAR but LCHARG=False (no CHGCAR saved)",
    ),
    # K-points related
    (
        lambda i: i.get('ismear', 0) == -5 and i.get('kspacing', 0) > 0.05,
        "ISMEAR=-5 (tetrahedron) requires dense k-mesh (kspacing < 0.05)",
    ),
    # Accuracy settings
    (
        lambda i: (
            i.get('lreal', False) == 'Auto'
            and str(i.get('prec', 'Normal')).lower() == 'accurate'
        ),
        "LREAL=Auto with PREC=Accurate: consider LREAL=False for small cells",
    ),
    (
        lambda i: i.get('ediff', 1e-4) > 1e-5 and i.get('ibrion', -1) == 2,
        "EDIFF > 1E-5 with IBRION=2: consider tightening EDIFF for accurate forces",
    ),
    # Spin settings
    (
        lambda i: i.get('ispin', 1) == 1 and i.get('magmom') is not None,
        "MAGMOM specified but ISPIN=1 (non-spin-polarized): MAGMOM will be ignored",
    ),
    (
        lambda i: i.get('ispin', 1) == 2 and i.get('magmom') is None,
        "ISPIN=2 without MAGMOM: consider specifying initial magnetic moments",
    ),
    # DFT+U settings
    (
        lambda i: i.get('ldau', False) is True and i.get('ldautype') is None,
        "LDAU=True but LDAUTYPE not specified",
    ),
    # Electronic convergence
    (
        lambda i: i.get('algo', 'Normal') == 'Fast' and i.get('ismear', 0) == -5,
        "ALGO=Fast with ISMEAR=-5: may have convergence issues",
    ),
    # Hybrid functional settings
    (
        lambda i: (
            i.get('lhfcalc', False) is True
            and i.get('algo', 'Normal') not in ('Damped', 'All', 'Normal')
        ),
        "LHFCALC=True (hybrid): use ALGO=Damped or All for better convergence",
    ),
    # LOPTICS requires NBANDS
    (
        lambda i: i.get('loptics', False) is True and i.get('nbands') is None,
        "LOPTICS=True: consider explicitly setting NBANDS for optical calculations",
    ),
]


def validate_incar(
    incar: Dict[str, Any],
    structure_natoms: Optional[int] = None,
) -> ValidationResult:
    """Validate INCAR parameters for common mistakes.

    Args:
        incar: INCAR parameters dict (lowercase keys as aiida-vasp uses)
        structure_natoms: Number of atoms (enables additional checks)

    Returns:
        ValidationResult with errors and warnings

    Example:
        >>> result = validate_incar({'encut': 520, 'ibrion': 2, 'nsw': 0})
        >>> print(result.warnings[0])
        'IBRION >= 0 requires NSW > 0 for ionic relaxation'
    """
    result = ValidationResult()

    # Normalize to lowercase
    incar_lower = {k.lower(): v for k, v in incar.items()}

    # Check rules
    for check_fn, message in INCAR_RULES:
        try:
            if check_fn(incar_lower):
                result.warnings.append(message)
        except Exception:
            pass  # Skip if rule can't be evaluated

    # Physical plausibility checks
    encut = incar_lower.get('encut', 0)
    if encut > 0 and encut < 200:
        result.warnings.append(f"ENCUT={encut} is very low (typical: 400-600 eV)")
    if encut > 800:
        result.warnings.append(
            f"ENCUT={encut} is very high - verify this is needed"
        )

    sigma = incar_lower.get('sigma', 0.2)
    ismear = incar_lower.get('ismear', 0)
    if ismear == 0 and sigma > 0.1:
        result.warnings.append(
            f"ISMEAR=0 with SIGMA={sigma}: consider SIGMA=0.05 for insulators"
        )
    if ismear in (-4, -5) and sigma != 0:
        result.info['sigma_note'] = "SIGMA is ignored for tetrahedron methods"

    # NCORE/NPAR checks if natoms available
    if structure_natoms:
        ncore = incar_lower.get('ncore')
        npar = incar_lower.get('npar')
        if ncore and npar:
            result.warnings.append(
                "Both NCORE and NPAR specified: VASP will use NCORE, NPAR ignored"
            )

        # Store useful info
        result.info['natoms'] = structure_natoms

    return result


def validate_builder_inputs(builder_inputs: Dict[str, Any]) -> ValidationResult:
    """Validate complete builder_inputs dict.

    Checks:
    - Required keys present
    - INCAR consistency
    - Resources are reasonable
    - K-points spacing is appropriate

    Args:
        builder_inputs: Standard PS-TEROS builder_inputs dictionary

    Returns:
        ValidationResult with errors and warnings

    Example:
        >>> result = validate_builder_inputs({
        ...     'parameters': {
        ...         'incar': {'encut': 520, 'ismear': 0, 'sigma': 0.05}
        ...     },
        ...     'options': {
        ...         'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4}
        ...     },
        ...     'kpoints_spacing': 0.03,
        ... })
        >>> result.is_valid
        True
    """
    result = ValidationResult()

    # Check required keys
    if 'parameters' not in builder_inputs:
        result.errors.append("Missing 'parameters' in builder_inputs")
        return result

    params = builder_inputs.get('parameters', {})
    if 'incar' not in params:
        result.errors.append("Missing 'parameters.incar' in builder_inputs")
        return result

    # Validate INCAR
    incar_result = validate_incar(params['incar'])
    result.warnings.extend(incar_result.warnings)
    result.errors.extend(incar_result.errors)
    result.info.update(incar_result.info)

    # Check options
    options = builder_inputs.get('options', {})
    resources = options.get('resources', {})

    if not resources:
        result.warnings.append("No 'options.resources' specified")
    else:
        # Check for common resource mistakes
        num_machines = resources.get('num_machines', 0)
        num_cores = resources.get('num_cores_per_machine')
        num_mpiprocs = resources.get('num_mpiprocs_per_machine')

        if num_cores and num_mpiprocs:
            result.warnings.append(
                "Both num_cores_per_machine and num_mpiprocs_per_machine specified; "
                "prefer num_mpiprocs_per_machine for VASP"
            )

        if num_machines > 1 and not num_mpiprocs:
            result.warnings.append(
                f"num_machines={num_machines} but num_mpiprocs_per_machine not set"
            )

    # K-points
    kpoints_spacing = builder_inputs.get('kpoints_spacing')
    if kpoints_spacing:
        if kpoints_spacing > 0.1:
            result.warnings.append(
                f"kpoints_spacing={kpoints_spacing} is very coarse (typical: 0.03-0.05)"
            )
        elif kpoints_spacing < 0.01:
            result.warnings.append(
                f"kpoints_spacing={kpoints_spacing} is very fine (may be slow)"
            )
        result.info['kpoints_spacing'] = kpoints_spacing

    # Potential family
    potential_family = builder_inputs.get('potential_family')
    if potential_family:
        result.info['potential_family'] = potential_family

    return result


def estimate_kpoints_mesh(
    lattice_vectors: List[List[float]],
    spacing: float,
) -> Tuple[int, int, int]:
    """Estimate k-points mesh from lattice vectors and spacing.

    Pure Python calculation - works without AiiDA.

    Args:
        lattice_vectors: 3x3 list of lattice vectors in Angstrom
        spacing: K-points spacing in Angstrom^-1

    Returns:
        (kx, ky, kz) mesh as tuple

    Example:
        >>> # Simple cubic with a=4 A
        >>> lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        >>> mesh = estimate_kpoints_mesh(lattice, spacing=0.03)
        >>> print(mesh)  # Dense mesh for small cell
        (53, 53, 53)
    """
    import numpy as np

    lattice = np.array(lattice_vectors)
    reciprocal = 2 * np.pi * np.linalg.inv(lattice).T

    # Calculate mesh based on reciprocal lattice vector lengths
    mesh = []
    for i in range(3):
        rec_length = np.linalg.norm(reciprocal[i])
        k = max(1, int(np.ceil(rec_length / spacing)))
        mesh.append(k)

    return tuple(mesh)


def validate_structure_for_vasp(
    structure,
    check_min_distance: bool = True,
    min_distance_threshold: float = 0.5,
) -> ValidationResult:
    """Validate structure for VASP calculation.

    Works with pymatgen Structure or AiiDA StructureData.

    Args:
        structure: pymatgen Structure or AiiDA StructureData
        check_min_distance: Whether to check for short interatomic distances
        min_distance_threshold: Minimum allowed distance in Angstrom

    Returns:
        ValidationResult with warnings about potential issues
    """
    result = ValidationResult()

    # Get pymatgen structure if needed
    if hasattr(structure, 'get_pymatgen_structure'):
        pmg_structure = structure.get_pymatgen_structure()
    else:
        pmg_structure = structure

    # Basic info
    result.info['formula'] = pmg_structure.composition.reduced_formula
    result.info['natoms'] = len(pmg_structure)
    result.info['elements'] = list(set(str(s) for s in pmg_structure.species))

    # Check for very small/large cells
    volume = pmg_structure.volume
    natoms = len(pmg_structure)
    volume_per_atom = volume / natoms

    if volume_per_atom < 5:
        result.warnings.append(
            f"Very small volume per atom ({volume_per_atom:.1f} A^3/atom)"
        )
    if volume_per_atom > 500:
        result.warnings.append(
            f"Very large volume per atom ({volume_per_atom:.1f} A^3/atom) - "
            "may need vacuum correction"
        )

    # Check minimum distance
    if check_min_distance:
        try:
            min_dist = min(
                pmg_structure.get_neighbors(site, r=3.0)[0].distance
                for site in pmg_structure
                if pmg_structure.get_neighbors(site, r=3.0)
            )
            if min_dist < min_distance_threshold:
                result.errors.append(
                    f"Minimum interatomic distance is {min_dist:.2f} A "
                    f"(threshold: {min_distance_threshold} A)"
                )
            result.info['min_distance'] = min_dist
        except (ValueError, IndexError):
            # Skip if we can't calculate distance (single atom, etc.)
            pass

    # Check for fractional occupancies
    if any(site.species.num_atoms != 1 for site in pmg_structure):
        result.warnings.append(
            "Structure contains fractional occupancies - "
            "VASP requires integer occupancies"
        )

    return result
