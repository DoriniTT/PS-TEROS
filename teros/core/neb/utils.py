"""Utility functions for NEB calculations.

Pure Python utilities for parameter preparation, structure validation,
and result formatting.
"""

import typing as t

from aiida import orm


# Default NEB INCAR parameters
DEFAULT_NEB_INCAR = {
    'ibrion': 3,          # BFGS optimizer for NEB
    'potim': 0.0,         # Use IOPT optimizer
    'spring': -5.0,       # Variable spring constant (negative = nudged)
    'iopt': 1,            # LBFGS optimizer (VTST)
    'nsw': 500,           # Max ionic steps
    'ediffg': -0.05,      # Force convergence (negative = force-based)
    'lreal': 'Auto',      # Real-space projection
    'lwave': False,       # Don't write WAVECAR
    'lcharg': False,      # Don't write CHGCAR
}


def get_neb_incar_parameters(
    n_images: int,
    climb: bool = True,
    spring_constant: float = -5.0,
    neb_optimizer: int = 1,
    force_convergence: float = -0.05,
    max_steps: int = 500,
) -> dict:
    """
    Generate NEB-specific INCAR parameters.

    Creates a dictionary of INCAR parameters optimized for NEB calculations
    using VASP with the VTST extensions.

    Args:
        n_images: Number of intermediate NEB images
        climb: Enable climbing image NEB (CI-NEB) for accurate saddle point
        spring_constant: NEB spring constant (negative = variable nudging)
        neb_optimizer: VTST optimizer (1=LBFGS, 3=FIRE, 7=fire)
        force_convergence: EDIFFG value (negative = force-based, eV/Å)
        max_steps: Maximum number of ionic steps (NSW)

    Returns:
        dict with NEB INCAR parameters suitable for deep_merge_dicts

    Example:
        >>> params = get_neb_incar_parameters(n_images=5, climb=True)
        >>> print(params)
        {'images': 5, 'lclimb': True, 'ibrion': 3, ...}

    Notes:
        - IBRION=3: Required for NEB (BFGS/quasi-Newton)
        - IOPT=1: LBFGS is generally most robust
        - IOPT=3: FIRE can be faster for well-behaved systems
        - Spring constant: Negative values give variable spring strength
        - For CI-NEB, the highest-energy image is allowed to climb
          to find the true saddle point
    """
    params = {
        'images': n_images,
        'lclimb': climb,
        'ibrion': 3,
        'potim': 0.0,
        'spring': spring_constant,
        'iopt': neb_optimizer,
        'nsw': max_steps,
        'ediffg': force_convergence,
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }
    return params


def validate_neb_structures(
    initial: orm.StructureData,
    final: orm.StructureData,
) -> None:
    """
    Validate that initial and final structures are compatible for NEB.

    Checks that:
    1. Both structures have the same number of atoms
    2. Both structures have the same chemical composition
    3. Both structures have compatible cell vectors

    Args:
        initial: Initial structure (reactant)
        final: Final structure (product)

    Raises:
        ValueError: If structures are incompatible for NEB
    """
    # Get pymatgen structures for comparison
    pmg_initial = initial.get_pymatgen_structure()
    pmg_final = final.get_pymatgen_structure()

    # Check atom count
    if len(pmg_initial) != len(pmg_final):
        raise ValueError(
            f"Initial and final structures have different atom counts: "
            f"{len(pmg_initial)} vs {len(pmg_final)}. "
            f"NEB requires the same number of atoms in both endpoints."
        )

    # Check composition
    comp_initial = pmg_initial.composition.reduced_formula
    comp_final = pmg_final.composition.reduced_formula
    if comp_initial != comp_final:
        raise ValueError(
            f"Initial and final structures have different compositions: "
            f"{comp_initial} vs {comp_final}. "
            f"NEB requires identical composition in both endpoints."
        )

    # Check cell vectors (allow small tolerance)
    import numpy as np
    lattice_initial = pmg_initial.lattice.matrix
    lattice_final = pmg_final.lattice.matrix

    # Allow 1% relative tolerance in cell parameters
    rel_diff = np.abs(lattice_initial - lattice_final) / (np.abs(lattice_initial) + 1e-10)
    max_diff = np.max(rel_diff)

    if max_diff > 0.01:  # 1% tolerance
        raise ValueError(
            f"Initial and final structures have significantly different cells. "
            f"Maximum relative difference: {max_diff * 100:.1f}%. "
            f"For NEB calculations, use fixed-cell calculations (ISIF=2) "
            f"or ensure both endpoints share the same cell."
        )


def compute_reaction_coordinate(
    structures: t.List[orm.StructureData],
) -> t.List[float]:
    """
    Calculate cumulative reaction coordinate along the NEB path.

    Computes the cumulative distance traveled along the minimum energy
    path by summing atomic displacements between consecutive images.

    Args:
        structures: List of StructureData along the NEB path
                   (including initial and final)

    Returns:
        List of cumulative distances (Å), normalized to [0, 1]

    Example:
        >>> coords = compute_reaction_coordinate(neb_images)
        >>> # coords = [0.0, 0.15, 0.35, 0.55, 0.78, 1.0]
    """
    import numpy as np

    if len(structures) < 2:
        return [0.0]

    cumulative_dist = [0.0]
    total_dist = 0.0

    for i in range(1, len(structures)):
        # Get positions from consecutive images
        pos_prev = structures[i - 1].get_ase().get_positions()
        pos_curr = structures[i].get_ase().get_positions()

        # Calculate displacement with minimum image convention
        # (for periodic boundary conditions)
        cell = structures[i].get_ase().get_cell()
        diff = pos_curr - pos_prev

        # Apply minimum image convention
        # This handles atoms that cross periodic boundaries
        frac_prev = np.linalg.solve(cell.T, pos_prev.T).T
        frac_curr = np.linalg.solve(cell.T, pos_curr.T).T
        frac_diff = frac_curr - frac_prev
        frac_diff = frac_diff - np.round(frac_diff)  # Apply MIC
        diff = np.dot(frac_diff, cell)

        # Compute RMS displacement for this step
        rms_displacement = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        total_dist += rms_displacement
        cumulative_dist.append(total_dist)

    # Normalize to [0, 1]
    if total_dist > 0:
        cumulative_dist = [d / total_dist for d in cumulative_dist]

    return cumulative_dist


def format_neb_results(results: dict) -> str:
    """
    Format NEB results into a human-readable string.

    Args:
        results: Dictionary from get_neb_results()

    Returns:
        Formatted string with NEB results summary
    """
    lines = []
    lines.append("=" * 60)
    lines.append("NEB CALCULATION RESULTS")
    lines.append("=" * 60)

    if 'summary' in results and results['summary']:
        summary = results['summary']

        # Barrier information
        lines.append("\nActivation Barriers:")
        lines.append("-" * 40)

        if 'forward_barrier_eV' in summary:
            forward_kj = summary.get('forward_barrier_kJ_mol', 0)
            reverse_kj = summary.get('reverse_barrier_kJ_mol', 0)
            reaction_kj = summary.get('reaction_energy_kJ_mol', 0)
            lines.append(
                f"  Forward barrier:  {summary['forward_barrier_eV']:.4f} eV "
                f"({forward_kj:.2f} kJ/mol)"
            )
            lines.append(
                f"  Reverse barrier:  {summary['reverse_barrier_eV']:.4f} eV "
                f"({reverse_kj:.2f} kJ/mol)"
            )
            lines.append(
                f"  Reaction energy:  {summary['reaction_energy_eV']:.4f} eV "
                f"({reaction_kj:.2f} kJ/mol)"
            )

            if summary.get('is_exothermic'):
                lines.append("  Reaction type:    Exothermic")
            else:
                lines.append("  Reaction type:    Endothermic")

        # Saddle point information
        if 'saddle_point_index' in summary:
            lines.append("\nSaddle Point:")
            lines.append("-" * 40)
            lines.append(f"  Image index:      {summary['saddle_point_index']}")
            lines.append(f"  Energy:           {summary.get('saddle_point_energy_eV', 'N/A'):.4f} eV")

        # Energy profile
        if 'energies_list_eV' in summary:
            lines.append("\nEnergy Profile:")
            lines.append("-" * 40)
            energies = summary['energies_list_eV']
            e_ref = energies[0]  # Reference to initial energy
            for i, e in enumerate(energies):
                rel_e = e - e_ref
                lines.append(f"  Image {i:02d}:  {e:.4f} eV  (ΔE = {rel_e:+.4f} eV)")

        # Calculation parameters
        lines.append("\nCalculation Parameters:")
        lines.append("-" * 40)
        lines.append(f"  Intermediate images: {summary.get('n_intermediate_images', 'N/A')}")
        lines.append(f"  Interpolation:       {summary.get('interpolation_method', 'N/A')}")
        lines.append(f"  Climbing image:      {summary.get('climbing_image', 'N/A')}")

    else:
        lines.append("\nNo summary data available.")

    lines.append("=" * 60)

    return '\n'.join(lines)


def estimate_neb_resources(
    n_images: int,
    n_atoms: int,
    base_time_hours: float = 24.0,
    base_cores: int = 40,
) -> dict:
    """
    Estimate computational resources for NEB calculation.

    Provides rough estimates for walltime and number of nodes based on
    system size and number of images. These are conservative estimates.

    Args:
        n_images: Number of intermediate NEB images
        n_atoms: Number of atoms in the unit cell
        base_time_hours: Base walltime estimate for simple systems (hours)
        base_cores: Base number of CPU cores

    Returns:
        dict with:
            - estimated_walltime_hours: Suggested walltime
            - estimated_nodes: Suggested number of nodes
            - cores_per_image: Recommended cores per image
            - notes: List of recommendations

    Note:
        These are rough estimates. Actual requirements depend heavily on:
        - System chemistry and convergence behavior
        - Distance between initial and final states
        - Quality of initial interpolation
        - Desired force convergence
    """
    notes = []

    # Scale walltime with system size and number of images
    # NEB scales roughly linearly with n_images
    size_factor = (n_atoms / 50) ** 1.5  # Nonlinear scaling with atoms
    image_factor = n_images / 5  # Linear scaling with images

    estimated_hours = base_time_hours * size_factor * image_factor

    # For NEB, it's often efficient to run with one node per image
    # or a fraction thereof
    if n_images <= 3:
        estimated_nodes = 1
        notes.append("Small NEB: single node should suffice")
    elif n_images <= 7:
        estimated_nodes = min(n_images, 3)
        notes.append(f"Medium NEB: {estimated_nodes} nodes recommended")
    else:
        estimated_nodes = min(n_images, 5)
        notes.append(f"Large NEB: {estimated_nodes} nodes recommended")

    cores_per_image = base_cores // max(1, min(n_images, estimated_nodes))

    # CI-NEB typically needs more iterations
    notes.append("CI-NEB may require 50-100% more walltime than regular NEB")

    # IDPP recommendation
    notes.append("IDPP interpolation recommended for reactive pathways")

    return {
        'estimated_walltime_hours': estimated_hours,
        'estimated_nodes': estimated_nodes,
        'cores_per_image': cores_per_image,
        'notes': notes,
    }


def get_endpoint_relax_incar() -> dict:
    """
    Get default INCAR parameters for endpoint relaxation.

    Returns INCAR settings appropriate for relaxing initial and final
    structures before NEB calculation.

    Returns:
        dict with INCAR parameters for ionic relaxation (ISIF=2)
    """
    return {
        'ibrion': 2,      # Conjugate gradient
        'isif': 2,        # Relax ions only, fixed cell
        'nsw': 100,       # Max ionic steps
        'ediffg': -0.02,  # Force convergence for endpoints
        'prec': 'Accurate',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    }


def get_stage2_parameters(
    stage1_params: dict,
    tighter_convergence: bool = True,
) -> dict:
    """
    Modify NEB parameters for Stage 2 (CI-NEB) calculation.

    Takes Stage 1 parameters and adjusts them for the climbing image
    refinement stage.

    Args:
        stage1_params: INCAR parameters from Stage 1
        tighter_convergence: If True, use tighter force convergence

    Returns:
        Modified parameters for Stage 2 CI-NEB
    """
    params = stage1_params.copy()

    # Enable climbing image
    params['lclimb'] = True

    # Optionally tighten convergence for more accurate TS
    if tighter_convergence:
        current_ediffg = params.get('ediffg', -0.05)
        if current_ediffg < 0:
            # Force-based: halve the convergence threshold
            params['ediffg'] = current_ediffg * 0.5
        else:
            # Energy-based: halve it
            params['ediffg'] = current_ediffg * 0.5

    return params
