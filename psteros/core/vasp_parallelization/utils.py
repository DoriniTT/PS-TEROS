"""Utility functions for VASP parallelization benchmarking.

This module provides functions to:
- Generate sensible NCORE/KPAR values based on available resources
- Parse OUTCAR timing and memory information
- Estimate k-points count from structure and spacing
- Prepare INCAR parameters for benchmark calculations
"""

from __future__ import annotations

import math
import re
import typing as t

import numpy as np


def get_divisors(n: int) -> list[int]:
    """Get all divisors of n in ascending order.

    Args:
        n: Positive integer to find divisors for.

    Returns:
        List of all divisors of n, sorted ascending.

    Example:
        >>> get_divisors(12)
        [1, 2, 3, 4, 6, 12]
    """
    if n <= 0:
        return []
    divisors = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.append(i)
            if i != n // i:
                divisors.append(n // i)
    return sorted(divisors)


def generate_ncore_values(
    num_procs: int,
    max_ncore: int | None = None,
    include_one: bool = True,
) -> list[int]:
    """Generate sensible NCORE values for benchmarking.

    NCORE determines the number of cores that work on a single orbital.
    It should be a divisor of the number of MPI processes and typically
    a power of 2 or small multiple for optimal performance.

    Args:
        num_procs: Total number of MPI processes available.
        max_ncore: Maximum NCORE value to consider. Defaults to num_procs.
        include_one: Whether to include NCORE=1 in the list.

    Returns:
        List of NCORE values to test, sorted ascending.

    Example:
        >>> generate_ncore_values(16)
        [1, 2, 4, 8, 16]
        >>> generate_ncore_values(24)
        [1, 2, 4, 6, 8, 12, 24]
    """
    if num_procs <= 0:
        return [1]

    if max_ncore is None:
        max_ncore = num_procs

    # Get all divisors of num_procs
    divisors = get_divisors(num_procs)

    # Filter by max_ncore
    ncore_values = [d for d in divisors if d <= max_ncore]

    # Optionally remove 1 (often inefficient)
    if not include_one and 1 in ncore_values and len(ncore_values) > 1:
        ncore_values.remove(1)

    return ncore_values


def generate_kpar_values(
    num_kpoints: int,
    num_procs: int,
    max_kpar: int | None = None,
) -> list[int]:
    """Generate sensible KPAR values for benchmarking.

    KPAR determines the number of k-point groups for parallelization.
    It must divide both the number of k-points and the number of MPI processes.
    Effective when num_kpoints >= num_procs.

    Args:
        num_kpoints: Total number of k-points in the calculation.
        num_procs: Total number of MPI processes available.
        max_kpar: Maximum KPAR value to consider. Defaults to min(num_kpoints, num_procs).

    Returns:
        List of KPAR values to test, sorted ascending.

    Example:
        >>> generate_kpar_values(num_kpoints=16, num_procs=8)
        [1, 2, 4, 8]
        >>> generate_kpar_values(num_kpoints=10, num_procs=8)
        [1, 2]
    """
    if num_kpoints <= 0 or num_procs <= 0:
        return [1]

    if max_kpar is None:
        max_kpar = min(num_kpoints, num_procs)

    # Get divisors of both
    kpoint_divisors = set(get_divisors(num_kpoints))
    proc_divisors = set(get_divisors(num_procs))

    # KPAR must divide both
    common_divisors = kpoint_divisors.intersection(proc_divisors)

    # Filter by max_kpar
    kpar_values = sorted([d for d in common_divisors if d <= max_kpar])

    if not kpar_values:
        return [1]

    return kpar_values


def generate_benchmark_combinations(
    num_procs: int,
    num_kpoints: int | None = None,
    ncore_values: list[int] | None = None,
    kpar_values: list[int] | None = None,
    max_combinations: int = 20,
) -> list[dict[str, int]]:
    """Generate NCORE/KPAR combinations for benchmarking.

    Creates a list of parameter combinations to test. If explicit values
    are not provided, sensible defaults are auto-generated based on
    available resources.

    Args:
        num_procs: Total number of MPI processes available.
        num_kpoints: Number of k-points (used for KPAR generation if kpar_values not provided).
        ncore_values: Explicit list of NCORE values to test. Auto-generated if None.
        kpar_values: Explicit list of KPAR values to test. Auto-generated if None.
        max_combinations: Maximum number of combinations to return.

    Returns:
        List of dicts with 'ncore' and 'kpar' keys, labeled with 'label'.

    Example:
        >>> generate_benchmark_combinations(num_procs=8, num_kpoints=4)
        [
            {'label': 'ncore1_kpar1', 'ncore': 1, 'kpar': 1},
            {'label': 'ncore2_kpar1', 'ncore': 2, 'kpar': 1},
            {'label': 'ncore4_kpar1', 'ncore': 4, 'kpar': 1},
            ...
        ]
    """
    # Auto-generate NCORE values if not provided
    if ncore_values is None:
        ncore_values = generate_ncore_values(num_procs)

    # Auto-generate KPAR values if not provided
    if kpar_values is None:
        if num_kpoints is not None and num_kpoints > 1:
            kpar_values = generate_kpar_values(num_kpoints, num_procs)
        else:
            kpar_values = [1]

    # Generate all combinations
    combinations = []
    for ncore in ncore_values:
        for kpar in kpar_values:
            # Validate: NCORE * KPAR should not exceed num_procs
            # (though VASP handles this, we skip clearly inefficient combos)
            if ncore * kpar <= num_procs:
                label = f"ncore{ncore}_kpar{kpar}"
                combinations.append(
                    {
                        "label": label,
                        "ncore": ncore,
                        "kpar": kpar,
                    }
                )

    # Limit number of combinations
    if len(combinations) > max_combinations:
        # Prioritize: keep varied NCORE values, fewer KPAR variations
        # Sort by NCORE (primary), KPAR (secondary)
        combinations.sort(key=lambda x: (x["ncore"], x["kpar"]))
        # Sample evenly
        step = len(combinations) / max_combinations
        indices = [int(i * step) for i in range(max_combinations)]
        combinations = [combinations[i] for i in indices]

    return combinations


def estimate_kpoints_count(
    structure,
    kpoints_spacing: float = 0.03,
) -> int:
    """Estimate the number of k-points from structure and spacing.

    Uses the reciprocal lattice vectors to estimate k-point mesh dimensions,
    then calculates the total number of k-points in the IBZ (approximated).

    Args:
        structure: AiiDA StructureData or pymatgen Structure.
        kpoints_spacing: K-points spacing in Å⁻¹.

    Returns:
        Estimated number of k-points.

    Note:
        This is an approximation. Actual k-points depend on symmetry
        and the specific Monkhorst-Pack mesh generated by VASP.
    """
    # Get cell from structure
    if hasattr(structure, "get_pymatgen"):
        # AiiDA StructureData
        pmg_structure = structure.get_pymatgen()
        cell = pmg_structure.lattice.matrix
    elif hasattr(structure, "lattice"):
        # pymatgen Structure
        cell = structure.lattice.matrix
    elif hasattr(structure, "cell"):
        # ASE Atoms
        cell = structure.cell
    else:
        # Assume it's a numpy array
        cell = np.array(structure)

    # Calculate reciprocal lattice vectors
    # b = 2π * (a2 × a3) / (a1 · (a2 × a3)) etc.
    cell = np.array(cell)
    volume = np.abs(np.dot(cell[0], np.cross(cell[1], cell[2])))
    reciprocal = (
        2
        * np.pi
        * np.array(
            [
                np.cross(cell[1], cell[2]) / volume,
                np.cross(cell[2], cell[0]) / volume,
                np.cross(cell[0], cell[1]) / volume,
            ]
        )
    )

    # Calculate k-point mesh dimensions
    # Each dimension: ceil(|b_i| / kpoints_spacing)
    b_lengths = np.linalg.norm(reciprocal, axis=1)
    mesh = np.ceil(b_lengths / kpoints_spacing).astype(int)

    # Total k-points in full mesh
    total_kpoints = int(np.prod(mesh))

    # Approximate IBZ reduction (rough factor of 2-8 depending on symmetry)
    # For a conservative estimate, use factor of 2
    estimated_ibz = max(1, total_kpoints // 2)

    return estimated_ibz


def parse_outcar_timing(outcar_content: str) -> dict[str, t.Any]:
    """Parse timing and memory information from VASP OUTCAR.

    Extracts:
    - Total elapsed time
    - CPU time
    - Maximum memory usage
    - Per-loop timing

    Args:
        outcar_content: String content of OUTCAR file.

    Returns:
        Dict with timing/memory metrics:
        - 'elapsed_time': Total elapsed time in seconds
        - 'cpu_time': Total CPU time in seconds
        - 'memory_kb': Maximum memory used in KB
        - 'loop_times': List of (cpu, real) times per electronic step
        - 'parsing_success': Whether parsing succeeded
    """
    result = {
        "elapsed_time": None,
        "cpu_time": None,
        "memory_kb": None,
        "loop_times": [],
        "parsing_success": False,
    }

    # Pattern for total elapsed time
    # "Elapsed time (sec):    12.345"
    elapsed_match = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", outcar_content)
    if elapsed_match:
        result["elapsed_time"] = float(elapsed_match.group(1))
        result["parsing_success"] = True

    # Pattern for total CPU time (from General timing)
    # "Total CPU time used (sec):    24.567"
    cpu_match = re.search(r"Total CPU time used \(sec\):\s*([\d.]+)", outcar_content)
    if cpu_match:
        result["cpu_time"] = float(cpu_match.group(1))

    # Pattern for memory usage
    # "Maximum memory used (kb):     123456"
    # or "maximum memory used =     123456 kBytes"
    memory_match = re.search(r"[Mm]aximum memory used[^:]*:\s*(\d+)", outcar_content)
    if memory_match:
        result["memory_kb"] = int(memory_match.group(1))

    # Pattern for per-loop timing
    # "LOOP:  cpu time   2.34: real time    0.65"
    loop_pattern = re.compile(r"LOOP:\s+cpu time\s+([\d.]+):\s+real time\s+([\d.]+)")
    for match in loop_pattern.finditer(outcar_content):
        cpu = float(match.group(1))
        real = float(match.group(2))
        result["loop_times"].append({"cpu": cpu, "real": real})

    # Alternative: parse from timing section at end
    # If elapsed_time not found, try to get from LOOP times
    if result["elapsed_time"] is None and result["loop_times"]:
        # Sum of real times as approximation
        result["elapsed_time"] = sum(lt["real"] for lt in result["loop_times"])
        result["parsing_success"] = True

    return result


def prepare_benchmark_incar(
    base_incar: dict[str, t.Any] | None = None,
    ncore: int = 1,
    kpar: int = 1,
    nelm: int = 3,
    force_short_scf: bool = True,
) -> dict[str, t.Any]:
    """Prepare INCAR parameters for benchmark calculation.

    Sets parallelization parameters and optionally limits the calculation
    to a short SCF run for fast benchmarking.

    Note:
        All INCAR keys are converted to lowercase as required by AiiDA-VASP.

    Args:
        base_incar: Base INCAR parameters to use. Will be updated (not mutated).
        ncore: NCORE value for parallelization.
        kpar: KPAR value for parallelization.
        nelm: Maximum number of electronic steps (for short SCF).
        force_short_scf: If True, force NSW=0 and set NELM for quick benchmark.

    Returns:
        Updated INCAR dict with parallelization and benchmark settings.
        All keys are lowercase.

    Example:
        >>> prepare_benchmark_incar(
        ...     base_incar={'encut': 400, 'ismear': 0},
        ...     ncore=4, kpar=1, nelm=3
        ... )
        {'encut': 400, 'ismear': 0, 'ncore': 4, 'kpar': 1, 'nelm': 3, 'nsw': 0, ...}
    """
    # Convert all keys to lowercase (AiiDA-VASP requirement)
    if base_incar:
        incar = {k.lower(): v for k, v in base_incar.items()}
    else:
        incar = {}

    # Set parallelization parameters (lowercase)
    incar["ncore"] = ncore
    incar["kpar"] = kpar

    # Remove conflicting NPAR if present (NCORE takes precedence)
    incar.pop("npar", None)

    if force_short_scf:
        # Limit to short SCF calculation
        incar["nelm"] = nelm
        incar["nsw"] = 0  # No ionic steps
        incar["ibrion"] = -1  # No ionic relaxation

        # Ensure reasonable defaults for quick benchmark
        incar.setdefault("algo", "Fast")
        incar.setdefault("prec", "Normal")  # Not 'Accurate' for speed
        incar.setdefault("lwave", False)
        incar.setdefault("lcharg", False)
        incar.setdefault("lreal", "Auto")

    return incar


def rank_benchmark_results(
    results: list[dict[str, t.Any]],
    weight_time: float = 0.7,
    weight_memory: float = 0.3,
) -> list[dict[str, t.Any]]:
    """Rank benchmark results by combined time and memory score.

    Normalizes time and memory values and computes weighted score.
    Lower scores are better.

    Args:
        results: List of dicts with 'elapsed_time' and 'memory_kb' keys.
        weight_time: Weight for elapsed time in scoring (0-1).
        weight_memory: Weight for memory in scoring (0-1).

    Returns:
        Results sorted by score (best first) with 'score' and 'rank' added.
    """
    # Filter out failed results
    valid_results = [r for r in results if r.get("elapsed_time") is not None]

    if not valid_results:
        return results

    # Normalize values
    times = [r["elapsed_time"] for r in valid_results]
    min_time, max_time = min(times), max(times)
    time_range = max_time - min_time if max_time > min_time else 1.0

    memories = [r.get("memory_kb", 0) or 0 for r in valid_results]
    if any(m > 0 for m in memories):
        min_mem, max_mem = min(m for m in memories if m > 0), max(memories)
        mem_range = max_mem - min_mem if max_mem > min_mem else 1.0
    else:
        min_mem, mem_range = 0, 1.0

    # Calculate scores
    for r in valid_results:
        time_norm = (r["elapsed_time"] - min_time) / time_range
        mem = r.get("memory_kb", 0) or 0
        mem_norm = (mem - min_mem) / mem_range if mem > 0 else 0.5

        r["score"] = weight_time * time_norm + weight_memory * mem_norm

    # Sort by score
    valid_results.sort(key=lambda x: x["score"])

    # Add ranks
    for i, r in enumerate(valid_results):
        r["rank"] = i + 1

    return valid_results
