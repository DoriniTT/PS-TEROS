#!/usr/bin/env python
"""Test script: VASP parallelization benchmark on bohr cluster.

This script tests the vasp_parallelization module by running short benchmark
calculations on the bohr cluster with VASP 6.4.3.

Usage:
    1. Activate AiiDA environment and ensure daemon is running:
       source ~/envs/aiida/bin/activate
       verdi daemon restart

    2. Run this script:
       python test_benchmark_bohr.py

    3. Monitor progress:
       verdi process show <PK>
       verdi process report <PK>

    4. After completion, view results:
       python -c "from teros.core.vasp_parallelization import print_benchmark_summary; print_benchmark_summary(<PK>)"

Note:
    This test uses a small Si structure (2 atoms) with minimal k-points
    and only 3 electronic steps for fast benchmarking.

    The workflow automatically configures VaspWorkChain to run only ONE
    VaspCalculation per benchmark (via max_iterations=1 and handler_overrides)
    to prevent restarts due to non-convergence with limited NELM.
"""

from aiida import load_profile, orm
from pymatgen.core import Lattice, Structure

load_profile()

from teros.core.vasp_parallelization import (
    build_parallelization_benchmark_workgraph,
    generate_benchmark_combinations,
)

# =============================================================================
# Configuration for bohr cluster
# =============================================================================

# VASP code on bohr (40 cores per node)
CODE_LABEL = "VASP-6.4.3@bohr"
POTENTIAL_FAMILY = "PBE"

# Number of MPI processes to test
# bohr has 40 cores per node, we'll test with a subset
NUM_PROCS = 40

# Scheduler options for PBS on bohr
OPTIONS = {
    "resources": {
        "num_machines": 1,
        "num_cores_per_machine": NUM_PROCS,
    },
    # PBS options for bohr
    "queue_name": "par40",
}

# Benchmark settings
NELM = 3  # Number of electronic steps (fast benchmark)
MAX_CONCURRENT_JOBS = 1  # Run sequentially for consistent timing

# Explicit NCORE values to test (subset for faster testing)
# For 8 procs: test 1, 2, 4, 8
NCORE_VALUES = [1, 2, 4]

# KPAR values to test
KPAR_VALUES = [1, 2]  # Keep simple for initial test

# Base INCAR parameters (lowercase as required by AiiDA-VASP)
BASE_INCAR = {
    "encut": 300,  # Lower for faster benchmark
    "ismear": 0,
    "sigma": 0.05,
    "ediff": 1e-4,  # Looser convergence for speed
    "prec": "Normal",  # Not Accurate, for speed
}

# K-points spacing (coarser = fewer k-points = faster)
KPOINTS_SPACING = 0.1


# =============================================================================
# Create test structure
# =============================================================================


def create_si_diamond() -> orm.StructureData:
    """Create a simple Si diamond structure (2 atoms)."""
    lattice = Lattice.cubic(5.43)
    si_structure = Structure(
        lattice,
        ["Si", "Si"],
        [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )
    structure = orm.StructureData(pymatgen=si_structure)
    structure.label = "Si-diamond"
    structure.description = "Silicon diamond structure for parallelization benchmark"
    return structure


# =============================================================================
# Main
# =============================================================================


def main():
    """Run the parallelization benchmark test."""
    print("=" * 70)
    print("VASP Parallelization Benchmark Test - bohr cluster")
    print("=" * 70)

    # Create structure
    print("\n1. Creating test structure...")
    structure = create_si_diamond()
    print(f"   Structure: {structure.get_formula()} ({len(structure.sites)} atoms)")

    # Preview combinations
    print("\n2. Benchmark combinations to test:")
    combinations = generate_benchmark_combinations(
        num_procs=NUM_PROCS,
        num_kpoints=1,  # Gamma-only expected with coarse spacing
        ncore_values=NCORE_VALUES,
        kpar_values=KPAR_VALUES,
    )
    for combo in combinations:
        print(f"   - {combo['label']}: NCORE={combo['ncore']}, KPAR={combo['kpar']}")
    print(f"   Total: {len(combinations)} benchmark runs")

    # Build WorkGraph
    print("\n3. Building benchmark WorkGraph...")
    wg = build_parallelization_benchmark_workgraph(
        structure=structure,
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        potential_mapping={"Si": "Si"},
        num_procs=NUM_PROCS,
        # Explicit NCORE/KPAR values
        ncore_values=NCORE_VALUES,
        kpar_values=KPAR_VALUES,
        # Benchmark configuration
        nelm=NELM,
        base_incar=BASE_INCAR,
        options=OPTIONS,
        kpoints_spacing=KPOINTS_SPACING,
        max_concurrent_jobs=MAX_CONCURRENT_JOBS,
        clean_workdir=True,  # Clean up after benchmark
        # Ranking weights
        weight_time=0.7,
        weight_memory=0.3,
        # WorkGraph name
        name="Si_Parallelization_Benchmark_Bohr",
    )

    # Submit
    print("\n4. Submitting workflow...")
    wg.submit()

    print(f"\n   Submitted WorkGraph: PK={wg.pk}")
    print(f"   Name: {wg.name}")

    # Instructions
    print("\n" + "=" * 70)
    print("Next steps:")
    print("=" * 70)

    print(f"\n1. Monitor progress:")
    print(f"   verdi process show {wg.pk}")
    print(f"   verdi process report {wg.pk}")

    print(f"\n2. Check daemon logs if stuck:")
    print(f"   verdi daemon logshow")

    print(f"\n3. After completion, view results:")
    print(
        f'   python -c "from teros.core.vasp_parallelization import print_benchmark_summary; print_benchmark_summary({wg.pk})"'
    )

    print(f"\n4. Or get results programmatically:")
    print(f"   from teros.core.vasp_parallelization import get_benchmark_results")
    print(f"   results = get_benchmark_results({wg.pk})")
    print(f"   print(results['recommended'])")

    print("\n" + "=" * 70)

    return wg.pk


if __name__ == "__main__":
    pk = main()
