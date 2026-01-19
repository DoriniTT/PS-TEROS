#!/usr/bin/env python
"""Example: Benchmark VASP parallelization parameters.

This script demonstrates how to use the vasp_parallelization module to find
optimal NCORE/KPAR settings for a given structure and computational resources.

Usage:
    1. Activate your AiiDA environment:
       source ~/envs/aiida/bin/activate
       verdi daemon restart

    2. Run the script:
       python run_benchmark.py

    3. Monitor progress:
       verdi process show <PK>
       verdi process report <PK>

    4. After completion, view results:
       python -c "from teros.core.vasp_parallelization import print_benchmark_summary; print_benchmark_summary(<PK>)"

Notes:
    - Adjust code_label, potential_family, and options for your cluster
    - For quick testing, reduce num_procs or limit ncore_values/kpar_values
    - The benchmark uses short SCF (NELM=3) for fast execution
"""

from aiida import orm
from pymatgen.core import Structure, Lattice

from teros.core.vasp_parallelization import (
    build_parallelization_benchmark_workgraph,
    generate_benchmark_combinations,
    print_benchmark_summary,
)

# =============================================================================
# Configuration
# =============================================================================

# VASP code configuration
# For obelix cluster (from home computer):
CODE_LABEL = "VASP-6.5.1-idefix@obelix"
POTENTIAL_FAMILY = "PBE"
NUM_PROCS = 4  # Number of MPI processes (PROCESS_MPI for hybrid MPI+OpenMP)

# Scheduler options for PBS on obelix
OPTIONS = {
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": NUM_PROCS,
    },
    "custom_scheduler_commands": """#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N VASPBenchmark""",
}

# Benchmark settings
NELM = 3  # Number of electronic steps (1-3 recommended for benchmarks)
MAX_CONCURRENT_JOBS = 2  # Run 2 benchmarks in parallel

# Base INCAR parameters (without parallelization - those are added by benchmark)
BASE_INCAR = {
    "ENCUT": 400,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "EDIFF": 1e-4,  # Looser convergence for benchmark speed
}

# K-points spacing
KPOINTS_SPACING = 0.04  # Coarser for faster benchmarks


# =============================================================================
# Create test structure
# =============================================================================


def create_test_structure() -> orm.StructureData:
    """Create a simple Si diamond structure for testing."""
    # Si diamond structure (2 atoms)
    lattice = Lattice.cubic(5.43)
    si_structure = Structure(
        lattice,
        ["Si", "Si"],
        [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )

    # Convert to AiiDA StructureData
    structure = orm.StructureData(pymatgen=si_structure)
    structure.label = "Si-diamond"
    structure.description = "Silicon diamond structure for parallelization benchmark"

    return structure


# =============================================================================
# Main
# =============================================================================


def main():
    """Run the parallelization benchmark."""
    print("=" * 60)
    print("VASP Parallelization Benchmark")
    print("=" * 60)

    # Create structure
    print("\n1. Creating test structure...")
    structure = create_test_structure()
    print(f"   Structure: {structure.get_formula()} ({len(structure.sites)} atoms)")

    # Preview combinations that will be tested
    print("\n2. Benchmark combinations to test:")
    combinations = generate_benchmark_combinations(
        num_procs=NUM_PROCS,
        num_kpoints=None,  # Will be auto-estimated
    )
    for combo in combinations:
        print(f"   - {combo['label']}: NCORE={combo['ncore']}, KPAR={combo['kpar']}")

    # Build WorkGraph
    print("\n3. Building benchmark WorkGraph...")
    wg = build_parallelization_benchmark_workgraph(
        structure=structure,
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        potential_mapping={"Si": "Si"},
        num_procs=NUM_PROCS,
        # Benchmark configuration
        nelm=NELM,
        base_incar=BASE_INCAR,
        options=OPTIONS,
        kpoints_spacing=KPOINTS_SPACING,
        max_concurrent_jobs=MAX_CONCURRENT_JOBS,
        # Ranking weights
        weight_time=0.7,
        weight_memory=0.3,
        # WorkGraph name
        name="Si_Parallelization_Benchmark",
    )

    # Submit
    print("\n4. Submitting workflow...")
    wg.submit()

    print(f"\n   Submitted WorkGraph: PK={wg.pk}")
    print(f"   Name: {wg.name}")
    print(f"   Number of benchmark tasks: {len(combinations)}")

    # Instructions
    print("\n" + "=" * 60)
    print("Next steps:")
    print("=" * 60)
    print(f"\n1. Monitor progress:")
    print(f"   verdi process show {wg.pk}")
    print(f"   verdi process report {wg.pk}")

    print(f"\n2. After completion, view results:")
    print(
        f'   python -c "from teros.core.vasp_parallelization import print_benchmark_summary; print_benchmark_summary({wg.pk})"'
    )

    print(f"\n3. Or in Python:")
    print(f"   from teros.core.vasp_parallelization import get_benchmark_results")
    print(f"   results = get_benchmark_results({wg.pk})")
    print(f"   print(results['recommended'])")

    return wg.pk


if __name__ == "__main__":
    pk = main()
