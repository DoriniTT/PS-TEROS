#!/usr/bin/env python
"""
PdIn Intermetallic Surface Energy Calculation on Obelix Cluster

This example demonstrates:
1. Pre-flight feasibility analysis for stoichiometric+symmetric surfaces
2. Surface energy calculation (stoichiometric+symmetric surfaces only)
3. Wulff shape construction from calculated surface energies

Material: PdIn (B2 CsCl-type intermetallic)
Cluster: Obelix (Skylake nodes, hybrid MPI+OpenMP)
Code: VASP 6.5.1

The B2 structure has Pd at corners and In at body center, making it an
ideal test case for the stoichiometric finder since different Miller
indices can have varying stoichiometry behavior.

Note: The surface_energy module ONLY works with stoichiometric+symmetric
surfaces. For non-stoichiometric surfaces, use teros.core.thermodynamics.

Usage:
    source ~/envs/aiida/bin/activate
    python pdin_surface_energy_obelix.py

After completion:
    verdi process show <PK>
    python -c "from teros.core.surface_energy import get_wulff_shape_summary; ..."
"""

import sys
from aiida import load_profile, orm
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor

from teros.core.surface_energy import (
    build_metal_surface_energy_workgraph,
    analyze_miller_feasibility,
    get_feasibility_summary,
    NoStoichiometricSymmetricSurfaceError,
)


# =============================================================================
# CONFIGURATION
# =============================================================================

# AiiDA profile
AIIDA_PROFILE = 'presto'

# Obelix cluster configuration
CODE_LABEL = 'VASP-6.5.1-idefix@obelix'
POTENTIAL_FAMILY = 'PBE'

# Obelix uses hybrid MPI+OpenMP: 4 MPI processes, rest OpenMP threads
OBELIX_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,  # PROCESS_MPI=4
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N PdIn_surf''',
}

# VASP parameters optimized for metals
BULK_PARAMETERS = {
    'prec': 'Accurate',
    'encut': 520,           # High cutoff for accuracy
    'ediff': 1e-6,          # Tight electronic convergence
    'ismear': 1,            # Methfessel-Paxton for metals
    'sigma': 0.2,           # Smearing width
    'ibrion': 2,            # Conjugate gradient relaxation
    'isif': 3,              # Full cell + ion relaxation for bulk
    'nsw': 100,             # Max ionic steps
    'ediffg': -0.01,        # Force convergence (eV/A)
    'algo': 'Normal',
    'lreal': 'Auto',
    'lwave': False,
    'lcharg': False,
    'ncore': 22,            # Parallelization: 88/4 = 22 cores per MPI rank
}

# Slab parameters (same as bulk but fix cell shape)
SLAB_PARAMETERS = BULK_PARAMETERS.copy()
SLAB_PARAMETERS['isif'] = 2  # Fix cell, relax ions only

# Miller indices to study
MILLER_INDICES = [[1, 1, 0], [1, 0, 0], [1, 1, 1]]

# Slab geometry
MIN_SLAB_THICKNESS = 20.0   # Angstroms
MIN_VACUUM_THICKNESS = 20.0  # Angstroms


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def create_pdin_structure():
    """
    Create B2 (CsCl-type) PdIn structure.

    B2 structure:
    - Pd at (0, 0, 0) - corners
    - In at (0.5, 0.5, 0.5) - body center

    Returns:
        tuple: (pymatgen Structure, AiiDA StructureData)
    """
    # Experimental lattice constant for PdIn is ~3.25 A
    lattice = Lattice.cubic(3.25)
    pmg_structure = Structure(
        lattice,
        ["Pd", "In"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    )

    # Convert to AiiDA StructureData
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    aiida_structure = orm.StructureData(ase=ase_atoms)

    return pmg_structure, aiida_structure


def run_feasibility_analysis(pmg_structure):
    """
    Pre-flight analysis: check which Miller indices have valid surfaces.

    This analysis runs BEFORE any DFT calculations and helps identify:
    - Orientations suitable for this module (stoichiometric+symmetric)
    - Orientations that need teros.core.thermodynamics instead
    """
    print("\n" + "-" * 70)
    print("PRE-FLIGHT FEASIBILITY ANALYSIS")
    print("-" * 70)

    miller_tuples = [tuple(m) for m in MILLER_INDICES]

    reports = analyze_miller_feasibility(
        pmg_structure,
        miller_indices=miller_tuples,
        min_slab_thickness=MIN_SLAB_THICKNESS,
    )

    # Print detailed summary
    print(get_feasibility_summary(reports))

    # Summary table
    print("\nQuick Summary:")
    print("  Miller Index | Stoich | Sym | Both | Status")
    print("  " + "-" * 50)

    all_feasible = True
    for miller, report in reports.items():
        status = "OK" if report.has_valid_surfaces else "USE THERMODYNAMICS"
        if not report.has_valid_surfaces:
            all_feasible = False
        print(f"  {str(miller):14s} | {report.n_stoichiometric:6d} | {report.n_symmetric:3d} | "
              f"{report.n_both:4d} | {status}")

    print("  " + "-" * 50)

    if not all_feasible:
        print("\n  WARNING: Some orientations have no stoichiometric+symmetric surfaces!")
        print("  These orientations will raise NoStoichiometricSymmetricSurfaceError.")
        print("  Remove them from MILLER_INDICES or use teros.core.thermodynamics.")

    return reports, all_feasible


def build_and_submit_workflow(aiida_structure):
    """Build and submit the surface energy WorkGraph."""

    print("\n" + "-" * 70)
    print("BUILDING WORKGRAPH")
    print("-" * 70)

    # Build workflow
    # Note: require_stoichiometric_symmetric=True is the DEFAULT
    # The module only works with stoichiometric+symmetric surfaces
    wg = build_metal_surface_energy_workgraph(
        # Structure
        bulk_structure=aiida_structure,

        # Code configuration
        code_label=CODE_LABEL,
        potential_family=POTENTIAL_FAMILY,
        potential_mapping={'Pd': 'Pd', 'In': 'In'},
        kpoints_spacing=0.03,  # Dense k-mesh for metals
        clean_workdir=False,

        # Bulk calculation
        bulk_parameters=BULK_PARAMETERS,
        bulk_options=OBELIX_OPTIONS,

        # Slab generation
        miller_indices=MILLER_INDICES,
        min_slab_thickness=MIN_SLAB_THICKNESS,
        min_vacuum_thickness=MIN_VACUUM_THICKNESS,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,

        # Stoichiometric finder strategies (default: all strategies)
        stoichiometric_strategies=['filter_first', 'symmetrize_check', 'thickness_scan'],
        stoichiometric_max_thickness=30.0,

        # Slab relaxation
        slab_parameters=SLAB_PARAMETERS,
        slab_options=OBELIX_OPTIONS,
        slab_kpoints_spacing=0.04,  # Slightly coarser for slabs

        # Concurrency
        max_concurrent_jobs=4,

        # Workflow name
        name='PdIn_surface_energy',
    )

    print(f"\n  WorkGraph '{wg.name}' built successfully")
    print("  Only stoichiometric+symmetric slabs will be generated")
    print("  Non-stoichiometric terminations are filtered BEFORE DFT")

    # Submit
    print("\n  Submitting to AiiDA daemon...")
    wg.submit(wait=False)
    print(f"  Submitted with PK: {wg.pk}")

    return wg


def print_post_submission_info(wg):
    """Print information for monitoring and analyzing results."""

    print("\n" + "=" * 70)
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print("=" * 70)

    print(f"""
WorkGraph PK: {wg.pk}

Monitor Progress:
  verdi process status {wg.pk}
  verdi process show {wg.pk}
  verdi process report {wg.pk}

Expected Outputs:
  - bulk_energy: Total bulk energy (eV)
  - bulk_energy_per_atom: Bulk energy per atom (eV)
  - bulk_structure: Relaxed bulk structure
  - surface_energies: Dict with gamma for each (hkl, termination)
  - wulff_shape: Wulff shape analysis (shape_factor, dominant_facet, etc.)

After Completion - View Results:

  from aiida import orm
  from teros.core.surface_energy import get_wulff_shape_summary, visualize_wulff_shape

  wg = orm.load_node({wg.pk})

  # Surface energies
  se = wg.outputs.surface_energies.get_dict()
  for hkl, terms in se.items():
      for term, data in terms.items():
          print(f"{{hkl}} {{term}}: {{data['gamma_J_m2']:.3f}} J/m2")

  # Wulff shape summary
  print(get_wulff_shape_summary(wg.outputs.wulff_shape))

  # 3D visualization (requires matplotlib)
  visualize_wulff_shape(
      wg.outputs.bulk_structure,
      wg.outputs.surface_energies,
      save_path='pdin_wulff.png'
  )

Expected Surface Energies (approximate):
  PdIn is a B2 intermetallic with anisotropic surface energies.
  Typical values: 1.0 - 2.0 J/m2 depending on orientation.

Module Behavior:
  - Only stoichiometric+symmetric slabs are generated
  - Non-stoichiometric terminations filtered BEFORE DFT
  - Saves computational resources
  - If no valid surface exists: NoStoichiometricSymmetricSurfaceError raised
  - For non-stoichiometric surfaces: use teros.core.thermodynamics
""")

    print("=" * 70 + "\n")


def main():
    """Main entry point."""

    print("\n" + "=" * 70)
    print("PdIn INTERMETALLIC SURFACE ENERGY CALCULATION")
    print("Cluster: Obelix (Skylake, hybrid MPI+OpenMP)")
    print("=" * 70)

    # 1. Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile=AIIDA_PROFILE)
    print(f"   Profile '{AIIDA_PROFILE}' loaded")

    # 2. Create structure
    print("\n2. Creating PdIn B2 structure...")
    pmg_structure, aiida_structure = create_pdin_structure()
    print(f"   Formula: {pmg_structure.composition.reduced_formula}")
    print(f"   Space group: Pm-3m (B2/CsCl-type)")
    print(f"   Lattice constant: {pmg_structure.lattice.a:.3f} A")
    print(f"   Atoms: Pd @ (0,0,0), In @ (0.5,0.5,0.5)")

    # 3. Pre-flight feasibility analysis
    print("\n3. Running pre-flight feasibility analysis...")
    reports, all_feasible = run_feasibility_analysis(pmg_structure)

    # 4. Build and submit workflow
    print("\n4. Building and submitting WorkGraph...")
    try:
        wg = build_and_submit_workflow(aiida_structure)
    except NoStoichiometricSymmetricSurfaceError as e:
        print(f"\n  ERROR: {e}")
        print("\n  Some Miller indices have no valid stoichiometric+symmetric surfaces.")
        print("  Options:")
        print("    1. Remove problematic Miller indices from MILLER_INDICES")
        print("    2. Use teros.core.thermodynamics for those orientations")
        sys.exit(1)

    # 5. Print post-submission information
    print_post_submission_info(wg)

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
