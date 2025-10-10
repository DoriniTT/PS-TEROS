#!/usr/bin/env python
"""
Full working example: Relaxation energy calculation for Ag2O (100) surface.

This script demonstrates the new relaxation energy feature:
1. Relaxes bulk Ag2O and reference structures (Ag, O2)
2. Generates slab structures for (100) orientation
3. Performs SCF calculation on unrelaxed slabs (NSW=0, IBRION=-1) - AUTOMATIC
4. Performs relaxation calculation on slabs - USER SPECIFIED
5. Calculates relaxation energy: E_relax = E_relaxed - E_unrelaxed - AUTOMATIC
6. Displays results

The relaxation energy calculation is AUTOMATIC when relax_slabs=True!

Usage:
    source ~/envs/psteros/bin/activate && python ag2o_100_relaxation_energy.py
"""

import sys
import os

# Ensure we're using the feature branch code
feature_branch_path = '/home/thiagotd/git/worktree/PS-TEROS/feature-relax-energy'
if feature_branch_path not in sys.path:
    sys.path.insert(0, feature_branch_path)

from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph


def main():
    """Main function to run Ag2O (100) relaxation energy calculation."""

    # Load AiiDA profile
    print("\n" + "=" * 80)
    print("RELAXATION ENERGY CALCULATION: Ag2O (100) SURFACE")
    print("=" * 80)
    print("\nLoading AiiDA profile...")
    load_profile()

    # Define structures directory
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"

    # Define calculation parameters
    code_label = "VASP-VTST-6.4.3@bohr"
    potential_family = "PBE"

    # ===== BULK RELAXATION PARAMETERS =====
    bulk_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 3,  # Full relaxation
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    bulk_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== REFERENCE STRUCTURE PARAMETERS =====
    metal_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 1,
        "SIGMA": 0.2,
        "IBRION": 2,
        "ISIF": 3,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    metal_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # For binary oxide, nonmetal is not used but required by framework
    nonmetal_parameters = metal_parameters.copy()
    nonmetal_options = metal_options.copy()

    oxygen_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.01,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": False,
        "LWAVE": False,
        "LCHARG": False,
    }

    oxygen_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== SLAB GENERATION PARAMETERS =====
    miller_indices = [1, 0, 0]  # (100) surface
    min_slab_thickness = 10.0   # Angstroms (small for testing)
    min_vacuum_thickness = 15.0  # Angstroms

    # ===== SLAB RELAXATION PARAMETERS =====
    # These are used for the RELAXATION calculation
    # The SCF calculation (NSW=0, IBRION=-1) is set automatically
    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,      # Conjugate gradient
        "ISIF": 2,        # Relax atoms only, keep cell fixed
        "NSW": 50,        # Reduced for testing (use 100-200 for production)
        "EDIFFG": -0.05,  # Force convergence
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== POTENTIAL MAPPINGS =====
    bulk_potential_mapping = {'Ag': 'Ag', 'O': 'O'}
    metal_potential_mapping = {'Ag': 'Ag'}
    oxygen_potential_mapping = {'O': 'O'}
    slab_potential_mapping = {'Ag': 'Ag', 'O': 'O'}

    # ===== PRINT WORKFLOW INFO =====
    print("\n" + "-" * 80)
    print("WORKFLOW CONFIGURATION")
    print("-" * 80)
    print(f"Structure directory: {structures_dir}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")
    print(f"\nBulk: ag2o.cif (Ag2O)")
    print(f"Metal reference: Ag.cif")
    print(f"Oxygen reference: O2.cif")
    print(f"\nSlab generation:")
    print(f"  Miller indices: {miller_indices}")
    print(f"  Min slab thickness: {min_slab_thickness} Å")
    print(f"  Min vacuum: {min_vacuum_thickness} Å")
    print(f"\nSlab relaxation: ENABLED")
    print(f"  ISIF=2 (relax atoms, fix cell)")
    print(f"  NSW={slab_parameters['NSW']}")
    print(f"  EDIFFG={slab_parameters['EDIFFG']}")

    print("\n" + "-" * 80)
    print("RELAXATION ENERGY CALCULATION")
    print("-" * 80)
    print("The workflow will automatically:")
    print("  1. Generate slab structures from relaxed bulk")
    print("  2. Run SCF on unrelaxed slabs (NSW=0, IBRION=-1)")
    print("  3. Run relaxation on slabs (your parameters)")
    print("  4. Calculate E_relax = E_relaxed - E_unrelaxed")
    print("\nAll calculations run in PARALLEL!")

    # ===== BUILD WORKGRAPH =====
    print("\n" + "=" * 80)
    print("Building WorkGraph...")
    print("=" * 80)

    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name="ag2o.cif",
        metal_name="Ag.cif",
        nonmetal_name="Ag.cif",  # Dummy for binary oxide
        oxygen_name="O2.cif",
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping,
        metal_potential_mapping=metal_potential_mapping,
        nonmetal_potential_mapping=metal_potential_mapping,  # Dummy
        oxygen_potential_mapping=oxygen_potential_mapping,
        kpoints_spacing=0.3,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=0.3,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        relax_slabs=True,  # Enable slab relaxation
        compute_relaxation_energy=True,  # ← Enable relaxation energy calculation!
        clean_workdir=False,
        name='Ag2O_100_RelaxationEnergy',
    )

    print(f"\n✓ WorkGraph '{wg.name}' created successfully!")
    print(f"  Tasks: {len(wg.tasks)}")
    print(f"  Outputs: {len(wg.outputs)}")

    # ===== SUBMIT WORKGRAPH =====
    print("\n" + "=" * 80)
    print("Submitting WorkGraph...")
    print("=" * 80)

    wg.submit(wait=False)
    
    print(f"\n✓ WorkGraph submitted successfully!")
    print(f"\n{'='*80}")
    print(f"WORKGRAPH PK: {wg.pk}")
    print(f"{'='*80}")

    # ===== MONITORING INSTRUCTIONS =====
    print("\n" + "-" * 80)
    print("MONITORING")
    print("-" * 80)
    print(f"\nCheck status:")
    print(f"  verdi process show {wg.pk}")
    print(f"\nWatch progress:")
    print(f"  watch -n 5 'verdi process show {wg.pk}'")
    print(f"\nCheck daemon:")
    print(f"  verdi daemon status")

    # ===== RESULTS INSTRUCTIONS =====
    print("\n" + "-" * 80)
    print("ACCESSING RESULTS (after completion)")
    print("-" * 80)
    print(f"\nFrom Python:")
    print(f"  from aiida import orm")
    print(f"  node = orm.load_node({wg.pk})")
    print(f"")
    print(f"  # Available outputs:")
    print(f"  # - node.outputs.slab_structures         (unrelaxed slabs)")
    print(f"  # - node.outputs.unrelaxed_slab_energies (SCF energies)")
    print(f"  # - node.outputs.relaxed_slabs           (relaxed structures)")
    print(f"  # - node.outputs.slab_energies           (relaxed energies)")
    print(f"  # - node.outputs.relaxation_energies     (E_relax = E_relaxed - E_unrelaxed)")
    print(f"")
    print(f"  # Print relaxation energies:")
    print(f"  print('\\nRelaxation Energies:')")
    print(f"  print('-' * 40)")
    print(f"  for label, energy in node.outputs.relaxation_energies.items():")
    print(f"      print(f'  {{label}}: {{energy.value:+.4f}} eV')")

    print("\n" + "=" * 80)
    print("Workflow submitted! Check status with commands above.")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
