#!/home/thiagotd/envs/aiida/bin/python
"""
Complete example testing ALL features of PS-TEROS WorkGraph for Ag2O.

This script demonstrates the full capabilities of the build_core_workgraph function,
including:
1. Bulk structure relaxation
2. Reference structures relaxation (metal, oxygen)
3. Formation enthalpy calculation
4. Slab generation from bulk structure
5. Slab relaxation (with unrelaxed SCF calculations)
6. Relaxation energy calculation (E_relaxed - E_unrelaxed)
7. Cleavage energy calculation for complementary slabs
8. Surface thermodynamics with chemical potential sampling

This tests all boolean flags and calculation modes:
- compute_relaxation_energy=True (default)
- compute_cleavage=True (default)
- compute_thermodynamics=True (default)

Material: Ag2O (binary oxide)
- Bulk: Ag2O (cuprite structure)
- Metal reference: Ag
- Oxygen reference: O2
- Slabs: (1,1,1) Miller index with multiple terminations

Usage:
    source ~/envs/aiida/bin/activate && python complete_ag2o_example.py
"""

import sys
import os
from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph


def main():
    """Main function to run the complete workflow testing all features."""

    # Load AiiDA profile
    print("=" * 80)
    print("PS-TEROS COMPLETE EXAMPLE - Testing ALL Features (Ag2O)")
    print("=" * 80)
    print("\nLoading AiiDA profile...")
    load_profile(profile='psteros')

    # ===== IMPORTANT: SET THIS TO YOUR PREVIOUS WORKGRAPH PK IN CASE YOU WANT TO RESTART =====
    # This is the PK of a previous PS-TEROS run that you want to restart from
    PREVIOUS_RUN_PK = 27054  # Replace with your actual PK

    # Define structure path (relative to this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    # Structure filenames
    bulk_filename = 'ag2o.cif'
    metal_filename = 'Ag.cif'
    oxygen_filename = 'O2.cif'

    # Create a dictionary of slab structures
    # Keys should follow the pattern "term_0", "term_1", etc.
    slabs_dir = "/home/thiagotd/git/worktree/PS-TEROS/feature-relax-energy/examples/complete/input_structures/ag2o"
    input_slabs = {}
    
    # Example: Load slab structures from CIF/POSCAR files
    # You can load as many slabs as you have files
    slab_files = [
        "slab_term_0_bin.cif",  # Replace with your actual filenames
        "slab_term_1_bin.cif",
    ]
    
    for idx, slab_file in enumerate(slab_files):
        try:
            slab_path = f"{slabs_dir}/{slab_file}"
            atoms = read(slab_path)
            structure = orm.StructureData(ase=atoms)
            # Store the structure so it can be used in the workflow
            structure.store()
            input_slabs[f"term_{idx}"] = structure
            print(f"  ✓ Loaded {slab_file} as term_{idx}")
        except FileNotFoundError:
            print(f"  ✗ Warning: {slab_file} not found, skipping...")
    
    if not input_slabs:
        print("\n✗ Error: No slab structures loaded!")
        print(f"  Please create slab structure files in: {slabs_dir}")
        print(f"  Expected files: {', '.join(slab_files)}")
        return None
    
    print(f"\n✓ Successfully loaded {len(input_slabs)} slab structures")
    print(f"✓ All structures stored in AiiDA database")

    # Define calculation parameters
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    print(f"\nStructures directory: {structures_dir}")
    print(f"  - Bulk: {bulk_filename}")
    print(f"  - Metal: {metal_filename}")
    print(f"  - Oxygen: {oxygen_filename}")

    # ===== BULK PARAMETERS (Ag2O) =====
    print("\n" + "=" * 80)
    print("BULK RELAXATION PARAMETERS (Ag2O)")
    print("=" * 80)

    bulk_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,  # Gaussian smearing for semiconductors
        'SIGMA': 0.05,
        'IBRION': 2,  # Conjugate gradient
        'ISIF': 3,    # Relax cell shape, volume, and atoms
        'NSW': 100,   # Max ionic steps
        'EDIFFG': -0.01,  # Force convergence (eV/Angstrom)
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    bulk_options = {
        'resources': {
            'num_machines': 1,
            "num_cores_per_machine": 40,
        },
        'queue_name': 'par40',
    }

    bulk_potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
    }

    print(f"  INCAR parameters: {len(bulk_parameters)} tags")
    print(f"  Resources: {bulk_options['resources']['num_cores_per_machine']} cores")

    # ===== METAL REFERENCE PARAMETERS (Ag) =====
    print("\n" + "=" * 80)
    print("METAL REFERENCE PARAMETERS (Ag)")
    print("=" * 80)

    metal_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 1,  # Methfessel-Paxton for metals
        'SIGMA': 0.2,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    metal_options = {
        'resources': {
            'num_machines': 1,
            "num_cores_per_machine": 40,
        },
        'queue_name': 'par40',
    }

    metal_potential_mapping = {
        'Ag': 'Ag',
    }

    print(f"  ISMEAR: {metal_parameters['ISMEAR']} (Methfessel-Paxton for metals)")
    print(f"  SIGMA: {metal_parameters['SIGMA']}")

    # ===== OXYGEN REFERENCE PARAMETERS (O2) =====
    print("\n" + "=" * 80)
    print("OXYGEN REFERENCE PARAMETERS (O2)")
    print("=" * 80)

    oxygen_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,  # Relax atoms only, keep cell fixed for O2 molecule
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': False,  # No LREAL for small systems
        'LWAVE': False,
        'LCHARG': False,
    }

    oxygen_options = {
        'resources': {
            'num_machines': 1,
            "num_cores_per_machine": 40,
        },
        'queue_name': 'par40',
    }

    oxygen_potential_mapping = {
        'O': 'O',
    }

    print(f"  Molecule: O2")
    print(f"  ISIF: {oxygen_parameters['ISIF']} (relax atoms only)")
    print(f"  LREAL: {oxygen_parameters['LREAL']} (exact for small systems)")

    # ===== SLAB PARAMETERS =====
    print("\n" + "=" * 80)
    print("SLAB GENERATION AND RELAXATION PARAMETERS")
    print("=" * 80)

    miller_indices = [1, 1, 1]  # (111) surface
    min_slab_thickness = 15.0   # Angstroms
    min_vacuum_thickness = 15.0  # Angstroms

    # Slab relaxation parameters (different from bulk - fix bottom layers)
    slab_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'NELM': 100,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.1,
        'IBRION': 2,
        'ISIF': 2,    # Relax atoms only, keep cell fixed for slabs
        'NSW': 100,
        'EDIFFG': -0.09,  # Slightly relaxed convergence for slabs
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': True,
        'LCHARG': True,
        'LASPH': True,
        #'IDIPOL': 3,  # Dipole correction in z-direction for slabs
        #'LDIPOL': True,
    }

    slab_options = {
        'resources': {
            'num_machines': 1,
            "num_cores_per_machine": 40,
        },
        'queue_name': 'par40',
    }

    print(f"  Miller indices: {miller_indices}")
    print(f"  Min slab thickness: {min_slab_thickness} Å")
    print(f"  Min vacuum thickness: {min_vacuum_thickness} Å")
    print(f"  ISIF: {slab_parameters['ISIF']} (atoms only)")
    #print(f"  IDIPOL: {slab_parameters['IDIPOL']} (dipole correction)")
    #print(f"  LDIPOL: {slab_parameters['LDIPOL']}")

    # ===== CALCULATION FLAGS =====
    print("\n" + "=" * 80)
    print("CALCULATION FLAGS (All calculations enabled)")
    print("=" * 80)

    relax_slabs = True
    compute_relaxation_energy = True
    compute_cleavage = True
    compute_thermodynamics = True
    thermodynamics_sampling = 50  # Grid points for chemical potential sampling

    print(f"  ✓ relax_slabs: {relax_slabs}")
    print(f"  ✓ compute_relaxation_energy: {compute_relaxation_energy}")
    print(f"      → Calculate E_relax - E_unrelaxed for each slab")
    print(f"  ✓ compute_cleavage: {compute_cleavage}")
    print(f"      → Calculate cleavage energies for complementary terminations")
    print(f"  ✓ compute_thermodynamics: {compute_thermodynamics}")
    print(f"      → Calculate surface energies with chemical potential sampling")
    print(f"      → Sampling grid: {thermodynamics_sampling} points")

    # ===== CREATE WORKGRAPH =====
    print("\n" + "=" * 80)
    print("CREATING WORKGRAPH")
    print("=" * 80)

    wg = build_core_workgraph(
        # Required parameters
        structures_dir=structures_dir,
        bulk_name=bulk_filename,

        # Code and potentials
        code_label=code_label,
        potential_family=potential_family,

        # Bulk parameters
        bulk_potential_mapping=bulk_potential_mapping,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,

        # Reference structures (for formation enthalpy)
        metal_name=metal_filename,
        metal_potential_mapping=metal_potential_mapping,
        metal_parameters=metal_parameters,
        metal_options=metal_options,

        # No nonmetal reference for binary oxide (Ag2O)
        # nonmetal_name=None,

        oxygen_name=oxygen_filename,
        oxygen_potential_mapping=oxygen_potential_mapping,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,

        # Slab generation
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,

        # Slab relaxation
        restart_from_node=PREVIOUS_RUN_PK,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        #input_slabs=input_slabs, # Predefined slab structures
        relax_slabs=relax_slabs,

        # Calculation flags (all enabled by default, but specified for clarity)
        compute_relaxation_energy=compute_relaxation_energy,
        compute_cleavage=compute_cleavage,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,

        # Other settings
        kpoints_spacing=0.4,
        clean_workdir=False,

        name='Ag2O_Complete_Workflow',
    )

    print("  ✓ WorkGraph created successfully")

    # ===== WORKGRAPH INFORMATION =====
    print("\n" + "=" * 80)
    print("WORKGRAPH STRUCTURE")
    print("=" * 80)

    print(f"\nWorkGraph name: {wg.name}")
    print(f"Total number of tasks: {len(wg.tasks)}")

    print("\nExpected workflow steps:")
    print("  1. Parallel relaxation of:")
    print("     - Bulk structure (Ag2O)")
    print("     - Metal reference (Ag)")
    print("     - Oxygen reference (O2)")
    print("  2. Formation enthalpy calculation")
    print("  3. Slab generation from relaxed bulk")
    print("  4. For each slab termination:")
    print("     a) SCF calculation (unrelaxed)")
    print("     b) Full relaxation")
    print("     c) Relaxation energy (E_relaxed - E_unrelaxed)")
    print("  5. Cleavage energy calculation")
    print("  6. Surface thermodynamics:")
    print("     - Oxide type identification")
    print("     - Chemical potential sampling")
    print("     - Surface energy calculation for each slab")

    # ===== EXPORT VISUALIZATION =====
    print("\n" + "=" * 80)
    print("EXPORTING VISUALIZATION")
    print("=" * 80)

    try:
        html_file = 'ag2o_complete_workgraph.html'
        wg.to_html(html_file)
        print(f"  ✓ WorkGraph visualization saved to: {html_file}")
        print(f"    Open this file in a browser to see the full workflow structure")
    except Exception as e:
        print(f"  ✗ Could not generate HTML visualization: {e}")

    # ===== SUBMIT WORKGRAPH =====
    print("\n" + "=" * 80)
    print("SUBMITTING WORKGRAPH")
    print("=" * 80)
    print("\nThis will submit the workflow to the AiiDA daemon.")
    print("The calculation will run on the cluster and may take several hours.")
    print("\nSubmitting...")

    wg.submit(wait=False)

    print(f"\n✓ WorkGraph submitted successfully!")
    print(f"  WorkGraph PK: {wg.pk}")

    # ===== MONITORING INSTRUCTIONS =====
    print("\n" + "=" * 80)
    print("MONITORING THE WORKFLOW")
    print("=" * 80)

    print(f"\nTo check status:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")

    print(f"\nTo monitor all running processes:")
    print(f"  verdi process list")
    print(f"  verdi process list -a -p1  # More detailed")

    print(f"\nTo check specific outputs after completion:")
    print(f"  verdi process show {wg.pk}")
    print(f"  # Outputs will include:")
    print(f"  #   - bulk_energy, metal_energy, oxygen_energy")
    print(f"  #   - formation_enthalpy")
    print(f"  #   - slab_structures (all generated terminations)")
    print(f"  #   - slab_energies (relaxed)")
    print(f"  #   - unrelaxed_slab_energies")
    print(f"  #   - relaxation_energies")
    print(f"  #   - cleavage_energies")
    print(f"  #   - surface_energies")

    # ===== EXPECTED OUTPUTS =====
    print("\n" + "=" * 80)
    print("EXPECTED OUTPUTS")
    print("=" * 80)

    print("\nUpon successful completion, the workflow will produce:")
    print("\n1. Bulk and Reference Energies:")
    print("   - bulk_energy: Total energy of relaxed Ag2O")
    print("   - metal_energy: Total energy of relaxed Ag")
    print("   - oxygen_energy: Total energy of relaxed O2")

    print("\n2. Formation Enthalpy:")
    print("   - formation_enthalpy: ΔH_f of Ag2O in eV/formula unit")

    print("\n3. Slab Structures:")
    print("   - slab_structures: All generated (111) terminations")

    print("\n4. Slab Energies:")
    print("   - unrelaxed_slab_energies: SCF energies")
    print("   - slab_energies: Relaxed energies")
    print("   - relaxation_energies: E_relax - E_unrelaxed for each slab")

    print("\n5. Cleavage Energies:")
    print("   - cleavage_energies: For complementary termination pairs")

    print("\n6. Surface Thermodynamics:")
    print("   - surface_energies: γ(μ_O) for each termination")
    print("   - Chemical potential range from metal-rich to oxygen-rich")

    print("\n" + "=" * 80)
    print("WORKFLOW SUBMITTED - Check status with verdi commands above")
    print("=" * 80 + "\n")

    return wg


if __name__ == '__main__':
    """
    Run the complete workflow testing all PS-TEROS features for Ag2O.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    5. Run this script:
       source ~/envs/aiida/bin/activate && python complete_ag2o_example.py

    This example will test ALL features:
    - Bulk and reference relaxations
    - Formation enthalpy calculation
    - Slab generation
    - Slab relaxation with unrelaxed SCF
    - Relaxation energy calculation
    - Cleavage energy calculation
    - Surface thermodynamics with chemical potential sampling
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
