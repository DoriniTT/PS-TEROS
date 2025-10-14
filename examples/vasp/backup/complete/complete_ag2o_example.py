#!/home/thiagotd/envs/aiida/bin/python
"""
Complete example testing ALL features of PS-TEROS WorkGraph for Ag2O.

This script demonstrates the full capabilities of the build_core_workgraph function,
including:
1. Bulk structure relaxation
2. Reference structures relaxation (metal, oxygen)
3. Formation enthalpy calculation
4. Electronic properties (DOS and band structure) for bulk
5. Slab generation from bulk structure
6. Slab relaxation (with unrelaxed SCF calculations)
7. Relaxation energy calculation (E_relaxed - E_unrelaxed)
8. Cleavage energy calculation for complementary slabs
9. Surface thermodynamics with chemical potential sampling
10. Electronic properties (DOS and band structure) for selected slabs (NEW!)

This tests all boolean flags and calculation modes:
- compute_relaxation_energy=True (default)
- compute_cleavage=True (default)
- compute_thermodynamics=True (default)
- compute_electronic_properties_bulk=True
- compute_electronic_properties_slabs=True (NEW!)

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
from teros.core.builders import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults,
)


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
    #slabs_dir = "/home/thiagotd/git/worktree/PS-TEROS/feature-relax-energy/examples/complete/input_structures/ag2o"
    #input_slabs = {}
    
    ## Example: Load slab structures from CIF/POSCAR files
    ## You can load as many slabs as you have files
    #slab_files = [
    #    "slab_term_0_bin.cif",  # Replace with your actual filenames
    #    "slab_term_1_bin.cif",
    #]
    
    #for idx, slab_file in enumerate(slab_files):
    #    try:
    #        slab_path = f"{slabs_dir}/{slab_file}"
    #        atoms = read(slab_path)
    #        structure = orm.StructureData(ase=atoms)
    #        # Store the structure so it can be used in the workflow
    #        structure.store()
    #        input_slabs[f"term_{idx}"] = structure
    #        print(f"  ✓ Loaded {slab_file} as term_{idx}")
    #    except FileNotFoundError:
    #        print(f"  ✗ Warning: {slab_file} not found, skipping...")
    
    #if not input_slabs:
    #    print("\n✗ Error: No slab structures loaded!")
    #    print(f"  Please create slab structure files in: {slabs_dir}")
    #    print(f"  Expected files: {', '.join(slab_files)}")
    #    return None
    
    #print(f"\n✓ Successfully loaded {len(input_slabs)} slab structures")
    #print(f"✓ All structures stored in AiiDA database")

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

    # ===== ELECTRONIC PROPERTIES PARAMETERS (NEW!) =====
    print("\n" + "=" * 80)
    print("ELECTRONIC PROPERTIES PARAMETERS (DOS & Bands)")
    print("=" * 80)

    # Get electronic properties defaults
    ep_defaults = get_electronic_properties_defaults(
        energy_cutoff=bulk_parameters['ENCUT'],  # Match bulk ENCUT
        electronic_convergence=1e-5,
        ncore=4,
        ispin=2,  # Spin-polarized for Ag2O
        lasph=True,
        lreal="Auto",
        kpoints_mesh_density=0.3,  # SCF k-mesh density
        band_kpoints_distance=0.2,  # Band path density
        dos_kpoints_distance=0.2,  # DOS k-mesh density
        line_density=0.2,  # Points along high-symmetry lines
        nedos=2000,  # DOS grid points
        sigma_bands=0.01,  # Smearing for bands (eV)
        symprec=1e-4,  # Symmetry precision
        band_mode="seekpath-aiida",  # Use seekpath for band paths
    )

    compute_electronic_properties_bulk = True  # Enable DOS and bands

    print(f"  Band mode: {ep_defaults['band_settings']['band_mode']}")
    print(f"  Band k-points distance: {ep_defaults['band_settings']['band_kpoints_distance']}")
    print(f"  DOS k-points distance: {ep_defaults['band_settings']['dos_kpoints_distance']}")
    print(f"  Line density: {ep_defaults['band_settings']['line_density']}")
    print(f"  NEDOS (DOS grid points): {ep_defaults['dos']['NEDOS']}")
    print(f"  SCF k-mesh density: {ep_defaults['scf_kpoints_distance']}")
    print(f"\n  ✓ compute_electronic_properties_bulk: {compute_electronic_properties_bulk}")
    print(f"      → Calculate DOS and band structure for relaxed bulk")

    # ===== SLAB ELECTRONIC PROPERTIES PARAMETERS (NEW!) =====
    print("\n" + "=" * 80)
    print("SLAB ELECTRONIC PROPERTIES PARAMETERS (DOS & Bands)")
    print("=" * 80)

    # Get slab electronic properties defaults (denser k-point sampling for 2D systems)
    slab_ep_defaults = get_slab_electronic_properties_defaults(
        energy_cutoff=slab_parameters['ENCUT'],  # Match slab ENCUT
        electronic_convergence=1e-5,
        ncore=4,
        ispin=2,  # Spin-polarized for Ag2O
        lasph=True,
        lreal="Auto",
        kpoints_mesh_density=0.25,  # Denser than bulk for 2D
        band_kpoints_distance=0.15,  # Denser path sampling
        dos_kpoints_distance=0.2,
        line_density=0.15,  # More points along paths
        nedos=2000,
        sigma_bands=0.01,
        symprec=1e-4,
        band_mode="seekpath-aiida",
    )

    # Enable slab electronic properties calculation
    compute_electronic_properties_slabs = True

    # Define which slabs to calculate electronic properties for
    # You can specify per-slab parameter overrides here
    # For this example, we'll calculate for term_0 and term_1 with different settings
    slab_electronic_properties = {
        'term_0': {
            'bands_parameters': slab_ep_defaults,
            'bands_options': {
                'resources': {
                    'num_machines': 1,
                    'num_cores_per_machine': 40,
                },
                'queue_name': 'par40',
            },
            'band_settings': slab_ep_defaults['band_settings'],
        },
        'term_1': {
            'bands_parameters': slab_ep_defaults,
            'bands_options': {
                'resources': {
                    'num_machines': 1,
                    'num_cores_per_machine': 40,
                },
                'queue_name': 'par40',
            },
            'band_settings': slab_ep_defaults['band_settings'],
        },
        # You can add more terminations here or use custom parameters per slab:
        # 'term_2': {
        #     'bands_parameters': custom_params,
        #     'bands_options': high_memory_options,
        #     'band_settings': custom_settings,
        # },
    }

    print(f"  Band mode: {slab_ep_defaults['band_settings']['band_mode']}")
    print(f"  Band k-points distance: {slab_ep_defaults['band_settings']['band_kpoints_distance']}")
    print(f"  DOS k-points distance: {slab_ep_defaults['band_settings']['dos_kpoints_distance']}")
    print(f"  Line density: {slab_ep_defaults['band_settings']['line_density']}")
    print(f"  NEDOS (DOS grid points): {slab_ep_defaults['dos']['NEDOS']}")
    print(f"  SCF k-mesh density: {slab_ep_defaults['scf_kpoints_distance']}")
    print(f"\n  ✓ compute_electronic_properties_slabs: {compute_electronic_properties_slabs}")
    print(f"      → Calculate DOS and band structure for selected slabs")
    print(f"      → Selected terminations: {list(slab_electronic_properties.keys())}")
    print(f"      → Note: Denser k-point sampling tuned for 2D slab systems")

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
        #restart_from_node=PREVIOUS_RUN_PK,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        #input_slabs=input_slabs, # Predefined slab structures
        relax_slabs=relax_slabs,

        # Calculation flags (all enabled by default, but specified for clarity)
        compute_relaxation_energy=compute_relaxation_energy,
        compute_cleavage=compute_cleavage,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,

        # Electronic properties (NEW!)
        compute_electronic_properties_bulk=compute_electronic_properties_bulk,
        bands_parameters=ep_defaults,
        band_settings=ep_defaults['band_settings'],
        bands_options=bulk_options,  # Use same resources as bulk

        # Slab electronic properties (NEW!)
        compute_electronic_properties_slabs=compute_electronic_properties_slabs,
        slab_electronic_properties=slab_electronic_properties,
        slab_bands_parameters=slab_ep_defaults,
        slab_band_settings=slab_ep_defaults['band_settings'],
        slab_bands_options=slab_options,  # Use same resources as slab relaxation

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
    print("  3. Electronic properties for bulk:")
    print("     a) SCF calculation (LWAVE=True, LCHARG=True)")
    print("     b) Band structure along high-symmetry paths")
    print("     c) Density of states (DOS) with tetrahedron method")
    print("  4. Slab generation from relaxed bulk")
    print("  5. For each slab termination:")
    print("     a) SCF calculation (unrelaxed)")
    print("     b) Full relaxation")
    print("     c) Relaxation energy (E_relaxed - E_unrelaxed)")
    print("  6. Cleavage energy calculation")
    print("  7. Surface thermodynamics:")
    print("     - Oxide type identification")
    print("     - Chemical potential sampling")
    print("     - Surface energy calculation for each slab")
    print("  8. Electronic properties for selected slabs (NEW!):")
    print("     a) SCF calculation (LWAVE=True, LCHARG=True)")
    print("     b) Band structure along high-symmetry paths")
    print("     c) Density of states (DOS)")
    print(f"     d) For terminations: {list(slab_electronic_properties.keys())}")

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
    print(f"  #   - bulk_bands, bulk_dos, bulk_primitive_structure, bulk_seekpath_parameters")
    print(f"  #   - slab_structures (all generated terminations)")
    print(f"  #   - slab_energies (relaxed)")
    print(f"  #   - unrelaxed_slab_energies")
    print(f"  #   - relaxation_energies")
    print(f"  #   - cleavage_energies")
    print(f"  #   - surface_energies")
    print(f"  #   - slab_bands, slab_dos, slab_primitive_structures, slab_seekpath_parameters (NEW!)")

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

    print("\n3. Electronic Properties:")
    print("   - bulk_bands: Band structure along high-symmetry paths")
    print("   - bulk_dos: Density of states")
    print("   - bulk_primitive_structure: Primitive cell used for band calculation")
    print("   - bulk_seekpath_parameters: Seekpath symmetry information")

    print("\n4. Slab Structures:")
    print("   - slab_structures: All generated (111) terminations")

    print("\n5. Slab Energies:")
    print("   - unrelaxed_slab_energies: SCF energies")
    print("   - slab_energies: Relaxed energies")
    print("   - relaxation_energies: E_relax - E_unrelaxed for each slab")

    print("\n6. Cleavage Energies:")
    print("   - cleavage_energies: For complementary termination pairs")

    print("\n7. Surface Thermodynamics:")
    print("   - surface_energies: γ(μ_O) for each termination")
    print("   - Chemical potential range from metal-rich to oxygen-rich")

    print("\n8. Slab Electronic Properties (NEW!):")
    print("   - slab_bands: Band structures for selected slabs")
    print("   - slab_dos: Density of states for selected slabs")
    print("   - slab_primitive_structures: Primitive cells for each slab")
    print("   - slab_seekpath_parameters: Seekpath info for each slab")
    print(f"   - Available for: {list(slab_electronic_properties.keys())}")

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
    - Electronic properties (DOS and bands) for bulk
    - Slab generation
    - Slab relaxation with unrelaxed SCF
    - Relaxation energy calculation
    - Cleavage energy calculation
    - Surface thermodynamics with chemical potential sampling
    - Electronic properties (DOS and bands) for selected slabs (NEW!)
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
