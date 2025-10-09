#!/home/thiagotd/envs/psteros/bin/python
"""
Example script for generating and relaxing slab structures from a relaxed BINARY oxide using PS-TEROS.

This script demonstrates the binary oxide workflow (Ag2O) with cleavage energy calculations:
1. Load the AiiDA profile
2. Set up parameters for bulk relaxation
3. Set up slab generation parameters (Miller indices, thickness, vacuum, etc.)
4. Set up slab relaxation parameters (VASP INCAR, scheduler options)
5. Create and run a WorkGraph that:
   - Relaxes the bulk structure
   - Generates all slab terminations
   - Relaxes all slabs in parallel with VASP
   - Computes surface energies γ(Δμ_O) for binary oxide
   - Computes cleavage energies for complementary slab pairs
6. Access the relaxed slab structures, energies, thermodynamics, and cleavage energies

Example: Ag2O (100) surface with full slab relaxation, thermodynamics, and cleavage energy
    - Relaxes Ag₂O bulk structure (BINARY OXIDE: Ag2O)
    - Generates all unique terminations for the (100) orientation
    - Relaxes each slab termination in parallel
    - Computes γ(Δμ_O) surface energies as function of oxygen chemical potential
    - Computes cleavage energies Ec for complementary termination pairs
    - Outputs relaxed slabs, energies, thermodynamics, and cleavage data

Usage:
    source ~/envs/psteros/bin/activate && python slabs_relax_ag2o_cleavage.py
"""

from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph_with_map


def main():
    """Main function to run the slab generation and relaxation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

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
        "ISIF": 3,
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
    # For Ag2O (binary oxide), we only need metal (Ag) and oxygen (O2) references
    # The nonmetal reference will be the same as metal (dummy for binary oxides)
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

    # For binary oxide, nonmetal is not physically meaningful but required by framework
    # We use same parameters as metal (won't affect binary thermodynamics calculation)
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
    min_slab_thickness = 10.0  # Angstroms
    min_vacuum_thickness = 15.0  # Angstroms

    # Slab generation options
    lll_reduce = True
    center_slab = True
    symmetrize = True
    primitive = True
    in_unit_planes = False
    max_normal_search = None

    # ===== THERMODYNAMICS PARAMETERS =====
    compute_thermodynamics = True  # Enable surface energy calculations
    thermodynamics_sampling = 100  # Grid resolution for chemical potential sampling

    # ===== CLEAVAGE ENERGY PARAMETERS =====
    compute_cleavage = True  # Enable cleavage energy calculations for complementary slab pairs

    # ===== SLAB RELAXATION PARAMETERS =====
    # Enable slab relaxation
    relax_slabs = True

    # VASP parameters for slab relaxation
    # For slabs, we typically want to:
    # - Fix bottom layers (not implemented here, but can be added via SELECTIVE DYNAMICS)
    # - Relax only atomic positions, not cell (ISIF=2)
    # - Use dipole corrections for asymmetric slabs (IDIPOL, LDIPOL)
    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,  # Conjugate gradient
        "ISIF": 2,  # Relax atoms only, keep cell fixed
        "NSW": 100,  # Max ionic steps
        "EDIFFG": -0.1,  # Tighter force convergence for surfaces
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
        # Optional: Add dipole corrections for asymmetric slabs
        # "IDIPOL": 3,      # Dipole correction along z-axis (c-axis)
        # "LDIPOL": True,   # Turn on dipole corrections
    }

    # Scheduler options for slab calculations
    # Slabs are typically larger than bulk, may need more resources
    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # Use the same potential mapping and kpoints as bulk
    # But you can customize if needed
    slab_potential_mapping = {"Ag": "Ag", "Ag": "Ag", "O": "O"}
    slab_kpoints_spacing = 0.3  # Can be different from bulk if needed

    # ===== PRINT WORKFLOW INFO =====
    print(f"\n{'='*80}")
    print(f"SLAB GENERATION + RELAXATION WORKFLOW FOR Ag2O (BINARY OXIDE)")
    print(f"{'='*80}")
    print(f"\nStructures directory: {structures_dir}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")

    print(f"\n{'-'*80}")
    print(f"BULK RELAXATION")
    print(f"{'-'*80}")
    print(f"  Structure: ag2o.cif")
    print(f"  ISIF=3 (full relaxation: cell + atoms)")

    print(f"\n{'-'*80}")
    print(f"SLAB GENERATION")
    print(f"{'-'*80}")
    print(f"  Miller indices: {miller_indices}")
    print(f"  Min slab thickness: {min_slab_thickness} Å")
    print(f"  Min vacuum thickness: {min_vacuum_thickness} Å")
    print(f"  Symmetrize: {symmetrize}")
    print(f"  LLL reduce: {lll_reduce}")

    print(f"\n{'-'*80}")
    print(f"SLAB RELAXATION")
    print(f"{'-'*80}")
    print(f"  Enabled: {relax_slabs}")
    print(f"  ISIF=2 (relax atoms only, fixed cell)")
    print(f"  EDIFFG={slab_parameters['EDIFFG']} (tighter convergence)")
    print(f"  All slabs will be relaxed in parallel")

    print(f"\n{'-'*80}")
    print(f"WORKFLOW STEPS")
    print(f"{'-'*80}")
    print(f"  1. Relax bulk Ag2O structure")
    print(f"  2. Relax reference structures (Ag, P, O2) in parallel")
    print(f"  3. Calculate formation enthalpy")
    print(
        f"  4. Generate all slab terminations for ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]}) orientation"
    )
    print(f"  5. Relax ALL slab terminations in parallel with VASP")
    print(f"  6. Extract energies for each relaxed slab")
    print(f"  7. Compute surface energies γ(Δμ_O) for each termination")
    print(f"  8. Compute cleavage energies Ec for complementary termination pairs")

    # ===== CREATE WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph...")
    print(f"{'='*80}")

    wg = build_core_workgraph_with_map(
        structures_dir=structures_dir,
        bulk_name="ag2o.cif",
        metal_name="Ag.cif",
        nonmetal_name="Ag.cif",
        oxygen_name="O2.cif",
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping={"Ag": "Ag", "Ag": "Ag", "O": "O"},
        metal_potential_mapping={"Ag": "Ag"},
        nonmetal_potential_mapping={"Ag": "Ag"},
        oxygen_potential_mapping={"O": "O"},
        kpoints_spacing=0.3,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=True,
        # Slab generation
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        # Slab relaxation
        relax_slabs=relax_slabs,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        # Thermodynamics
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        # Cleavage energy
        compute_cleavage=compute_cleavage,
        name=f"Ag2O_SlabsRelax_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}",
    )

    # Optional: Export to HTML
    try:
        html_file = f"ag2o_slabs_relax_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}.html"
        wg.to_html(html_file)
        print(f"\n✓ WorkGraph visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML visualization: {e}")

    # ===== SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Submitting WorkGraph...")
    print(f"{'='*80}")
    wg.submit(wait=False)

    print(f"\n✓ WorkGraph submitted successfully!")
    print(f"\nWorkGraph PK: {wg.pk}")

    print(f"\n{'-'*80}")
    print(f"MONITORING")
    print(f"{'-'*80}")
    print(f"  Check status:     verdi process show {wg.pk}")
    print(f"  Check report:     verdi process report {wg.pk}")
    print(f"  List processes:   verdi process list -a -p1")

    print(f"\n{'-'*80}")
    print(f"EXPECTED OUTPUTS")
    print(f"{'-'*80}")
    print(f"  Formation enthalpy:")
    print(f"    - formation_enthalpy (Dict with ΔH_f)")
    print(f"\n  Unrelaxed slabs:")
    print(f"    - slab_structures.term_0, term_1, ... (StructureData)")
    print(f"\n  Relaxed slabs:")
    print(f"    - relaxed_slabs.term_0, term_1, ... (StructureData)")
    print(f"    - slab_energies.term_0, term_1, ... (Float)")
    print(f"\n  Surface energies (thermodynamics):")
    print(f"    - surface_energies.term_0, term_1, ... (Dict with γ(Δμ_O) - binary oxide)")
    print(f"\n  Cleavage energies:")
    print(f"    - cleavage_energies.pair_0_N, pair_1_N-1, ... (Dict with Ec for complementary pairs)")

    print(f"\n{'-'*80}")
    print(f"ACCESSING RESULTS (after completion)")
    print(f"{'-'*80}")
    print(f"  In Python:")
    print(f"    from aiida import load_node")
    print(f"    wg = load_node({wg.pk})")
    print(f"    ")
    print(f"    # Get unrelaxed slabs")
    print(f"    unrelaxed = wg.outputs.slab_structures")
    print(f"    print(list(unrelaxed.keys()))")
    print(f"    ")
    print(f"    # Get relaxed slabs")
    print(f"    relaxed = wg.outputs.relaxed_slabs")
    print(f"    print(list(relaxed.keys()))")
    print(f"    ")
    print(f"    # Get slab energies")
    print(f"    energies = wg.outputs.slab_energies")
    print(f"    for term_id in energies.keys():")
    print(f"        energy = energies[term_id].value")
    print(f"        print(f'{{term_id}}: {{energy}} eV')")
    print(f"    ")
    print(f"    # Get cleavage energies")
    print(f"    cleavages = wg.outputs.cleavage_energies")
    print(f"    for pair_id in cleavages.keys():")
    print(f"        data = cleavages[pair_id].get_dict()")
    print(f"        ec = data['cleavage_energy_eV_A2']")
    print(f"        print(f'{{pair_id}}: {{ec:.4f}} eV/Å²')")
    print(f"    ")
    print(f"    # Export relaxed slab to file")
    print(f"    relaxed_term_0 = relaxed['term_0']")
    print(f"    atoms = relaxed_term_0.get_ase()")
    print(f"    atoms.write('relaxed_term_0.cif')")

    print(f"\n{'='*80}\n")

    return wg


if __name__ == "__main__":
    """
    Run the slab generation and relaxation workflow.

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
       source ~/envs/aiida/bin/activate && python slabs_relax.py
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\n{'='*80}")
        print(f"ERROR")
        print(f"{'='*80}")
        print(f"{e}")
        print(f"\nFull traceback:")
        import traceback

        traceback.print_exc()
        import sys

        sys.exit(1)
