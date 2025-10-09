#!/home/thiagotd/envs/psteros/bin/python
"""
Example script for relaxing slab structures from user-provided input slabs using PS-TEROS.

This script demonstrates how to:
1. Load the AiiDA profile
2. Load pre-generated slab structures from files
3. Set up parameters for bulk and reference relaxation
4. Set up slab relaxation parameters (VASP INCAR, scheduler options)
5. Create and run a WorkGraph that:
   - Relaxes the bulk structure
   - Relaxes reference structures (metal, oxygen) - NO nonmetal for binary oxides
   - Calculates formation enthalpy
   - Relaxes user-provided slabs in parallel with VASP (skips slab generation)
6. Access the relaxed slab structures and energies

Example: Ag2O (100) surface with user-provided slab structures
    - Relaxes Ag₂O bulk structure and references
    - Binary oxide: only metal (Ag) and oxygen (O2) references needed
    - Uses pre-generated slab structures from files
    - Relaxes each provided slab in parallel
    - Outputs relaxed slabs and their energies

Usage:
    source ~/envs/psteros/bin/activate && python slabs_input_relax_ag2o.py
"""

from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map


def main():
    """Main function to run the slab relaxation workflow with input slabs."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define structures directory
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"
    
    # Define directory containing pre-generated slab structures
    slabs_dir = "/home/thiagotd/git/PS-TEROS/examples/slabs/input_structures"

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
    # For binary oxides like Ag2O, we only need metal (Ag) and oxygen (O2)
    # NO nonmetal parameters needed!
    
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

    # ===== SLAB RELAXATION PARAMETERS =====
    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 2,  # Relax ions only, keep cell fixed
        "NSW": 100,
        "EDIFFG": -0.1,
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

    slab_potential_mapping = {"Ag": "Ag", "O": "O"}
    slab_kpoints_spacing = 0.3

    # ===== LOAD PRE-GENERATED SLABS FROM FILES =====
    print("\nLoading pre-generated slab structures from files...")
    
    # Create a dictionary of slab structures
    # Keys should follow the pattern "term_0", "term_1", etc.
    input_slabs = {}
    
    # Example: Load slab structures from CIF/POSCAR files
    # You can load as many slabs as you have files
    slab_files = [
        "slab_term_0_bin.cif",  # Replace with your actual filenames
        #"slab_term_1.cif",
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

    # ===== WORKFLOW DESCRIPTION =====
    print(f"\n{'='*80}")
    print("PS-TEROS WORKFLOW: User-Provided Slab Relaxation (Binary Oxide)")
    print(f"{'='*80}")
    print(f"\nThis workflow will:")
    print(f"  1. Relax bulk Ag₂O structure")
    print(f"  2. Relax reference structures (Ag metal, O₂ molecule)")
    print(f"     NOTE: No nonmetal reference needed for binary oxides!")
    print(f"  3. Calculate formation enthalpy")
    print(f"  4. Relax {len(input_slabs)} user-provided slab structures in parallel")
    print(f"  5. Extract energies for each relaxed slab")
    print(f"  6. Calculate surface energies as function of chemical potential")


    # ===== CREATE WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph...")
    print(f"{'='*80}")

    wg = build_core_workgraph_with_map(
        structures_dir=structures_dir,
        bulk_name="ag2o.cif",
        metal_name="Ag.cif",
        # NO nonmetal_name for binary oxides!
        oxygen_name="O2.cif",
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping={"Ag": "Ag", "O": "O"},
        metal_potential_mapping={"Ag": "Ag"},
        # NO nonmetal_potential_mapping for binary oxides!
        oxygen_potential_mapping={"O": "O"},
        kpoints_spacing=0.3,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        # NO nonmetal_parameters or nonmetal_options for binary oxides!
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=True,
        # Input slabs parameters (no slab generation needed!)
        input_slabs=input_slabs,
        # Slab relaxation
        relax_slabs=True,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        # Thermodynamics calculation
        compute_thermodynamics=True,
        thermodynamics_sampling=100,
        name="Ag2O_InputSlabs_Relax",
    )

    # Optional: Export to HTML
    try:
        html_file = "ag2o_input_slabs_relax.html"
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
    print(f"\n  Input slabs (passed through):")
    print(f"    - slab_structures.term_0, term_1, ... (StructureData)")
    print(f"\n  Relaxed slabs:")
    print(f"    - relaxed_slabs.term_0, term_1, ... (StructureData)")
    print(f"    - slab_energies.term_0, term_1, ... (Float)")
    print(f"\n  Surface energies:")
    print(f"    - surface_energies.term_0, term_1, ... (Dict with γ(Δμ))")

    print(f"\n{'-'*80}")
    print(f"ACCESSING RESULTS (after completion)")
    print(f"{'-'*80}")
    print(f"  In Python:")
    print(f"    from aiida import load_node")
    print(f"    wg = load_node({wg.pk})")
    print(f"    ")
    print(f"    # Get input slabs")
    print(f"    input_slabs = wg.outputs.slab_structures")
    print(f"    print(list(input_slabs.keys()))")
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
    print(f"    # Get surface energies")
    print(f"    surface_energies = wg.outputs.surface_energies")
    print(f"    for term_id in surface_energies.keys():")
    print(f"        data = surface_energies[term_id].get_dict()")
    print(f"        print(f'{{term_id}}: γ = {{data[\"gamma_array\"]}} J/m²')")
    print(f"    ")
    print(f"    # Export relaxed slab to file")
    print(f"    relaxed_term_0 = relaxed['term_0']")
    print(f"    atoms = relaxed_term_0.get_ase()")
    print(f"    atoms.write('relaxed_term_0.cif')")

    print(f"\n{'='*80}\n")

    return wg


if __name__ == "__main__":
    """
    Run the slab relaxation workflow with user-provided slabs.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Create your slab structure files in:
       /home/thiagotd/git/PS-TEROS/examples/slabs/input_structures/
       
       Example formats: CIF, POSCAR, xyz, etc.
       Naming: slab_term_0.cif, slab_term_1.cif, ...

    5. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    6. Run this script:
       source ~/envs/psteros/bin/activate && python slabs_input_relax.py
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
