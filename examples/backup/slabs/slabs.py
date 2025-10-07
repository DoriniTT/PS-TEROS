#!/usr/bin/env python
"""
Example script for generating slab structures from a relaxed bulk oxide using PS-TEROS.

This script demonstrates how to:
1. Load the AiiDA profile
2. Set up parameters for bulk relaxation
3. Set up slab generation parameters (Miller indices, thickness, vacuum, etc.)
4. Create and run a WorkGraph that relaxes the bulk and generates slabs
5. Access the generated slab structures

Example: Ag3PO4 (100) surface
    - Relaxes Ag3PO4 bulk structure
    - Generates all unique terminations for the (100) orientation
    - Each termination is output as a separate StructureData node

Usage:
    source ~/envs/aiida/bin/activate && python slabs.py
"""

from aiida import load_profile, orm
from teros.workgraph import build_core_workgraph


def main():
    """Main function to run the slab generation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile()

    # Define structures directory
    structures_dir = "/home/thiagotd/git/PS-TEROS/teros/structures"

    # Define calculation parameters
    code_label = "VASP-VTST-6.4.3@bohr"
    potential_family = "PBE"

    # ===== BULK RELAXATION PARAMETERS =====
    # VASP parameters for bulk relaxation (Ag3PO4)
    bulk_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,  # Conjugate gradient
        "ISIF": 3,  # Relax cell shape, volume, and atoms
        "NSW": 100,  # Max ionic steps
        "EDIFFG": -0.1,  # Force convergence
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    # Scheduler options for bulk
    bulk_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== REFERENCE STRUCTURE PARAMETERS =====
    # Metal reference (Ag)
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

    # Nonmetal reference (P)
    nonmetal_parameters = {
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

    nonmetal_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # Oxygen reference (O2)
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

    # Optional slab parameters
    lll_reduce = True  # Apply LLL reduction to get more orthogonal cell
    center_slab = True  # Center slab in the c direction
    symmetrize = (
        True  # Generate all unique terminations (False) or only symmetric ones (True)
    )
    primitive = True  # Use primitive cell for slab generation
    in_unit_planes = False  # Restrict to unit planes
    max_normal_search = None  # Max search distance for surface normal

    print(f"\n{'='*70}")
    print(f"SLAB GENERATION WORKFLOW FOR Ag3PO4")
    print(f"{'='*70}")
    print(f"\nStructures directory: {structures_dir}")
    print(f"Code: {code_label}")
    print(f"Potential family: {potential_family}")

    print(f"\n{'-'*70}")
    print(f"BULK RELAXATION")
    print(f"{'-'*70}")
    print(f"  Structure: ag3po4.cif")
    print(f"  Will relax cell shape, volume, and atomic positions (ISIF=3)")

    print(f"\n{'-'*70}")
    print(f"REFERENCE CALCULATIONS (for formation enthalpy)")
    print(f"{'-'*70}")
    print(f"  Metal: Ag.cif")
    print(f"  Nonmetal: P.cif")
    print(f"  Oxygen: O2.cif")

    print(f"\n{'-'*70}")
    print(f"SLAB GENERATION PARAMETERS")
    print(f"{'-'*70}")
    print(f"  Miller indices: {miller_indices}")
    print(f"  Minimum slab thickness: {min_slab_thickness} Å")
    print(f"  Minimum vacuum thickness: {min_vacuum_thickness} Å")
    print(f"  LLL reduce: {lll_reduce}")
    print(f"  Center slab: {center_slab}")
    print(f"  Symmetrize: {symmetrize}")
    print(f"  Use primitive cell: {primitive}")
    print(f"  In unit planes: {in_unit_planes}")

    print(f"\n{'-'*70}")
    print(f"WORKFLOW DETAILS")
    print(f"{'-'*70}")
    print(f"  1. Relaxes bulk Ag3PO4 structure")
    print(f"  2. Relaxes reference structures (Ag, P, O2) in parallel")
    print(f"  3. Calculates formation enthalpy")
    print(
        f"  4. Generates all slab terminations for ({miller_indices[0]}{miller_indices[1]}{miller_indices[2]}) orientation"
    )
    print(f"  5. Each termination is exported as orthogonal, c-axis aligned structure")

    # Create the WorkGraph
    print(f"\n{'='*70}")
    print("Creating WorkGraph...")
    print(f"{'='*70}")

    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name="ag3po4.cif",
        metal_name="Ag.cif",
        nonmetal_name="P.cif",
        oxygen_name="O2.cif",
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
        metal_potential_mapping={"Ag": "Ag"},
        nonmetal_potential_mapping={"P": "P"},
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
        # Slab generation parameters
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        name=f"Ag3PO4_Slabs_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}",
    )

    # Optional: Export WorkGraph to HTML for visualization
    try:
        html_file = f"ag3po4_slabs_{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}.html"
        wg.to_html(html_file)
        print(f"\n✓ WorkGraph visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML visualization: {e}")

    # Submit the WorkGraph
    print(f"\n{'='*70}")
    print("Submitting WorkGraph...")
    print(f"{'='*70}")
    wg.submit(wait=False)

    print(f"\n✓ WorkGraph submitted successfully!")
    print(f"\nWorkGraph PK: {wg.pk}")

    print(f"\n{'-'*70}")
    print(f"MONITORING COMMANDS")
    print(f"{'-'*70}")
    print(f"  Check status:     verdi process show {wg.pk}")
    print(f"  Check report:     verdi process report {wg.pk}")
    print(f"  List processes:   verdi process list")
    print(f"  Watch progress:   verdi process list -a -p1")

    print(f"\n{'-'*70}")
    print(f"EXPECTED OUTPUTS")
    print(f"{'-'*70}")
    print(f"  Formation enthalpy:")
    print(f"    - bulk_energy, metal_energy, nonmetal_energy, oxygen_energy")
    print(
        f"    - bulk_structure, metal_structure, nonmetal_structure, oxygen_structure"
    )
    print(f"    - formation_enthalpy (Dict with ΔH_f in eV and kJ/mol)")
    print(f"\n  Slab structures:")
    print(f"    - slab_structures.term_0  (first termination)")
    print(f"    - slab_structures.term_1  (second termination)")
    print(f"    - slab_structures.term_N  (Nth termination)")
    print(f"    - Each is a StructureData node with orthogonal c-axis aligned slab")

    print(f"\n{'-'*70}")
    print(f"ACCESSING RESULTS (after completion)")
    print(f"{'-'*70}")
    print(f"  In Python:")
    print(f"    from aiida import load_node")
    print(f"    wg = load_node({wg.pk})")
    print(f"    ")
    print(f"    # Get slab structures")
    print(f"    slabs = wg.outputs.slab_structures")
    print(f"    print(list(slabs.keys()))  # ['term_0', 'term_1', ...]")
    print(f"    ")
    print(f"    # Access individual termination")
    print(f"    term_0 = slabs['term_0']")
    print(f"    atoms = term_0.get_ase()  # Get as ASE Atoms object")
    print(f"    atoms.write('term_0.cif')  # Export to file")
    print(f"    ")
    print(f"    # Get formation enthalpy")
    print(f"    hf = wg.outputs.formation_enthalpy.get_dict()")
    print(f"    print(f'ΔH_f = {{hf[\"formation_enthalpy_ev\"]}} eV/f.u.')")

    print(f"\n{'='*70}\n")

    return wg


if __name__ == "__main__":
    """
    Run the slab generation workflow.

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
       source ~/envs/aiida/bin/activate && python slabs.py
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\n{'='*70}")
        print(f"ERROR")
        print(f"{'='*70}")
        print(f"{e}")
        print(f"\nFull traceback:")
        import traceback

        traceback.print_exc()
        import sys

        sys.exit(1)
