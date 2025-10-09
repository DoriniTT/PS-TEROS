#!/home/thiagotd/envs/psteros/bin/python
"""
Example using default_ag3po4_builders for automatic slab generation and relaxation.

This script demonstrates how to use the default builders module WITHOUT providing
input slabs. The slabs will be automatically generated using the Miller indices,
slab thickness, and vacuum thickness specified in the defaults.

Usage:
    source ~/envs/psteros/bin/activate && python slabs_autogen_ag3po4.py
"""

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph_with_map
from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults

def main():
    """Main function to run the slab auto-generation and relaxation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define directories
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"

    # ===== GET DEFAULT BUILDERS FOR Ag3PO4 =====
    # This includes slab generation parameters!
    defaults = get_ag3po4_defaults(
        structures_dir=structures_dir,
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
    )

    # ===== OPTIONAL: OVERRIDE SLAB GENERATION PARAMETERS =====
    # Note: Only ONE Miller index can be specified at a time
    # Example: Generate (110) surface instead of (100)
    # defaults['miller_indices'] = [1, 1, 0]
    
    # Example: Use thicker slabs
    # defaults['min_slab_thickness'] = 15.0
    
    # Example: Use more vacuum
    # defaults['min_vacuum_thickness'] = 20.0

    # ===== OPTIONAL: OVERRIDE OTHER PARAMETERS =====
    # Example: Change ENCUT for bulk calculation
    # defaults['bulk_parameters']['ENCUT'] = 600
    
    # Example: Use more machines for slab calculations
    # defaults['slab_options']['resources']['num_machines'] = 2

    print(f"\n{'='*80}")
    print("DEFAULT SLAB GENERATION PARAMETERS")
    print(f"{'='*80}")
    print(f"Miller index: {defaults['miller_indices']}")
    print(f"Min slab thickness: {defaults['min_slab_thickness']} Å")
    print(f"Min vacuum thickness: {defaults['min_vacuum_thickness']} Å")
    print(f"\nNote: Currently only ONE surface can be generated at a time.")
    print(f"Slab will be automatically generated for this surface!")
    print(f"\nSlabs will be automatically generated for these surfaces!")

    # ===== CREATE AND SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph with automatic slab generation...")
    print(f"{'='*80}")

    wg = build_core_workgraph_with_map(
        **defaults,  # Unpack all default parameters (NO input_slabs provided!)
        name="Ag3PO4_AutoGen_Slabs",
    )

    # Optional: Export to HTML
    try:
        html_file = "ag3po4_autogen_slabs.html"
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
    print(f"\n  Generated slabs:")
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
    print(f"    # Get generated slabs")
    print(f"    generated_slabs = wg.outputs.slab_structures")
    print(f"    print(list(generated_slabs.keys()))")
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
    Run the automatic slab generation and relaxation workflow.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Ensure structure files exist:
       ls /home/thiagotd/git/PS-TEROS/examples/structures/
       Should contain: ag3po4.cif, Ag.cif, P.cif, O2.cif

    5. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    6. Run this script:
       source ~/envs/psteros/bin/activate && python slabs_autogen_ag3po4.py
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
