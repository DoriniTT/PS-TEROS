#!/home/thiagotd/envs/psteros/bin/python
"""
Example using default_ag2o_builders for automatic slab generation and relaxation.

This script demonstrates the BINARY OXIDE workflow for Ag2O:
1. Uses default builders to get all parameters
2. Automatically generates slab terminations for specified Miller index
3. Relaxes all slabs in parallel
4. Computes surface energies γ(Δμ_O) for binary oxide

Note: For binary oxides (Ag2O), the nonmetal reference is not physically meaningful
but is required by the framework. It uses the same parameters as the metal reference.

Usage:
    source ~/envs/psteros/bin/activate && python slabs_autogen_ag2o.py
"""

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph_with_map
from teros.core.builders.default_ag2o_builders import get_ag2o_defaults


def main():
    """Main function to run the Ag2O slab auto-generation and relaxation workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define directories
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"

    # ===== GET DEFAULT BUILDERS FOR Ag2O (BINARY OXIDE) =====
    # This includes ALL parameters: bulk, references, slab generation, etc.
    defaults = get_ag2o_defaults(
        structures_dir=structures_dir,
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
    )

    # ===== OPTIONAL: OVERRIDE SLAB GENERATION PARAMETERS =====
    # Note: Only ONE Miller index can be specified at a time
    
    # Example: Generate (110) surface instead of (100)
    # defaults['miller_indices'] = [1, 1, 0]
    
    # Example: Generate (111) surface
    # defaults['miller_indices'] = [1, 1, 1]
    
    # Example: Use thicker slabs
    # defaults['min_slab_thickness'] = 15.0
    
    # Example: Use more vacuum
    # defaults['min_vacuum_thickness'] = 20.0
    
    # Example: Disable symmetrization
    # defaults['symmetrize'] = False

    # ===== OPTIONAL: OVERRIDE OTHER PARAMETERS =====
    # Example: Change ENCUT for bulk calculation
    # defaults['bulk_parameters']['ENCUT'] = 600
    
    # Example: Use more machines for slab calculations
    # defaults['slab_options']['resources']['num_machines'] = 2
    
    # Example: Change thermodynamics sampling
    # defaults['thermodynamics_sampling'] = 200

    print(f"\n{'='*80}")
    print("DEFAULT SLAB GENERATION PARAMETERS")
    print(f"{'='*80}")
    print(f"Miller index: {defaults['miller_indices']}")
    print(f"Min slab thickness: {defaults['min_slab_thickness']} Å")
    print(f"Min vacuum thickness: {defaults['min_vacuum_thickness']} Å")
    print(f"Symmetrize: {defaults['symmetrize']}")
    print(f"Primitive: {defaults['primitive']}")
    print(f"\nNote: This is a BINARY OXIDE system (Ag2O)")
    print(f"      Surface energies will be γ(Δμ_O)")
    print(f"      Slab will be automatically generated for this surface!")

    # ===== CREATE AND SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph with automatic slab generation...")
    print(f"{'='*80}")

    miller_str = f"{defaults['miller_indices'][0]}{defaults['miller_indices'][1]}{defaults['miller_indices'][2]}"
    
    wg = build_core_workgraph_with_map(
        **defaults,  # Unpack all default parameters (NO input_slabs provided!)
        name=f"Ag2O_AutoGen_{miller_str}",
    )

    # Optional: Export to HTML
    try:
        html_file = f"ag2o_autogen_{miller_str}.html"
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
    print(f"EXPECTED OUTPUTS (BINARY OXIDE)")
    print(f"{'-'*80}")
    print(f"  Formation enthalpy:")
    print(f"    - formation_enthalpy (Dict with ΔH_f for Ag2O)")
    print(f"\n  Generated slabs:")
    print(f"    - slab_structures.term_0, term_1, ... (StructureData)")
    print(f"    Each termination will be generated automatically")
    print(f"\n  Relaxed slabs:")
    print(f"    - relaxed_slabs.term_0, term_1, ... (StructureData)")
    print(f"    - slab_energies.term_0, term_1, ... (Float)")
    print(f"\n  Surface energies (binary oxide thermodynamics):")
    print(f"    - surface_energies.term_0, term_1, ... (Dict with γ(Δμ_O))")
    print(f"    Surface energy as function of oxygen chemical potential")

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
    print(f"    # Get surface energies (binary oxide)")
    print(f"    surface_energies = wg.outputs.surface_energies")
    print(f"    for term_id in surface_energies.keys():")
    print(f"        data = surface_energies[term_id].get_dict()")
    print(f"        print(f'{{term_id}}: γ(Δμ_O) = {{data[\"gamma_array\"]}} J/m²')")
    print(f"    ")
    print(f"    # Export relaxed slab to file")
    print(f"    relaxed_term_0 = relaxed['term_0']")
    print(f"    atoms = relaxed_term_0.get_ase()")
    print(f"    atoms.write('ag2o_relaxed_term_0.cif')")

    print(f"\n{'='*80}")
    print(f"IMPORTANT NOTES FOR BINARY OXIDES")
    print(f"{'='*80}")
    print(f"  • Ag2O is a BINARY oxide (only metal and oxygen)")
    print(f"  • Nonmetal reference uses Ag parameters (dummy)")
    print(f"  • Surface energies are γ(Δμ_O) - function of oxygen chemical potential")
    print(f"  • Formation enthalpy is calculated for Ag2O")
    print(f"\n{'='*80}\n")

    return wg


if __name__ == "__main__":
    """
    Run the automatic slab generation and relaxation workflow for Ag2O.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Ensure structure files exist:
       ls /home/thiagotd/git/PS-TEROS/examples/structures/
       Should contain: ag2o.cif, Ag.cif, O2.cif

    5. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    6. Run this script:
       source ~/envs/psteros/bin/activate && python slabs_autogen_ag2o.py
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
