#!/home/thiagotd/envs/psteros/bin/python
"""
Simplified example using default_builders module for Ag3PO4 slab relaxation.

This script demonstrates how to use the teros.default_builders module to
reduce boilerplate code. All default VASP parameters are stored in the module
and can be easily overridden as needed.

Usage:
    source ~/envs/psteros/bin/activate && python slabs_input_relax_ag3po4_simplified.py
"""

from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map
from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults

def main():
    """Main function to run the slab relaxation workflow with input slabs."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define directories
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"
    slabs_dir = "/home/thiagotd/git/PS-TEROS/examples/slabs/input_structures"

    # ===== LOAD INPUT SLABS FROM FILES =====
    print(f"\nLoading input slabs from: {slabs_dir}")
    import os
    from pathlib import Path

    slab_files = sorted(Path(slabs_dir).glob("slab_term_*.cif"))
    if not slab_files:
        raise FileNotFoundError(
            f"No slab files found in {slabs_dir}. "
            "Expected files like: slab_term_0.cif, slab_term_1.cif, ..."
        )

    input_slabs = {}
    for slab_file in slab_files:
        term_id = slab_file.stem.replace("slab_", "")
        atoms = read(str(slab_file))
        input_slabs[term_id] = orm.StructureData(ase=atoms)
        print(f"  ✓ Loaded {slab_file.name} as '{term_id}'")

    print(f"\nTotal slabs loaded: {len(input_slabs)}")

    # ===== GET DEFAULT BUILDERS FOR Ag3PO4 =====
    # This replaces ~200 lines of parameter definitions!
    defaults = get_ag3po4_defaults(
        structures_dir=structures_dir,
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
    )

    # ===== OPTIONAL: OVERRIDE SPECIFIC PARAMETERS =====
    # Example: Change ENCUT for bulk calculation
    # defaults['bulk_parameters']['ENCUT'] = 600
    
    # Example: Use more machines for slab calculations
    # defaults['slab_options']['resources']['num_machines'] = 2
    
    # Example: Change convergence criteria
    # defaults['slab_parameters']['EDIFFG'] = -0.01

    # ===== CREATE AND SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph with default builders...")
    print(f"{'='*80}")

    wg = build_core_workgraph_with_map(
        **defaults,  # Unpack all default parameters
        #input_slabs=input_slabs,  # Add run-specific parameters
        name="Ag3PO4_InputSlabs_Relax_Simplified",
    )

    # Optional: Export to HTML
    try:
        html_file = "ag3po4_input_slabs_relax_simplified.html"
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
       source ~/envs/psteros/bin/activate && python slabs_input_relax_ag3po4_simplified.py
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
