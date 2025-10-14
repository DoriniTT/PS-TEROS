#!/usr/bin/env python
"""
Comparison example showing both automatic slab generation and user-provided slabs.

This script demonstrates:
1. How to use the automatic slab generation mode (traditional)
2. How to use the user-provided slabs mode (new feature)
3. The differences in parameter requirements for each mode

This is a demonstration script that shows the API differences.
It does NOT submit calculations - use for reference only.

Usage:
    python compare_modes.py
"""

from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map


def show_automatic_generation_mode():
    """Demonstrate automatic slab generation mode (traditional approach)."""
    
    print("\n" + "="*80)
    print("MODE 1: AUTOMATIC SLAB GENERATION (Traditional)")
    print("="*80)
    
    print("\nIn this mode, PS-TEROS generates slabs from the relaxed bulk structure.")
    print("You must provide slab generation parameters.")
    
    print("\nRequired parameters:")
    print("  - miller_indices: [1, 0, 0]")
    print("  - min_slab_thickness: 10.0")
    print("  - min_vacuum_thickness: 15.0")
    
    print("\nOptional generation parameters:")
    print("  - lll_reduce: True/False")
    print("  - center_slab: True/False")
    print("  - symmetrize: True/False")
    print("  - primitive: True/False")
    
    print("\nExample code:")
    print("-" * 80)
    print("""
wg = build_core_workgraph_with_map(
    structures_dir="/path/to/structures",
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    # ... potential mappings, parameters, options ...
    
    # REQUIRED: Slab generation parameters
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    
    # OPTIONAL: Slab generation controls
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
    
    # Slab relaxation
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options=slab_options,
    
    name="AutoGeneration_Example",
)
    """)
    print("-" * 80)
    
    print("\nWorkflow steps:")
    print("  1. Relax bulk structure")
    print("  2. Generate slabs from relaxed bulk using Pymatgen")
    print("  3. Relax all generated slabs in parallel")
    print("  4. Output relaxed slabs and energies")


def show_user_provided_mode():
    """Demonstrate user-provided slabs mode (new feature)."""
    
    print("\n" + "="*80)
    print("MODE 2: USER-PROVIDED SLABS (New Feature)")
    print("="*80)
    
    print("\nIn this mode, you provide pre-generated slab structures.")
    print("PS-TEROS skips slab generation and uses your structures directly.")
    
    print("\nRequired setup:")
    print("  1. Create/obtain slab structure files (CIF, POSCAR, etc.)")
    print("  2. Load them into a dictionary with AiiDA StructureData")
    print("  3. Pass the dictionary as 'input_slabs' parameter")
    
    print("\nNOT required:")
    print("  - miller_indices (ignored)")
    print("  - min_slab_thickness (ignored)")
    print("  - min_vacuum_thickness (ignored)")
    print("  - lll_reduce, center_slab, symmetrize, etc. (ignored)")
    
    print("\nExample code:")
    print("-" * 80)
    print("""
# Step 1: Load your pre-generated slabs
from ase.io import read

input_slabs = {}
slab_files = ["slab_0.cif", "slab_1.cif", "slab_2.cif"]

for idx, slab_file in enumerate(slab_files):
    atoms = read(f"path/to/slabs/{slab_file}")
    input_slabs[f"term_{idx}"] = orm.StructureData(ase=atoms)

# Step 2: Build workgraph with input_slabs
wg = build_core_workgraph_with_map(
    structures_dir="/path/to/structures",
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    # ... potential mappings, parameters, options ...
    
    # NEW: Provide your slabs
    input_slabs=input_slabs,
    
    # NOT NEEDED: Generation parameters (will be ignored if provided)
    # miller_indices=[1, 0, 0],      # <-- Not required
    # min_slab_thickness=10.0,       # <-- Not required
    # min_vacuum_thickness=15.0,     # <-- Not required
    
    # Slab relaxation (still works as before)
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options=slab_options,
    
    name="UserProvided_Example",
)
    """)
    print("-" * 80)
    
    print("\nWorkflow steps:")
    print("  1. Relax bulk structure")
    print("  2. Use your provided slab structures (skip generation)")
    print("  3. Relax all provided slabs in parallel")
    print("  4. Output relaxed slabs and energies")


def show_use_cases():
    """Show when to use each mode."""
    
    print("\n" + "="*80)
    print("WHEN TO USE EACH MODE")
    print("="*80)
    
    print("\nUse AUTOMATIC GENERATION when:")
    print("  ✓ You want standard surface terminations")
    print("  ✓ You're exploring multiple Miller indices")
    print("  ✓ You want all symmetrically distinct terminations")
    print("  ✓ You trust Pymatgen's slab generation")
    print("  ✓ You're doing systematic surface studies")
    
    print("\nUse USER-PROVIDED SLABS when:")
    print("  ✓ You have specific slab structures from literature")
    print("  ✓ You need custom surface reconstructions")
    print("  ✓ You want to add adsorbates or defects")
    print("  ✓ You've manually edited slab structures")
    print("  ✓ You're using specialized slab generation tools")
    print("  ✓ You need exact reproducibility from external sources")
    print("  ✓ You want to test specific surface configurations")


def show_benefits():
    """Show benefits of the new feature."""
    
    print("\n" + "="*80)
    print("BENEFITS OF USER-PROVIDED SLABS")
    print("="*80)
    
    print("\n1. FLEXIBILITY")
    print("   - Use any slab generation method you prefer")
    print("   - Create custom surface modifications")
    print("   - Reproduce structures from publications")
    
    print("\n2. CONTROL")
    print("   - Exact control over surface terminations")
    print("   - Add specific adsorbates or defects")
    print("   - Use non-standard surface orientations")
    
    print("\n3. EFFICIENCY")
    print("   - Skip generation if you already have structures")
    print("   - Faster workflow setup")
    print("   - Direct use of pre-optimized geometries")
    
    print("\n4. COMPATIBILITY")
    print("   - Works with any structure file format (via ASE)")
    print("   - Integrates with external tools")
    print("   - Maintains same output structure as automatic mode")


def main():
    """Main comparison function."""
    
    print("\n" + "="*80)
    print("PS-TEROS: SLAB INPUT MODES COMPARISON")
    print("="*80)
    print("\nThis script demonstrates the two ways to work with slabs in PS-TEROS:")
    print("  1. Automatic generation from bulk structure")
    print("  2. User-provided pre-generated slabs (NEW)")
    
    show_automatic_generation_mode()
    show_user_provided_mode()
    show_use_cases()
    show_benefits()
    
    print("\n" + "="*80)
    print("OUTPUTS (Both modes produce the same structure)")
    print("="*80)
    print("""
After workflow completion, access results the same way:

    from aiida import load_node
    
    wg = load_node(PK)
    
    # Slab structures (generated OR provided)
    slabs = wg.outputs.slab_structures  # {'term_0': ..., 'term_1': ...}
    
    # Relaxed structures (if relax_slabs=True)
    relaxed = wg.outputs.relaxed_slabs  # {'term_0': ..., 'term_1': ...}
    
    # Energies (if relax_slabs=True)
    energies = wg.outputs.slab_energies  # {'term_0': ..., 'term_1': ...}
    """)
    
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("\nFor working examples, see:")
    print("  - examples/backup/slabs/slabs_relax.py")
    print("    (automatic generation)")
    print("  - examples/slabs/slabs_input_relax.py")
    print("    (user-provided slabs)")
    print("\nFor detailed documentation, see:")
    print("  - examples/slabs/README_INPUT_SLABS.md")
    
    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    main()
