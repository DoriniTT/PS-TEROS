import sys
from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.visualize import view
from ase.io import write

def generate_slabs_from_cif(
    cif_path,
    miller_indices=(1, 0, 0),
    min_slab_thickness=8.0,
    min_vacuum_thickness=15.0,
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
    max_normal_search=None,
    in_unit_planes=False,
    max_terminations=None,
    supercell=(1, 1, 1)
):
    # Load bulk structure from CIF
    bulk_structure = Structure.from_file(str(cif_path))

    # Optionally find primitive cell
    if primitive:
        analyzer = SpacegroupAnalyzer(bulk_structure)
        bulk_structure = analyzer.get_primitive_standard_structure()

    # Create slab generator
    slab_gen = SlabGenerator(
        bulk_structure,
        miller_indices,
        min_slab_thickness,
        min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        max_normal_search=max_normal_search,
        in_unit_planes=in_unit_planes
    )

    # Generate slabs with different terminations
    slabs = slab_gen.get_slabs(symmetrize=symmetrize)

    # Limit the number of terminations if specified
    if max_terminations and len(slabs) > max_terminations:
        slabs = slabs[:max_terminations]

    # Only take the first two slabs for testing
    #slabs = slabs[:1]

    # Convert all slabs to orthogonal cells and then to ASE Atoms
    ase_slabs = []
    for i, slab in enumerate(slabs):
        ortho_slab = slab.get_orthogonal_c_slab()
        super_slab = ortho_slab.make_supercell(supercell)
        ase_atoms = AseAtomsAdaptor().get_atoms(super_slab)
        ase_slabs.append((f"s_{i}", ase_atoms))
    return ase_slabs

if __name__ == "__main__":
    # Usage: python test_slab_generation.py /path/to/structure.cif [output_folder]
    if len(sys.argv) < 2:
        print("Usage: python test_slab_generation.py /path/to/structure.cif [output_folder]")
        sys.exit(1)

    cif_path = Path(sys.argv[1])
    output_folder = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")

    output_folder.mkdir(parents=True, exist_ok=True)
    slabs = generate_slabs_from_cif(cif_path)

    for name, ase_atoms in slabs:
        print(f"Visualizing and saving: {name}")
        # Visualize using ASE GUI (close window to continue)
        view(ase_atoms)
        # Save to VASP POSCAR format for inspection
        out_file = output_folder / f"{name}.vasp"
        write(str(out_file), ase_atoms, format="vasp")
        print(f"Saved {out_file}")
