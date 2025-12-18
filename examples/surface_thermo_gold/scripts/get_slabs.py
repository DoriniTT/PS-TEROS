from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read, write
from ase.visualize import view

primitive_structure_ase = read("au.cif")

# Convert to a pymatgen structure
primitive_structure = AseAtomsAdaptor().get_structure(primitive_structure_ase)
# Load the primitive structure from the .vasp file
#primitive_structure = Poscar.from_file("ag2cro4.cif").structure

# Define the Miller indices for several low index surface orientations
miller_list = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]

# Define the minimum and maximum slab thickness
min_slab_thickness = 18 # in Angstroms
vacuum = 20  # in Angstroms

for miller_indices in miller_list:
    # Generate the slabs
    slab_generator = SlabGenerator(primitive_structure, miller_indices, min_slab_thickness, vacuum, lll_reduce=True, center_slab=True)
    slabs = slab_generator.get_slabs(symmetrize=True)
    
    h, k, l = miller_indices
    print(f"Generating slabs for surface ({h}, {k}, {l})")

    for n, slab in enumerate(slabs):
        slab = slab.get_orthogonal_c_slab()
        ase_atoms = AseAtomsAdaptor().get_atoms(slab) 
        filename = f'structure_{h}{k}{l}_{n+1}.vasp'
        write(filename, ase_atoms)
        print(f"  Saved {filename}")