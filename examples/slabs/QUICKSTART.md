# Quick Start: User-Provided Slabs

## 5-Minute Guide to Using Your Own Slab Structures

### Step 1: Prepare Your Slab Files

Save your slab structures in any format ASE supports (CIF, POSCAR, xyz, etc.):

```
my_slabs/
├── slab_term_0.cif
├── slab_term_1.cif
└── slab_term_2.cif
```

### Step 2: Load Structures in Python

```python
from aiida import load_profile, orm
from ase.io import read

load_profile()

# Load your slabs
input_slabs = {}
slab_files = ["slab_term_0.cif", "slab_term_1.cif", "slab_term_2.cif"]

for idx, filename in enumerate(slab_files):
    atoms = read(f"my_slabs/{filename}")
    structure = orm.StructureData(ase=atoms)
    structure.store()  # IMPORTANT: Store before using!
    input_slabs[f"term_{idx}"] = structure
```

### Step 3: Build and Submit Workflow

```python
from teros.core.workgraph import build_core_workgraph_with_map

wg = build_core_workgraph_with_map(
    structures_dir="/path/to/bulk/structures",
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"P": "P"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters=bulk_params,
    bulk_options=bulk_opts,
    metal_parameters=metal_params,
    metal_options=metal_opts,
    nonmetal_parameters=nonmetal_params,
    nonmetal_options=nonmetal_opts,
    oxygen_parameters=oxygen_params,
    oxygen_options=oxygen_opts,
    
    # KEY: Provide your slabs here!
    input_slabs=input_slabs,
    
    # Relax them
    relax_slabs=True,
    slab_parameters=slab_params,
    slab_options=slab_opts,
    slab_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    
    name="MyProject_CustomSlabs",
)

wg.submit(wait=False)
print(f"Submitted! PK: {wg.pk}")
```

### Step 4: Get Results

```python
from aiida import load_node

wg = load_node(PK)  # Replace PK with your workflow PK

# Access relaxed slabs
relaxed = wg.outputs.relaxed_slabs
energies = wg.outputs.slab_energies

for key in relaxed.keys():
    structure = relaxed[key]
    energy = energies[key].value
    print(f"{key}: {energy:.3f} eV")
    
    # Export to file
    atoms = structure.get_ase()
    atoms.write(f"relaxed_{key}.cif")
```

## That's It!

You've just run PS-TEROS with your own custom slab structures! 🎉

## Key Points

✅ **No slab generation parameters needed** when using `input_slabs`
✅ **Works with any file format** ASE can read
✅ **Same outputs** as automatic generation mode
✅ **Full backward compatibility** - existing scripts still work

## Next Steps

- See `slabs_input_relax.py` for a complete example
- Read `README_INPUT_SLABS.md` for detailed documentation
- Run `compare_modes.py` to see both modes side-by-side

## Common Patterns

### Pattern 1: Load from Directory

```python
import glob
from pathlib import Path

input_slabs = {}
for idx, filepath in enumerate(sorted(glob.glob("slabs/*.cif"))):
    atoms = read(filepath)
    structure = orm.StructureData(ase=atoms)
    structure.store()
    input_slabs[f"term_{idx}"] = structure
```

### Pattern 2: Descriptive Keys

```python
input_slabs = {}
for name, file in [("oxygen_rich", "o_rich.cif"), 
                    ("metal_rich", "m_rich.cif"),
                    ("stoichiometric", "stoich.cif")]:
    structure = orm.StructureData(ase=read(file))
    structure.store()
    input_slabs[name] = structure
```

### Pattern 3: Mixed Sources

```python
from pymatgen.core import Structure

input_slabs = {}

# From file
struct1 = orm.StructureData(ase=read("slab_0.cif"))
struct1.store()
input_slabs["term_0"] = struct1

# From Pymatgen
pmg_struct = Structure.from_file("slab_1.json")
struct2 = orm.StructureData(pymatgen=pmg_struct)
struct2.store()
input_slabs["term_1"] = struct2

# From ASE Atoms
from ase import Atoms
atoms = Atoms(...)  # Your custom ASE structure
struct3 = orm.StructureData(ase=atoms)
struct3.store()
input_slabs["term_2"] = struct3
```

## Troubleshooting

**Q: Import error when loading structure?**
```python
# Make sure ASE is installed and file exists
try:
    atoms = read("slab.cif")
except Exception as e:
    print(f"Error: {e}")
```

**Q: How do I check my structures before running?**
```python
# Visualize with ASE
from ase.visualize import view
atoms = read("slab.cif")
view(atoms)

# Or export to different format
atoms.write("slab.xyz")
```

**Q: Can I use different naming conventions?**
Yes! Use any keys you want:
```python
input_slabs = {
    "my_custom_name_1": structure1,
    "another_name": structure2,
}
```

## Complete Minimal Example

```python
#!/usr/bin/env python
from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map

load_profile()

# Load slabs and store them
input_slabs = {}
for idx, filename in enumerate(["slab_0.cif", "slab_1.cif"]):
    structure = orm.StructureData(ase=read(filename))
    structure.store()  # Must store!
    input_slabs[f"term_{idx}"] = structure

# Build workflow
wg = build_core_workgraph_with_map(
    structures_dir="structures/",
    bulk_name="bulk.cif",
    metal_name="metal.cif",
    nonmetal_name="nonmetal.cif",
    oxygen_name="O2.cif",
    code_label="vasp@localhost",
    potential_family="PBE",
    bulk_potential_mapping={"X": "X", "Y": "Y", "O": "O"},
    metal_potential_mapping={"X": "X"},
    nonmetal_potential_mapping={"Y": "Y"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters={...},
    bulk_options={...},
    metal_parameters={...},
    metal_options={...},
    nonmetal_parameters={...},
    nonmetal_options={...},
    oxygen_parameters={...},
    oxygen_options={...},
    input_slabs=input_slabs,  # Your slabs!
    relax_slabs=True,
    slab_parameters={...},
    slab_options={...},
)

wg.submit()
```

Done! 🚀
