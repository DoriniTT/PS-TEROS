# Slab Examples Directory

This directory contains examples for working with surface slab structures in PS-TEROS.

## Quick Navigation

### 🚀 New User? Start Here!
- **[QUICKSTART.md](QUICKSTART.md)** - 5-minute guide to get started

### 📚 Full Documentation
- **[README_INPUT_SLABS.md](README_INPUT_SLABS.md)** - Comprehensive guide for user-provided slabs
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Technical implementation details

### 💻 Example Scripts

#### User-Provided Slabs (New Feature)
- **[slabs_input_relax.py](slabs_input_relax.py)** - Use your own pre-generated slab structures
  - Load slabs from files (CIF, POSCAR, etc.)
  - Skip automatic slab generation
  - Full flexibility and control

#### Automatic Generation (Original)
- **[slabs_relax.py](slabs_relax.py)** - Generate slabs automatically from bulk
  - Uses Pymatgen's SlabGenerator
  - Specify Miller indices and thickness
  - Standard approach for systematic studies

#### Comparison
- **[compare_modes.py](compare_modes.py)** - Compare both modes side-by-side
  - Educational demonstration
  - Shows API differences
  - Explains when to use each mode

## Two Ways to Work with Slabs

### Mode 1: Automatic Generation (Traditional)
```python
wg = build_core_workgraph_with_map(
    ...,
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
)
```

### Mode 2: User-Provided Slabs (New)
```python
# Load your slabs
input_slabs = {}
for idx, file in enumerate(["slab_0.cif", "slab_1.cif"]):
    structure = orm.StructureData(ase=read(file))
    structure.store()  # Must be stored!
    input_slabs[f"term_{idx}"] = structure

wg = build_core_workgraph_with_map(
    ...,
    input_slabs=input_slabs,  # Provide your slabs
    relax_slabs=True,
)
```

## Directory Structure

```
examples/slabs/
├── README.md                      # This file
├── QUICKSTART.md                  # Quick start guide
├── README_INPUT_SLABS.md          # Detailed user guide
├── IMPLEMENTATION_SUMMARY.md      # Technical documentation
│
├── slabs_input_relax.py          # Example: User-provided slabs
├── slabs_relax.py                # Example: Automatic generation
├── compare_modes.py              # Comparison demonstration
│
└── input_structures/             # Directory for your slab files
    └── README.md                 # Instructions
```

## When to Use Each Mode

### Use Automatic Generation When:
- ✅ You want standard surface terminations
- ✅ You're exploring multiple Miller indices
- ✅ You want all symmetrically distinct terminations
- ✅ You're doing systematic surface studies

### Use User-Provided Slabs When:
- ✅ You have specific structures from literature
- ✅ You need custom surface reconstructions
- ✅ You want to add adsorbates or defects
- ✅ You've manually edited slab structures
- ✅ You're using specialized generation tools

## Getting Started

1. **New to PS-TEROS slabs?**
   - Start with [QUICKSTART.md](QUICKSTART.md)
   - Try the automatic generation example first: `slabs_relax.py`

2. **Want to use your own slabs?**
   - Read [README_INPUT_SLABS.md](README_INPUT_SLABS.md)
   - Check the example: `slabs_input_relax.py`
   - Place your slab files in `input_structures/`

3. **Not sure which mode to use?**
   - Run `python compare_modes.py` to see both modes
   - Read the comparison guide in output

## Common Workflows

### Workflow 1: Systematic Surface Study
```bash
# Use automatic generation
python slabs_relax.py
# Generates all terminations for specified Miller index
# Relaxes all slabs in parallel
```

### Workflow 2: Specific Surface Termination
```bash
# Use your own slab structure
# 1. Create/obtain slab structure file
# 2. Place in input_structures/
# 3. Run:
python slabs_input_relax.py
```

### Workflow 3: Literature Reproduction
```bash
# Use structure from paper
# 1. Download/create exact structure from publication
# 2. Load in script as input_slabs
# 3. Submit workflow
```

## Output Structure

Both modes produce the same outputs:

```python
from aiida import load_node

wg = load_node(PK)

# Slab structures (generated or provided)
slabs = wg.outputs.slab_structures
# → {'term_0': StructureData, 'term_1': StructureData, ...}

# Relaxed slabs (if relax_slabs=True)
relaxed = wg.outputs.relaxed_slabs
# → {'term_0': StructureData, 'term_1': StructureData, ...}

# Energies (if relax_slabs=True)
energies = wg.outputs.slab_energies
# → {'term_0': Float, 'term_1': Float, ...}
```

## Requirements

- AiiDA profile configured
- VASP code set up in AiiDA
- Pseudopotentials installed
- ASE for structure file reading

## Tips

1. **File Formats**: Use CIF for easy visualization and sharing
2. **Naming**: Use descriptive names like `ag3po4_100_term_0.cif`
3. **Validation**: Check structures before submitting workflows
4. **Vacuum**: Ensure adequate vacuum spacing (10-15 Å typical)
5. **Starting Point**: Use automatic generation to create initial structures

## Troubleshooting

**Q: My slabs aren't loading?**
- Check file path and format
- Verify ASE can read the file: `ase.io.read("file.cif")`

**Q: Which mode should I use?**
- Use automatic for standard studies
- Use user-provided for custom configurations

**Q: Can I modify auto-generated slabs?**
- Yes! Generate first, export, modify, then use as input

**Q: How do I visualize slabs?**
```python
from ase.io import read
from ase.visualize import view
atoms = read("slab.cif")
view(atoms)
```

## Support

- Full documentation: [README_INPUT_SLABS.md](README_INPUT_SLABS.md)
- Technical details: [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)
- Quick reference: [QUICKSTART.md](QUICKSTART.md)
- Main PS-TEROS docs: `../../docs/`

## Version

This feature was added in PS-TEROS v0.2.0

See `../../CHANGE.md` for full changelog.
