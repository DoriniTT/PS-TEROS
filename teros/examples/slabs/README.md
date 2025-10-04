# Slab Generation Example

This example demonstrates how to use PS-TEROS to generate surface slab structures from a relaxed bulk oxide.

## What This Example Does

1. **Relaxes bulk structure** - Performs a full DFT relaxation of the Ag₃PO₄ bulk structure
2. **Calculates formation enthalpy** - Relaxes reference structures (Ag, P, O₂) and calculates ΔH_f
3. **Generates slab structures** - Creates all unique terminations for the (100) orientation
4. **Exports orthogonal slabs** - Each termination is output as a c-axis orthogonal structure ready for surface calculations

## Quick Start

```bash
# 1. Set up AiiDA profile
verdi profile set-default psteros
verdi status

# 2. Start the daemon if needed
verdi daemon start

# 3. Run the example
source ~/envs/aiida/bin/activate
cd /home/thiagotd/git/PS-TEROS/teros/examples/slabs
python slabs.py
```

## Key Parameters

### Slab Generation

- **miller_indices**: `[1, 0, 0]` - Surface orientation (hkl)
- **min_slab_thickness**: `10.0` Å - Minimum thickness of the slab
- **min_vacuum_thickness**: `15.0` Å - Minimum vacuum between periodic images
- **symmetrize**: `False` - Generate all unique terminations (not just symmetric ones)
- **primitive**: `True` - Use primitive cell before generating slabs
- **lll_reduce**: `False` - Apply LLL reduction for more orthogonal cells
- **center_slab**: `True` - Center the slab in the c direction

### Customization

To generate slabs for a different orientation, modify:

```python
miller_indices = [1, 1, 0]  # (110) surface
```

To change the slab size:

```python
min_slab_thickness = 15.0   # Thicker slab
min_vacuum_thickness = 20.0  # More vacuum
```

## Output Structure

After the workflow completes, you can access:

### Slab Structures

```python
from aiida import load_node

wg = load_node(PK)  # Replace PK with your WorkGraph PK

# List all generated terminations
slabs = wg.outputs.slab_structures
print(list(slabs.keys()))  # ['term_0', 'term_1', 'term_2', ...]

# Access individual termination
term_0 = slabs['term_0']
atoms = term_0.get_ase()

# Export to file
atoms.write('ag3po4_100_term_0.cif')
atoms.write('ag3po4_100_term_0.vasp')
```

### Formation Enthalpy

```python
# Get formation enthalpy results
hf = wg.outputs.formation_enthalpy.get_dict()
print(f"ΔH_f = {hf['formation_enthalpy_ev']:.3f} eV/f.u.")
print(f"ΔH_f = {hf['formation_enthalpy_kjmol']:.3f} kJ/mol")
```

## Workflow Details

The workflow executes the following tasks in parallel:

1. **Bulk relaxation** → VASP relaxation with ISIF=3
2. **Metal relaxation** → Ag with metallic smearing (ISMEAR=1)
3. **Nonmetal relaxation** → P with Gaussian smearing (ISMEAR=0)
4. **Oxygen relaxation** → O₂ molecule with ISIF=2

Then sequentially:

5. **Formation enthalpy calculation** → Uses results from steps 1-4
6. **Slab generation** → Uses relaxed bulk from step 1
   - Converts to primitive cell (if primitive=True)
   - Generates all unique terminations using Pymatgen's SlabGenerator
   - Converts each slab to orthogonal cell with c-axis perpendicular to surface
   - Outputs each as a separate StructureData node

## Monitoring

```bash
# Check workflow status
verdi process show <PK>

# Check detailed report
verdi process report <PK>

# List all processes
verdi process list -a -p1

# Watch in real-time
watch -n 5 'verdi process list -a -p1'
```

## Typical Output

For Ag₃PO₄ (100) with `symmetrize=False`, you might get:
- **term_0**: Ag-terminated surface
- **term_1**: P-terminated surface
- **term_2**: O-terminated surface
- **term_3**: Mixed termination

Each termination is automatically:
- ✓ Orthogonal with c-axis perpendicular to surface
- ✓ Centered in the c direction
- ✓ Ready for surface energy calculations

## Next Steps

After generating slabs, you can:

1. **Visualize structures** - Use ASE, VESTA, or Materials Studio
2. **Relax slab structures** - Perform VASP relaxations with fixed bottom layers
3. **Calculate surface energies** - Using the formation enthalpy as reference
4. **Screen different orientations** - Run for (110), (111), etc.

## File Structure

```
examples/slabs/
├── README.md          # This file
└── slabs.py          # Main example script
```

## Troubleshooting

**Issue**: No slabs generated
- Check that the bulk structure relaxed successfully
- Verify Miller indices are valid for your structure
- Try increasing `max_normal_search` parameter

**Issue**: Too many terminations
- Set `symmetrize=True` to get only symmetrically distinct terminations
- Increase `min_slab_thickness` to reduce number of possible cuts

**Issue**: Slabs not orthogonal enough
- Set `lll_reduce=True` to apply LLL reduction algorithm
- Try different Miller indices

## References

- [Pymatgen SlabGenerator Documentation](https://pymatgen.org/pymatgen.core.surface.html)
- [AiiDA-WorkGraph Documentation](https://aiida-workgraph.readthedocs.io/)
- [PS-TEROS Documentation](../../CLAUDE.md)
