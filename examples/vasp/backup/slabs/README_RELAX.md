# Slab Generation + Relaxation Example

This example demonstrates how to use PS-TEROS to generate surface slab structures from a relaxed bulk oxide **and relax them with VASP in parallel**.

## What This Example Does

1. **Relaxes bulk structure** - Performs a full DFT relaxation of the Ag₃PO₄ bulk structure
2. **Calculates formation enthalpy** - Relaxes reference structures (Ag, P, O₂) and calculates ΔH_f
3. **Generates slab structures** - Creates all unique terminations for the (100) orientation
4. **Relaxes all slabs in parallel** - Each slab termination is relaxed with VASP simultaneously
5. **Outputs relaxed structures and energies** - All relaxed slabs and their final energies

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
python slabs_relax.py
```

## Architecture

### Workflow Structure

The workflow uses a **nested @task.graph** pattern to handle dynamic slab relaxation:

1. **Main workflow** (`formation_workgraph`):
   - Relaxes bulk and references
   - Generates slabs with `get_slabs`
   - Calls `relax_all_slabs` sub-workflow

2. **Slab relaxation sub-workflow** (`relax_all_slabs`):
   - Takes dynamic slabs dict from `get_slabs`
   - Loops over each slab (term_0, term_1, ...)
   - Creates VASP relaxation task for each
   - Returns dynamic namespaces: `relaxed_slabs` and `slab_energies`

```python
@task.graph
def relax_all_slabs(
    slabs: Annotated[dict, spec.dynamic(Atoms)],
    ...
) -> Annotated[dict, spec.namespace(
    relaxed_slabs=spec.dynamic(orm.StructureData),
    slab_energies=spec.dynamic(orm.Float)
)]:
    ...
    for slab_id, slab_atoms in slabs.items():
        # Create VASP task for this slab
        vasp_task = VaspTask(...)
        energy = extract_total_energy(energies=vasp_task.misc)

        relaxed_slabs[slab_id] = vasp_task.structure
        slab_energies[slab_id] = energy.result

    return {
        'relaxed_slabs': relaxed_slabs,
        'slab_energies': slab_energies,
    }
```

### Key Parameters

#### Slab Generation
- **miller_indices**: `[1, 0, 0]` - Surface orientation (hkl)
- **min_slab_thickness**: `10.0` Å - Minimum thickness of the slab
- **min_vacuum_thickness**: `15.0` Å - Minimum vacuum between periodic images
- **symmetrize**: `True` - Generate only symmetrically distinct terminations
- **lll_reduce**: `True` - Apply LLL reduction for more orthogonal cells

#### Slab Relaxation
- **relax_slabs**: `True` - Enable slab relaxation
- **slab_parameters**: VASP INCAR parameters
  - `ISIF=2` - Relax atoms only, keep cell fixed
  - `EDIFFG=-0.02` - Tighter force convergence for surfaces
  - Optional: `IDIPOL=3`, `LDIPOL=True` for dipole corrections
- **slab_options**: Scheduler options (can be different from bulk)
- **slab_potential_mapping**: Element to potential mapping (defaults to bulk mapping)
- **slab_kpoints_spacing**: K-points spacing (defaults to bulk spacing)

### Customization

#### Different Orientation

```python
miller_indices = [1, 1, 0]  # (110) surface
```

#### Slab Size

```python
min_slab_thickness = 15.0   # Thicker slab
min_vacuum_thickness = 20.0  # More vacuum
```

#### Dipole Corrections

For asymmetric slabs, add dipole corrections:

```python
slab_parameters = {
    ...
    "IDIPOL": 3,      # Dipole correction along z-axis (c-axis)
    "LDIPOL": True,   # Turn on dipole corrections
}
```

#### Fixed Bottom Layers

To fix bottom layers (common for surface calculations), you would need to modify the structure with selective dynamics before relaxation. This can be added to the workflow if needed.

## Output Structure

After the workflow completes:

### Unrelaxed Slabs

```python
from aiida import load_node

wg = load_node(PK)

# List all generated terminations
slabs = wg.outputs.slab_structures
print(list(slabs.keys()))  # ['term_0', 'term_1', ...]

# Access individual slab
term_0 = slabs['term_0']
atoms = term_0.get_ase()
atoms.write('unrelaxed_term_0.cif')
```

### Relaxed Slabs

```python
# List all relaxed terminations
relaxed = wg.outputs.relaxed_slabs
print(list(relaxed.keys()))  # ['term_0', 'term_1', ...]

# Access relaxed slab
relaxed_term_0 = relaxed['term_0']
atoms = relaxed_term_0.get_ase()
atoms.write('relaxed_term_0.cif')
```

### Slab Energies

```python
# Get all slab energies
energies = wg.outputs.slab_energies

for term_id in energies.keys():
    energy = energies[term_id].value
    print(f'{term_id}: {energy:.6f} eV')

# Output:
# term_0: -285.123456 eV
# term_1: -284.987654 eV
```

### Formation Enthalpy

```python
# Get formation enthalpy results
hf = wg.outputs.formation_enthalpy.get_dict()
print(f"ΔH_f = {hf['formation_enthalpy_ev']:.3f} eV/f.u.")
print(f"ΔH_f = {hf['formation_enthalpy_kjmol']:.3f} kJ/mol")
```

## Workflow Details

The workflow executes the following tasks:

### Phase 1: Parallel Relaxations (tasks 1-4)
1. **Bulk relaxation** → VASP with ISIF=3
2. **Metal relaxation** → Ag with metallic smearing
3. **Nonmetal relaxation** → P with Gaussian smearing
4. **Oxygen relaxation** → O₂ molecule with ISIF=2

### Phase 2: Sequential (tasks 5-6)
5. **Formation enthalpy** → Uses results from Phase 1
6. **Slab generation** → Creates all terminations from relaxed bulk

### Phase 3: Parallel Slab Relaxations (task 7)
7. **Slab relaxations** → All slabs relaxed in parallel
   - Creates N VASP tasks (one per termination)
   - Each task relaxes its slab with ISIF=2
   - Extracts energies from all tasks

## Task Graph

```
┌─ bulk_vasp ─ extract_energy ─┐
├─ metal_vasp ─ extract_energy ─┤
├─ nonmetal_vasp ─ extract_energy ─┼─ formation_hf
├─ oxygen_vasp ─ extract_energy ─┘
│
└─ bulk_vasp ─ get_slabs ─ relax_all_slabs
                            ├─ term_0: vasp ─ extract_energy
                            ├─ term_1: vasp ─ extract_energy
                            └─ term_N: vasp ─ extract_energy
```

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

For Ag₃PO₄ (100) with `symmetrize=True`, you might get:
- **term_0**: Ag-terminated (lowest energy)
- **term_1**: P-terminated (higher energy)

Each termination is:
- ✓ Orthogonal with c-axis perpendicular to surface
- ✓ Relaxed with ISIF=2 (atoms only)
- ✓ Converged to EDIFFG=-0.02 eV/Å
- ✓ Ready for surface energy calculations

## Surface Energy Calculation

After relaxation, you can calculate surface energies:

```python
# Get bulk and slab energies
E_bulk = wg.outputs.bulk_energy.value
n_bulk_atoms = len(wg.outputs.bulk_structure.get_ase())

for term_id in wg.outputs.relaxed_slabs.keys():
    E_slab = wg.outputs.slab_energies[term_id].value
    slab_struct = wg.outputs.relaxed_slabs[term_id]
    n_slab_atoms = len(slab_struct.get_ase())

    # Surface area (assuming orthogonal cell)
    cell = slab_struct.get_ase().get_cell()
    area = cell[0][0] * cell[1][1]  # A x B for c-axis orthogonal

    # Surface energy (per surface, there are 2 surfaces)
    E_surf = (E_slab - (n_slab_atoms / n_bulk_atoms) * E_bulk) / (2 * area)

    print(f'{term_id}: {E_surf:.4f} eV/Å²')
```

## Advanced: Selective Dynamics

To implement fixed bottom layers:

```python
from ase.constraints import FixAtoms

# In a custom task
@task
def add_constraints(structure: orm.StructureData, fix_layers: int) -> orm.StructureData:
    atoms = structure.get_ase()

    # Fix bottom N layers based on z-coordinates
    z_coords = atoms.get_positions()[:, 2]
    z_min = z_coords.min()
    layer_height = 2.5  # Angstroms per layer

    fix_indices = [i for i, z in enumerate(z_coords)
                   if z < z_min + fix_layers * layer_height]

    atoms.set_constraint(FixAtoms(indices=fix_indices))
    return orm.StructureData(ase=atoms)
```

## Next Steps

After relaxing slabs, you can:

1. **Calculate surface energies** - Use energies and areas
2. **Compare terminations** - Find most stable surface
3. **Run adsorption studies** - Add molecules to relaxed slabs
4. **Calculate work functions** - From VASP LOCPOT
5. **Screen orientations** - Run for (110), (111), etc.

## File Structure

```
examples/slabs/
├── README.md              # Basic slab generation
├── README_RELAX.md        # This file (slab relaxation)
├── slabs.py              # Basic example (no relaxation)
├── slabs_relax.py        # Full example (with relaxation)
├── test_slabs.py         # Test slab generation
├── test_workflow.py      # Test workflow build
└── test_slab_relax.py    # Test slab relaxation build
```

## Troubleshooting

**Issue**: Slabs not relaxing
- Check that `relax_slabs=True` is set
- Verify slab_parameters are provided or bulk_parameters exist
- Check daemon is running: `verdi daemon status`

**Issue**: Some slabs fail to converge
- Increase NSW (max ionic steps)
- Loosen EDIFFG convergence criterion
- Check for imaginary frequencies (unstable surface)
- Add dipole corrections for polar surfaces

**Issue**: Memory issues
- Reduce number of cores: `num_cores_per_machine`
- Process slabs sequentially instead of parallel (modify workflow)
- Reduce slab thickness or vacuum

## References

- [Pymatgen SlabGenerator](https://pymatgen.org/pymatgen.core.surface.html)
- [AiiDA-WorkGraph Dynamic Namespaces](https://aiida-workgraph.readthedocs.io/)
- [VASP Surface Calculations](https://www.vasp.at/wiki/index.php/Category:Surface_science)
- [PS-TEROS Documentation](../../CLAUDE.md)
