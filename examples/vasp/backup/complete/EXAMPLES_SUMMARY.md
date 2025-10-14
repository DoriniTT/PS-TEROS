# PS-TEROS Examples Summary

This document provides an overview of the available examples for the PS-TEROS workflow system after the recent updates that removed the "bulk only mode" and replaced it with explicit boolean flags for controlling calculations.

## Overview of Changes

The PS-TEROS `build_core_workgraph` function has been updated to:

1. **Remove implicit "modes"** - No more MODE 1 (bulk-only) vs MODE 2 (formation enthalpy)
2. **Add explicit control flags** - All calculations controlled by clear boolean parameters
3. **Set sensible defaults** - Everything calculates by default when inputs are provided
4. **Minimize required inputs** - Many parameters have reasonable default values

## Available Examples

### 1. Simple Example: Bulk Relaxation Only

**Location**: `examples/relaxation/relaxation.py`

**Purpose**: Minimal example showing bulk relaxation without reference structures.

**Features Tested**:
- ✓ Bulk structure relaxation
- ✗ Reference relaxations
- ✗ Formation enthalpy
- ✗ Slab generation/relaxation
- ✗ Relaxation energies
- ✗ Cleavage energies
- ✗ Surface thermodynamics

**Material**: Ag₃PO₄

**Required Inputs**:
```python
- structures_dir
- bulk_name
- bulk_potential_mapping
- bulk_parameters
- bulk_options
```

**Optional Inputs**: All others (with defaults)

**Estimated Runtime**: ~1 hour

**Key Points**:
- Shows minimal usage pattern
- No metal_name/oxygen_name provided → thermodynamics skipped automatically
- Demonstrates default values for flags

**Usage**:
```bash
cd examples/relaxation
source ~/envs/aiida/bin/activate
python relaxation.py
```

---

### 2. Complete Example: Full Workflow Testing

**Location**: `examples/complete/complete_example.py`

**Purpose**: Comprehensive example testing **ALL** features of PS-TEROS.

**Features Tested**:
- ✓ Bulk structure relaxation
- ✓ Reference relaxations (metal, nonmetal, oxygen)
- ✓ Formation enthalpy calculation
- ✓ Slab generation from bulk
- ✓ Slab relaxation with unrelaxed SCF
- ✓ Relaxation energies (E_relaxed - E_unrelaxed)
- ✓ Cleavage energies for complementary slabs
- ✓ Surface thermodynamics with chemical potential sampling

**Material**: Ag₃PO₄ (ternary oxide)

**Required Inputs**:
```python
- structures_dir
- bulk_name, metal_name, nonmetal_name, oxygen_name
- bulk_potential_mapping, metal_potential_mapping, etc.
- bulk_parameters, metal_parameters, etc.
- bulk_options, metal_options, etc.
- miller_indices (for slab generation)
```

**Optional Inputs**:
- min_slab_thickness (default: 15.0 Å)
- min_vacuum_thickness (default: 15.0 Å)
- compute_relaxation_energy (default: True)
- compute_cleavage (default: True)
- compute_thermodynamics (default: True)
- thermodynamics_sampling (default: 100)

**Estimated Runtime**: Several hours to a day

**Key Points**:
- Production-quality example
- Demonstrates all calculation flags
- Shows proper VASP parameters for different material types
- Includes detailed monitoring instructions
- Creates HTML visualization of workflow

**Usage**:
```bash
cd examples/complete
source ~/envs/aiida/bin/activate
python complete_example.py
```

---

## Comparison Table

| Feature | Simple Example | Complete Example |
|---------|---------------|------------------|
| **Purpose** | Minimal usage | Full capabilities |
| **Bulk relaxation** | ✓ | ✓ |
| **Reference relaxations** | ✗ | ✓ (Metal, Nonmetal, O₂) |
| **Formation enthalpy** | ✗ | ✓ |
| **Slab generation** | ✗ | ✓ |
| **Slab relaxation** | ✗ | ✓ |
| **Unrelaxed SCF** | ✗ | ✓ |
| **Relaxation energies** | ✗ | ✓ |
| **Cleavage energies** | ✗ | ✓ |
| **Surface thermodynamics** | ✗ | ✓ |
| **Chemical potential sampling** | ✗ | ✓ |
| **Material** | Ag₃PO₄ | Ag₃PO₄ |
| **Lines of code** | ~150 | ~530 |
| **Runtime estimate** | ~1 hour | Several hours |
| **Best for** | Testing, learning | Production research |

---

## Control Flags Reference

All examples can control calculations using these boolean flags:

### `relax_slabs` (default: False)
- Controls whether to perform slab relaxations
- If False, only bulk/reference relaxations are done
- Required for: relaxation_energy, cleavage, thermodynamics

### `compute_relaxation_energy` (default: True)
- Calculate energy difference between unrelaxed and relaxed slabs
- Requires: relax_slabs=True
- Outputs: relaxation_energies namespace

### `compute_cleavage` (default: True)
- Calculate cleavage energies for complementary slab pairs
- Requires: relax_slabs=True
- Outputs: cleavage_energies namespace

### `compute_thermodynamics` (default: True)
- Calculate surface energies with chemical potential sampling
- Requires: relax_slabs=True, metal_name, oxygen_name
- Outputs: surface_energies namespace
- Controls: thermodynamics_sampling parameter (grid resolution)

---

## Default Values

The following parameters have sensible defaults:

| Parameter | Default Value | Description |
|-----------|--------------|-------------|
| `code_label` | `'VASP-VTST-6.4.3@bohr'` | VASP code in AiiDA |
| `potential_family` | `'PBE'` | Pseudopotential family |
| `kpoints_spacing` | `0.3` | K-points spacing (Å⁻¹ × 2π) |
| `min_slab_thickness` | `15.0` | Minimum slab thickness (Å) |
| `min_vacuum_thickness` | `15.0` | Minimum vacuum thickness (Å) |
| `clean_workdir` | `False` | Clean work directory after completion |
| `relax_slabs` | `False` | Perform slab relaxations |
| `compute_relaxation_energy` | `True` | Calculate relaxation energies |
| `compute_cleavage` | `True` | Calculate cleavage energies |
| `compute_thermodynamics` | `True` | Calculate surface energies |
| `thermodynamics_sampling` | `100` | Grid points for μ sampling |
| `lll_reduce` | `False` | LLL cell reduction |
| `center_slab` | `True` | Center slab in c direction |
| `symmetrize` | `False` | Symmetrize terminations |
| `primitive` | `True` | Find primitive cell |

---

## Which Example Should I Use?

### Use the **Simple Example** if you want to:
- Learn the basics of PS-TEROS
- Test your AiiDA/VASP setup
- Only need bulk relaxation
- Want a quick calculation for testing
- Understand minimal required inputs

### Use the **Complete Example** if you want to:
- Perform full surface energy calculations
- Calculate formation enthalpies
- Study surface thermodynamics
- Compare different terminations
- Perform production research calculations
- Understand all PS-TEROS capabilities

---

## Example Directory Structure

```
examples/
├── EXAMPLES_SUMMARY.md          # This file
├── structures/                  # Shared structure files
│   ├── ag3po4.cif
│   ├── Ag.cif
│   ├── P.cif
│   └── O2.cif
├── relaxation/                  # Simple example
│   ├── relaxation.py
│   └── structures/
│       └── ag3po4.cif
└── complete/                    # Complete example
    ├── README.md                # Detailed documentation
    ├── complete_example.py      # Main script
    └── structures/              # All structure files
        ├── ag3po4.cif
        ├── Ag.cif
        ├── P.cif
        └── O2.cif
```

---

## Output Structure

Both examples produce outputs in the same format, but the complete example produces more:

### Simple Example Outputs:
```
bulk_energy: Float
bulk_structure: StructureData
metal_energy: Float (placeholder)
metal_structure: StructureData (placeholder)
... (other placeholders)
```

### Complete Example Outputs:
```
# Bulk and References
bulk_energy: Float
bulk_structure: StructureData
metal_energy: Float
metal_structure: StructureData
nonmetal_energy: Float
nonmetal_structure: StructureData
oxygen_energy: Float
oxygen_structure: StructureData

# Formation Enthalpy
formation_enthalpy: Dict
    - delta_hf: Float (eV/formula unit)
    - stoichiometry: Dict
    - energies: Dict

# Slabs
slab_structures: namespace
    - term_0: StructureData
    - term_1: StructureData
    - ...

# Slab Energies
unrelaxed_slab_energies: namespace
slab_energies: namespace
relaxed_slabs: namespace

# Derived Properties
relaxation_energies: namespace
    - term_0: Float (eV)
    - ...

cleavage_energies: namespace
    - pair_0: Dict
        - E_cleavage: Float (J/m²)
        - terminations: List[str]
    - ...

surface_energies: namespace
    - term_0: Dict
        - gamma: Array[n_O, n_P] (J/m²)
        - mu_O_range: Array
        - mu_P_range: Array
        - stable_region: Array[bool]
    - ...
```

---

## Monitoring Commands

For both examples, use these commands to monitor progress:

```bash
# Check overall workflow status
verdi process show <PK>

# Check detailed report
verdi process report <PK>

# List all running processes
verdi process list

# Watch processes in real-time
watch -n 5 'verdi process list'

# Check specific VASP calculation output
verdi calcjob outputcat <CALCJOB_PK>
```

---

## Next Steps

After running these examples:

1. **Modify parameters** - Adjust VASP parameters for your system
2. **Change materials** - Use your own structure files
3. **Control calculations** - Toggle flags to compute only what you need
4. **Scale up** - Use the patterns for production calculations
5. **Extend** - Build upon these examples for custom workflows

---

## Additional Resources

- **PS-TEROS Documentation**: [Link to docs]
- **AiiDA Documentation**: https://aiida.readthedocs.io/
- **WorkGraph Documentation**: https://aiida-workgraph.readthedocs.io/

---

## Troubleshooting

### Python cache issues
```bash
find . -type d -name __pycache__ -exec rm -rf {} +
find . -name "*.pyc" -delete
verdi daemon restart
```

### Check installation
```bash
verdi profile show
verdi status
verdi code list
```

### Failed calculations
```bash
verdi process show <PK>
verdi calcjob res <CALCJOB_PK>
verdi calcjob outputcat <CALCJOB_PK>
```

---

## Contributing

If you create additional examples, please:
1. Document them clearly
2. Test thoroughly
3. Add to this summary
4. Follow the established naming conventions

---

*Last updated: 2025-10-11*
