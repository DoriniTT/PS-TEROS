# PS-TEROS: Predicting Stability of TERminations of Oxide Surfaces

PS-TEROS is a Python package for calculating surface Gibbs free energy using ab initio atomistic thermodynamics with the AiiDA-WorkGraph framework. It automates the calculation of surface energies for binary and ternary oxides under varying thermodynamic conditions.

## Features

- **Automated Surface Energy Calculations**: Complete workflow for surface Gibbs free energy calculations using ab initio atomistic thermodynamics
- **DFT Code Support**: Currently supports VASP via AiiDA plugins (CP2K and Quantum ESPRESSO support planned for future releases)
- **Flexible Slab Generation**: Automatic generation or manual input of surface terminations
- **Cleavage Energy Calculations**: Built-in module for computing cleavage energies of complementary surface pairs
- **Restart Capability**: Resume incomplete calculations from previous runs using remote data
- **Default Builders**: Pre-configured parameter sets for rapid workflow setup
- **Full Provenance Tracking**: Complete data provenance via AiiDA infrastructure
- **Parallel Execution**: Efficient parallel processing of multiple terminations and references

> **Note**: The current version only supports VASP. CP2K and Quantum ESPRESSO are planned as future updates. If you need PS-TEROS code for these versions, please use the [legacy version](legacy/).

## Installation

### Requirements

- Python 3.9+
- AiiDA Core (latest version)
- AiiDA-WorkGraph (latest version)
- pymatgen
- numpy
- ase
- VASP with the AiiDA-VASP plugin (currently the only supported DFT code)

### Setup

```bash
# Clone the repository
git clone git@github.com:DoriniTT/PS-TEROS.git
cd PS-TEROS

# Install in development mode
pip install -e .

# Verify AiiDA profile is configured
verdi profile list
verdi status

# Start the AiiDA daemon if needed
verdi daemon start
```

**Note**: After any code modifications, restart the daemon:
```bash
verdi daemon restart
```

## Quick Start

### Basic Workflow with Presets (Recommended)

The simplest way to run PS-TEROS is using workflow presets:

```python
from teros.core.workgraph import build_core_workgraph
from aiida.engine import submit

# Create workgraph using preset
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # One line activates full workflow!

    # Structure files
    structures_dir="/path/to/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    oxygen_name="O2.cif",

    # Code configuration
    code_label="vasp@localhost",
    potential_family="PBE",
    kpoints_spacing=0.4,

    # Potential mappings
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},

    # VASP parameters (simplified)
    bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520},

    # Computational resources
    bulk_options={'resources': {'num_machines': 1}},
    metal_options={'resources': {'num_machines': 1}},
    oxygen_options={'resources': {'num_machines': 1}},

    # Slab generation
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,

    name="Ag2O_100_surface"
)

# Submit to AiiDA
node = submit(wg)
print(f"Submitted PS-TEROS workflow: PK={node.pk}")
```

**What this does:**
- Relaxes bulk and reference structures
- Calculates formation enthalpy
- Generates and relaxes (100) surface slabs
- Calculates surface energies vs. chemical potential

**Optional features** (disabled by default, add these flags to enable):
```python
compute_cleavage=True,              # Enable cleavage energy calculations
compute_relaxation_energy=True,      # Enable relaxation energy calculations
```

### Workflow with Manual Terminations

For precise control over surface structures:

```python
from aiida.orm import load_node, StructureData
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map

# Load pre-generated slab structures
slabs = {
    'term_0': StructureData(ase=read('/path/to/slab_0.vasp')),
    'term_1': StructureData(ase=read('/path/to/slab_1.vasp')),
    'term_2': StructureData(ase=read('/path/to/slab_2.vasp')),
}

# Create workgraph with manual slabs
wg = build_core_workgraph_with_map(
    **defaults,
    input_slabs=slabs,  # Provide pre-generated slabs
    # Miller indices not required when using input_slabs
    relax_slabs=True,
    compute_thermodynamics=True,
    name="Custom_Terminations"
)

node = submit(wg)
```

### Restart from Previous Calculation

Resume calculations that didn't complete:

```python
# Restart from a previous PS-TEROS run
wg = build_core_workgraph(
    **defaults,
    restart_from_node=12345,  # PK of previous calculation
    relax_slabs=True,
    name="Ag3PO4_restart"
)

node = submit(wg)
```

The restart feature automatically:
- Extracts slab structures from the previous run
- Loads RemoteData folders containing WAVECAR and CONTCAR
- Continues relaxations from the last ionic step
- Saves computational time and improves convergence

### Cleavage Energy Calculation

Calculate cleavage energies for complementary surface pairs:

```python
wg = build_core_workgraph_with_map(
    **defaults,
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    oxygen_name="O2.cif",
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    compute_cleavage=True,  # Enable cleavage calculation
    name="Ag2O_cleavage"
)

node = submit(wg)

# Access results after completion
results = load_node(node.pk).outputs
cleavage_energies = results.cleavage_energies
```

## Three-Tier Workflow System

PS-TEROS uses a powerful three-tier configuration system that balances simplicity with flexibility:

### Tier 1: Named Workflow Presets (Recommended)

Use a single parameter to activate complete workflows:

```python
workflow_preset='surface_thermodynamics'  # Activates full thermodynamics workflow
```

**Available presets:**
- `surface_thermodynamics` (default) - Complete surface energy calculations
- `surface_thermodynamics_unrelaxed` - Quick screening with unrelaxed slabs
- `cleavage_only` - Cleavage energy calculations
- `relaxation_energy_only` - Surface reconstruction analysis
- `bulk_only` - Bulk optimization only
- `formation_enthalpy_only` - Formation enthalpy without surfaces
- `electronic_structure_bulk_only` - DOS/bands for bulk
- `electronic_structure_slabs_only` - DOS/bands for slabs
- `electronic_structure_bulk_and_slabs` - DOS/bands for both
- `aimd_only` - Molecular dynamics simulations
- `comprehensive` - Everything enabled

**List available presets:**
```python
from teros.core import list_workflow_presets
list_workflow_presets()
```

### Tier 2: Individual Component Flags (Fine-Grained Control)

Override preset defaults with specific flags:

```python
workflow_preset='surface_thermodynamics',
compute_cleavage=True,              # Add cleavage energies
compute_relaxation_energy=False,    # Skip relaxation energies
```

**Available flags:**
- `relax_slabs` - Enable/disable slab relaxation
- `compute_thermodynamics` - Enable/disable surface energy calculations
- `compute_cleavage` - Enable/disable cleavage energy calculations
- `compute_relaxation_energy` - Enable/disable relaxation energy calculations
- `compute_electronic_properties_bulk` - Enable/disable bulk DOS/bands
- `compute_electronic_properties_slabs` - Enable/disable slab DOS/bands
- `run_aimd` - Enable/disable AIMD simulations

### Tier 3: Automatic Dependency Resolution

The system automatically validates your configuration and warns about conflicts:

```python
# This will warn you:
workflow_preset='surface_thermodynamics',
relax_slabs=False,  # Breaks dependency!
compute_cleavage=True,  # Needs relaxed slabs
```

**For detailed documentation, see:**
- [Workflow System Explained](docs/WORKFLOW_SYSTEM_EXPLAINED.md) - Understanding the three tiers
- [Workflow Presets Guide](docs/WORKFLOW_PRESETS_GUIDE.md) - Complete preset reference
- [Workflow Presets Examples](docs/WORKFLOW_PRESETS_EXAMPLES.md) - Runnable examples
- [Migration Guide](docs/WORKFLOW_MIGRATION_GUIDE.md) - Updating existing scripts

---

## Workflow Architecture

PS-TEROS workflows follow this execution pattern:

1. **Bulk and Reference Relaxations** (parallel)
   - Bulk oxide structure optimization
   - Elemental references (metal, nonmetal)
   - Oxygen molecule calculation

2. **Formation Enthalpy Calculation**
   - Compute bulk formation energy from relaxed structures

3. **Slab Generation or Input**
   - Automatic: Generate terminations using pymatgen SlabGenerator
   - Manual: Use provided StructureData objects

4. **Slab Relaxation** (parallel, optional)
   - Relax all surface terminations
   - Optional restart from previous RemoteData

5. **Analysis** (parallel)
   - **Thermodynamics**: Surface Gibbs free energy vs chemical potentials
   - **Cleavage Energy**: Energy cost of cleaving complementary surfaces

## Core Modules

### `teros.core.workgraph`
Main workflow construction functions:
- `build_core_workgraph()`: Primary workflow builder
- `build_core_workgraph_with_map()`: Simplified interface with structure file mapping
- `core_workgraph()`: Low-level workflow graph construction

### `teros.core.slabs`
Slab generation and relaxation:
- `get_slabs()`: Generate surface terminations from bulk structure
- `relax_slabs_scatter()`: Parallel slab relaxation
- `extract_restart_folders_from_node()`: Extract RemoteData for restart

### `teros.core.thermodynamics`
Surface energy calculations:
- Compute surface Gibbs free energy as function of chemical potentials
- Generate stability diagrams for binary and ternary oxides
- Calculate allowed chemical potential ranges

### `teros.core.cleavage`
Cleavage energy calculations:
- `calculate_cleavage_energy()`: Energy for single complementary pair
- `compute_cleavage_energies_scatter()`: Parallel calculation for all pairs
- Automatic pairing of complementary terminations

### `teros.core.builders`
Default parameter sets:
- `get_ag2o_defaults()`: Silver oxide system
- `get_ag3po4_defaults()`: Silver phosphate system
- Easily customizable for other systems

## Key Parameters

### Structure Parameters
- `bulk_name`: Bulk oxide structure file
- `metal_name`: Metal reference structure
- `nonmetal_name`: Nonmetal reference (ternary oxides only)
- `oxygen_name`: O2 molecule structure

### Slab Generation Parameters
- `miller_indices`: Miller indices for surface plane (e.g., [1, 0, 0])
- `min_slab_thickness`: Minimum slab thickness in Angstroms
- `min_vacuum_thickness`: Minimum vacuum gap in Angstroms
- `input_slabs`: Dictionary of pre-generated slabs (optional)

### Calculation Control
- `workflow_preset`: Named workflow preset (e.g., 'surface_thermodynamics') - **Recommended**
- `relax_slabs`: Enable/disable slab relaxation
- `compute_thermodynamics`: Enable/disable surface energy calculations
- `compute_cleavage`: Enable/disable cleavage energy calculations (optional in most presets)
- `compute_relaxation_energy`: Enable/disable relaxation energy calculations (optional in most presets)
- `compute_electronic_properties_bulk`: Enable/disable bulk DOS/bands
- `compute_electronic_properties_slabs`: Enable/disable slab DOS/bands
- `run_aimd`: Enable/disable AIMD simulations
- `restart_from_node`: PK of previous calculation for restart (optional)

**Note:** When using `workflow_preset`, most flags are set automatically. Use individual flags only to override preset defaults.

### DFT Parameters
- `code_label`: AiiDA code label (e.g., "vasp@cluster")
- `potential_family`: Pseudopotential family
- `*_parameters`: INCAR/input parameters for each calculation type
- `*_kpoints_distance`: k-point density
- `*_options`: Computational resources

## Examples

The `examples/` directory contains comprehensive demonstrations:

- **`examples/default_builders/`**: Using default parameter sets
- **`examples/restart/`**: Restart functionality and patterns
- **`examples/cleavage/`**: Cleavage energy calculations
- **`examples/slabs/`**: Manual slab input workflows
- **`examples/vasp/`**: VASP-specific examples

Each example includes detailed documentation in accompanying `.md` files.

## Analyzing Results

```python
from aiida.orm import load_node

# Load completed workflow
node = load_node(YOUR_PK)

# Access outputs
results = node.outputs

# Bulk and reference energies
bulk_energy = results.bulk_energy
metal_energy = results.metal_energy
oxygen_energy = results.oxygen_energy

# Formation enthalpy
formation_enthalpy = results.formation_enthalpy

# Slab structures and energies (if relaxed)
slab_structures = results.slab_structures
slab_energies = results.slab_energies

# Thermodynamics data (if computed)
thermo_data = results.thermodynamics

# Cleavage energies (if computed)
cleavage_data = results.cleavage_energies

# Remote folders (for restart)
remote_folders = results.slab_remote
```

Use `verdi` commands for workflow inspection:
```bash
# Check workflow status
verdi process show <PK>

# View workflow report
verdi process report <PK>

# List all child processes
verdi process list -a -p <PK>

# Get workflow outputs
verdi node show <PK>
```

## Troubleshooting

### Daemon Issues
```bash
# Check daemon status
verdi daemon status

# Restart daemon after code changes
verdi daemon restart

# View daemon logs
verdi daemon logshow
```

### Workflow Debugging
```bash
# Check failed processes
verdi process list -S excepted

# Get detailed error report
verdi process report <PK>

# Inspect specific calculation
verdi calcjob outputcat <CALCJOB_PK>
```

### Common Issues

1. **Slab generation fails**: Check structure file format and miller indices
2. **Relaxation doesn't converge**: Adjust VASP parameters (EDIFF, NELM, etc.) or use restart feature
3. **Memory errors**: Increase resources in `*_options` parameters
4. **Code not found**: Verify AiiDA code setup with `verdi code list`

## Development

### Project Structure
```
PS-TEROS/
├── teros/
│   ├── core/           # Main workflow modules
│   │   ├── workgraph.py
│   │   ├── slabs.py
│   │   ├── thermodynamics.py
│   │   ├── cleavage.py
│   │   └── builders/   # Default parameter sets
│   └── experimental/   # New features under development
├── examples/           # Example workflows
├── docs/              # Documentation
├── CHANGE.md          # Version history
└── README.md          # This file
```

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Update documentation
5. Submit a pull request

## Citation

If you use PS-TEROS in your research, please cite:
```
[Citation information to be added]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For questions and issues:
- GitHub Issues: https://github.com/DoriniTT/PS-TEROS/issues
- Documentation: See `docs/` directory
- Examples: See `examples/` directory

## Acknowledgments

PS-TEROS is built on the AiiDA and AiiDA-WorkGraph frameworks. For more information:
- AiiDA: https://www.aiida.net
- AiiDA-WorkGraph: https://aiida-workgraph.readthedocs.io
