# PS-TEROS: Predicting Stability of TERminations of Oxide Surfaces

PS-TEROS is a Python package for calculating surface Gibbs free energy using ab initio atomistic thermodynamics with the AiiDA-WorkGraph framework. It automates the calculation of surface energies for binary and ternary oxides under varying thermodynamic conditions.

## Features

- **Automated Surface Energy Calculations**: Complete workflow for surface Gibbs free energy calculations using ab initio atomistic thermodynamics
- **Multiple DFT Code Support**: Compatible with VASP, Quantum ESPRESSO, and CP2K via AiiDA plugins
- **Flexible Slab Generation**: Automatic generation or manual input of surface terminations
- **Cleavage Energy Calculations**: Built-in module for computing cleavage energies of complementary surface pairs
- **Restart Capability**: Resume incomplete calculations from previous runs using remote data
- **Default Builders**: Pre-configured parameter sets for rapid workflow setup
- **Full Provenance Tracking**: Complete data provenance via AiiDA infrastructure
- **Parallel Execution**: Efficient parallel processing of multiple terminations and references

## Installation

### Requirements

- Python 3.9+
- AiiDA Core (latest version)
- AiiDA-WorkGraph (latest version)
- pymatgen
- numpy
- ase
- A supported DFT code (VASP, Quantum ESPRESSO, or CP2K) with corresponding AiiDA plugin

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

### Basic Workflow with Default Builders

The simplest way to run PS-TEROS is using default builders:

```python
from teros.core.builders import get_ag3po4_defaults
from teros.core.workgraph import build_core_workgraph_with_map
from aiida.engine import submit

# Get default parameters for Ag3PO4 system
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="vasp@localhost",
    potential_family="PBE.54",
)

# Create and submit workgraph
wg = build_core_workgraph_with_map(
    **defaults,
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    compute_thermodynamics=True,
    name="Ag3PO4_100_surface"
)

# Submit to AiiDA
node = submit(wg)
print(f"Submitted PS-TEROS workflow: PK={node.pk}")
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
- `relax_slabs`: Enable/disable slab relaxation (default: False)
- `compute_thermodynamics`: Enable surface energy calculations (default: False)
- `compute_cleavage`: Enable cleavage energy calculations (default: False)
- `restart_from_node`: PK of previous calculation for restart (optional)

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
