# PS-TEROS: Predicting Stability of TERminations of Oxide Surfaces

PS-TEROS is a Python package for automated calculation of oxide surface energies using ab initio atomistic thermodynamics. Built on the AiiDA-WorkGraph framework, it provides a complete workflow from bulk relaxation through slab generation to surface energy calculations.

## Features

- **Automated Workflow**: Complete end-to-end workflow for surface energy calculations of ternary oxides
- **Parallel Execution**: Bulk and reference relaxations run in parallel; all slab terminations relax simultaneously
- **Dynamic Slab Generation**: Automatically generates and processes all unique surface terminations
- **Provenance Tracking**: Full AiiDA integration for complete computational provenance
- **Modular Architecture**: Clean separation of concerns with reusable task components
- **VASP Integration**: Native support for VASP calculations via aiida-vasp plugin

## How PS-TEROS Works

### Workflow Architecture

PS-TEROS uses a **graph-based workflow** built with AiiDA-WorkGraph. The workflow consists of several stages that execute in a coordinated manner:

#### 1. Structure Loading & Parallel Relaxations
The workflow begins by loading all input structures and relaxing them in parallel:
- **Bulk oxide** (e.g., Ag₃PO₄) - relaxed with full cell optimization (ISIF=3)
- **Metal reference** (e.g., Ag) - for elemental chemical potential
- **Nonmetal reference** (e.g., P) - for elemental chemical potential
- **Oxygen reference** (O₂) - for oxygen chemical potential

All four relaxations execute **simultaneously** since they have no dependencies.

#### 2. Formation Enthalpy Calculation
Once all relaxations complete, the workflow calculates the formation enthalpy:
```
ΔH_f = E_bulk - (n_metal × E_metal/atom + n_nonmetal × E_nonmetal/atom + n_O × E_O2/atom)
```
This provides the thermodynamic stability of the bulk oxide relative to its elements.

#### 3. Slab Generation
Using the relaxed bulk structure, PS-TEROS automatically generates all unique surface terminations:
- Converts bulk to primitive cell (optional)
- Generates slabs for specified Miller indices (e.g., (100), (110), (111))
- Creates orthogonal supercells with appropriate vacuum spacing
- Identifies all symmetrically distinct terminations

The number of terminations is **dynamically determined** - the workflow adapts to however many unique surfaces exist.

#### 4. Parallel Slab Relaxations (Map Zone)
When enabled (`relax_slabs=True`), PS-TEROS relaxes all generated slabs in parallel using the **Map context manager**:
- Each termination is relaxed independently with VASP
- Atomic positions optimized (ISIF=2), cell kept fixed
- All relaxations run **simultaneously** for maximum efficiency
- Energies extracted automatically for surface energy calculations

### Key Components

**Core Modules** (`teros/CORE/`):
- `workgraph.py`: Main workflow definitions
  - `core_workgraph`: Graph task using `@task.graph` decorator
  - `build_core_workgraph_with_map`: Recommended builder with Map zone for slab relaxations
  - `load_structure`, `extract_total_energy`: Utility tasks

- `modules/slabs.py`: Slab generation using Pymatgen
  - `get_slabs`: Creates all terminations with dynamic namespace output

- `modules/hf.py`: Thermodynamics calculations
  - `calculate_formation_enthalpy`: Formation enthalpy from DFT energies

**Task Types**:
- `@task`: Wraps Python functions as workflow tasks
- `@task.calcfunction`: AiiDA calcfunction for provenance
- `@task.graph`: Builds sub-workflows with Python logic
- `task(WorkChain)`: Wraps AiiDA WorkChains (e.g., VASP) as tasks

**Dynamic Namespaces**:
PS-TEROS uses dynamic namespaces to handle variable numbers of slabs - the workflow automatically adapts whether you have 2 or 20 terminations.

## Installation

```bash
# Clone the repository
git clone git@github.com:DoriniTT/PS-TEROS.git
cd PS-TEROS

# Install in development mode
pip install -e .
```

**Requirements:**
- Python 3.9+
- AiiDA Core (configured with a profile)
- AiiDA-WorkGraph (version 1.0.0b3 or compatible)
- aiida-vasp plugin
- Pymatgen
- ASE (Atomic Simulation Environment)

**Important:** PS-TEROS currently uses AiiDA-WorkGraph 1.0 beta. Install compatible version:
```bash
pip install aiida-workgraph>=1.0.0b3
```

## Quick Start

### Example: Ag₃PO₄ (100) Surface Relaxation

```python
from aiida import load_profile
from teros.workgraph import build_core_workgraph_with_map

# Load AiiDA profile
load_profile()

# Define calculation parameters
structures_dir = "/path/to/structures"
code_label = "VASP-6.4.3@cluster"
potential_family = "PBE"

# VASP parameters for bulk relaxation (full optimization)
bulk_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,  # Relax cell + atoms
    "NSW": 100,
    "EDIFFG": -0.1,
}

bulk_options = {
    "resources": {"num_machines": 1, "num_cores_per_machine": 40},
    "queue_name": "parallel",
}

# VASP parameters for slab relaxation (atoms only)
slab_parameters = {
    **bulk_parameters,
    "ISIF": 2,  # Relax atoms only, fix cell
    "EDIFFG": -0.02,  # Tighter convergence for surfaces
}

# Build the workflow
wg = build_core_workgraph_with_map(
    structures_dir=structures_dir,
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label=code_label,
    potential_family=potential_family,
    bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"P": "P"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters=bulk_parameters,
    bulk_options=bulk_options,
    # Reference parameters (metal, nonmetal, O2)
    metal_parameters={...},
    metal_options={...},
    nonmetal_parameters={...},
    nonmetal_options={...},
    oxygen_parameters={...},
    oxygen_options={...},
    # Slab generation
    miller_indices=[1, 0, 0],  # (100) surface
    min_slab_thickness=10.0,  # Angstroms
    min_vacuum_thickness=15.0,  # Angstroms
    symmetrize=True,
    # Slab relaxation
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options=bulk_options,
    name="Ag3PO4_100_Surface",
)

# Submit the workflow
wg.submit(wait=False)
print(f"Submitted WorkGraph: {wg.pk}")
```

### Monitoring

```bash
# Check workflow status
verdi process show <PK>

# Monitor all processes
verdi process list -a -p1

# View workflow report
verdi process report <PK>
```

### Accessing Results

```python
from aiida import load_node

# Load completed workflow
wg = load_node(<PK>)

# Formation enthalpy
hf = wg.outputs.formation_enthalpy.get_dict()
print(f"ΔH_f = {hf['formation_enthalpy_ev']:.3f} eV/formula unit")

# Generated slabs (unrelaxed)
slab_structures = wg.outputs.slab_structures
print(f"Generated {len(slab_structures)} terminations")

# Relaxed slabs
relaxed_slabs = wg.outputs.relaxed_slabs
for term_id, structure in relaxed_slabs.items():
    energy = wg.outputs.slab_energies[term_id].value
    print(f"{term_id}: E = {energy:.3f} eV")

    # Export to file
    atoms = structure.get_ase()
    atoms.write(f"relaxed_{term_id}.cif")
```

For complete examples with all parameters, see `teros/examples/slabs/slabs_relax.py`.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Workflow Visualization

PS-TEROS workflows can be exported to interactive HTML visualizations:

```python
# After building the workflow
wg.to_html("my_workflow.html")
```

This creates a visual graph showing:
- All tasks and their connections
- Data flow between tasks
- Parallel execution paths
- Task status (when viewing a running/completed workflow)

## Advanced Features

### Custom VASP Parameters per Component

Each component (bulk, metal, nonmetal, oxygen, slabs) can have independent VASP parameters:

```python
# Different ISMEAR for metals vs. insulators
metal_parameters = {
    "ISMEAR": 1,  # Methfessel-Paxton for metals
    "SIGMA": 0.2,
}

bulk_parameters = {
    "ISMEAR": 0,  # Gaussian smearing for insulators
    "SIGMA": 0.05,
}
```

### Slab Generation Options

Fine-tune slab generation behavior:

```python
wg = build_core_workgraph_with_map(
    # ... other parameters ...
    miller_indices=[1, 1, 0],  # Different surface orientation
    min_slab_thickness=15.0,  # Thicker slabs
    min_vacuum_thickness=20.0,  # More vacuum
    symmetrize=True,  # Only symmetrically distinct terminations
    primitive=True,  # Use primitive cell (smaller slabs)
    lll_reduce=True,  # LLL-reduced cell vectors
    center_slab=True,  # Center slab in cell
)
```

### Two Workflow Implementations

PS-TEROS provides two equivalent workflow builders:

1. **`build_core_workgraph_with_map`** (Recommended)
   - Uses Map context manager for slab relaxations
   - Cleaner provenance graph
   - Better visualization
   - Full AiiDA-WorkGraph 1.0 features

2. **`build_core_workgraph`**
   - Uses `@task.graph` decorator with two-pass approach
   - Compatible with older AiiDA-WorkGraph versions
   - Functionally equivalent results

Both produce identical scientific results. Use `build_core_workgraph_with_map` for new workflows.

## Project Structure

```
PS-TEROS/
├── teros/
│   ├── CORE/
│   │   ├── workgraph.py          # Main workflow logic
│   │   └── modules/
│   │       ├── slabs.py          # Slab generation
│   │       ├── hf.py             # Formation enthalpy
│   │       └── helper_functions.py
│   ├── examples/
│   │   └── slabs/
│   │       └── slabs_relax.py    # Complete example script
│   ├── structures/               # Example structure files
│   └── workgraph.py              # Re-exports for compatibility
├── README.md
└── setup.py
```

## Documentation

Comprehensive documentation is available in the `teros/CLAUDE.md` file, including:
- AiiDA-WorkGraph concepts (Tasks, Sockets, Namespaces)
- Graph tasks vs. Context managers
- Map zone implementation details
- Quick reference commands

For AiiDA and AiiDA-WorkGraph documentation:
- [AiiDA Core Documentation](https://aiida.readthedocs.io/)
- [AiiDA-WorkGraph Documentation](https://aiida-workgraph.readthedocs.io/)
- [Materials Science Tutorial](https://aiida-workgraph.readthedocs.io/en/latest/tutorial/autogen/materials_science_ase.html)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use PS-TEROS in your research, please cite:
```
[Citation information will be added upon publication]
```
