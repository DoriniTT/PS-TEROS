# PS-TEROS: Predicting Stability of TERminations of Oxide Surfaces

PS-TEROS is a Python package for calculating surface Gibbs free energy using ab initio atomistic thermodynamics with the AiiDA Workgraph framework.

## Features

- Automated calculation of surface energies for binary and ternary oxides
- Support for multiple DFT codes (VASP, Quantum ESPRESSO, CP2K)
- Modular design with customizable workflow components
- Integration with AiiDA for provenance tracking and workflow management

## Installation

```bash
# Clone the repository
git clone git@github.com:DoriniTT/PS-TEROS.git
cd PS-TEROS

# Install the package in development mode
pip install -e .

**Requirements:**
- Python 3.9+
- AiiDA Core (ensure it's configured)
- AiiDA WorkGraph
- pymatgen
- numpy
(and a supported DFT code + AiiDA plugin for actual calculations)
```

## Usage

### Basic Workflow

To create and submit a PS-TEROS workflow:

```python
from teros import create_teros_workgraph
from aiida.engine import submit # Changed import
# from aiida.orm import StructureData # Example for creating structures
# from aiida_vasp.workflows.relax import RelaxWorkChain # Example DFT workchain

# 1. Set up your DFT builders (code-specific)
# Example: builder_bulk for the bulk material's relaxation
# builder_bulk = RelaxWorkChain.get_builder()
# builder_bulk.structure = your_bulk_structure_data
# ... (configure other settings: k-points, pseudopotentials, parameters, etc.)

# Example: builder_slab for the slab model's relaxation
# builder_slab = RelaxWorkChain.get_builder()
# builder_slab.structure = your_slab_structure_data
# ... (configure settings for slab)

# Example: reference_builders for elemental/gas-phase references
# builder_o2 = RelaxWorkChain.get_builder() 
# ... (configure for O2 calculation)
# builder_ag = RelaxWorkChain.get_builder()
# ... (configure for Ag bulk calculation)
# references = {'O': builder_o2, 'Ag': builder_ag} 
    
# Placeholder variables - replace with your actual configured builders
dft_workchain_class = None # e.g., RelaxWorkChain
builder_bulk = None
builder_slab = None
references = {}

# 2. Create the workgraph
wg = create_teros_workgraph(
    dft_workchain=dft_workchain_class, # Pass the actual workchain class
    builder_bulk=builder_bulk,
    builder_slab=builder_slab,
    reference_builders=references,
    code="VASP" # Specify your DFT code: "VASP", "QE", or "CP2K"
)

# 3. Submit the workgraph
node = submit(wg) # Changed submission call
print(f"Submitted Teros workgraph with PK: {node.pk}")
```

For more detailed examples, check the `examples` directory.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Documentation

The documentation is available online at [Read the Docs](https://aiida-teros.readthedocs.io/).

To build the documentation locally:

```bash
# Install documentation dependencies
pip install -e .[docs]

# Build the documentation
cd docs
make html

# Open the documentation in your browser
open build/html/index.html  # On Linux use: xdg-open build/html/index.html
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
