# TEROS: Thermodynamics of Oxide Surfaces

Teros is a Python package for calculating surface Gibbs free energy using ab initio atomistic thermodynamics with the AiiDA Workgraph framework.

## Features

- Automated calculation of surface energies for binary and ternary oxides
- Support for multiple DFT codes (VASP, Quantum ESPRESSO, CP2K)
- Modular design with customizable workflow components
- Integration with AiiDA for provenance tracking and workflow management

## Installation

```bash
# Clone the repository
git clone git@github.com:DoriniTT/aiida-teros.git
cd aiida-teros

# Install the package in development mode
pip install -e .
```

## Usage

### Basic Workflow

To create and submit a Teros workflow:

```python
from teros import create_teros_workgraph
from aiida.orm import load_node

# Set up your DFT builders
# ...

# Create the workgraph
wg = create_teros_workgraph(
    dft_workchain=MyDFTWorkChain,
    builder_bulk=builder_bulk,
    builder_slab=builder_slab,
    reference_builders=references,
    code="VASP"
)

# Submit the workgraph
node = wg.submit()
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
