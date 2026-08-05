# PS-TEROS

PS-TEROS is a Python library for organizing surface thermodynamics calculations
with AiiDA. It helps you turn starting structures, a calculation recipe, and
an execution policy into an inspectable calculation graph. When the graph runs,
AiiDA records the inputs, outputs, and calculation steps that produced the
result.

The library provides building blocks; it does not choose a material model,
convergence settings, compute resource, or scientific conclusion for you.
Those choices belong to the project that uses it.

## Start here

- **Install the library:** [installation guide](docs/source/installation.rst)
- **Build one graph without submitting work:** [first tutorial](docs/source/tutorial.rst)
- **Understand the pieces of a calculation:** [core concepts](docs/source/concepts.rst)
- **Prepare a relaxation followed by a static QE calculation:**
  [QE guide](docs/source/qe-first-workflow.rst)
- **Understand the SnO2 surface-energy model:**
  [model explanation](docs/source/examples.rst)
- **Look up public names and inputs:** [reference](docs/source/api.rst)

## A safe first interaction

The structure helpers can be used before you connect to a scheduler or submit a
calculation. This example creates an in-memory starting model for one rutile
SnO2(110) surface; it does not use an AiiDA profile or external compute time.

```python
import psteros

slab, identity = psteros.sno2_110_slab(
    termination="o",
    triple_layers=9,
    vacuum_angstrom=20.0,
)

print(identity.termination)  # o
print(slab.composition.reduced_formula)  # SnO2
```

A starting structure is not yet a converged result. The SnO2 model page explains
what a slab and a termination are, which calculations provide the needed
energies, and which assumptions must be checked before interpreting a surface
energy.

## Installation from a source checkout

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
python -c "import psteros; print(psteros.__version__)"
```

Before building a calculation graph, configure an AiiDA profile, a
`quantumespresso.pw` code, and an `aiida-pseudo` family for your own computing
environment. The [installation guide](docs/source/installation.rst) explains
what each of those pieces contributes and points to AiiDA's setup guides.

The main supported workflow uses Quantum ESPRESSO. The installation guide also
describes the optional VASP adapter for existing studies.
