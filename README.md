# PS-TEROS

PS-TEROS helps surface-science projects prepare and track related
electronic-structure calculations in AiiDA. You still choose and validate the
physical model, numerical settings, computing resources, and scientific
interpretation.

The current API retains legacy deployment-specific execution defaults for
compatibility. New calculations should pass their own execution settings
explicitly; those defaults are not portable recommendations.

## Start here

- **Install the library:** [installation guide](docs/source/installation.rst)
- **Build one graph without submitting work:** [first tutorial](docs/source/tutorial.rst)
- **Understand the pieces of a calculation:** [core concepts](docs/source/concepts.rst)
- **Understand the SnO2 surface-energy model:**
  [model explanation](docs/source/examples.rst)
- **Prepare a relaxation followed by a static QE calculation:**
  [QE guide](docs/source/qe-first-workflow.rst)
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

## Install

Follow the [installation guide](docs/source/installation.rst) for the source
checkout, a clean virtual environment, package verification, and the AiiDA
computer, code, and pseudopotential setup needed before the tutorial.
