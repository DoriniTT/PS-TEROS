.. _tutorial:

========
Tutorial
========

Basic Usage
----------

This tutorial will guide you through setting up and running a basic TEROS workflow.

Setting up AiiDA
~~~~~~~~~~~~~~~

Before using TEROS, make sure you have AiiDA properly set up:

1. Initialize your AiiDA profile if you haven't already:

   .. code-block:: bash

       verdi quicksetup

2. Configure your computer and codes in AiiDA:

   .. code-block:: bash

       verdi computer setup
       verdi computer configure
       verdi code create

Running a Simple Workflow
~~~~~~~~~~~~~~~~~~~~~~~

Here's a simple example of how to set up and run a TEROS workflow with VASP:

.. code-block:: python

    from teros import create_teros_workgraph
    from aiida.engine import submit
    from aiida.orm import load_node

    # Prepare your DFT builders
    # This is code-specific, here is a placeholder example
    from aiida_vasp.workflows.relax import RelaxWorkChain
    
    # Set up bulk calculation builder
    builder_bulk = RelaxWorkChain.get_builder()
    # ... configure builder_bulk settings ...
    
    # Set up slab calculation builder
    builder_slab = RelaxWorkChain.get_builder()
    # ... configure builder_slab settings ...
    
    # Set up reference calculations (e.g., O2, pure metals)
    references = {
        'O': builder_o2,
        'Ag': builder_ag,
        # Add more references as needed
    }
    
    # Create the workgraph
    wg = create_teros_workgraph(
        dft_workchain=RelaxWorkChain,
        builder_bulk=builder_bulk,
        builder_slab=builder_slab,
        reference_builders=references,
        code="VASP"
    )
    
    # Submit the workgraph
    node = submit(wg)
    print(f"Submitted workgraph with pk: {node.pk}")

Analyzing Results
~~~~~~~~~~~~~~~

Once your workflow has completed, you can analyze the results:

.. code-block:: python

    from aiida.orm import load_node
    
    # Load the completed workflow node
    wg_node = load_node(pk)  # Replace with your workflow PK
    
    # Access the results
    surface_energies = wg_node.outputs.surface_energies
    
    # Print or process the results
    for surface, energy in surface_energies.items():
        print(f"Surface: {surface}, Energy: {energy} J/m²")

Visualization
~~~~~~~~~~~~

TEROS also provides utilities for visualizing your results in the `teros.utils.plots` module:

.. code-block:: python

    from teros.utils.plots import plot_surface_energy_vs_potential
    
    # Generate a surface energy plot as a function of chemical potential
    figure = plot_surface_energy_vs_potential(wg_node)
    figure.savefig('surface_energy_plot.png')
