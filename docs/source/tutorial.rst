.. _tutorial:

========
Tutorial
========

Basic Usage
-----------

This tutorial will guide you through setting up and running a basic TEROS workflow.

Setting up AiiDA
~~~~~~~~~~~~~~~~

Before using TEROS, make sure you have AiiDA properly set up. 
AiiDA's `verdi presto` command is recommended for quickly setting up a new profile with necessary services.

1. Start and initialize your AiiDA profile if you haven't already:

   .. code-block:: bash

       verdi presto start

   This command will guide you through creating a new profile (if one doesn't exist) or start an existing one, along with its associated database and message broker services.

2. Configure your computer(s) and code(s) in AiiDA:

   .. code-block:: bash

       verdi computer setup
       verdi computer configure
       verdi code create

Running a Simple Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    
    # Example: Setting up reference calculation builders for VASP
    # Oxygen reference (assuming O2 molecule calculation)
    builder_o2_relax = RelaxWorkChain.get_builder()
    # ... configure builder_o2_relax for an O2 molecule ...
    # Metal reference (assuming Ag bulk calculation)
    builder_ag_bulk_relax = RelaxWorkChain.get_builder()
    # ... configure builder_ag_bulk_relax for Ag bulk ...

    references = {
        'O': builder_o2_relax,  # Builder for 1/2 O2 energy
        'Ag': builder_ag_bulk_relax, # Builder for Ag bulk energy
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
~~~~~~~~~~~~~~~~~~

Once your workflow has completed, you can analyze the results:

.. code-block:: python

    from aiida.orm import load_node
    
    # Load the completed workflow node
    wg_node = load_node(pk)  # Replace with your workflow PK
    
    # Access the results
    surface_energies = wg_node.outputs.surface_energies
    # surface_energies is a dictionary where keys are surface identifiers (e.g., strings)
    # and values are the calculated surface energies (e.g., floats in J/m²).
    
    # Print or process the results
    for surface, energy in surface_energies.items():
        print(f"Surface: {surface}, Energy: {energy} J/m²")

Visualization
~~~~~~~~~~~~~

Utilities for visualizing results, such as plotting surface energy diagrams, are currently under development and will be available in an upcoming release. These tools will help in analyzing the output of TEROS workflows.

Stay tuned for updates on this feature!
