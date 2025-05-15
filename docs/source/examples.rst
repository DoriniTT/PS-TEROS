.. _examples:

========
Examples
========

This section provides examples of how to use TEROS with different DFT codes and for different materials.

VASP Example
-----------

.. code-block:: python

    """
    Example for running a TEROS workflow with VASP
    """
    from teros import create_teros_workgraph
    from aiida.engine import submit
    
    # Import VASP-specific components
    from aiida_vasp.workflows.relax import RelaxWorkChain
    
    # ... setup VASP builders ...
    
    # Create and submit the workflow
    wg = create_teros_workgraph(
        dft_workchain=RelaxWorkChain,
        builder_bulk=builder_bulk,
        builder_slab=builder_slab,
        reference_builders=references,
        code="VASP"
    )
    
    node = submit(wg)

Quantum ESPRESSO Example
----------------------

.. code-block:: python

    """
    Example for running a TEROS workflow with Quantum ESPRESSO
    """
    from teros import create_teros_workgraph
    from aiida.engine import submit
    
    # Import QE-specific components
    from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain
    
    # ... setup QE builders ...
    
    # Create and submit the workflow
    wg = create_teros_workgraph(
        dft_workchain=PwRelaxWorkChain,
        builder_bulk=builder_bulk,
        builder_slab=builder_slab,
        reference_builders=references,
        code="QE"
    )
    
    node = submit(wg)

CP2K Example
-----------

.. code-block:: python

    """
    Example for running a TEROS workflow with CP2K
    """
    from teros import create_teros_workgraph
    from aiida.engine import submit
    
    # Import CP2K-specific components
    from aiida_cp2k.workflows import Cp2kRelaxWorkChain
    
    # ... setup CP2K builders ...
    
    # Create and submit the workflow
    wg = create_teros_workgraph(
        dft_workchain=Cp2kRelaxWorkChain,
        builder_bulk=builder_bulk,
        builder_slab=builder_slab,
        reference_builders=references,
        code="CP2K"
    )
    
    node = submit(wg)

Full examples are available in the `teros/examples/` directory, including:

* `run_workflow_vasp.py`
* `run_workflow_qe.py`
* `run_workflow_cp2k.py`
