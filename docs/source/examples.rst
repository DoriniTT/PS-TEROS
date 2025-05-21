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
    
    # --- Placeholder: Setup VASP builders ---
    # This is a simplified illustration. Refer to full examples for details.
    
    # Structure for bulk Ag2O (example)
    # bulk_structure = StructureData() # Replace with actual structure loading/definition
    # ParameterData = DataFactory('core.dict') # AiiDA's Dict class

    # Builder for bulk calculation
    # builder_bulk = RelaxWorkChain.get_builder()
    # builder_bulk.structure = bulk_structure
    # builder_bulk.relax.algo = "cg" 
    # builder_bulk.kpoints_distance = 0.2
    # Some VASP specific settings might be passed via builder_bulk.vasp.parameters = ParameterData(dict={...})
    # ... other VASP specific settings for bulk ...

    # Structure for Ag2O slab (example)
    # slab_structure = StructureData() # Replace with actual slab structure

    # Builder for slab calculation
    # builder_slab = RelaxWorkChain.get_builder()
    # builder_slab.structure = slab_structure
    # ... similar settings as builder_bulk, adjusted for slabs ...
    
    # Reference calculations (e.g., O2, Ag)
    # builder_o2 = RelaxWorkChain.get_builder()
    # ... setup for O2 calculation ...
    # builder_ag = RelaxWorkChain.get_builder()
    # ... setup for Ag bulk calculation ...
    # references = {'O': builder_o2, 'Ag': builder_ag}
    # --- End Placeholder ---
    
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
    
    # --- Placeholder: Setup QE builders ---
    # This is a simplified illustration. Refer to full examples for details.
    
    # Structure for bulk (example)
    # bulk_structure = StructureData() # Replace with actual structure

    # Builder for bulk calculation
    # builder_bulk = PwRelaxWorkChain.get_builder()
    # builder_bulk.structure = bulk_structure
    # builder_bulk.base.pw.parameters = {'CONTROL': {'calculation': 'relax'}} # Example parameters
    # builder_bulk.base.kpoints_distance = 0.2
    # ... other QE specific settings ...

    # Structure for slab (example)
    # slab_structure = StructureData() # Replace with actual slab structure

    # Builder for slab calculation
    # builder_slab = PwRelaxWorkChain.get_builder()
    # builder_slab.structure = slab_structure
    # ... similar settings, adjusted for slabs ...

    # Reference calculations
    # references = {'O': builder_o2_qe, 'Elem': builder_elem_qe} 
    # ... setup for reference builders ...
    # --- End Placeholder ---
    
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
    
    # --- Placeholder: Setup CP2K builders ---
    # This is a simplified illustration. Refer to full examples for details.

    # Structure for bulk (example)
    # bulk_structure = StructureData() # Replace with actual structure
    
    # Builder for bulk calculation
    # builder_bulk = Cp2kRelaxWorkChain.get_builder()
    # builder_bulk.structure = bulk_structure
    # builder_bulk.cp2k.parameters = { 'FORCE_EVAL': { 'DFT': { } } } # Example parameters
    # ... other CP2K specific settings ...

    # Structure for slab (example)
    # slab_structure = StructureData() # Replace with actual slab structure

    # Builder for slab calculation
    # builder_slab = Cp2kRelaxWorkChain.get_builder()
    # builder_slab.structure = slab_structure
    # ... similar settings, adjusted for slabs ...

    # Reference calculations
    # references = {'O': builder_o2_cp2k, 'Elem': builder_elem_cp2k}
    # ... setup for reference builders ...
    # --- End Placeholder ---
    
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
