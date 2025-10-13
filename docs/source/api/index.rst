==============
API Reference
==============

Detailed API documentation for PS-TEROS modules. This section is auto-generated from Python docstrings and provides complete parameter specifications.

Core Modules
============

**Main Workflow Builder**

* :doc:`workgraph` - ``build_core_workgraph()`` and workflow construction

**Feature Modules**

* :doc:`slabs` - Slab generation, relaxation, and energy calculations
* :doc:`thermodynamics` - Surface energy and chemical potential calculations
* :doc:`cleavage` - Cleavage energy for complementary terminations
* :doc:`aimd` - Ab initio molecular dynamics on slabs

**Utilities**

* :doc:`builders` - Pre-configured parameter builders (electronic properties, defaults)

Module Documentation
====================

.. toctree::
   :maxdepth: 2

   workgraph
   slabs
   thermodynamics
   cleavage
   aimd
   builders

Quick Reference: Parameter Inheritance
=======================================

PS-TEROS uses consistent parameter naming across modules:

**Structure Parameters**

* ``structures_dir``: Directory containing input structures
* ``bulk_name``, ``metal_name``, ``nonmetal_name``, ``oxygen_name``: Structure filenames

**Calculation Parameters**

* ``*_parameters``: VASP INCAR dictionaries (e.g., ``bulk_parameters``, ``slab_parameters``)
* ``*_options``: Scheduler options (resources, queue, time limits)
* ``*_potential_mapping``: Element to pseudopotential mapping

**Feature Flags**

* ``relax_slabs``: Enable slab relaxation (default: True)
* ``compute_relaxation_energy``: Calculate E_relax - E_unrelaxed (default: False)
* ``compute_cleavage``: Calculate cleavage energies (default: False)
* ``compute_thermodynamics``: Calculate surface energies (default: False)
* ``compute_electronic_properties_bulk``: DOS/bands for bulk (default: False)
* ``compute_electronic_properties_slabs``: DOS/bands for slabs (default: False)

Common Patterns
===============

**Boolean Flags**: Most features are controlled by boolean flags that default to False (opt-in)

**Parameter Dictionaries**: INCAR parameters are passed as Python dicts, converted internally to AiiDA Dict nodes

**Scatter-Gather**: All parallel operations use scatter-gather pattern for consistency

**Restart Support**: Most calculations accept ``restart_folder`` or ``restart_from_node`` for continuation

Navigation Tips
===============

**Looking for a specific function?** Use the :ref:`genindex` or search box

**Want examples?** See :doc:`Workflows </workflows/index>` for complete working code

**Need quick solutions?** Check :doc:`How-To Guides </how-to/index>` for common tasks
