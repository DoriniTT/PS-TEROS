.. PS-TEROS documentation master file

===============================================
PS-TEROS: Predicting Stability of TERminations of Oxide Surfaces
===============================================

PS-TEROS (Predicting Stability of TERminations of Oxide Surfaces) is a Python package for calculating surface energies and electronic properties of oxide materials using ab initio atomistic thermodynamics with the AiiDA-WorkGraph framework.

Key Features
============

* **Automated workflows**: End-to-end calculations from structure input to surface thermodynamics
* **Comprehensive analysis**: Surface energies, cleavage energies, formation enthalpies, DOS, band structures
* **AIMD support**: Ab initio molecular dynamics on slab structures with sequential staging
* **Material-agnostic**: Works for binary, ternary oxides, and more complex systems
* **AiiDA integration**: Full provenance tracking, workflow management, and restart capabilities
* **Parallel execution**: Scatter-gather patterns for efficient computation

.. important::
   **DFT Code Support**: The current version only supports VASP. CP2K and Quantum ESPRESSO are planned as future updates. If you need PS-TEROS code for these versions, please use the legacy version in the ``legacy/`` directory.

Quick Links
===========

* **New to PS-TEROS?** Start with the :doc:`Beginner Workflow <workflows/beginner-surface-energy>`
* **Need specific features?** Check the :doc:`How-To Guides <how-to/index>`
* **API reference?** See the :doc:`API Documentation <api/index>`
* **What's new?** Read the :doc:`Changelog <changelog>`

Documentation Structure
=======================

This documentation is organized to help computational materials scientists get productive quickly:

1. **Workflows**: Complete examples from simple to advanced, organized by complexity
2. **How-To Guides**: Task-focused solutions for specific problems
3. **API Reference**: Detailed function and parameter documentation
4. **Theory**: Mathematical background and implementation details

Getting Started
===============

.. toctree::
   :maxdepth: 2
   :caption: Installation & Setup

   installation

.. toctree::
   :maxdepth: 2
   :caption: Workflows

   workflows/index
   workflows/beginner-surface-energy
   workflows/intermediate-with-features
   workflows/advanced-complete
   workflows/aimd-molecular-dynamics

.. toctree::
   :maxdepth: 2
   :caption: How-To Guides

   how-to/index
   how-to/restart-calculations
   how-to/custom-slabs
   how-to/aimd-stages
   how-to/electronic-properties

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index
   api/workgraph
   api/slabs
   api/thermodynamics
   api/cleavage
   api/aimd
   api/builders

.. toctree::
   :maxdepth: 2
   :caption: Additional Information

   theory
   changelog
   contributing
   authors

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
