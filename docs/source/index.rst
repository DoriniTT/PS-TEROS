======================
PS-TEROS documentation
======================

PS-TEROS helps you organize surface thermodynamics calculations with AiiDA. The
Python package is imported as ``psteros``. It keeps the starting structures,
calculation inputs, and execution choices visible as you build a graph. It does
not replace the scientific checks needed to choose a model or interpret a
result.

The pages below follow the questions readers usually ask first: how to install
the library, how to build a safe first graph, how the pieces fit together, and
where to find the supported public API.

Choose the path that matches what you need now:

* **New to PS-TEROS?** :doc:`Install the package <installation>`, then
  :doc:`build one unsubmitted graph <tutorial>`.
* **Learning the model?** Read :doc:`how a calculation fits together <concepts>`
  before the :doc:`SnO2 surface-energy model <examples>`.
* **Preparing QE work?** Use the :doc:`relaxation-to-static guide
  <qe-first-workflow>` after the tutorial.
* **Looking up an input?** Open the :doc:`API reference <api>`.

.. toctree::
   :maxdepth: 1
   :caption: Start here

   installation
   tutorial

.. toctree::
   :maxdepth: 1
   :caption: Learn the model

   concepts
   examples

.. toctree::
   :maxdepth: 1
   :caption: Prepare a calculation

   qe-first-workflow

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api
