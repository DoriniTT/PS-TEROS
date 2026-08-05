.. _tutorial:

========
Tutorial
========

This tutorial builds the supported psteros workflow: a QE calculation through
``aiida-quantumespresso``.  It creates a graph but does not submit it, so it
is safe to run while learning the API.

Prerequisites
-------------

Create an AiiDA profile, register a ``quantumespresso.pw`` code and install
an ``aiida-pseudo`` family.  The maintained Bohr route uses
``QE-7.6-PW-GPU-A100@bohr`` and ``SSSP/1.3/PBE/precision``.  Use the pinned
environment described in :doc:`installation`.

Build a QE graph
----------------

.. code-block:: python

   import psteros

   parameters = {
       "CONTROL": {"calculation": "scf", "tstress": True, "tprnfor": True},
       "SYSTEM": {"ecutwfc": 75.0, "ecutrho": 600.0, "occupations": "fixed"},
       "ELECTRONS": {"conv_thr": 1.0e-10},
   }
   recipe = psteros.SurfaceWorkflowConfig(
       backend="qe",
       calculation=psteros.QeCalculationConfig(
           code_label="QE-7.6-PW-GPU-A100@bohr",
           pseudo_family="SSSP/1.3/PBE/precision",
           parameters=parameters,
           kpoints_distance=0.20,
       ),
       execution=psteros.ExecutionPolicy(),
       name="sno2_bulk_scf",
   )
   graph = psteros.build_surface_workgraph(
       {"bulk_sno2": psteros.rutile_sno2_bulk()}, recipe
   )

``ExecutionPolicy`` enforces one active calculation and encodes the Bohr
``gpu_a100`` queue.  Submit only from a project controller which owns the
global queue policy:

.. code-block:: python

   node = graph.submit()
   print(node.pk)

For a surface calculation, use ``psteros.sno2_110_slab`` to make an O, SnO,
or Sn2O symmetric (110) starting model.  Relaxation followed by a final SCF
is constructed with ``build_qe_relax_static_workgraph``.  See
:doc:`qe-first-workflow` for the full method and analysis convention.
