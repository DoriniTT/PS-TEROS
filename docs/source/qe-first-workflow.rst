.. _qe-first-workflow:

========================================
Prepare a QE relaxation-to-static graph
========================================

Use this guide after the :doc:`first tutorial <tutorial>` when you want a final
static Quantum ESPRESSO calculation to use the geometry from a preceding
relaxation. The guide builds the graph only; it does not submit work.

Why use two calculation stages?
-------------------------------

A relaxation changes atomic positions until the chosen force criterion is met.
A subsequent static self-consistent-field (SCF) calculation evaluates the final
energy on that relaxed geometry. Keeping the stages separate makes the handoff
visible in the AiiDA provenance record and avoids rebuilding an intermediate
structure in client code.

Before you start
----------------

You need an active AiiDA profile, a registered ``quantumespresso.pw`` code, and
an ``aiida-pseudo`` family. The labels and scheduler settings below are
placeholders. Replace them with values from your own environment, and validate
the numerical settings for your material before submitting a calculation.

Choose one execution policy
---------------------------

Both stages must use the same execution policy. The builder turns its queue,
resource, wall-time, and MPI choices into AiiDA task metadata. The registered
code in each recipe selects the actual AiiDA computer. ``ExecutionPolicy`` also
has a descriptive ``computer`` field; keep it consistent with ``code_label``
because the builder does not cross-check them.

.. code-block:: python

   import psteros

   execution = psteros.ExecutionPolicy(
       computer="your-computer",
       queue="your-scheduler-queue",
       resources={
           "num_machines": 1,
           "num_mpiprocs_per_machine": 1,
       },
       max_wallclock_seconds=86_400,
       with_mpi=True,
       max_concurrent_jobs=1,
   )

``build_qe_relax_static_workgraph`` currently permits one active calculation in
a graph, so ``max_concurrent_jobs`` must remain ``1``. That limit belongs to
not a recommendation about how many jobs your cluster can run.

.. warning::

   ``ExecutionPolicy`` currently emits PBS-style custom scheduler directives.
   If your AiiDA computer uses another scheduler, adapt and validate that backend
   before enabling submission; changing the placeholder queue name is not
   sufficient.

Define the two recipes
----------------------

The helper below keeps the shared code, pseudopotential family, k-point
spacing, and execution policy in one place. The relaxation-specific force
criterion belongs in ``CONTROL``; psteros validates that placement before a job
is created.

.. code-block:: python

   relax_parameters = {
       "CONTROL": {
           "calculation": "relax",
           "forc_conv_thr": 1.0e-3,
           "tstress": True,
           "tprnfor": True,
       },
       "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
       "ELECTRONS": {"conv_thr": 1.0e-8},
       "IONS": {"ion_dynamics": "bfgs"},
   }
   static_parameters = {
       "CONTROL": {"calculation": "scf", "tstress": True, "tprnfor": True},
       "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
       "ELECTRONS": {"conv_thr": 1.0e-8},
   }

   def qe_recipe(name, parameters):
       return psteros.SurfaceWorkflowConfig(
           backend="qe",
           calculation=psteros.QeCalculationConfig(
               code_label="your-qe-code@your-computer",
               pseudo_family="your-pseudo-family",
               parameters=parameters,
               kpoints_distance=0.20,
           ),
           execution=execution,
           name=name,
       )

   relax = qe_recipe("sno2_110_relax", relax_parameters)
   static = qe_recipe("sno2_110_static", static_parameters)

Connect the stages
------------------

Pass a labelled starting structure and the two recipes to the dedicated
builder. The static task receives the relaxed ``StructureData`` output from the
first task.

.. code-block:: python

   slab, _ = psteros.sno2_110_slab(
       termination="o",
       triple_layers=9,
       vacuum_angstrom=20.0,
   )

   graph = psteros.build_qe_relax_static_workgraph(
       {"sno2_110_o": slab},
       relax,
       static,
   )

As in the tutorial, this creates an unsubmitted graph because ``submit`` was
not set to ``True``. Inspect the code label, pseudopotential family, input
parameters, and scheduler metadata before you request compute time.

Before submission
-----------------

Check the items that the library cannot choose for you:

* the starting structure and any atoms you intend to constrain;
* convergence of cutoffs, k-point sampling, vacuum, and slab thickness;
* compatibility of the pseudopotentials and numerical settings across every
  calculation you plan to compare; and
* the queue, resources, wall-time, and MPI settings accepted by your AiiDA
  computer and scheduler.

The :doc:`SnO2 surface-energy model <examples>` explains how compatible bulk,
slab, and oxygen-reference energies are combined after the calculations finish.
