.. _tutorial:

=========================
Your first PS-TEROS graph
=========================

Here you will build one Quantum ESPRESSO graph from a bulk SnO2 structure. The
graph stays unsubmitted, so you can inspect how the pieces fit together without
starting a calculation.

By the end, you will know how to:

* describe the calculation inputs with ``QeCalculationConfig``;
* identify the registered code that selects the AiiDA computer;
* make the queue and resource choices explicit with ``ExecutionPolicy``; and
* build an unsubmitted graph from a labelled structure.

Before you start
----------------

Complete the :doc:`installation steps <installation>` first. To build the graph,
you also need an active AiiDA profile, a registered ``quantumespresso.pw`` code,
and an ``aiida-pseudo`` family. The identifiers in the example below are
placeholders. Replace them with values from your own AiiDA profile.

Build the recipe
----------------

A **recipe** collects the settings that should be shared by a calculation. The
example uses a bulk rutile SnO2 starting structure because it is small and
makes the shape of the API easy to see. The numerical values only make this
example concrete; converge them for the material and property you study.

.. code-block:: python

   import psteros

   qe_parameters = {
       "CONTROL": {"calculation": "relax"},
       "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
       "ELECTRONS": {"conv_thr": 1.0e-8},
       "IONS": {"ion_dynamics": "bfgs"},
   }

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

   recipe = psteros.SurfaceWorkflowConfig(
       backend="qe",
       calculation=psteros.QeCalculationConfig(
           code_label="your-qe-code@your-computer",
           pseudo_family="your-pseudo-family",
           parameters=qe_parameters,
           kpoints_distance=0.20,
       ),
       execution=execution,
       name="sno2_bulk_relax",
   )

   graph = psteros.build_surface_workgraph(
       {"bulk_sno2": psteros.rutile_sno2_bulk()},
       recipe,
   )

   print(graph.name)
   print(graph.max_number_jobs)

You should see:

.. code-block:: text

   sno2_bulk_relax
   1

The first line confirms the graph name from the recipe. The second is the
number of calculations the graph may run at once.

What the code did
-----------------

``QeCalculationConfig`` describes the electronic-structure calculation and
identifies the registered code; that code selects the actual AiiDA computer.
``ExecutionPolicy`` supplies the queue, resource, wall-time, and MPI metadata.
Its ``computer`` field is descriptive in the current API, so keep it consistent
with the computer in ``code_label``. The builder does not cross-check them.
``build_surface_workgraph`` currently keeps one active calculation in a graph,
so ``max_concurrent_jobs`` must be ``1``. The builder then connects
the labelled structure and the recipe into an AiiDA WorkGraph.

The call above leaves ``submit`` at its default value, ``False``. It builds a
graph in your local AiiDA environment but does not request scheduler resources
or start Quantum ESPRESSO.

What you have now
-----------------

You have an inspectable graph for one labelled bulk calculation. Before turning
it into a submitted calculation, check that the code label, pseudopotential
family, numerical settings, queue, and resource request are appropriate for
your project.

The :doc:`QE guide <qe-first-workflow>` shows how to connect a geometry
relaxation to a final static calculation. Read :doc:`how a PS-TEROS calculation
fits together <concepts>` for the roles of the structures, recipe, graph, and
analysis helpers.
