.. _tutorial:

=========================
Your first PS-TEROS graph
=========================

Here you will build one Quantum ESPRESSO graph from a bulk SnO2 structure. The
graph stays unsubmitted, so you can inspect how the pieces fit together without
starting a calculation.

By the end, you will know how to:

* separate the calculation inputs from the scheduler settings;
* combine them into a recipe for one labelled structure; and
* confirm that the resulting graph remains unsubmitted.

Before you start
----------------

Complete the :doc:`installation steps <installation>` first. To build the graph,
you also need an active AiiDA profile, a registered ``quantumespresso.pw`` code,
and an ``aiida-pseudo`` family. The identifiers below are placeholders. Replace
them with values from your own AiiDA profile.

Choose the calculation settings
-------------------------------

A calculation configuration describes what Quantum ESPRESSO should calculate.
This first graph contains one bulk relaxation. The numerical settings make the
example concrete; choose converged cutoffs and k-point spacing for your own
system.

Quantum ESPRESSO reads the namelist values in its native units. In this example,
``ecutwfc``, ``ecutrho``, and ``conv_thr`` are in Ry, while
``kpoints_distance`` is in Å⁻¹.

.. code-block:: python

   import psteros

   qe_parameters = {
       "CONTROL": {"calculation": "relax"},
       "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
       "ELECTRONS": {"conv_thr": 1.0e-8},
       "IONS": {"ion_dynamics": "bfgs"},
   }

   calculation = psteros.QeCalculationConfig(
       code_label="your-qe-code@your-computer",
       pseudo_family="your-pseudo-family",
       parameters=qe_parameters,
       kpoints_distance=0.20,
   )

The registered code selects the actual AiiDA computer. The pseudopotential
family supplies the element-specific potentials, and the parameter mapping is
passed to ``PwBaseWorkChain``.

Choose the execution settings
-----------------------------

An execution policy supplies the scheduler queue, resource request, wall time,
and MPI choice. Set every field for your own AiiDA environment. The current API
retains legacy deployment-specific defaults for compatibility; do not rely on
them for new calculations.

.. code-block:: python

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

The ``computer`` field is descriptive in the current API. Keep it consistent
with the computer in ``code_label`` because the builder does not cross-check
them. ``ExecutionPolicy`` currently emits PBS-style scheduler directives; a
different scheduler needs a validated backend adaptation before submission.

Assemble and inspect the graph
------------------------------

A **recipe** combines the calculation and execution settings. Give the recipe a
stable name, then pass it and one labelled structure to the graph builder.

.. code-block:: python

   recipe = psteros.SurfaceWorkflowConfig(
       backend="qe",
       calculation=calculation,
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
number of calculations the graph may run at once. PS-TEROS currently keeps one
active calculation in each graph; this graph-local limit does not restrict how
many unrelated jobs your cluster can run.

The call leaves ``submit`` at its default value, ``False``. It creates the graph
in your local AiiDA environment but does not request scheduler resources or
start Quantum ESPRESSO.

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
