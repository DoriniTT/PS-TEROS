.. _tutorial:

=========================
Your first PS-TEROS graph
=========================

Here you will build one VASP graph from a bulk SnO2 structure. The
graph stays unsubmitted, so you can inspect how the pieces fit together without
starting a calculation.

By the end, you will know how to:

* separate the calculation inputs from the scheduler settings;
* combine them into a recipe for one labelled structure; and
* confirm that the resulting graph remains unsubmitted.

Before you start
----------------

Complete the :doc:`installation steps <installation>` first. To build the graph,
you also need an active AiiDA profile, a registered ``vasp.vasp`` code, and an
uploaded POTCAR family. The identifiers below are placeholders. Replace them
with values from your own AiiDA profile.

Choose the calculation settings
-------------------------------

A calculation configuration describes what VASP should calculate. This first
graph contains one bulk relaxation. The numerical settings make the example
concrete; choose a converged cutoff and k-point spacing for your own system.

The INCAR is a flat mapping of VASP tags in VASP's units (``ENCUT`` in eV,
``EDIFF`` in eV); ``kpoints_spacing`` is in Å⁻¹, as in aiida-vasp.

.. code-block:: python

   import psteros

   incar = {
       "ENCUT": 520,
       "PREC": "Accurate",
       "EDIFF": 1.0e-6,
       "ISMEAR": 0,
       "SIGMA": 0.05,
       "IBRION": 2,
       "NSW": 100,
       "ISIF": 3,
       "EDIFFG": -0.01,
   }

   calculation = psteros.VaspCalculationConfig(
       code_label="your-vasp-code@your-computer",
       incar=incar,
       potential_family="your-potcar-family",
       potential_mapping={"Sn": "Sn_d"},
       kpoints_spacing=0.20,
   )

The registered code selects the actual AiiDA computer. The POTCAR family
supplies the potentials; ``potential_mapping`` chooses a POTCAR per element
(here ``Sn_d``), and elements that are not listed, like O, use the POTCAR of
the same name. psteros passes the INCAR to aiida-vasp's ``VaspWorkChain``.

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
them. The queue is passed to AiiDA as ``queue_name``, which your computer's
scheduler plugin writes in its own syntax (PBS, SLURM, and so on).

Assemble and inspect the graph
------------------------------

A **recipe** combines the calculation and execution settings. Give the recipe a
stable name, then pass it and one labelled structure to the graph builder.

.. code-block:: python

   recipe = psteros.SurfaceWorkflowConfig(
       backend="vasp",
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
start VASP.

What you have now
-----------------

You have an inspectable graph for one labelled bulk calculation. Before turning
it into a submitted calculation, check that the code label, POTCAR family,
numerical settings, queue, and resource request are appropriate for your
project.

The :doc:`VASP guide <vasp-workflow>` shows how to connect a geometry
relaxation to a final static calculation. Read :doc:`how a PS-TEROS calculation
fits together <concepts>` for the roles of the structures, recipe, graph, and
analysis helpers.
