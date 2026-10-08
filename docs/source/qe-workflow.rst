.. _qe-workflow:

========================================
Prepare a QE relaxation-to-static graph
========================================

VASP is the central engine of psteros (see :doc:`vasp-workflow`). The same
recipes, graph builders and analysis also run Quantum ESPRESSO, installed with
``pip install '.[qe]'``. Use this guide when you want a final static Quantum
ESPRESSO calculation to use the geometry from a preceding relaxation. The guide
builds the graph only; it does not submit work.

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

Both stages must use the same execution policy, with the settings of your own
computer and scheduler. The builder turns its queue, resource, wall-time, and MPI choices into AiiDA task
metadata. The registered code in each recipe selects the actual AiiDA computer.
``ExecutionPolicy`` also has a descriptive ``computer`` field; keep it consistent
with ``code_label`` because the builder does not cross-check them.

.. code-block:: python

   import psteros

   execution = psteros.ExecutionPolicy(
       computer="my-cluster",
       queue="my-queue",
       resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
       max_wallclock_seconds=86_400,
       max_concurrent_jobs=1,
   )

``max_concurrent_jobs`` limits how many calculations of the graph run at once;
1 runs them one after the other, ``None`` removes the limit. Options you leave
out, such as ``queue`` or ``max_wallclock_seconds``, are not sent to the
scheduler.

.. note::

   The ``queue`` is passed to AiiDA as ``queue_name``, so the scheduler plugin
   of your AiiDA computer writes it in its own syntax (``#PBS -q``,
   ``#SBATCH --partition``). The accepted ``resources`` keys still depend on
   that plugin; PBS Pro, for instance, turns ``num_cores_per_machine`` into
   ``ncpus``. The job limit applies within one graph: two graphs submitted
   together can run twice as many jobs, which matters on queues with a
   per-user limit.

Define the two recipes
----------------------

The helper below keeps the shared code, pseudopotential family, k-point
spacing, and execution policy in one place. The relaxation-specific force
criterion belongs in ``CONTROL``; psteros validates that placement before a job
is created. The cutoffs and ``conv_thr`` below use Ry, ``forc_conv_thr`` uses
Ry/bohr, and ``kpoints_distance`` uses Å⁻¹.

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
       "CONTROL": {
           "calculation": "scf",
           "tstress": True,
           "tprnfor": True,
       },
       "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
       "ELECTRONS": {"conv_thr": 1.0e-8},
   }

   def qe_recipe(name, parameters):
       return psteros.SurfaceWorkflowConfig(
           backend="qe",
           calculation=psteros.QeCalculationConfig(
               code_label="pw@my-cluster",
               pseudo_family="SSSP/1.3/PBE/efficiency",
               parameters=parameters,
               kpoints_distance=0.20,
           ),
           execution=execution,
           name=name,
       )

   relax = qe_recipe("sno2_110", relax_parameters)
   static = qe_recipe("sno2_110_static", static_parameters)

Connect the stages
------------------

Pass a labelled starting structure and the two recipes to the builder (the
same ``build_relax_static_workgraph`` as for VASP;
``build_qe_relax_static_workgraph`` also accepts only QE recipes). The static
task receives the relaxed ``StructureData`` output from the first task.

.. code-block:: python

   slab, _ = psteros.sno2_110_slab(
       termination="o",
       triple_layers=9,
       vacuum_angstrom=20.0,
   )

   graph = psteros.build_relax_static_workgraph(
       {"sno2_110_o": slab},
       relax,
       static,
   )

   print(graph.name)
   print(graph.max_number_jobs)

You should see:

.. code-block:: text

   sno2_110_relax_static
   1

The name confirms that the builder connected the relaxation and static stages;
the second line is the graph's concurrency limit. This remains an
unsubmitted graph because ``submit`` was not set to ``True``. Check the code
label, pseudopotential family, input parameters, and scheduler metadata before
you request compute time.

Before submission
-----------------

Check the items that the library cannot choose for you:

* the starting structure and any atoms you intend to constrain;
* convergence of cutoffs, k-point sampling, vacuum, and slab thickness;
* compatibility of the pseudopotentials and numerical settings across every
  calculation you plan to compare; and
* the queue, resources, wall-time, and MPI settings accepted by your AiiDA
  computer and scheduler.

A ``vc-relax`` often ends with exit status 501: the ionic cycle converged, but
the final SCF, recomputed with the plane-wave basis of the new cell, exceeds
the stress threshold. The relaxation stage accepts that structure and the
static SCF then evaluates its energy, so the graph continues. The static stage
is also where to read the energy for thermodynamics; the relaxation energy
belongs to the basis of the starting cell. ``psteros.read_qe_results(graph,
labels)`` reads the static energies and relaxed structures of a finished graph.

Atoms are fixed as with VASP, with ``CalculationOverride(fixed_sites=...)``;
psteros writes the ``FIXED_COORDS`` setting of aiida-quantumespresso.

The :doc:`SnO2 surface-energy model <examples>` explains how compatible bulk,
slab, and oxygen-reference energies are combined after the calculations finish,
and :doc:`phase-diagram` turns them into a surface phase diagram.
