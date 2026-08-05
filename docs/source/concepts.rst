.. _concepts:

========================================
How a PS-TEROS calculation fits together
========================================

A surface-energy result is more than a number returned by an
electronic-structure code. It depends on the structure that was modelled, the
method used to calculate its energy, and the references used to compare that
energy. psteros keeps those parts separate so you can inspect them before and
after a run.

Four pieces, one calculation story
----------------------------------

**Starting structures**
   A structure is the geometry given to a calculation. psteros can create common
   starting models, such as bulk rutile SnO2 and symmetric SnO2(110) slabs. A
   starting model is a hypothesis to relax and test; it is not a converged
   result.

**Calculation recipe**
   A ``SurfaceWorkflowConfig`` groups the settings that describe a calculation:
   the backend, code label, pseudopotential family, Quantum ESPRESSO namelists,
   and any per-structure overrides. It answers, “What calculation do I want to
   run?”

**Execution policy**
   An ``ExecutionPolicy`` carries scheduler queue, resource, wall-time, and MPI
   choices. The registered code in the calculation recipe selects the actual
   AiiDA computer. The policy's ``computer`` field is descriptive in the current
   API and should match that code, but the builder does not cross-check them.
   These values belong to your own AiiDA environment and should be supplied
   explicitly in public examples.

**WorkGraph**
   A WorkGraph connects structures and recipes into an ordered set of AiiDA
   tasks. In psteros, a graph is currently limited to one active calculation at
   a time. Building a graph is separate from submitting it, which gives you a
   chance to inspect the plan before it uses compute time.

What AiiDA records
------------------

AiiDA records the inputs, outputs, and links between the calculation steps. This
record is called **provenance**. In practice, provenance lets you trace a final
energy back to the structure, code, pseudopotentials, parameters, and scheduler
metadata that produced it.

psteros adds a small layer on top of that record: it supplies typed recipes,
structure helpers, and a way to connect common surface-calculation stages. It
does not replace an AiiDA profile, a scheduler, or the scientific judgement
needed to choose converged settings.

From structures to surface energy
---------------------------------

A typical surface study has three stages:

1. Build or import bulk and slab structures.
2. Relax the structures and calculate compatible energies with the same
   electronic-structure method.
3. Pass those energies, atom counts, and surface area to a pure analysis helper
   such as ``surface_energy_oxide_equilibrium``.

The analysis helper does not run Quantum ESPRESSO and cannot decide whether the
input calculations are scientifically compatible. Its job is to make the
thermodynamic expression explicit and return the resulting surface energy in
both eV/Å² and J/m².

The :doc:`SnO2 surface-energy model <examples>` explains this sequence with a
symmetric slab, defines the terms used in the equation, and shows what to
record before interpreting a result. If you already know the concepts and want
to prepare a two-stage calculation, continue to the :doc:`QE guide
<qe-first-workflow>`.
