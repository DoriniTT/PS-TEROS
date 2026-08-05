.. _api:

=========
Reference
=========

Use this page as a map of the supported top-level ``psteros`` API. Start with
the :doc:`tutorial <tutorial>` to learn the flow of a calculation or the
:doc:`concepts page <concepts>` to understand the roles of the objects below.

Structures
----------

``rutile_sno2_bulk()``
   Return a starting bulk rutile SnO2 structure.

``sno2_110_slab(termination=..., triple_layers=..., vacuum_angstrom=...)``
   Build a symmetric 1x1 rutile SnO2(110) starting slab and return the slab
   together with a ``SlabIdentity`` record. Valid terminations are ``"o"``,
   ``"sno"``, and ``"sn2o"``. See the :doc:`worked example <examples>` for the
   physical meaning of these labels.

``alpha_sn_bulk()``, ``litharge_sno_bulk()``, ``triplet_o2_cell()``
   Return reference starting structures used by the SnO2 example. These helpers
   make geometry; they do not run an electronic-structure calculation.

Configuration
-------------

``QeCalculationConfig``
   Describe a Quantum ESPRESSO calculation with a code label, pseudopotential
   family, QE namelists, k-point distance, and optional iteration/cleanup
   choices. ``CONTROL``, ``SYSTEM``, and ``ELECTRONS`` are required namelists.

``VaspCalculationConfig``
   Describe the VASP configuration retained for established VASP studies.

``ExecutionPolicy``
   Describe the AiiDA computer, queue, resources, wall-time, MPI choice, and
   graph-local concurrency bound. Supply deployment-specific values explicitly.
   The public builder currently requires ``max_concurrent_jobs=1`` and emits
   PBS-style custom scheduler directives; another scheduler needs a validated
   backend adaptation before submission.

``SurfaceWorkflowConfig``
   Combine a backend-specific calculation configuration with an execution
   policy. Use ``role_overrides`` when a labelled structure needs a deliberate
   change from the shared recipe.

``CalculationOverride`` and ``qe_fixed_coordinate_flags()``
   Describe a per-structure change or create the QE fixed-coordinate flags for
   a constrained slab relaxation. In the AiiDA Quantum ESPRESSO plugin,
   ``True`` means a coordinate is fixed.

Graph builders
--------------

``build_surface_workgraph(structures, config, submit=False)``
   Build one serial graph from labelled structures and a shared recipe. Set
   ``submit=True`` only when you are ready for the graph to submit work through
   your configured AiiDA environment.

``build_qe_relax_static_workgraph(structures, relaxation, static, submit=False)``
   Build a QE graph in which each static task consumes the relaxed structure
   emitted by its preceding relaxation. The two recipes must share the same
   execution policy.

Thermodynamic analysis
----------------------

``surface_energy_oxide_equilibrium(...)``
   Return a ``SurfaceEnergyPoint`` for an MO2 slab in equilibrium with its bulk
   oxide. The point exposes ``gamma_ev_per_angstrom2`` and the converted
   ``gamma_j_per_m2`` property. See the :doc:`worked example <examples>` for
   the equation, input meanings, and assumptions.

``surface_energy_elemental(...)``
   Return a surface energy from explicit elemental chemical potentials.

``stable_termination(points_by_label)``
   Select the lowest-energy labelled point at every shared chemical-potential
   value. The helper intentionally requires a common input grid rather than
   interpolating values silently.

``EV_PER_ANGSTROM2_TO_J_PER_M2``
   Conversion factor used by ``SurfaceEnergyPoint.gamma_j_per_m2``.

Compatibility boundary
----------------------

The pre-1.0 VASP builders remain under ``psteros.compat`` for existing
projects. New work should use the typed top-level API described here. The
compatibility layer is not part of the current tutorial path.
