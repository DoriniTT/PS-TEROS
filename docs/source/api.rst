.. _api:

=========
Reference
=========

Use this page to look up the supported top-level ``psteros`` API. For a guided
introduction, start with the :doc:`tutorial <tutorial>`. The :doc:`concepts page
<concepts>` explains how the objects fit together.

The signatures below show the arguments readers should supply. Deployment
values such as computer, queue, code, and pseudopotential labels are deliberately
shown as required placeholders rather than inherited defaults.

Structures
----------

``rutile_sno2_bulk(*, a=4.737, c=3.186)``
   Return a pymatgen ``Structure`` for bulk rutile SnO2. ``a`` and ``c`` are the
   tetragonal lattice parameters in Å.

``sno2_110_slab(*, termination="o", triple_layers=9, vacuum_angstrom=20.0, a=4.737, c=3.186)``
   Return ``(slab, identity)`` for a symmetric 1x1 rutile SnO2(110) starting
   slab. ``termination`` accepts ``"o"``, ``"sno"``, or ``"sn2o"``.
   ``triple_layers`` must be an odd integer of at least three.
   ``vacuum_angstrom`` must be positive. ``a`` and ``c`` set the bulk lattice
   used to construct the slab.

``SlabIdentity(miller_index, termination, triple_layers, vacuum_angstrom)``
   Immutable record returned beside a slab. Its fields preserve the Miller
   index, termination label, layer count, and requested vacuum width.

``alpha_sn_bulk(*, a=6.489)``
   Return a diamond-cubic alpha-Sn starting structure with lattice parameter
   ``a`` in Å.

``litharge_sno_bulk(*, a=3.802, c=4.838, z_sn=0.2367)``
   Return a tetragonal litharge SnO starting structure. ``z_sn`` is the
   fractional tin coordinate along the c axis.

``triplet_o2_cell(*, cell_length=18.0, bond_length=1.208)``
   Return an O2 molecule in a cubic periodic cell. Lengths are in Å.
   Spin polarization is a calculation input, not a property of this structure.

Calculation configuration
-------------------------

``QeCalculationConfig(code_label, pseudo_family, parameters, kpoints_distance=0.20, max_iterations=1, clean_workdir=False)``
   Describe inputs shared by Quantum ESPRESSO ``PwBaseWorkChain`` calculations.

   ``code_label``
      Full AiiDA label of a registered ``quantumespresso.pw`` code.

   ``pseudo_family``
      Label of an ``aiida-pseudo`` family.

   ``parameters``
      Mapping of QE namelists. ``CONTROL``, ``SYSTEM``, and ``ELECTRONS`` are
      required. Relaxation controls such as ``forc_conv_thr``,
      ``etot_conv_thr``, and ``nstep`` belong in ``CONTROL``.

   ``kpoints_distance``
      Positive reciprocal-space distance passed to ``PwBaseWorkChain``.

   ``max_iterations``
      Positive maximum number of work-chain restart iterations.

   ``clean_workdir``
      Whether the AiiDA work chain should clean its remote working directory.

``VaspCalculationConfig(code_label, incar, potential_family="PBE", potential_mapping={}, kpoints_spacing=0.20, clean_workdir=False)``
   Describe the VASP configuration retained for established VASP studies.
   ``code_label`` identifies the AiiDA code, ``incar`` stores INCAR settings,
   and the potential fields select the family and optional per-element mapping.
   ``kpoints_spacing`` must be positive.

``ExecutionPolicy(computer=..., queue=..., resources=..., max_concurrent_jobs=1, max_wallclock_seconds=86400, with_mpi=True)``
   Supply scheduler queue, resource, wall-time, and MPI choices. The registered
   code in ``QeCalculationConfig`` or ``VaspCalculationConfig`` selects the
   actual AiiDA computer. The policy's ``computer`` field is descriptive in the
   current API; keep it consistent with the computer in ``code_label`` because
   the builder does not cross-check them. ``resources`` is a scheduler resource
   mapping accepted by AiiDA, and the wall time is in seconds.

   The graph builder currently requires ``max_concurrent_jobs=1``.
   ``scheduler_options()`` emits PBS-style custom scheduler directives. A
   computer using another scheduler needs a validated backend adaptation before
   submission.

``CalculationOverride(parameters=None, kpoints_distance=None, settings={}, metadata={})``
   Describe a deliberate change for one labelled structure. For QE,
   ``parameters`` deep-merges namelist values into the shared recipe.
   ``settings`` and ``metadata`` are passed to the backend task. A supplied
   k-point distance must be positive.

``SurfaceWorkflowConfig(backend, calculation, execution, name="psteros_surface", role_overrides={})``
   Combine one calculation configuration with an execution policy.
   ``backend`` is ``"qe"`` or ``"vasp"`` and must match the calculation type.
   ``name`` may contain letters, numbers, hyphens, and underscores.
   ``role_overrides`` maps structure labels to ``CalculationOverride`` objects.

``qe_fixed_coordinate_flags(number_of_sites, fixed_site_indices)``
   Return one three-boolean ``FIXED_COORDS`` row per site. ``True`` means that
   the corresponding Cartesian coordinate is fixed in the AiiDA Quantum
   ESPRESSO plugin. Site indices are zero-based and must lie inside the
   structure.

Graph builders
--------------

``build_surface_workgraph(structures, config, *, submit=False)``
   Build one AiiDA WorkGraph from a mapping of stable labels to pymatgen
   structures, AiiDA ``StructureData`` nodes, or node PKs. Labels must be
   nonempty and may contain letters, numbers, hyphens, and underscores. Each
   label becomes one backend task. The function returns the graph; with
   ``submit=True`` it also calls ``graph.submit()``.

``build_qe_relax_static_workgraph(structures, relaxation, static, *, submit=False)``
   Build a QE graph in which each static task consumes the relaxed structure
   from its preceding relaxation. Both configurations must use the QE backend,
   contain ``QeCalculationConfig`` objects, and share the same execution policy.
   The function returns the graph; with ``submit=True`` it also submits it.

Thermodynamic analysis
----------------------

``SurfaceEnergyPoint(delta_mu_oxygen_ev, gamma_ev_per_angstrom2)``
   Immutable value at one oxygen chemical-potential offset. The
   ``gamma_j_per_m2`` property converts the stored surface energy to J/m².

``surface_energy_oxide_equilibrium(*, slab_energy_ev, n_metal, n_oxygen, bulk_formula_energy_ev, oxygen_reference_energy_ev, delta_mu_oxygen_ev, surface_area_angstrom2, surfaces=2)``
   Return a ``SurfaceEnergyPoint`` for an MO2 slab in equilibrium with its bulk
   oxide. Energies are in eV, ``surface_area_angstrom2`` is the area of one
   exposed face, and ``surfaces`` is the number of equivalent faces represented
   by the slab. See the :doc:`SnO2 model <examples>` for the equation and the
   meaning of each term.

``surface_energy_elemental(*, slab_energy_ev, stoichiometry, chemical_potentials_ev, surface_area_angstrom2, surfaces=2)``
   Return the surface energy in eV/Å² from an explicit element-count
   mapping and matching elemental chemical potentials in eV.

``stable_termination(points_by_label)``
   Compare labelled ``SurfaceEnergyPoint`` sequences on a common chemical-
   potential grid. Return ``(delta_mu_oxygen_ev, label,
   gamma_ev_per_angstrom2)`` tuples for the lowest point at each grid value.
   The function raises an error instead of interpolating mismatched grids.

``EV_PER_ANGSTROM2_TO_J_PER_M2``
   Conversion factor ``16.02176634`` used by
   ``SurfaceEnergyPoint.gamma_j_per_m2``.

Compatibility boundary
----------------------

The pre-1.0 VASP builders remain under ``psteros.compat`` for existing
projects. New work should use the typed top-level API described here. The
compatibility layer is not part of the current tutorial path.
