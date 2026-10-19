.. _api:

=========
Reference
=========

Use this page to look up the supported top-level ``psteros`` API. For a guided
introduction, start with the :doc:`tutorial <tutorial>`. The :doc:`concepts page
<concepts>` explains how the objects fit together.

The signatures below emphasize the arguments portable calculations should
supply. Deployment values appear as placeholders instead of copying settings
from one computing environment.

.. note::

   psteros assumes nothing about your computer: code labels, POTCAR or
   pseudopotential families, queues and resources are always yours to give.
   Scheduler options you leave out are not sent, so the AiiDA computer and the
   scheduler apply their own defaults.

Structures
----------

.. _api-rutile-sno2-bulk:

.. index:: rutile_sno2_bulk

``rutile_sno2_bulk(*, a=4.737, c=3.186)``
   Return a pymatgen ``Structure`` for bulk rutile SnO2. ``a`` and ``c`` are the
   tetragonal lattice parameters in Å.

.. _api-sno2-110-slab:

.. index:: sno2_110_slab

``sno2_110_slab(*, termination="o", triple_layers=9, vacuum_angstrom=20.0, a=4.737, c=3.186)``
   Return ``(slab, identity)`` for a symmetric 1x1 rutile SnO2(110) starting
   slab. ``termination`` is ``"o"``, ``"sno"``, or ``"sn2o"`` and
   ``triple_layers`` an odd integer of at least three; other values raise an
   error. ``vacuum_angstrom`` is the positive minimum vacuum requested from
   pymatgen; ``a`` and ``c`` set the bulk lattice used to construct the slab.
   Pass the relaxed bulk lattice so that slab and bulk energies refer to the
   same cell.

.. _api-slab-identity:

.. index:: SlabIdentity

``SlabIdentity(miller_index, termination, triple_layers, vacuum_angstrom)``
   Immutable record returned beside a slab. Its fields preserve the Miller
   index, termination label, layer count, and requested vacuum width.

.. _api-alpha-sn-bulk:

.. index:: alpha_sn_bulk

``alpha_sn_bulk(*, a=6.489)``
   Return a diamond-cubic alpha-Sn starting structure with lattice parameter
   ``a`` in Å.

.. _api-litharge-sno-bulk:

.. index:: litharge_sno_bulk

``litharge_sno_bulk(*, a=3.802, c=4.838, z_sn=0.2367)``
   Return a tetragonal litharge SnO starting structure. ``z_sn`` is the
   fractional tin coordinate along the c axis.

.. _api-triplet-o2-cell:

.. index:: triplet_o2_cell

``triplet_o2_cell(*, cell_length=18.0, bond_length=1.208)``
   Return an O2 molecule in a cubic periodic cell. Lengths are in Å.
   Spin polarization is a calculation input, not a property of this structure.

Calculation configuration
-------------------------

.. _api-vasp-calculation-config:

.. index:: VaspCalculationConfig

``VaspCalculationConfig(code_label, incar=None, potential_family="PBE", potential_mapping={}, kpoints_spacing=0.20, clean_workdir=False, max_iterations=None, precision=None)``
   Describe inputs shared by aiida-vasp ``VaspWorkChain`` calculations, the
   central psteros engine.

   ``code_label``
      Full AiiDA label of a registered ``vasp.vasp`` code.

   ``incar``
      INCAR tags in VASP's units, flat (``{"ENCUT": 520, "IBRION": 2, ...}``) or
      already in aiida-vasp's namespaces (``{"incar": {...}, "dynamics": {...}}``).
      psteros passes the tags to aiida-vasp under ``"incar"`` either way; a flat
      mapping that holds nested mappings is rejected.

   ``potential_family`` and ``potential_mapping``
      The uploaded POTCAR family and the POTCAR per element or kind
      (``{"Sn": "Sn_d"}``). Elements that are not listed use the POTCAR of
      the same name; kinds that are not elements, such as pseudo-hydrogen
      ``H0p75``, must be listed.

   ``kpoints_spacing``
      Positive k-point spacing in Å⁻¹ that includes the 2π, as VASP's
      ``KSPACING`` does: the Γ-centred mesh has ``ceil(|b_i| / kpoints_spacing)``
      points along each reciprocal vector (6×6×6 for GaAs, ``a = 5.75`` Å, at
      0.2 Å⁻¹). The VASP adapter converts it for aiida-vasp, which multiplies
      its own spacing input by 2π.

   ``max_iterations``
      Positive maximum number of work-chain restarts; aiida-vasp's default
      when ``None``.

   ``clean_workdir``
      Whether the AiiDA work chain should clean its remote working directory.

   ``precision``
      ``"fast"``, ``"moderate"`` or ``"precise"`` (:ref:`vasp-precision`): the
      level fills in the cutoff, the k-point spacing, the convergence tags and
      the POTCAR of the elements that ``potential_mapping`` does not name, and
      what you write explicitly wins. ``None`` adds nothing, as before. Leaving
      ``incar`` out selects ``"moderate"``. Unknown names raise a
      ``ValueError`` that lists the three levels.

``PrecisionLevel`` and ``PRECISION_LEVELS``
   ``PRECISION_LEVELS`` maps ``"fast"``, ``"moderate"`` and ``"precise"`` to
   :class:`PrecisionLevel` objects (``encut_ev``, ``kpoints_spacing`` in Å⁻¹ with
   the 2π, ``kpoints_spacing_aiida_vasp``, ``ediff_ev``,
   ``force_cutoff_ev_per_angstrom``, ``ediffg_ev_per_angstrom``) so that you can
   read what a level sets: ``psteros.PRECISION_LEVELS["moderate"].encut_ev``.

.. _api-qe-calculation-config:

.. index:: QeCalculationConfig

``QeCalculationConfig(code_label, pseudo_family, parameters, kpoints_distance=0.20, max_iterations=1, clean_workdir=False)``
   Describe inputs shared by Quantum ESPRESSO ``PwBaseWorkChain`` calculations
   (``pip install '.[qe]'``).

   ``code_label``
      Full AiiDA label of a registered ``quantumespresso.pw`` code.

   ``pseudo_family``
      Label of an ``aiida-pseudo`` family.

   ``parameters``
      Mapping of QE namelists. ``CONTROL``, ``SYSTEM``, and ``ELECTRONS`` are
      required. Relaxation controls such as ``forc_conv_thr``,
      ``etot_conv_thr``, and ``nstep`` belong in ``CONTROL``. Values are passed
      to Quantum ESPRESSO without unit conversion; use its native units. The
      tutorial examples use Ry for cutoffs and energy thresholds and Ry/bohr
      for force thresholds.

   ``kpoints_distance``
      Positive reciprocal-space distance in Å⁻¹ passed to
      ``PwBaseWorkChain``.

   ``max_iterations``
      Positive maximum number of work-chain restart iterations.

   ``clean_workdir``
      Whether the AiiDA work chain should clean its remote working directory.

.. _api-cp2k-calculation-config:

.. index:: Cp2kCalculationConfig

``Cp2kCalculationConfig(code_label, parameters, basis_sets, potentials, basis_set_file="BASIS_MOLOPT", potential_file="GTH_POTENTIALS", kpoints_spacing=0.20, input_files={}, max_iterations=3, clean_workdir=False, vacuum_gap_angstrom=8.0)``
   Describe inputs shared by aiida-cp2k ``Cp2kBaseWorkChain`` calculations
   (``pip install '.[cp2k]'``); see :doc:`cp2k-workflow`.

   ``code_label``
      Full AiiDA label of a registered ``cp2k`` code.

   ``parameters``
      The CP2K input as aiida-cp2k's nested mapping, with upper-case section
      names; ``FORCE_EVAL/DFT`` is required. ``GLOBAL/RUN_TYPE`` is the run
      type (``ENERGY_FORCE`` when left out). The structure (``SUBSYS/CELL``,
      ``COORD``, ``TOPOLOGY``) and ``DFT/KPOINTS`` must not be given: psteros
      writes them. psteros adds one ``KIND`` per kind (unless ``SUBSYS/KIND``
      is given), ``PRINT_LEVEL MEDIUM``, and for ``ENERGY_FORCE`` the force
      file ``psteros_forces.dat`` and the density-of-states file
      ``psteros.dos``.

   ``basis_sets``, ``potentials``
      Element -> basis-set name and element -> pseudopotential name; both
      must name the same elements.

   ``basis_set_file``, ``potential_file``
      Files CP2K reads them from: its data directory, or ``input_files``
      (file name -> local path), which are copied into every job.

   ``kpoints_spacing``
      Positive spacing in Å⁻¹ with the 2π: a Monkhorst-Pack mesh with
      ``ceil(|b_i| / kpoints_spacing)`` points. ``None`` or a 1×1×1 mesh runs
      at Γ only, without a ``KPOINTS`` section.

   ``vacuum_gap_angstrom``
      A lattice vector crossed by an empty gap wider than this (the normal of
      a slab, every direction of a molecule in a box) gets one k-point.
      ``None`` switches this off.

   ``max_iterations``, ``clean_workdir``
      Restarts of the work chain (aiida-cp2k's restart handler is on) and
      cleaning of the remote folder.

.. _api-execution-policy:

.. index:: ExecutionPolicy

``ExecutionPolicy(computer=None, queue=None, max_concurrent_jobs=1, resources={"num_machines": 1}, max_wallclock_seconds=None, with_mpi=True, prepend_text="", extra_options={}, account=None, custom_scheduler_commands="")``
   Scheduler settings of your computer. The registered code in
   ``VaspCalculationConfig`` or ``QeCalculationConfig`` selects the AiiDA
   computer that runs the job; ``computer`` only records its name. ``queue``
   is passed as AiiDA's ``queue_name``, which the scheduler plugin writes in
   its own syntax (``#PBS -q``, ``#SBATCH --partition``, ...). ``resources`` is
   an AiiDA resource mapping, e.g. ``{"num_machines": 1,
   "num_mpiprocs_per_machine": 32}``. ``max_wallclock_seconds`` (seconds),
   ``account`` and ``custom_scheduler_commands`` are passed when given.
   ``prepend_text`` adds shell lines to the job script before the executable,
   such as module loads (pass ``with_mpi=False`` for a code whose wrapper
   launches MPI itself). ``extra_options`` adds or replaces AiiDA scheduler
   options, for example ``{"import_sys_environment": False}`` on clusters whose
   login environment breaks the scheduler parser; it is applied last.

   ``max_concurrent_jobs`` is the largest number of calculations of one graph
   running at once (``None``: no limit); the default of 1 runs them one after
   the other. ``scheduler_options()`` returns the AiiDA metadata options; an
   option that was not given is left out, and the computer's or scheduler's
   default applies.

.. _api-calculation-override:

.. index:: CalculationOverride

``CalculationOverride(parameters=None, kpoints_distance=None, settings={}, metadata={}, fixed_sites=())``
   Describe a deliberate change for one labelled structure. For VASP,
   ``parameters={"INCAR": {...}}`` updates the INCAR of the shared recipe; for
   QE, ``parameters`` deep-merges namelist values; for CP2K it deep-merges
   input sections. ``kpoints_distance`` replaces the recipe's k-point distance
   (QE ``kpoints_distance``, VASP and CP2K ``kpoints_spacing``) and must be
   positive. ``fixed_sites`` holds the zero-based indices of sites kept fixed
   during a relaxation: VASP selective dynamics, QE ``FIXED_COORDS``, CP2K
   ``FIXED_ATOMS``. ``settings`` and ``metadata`` are passed to
   the backend task.

.. _api-surface-workflow-config:

.. index:: SurfaceWorkflowConfig

``SurfaceWorkflowConfig(backend, calculation, execution, name="psteros_surface", role_overrides={})``
   Combine one calculation configuration with an execution policy. Use
   ``backend="vasp"`` with ``VaspCalculationConfig``, ``backend="qe"`` with
   ``QeCalculationConfig`` or ``backend="cp2k"`` with
   ``Cp2kCalculationConfig``; any other backend string or a mismatched
   configuration raises an error. CP2K recipes go to the CP2K builders. ``name`` may contain letters, numbers,
   hyphens, and underscores. ``role_overrides`` maps structure labels to
   ``CalculationOverride`` objects.

.. _api-central-sites:

.. index:: central_sites

``central_sites(structure, half_width=1.5)``
   Indices of the sites within ``half_width`` Å of the mid-plane of a slab,
   for ``CalculationOverride(fixed_sites=...)``.

.. _api-qe-fixed-coordinate-flags:

.. index:: qe_fixed_coordinate_flags

``qe_fixed_coordinate_flags(number_of_sites, fixed_site_indices)``
   Return one three-boolean ``FIXED_COORDS`` row per site. ``True`` means that
   the corresponding Cartesian coordinate is fixed in the AiiDA Quantum
   ESPRESSO plugin. Site indices are zero-based and must lie inside the
   structure.

Graph builders
--------------

.. _api-build-surface-workgraph:

.. index:: build_surface_workgraph

``build_surface_workgraph(structures, config, *, submit=False)``
   Build one AiiDA WorkGraph from a mapping of stable labels to pymatgen
   structures, AiiDA ``StructureData`` nodes, or node PKs. Labels must be
   nonempty and may contain letters, numbers, hyphens, and underscores. Each
   label becomes one backend task. The function returns the graph; with
   ``submit=True`` it also calls ``graph.submit()``.

.. _api-build-relax-static-workgraph:

.. index:: build_relax_static_workgraph

``build_relax_static_workgraph(structures, relaxation, static, *, submit=False)``
   Build a graph in which each static task consumes the relaxed structure
   from its preceding relaxation, with VASP or QE. Both configurations must
   use the same backend and share the same execution policy. The function
   returns the graph; with ``submit=True`` it also submits it.

   With VASP, every relaxation INCAR (recipe plus override) must move the ions
   (``IBRION >= 0`` and ``NSW > 0``) and every static INCAR must not. Graph
   outputs are ``<label>_relaxed_structure``, ``<label>_relax_misc``,
   ``<label>_static_misc`` and the matching ``_retrieved`` folders;
   ``read_vasp_results(graph, labels)`` reads the static energies and relaxed
   structures. ``build_surface_workgraph`` names them ``<label>_misc``,
   ``<label>_structure`` and ``<label>_retrieved``.

   With QE, relaxations run as ``PwRelaxStageWorkChain``, a
   ``PwBaseWorkChain`` that accepts a ``vc-relax`` whose final SCF exceeded
   the force or stress thresholds (exit status 501), so that the static SCF
   still runs on the relaxed structure; psteros must therefore be installed in
   the Python environment of the AiiDA daemon. Graph outputs are
   ``<label>_relaxed_structure``, ``<label>_relax_parameters``,
   ``<label>_static_parameters`` and the ``_retrieved`` folders (with
   ``build_surface_workgraph``: ``<label>_parameters``, ``<label>_structure``
   and ``<label>_retrieved``); ``read_qe_results`` reads them.

   aiida-workgraph attaches graph-level outputs only when every task finished
   successfully; otherwise read them from the work chain called by each task.

.. _api-build-qe-relax-static-workgraph:

.. index:: build_qe_relax_static_workgraph

``build_qe_relax_static_workgraph(structures, relaxation, static, *, submit=False)``
   ``build_relax_static_workgraph`` restricted to ``QeCalculationConfig``
   recipes.

``build_surface_workgraph``, ``build_relax_static_workgraph`` and
``build_vibrations_workgraph`` raise ``ValueError`` for a CP2K recipe.

.. _api-build-cp2k-surface-workgraph:

.. index:: build_cp2k_surface_workgraph

``build_cp2k_surface_workgraph(structures, config, *, submit=False)``
   One CP2K calculation per label, with the recipe's run type
   (``ENERGY_FORCE`` by default). Graph outputs: ``<label>_parameters``,
   ``<label>_structure`` (when CP2K writes one) and ``<label>_retrieved``.

.. _api-build-cp2k-relax-static-workgraph:

.. index:: build_cp2k_relax_static_workgraph

``build_cp2k_relax_static_workgraph(structures, relaxation, static, *, submit=False)``
   The CP2K counterpart of ``build_relax_static_workgraph``. Every relaxation
   (recipe plus override) needs ``GLOBAL/RUN_TYPE`` ``GEO_OPT`` or
   ``CELL_OPT`` and every static calculation ``ENERGY_FORCE`` or ``ENERGY``;
   both recipes share one execution policy. Graph outputs:
   ``<label>_relaxed_structure``, ``<label>_relax_parameters``,
   ``<label>_static_parameters`` and the ``_retrieved`` folders of both
   stages.

.. _api-read-cp2k-results:

.. index:: read_cp2k_results, read_cp2k_forces, read_cp2k_dos, Cp2kDos

``read_cp2k_results(graph, labels)``, ``read_cp2k_forces(graph, labels)``, ``read_cp2k_dos(graph, labels)``
   Read a finished CP2K graph (``orm.load_node(pk)``): energies in eV (the
   static stage when there is one) and structures; forces in eV/Å from the
   force file of each ``ENERGY_FORCE`` run (unit taken from its header); and
   the density of states as ``Cp2kDos``, whose
   ``partially_occupied_bins(tolerance=0.01)`` counts energy bins filled
   between 1% and 99%: more than zero means partly filled states at the
   Fermi level, a metal.

Thermodynamic analysis
----------------------

.. _api-surface-energy-point:

.. index:: SurfaceEnergyPoint

``SurfaceEnergyPoint(delta_mu_oxygen_ev, gamma_ev_per_angstrom2)``
   Immutable value at one chemical-potential offset (of oxygen for an oxide;
   ``delta_mu_ev`` is the element-neutral name). The ``gamma_j_per_m2``
   property converts the stored surface energy to J/m².

.. _api-surface-energy-binary-equilibrium:

.. index:: surface_energy_binary_equilibrium

``surface_energy_binary_equilibrium(*, slab_energy_ev, n_other, n_variable, bulk_formula_energy_ev, variable_reference_energy_ev, delta_mu_ev, surface_area_angstrom2, surfaces=2, formula_unit=(1, 1), reservoir_correction_ev=0.0)``
   Return a ``SurfaceEnergyPoint`` for an A\ :sub:`x`\ B\ :sub:`y` slab in
   equilibrium with its bulk, with μ\ :sub:`B` = reference + Δμ. B is the
   axis element and A the other one; ``formula_unit`` is ``(x, y)``.
   ``reservoir_correction_ev`` is subtracted from the slab energy for any
   further reservoir terms, such as pseudo-hydrogen on a passivated bottom.
   ``surface_energy_oxide_equilibrium`` is the oxide special case.

.. _api-surface-energy-oxide-equilibrium:

.. index:: surface_energy_oxide_equilibrium

``surface_energy_oxide_equilibrium(*, slab_energy_ev, n_metal, n_oxygen, bulk_formula_energy_ev, oxygen_reference_energy_ev, delta_mu_oxygen_ev, surface_area_angstrom2, surfaces=2, formula_unit=(1, 2))``
   Return a ``SurfaceEnergyPoint`` for an M\ :sub:`x`\ O\ :sub:`y` slab in
   equilibrium with its bulk oxide. ``formula_unit`` is ``(x, y)`` for the
   formula unit whose energy is ``bulk_formula_energy_ev``: ``(1, 2)`` for MO2
   (the default), ``(1, 1)`` for MO, ``(2, 3)`` for M2O3. Energies are in eV,
   ``surface_area_angstrom2`` is the area of one exposed face, and ``surfaces``
   is the number of equivalent faces represented by the slab. See the
   :doc:`SnO2 model <examples>` for the equation and the meaning of each term.

.. _api-surface-energy-elemental:

.. index:: surface_energy_elemental

``surface_energy_elemental(*, slab_energy_ev, stoichiometry, chemical_potentials_ev, surface_area_angstrom2, surfaces=2)``
   Return the surface energy in eV/Å² from an explicit element-count
   mapping and matching elemental chemical potentials in eV.

.. _api-stable-termination:

.. index:: stable_termination

``stable_termination(points_by_label)``
   Compare labelled ``SurfaceEnergyPoint`` sequences on a common chemical-
   potential grid. Return ``(delta_mu_oxygen_ev, label,
   gamma_ev_per_angstrom2)`` tuples for the lowest point at each grid value.
   The function raises an error instead of interpolating mismatched grids.

Surface phase diagrams
----------------------

See :doc:`phase-diagram` for a worked introduction.

.. _api-binary-oxide-references:

.. index:: BinaryOxideReferences

``BinaryOxideReferences(bulk_energy_ev, bulk_composition, oxygen_molecule_energy_ev, metal_energy_per_atom_ev=None)``
   Bulk oxide and gas references for a binary oxide M\ :sub:`x`\ O\ :sub:`y`.
   ``bulk_energy_ev`` is the energy of the calculated bulk cell and
   ``bulk_composition`` its element counts, as a mapping or a pymatgen
   ``Composition``; the cell is reduced to one formula unit automatically.
   ``oxygen_molecule_energy_ev`` is E(O2) of a triplet calculation. With
   ``metal_energy_per_atom_ev`` the object also provides
   ``formation_enthalpy_ev`` (per formula unit) and ``oxygen_poor_limit_ev``
   (Δμ\ :sub:`O` = ΔH\ :sub:`f`/y); a non-negative formation enthalpy raises
   an error because no stability window exists. ``metal``, ``formula`` and
   ``formula_unit`` describe the oxide.

.. _api-binary-references:

.. index:: BinaryReferences

``BinaryReferences(bulk_energy_ev, bulk_composition, reference_energies_per_atom_ev, variable=None, reservoir_labels={})``
   References of any binary compound A\ :sub:`x`\ B\ :sub:`y`.
   ``reference_energies_per_atom_ev`` maps each element to its reference
   energy per atom (elemental solid, or half a molecule); the axis element
   ``variable`` (default: the more electronegative one) must be present.
   With the other element's reference too, ``formation_enthalpy_ev`` and
   ``poor_limit_ev`` (Δμ\ :sub:`B` = ΔH\ :sub:`f`/y) are known.
   ``chemical_potentials_ev(delta_mu_ev)`` returns μ\ :sub:`A` and
   μ\ :sub:`B`; ``other``, ``formula`` and ``formula_unit`` describe the
   compound, and ``reservoir_labels`` names the references in figures.

.. _api-slab-termination:

.. index:: SlabTermination

``SlabTermination(label, slab_energy_ev, composition, surface_area_angstrom2, surfaces=2, pseudo_hydrogen={}, bottom_fingerprint=None, face=None)``
   One slab termination. ``surface_area_angstrom2`` is the area of one
   exposed face. ``SlabTermination.from_structure(label, slab_energy_ev,
   structure, *, surfaces=2)`` reads the composition and the area of the plane
   of the first two lattice vectors from a pymatgen structure or an AiiDA
   ``StructureData``. A polar slab with a passivated bottom has
   ``surfaces=1``, its pseudo-hydrogen counts, bottom fingerprint and face;
   ``SlabTermination.from_polar(slab, slab_energy_ev)`` fills them from a
   ``PolarTermination``.

.. _api-surface-phase-diagram:

.. index:: surface_phase_diagram

``surface_phase_diagram(terminations, references, *, delta_mu_range=None, points=201, pseudo_hydrogen=None)``
   Evaluate γ(Δμ) of every termination on a common grid and return a
   ``SurfacePhaseDiagram``; ``references`` is a ``BinaryOxideReferences`` or a
   ``BinaryReferences`` and Δμ is that of its axis element. The default range
   is the stability window, from the poor limit (requires the reference of
   the other element) to Δμ = 0. A wider ``delta_mu_range`` is allowed; the
   window limits are then added to the grid. Terminations must contain only
   the two elements of the compound, with unique labels. Passivated polar
   slabs need ``pseudo_hydrogen`` (a ``PseudoHydrogenReferences``) and, per
   face, one shared bottom.

.. _api-surface-phase-diagram-class:

.. index:: SurfacePhaseDiagram

``SurfacePhaseDiagram``
   The evaluated diagram. ``curves`` maps each label to its
   ``SurfaceEnergyPoint`` tuple, ``stable`` gives the lowest termination at
   each grid point, and ``transitions`` holds the exact
   ``(delta_mu_oxygen_ev, stable_below, stable_above)`` crossings of the lower
   envelope.

   ``plot(path, *, units="J/m2", title=None, dpi=200)``
      Draw γ(Δμ\ :sub:`O`) with the stability window, the transitions and a
      stable-termination strip, and save it; the format follows the file
      suffix (``.png``, ``.pdf``, ``.svg``). ``figure(...)`` returns the
      matplotlib ``Figure`` instead, for further editing.

   ``to_csv(path, *, units="J/m2")``
      Write one row per Δμ value with the columns ``delta_mu_<B>_eV``
      (``delta_mu_O_eV`` for an oxide),
      ``gamma_<label>_Jm2`` (``_eVA2`` with ``units="eV/A2"``) for every
      termination, ``stable_termination`` and ``in_stability_window``.

.. _api-ternary-oxide-references:

.. index:: TernaryOxideReferences, CompetingPhase

``TernaryOxideReferences(bulk_energy_ev, bulk_composition, oxygen_molecule_energy_ev, element_energies_per_atom_ev, competing_phases=(), independent=None)``
   References of a ternary oxide A\ :sub:`x`\ B\ :sub:`y`\ O\ :sub:`z`.
   ``element_energies_per_atom_ev`` maps both non-oxygen elements to their
   reference energy per atom (Δμ = 0). ``competing_phases`` holds
   ``CompetingPhase(label, energy_ev, composition)`` objects, each the energy
   and element counts of a calculated cell. ``independent`` is the element on
   the horizontal axis (alphabetically first by default); the other,
   ``eliminated``, follows from bulk equilibrium through
   ``delta_mu_eliminated_ev(dmu_A, dmu_O)``. The object provides
   ``formation_enthalpy_ev``, ``stability_region`` (polygon vertices in
   (Δμ\ :sub:`A`, Δμ\ :sub:`O`)), ``stability_boundaries`` (the phase limiting
   each edge) and ``in_stability_region(dmu_A, dmu_O)``. It raises an error
   when the oxide is unstable against its elements or the competing phases.

.. _api-ternary-references:

.. index:: TernaryReferences

``TernaryReferences(bulk_energy_ev, bulk_composition, reference_energies_per_atom_ev, competing_phases=(), independent=None, vertical=None, reservoir_labels={})``
   References of any ternary compound A\ :sub:`x`\ B\ :sub:`y`\ C\ :sub:`z`,
   with the same interface as ``TernaryOxideReferences``.
   ``reference_energies_per_atom_ev`` maps all three elements to their
   reference energy per atom. ``vertical`` (C, default: the most
   electronegative element) is on the vertical axis and ``independent``
   (default: the alphabetically first of the other two) on the horizontal one.
   ``reservoir_labels`` names a reference in the figure.

.. _api-ternary-surface-phase-diagram:

.. index:: ternary_surface_phase_diagram, TernarySurfacePhaseDiagram

``ternary_surface_phase_diagram(terminations, references, *, points=101)``
   Evaluate every ``SlabTermination`` over the stability region and return a
   ``TernarySurfacePhaseDiagram``. ``planes`` maps each label to
   ``(gamma_0, d gamma/d dmu_A, d gamma/d dmu_O)`` in eV/Å², ``regions`` to
   the exact polygon where it is the most stable termination (empty if it
   never is). ``gamma_ev_per_angstrom2(label, dmu_A, dmu_O)`` and
   ``stable_termination(dmu_A, dmu_O)`` evaluate single points. ``plot(path,
   *, title=None, dpi=200)`` saves the region map (``figure()`` returns the
   matplotlib ``Figure``), and ``to_csv(path, *, units="J/m2")`` writes one row
   per point of the ``points`` × ``points`` grid with ``delta_mu_<A>_eV``,
   ``delta_mu_<C>_eV`` (``delta_mu_O_eV`` for an oxide), ``delta_mu_<B>_eV``,
   ``gamma_<label>_Jm2`` columns,
   ``stable_termination`` and ``in_stability_region``.

Charge-neutral surfaces
-----------------------

The method and a worked example are in ``docs/CHARGE_NEUTRAL_TERMINATIONS.md``.

.. _api-charge-neutral:

.. index:: find_charge_neutral_terminations, TerminationSet, Termination

``find_charge_neutral_terminations(bulk, miller_index, min_slab_thickness, min_vacuum_thickness=15.0, *, oxidation_states=None, unit_bonds=None, supercell=None)``
   Symmetric slabs with zero net formal charge, repaired by removing
   symmetry-related surface units where needed. Returns a ``TerminationSet``
   (printable table, ``write(directory)``, ``plot(path)``) of
   ``Termination`` objects (``structure``, ``formula``, ``thickness``,
   ``is_stoichiometric``, ``origin``, ``to_dict()``). Polar directions raise
   ``NoChargeNeutralTerminationError``.

.. index:: ChargeNeutralSurfaceStudy, read_qe_results

``ChargeNeutralSurfaceStudy(bulk, miller_indices, references, *, competing_phases={}, min_slab_thickness=12.0, vacuum=15.0, stoichiometric_only=False, oxidation_states=None, unit_bonds=None, supercell=None, termination_options={})``
   The bulk, elemental references, competing phases and every charge-neutral
   slab as labelled ``structures`` (with ``roles``) for
   ``build_surface_workgraph`` or ``build_qe_relax_static_workgraph``;
   ``vasp_overrides(stage)``, ``potential_mapping(base)`` and
   ``qe_overrides(stage)`` give their settings, and
   ``analyse(energies_ev, relaxed_structures=None)`` returns a
   ``ChargeNeutralStudyResult`` with the references and the binary or
   ternary phase diagram. ``read_vasp_results(graph, labels)`` (VASP) and
   ``read_qe_results(graph, labels)`` (QE) read energies and relaxed
   structures from a finished WorkGraph.

Polar surfaces
--------------

The method and a worked example are in ``docs/POLAR_SURFACES.md``.

.. _api-polar:

.. index:: find_polar_terminations, PolarTerminationSet, PolarTermination

``find_polar_terminations(bulk, miller_index, *, bilayers=None, layers=None, vacuum=15.0, oxidation_states=None, electron_counting=True, include_ideal=True, supercell=None, hydrogen_bond_lengths=None, max_variants=20)``
   Slabs of one (hkl) face of a tetrahedrally bonded compound, all on one
   pseudo-hydrogen passivated bottom and in one surface cell: the ideal top
   and the tops that satisfy electron counting. Returns a
   ``PolarTerminationSet`` (printable table, ``write(directory)``,
   ``plot(path)``, ``bottom_fingerprint``) of ``PolarTermination`` objects
   (``structure`` with ``kind_name``, ``pseudo_hydrogen`` and ``bottom`` site
   properties, ``composition`` without pseudo-H, ``pseudo_hydrogen_counts``,
   ``bottom_indices``, ``electron_counting``, ``origin``).

.. index:: pseudo_hydrogens, PseudoHydrogen, pseudo_hydrogen_charge

``pseudo_hydrogens(bulk, oxidation_states=None)``
   ``{element: PseudoHydrogen}`` with ``charge`` (2 − Z/4),
   ``formal_charge``, ``kind_name`` and ``vasp_potential``.
   ``pseudo_hydrogen_charge(element)`` returns the charge alone.

.. index:: PseudoHydrogenReferences, pseudo_molecule, tetrahedral_cluster, fit_cluster_pseudo_chemical_potentials

``PseudoHydrogenReferences.from_pseudo_molecules(energies_ev)``
   Pseudo chemical potentials μ̂ = κ − μ\ :sub:`X`/4 from the energies of the
   pseudo-molecules built by ``pseudo_molecule(bulk, element)``. For the
   cluster method, ``tetrahedral_cluster(bulk, outer, size)`` builds the
   clusters and ``fit_cluster_pseudo_chemical_potentials(outer, energies_ev,
   mu_outer_ev)`` fits Eq. 9 of the Sci. Rep. paper; its ``reference`` goes
   into a ``PseudoHydrogenReferences``.

.. index:: check_bottoms, polar_slab_terminations, eq7_check, nonpolar_check

``check_bottoms(terminations, relaxed_structures, *, reference=None, rmsd_tolerance=0.02, max_tolerance=0.05)``
   Compare the relaxed bottom of each slab with the reference slab and return
   a ``BottomCheckReport`` (``passed``, ``failed``, printable table).
   ``polar_slab_terminations(terminations, energies_ev, relaxed_structures)``
   returns the ``SlabTermination`` objects of the slabs that pass, with the
   report. ``eq7_check`` and ``nonpolar_check`` report the self-consistency of
   the pseudo-hydrogen energies in meV/Å².

.. index:: PolarSurfaceStudy, read_vasp_results

``PolarSurfaceStudy(bulk, faces, references, *, bilayers=9, pseudo_hydrogen_method="molecules", cluster_sizes=(2, 3, 8, 9), eq7_check=True, nonpolar_check=None, electron_counting=True)``
   Every structure of a VASP polar-surface calculation set (``structures``,
   ``roles``), its POTCAR mapping (``potential_mapping(base)``) and INCAR
   overrides (``vasp_overrides(stage)``), and ``analyse(energies_ev,
   relaxed_structures)``, which returns the references, the pseudo chemical
   potentials, the phase diagram, the bottom checks and the consistency
   checks. ``read_vasp_results(graph, labels)`` reads energies and relaxed
   structures from a finished VASP WorkGraph (``build_surface_workgraph`` or
   ``build_relax_static_workgraph``).

Vibrational contributions
-------------------------

See :doc:`vibrations` for the model and a worked example.

.. _api-harmonic-vibrations:

.. index:: HarmonicVibrations

``HarmonicVibrations(frequencies_cm1, composition={}, frozen_composition={}, supercell_size=1, imaginary_modes="raise", low_frequency_cutoff_cm1=None)``
   Harmonic modes of one structure (cm^-1, imaginary ones negative, zero
   modes removed). ``zero_point_energy_ev``, ``free_energy_ev(T)``,
   ``internal_energy_ev(T)`` and ``entropy_ev_per_k(T)`` are per input cell
   (a supercell is divided by ``supercell_size``). ``to_dict()`` and
   ``from_dict(data)`` convert to and from a plain mapping or AiiDA ``Dict``.

.. index:: solid_free_energy_ev, molecule_reference_energy_ev

``solid_free_energy_ev(energy_ev, vibrations, temperature_k, *, bulk=None)``
   ``E + F_vib(T)`` of a slab or bulk, to use in place of its total energy. A
   frozen region (``vibrations.frozen_composition``) is counted as bulk with
   the ``bulk`` vibrations and must have the bulk stoichiometry.
``molecule_reference_energy_ev(energy_ev, vibrations)``
   ``E + ZPE`` of a gas molecule that fixes ``Delta mu = 0`` (e.g. O2).
``free_energies_with_vibrations(energies_ev, vibrations, temperature_k, *, bulk_label, molecule_labels=())``
   Every energy of a diagram as a free energy, by label: ``molecule_labels``
   get ``E + ZPE``, every other label ``E + F_vib(T)`` with the vibrations of
   ``bulk_label`` for a frozen region. Every label needs its vibrations.

.. index:: harmonic_vibrations_from_forces

``harmonic_vibrations_from_forces(*, masses_amu, displaced_sites, displacement_angstrom, forces_plus, forces_minus, composition=None, frozen_composition=None, zero_modes=0, supercell_size=1, imaginary_modes="raise", low_frequency_cutoff_cm1=None)``
   Γ-point modes from central finite differences of forces (eV/Å); a pure
   numpy function, usable with forces from any code.

.. _api-build-vibrations-workgraph:

.. index:: build_vibrations_workgraph, VibrationsConfig, read_vibrations

``build_vibrations_workgraph(structures, static, vibrations=None, *, submit=False)``
   Γ-point finite-difference vibrations of relaxed structures with a VASP or
   QE static recipe: 6 static calculations per displaced site, at most
   ``static.execution.max_concurrent_jobs`` at once. The recipe's
   ``fixed_sites`` are not displaced. Outputs ``<label>_vibrations``.
``VibrationsConfig(displacement_angstrom=0.01, displaced_sites={}, supercells={}, molecules=())``
   Finite-difference step, explicit sites per label, Γ-only supercells per
   label (bulk) and the labels of gas molecules.
``read_vibrations(graph, labels, *, imaginary_modes="raise", low_frequency_cutoff_cm1=None)``
   ``{label: HarmonicVibrations}`` from a finished vibrations graph, or from
   the ``<label>_mlip_vibrations`` of an MLIP graph built with ``vibrations``.

.. _api-ev-per-angstrom2-to-j-per-m2:

.. index:: EV_PER_ANGSTROM2_TO_J_PER_M2

``EV_PER_ANGSTROM2_TO_J_PER_M2``
   Conversion factor ``16.02176634`` used by
   ``SurfaceEnergyPoint.gamma_j_per_m2``.

Reference systems and thermochemistry
-------------------------------------

See :doc:`reference-thermochemistry` for a worked introduction. Frequencies
are in cm\ :sup:`-1` (imaginary modes negative), energies in eV, temperatures
in K and pressures in bar.

.. _api-reference-blocks:

.. index:: Relax, Static, Vibrations

``Relax(name="relax", incar={}, structure_from=None)``, ``Static(name="static", incar={}, structure_from=None)``, ``Vibrations(name="vibrations", incar={}, structure_from=None, ibrion=None, potim=0.015, nfree=2)``
   Calculation blocks, run in order for every reference. ``structure_from``
   names an earlier block whose structure the block takes (by default the
   block before it). ``Static`` forces ``NSW = 0``. ``Vibrations`` runs VASP
   finite differences with ``IBRION`` 5 (gas) or 6 (solid) unless ``ibrion``
   is given, and sets ``POTIM``, ``NFREE`` and ``NSW = 1`` from its fields. Its INCAR
   defaults are ``ISIF = 2`` and, because VASP cannot change its k-point set
   under band parallelisation, ``ISYM = 0`` for a gas and ``NCORE = 1`` for a
   solid.

.. _api-reference-system:

.. index:: ReferenceSystem

``ReferenceSystem(structure, phase, override=None, block_overrides={}, supercell=(1, 1, 1), symmetry_number=None, spin=None)``
   One reference: ``phase`` is ``"gas"`` or ``"solid"``. ``override`` changes
   the recipe for all blocks of this reference and ``block_overrides`` for the
   named block; VASP INCAR tags go under ``"INCAR"`` and ``kpoints_distance``
   is the k-point spacing in the unit of ``VaspCalculationConfig.kpoints_spacing``
   (Å⁻¹ with the 2π; 5 Å⁻¹ gives a Γ-only mesh for a molecule in a box). A solid vibrates in ``supercell``; a
   gas needs ``symmetry_number`` and ``spin`` (total electron spin).

.. _api-build-vasp-reference-workgraph:

.. index:: build_vasp_reference_workgraph

``build_vasp_reference_workgraph(references, config, *, blocks=(Relax(), Static(), Vibrations()), submit=False)``
   Build a VASP graph running ``blocks`` for every labelled
   ``ReferenceSystem`` with the recipe of ``config`` (backend ``"vasp"``).
   Tasks are ``<label>_<block>_vasp``, ``_energy`` and ``_frequencies``;
   outputs ``<label>_<block>_energy``, ``_structure``, ``_frequencies``,
   ``_misc`` and ``_retrieved``. With ``submit=True`` it submits the graph and
   records the references on its node for the two readers below.

.. _api-reference-results:

.. index:: reference_results, reference_thermochemistry

``reference_results(pk)``, ``reference_thermochemistry(pk, *, energy_block=None, vibrations_block=None, imaginary_tolerance_cm1=0.0, corrections_ev=None)``
   ``reference_results`` returns ``{label: {block: result}}`` with ``state``,
   ``pk``, ``energy`` or ``frequencies``, ``structure`` and ``misc``, also for
   running or failed graphs. ``reference_thermochemistry`` returns
   ``{label: IdealGasMolecule | HarmonicSolid}`` from a finished graph.

.. _api-ideal-gas-molecule:

.. index:: IdealGasMolecule, HarmonicSolid, FreeEnergy

``IdealGasMolecule(electronic_energy_ev, frequencies_cm1, masses_amu, positions_angstrom, symmetry_number, spin, geometry=None, correction_ev=0.0, imaginary_tolerance_cm1=0.0)``
   Ideal gas, rigid rotor, harmonic oscillator. ``frequencies_cm1`` is all
   3N modes (translations and rotations are dropped) or only the vibrational
   ones. ``IdealGasMolecule.from_structure(structure, *, electronic_energy_ev,
   frequencies_cm1, symmetry_number, spin, ...)`` reads masses and positions.
   ``free_energy(temperature_k, pressure_bar=1.0)`` returns a ``FreeEnergy``.

``HarmonicSolid(electronic_energy_ev, atoms_in_cell, frequencies_cm1, atoms_in_supercell, correction_ev=0.0, imaginary_tolerance_cm1=0.0)``
   Harmonic crystal whose supercell modes are scaled to the cell of
   ``electronic_energy_ev``. ``free_energy(temperature_k)`` returns a
   ``FreeEnergy``.

``FreeEnergy``
   Immutable terms ``electronic_energy_ev``, ``zero_point_energy_ev``,
   ``thermal_enthalpy_ev``, ``entropy_ev_per_k`` and ``correction_ev``, with
   ``enthalpy_ev``, ``entropy_term_ev`` (−TS), ``free_energy_ev``,
   ``vibrational_free_energy_ev`` and ``as_dict()``.

.. _api-delta-mu-oxygen:

.. index:: delta_mu_oxygen_ev, oxygen_pressure_bar, free_energies, parse_vasp_frequencies_cm1

``delta_mu_oxygen_ev(oxygen, temperature_k, pressure_bar=1.0, *, include_zero_point=True)``, ``oxygen_pressure_bar(oxygen, temperature_k, delta_mu_oxygen, *, include_zero_point=True)``
   Convert between O\ :sub:`2` gas at (T, p) and Δμ\ :sub:`O` on the axis
   ``mu_O = E(O2)/2 + Delta mu_O`` (``E(O2)`` the bare DFT energy).

``free_energies(systems, temperature_k, pressure_bar=1.0)``
   ``{label: FreeEnergy}`` of a mapping of molecules and solids.

``parse_vasp_frequencies_cm1(outcar)``
   Frequencies of the last dynamical matrix in the text of a VASP OUTCAR.

Machine-learned potentials (the MLIP route)
-------------------------------------------

See :doc:`mlip` for the route and its one-source rule. Energies are in eV,
forces in eV/Å; nothing here needs AiiDA.

``LabelledStructure(label, structure, energy_ev, forces_ev_per_angstrom, role="train")``
   One DFT static calculation as a training example (``role``: ``"train"``,
   ``"valid"`` or ``"test"``). ``write_extxyz(frames, destination)`` writes
   frames with ``REF_energy`` and ``REF_forces`` as MACE reads them.

``fit_element_shifts(compositions, residuals_ev)``
   Least-squares per-element shifts of E_DFT − E_foundation; raises when the
   compositions do not separate the elements.

``MaceFineTuneConfig(foundation_model, foundation_model_sha256=None, foundation_head="omat_pbe", replay=True, replay_file=None, replay_min_atoms=2, replay_max_atoms=128, replay_validation_every=5, epochs=100, learning_rate=1e-3, ema_decay=0.999, energy_weight=10.0, forces_weight=10.0, huber_delta=0.01, clip_grad=1.0, seeds=(20260920, 20260921, 20260922), default_dtype="float64", lr_scheduler_gamma=1.0)``
   The fine-tune recipe; ``mace_training_options(config, *, name, seed, atomic_numbers)``
   gives its ``mace_run_train`` options.

``select_replay_structures(path, *, elements, min_atoms=2, max_atoms=128, validation_every=5)``
   Replay geometries (training, validation) and a record of the selection.

``prepare_finetune_job(directory, frames, config, *, name, replay_train=(), replay_valid=(), device="cuda", seeds=None)``, ``prepare_relax_job(directory, structures, models, config=MlipRelaxConfig(), *, fixed_sites=None, device="cuda")``
   Job folders with the inputs and the self-contained scripts ``finetune.py``
   and ``relax.py`` (they need only ase, torch and mace-torch).
``prepare_vibrations_job(directory, structures, model, displaced_sites, *, displacement_angstrom=0.01, default_dtype="float64", device="cuda")``
   The folder of a vibrations job: the model ``(path, head)`` gives the forces
   of every structure with each site of ``displaced_sites[label]`` moved by
   ``±d`` along x, y and z (``vibrations.py``; ``vibrations/<label>.npz``).

``MlipRelaxConfig(fmax_ev_per_angstrom=0.01, max_steps=2000, optimizer="FIRE", head="Default", default_dtype="float64")``, ``relax_with_calculator(structure, calculator, *, config=MlipRelaxConfig(), fixed_sites=())``
   Relax positions with any ASE calculator; returns ``MlipRelaxation``
   (``structure``, ``energy_ev``, ``max_force_ev_per_angstrom``, ``steps``,
   ``converged``).

``partially_occupied_count(occupations, *, tolerance=0.01)``, ``require_metallic(label, *, partially_occupied)``
   The metallic-bulk check: states filled between 1% and 99% of the fullest
   (for CP2K, ``Cp2kDos.partially_occupied_bins()``); ``require_metallic``
   raises ``ValueError`` when there are none.

``check_one_energy_source(sources)``
   Returns the single source of label -> ``"dft"``, ``"mlip"`` or
   ``"mlip_dft"``, or raises ``ValueError`` naming the labels of each source.

``MlipSurfaceStudy(bulk, references, slabs, fixed_half_width_angstrom=1.5, rattled=0, rattle_sigma_angstrom=0.05, rattle_seed=20261017)``
   The MLIP route for your own slabs of a binary compound (see
   :doc:`mlip` and ``docs/MLIP_SURFACES.md``): ``structures`` (``bulk``,
   ``ref_<element>``, the slabs), ``fixed_sites``, ``training_structures``,
   ``molecule_labels``, ``calculation_count(active_learning_rounds=1, dft_check=True, vibrations=False)``;
   ``build(static, fine_tune, relax=MlipRelaxConfig(), *, job, active_learning_rounds=1, dft_check=True, vibrations=None, allow_non_metallic=False, submit=False)``
   returns the graph of ``build_mlip_surface_workgraph`` (the bulk checked
   to be a metal unless ``allow_non_metallic``) and
   ``analyse(source, *, energy_source="mlip_dft", temperature_k=None, vibrations=None, **diagram_options)`` an
   ``MlipStudyResult`` (``references``, ``diagram``, ``terminations``,
   ``energy_source``, ``temperature_k``) from the finished graph or from
   energies typed by hand; with ``temperature_k``, from the free energies of
   the graph's vibrations or of ``vibrations``.
   ``dft_vibrations(source, static, vibrations=None, *, submit=False)`` builds
   DFT finite differences (QE, VASP or CP2K) at the model's geometries.

``MlipRoute(study, fixed_half_width_angstrom=1.5, rattled=0, rattle_sigma_angstrom=0.05, rattle_seed=20261017)``
   The MLIP route for a ``ChargeNeutralSurfaceStudy`` or a
   ``PolarSurfaceStudy``, whose slabs have the right electron count, so the
   bulk is not checked to be a metal: ``structures``, ``slabs``,
   ``fixed_sites`` (the centres of charge-neutral slabs; none for polar
   slabs), ``training_structures``, ``molecule_labels``,
   ``calculation_count(...)``; ``build(static, fine_tune, relax=MlipRelaxConfig(), *, job, active_learning_rounds=1, dft_check=True, vibrations=None, submit=False)``
   (polar: VASP, no vibrations); ``read(source, *, energy_source="mlip_dft")``
   gives ``(energies, relaxed_structures)``;
   ``analyse(source, *, energy_source="mlip_dft", relaxed_structures=None, temperature_k=None, vibrations=None, **options)``
   gives the study's result (``MlipStudyResult`` or ``PolarStudyResult``);
   ``dft_vibrations`` as above.

``build_mlip_surface_workgraph(structures, static, fine_tune, relax=MlipRelaxConfig(), *, job, metallic_bulk, fixed_sites=None, training_structures=None, active_learning_rounds=1, dft_check=True, vibrations=None, allow_non_metallic=False, submit=False)``
   The route as one WorkGraph (DFT statics with any backend, aiida-shell GPU
   jobs). Outputs per label: ``<label>_static_parameters``,
   ``<label>_static_retrieved``, ``<label>_mlip_relaxed_structure``,
   ``<label>_mlip_energy`` (eV) and, with ``dft_check``,
   ``<label>_mlip_dft_energy`` (eV), ``<label>_mlip_dft_parameters``,
   ``<label>_mlip_dft_retrieved``; per round ``r``:
   ``mlip_round<r>_finetune_summary``, ``mlip_round<r>_model_folder``,
   ``mlip_round<r>_relax_results``; with ``vibrations`` (a
   ``VibrationsConfig``) one more GPU job gives ``<label>_mlip_vibrations``
   and ``mlip_vibrations_results``. ``allow_non_metallic=True`` skips the
   metal check of ``metallic_bulk``, with a warning about the electron count
   of the slabs. ``MlipJobConfig(code_label, execution, device="cuda")``
   names the Python of the MACE environment and its scheduler settings
   (``with_mpi=False``).

``read_mlip_results(graph, labels, *, seed=None, round=None)``, ``read_mlip_energies(graph, labels, *, energy_source="mlip_dft")``
   ``MlipRelaxation`` per label (the deployed seed and the last finished round
   by default), and the energies (eV) of one source: ``"mlip_dft"`` or
   ``"mlip"``; both work while the graph runs, for the parts that finished.

``dft_label(backend, *, label, structure, parameters, files, role="train")``, ``dft_diagram_energy(backend, *, label, parameters)``, ``dft_partially_occupied(backend, *, label, files)``
   A training frame from a finished QE, VASP or CP2K static (the free
   energy, whose derivatives are the forces), the energy psteros uses for a
   diagram (VASP: sigma -> 0), and the partly filled states of the bulk.

``rattled_copies(structures, *, count=5, sigma_angstrom=0.05, seed=20261017)``
   ``{"<label>_r<i>": structure}`` with every site moved by a normal vector of
   ``sigma_angstrom`` per component: extra training statics.

``rutile_bulk(metal, *, a, c, u=0.3056)``, ``rutile_110_slab(*, metal, termination="o_bridge", triple_layers=3, vacuum_angstrom=15.0, a, c, u=0.3056, repeat_c=1, cus_bond_angstrom=1.70)``
   Rutile MO\ :sub:`2` and its symmetric (110) slabs with the terminations of
   ``RUTILE_110_TERMINATIONS``: ``"o_bridge"`` (stoichiometric), ``"o_cus"``,
   ``"metal"`` (no bridging O), ``"o_cus_half"`` and ``"o_bridge_half"`` (cell
   doubled along [001]) and ``"o_cus_no_bridge"``.

Compatibility boundary
----------------------

The former VASP builders remain under ``psteros.compat`` and ``psteros.core``
for existing projects. New work should use the typed top-level API described
here, which runs VASP through the same aiida-vasp work chain. The
compatibility layer is not part of the current tutorial path.
