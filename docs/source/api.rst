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

.. warning::

   ``ExecutionPolicy`` retains legacy deployment-specific defaults for
   compatibility, and omitting ``SurfaceWorkflowConfig.execution`` constructs
   that default policy. Always pass an explicit policy for new calculations;
   the retained defaults are not portable recommendations.

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

``VaspCalculationConfig(code_label, incar, potential_family="PBE", potential_mapping={}, kpoints_spacing=0.20, clean_workdir=False, max_iterations=None)``
   Describe inputs shared by aiida-vasp ``VaspWorkChain`` calculations, the
   central psteros engine.

   ``code_label``
      Full AiiDA label of a registered ``vasp.vasp`` code.

   ``incar``
      Flat mapping of INCAR tags (``{"ENCUT": 520, "IBRION": 2, ...}``) in
      VASP's units. psteros passes it to aiida-vasp in its ``incar``
      namespace; nested mappings are rejected.

   ``potential_family`` and ``potential_mapping``
      The uploaded POTCAR family and the POTCAR per element or kind
      (``{"Sn": "Sn_d"}``). Elements that are not listed use the POTCAR of
      the same name; kinds that are not elements, such as pseudo-hydrogen
      ``H0p75``, must be listed.

   ``kpoints_spacing``
      Positive k-point spacing in Å⁻¹ (with 2π, as in aiida-vasp).

   ``max_iterations``
      Positive maximum number of work-chain restarts; aiida-vasp's default
      when ``None``.

   ``clean_workdir``
      Whether the AiiDA work chain should clean its remote working directory.

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

.. _api-execution-policy:

.. index:: ExecutionPolicy

``ExecutionPolicy(computer=..., queue=..., resources=..., max_concurrent_jobs=1, max_wallclock_seconds=86400, with_mpi=True, prepend_text="")``
   Supply scheduler queue, resource, wall-time, and MPI choices. The registered
   code in ``VaspCalculationConfig`` or ``QeCalculationConfig`` selects the
   actual AiiDA computer. The policy's ``computer`` field is descriptive in the
   current API; keep it consistent with the computer in ``code_label`` because
   the builder does not cross-check them. ``resources`` is a scheduler resource
   mapping accepted by AiiDA, and the wall time is in seconds. ``prepend_text``
   adds shell lines to the job script before the executable, for example
   module loads or ``export QE_MPI_RANKS=88`` for a code whose wrapper launches
   MPI itself (then also pass ``with_mpi=False``).

   The graph builder currently requires ``max_concurrent_jobs=1``.
   ``scheduler_options()`` passes ``queue`` as AiiDA's ``queue_name`` option,
   which the AiiDA scheduler plugin renders in its own syntax (``#PBS -q`` for
   PBS Pro, ``#SBATCH --partition`` for SLURM). It also adds ``#PBS -j oe``,
   a comment to other schedulers.

.. _api-calculation-override:

.. index:: CalculationOverride

``CalculationOverride(parameters=None, kpoints_distance=None, settings={}, metadata={}, fixed_sites=())``
   Describe a deliberate change for one labelled structure. For VASP,
   ``parameters={"INCAR": {...}}`` updates the INCAR of the shared recipe; for
   QE, ``parameters`` deep-merges namelist values. ``kpoints_distance``
   replaces the k-point spacing and must be positive. ``fixed_sites`` holds the
   zero-based indices of sites kept fixed during a relaxation: VASP selective
   dynamics, QE ``FIXED_COORDS``. ``settings`` and ``metadata`` are passed to
   the backend task.

.. _api-surface-workflow-config:

.. index:: SurfaceWorkflowConfig

``SurfaceWorkflowConfig(backend, calculation, execution, name="psteros_surface", role_overrides={})``
   Combine one calculation configuration with an execution policy. Use
   ``backend="vasp"`` with ``VaspCalculationConfig`` or ``backend="qe"`` with
   ``QeCalculationConfig``; any other backend string or a mismatched
   configuration raises an error. ``name`` may contain letters, numbers,
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

.. _api-ev-per-angstrom2-to-j-per-m2:

.. index:: EV_PER_ANGSTROM2_TO_J_PER_M2

``EV_PER_ANGSTROM2_TO_J_PER_M2``
   Conversion factor ``16.02176634`` used by
   ``SurfaceEnergyPoint.gamma_j_per_m2``.

Compatibility boundary
----------------------

The former VASP builders remain under ``psteros.compat`` and ``psteros.core``
for existing projects. New work should use the typed top-level API described
here, which runs VASP through the same aiida-vasp work chain. The
compatibility layer is not part of the current tutorial path.
