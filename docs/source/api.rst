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

.. _api-qe-calculation-config:

.. index:: QeCalculationConfig

``QeCalculationConfig(code_label, pseudo_family, parameters, kpoints_distance=0.20, max_iterations=1, clean_workdir=False)``
   Describe inputs shared by Quantum ESPRESSO ``PwBaseWorkChain`` calculations.

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

.. _api-vasp-calculation-config:

.. index:: VaspCalculationConfig

``VaspCalculationConfig(code_label, incar, potential_family="PBE", potential_mapping={}, kpoints_spacing=0.20, clean_workdir=False)``
   Describe the VASP configuration retained for established VASP studies.
   ``code_label`` identifies the AiiDA code, ``incar`` stores INCAR settings,
   and the potential fields select the family and optional per-element mapping.
   ``kpoints_spacing`` is a positive reciprocal-space distance in aiida-vasp's
   units of 2π Å⁻¹: ``0.03`` is a spacing of about 0.19 Å⁻¹, and the default
   ``0.20`` (1.26 Å⁻¹) gives a Gamma-only mesh for most cells.
   ``incar`` may be flat (``{"encut": 520, ...}``) or already in aiida-vasp's
   namespaces (``{"incar": {...}, "dynamics": {...}}``); the VASP builders pass
   the tags to aiida-vasp under ``"incar"`` either way.

.. _api-execution-policy:

.. index:: ExecutionPolicy

``ExecutionPolicy(computer=..., queue=..., resources=..., max_concurrent_jobs=1, max_wallclock_seconds=86400, with_mpi=True, prepend_text="", extra_options={})``
   Supply scheduler queue, resource, wall-time, and MPI choices. The registered
   code in ``QeCalculationConfig`` or ``VaspCalculationConfig`` selects the
   actual AiiDA computer. The policy's ``computer`` field is descriptive in the
   current API; keep it consistent with the computer in ``code_label`` because
   the builder does not cross-check them. ``resources`` is a scheduler resource
   mapping accepted by AiiDA, and the wall time is in seconds. ``prepend_text``
   adds shell lines to the job script before the executable, for example
   module loads or ``export QE_MPI_RANKS=88`` for a code whose wrapper launches
   MPI itself (then also pass ``with_mpi=False``). ``extra_options`` adds or
   replaces AiiDA scheduler options, for example
   ``{"import_sys_environment": False}`` on clusters whose login environment
   breaks the scheduler parser.

   The graph builder currently requires ``max_concurrent_jobs=1``.
   ``scheduler_options()`` passes ``queue`` as AiiDA's ``queue_name`` option,
   which the AiiDA scheduler plugin renders in its own syntax (``#PBS -q`` for
   PBS Pro, ``#SBATCH --partition`` for SLURM). It also adds ``#PBS -j oe``,
   a comment to other schedulers.

.. _api-calculation-override:

.. index:: CalculationOverride

``CalculationOverride(parameters=None, kpoints_distance=None, settings={}, metadata={})``
   Describe a deliberate change for one labelled structure. For QE,
   ``parameters`` deep-merges namelist values into the shared recipe; for
   VASP, the tags under ``parameters["INCAR"]`` update the recipe INCAR.
   ``kpoints_distance`` replaces the recipe's k-point distance (QE
   ``kpoints_distance``, VASP ``kpoints_spacing``). ``settings`` and
   ``metadata`` are passed to the backend task. A supplied k-point distance
   must be positive.

.. _api-surface-workflow-config:

.. index:: SurfaceWorkflowConfig

``SurfaceWorkflowConfig(backend, calculation, execution, name="psteros_surface", role_overrides={})``
   Combine one calculation configuration with an execution policy. Use
   ``backend="qe"`` with ``QeCalculationConfig`` or ``backend="vasp"`` with
   ``VaspCalculationConfig``; any other backend string or a mismatched
   configuration raises an error. ``name`` may contain letters, numbers,
   hyphens, and underscores. ``role_overrides`` maps structure labels to
   ``CalculationOverride`` objects.

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

.. _api-build-qe-relax-static-workgraph:

.. index:: build_qe_relax_static_workgraph

``build_qe_relax_static_workgraph(structures, relaxation, static, *, submit=False)``
   Build a QE graph in which each static task consumes the relaxed structure
   from its preceding relaxation. Both configurations must use the QE backend,
   contain ``QeCalculationConfig`` objects, and share the same execution policy.
   The function returns the graph; with ``submit=True`` it also submits it.

   Relaxations run as ``PwRelaxStageWorkChain``, a ``PwBaseWorkChain`` that
   accepts a ``vc-relax`` whose final SCF exceeded the force or stress
   thresholds (exit status 501), so that the static SCF still runs on the
   relaxed structure. psteros must therefore be installed in the Python
   environment of the AiiDA daemon.

   Graph-level outputs are named ``<label>_relaxed_structure``,
   ``<label>_relax_parameters``, ``<label>_static_parameters`` and the
   matching ``_retrieved`` folders. ``build_surface_workgraph`` names them
   ``<label>_parameters``, ``<label>_structure`` and ``<label>_retrieved``.
   aiida-workgraph attaches graph-level outputs only when every task finished
   successfully; otherwise read them from the work chain called by each task.

Thermodynamic analysis
----------------------

.. _api-surface-energy-point:

.. index:: SurfaceEnergyPoint

``SurfaceEnergyPoint(delta_mu_oxygen_ev, gamma_ev_per_angstrom2)``
   Immutable value at one oxygen chemical-potential offset. The
   ``gamma_j_per_m2`` property converts the stored surface energy to J/m².

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

.. _api-slab-termination:

.. index:: SlabTermination

``SlabTermination(label, slab_energy_ev, composition, surface_area_angstrom2, surfaces=2)``
   One symmetric slab termination. ``surface_area_angstrom2`` is the area of
   one exposed face. ``SlabTermination.from_structure(label, slab_energy_ev,
   structure, *, surfaces=2)`` reads the composition and the area of the plane
   of the first two lattice vectors from a pymatgen structure or an AiiDA
   ``StructureData``.

.. _api-surface-phase-diagram:

.. index:: surface_phase_diagram

``surface_phase_diagram(terminations, references, *, delta_mu_range=None, points=201)``
   Evaluate γ(Δμ\ :sub:`O`) of every termination on a common grid and return a
   ``SurfacePhaseDiagram``. The default range is the stability window, from the
   O-poor limit (requires ``metal_energy_per_atom_ev``) to Δμ\ :sub:`O` = 0. A
   wider ``delta_mu_range`` is allowed; the window limits are then added to the
   grid. Terminations must contain only the metal and oxygen, with unique
   labels.

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
      Write one row per Δμ\ :sub:`O` value with the columns ``delta_mu_O_eV``,
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
   ``delta_mu_O_eV``, ``delta_mu_<B>_eV``, ``gamma_<label>_Jm2`` columns,
   ``stable_termination`` and ``in_stability_region``.

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
   defaults are ``ISIF = 2`` and ``NCORE = 1`` (VASP cannot change its k-point
   set under band parallelisation, which displaced cells need).

.. _api-reference-system:

.. index:: ReferenceSystem

``ReferenceSystem(structure, phase, override=None, block_overrides={}, supercell=(1, 1, 1), symmetry_number=None, spin=None)``
   One reference: ``phase`` is ``"gas"`` or ``"solid"``. ``override`` changes
   the recipe for all blocks of this reference and ``block_overrides`` for the
   named block; VASP INCAR tags go under ``"INCAR"`` and ``kpoints_distance``
   is the aiida-vasp ``kpoints_spacing``. A solid vibrates in ``supercell``; a
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

Compatibility boundary
----------------------

The pre-1.0 VASP builders remain under ``psteros.compat`` for existing
projects. New work should use the typed top-level API described here. The
compatibility layer is not part of the current tutorial path.
