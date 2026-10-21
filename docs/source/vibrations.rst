.. _vibrations:

=========================
Vibrational contributions
=========================

The phase diagrams of :doc:`phase-diagram` use 0 K total energies by default.
This page adds the vibrational free energy of the harmonic approximation, so
that the slab, bulk and reference energies become free energies at a chosen
temperature. The feature is opt-in: leave it out and every result is the
same as before.

The model
---------

For a solid, ``pV`` is negligible and the Gibbs free energy is

.. math::

   G(T) \simeq E_\text{DFT} + F_\text{vib}(T), \qquad
   F_\text{vib}(T) = \sum_i \left[ \frac{h\nu_i}{2}
   + k_B T \ln\left(1 - e^{-h\nu_i / k_B T}\right) \right]

The free energies replace the total energies in the existing references and
terminations. At a fixed temperature, γ stays linear in Δμ, so the binary,
ternary and polar phase diagrams are drawn by the same code.

**Slabs with a frozen centre.** Only the sites that were free during the
relaxation are displaced (a partial Hessian). The frozen centre is counted as
bulk:

.. math::

   F_\text{vib}(\text{slab}) = F_\text{vib}(\text{free sites}) + r\,F_\text{vib}(\text{bulk cell})

where :math:`r` is the number of bulk cells in the frozen region. The frozen
region must therefore have the bulk stoichiometry (whole formula units), which
psteros checks.

**Gas references.** A molecule enters through the energy that fixes
:math:`\Delta\mu = 0`, for example :math:`E(\mathrm{O_2})`. Only its zero-point
energy is added (:math:`E + \mathrm{ZPE}`): the thermal part of the gas
belongs to :math:`\Delta\mu(T, p)` itself.

**The poor limit.** With the vibrations of the bulk and of the other element's
reference (the metal for an oxide), the poor limit becomes
:math:`\Delta G_f(T) / y` instead of :math:`\Delta H_f / y`.

**Frequencies.** The modes come from Γ-point central finite differences of the
forces. A bulk is calculated in a Γ-only supercell, whose 3 acoustic modes at
Γ are removed. Phonons on a q-point mesh are the opt-in phonopy mode
(:ref:`vibrations-phonopy`). A molecule loses its 5 (linear) or 6 rotational and
translational zero modes.

Imaginary modes raise an error, because they mean the structure is not at a
minimum. ``imaginary_modes="drop"`` leaves them out, and
``low_frequency_cutoff_cm1`` raises very soft modes to a cutoff (they dominate
the entropy and are the least converged).

.. _vibrations-step-1:

Step 1: the vibrations graph
----------------------------

Run it after the relaxations, on the relaxed structures, with the static
recipe of the energies (converged, see below) and a tight electronic
convergence (VASP ``EDIFF <= 1e-7``, QE ``conv_thr <= 1e-10``). The ``fixed_sites`` of the
overrides select the sites that are *not* displaced.

**The recipe of the displaced calculations.** Force constants are more
sensitive to the recipe than total energies are. Use a plane-wave cutoff well
above the largest ``ENMAX`` of the POTCARs (QE: a converged ``ecutwfc`` and
``ecutrho``) and, with VASP, ``PREC = Accurate`` with ``ADDGRID = .TRUE.``. A
cheaper recipe that is good enough for energies can shift the modes by several
cm\ :sup:`-1`, by different amounts for surface and bulk atoms, and that error
does not cancel in :math:`\gamma_\text{vib}`. The same holds for a DFT check
of a model's modes (``dft_vibrations`` of the MLIP route); with phonopy,
``PhonopyConfig(plus_minus=True)`` adds central differences.

.. code-block:: python

   static = psteros.SurfaceWorkflowConfig(
       backend="vasp",
       calculation=psteros.VaspCalculationConfig(
           code_label="vasp@my-cluster",
           incar={"ENCUT": 520, "PREC": "Accurate", "ADDGRID": True, "EDIFF": 1e-7, "ISMEAR": 0, "SIGMA": 0.05,
                  "IBRION": -1, "NSW": 0},
           potential_mapping={"Sn": "Sn_d"},
       ),
       execution=psteros.ExecutionPolicy(max_concurrent_jobs=4),  # 1 by default
       role_overrides={
           "slab_o": psteros.CalculationOverride(fixed_sites=fixed_sites_of_the_relaxation),
           "o2": triplet_o2_override,
       },
   )
   graph = psteros.build_vibrations_workgraph(
       {"slab_o": relaxed_slab, "sno2_bulk": relaxed_bulk, "alpha_sn": relaxed_metal, "o2": relaxed_o2},
       static,
       psteros.VibrationsConfig(
           displacement_angstrom=0.01,
           supercells={"sno2_bulk": (2, 2, 2), "alpha_sn": (1, 1, 1)},
           molecules=("o2",),
       ),
       submit=True,
   )

Each displaced site costs 6 static calculations. ``ExecutionPolicy.max_concurrent_jobs``
(1 by default) sets how many run at once. QE calculations run with ``tprnfor``
switched on, and the relaxations' selective dynamics are not passed on to the
displaced calculations. Each label gives a ``<label>_vibrations`` output.

Step 2: free energies in the phase diagram
------------------------------------------

.. code-block:: python

   from aiida import orm

   modes = psteros.read_vibrations(orm.load_node(graph.pk), ["slab_o", "sno2_bulk", "alpha_sn", "o2"])
   T = 800.0  # K
   bulk = modes["sno2_bulk"]

   references = psteros.BinaryOxideReferences(
       bulk_energy_ev=psteros.solid_free_energy_ev(energy["sno2_bulk"], bulk, T),
       bulk_composition={"Sn": 2, "O": 4},
       oxygen_molecule_energy_ev=psteros.molecule_reference_energy_ev(energy["o2"], modes["o2"]),
       metal_energy_per_atom_ev=psteros.solid_free_energy_ev(energy["alpha_sn"], modes["alpha_sn"], T) / 8,
   )
   slab = psteros.SlabTermination.from_structure(
       "slab_o", psteros.solid_free_energy_ev(energy["slab_o"], modes["slab_o"], T, bulk=bulk), relaxed_slab,
   )
   diagram = psteros.surface_phase_diagram([slab], references)

Every free energy is per input cell: a bulk calculated in a supercell is
divided back to its cell. ``HarmonicVibrations`` also gives
``zero_point_energy_ev``, ``free_energy_ev(T)``, ``internal_energy_ev(T)`` and
``entropy_ev_per_k(T)`` for your own analysis, and it can be built directly
from frequencies calculated elsewhere:
``psteros.HarmonicVibrations(frequencies_cm1, composition=..., frozen_composition=...)``.

The complete VASP example is ``examples/vasp_surface_phase_diagram``
(``vibrations.py``, then ``phase_diagram.py --vib-pk ... --temperature ...``).

.. _vibrations-phonopy:

Phonons on a q-point mesh (phonopy)
-----------------------------------

A Γ-point calculation sees a crystal only at the q-points commensurate with
its cell. In the small cells of surface calculations this misses most of the
acoustic branches, which carry much of the vibrational entropy, so the bulk
:math:`F_\text{vib}(T)`, and with it the vibrational change of γ, is not
converged (see the SnO\ :sub:`2` example under :ref:`vibrations-limits`).
``VibrationsConfig(phonopy=psteros.PhonopyConfig(...))`` replaces the Γ-point
modes of every structure that is not a molecule by phonons on a q-point mesh
(``pip install "psteros[phonopy]"``; leave it out and every result is as
before):

* **Force constants.** phonopy displaces one atom at a time in a supercell
  (its symmetry-reduced set, sites of one element with different kind names
  kept apart); one static runs per displacement, plus one of the undisplaced
  supercell, whose forces are subtracted from the others. That subtraction
  removes the residual forces of a geometry that is not the minimum of the
  method that gives the forces (DFT forces at an MLIP geometry).
* **Free energy.** The force constants are interpolated onto a q-point mesh of
  :math:`N_q` points, and the free energy per input cell is the mesh average

  .. math::

     F_\text{vib}(T) = \frac{1}{N_q} \sum_{\mathbf{q}, \nu} \left[ \frac{h\nu_{\mathbf{q}\nu}}{2}
     + k_B T \ln\left(1 - e^{-h\nu_{\mathbf{q}\nu} / k_B T}\right) \right]

  The mesh frequencies go into a ``HarmonicVibrations`` with
  ``supercell_size`` :math:`= N_q`, so ``read_vibrations``,
  ``solid_free_energy_ev`` and ``free_energies_with_vibrations`` work as
  before. The mesh is shifted by half a step on every axis with more than one
  point: Γ, whose three acoustic modes are zero, is not sampled, and
  :math:`F_\text{vib}` converges much faster with the mesh than on a
  Γ-centred one (``shifted_mesh=False``, which leaves the 3 acoustic modes at
  Γ out).
* **Slabs vibrate as a whole** (full-slab phonons): every atom, in an in-plane
  supercell :math:`(n_a, n_b, 1)` with a 2D mesh :math:`(m_a, m_b, 1)`, the
  vacuum along :math:`c`. No region is counted as bulk; the frozen centre of a
  relaxation is only a constraint of that relaxation. For a stoichiometric
  slab of :math:`n` bulk cells

  .. math::

     \gamma_\text{vib}(T) = \frac{F_\text{vib}(\text{slab}) - n\,F_\text{vib}(\text{bulk})}{2A}

  which needs slabs thick enough for the centre to vibrate like the bulk:
  check it with one thicker slab.
* **Molecules** keep their Γ-point modes (5 or 6 zero modes removed).

.. code-block:: python

   vibrations = psteros.VibrationsConfig(
       displacement_angstrom=0.01,
       molecules=("o2",),
       phonopy=psteros.PhonopyConfig(
           supercell_length_angstrom=10.0,  # faces of every supercell at least 10 A apart
           mesh_length_angstrom=40.0,       # m_i = round(40 A |b_i|) q-points per axis
           slabs=("slab_o",),               # in-plane supercell, 2D mesh
       ),
   )
   graph = psteros.build_vibrations_workgraph(
       {"slab_o": relaxed_slab, "sno2_bulk": relaxed_bulk, "alpha_sn": relaxed_metal, "o2": relaxed_o2},
       static, vibrations, submit=True,
   )

The displacements are generated when the graph is built, so the graph holds
one static per displacement (``<label>_vib_p<i>`` and ``<label>_vib_perfect``)
and one calcfunction per label. ``PhonopyConfig.supercells`` and
``PhonopyConfig.meshes`` set them per label instead of the lengths (a slab's
with 1 along :math:`c`); ``plus_minus=True`` adds the opposite of every
displacement (central differences); ``symprec`` is phonopy's symmetry
tolerance. ``displaced_sites`` and the Γ-only ``supercells`` of
``VibrationsConfig`` cannot be combined with phonopy (it displaces every atom
of its own supercell). A slab listed in ``slabs`` needs a vacuum gap along
:math:`c`; a structure with one that is not listed gets a warning.

**Stability.** Every ``<label>_vibrations`` of phonopy also holds the
frequencies searched for imaginary modes on the mesh and at the q-points
commensurate with the supercell (where the force constants are exact, Γ
included): ``psteros.read_phonon_stability(graph, labels)`` returns a
``PhononStability`` per label (``stable``, ``imaginary``,
``lowest_frequency_cm1``, ``summary()``). An imaginary mode at some q means
the structure is not a minimum: a saddle point of the relaxation, or a
distortion that the cell of the calculation could not follow.

**Cost.** A supercell of 10 Å holds 100 to 300 atoms, and a slab of low
symmetry needs a few tens of displacements, so the DFT graph is expensive; on
the MLIP route the GPU job costs seconds. Converge ``supercell_length_angstrom``
and ``mesh_length_angstrom`` with the model before relying on a number, and
check the model against DFT in smaller supercells, with the tight recipe of
:ref:`vibrations-step-1` and ``plus_minus=True``.

**Without AiiDA.** ``phonopy_displacements(structure, supercell)`` gives the
displacements (``PhonopyDisplacements``: the undisplaced and displaced
supercells as pymatgen structures), and ``phonopy_modes(structure,
displacements, forces, reference_forces_ev_per_angstrom=..., mesh=...)`` the
modes (``PhononModes``: ``vibrations``, ``stability``) from forces computed in
any way. ``examples/vibrations_phonopy/cu100_emt.py`` does it with ASE's EMT
potential.

On the MLIP route
-----------------

The MLIP route (:doc:`mlip`) gives the same modes from its fine-tuned MACE
model, on the GPU and without DFT: ``build(..., vibrations=VibrationsConfig(...))``
of ``MlipRoute`` or ``MlipSurfaceStudy`` (or ``build_mlip_surface_workgraph``)
adds one job with the same finite differences at the MACE-relaxed geometries,
the sites held in place in the relaxation being the frozen region. Its
``<label>_mlip_vibrations`` outputs are read by ``read_vibrations``, and
``analyse(graph, temperature_k=T)`` turns every energy into its free energy
with ``free_energies_with_vibrations`` (solids ``E + F_vib(T)``, the frozen
region counted as bulk; gas references ``E + ZPE``). Because the model is
cheap, a bulk supercell much larger than the DFT one costs little. To check
the model's modes, ``dft_vibrations`` builds this page's graph on the MACE
geometries (with QE, VASP or CP2K). ``VibrationsConfig(phonopy=...)``
works on the route too: the daemon generates phonopy's displacements of each
relaxed supercell, one GPU job gives their forces, and the slabs of the study
vibrate as a whole (:ref:`vibrations-phonopy`).

Imaginary modes and the phase diagram
-------------------------------------

An imaginary mode means that the structure is not a minimum of the energy, so
its harmonic free energy is not defined. ``read_vibrations`` and
``HarmonicVibrations`` therefore raise by default and name the modes.
``imaginary_modes="drop"`` leaves them out of the free energy, which is a
choice with consequences: the modes carry a large part of the vibrational
free energy of a soft structure.

The SrTiO\ :sub:`3`\ (001) test of the Lovelace validation shows it. The
cubic perovskite (three imaginary modes at about -60 cm\ :sup:`-1`) and the
TiO\ :sub:`2`-terminated slab (four, down to -207 cm\ :sup:`-1`) are not
minima of the ideal 1x1 structures. With the modes dropped the
TiO\ :sub:`2`-terminated region of the ternary diagram shrinks from 15 % of the
stability region to nothing at 900 K; counting the same modes as real keeps it
at 16 to 20 %. Treat a diagram that depends on dropped modes as unreliable:
relax the structure further (for example along the imaginary mode, in a larger
cell) or compare both choices before drawing conclusions.

Other compounds
---------------

Nothing in the vibrations assumes two elements. For a ternary oxide give the
free energy, not the total energy, of every energy that enters
``TernaryOxideReferences``: the bulk, the energy per atom of each metal, the
oxygen molecule with its zero-point energy, and every competing phase, each with
its own vibrations graph. A slab with a frozen centre needs a frozen region of
whole formula units of the bulk (``solid_free_energy_ev`` refuses anything
else); choose its sites with ``CalculationOverride(fixed_sites=...)`` when
``central_sites`` does not give whole units. ``validation/ternary_lovelace`` is a
complete example for SrTiO\ :sub:`3`.

Running a long graph
--------------------

A displaced site costs 6 calculations, so a few hundred jobs are normal. What
decides whether such a graph finishes:

* **The queue's limit per user.** ``max_concurrent_jobs`` cannot exceed what the
  queue lets one user have: ``qstat -Qf <queue>`` shows ``max_queued`` and
  ``max_run`` per user, shared with all your other jobs and projects. A ``qsub``
  that is refused is retried five times and the calculation is then paused
  (``Pausing after failed transport task: submit_calculation``); play it again
  with ``verdi process play``.
* **One daemon worker with the default SQLite storage.** Two workers on a
  ``core.sqlite_dos`` profile produced ``database is locked`` errors and ran the
  same WorkGraph twice (a displacement calculated twice, and mode tasks started
  with retrieved folders missing). Use one worker (a PostgreSQL profile, which allows
  several, was not tried).
* **psteros installed where the daemon runs**, because the displaced structures
  and the modes are calculated by calcfunctions on the daemon.
* **A daemon host and an SSH connection that stay up.** After a restart of the
  host stop the daemon, run ``verdi process repair`` (it recreates the lost
  process tasks), start the daemon and ``verdi process play --all``. After a
  lost SSH connection restart the daemon too: a worker that was in the middle of
  the dead connection does not always recover by itself.
* **The state of the WorkGraph tasks, not only of the processes.** A task can be
  marked FAILED or SKIPPED (``WorkGraph.load(pk).tasks``) while its calculations
  finished with exit status 0. If a ``<label>_vibrations`` task did not run,
  collect the retrieved folders of its finished calculations and run the same
  ``harmonic_modes`` calcfunction again in a small graph;
  ``validation/ternary_lovelace/rerun_modes.py`` does that and keeps the provenance.

.. _vibrations-limits:

Limits
------


* Harmonic. Without phonopy, Γ-point only: bulk dispersion is sampled by the
  size of the supercell. A small supercell has no acoustic region and
  overestimates the free energy of the bulk against that of a slab's soft
  surface atoms: in the SnO\ :sub:`2`\ (110) test a 1x1x2 bulk gave a
  vibrational change of the surface energy of -26 meV/Å\ :sup:`2` at 900 K for
  the stoichiometric slab, a value that has to be checked against a larger
  supercell before it is used. The phonopy mode (:ref:`vibrations-phonopy`)
  samples a q-point mesh instead. No thermal expansion (quasi-harmonic) yet.
* The frozen region of a slab must be bulk-like. Polar slabs with a
  pseudo-hydrogen passivated bottom are not covered yet: their pseudo chemical
  potentials are 0 K energies, so leave them at 0 K.
* :math:`\Delta\mu` is still an axis, not a :math:`(T, p)` pair. The vibrations
  make γ depend on T at each :math:`\Delta\mu`.
