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
Γ are removed. A future implementation may add the bulk vibrations from
phonopy on a q-point mesh. A molecule loses its 5 (linear) or 6 rotational and
translational zero modes.

Imaginary modes raise an error, because they mean the structure is not at a
minimum. ``imaginary_modes="drop"`` leaves them out, and
``low_frequency_cutoff_cm1`` raises very soft modes to a cutoff (they dominate
the entropy and are the least converged).

Step 1: the vibrations graph
----------------------------

Run it after the relaxations, on the relaxed structures, with the static
recipe of the energies and a tight electronic convergence (VASP
``EDIFF <= 1e-7``, QE ``conv_thr <= 1e-10``). The ``fixed_sites`` of the
overrides select the sites that are *not* displaced:

.. code-block:: python

   static = psteros.SurfaceWorkflowConfig(
       backend="vasp",
       calculation=psteros.VaspCalculationConfig(
           code_label="vasp@my-cluster",
           incar={"ENCUT": 520, "EDIFF": 1e-7, "ISMEAR": 0, "SIGMA": 0.05, "IBRION": -1, "NSW": 0},
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

Limits
------

* Harmonic, Γ-point only. Bulk dispersion is sampled by the size of the
  supercell (no q-point mesh yet).
* The frozen region of a slab must be bulk-like. Polar slabs with a
  pseudo-hydrogen passivated bottom are not covered yet: their pseudo chemical
  potentials are 0 K energies, so leave them at 0 K.
* :math:`\Delta\mu` is still an axis, not a :math:`(T, p)` pair. The vibrations
  make γ depend on T at each :math:`\Delta\mu`.
