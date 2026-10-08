.. _reference-thermochemistry:

=========================================
Free energies of the reference systems
=========================================

The phase-diagram analysis uses 0 K DFT energies for the bulk oxide, the
metal and the O\ :sub:`2` molecule. This guide adds the vibrational (and, for
gases, translational, rotational and electronic) contributions of those
**reference systems**, so that

* Δμ\ :sub:`O` can be read as a temperature and an O\ :sub:`2` pressure; and
* the stability window of the oxide uses formation *free* energies at a
  chosen temperature.

The calculations run with VASP (aiida-vasp, ``vasp.v2.vasp``) as
:ref:`blocks <reference-blocks>`; the thermochemistry itself is pure Python
and works with frequencies from any source. Nothing here changes the existing
builders or the phase-diagram functions.

The model
---------

A **molecule** is an ideal gas with a rigid rotor and harmonic vibrations:

.. math::

   G(T, p) = E_\mathrm{DFT} + E_\mathrm{ZPE} + \int_0^T C_p\,\mathrm{d}T
             - T S(T, p) \qquad
   S = S_\mathrm{trans}(T, p) + S_\mathrm{rot} + S_\mathrm{vib} + k_B \ln(2S_e + 1)

It needs the vibrational frequencies, the geometry (from the structure), the
rotational symmetry number (2 for O\ :sub:`2`, H\ :sub:`2`, H\ :sub:`2`\ O)
and the electron spin S\ :sub:`e` (1 for triplet O\ :sub:`2`).

A **solid** is a harmonic crystal: :math:`F(T) = E_\mathrm{DFT} + E_\mathrm{ZPE}
+ U_\mathrm{vib}(T) - T S_\mathrm{vib}(T)`; :math:`pV` is neglected. Its
frequencies come from a **supercell**: the Gamma point of a small cell misses
most of the phonons, and an elemental metal with one atom per cell has only
the three zero acoustic modes there. The vibrational terms of the supercell
are scaled to the cell whose DFT energy is used.

Every result is a ``FreeEnergy`` that keeps the
terms apart (``electronic_energy_ev``, ``zero_point_energy_ev``,
``thermal_enthalpy_ev``, ``entropy_ev_per_k``, ``correction_ev``) and gives
``free_energy_ev``; ``as_dict()`` is ready for a table.

On the psteros oxygen axis :math:`\mu_\mathrm{O} = \tfrac12 E_\mathrm{DFT}(\mathrm{O_2})
+ \Delta\mu_\mathrm{O}`, O\ :sub:`2` gas at :math:`(T, p)` sits at

.. math::

   \Delta\mu_\mathrm{O}(T, p) = \tfrac12\left[G_\mathrm{O_2}(T, p) - E_\mathrm{DFT}(\mathrm{O_2})\right].

With ``include_zero_point=False`` this is the tabulated Δμ\ :sub:`O` of
Reuter and Scheffler (Phys. Rev. B 65, 035406, 2001), which psteros
reproduces within 0.01 eV.

.. _reference-blocks:

Run the reference calculations
------------------------------

Each reference is a ``ReferenceSystem``; every
reference runs the same **blocks** with one VASP recipe. The default blocks
are ``Relax`` → ``Static`` →
``Vibrations``; each block takes the structure of the
block before it (or of the block named in ``structure_from``).

.. code-block:: python

   import psteros

   recipe = psteros.SurfaceWorkflowConfig(
       backend="vasp",
       calculation=psteros.VaspCalculationConfig(
           code_label="VASP-6.4@cluster",
           incar={"encut": 520, "ediff": 1e-7, "ismear": 0, "sigma": 0.05,
                  "ibrion": 2, "nsw": 100, "ediffg": -0.005},
           potential_mapping={"Sn": "Sn_d", "O": "O"},
           kpoints_spacing=0.03,  # aiida-vasp units of 2*pi/A: about 0.19 1/A
       ),
       execution=psteros.ExecutionPolicy(computer="cluster", queue="standard",
                                         resources={"num_machines": 1, "num_mpiprocs_per_machine": 32}),
       name="sno2",
   )
   cell_relax = psteros.CalculationOverride(parameters={"INCAR": {"isif": 3}})
   references = {
       "o2": psteros.ReferenceSystem(
           psteros.triplet_o2_cell(), "gas", symmetry_number=2, spin=1.0,
           override=psteros.CalculationOverride(parameters={"INCAR": {"ispin": 2, "nupdown": 2}},
                                                kpoints_distance=5.0),  # Gamma only
       ),
       "sno2": psteros.ReferenceSystem(psteros.rutile_sno2_bulk(), "solid", supercell=(2, 2, 3),
                                       block_overrides={"relax": cell_relax}),
       "sn": psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid", supercell=(2, 2, 2),
                                     block_overrides={"relax": cell_relax}),
   }
   graph = psteros.build_vasp_reference_workgraph(references, recipe, submit=True)

The INCAR of a block is built in layers, each overriding the one before: the
recipe INCAR, the block defaults, the block's own ``incar``, the reference's
``override`` and its ``block_overrides[<block name>]``. The tags that make a
block what it is come last: ``NSW = 0`` for a static block; ``IBRION``
(5 for a gas, 6 for a solid), ``POTIM``, ``NFREE`` and ``NSW = 1`` for the
vibrations. The vibrations block also defaults to ``ISIF = 2`` and avoids band
parallelisation: displacing atoms lowers the symmetry, and VASP then stops with
"requested a change of the k-point set ... remove the tag NPAR" if the run uses
``NCORE > 1``. A gas gets ``ISYM = 0`` (it keeps the recipe's ``NCORE``; a
single rank per band would pad a small molecule to one band per rank), a solid
``NCORE = 1``. Keep ``NPAR`` out of the recipe INCAR. A different protocol is a different list of blocks:

.. code-block:: python

   blocks = (
       psteros.Relax(name="coarse", incar={"ediffg": -0.02}),
       psteros.Relax(name="fine"),
       psteros.Static(),
       psteros.Vibrations(structure_from="fine", nfree=4),
   )
   psteros.build_vasp_reference_workgraph(references, recipe, blocks=blocks, submit=True)

Frequencies are only as good as the minimum they are computed at: relax to
forces of a few meV/Å and converge the electrons tightly (``EDIFF`` ≤ 1e-7).
An imaginary mode left after the translations, rotations or acoustic modes
are removed raises an error instead of being dropped silently.

The graph names its tasks ``<label>_<block>_vasp``, ``<label>_<block>_energy``
and ``<label>_<block>_frequencies``, and its outputs
``<label>_<block>_energy`` (eV), ``<label>_<block>_structure`` (relaxations),
``<label>_<block>_frequencies`` (cm\ :sup:`-1`, imaginary as negative),
``<label>_<block>_misc`` and ``<label>_<block>_retrieved``.

Read the results
----------------

``reference_results()`` reads every block, also while
the graph runs or after one block failed.
``reference_thermochemistry()`` turns a finished graph
into thermochemistry objects, with the energy of the static block and the
frequencies of the vibrations block:

.. code-block:: python

   import psteros
   from aiida import load_profile

   load_profile()
   systems = psteros.reference_thermochemistry(graph_pk)       # {"o2": IdealGasMolecule, ...}
   for label, terms in psteros.free_energies(systems, 600.0, pressure_bar=0.2).items():
       print(label, terms.as_dict())

   psteros.delta_mu_oxygen_ev(systems["o2"], 600.0, 0.2)         # O2 at 600 K, 0.2 bar on the axis
   psteros.oxygen_pressure_bar(systems["o2"], 600.0, -1.0)       # O2 pressure where Delta mu_O = -1 eV

The same objects can be built from frequencies computed elsewhere, for
example with ``IdealGasMolecule.from_structure`` and a frequency list, or
with ``parse_vasp_frequencies_cm1()`` on an OUTCAR.

Use them in the phase diagram
-----------------------------

**Placing (T, p) on the axis** needs only the O\ :sub:`2` molecule:
``delta_mu_oxygen_ev`` and ``oxygen_pressure_bar`` convert between the
Δμ\ :sub:`O` axis of ``surface_phase_diagram()`` and an
O\ :sub:`2` temperature and pressure; the axis itself keeps the bare DFT
energy of O\ :sub:`2`.

**The stability window at a temperature** follows from the free energies of
the solids, keeping the bare O\ :sub:`2` energy as the axis origin:

.. code-block:: python

   G = psteros.free_energies(systems, 600.0)
   window = psteros.BinaryOxideReferences(
       bulk_energy_ev=G["sno2"].free_energy_ev,
       bulk_composition=bulk_composition,                       # composition of the relaxed SnO2 cell
       oxygen_molecule_energy_ev=systems["o2"].electronic_energy_ev,
       metal_energy_per_atom_ev=G["sn"].free_energy_ev / atoms_in_sn_cell,
   )
   window.oxygen_poor_limit_ev                                  # Delta G_f(600 K) / y on the same axis

.. warning::

   Keep the slab and bulk energies inside γ at the same level of theory.
   The bulk enters γ as :math:`\frac{N_\mathrm{M}}{x} E_\mathrm{bulk}`;
   replacing it by a free energy while the slabs keep their 0 K DFT energy
   adds the bulk's vibrational free energy to γ without the slab's, an error
   of tenths of an eV per formula unit. Until the slabs have vibrational
   free energies too, build the diagram with DFT energies, and use the free
   energies above for the T, p reading of the axis and the stability window.

What this leaves out
--------------------

Anharmonicity and thermal expansion (quasi-harmonic), the electronic entropy
of metals, configurational entropy, and empirical corrections such as the
O\ :sub:`2` binding correction. A correction can be passed explicitly
(``correction_ev`` or ``reference_thermochemistry(..., corrections_ev=...)``)
and is reported as its own term; none is applied by default.
