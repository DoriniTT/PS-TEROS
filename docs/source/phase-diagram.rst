.. _phase-diagram:

==============================
Build a surface phase diagram
==============================

A surface phase diagram shows, for every candidate termination, the surface
free energy γ as a function of a chemical potential (Δμ\ :sub:`O` for an
oxide, Δμ\ :sub:`As` for GaAs, ...), and which termination is the most stable
at each value. The oxide case is described first; :ref:`other binary
compounds <phase-diagram-binary>` follow the same steps. This guide turns
calculated energies into that diagram and writes it in two forms:

* a **figure** drawn by psteros, ready to inspect or publish; and
* a **CSV table** with every curve, for plotting in the tool of your choice.

The analysis is pure Python. It needs no AiiDA profile and works with energies
from any source; the :ref:`last section <phase-diagram-from-graphs>` shows how
to read them from psteros graphs.

Before you start
----------------

You need, from calculations with compatible settings:

* the energy of the relaxed **bulk oxide** cell and its composition;
* the energy of a **triplet O**\ :sub:`2` molecule;
* the energy per atom of the **elemental metal**, which sets the O-poor limit;
* for every **termination**, the static energy and structure of a relaxed,
  symmetric slab cut from the relaxed bulk lattice.

The model
---------

For a binary oxide M\ :sub:`x`\ O\ :sub:`y` in equilibrium with its bulk,
psteros evaluates

.. math::

   \gamma(\Delta\mu_\mathrm{O}) =
   \frac{E_\mathrm{slab} - \frac{N_\mathrm{M}}{x} E_\mathrm{bulk}
         - \left(N_\mathrm{O} - \frac{y}{x} N_\mathrm{M}\right) \mu_\mathrm{O}}{2A},
   \qquad
   \mu_\mathrm{O} = \tfrac{1}{2} E(\mathrm{O_2}) + \Delta\mu_\mathrm{O},

where :math:`E_\mathrm{bulk}` is the energy per formula unit and :math:`A` the
area of one face. Stoichiometric slabs give a horizontal line; oxygen-poor
slabs rise with Δμ\ :sub:`O`. The oxide is stable between

* the **O-poor limit** :math:`\Delta\mu_\mathrm{O} = \Delta H_f / y`, below
  which it decomposes into the metal, with
  :math:`\Delta H_f = E_\mathrm{bulk} - x E_\mathrm{M} - \frac{y}{2} E(\mathrm{O_2})`;
* the **O-rich limit** :math:`\Delta\mu_\mathrm{O} = 0`, beyond which O\ :sub:`2`
  would condense.

Describe the references
-----------------------

.. code-block:: python

   import psteros

   references = psteros.BinaryOxideReferences(
       bulk_energy_ev=-6711.2091,          # Sn2O4 conventional cell
       bulk_composition={"Sn": 2, "O": 4},  # or bulk_structure.composition
       oxygen_molecule_energy_ev=-1129.8110,
       metal_energy_per_atom_ev=-2220.9247,
   )
   print(references.formula, references.formation_enthalpy_ev, references.oxygen_poor_limit_ev)

The bulk cell is reduced to one formula unit automatically, so the same call
works for MO, MO\ :sub:`2`, M\ :sub:`2`\ O\ :sub:`3` and other binary oxides.
Without ``metal_energy_per_atom_ev`` the O-poor limit is unknown and you must
choose the Δμ\ :sub:`O` range yourself.

Describe the terminations
-------------------------

.. code-block:: python

   terminations = [
       psteros.SlabTermination.from_structure("O-bridge", -20130.8076, slab_o),
       psteros.SlabTermination.from_structure("-1 O/side", -18996.6211, slab_sno),
       psteros.SlabTermination.from_structure("-2 O/side", -17864.0290, slab_sn2o),
   ]

``from_structure`` reads the composition and the area of one face from a
pymatgen structure or an AiiDA ``StructureData``. With the numbers in hand,
``psteros.SlabTermination(label, energy, {"Sn": 6, "O": 12}, area)`` does the
same. Each slab is assumed symmetric, with two equivalent faces.

Evaluate the diagram
--------------------

.. code-block:: python

   diagram = psteros.surface_phase_diagram(terminations, references)
   for delta_mu, below, above in diagram.transitions:
       print(f"{below} -> {above} at {delta_mu:.3f} eV")

By default the grid spans the stability window with 201 points. Pass
``delta_mu_range=(-3.0, 0.5)`` to look beyond it; the window limits are then
added to the grid. The transitions are the exact crossings of the lowest
curves, not grid estimates.

Write the figure
----------------

.. code-block:: python

   diagram.plot("sno2_110_phase_diagram.png", title="SnO$_2$(110)")

The figure shows one line per termination, the O-poor and O-rich limits, the
transitions, and a strip that colours the most stable termination along
Δμ\ :sub:`O`; anything outside the stability window is shaded. The file suffix
selects the format (``.png``, ``.pdf``, ``.svg``) and ``units="eV/A2"``
changes the γ axis. ``diagram.figure()`` returns the matplotlib ``Figure``
instead, if you want to adjust it before saving.

Export the data
---------------

.. code-block:: python

   diagram.to_csv("sno2_110_phase_diagram.csv")

The table has one row per Δμ\ :sub:`O` value:

``delta_mu_O_eV``
   The oxygen chemical potential relative to ½E(O\ :sub:`2`), in eV.
``gamma_<label>_Jm2``
   γ of each termination in J/m² (``gamma_<label>_eVA2`` with
   ``units="eV/A2"``).
``stable_termination``
   The label with the lowest γ at that point.
``in_stability_window``
   ``False`` where the bulk oxide would decompose (below the O-poor limit or
   above Δμ\ :sub:`O` = 0).

It reads directly into a spreadsheet, gnuplot, or pandas:

.. code-block:: python

   import pandas as pd

   table = pd.read_csv("sno2_110_phase_diagram.csv")
   window = table[table.in_stability_window]

.. _phase-diagram-binary:

Other binary compounds
----------------------

Any binary compound A\ :sub:`x`\ B\ :sub:`y` works the same way with
``BinaryReferences``. Give the reference energy per atom of each element,
which fixes μ = reference + Δμ: the elemental solid (Ga, As, Zn, ...) or half
the molecule (O\ :sub:`2`, N\ :sub:`2`). The axis is Δμ of the more
electronegative element unless ``variable`` says otherwise:

.. code-block:: python

   references = psteros.BinaryReferences(
       bulk_energy_ev=e_gaas_cell,                  # relaxed GaAs bulk cell
       bulk_composition={"Ga": 4, "As": 4},
       reference_energies_per_atom_ev={"Ga": e_ga / n_ga, "As": e_as / n_as},
       reservoir_labels={"Ga": "Ga bulk", "As": "As bulk"},  # figure labels
   )
   diagram = psteros.surface_phase_diagram(terminations, references)
   print(references.variable, references.poor_limit_ev)  # As, Delta H_f / y

The model is the one above with B in place of O,
:math:`\mu_\mathrm{B} = \mu_\mathrm{B}^\mathrm{ref} + \Delta\mu_\mathrm{B}`, and
the window runs from the B-poor limit Δμ\ :sub:`B` = ΔH\ :sub:`f`/y to 0.
The CSV column is ``delta_mu_<B>_eV`` and the figure names the references.
For an oxide with an O\ :sub:`2` reference, ``BinaryReferences`` and
``BinaryOxideReferences`` give the same diagram.

.. _phase-diagram-from-graphs:

From psteros graphs
-------------------

``examples/qe_surface_phase_diagram`` in the repository runs the whole
campaign with Quantum ESPRESSO — bulk SnO\ :sub:`2`, α-Sn and O\ :sub:`2`
references, then three SnO\ :sub:`2`\ (110) terminations built on the relaxed
bulk lattice — and its ``phase_diagram.py`` reads the static energies and
relaxed structures from the finished graphs before calling the functions above.

Ternary compounds
-----------------

A ternary oxide A\ :sub:`x`\ B\ :sub:`y`\ O\ :sub:`z` has two independent
chemical potentials. Bulk equilibrium,

.. math::

   x\,\Delta\mu_\mathrm{A} + y\,\Delta\mu_\mathrm{B} + z\,\Delta\mu_\mathrm{O} = \Delta H_f,

eliminates Δμ\ :sub:`B`, so each termination has a surface energy that is a
plane over (Δμ\ :sub:`A`, Δμ\ :sub:`O`), and the diagram becomes a map. The
bulk is stable inside a polygon: no element precipitates (every
Δμ ≤ 0) and no competing phase P forms
(:math:`\sum_i n_i^P \Delta\mu_i \le \Delta H_f^P`). psteros computes that
polygon and the region where each termination is the most stable exactly.

.. code-block:: python

   references = psteros.TernaryOxideReferences(
       bulk_energy_ev=e_bulk,
       bulk_composition={"Sr": 2, "Ti": 2, "O": 6},
       oxygen_molecule_energy_ev=e_o2,
       element_energies_per_atom_ev={"Sr": e_sr, "Ti": e_ti},
       competing_phases=(
           psteros.CompetingPhase("SrO", e_sro, {"Sr": 4, "O": 4}),
           psteros.CompetingPhase("TiO2", e_tio2, {"Ti": 2, "O": 4}),
       ),
       independent="Sr",  # horizontal axis; Ti is eliminated
   )
   diagram = psteros.ternary_surface_phase_diagram(terminations, references)
   diagram.plot("srtio3_001_phase_diagram.png")
   diagram.to_csv("srtio3_001_phase_diagram.csv")

Competing phases are optional but matter: for SrTiO\ :sub:`3` the elements
alone allow a large triangle, while SrO and TiO\ :sub:`2` reduce it to the
narrow strip in which the perovskite actually exists. Give their energies with
the same settings as the other calculations; an oxide that is unstable against
them raises an error naming the phase.

The figure fills each termination's region inside the stability polygon and
names the phase that forms beyond each edge. The CSV has one row per grid point
with ``delta_mu_<A>_eV``, ``delta_mu_O_eV``, ``delta_mu_<B>_eV`` (from bulk
equilibrium), ``gamma_<label>_Jm2`` for every termination,
``stable_termination`` and ``in_stability_region``. The choice of
``independent`` changes the axes, not the physics: the same state has the same
γ in either representation. ``diagram.planes`` gives each plane's coefficients
and ``diagram.regions`` the region polygons.

Other ternary compounds use ``TernaryReferences`` with the reference energy per
atom of all three elements. The vertical axis is the most electronegative
element unless ``vertical`` says otherwise, and the CSV column becomes
``delta_mu_<C>_eV``:

.. code-block:: python

   references = psteros.TernaryReferences(
       bulk_energy_ev=e_bulk,
       bulk_composition={"Cu": 4, "In": 4, "S": 8},
       reference_energies_per_atom_ev={"Cu": e_cu, "In": e_in, "S": e_s},
       competing_phases=(psteros.CompetingPhase("CuS", e_cus, {"Cu": 2, "S": 2}),),
       reservoir_labels={"S": "S8"},   # names the S-rich edge in the figure
   )   # axes: Delta mu_Cu (horizontal), Delta mu_S (vertical); In is eliminated

For an oxide with an O\ :sub:`2` reference, ``TernaryReferences`` and
``TernaryOxideReferences`` give the same diagram.

What the model leaves out
-------------------------

The energies are 0 K total energies: vibrational and configurational
contributions are not included, and Δμ is not converted to a temperature and
pressure. Binary and ternary compounds are supported; slabs must be symmetric
so that both faces are the same termination, unless they are polar slabs with
a passivated bottom described below.
