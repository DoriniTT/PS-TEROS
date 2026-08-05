==============
Worked example
==============

Suppose you want to compare three SnO2(110) surfaces that expose different
numbers of oxygen atoms. Before comparing their energies, you need to know what
each structure represents and how the difference in composition enters the
thermodynamics. This walkthrough follows that small model from its starting
structures to the surface-energy expression.

Rutile SnO2(110) keeps the example compact enough to inspect; its settings are
not a recipe for every oxide.

You will learn what a slab and a surface termination are, why this example uses
a symmetric cell, and which calculated energies the analysis helper needs. The
page does not submit calculations or establish converged parameters for a new
material.

.. important::

   Treat the structures and numerical values below as an illustrative model.
   Relax the structures and converge the slab thickness, vacuum, k-point
   sampling, and basis or cutoff settings for the material and property you
   study.

The physical question
---------------------

A crystal has a repeating bulk structure. To study one of its exposed faces, we
make a finite piece of that crystal and add empty space around it. The finite
piece is called a **slab**; the empty region prevents periodic copies of the
slab from interacting too strongly.

The atoms that appear on the outside of a slab define its **surface
termination**. Cutting the same crystal plane in different places can leave
different atoms exposed, so the resulting terminations can have different
energies.

For this example, ``(110)`` is the Miller index that identifies the crystal
orientation. psteros labels three starting models ``o``, ``sno``, and ``sn2o``:

* ``o`` is the stoichiometric starting surface;
* ``sno`` removes one equivalent pair of oxygen atoms from the outer faces; and
* ``sn2o`` removes two equivalent pairs.

Removing matching atoms from the top and bottom faces keeps the slab
**symmetric**. Both exposed faces then represent the same termination, which is
why the surface-energy expression below divides by two surface areas.

Build the three starting structures
-----------------------------------

The following code creates the three in-memory structures. It does not use an
AiiDA profile or submit work.

.. code-block:: python

   import psteros

   structures = {
       label: psteros.sno2_110_slab(
           termination=label,
           triple_layers=9,
           vacuum_angstrom=20.0,
       )[0]
       for label in ("o", "sno", "sn2o")
   }

For this ``1x1`` starting-cell model, the formulas are ``Sn18O36``,
``Sn18O34``, and ``Sn18O32``. Nine triple layers and 20 Å of vacuum reproduce
this example only. A real study should report the checks used to choose its
own thickness and vacuum.

Choose compatible energy references
-----------------------------------

To compare slabs with different numbers of oxygen atoms, the analysis needs a
way to account for oxygen exchanged with a reservoir. In this model, the oxygen
chemical potential is written as

.. math::

   \mu_{\mathrm{O}} = \frac{E(\mathrm{O}_2)}{2} + \Delta\mu_{\mathrm{O}}.

Here :math:`E(\mathrm{O}_2)` is the energy of the oxygen reference calculation.
The offset :math:`\Delta\mu_{\mathrm{O}}` lets you explore oxygen-poor and
oxygen-rich conditions. Calculate the bulk, oxygen, and slab energies with a
compatible electronic-structure method; combining unrelated functionals,
pseudopotential families, or numerical settings makes the comparison unclear.

Calculate the surface energy
----------------------------

For a symmetric SnO2 slab in equilibrium with bulk SnO2, psteros evaluates

.. math::

   \gamma =
   \frac{
       E_{\mathrm{slab}}
       - N_{\mathrm{Sn}} E_{\mathrm{bulk}}
       - \left(N_{\mathrm{O}} - 2N_{\mathrm{Sn}}\right)\mu_{\mathrm{O}}
   }{2A}.

The expression answers, “How much energy is associated with one exposed
surface?” Its terms are:

* :math:`\gamma`, the surface energy;
* :math:`E_{\mathrm{slab}}`, the calculated energy of the slab;
* :math:`E_{\mathrm{bulk}}`, the energy per SnO2 formula unit in the bulk
  reference calculation;
* :math:`N_{\mathrm{Sn}}` and :math:`N_{\mathrm{O}}`, the numbers of tin and
  oxygen atoms in the slab;
* :math:`N_{\mathrm{O}} - 2N_{\mathrm{Sn}}`, the oxygen excess or deficit
  relative to bulk SnO2; and
* :math:`A`, the area of one exposed face. The denominator is :math:`2A`
  because the slab has two equivalent surfaces.

Use the analysis helper
-----------------------

After you have the compatible energies and the slab area, pass them to the
pure-Python helper. This calculation is local: it does not run Quantum
ESPRESSO or submit an AiiDA process.

.. code-block:: python

   point = psteros.surface_energy_oxide_equilibrium(
       slab_energy_ev=slab_energy,
       n_metal=18,
       n_oxygen=36,
       bulk_formula_energy_ev=bulk_energy_per_sno2,
       oxygen_reference_energy_ev=o2_energy,
       delta_mu_oxygen_ev=-1.0,
       surface_area_angstrom2=surface_area,
       surfaces=2,
   )

   print(point.gamma_ev_per_angstrom2)
   print(point.gamma_j_per_m2)

``point`` provides the surface energy in eV/Å² and J/m². The variables in this
snippet stand for values calculated and checked by your project; they are not
reference energies supplied by psteros.

Before you draw a conclusion
----------------------------

Record the structure, calculation, and analysis choices that make a result
reproducible:

* slab orientation, termination, thickness, vacuum, and relaxed geometry;
* bulk and oxygen-reference calculations used in the expression;
* convergence evidence for numerical settings; and
* the code, pseudopotentials, scheduler configuration, and AiiDA provenance
  records that produced the energies.

This is a model of clean, symmetric SnO2(110) terminations. It does not decide
which termination is stable for another material, nor does it include every
possible temperature, pressure, adsorbate, vibrational, or reconstruction
effect.

For the workflow that creates a relaxation followed by a static calculation,
see the :doc:`QE guide <qe-first-workflow>`. For the role of each object in the
larger calculation story, return to :doc:`the concepts page <concepts>`.
