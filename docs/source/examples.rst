========
Examples
========

The walkthrough below explains one PS-TEROS model in enough detail to make its
assumptions visible. It is an instructional starting point, not a convergence
protocol or a set of transferable parameters for every material.

Rutile SnO2(110) surface-thermodynamics walkthrough
----------------------------------------------------

This example studies one exposed face of rutile tin dioxide, written
``SnO2(110)``. The ``(110)`` is a Miller index: it identifies the orientation
of the crystal face being modelled.

What the terms mean
^^^^^^^^^^^^^^^^^^^

A surface calculation starts with a finite piece of a bulk crystal, called a
**slab**, separated from its periodic images by vacuum. The atoms exposed at
the outside of that slab define a **surface termination**. Cutting the same
bulk crystal at the same orientation can expose different atoms, so different
terminations can have different energies and stability ranges.

PS-TEROS provides three labelled starting models for this particular example:
``o``, ``sno``, and ``sn2o``. The ``o`` model is stoichiometric. The ``sno``
and ``sn2o`` models are obtained by removing one or two equivalent pairs of
oxygen atoms from the two outer faces. Removing matching atoms from both sides
keeps the slab **symmetric**: its top and bottom surfaces represent the same
termination. This matters below because the surface-energy expression then
accounts for two equivalent surfaces.

Build the starting structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

For this ``1x1`` starting-cell configuration, the formulas are ``Sn18O36``,
``Sn18O34``, and ``Sn18O32``. Nine triple layers and 20 Å of vacuum reproduce
this illustrative model; they are not convergence recommendations. Relax the
structures and test slab thickness, vacuum, k-point sampling, and basis or
cutoff settings for the material and property you are studying.

From calculated energies to surface energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For an oxide slab in equilibrium with bulk SnO2, the analysis uses an oxygen
chemical potential defined by

.. math::

   \mu_{\mathrm{O}} = \frac{E(\mathrm{O}_2)}{2} + \Delta\mu_{\mathrm{O}}.

Here :math:`E(\mathrm{O}_2)` is the energy of the chosen oxygen reference
calculation and :math:`\Delta\mu_{\mathrm{O}}` is an offset used to explore
oxygen-poor and oxygen-rich conditions. Keep that reference compatible with
the bulk and slab calculations: do not combine energies from unrelated
methods, pseudopotential families, or numerical settings.

For a symmetric slab, PS-TEROS evaluates

.. math::

   \gamma =
   \frac{
       E_{\mathrm{slab}}
       - N_{\mathrm{Sn}} E_{\mathrm{bulk}}
       - \left(N_{\mathrm{O}} - 2N_{\mathrm{Sn}}\right)\mu_{\mathrm{O}}
   }{2A}.

The terms have the following roles:

* :math:`\gamma` is the surface energy.
* :math:`E_{\mathrm{slab}}` is the calculated energy of the slab.
* :math:`E_{\mathrm{bulk}}` is the energy per SnO2 formula unit in the bulk
  reference calculation.
* :math:`N_{\mathrm{Sn}}` and :math:`N_{\mathrm{O}}` are the numbers of tin
  and oxygen atoms in the slab. The term
  :math:`N_{\mathrm{O}} - 2N_{\mathrm{Sn}}` measures oxygen excess or deficit
  relative to bulk SnO2 stoichiometry.
* :math:`A` is the area of one exposed face. The denominator is :math:`2A`
  because a symmetric slab has two equivalent surfaces.

Analyze a calculated slab
^^^^^^^^^^^^^^^^^^^^^^^^^

After obtaining compatible bulk, oxygen-reference, and slab energies, pass
those values to the pure-Python helper:

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

The helper returns the surface energy in both eV/Å² and J/m². The variable
names in the example stand for results that you calculated and validated for
your own project; they are not supplied reference energies.

Before interpreting a result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a defensible study, record and validate at least the following:

* the slab orientation, termination, thickness, vacuum, and relaxed geometry;
* the bulk and oxygen-reference calculations used in the expression;
* convergence evidence for the chosen numerical settings; and
* the code, pseudopotentials, scheduler configuration, and provenance records
  for the calculations that produced the energies.

This walkthrough explains the SnO2(110) model and its thermodynamic convention.
It does not establish which termination, numerical settings, or computational
resources are appropriate for a different material or research question.
