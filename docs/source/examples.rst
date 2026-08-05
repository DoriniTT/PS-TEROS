========
Examples
========

The maintained examples use QE through ``aiida-quantumespresso``.  VASP is
available only as the secondary compatibility backend; new examples should
use only the typed public ``psteros`` interface.

SnO2(110) termination structures
--------------------------------

.. code-block:: python

   import psteros

   structures = {
       termination: psteros.sno2_110_slab(
           termination=termination, triple_layers=9, vacuum_angstrom=20.0
       )[0]
       for termination in ("o", "sno", "sn2o")
   }

The 1x1 starting cells have formulas Sn18O36, Sn18O34, and Sn18O32.  They are
symmetrically terminated on both faces; a final relaxation remains the
authority for the calculated geometry.

Surface thermodynamics
----------------------

.. code-block:: python

   point = psteros.surface_energy_oxide_equilibrium(
       slab_energy_ev=slab_energy,
       n_metal=18,
       n_oxygen=36,
       bulk_formula_energy_ev=bulk_energy_per_sno2,
       oxygen_reference_energy_ev=o2_energy,
       delta_mu_oxygen_ev=-1.0,
       surface_area_angstrom2=surface_area,
   )
   print(point.gamma_j_per_m2)

The helper uses ``mu_O = E(O2)/2 + Delta mu_O`` and explicitly divides by two
surfaces.  It is intentionally independent of the workflow engine, allowing
retrieved AiiDA results to be analyzed reproducibly.
