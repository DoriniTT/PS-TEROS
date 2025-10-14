.. _theory:

======
Theory
======

Ab Initio Atomistic Thermodynamics
----------------------------------

The TEROS package implements ab initio atomistic thermodynamics for calculating surface energies of oxide materials. This approach combines density functional theory (DFT) calculations with thermodynamics to model surface properties under realistic conditions.

Surface Energy Calculation
--------------------------

The surface Gibbs free energy (γ) is calculated as:

.. math::

    \gamma(T, p) = \frac{1}{2A} \left[ G_\text{slab}(T,p,N_i) - \sum_i N_i \mu_i(T,p) \right]

Where:

* :math:`G_\text{slab}` is the Gibbs free energy of the slab
* :math:`A` is the surface area (the factor 2 accounts for the two surfaces of the slab)
* :math:`N_i` is the number of atoms of species :math:`i`
* :math:`\mu_i` is the chemical potential of species :math:`i`

Binary Oxides
-------------

For a binary oxide M₂O (e.g., Ag₂O), the surface energy can be expressed as:

.. math::

    \gamma(T, p_{\text{O}_2}) = \frac{1}{2A} \left[ G_\text{slab} - N_\text{M} \mu_\text{M} - N_\text{O} \mu_\text{O}(T,p_{\text{O}_2}) \right]

The oxygen chemical potential is typically referenced to gaseous O₂:

.. math::

    \mu_\text{O}(T,p_{\text{O}_2}) = \frac{1}{2} \left[ E_{\text{O}_2} + \Delta \mu_{\text{O}_2}(T,p_{\text{O}_2}) \right]

Ternary Oxides
--------------

For ternary oxides (e.g., Ag₃PO₄), the surface energy expression becomes more complex, involving chemical potentials of multiple species and considering formation energies of possible secondary phases.

Phase Diagrams
--------------

Surface stability can be represented in phase diagrams showing the most stable surface terminations as a function of chemical potentials, restricted by thermodynamic constraints to avoid precipitation of secondary phases.
