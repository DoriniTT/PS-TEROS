================================
Intermediate: Adding Features
================================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

This workflow builds on the :doc:`beginner workflow <beginner-surface-energy>` by enabling three optional features that provide deeper insights into surface properties.

**Prerequisites**: Complete the :doc:`beginner workflow <beginner-surface-energy>` first

**What's New**:

1. **Relaxation energy**: Quantify energetic stabilization from surface relaxation (E_relaxed - E_unrelaxed)
2. **Cleavage energy**: Calculate energy to split bulk into complementary surface pairs
3. **Electronic properties**: DOS and band structure for relaxed bulk structure

**Material**: Ag₂O (same as beginner)

**Additional Runtime**: +1-2 hours (total ~4-6 hours)

Why These Features?
===================

**Relaxation Energy**

* Measures how much energy is gained by allowing surface atoms to relax
* Helps understand which terminations have strong surface reconstruction
* Useful for comparing "frozen" vs. relaxed surface energies

**Formula**:

.. math::

   E_{relax} = E_{relaxed} - E_{unrelaxed}

Negative values indicate stabilization from relaxation.

**Cleavage Energy**

* Energy required to split bulk crystal into two complementary surfaces
* Different from surface energy (which references chemical potential)
* Useful for understanding mechanical surface creation

**Formula**:

.. math::

   E_c(i,j) = \frac{1}{2A} (E_i^{slab} + E_j^{slab} - n \cdot E_{bulk})

Where terminations i and j are complementary pairs.

**Electronic Properties**

* Band structure shows electronic band dispersion along high-symmetry paths
* DOS shows density of states for analyzing gap, band edges, and orbital character
* Essential for understanding semiconducting/insulating behavior

What You'll Learn
=================

* How to enable optional features with boolean flags
* How to configure electronic properties calculations
* How to interpret relaxation and cleavage energies
* How to visualize band structures and DOS

Code Changes from Beginner
===========================

Only three changes are needed to enable all new features:

1. Import the Electronic Properties Builder
-------------------------------------------

Add one import to your beginner script:

.. code-block:: python

    from teros.core.workgraph import build_core_workgraph
    from teros.core.builders import get_electronic_properties_defaults  # NEW

2. Configure Electronic Properties
-----------------------------------

Add this section before ``build_core_workgraph()``:

.. code-block:: python

    # Get electronic properties defaults for bulk
    ep_defaults = get_electronic_properties_defaults(
        energy_cutoff=bulk_parameters['ENCUT'],  # Match bulk ENCUT
        electronic_convergence=1e-5,
        ncore=4,
        ispin=2,  # Spin-polarized for Ag₂O
        lasph=True,
        lreal="Auto",
        kpoints_mesh_density=0.3,     # SCF k-mesh density
        band_kpoints_distance=0.2,    # Band path density
        dos_kpoints_distance=0.2,     # DOS mesh density
        line_density=0.2,             # Points along high-symmetry lines
        nedos=2000,                   # DOS grid points
        sigma_bands=0.01,             # Smearing for bands (eV)
        symprec=1e-4,                 # Symmetry precision
        band_mode="seekpath-aiida",   # Auto high-symmetry paths
    )

**Understanding the Parameters**:

* ``kpoints_mesh_density``: Density of k-point mesh for SCF and DOS (Å⁻¹)
* ``band_kpoints_distance``: Spacing along band structure path (smaller = denser)
* ``line_density``: Points per Å along high-symmetry lines
* ``nedos``: Number of DOS grid points (higher = smoother DOS)
* ``band_mode``: Uses seekpath to automatically find high-symmetry paths

3. Enable New Features
-----------------------

Modify the ``build_core_workgraph()`` call:

.. code-block:: python
   :emphasize-lines: 7-8,9-10,14-17

    wg = build_core_workgraph(
        # ... all beginner parameters stay the same ...

        # ===== NEW: Feature Flags =====
        relax_slabs=True,

        compute_relaxation_energy=True,           # NEW: Enable relaxation energy
        compute_cleavage=True,                    # NEW: Enable cleavage energy
        compute_thermodynamics=True,
        thermodynamics_sampling=10,

        # ===== NEW: Electronic Properties =====
        compute_electronic_properties_bulk=True,  # NEW: Enable DOS/bands
        bands_parameters=ep_defaults,             # NEW: Use builder defaults
        band_settings=ep_defaults['band_settings'], # NEW
        bands_options=bulk_options,               # NEW: Same resources as bulk

        # ... rest of parameters ...
        name='Ag2O_Intermediate_Workflow',
    )

That's it! Three simple changes enable all new features.

Understanding the New Workflow
===============================

The workflow now has additional steps:

.. code-block:: text

    Bulk + Metal + O₂ Relaxation (parallel)
         ↓
    Formation Enthalpy
         ↓
    ┌─────────────────────────────────────┐
    │ NEW: Bulk Electronic Properties     │
    │   1. SCF (LWAVE=True, LCHARG=True) │
    │   2. Band structure (non-SCF)       │
    │   3. DOS (non-SCF, tetrahedron)     │
    └─────────────────────────────────────┘
         ↓
    Generate Slabs
         ↓
    ┌─────────────────────────────────────┐
    │ NEW: Parallel Unrelaxed SCF         │
    │   (NSW=0, IBRION=-1)               │
    └─────────────────────────────────────┘
         ↓
    Parallel Slab Relaxation (existing)
         ↓
    ┌─────────────────────────────────────┐
    │ NEW: Calculate Relaxation Energies  │
    │   E_relax = E_relaxed - E_unrelaxed │
    └─────────────────────────────────────┘
         ↓
    ┌─────────────────────────────────────┐
    │ NEW: Calculate Cleavage Energies    │
    │   For complementary termination pairs│
    └─────────────────────────────────────┘
         ↓
    Surface Thermodynamics (existing)

**Key Additions**:

* Unrelaxed SCF runs in parallel with relaxations
* Electronic properties calculated after bulk relaxation
* Cleavage energies computed for complementary pairs

New Outputs
===========

Relaxation Energies
-------------------

Access relaxation energies:

.. code-block:: python

    from aiida import load_profile
    from aiida.orm import load_node

    load_profile()
    wg = load_node(<YOUR_PK>)

    # Get relaxation energies
    relaxation_energies = wg.outputs.relaxation_energies.get_dict()

    print("Relaxation energies (eV):")
    for term_label, energy_pk in relaxation_energies.items():
        energy = load_node(energy_pk).value
        print(f"  {term_label}: {energy:.3f} eV")

**Interpretation**:

* **Negative values**: Relaxation stabilizes the surface (common)
* **Large magnitude** (e.g., -5 eV): Strong surface reconstruction
* **Small magnitude** (e.g., -0.5 eV): Minimal relaxation

**Typical values for Ag₂O**: -1 to -3 eV per slab (depends on termination)

Cleavage Energies
-----------------

.. code-block:: python

    # Get cleavage energies
    cleavage_energies = wg.outputs.cleavage_energies.get_dict()

    print("\nCleavage energies:")
    for pair_label, data_pk in cleavage_energies.items():
        data = load_node(data_pk).get_dict()
        print(f"  {pair_label}:")
        print(f"    E_c = {data['cleavage_energy_eV_per_A2']:.4f} eV/Ų")
        print(f"        = {data['cleavage_energy_J_per_m2']:.4f} J/m²")
        print(f"    Pair: {data['term_i']} + {data['term_j']}")

**Understanding Pairs**:

* pymatgen pairs complementary terminations (e.g., term_0 + term_1)
* Symmetric slabs pair with themselves (term_i + term_i)
* Cleavage energy is always positive (energy required to create surfaces)

**Typical values**: 1-3 J/m² for oxide surfaces

Bulk Electronic Properties
--------------------------

**Band Structure**:

.. code-block:: python

    # Access band structure
    bulk_bands = wg.outputs.bulk_bands

    # Get seekpath information (high-symmetry points)
    seekpath_params = wg.outputs.bulk_seekpath_parameters.get_dict()

    print("\nHigh-symmetry path:")
    print(f"  Path: {' → '.join(seekpath_params['path'])}")
    print(f"  Point labels: {seekpath_params['point_coords'].keys()}")

**DOS**:

.. code-block:: python

    # Access DOS
    bulk_dos = wg.outputs.bulk_dos

    # DOS data is in AiiDA BandsData format
    # Can be plotted using aiida-quantumespresso tools or custom scripts

**Primitive Structure**:

.. code-block:: python

    # Seekpath standardizes to primitive cell
    primitive_structure = wg.outputs.bulk_primitive_structure

    print(f"Original bulk: {wg.outputs.bulk_structure.get_formula()}")
    print(f"Primitive cell: {primitive_structure.get_formula()}")

Visualization Examples
======================

Plot Relaxation Energy Distribution
------------------------------------

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np

    # Get relaxation energies
    relaxation_energies = wg.outputs.relaxation_energies.get_dict()

    terms = []
    energies = []
    for term_label, energy_pk in sorted(relaxation_energies.items()):
        terms.append(term_label)
        energies.append(load_node(energy_pk).value)

    # Bar plot
    plt.figure(figsize=(8, 5))
    plt.bar(terms, energies, color='steelblue', alpha=0.7)
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    plt.xlabel('Termination')
    plt.ylabel('Relaxation Energy (eV)')
    plt.title('Surface Relaxation Energies for Ag₂O (111)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('relaxation_energies.png', dpi=300)

Compare Surface vs. Cleavage Energies
--------------------------------------

.. code-block:: python

    import matplotlib.pyplot as plt

    # Get cleavage energies (J/m²)
    cleavage_energies = wg.outputs.cleavage_energies.get_dict()

    cleavage_dict = {}
    for pair_label, data_pk in cleavage_energies.items():
        data = load_node(data_pk).get_dict()
        # Average for pairs
        avg_energy = data['cleavage_energy_J_per_m2']
        cleavage_dict[data['term_i']] = avg_energy

    # Get surface energies at mid-range μ_O
    surface_energies = wg.outputs.surface_energies.get_dict()

    surface_dict = {}
    for term_label, energies_pk in surface_energies.items():
        energies = load_node(energies_pk).get_dict()
        # Take middle value
        mid_idx = len(energies['gamma']) // 2
        gamma_j_m2 = energies['gamma'][mid_idx] * 16.0217662  # eV/Ų → J/m²
        surface_dict[term_label] = gamma_j_m2

    # Plot comparison
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(surface_dict))
    width = 0.35

    ax.bar(x - width/2, list(surface_dict.values()), width,
           label='Surface Energy (mid μ_O)', color='coral', alpha=0.7)
    ax.bar(x + width/2, [cleavage_dict.get(t, 0) for t in surface_dict.keys()], width,
           label='Cleavage Energy', color='steelblue', alpha=0.7)

    ax.set_xlabel('Termination')
    ax.set_ylabel('Energy (J/m²)')
    ax.set_title('Surface vs. Cleavage Energies')
    ax.set_xticks(x)
    ax.set_xticklabels(surface_dict.keys())
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('surface_vs_cleavage.png', dpi=300)

Plot Band Structure (Basic)
----------------------------

.. code-block:: python

    from aiida.tools.data.array.bands import get_bands_data
    import matplotlib.pyplot as plt

    # Get band structure
    bands = wg.outputs.bulk_bands
    seekpath_params = wg.outputs.bulk_seekpath_parameters.get_dict()

    # Extract bands data
    bands_array, labels, k_coords = bands.get_bands()

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for band in bands_array.T:  # Transpose to iterate over bands
        ax.plot(k_coords, band, 'k-', linewidth=0.5)

    # Add high-symmetry points
    for i, label in enumerate(labels):
        if label:
            ax.axvline(x=k_coords[i], color='gray', linestyle='--', linewidth=0.5)
            ax.text(k_coords[i], ax.get_ylim()[1], label,
                   ha='center', va='bottom')

    ax.set_xlabel('k-path')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('Band Structure - Ag₂O Bulk')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(k_coords[0], k_coords[-1])

    plt.tight_layout()
    plt.savefig('band_structure.png', dpi=300)

.. note::
   For production-quality band structure plots, consider using specialized tools like `sumo <https://github.com/SMTG-UCL/sumo>`_ or the AiiDA-QuantumESPRESSO plotting utilities.

Understanding Output Patterns
==============================

**Relaxation Energy Patterns**:

* Polar terminations (charged): Large |E_relax| due to strong reconstruction
* Stoichiometric terminations: Smaller |E_relax|, less reconstruction
* Metal-rich terminations: Often larger relaxation than oxygen-rich

**Cleavage vs. Surface Energy**:

* **Cleavage energy**: Material property, independent of chemical potential
* **Surface energy**: Environment-dependent, varies with μ_O
* Cleavage energy ≈ 2× surface energy (for symmetric slabs at mid-range μ_O)

**Electronic Properties**:

* **Direct gap**: Valence band maximum and conduction band minimum at same k-point
* **Indirect gap**: VBM and CBM at different k-points
* **Ag₂O**: Indirect gap ~1.2 eV (experimental), often underestimated by PBE

Performance Notes
=================

**Computational Cost**:

* **Relaxation energy**: +30-50% time (unrelaxed SCF in parallel with relaxation)
* **Cleavage energy**: Negligible (post-processing only)
* **Electronic properties**: +50-100% time (band structure + DOS)

**Total time**: ~4-6 hours (vs. 2-4 hours for beginner)

**Memory Requirements**:

* Electronic properties require LWAVE=True and LCHARG=True for SCF
* Ensure sufficient disk space (~5-10 GB per calculation)

Troubleshooting
===============

Problem: Band structure looks wrong
------------------------------------

**Symptoms**: Flat bands, strange dispersion, no gap

**Solution**:

1. Check SCF converged properly (``verdi process report <SCF_PK>``)
2. Verify LWAVE=True and LCHARG=True in SCF stage
3. Check ICHARG=11 in band structure stage (non-self-consistent)
4. Ensure k-point path is reasonable (check seekpath output)

Problem: Relaxation energies are positive
------------------------------------------

**Symptoms**: E_relax > 0 (unrelaxed more stable than relaxed)

**Solution**:

1. This should NOT happen physically
2. Check that unrelaxed SCF used correct structure (from slab generation)
3. Verify relaxed calculation actually relaxed (NSW > 0, forces decreased)
4. Check for calculation failures in either unrelaxed or relaxed

Problem: No cleavage energies output
-------------------------------------

**Symptoms**: ``cleavage_energies`` output is empty or missing

**Solution**:

1. Check that pymatgen generated complementary terminations
2. Verify both terminations in pair were successfully relaxed
3. Some terminations may not have natural complements (unpaired terminations)

Next Steps
==========

You've now mastered intermediate PS-TEROS features:

✓ Relaxation and cleavage energy calculations
✓ Bulk electronic properties (DOS & bands)
✓ Advanced output interpretation
✓ Visualization techniques

**Ready for advanced features?**

* :doc:`Advanced Workflow <advanced-complete>` - Ternary oxides, slab electronic properties, complete production calculations
* :doc:`How-To Guide: Electronic Properties <../how-to/electronic-properties>` - Detailed configuration options
* :doc:`AIMD Workflow <aimd-molecular-dynamics>` - Finite-temperature molecular dynamics

**Need Help?**

* :doc:`API Reference </api/builders>` - Electronic properties builder documentation
* :doc:`Theory </theory>` - Mathematical background on surface/cleavage energies
* `GitHub Issues <https://github.com/your-repo/PS-TEROS/issues>`_ - Report bugs or ask questions

Complete Example Script
========================

.. code-block:: python

    #!/usr/bin/env python
    """
    Intermediate PS-TEROS workflow: Ag₂O with optional features enabled.

    Adds relaxation energy, cleavage energy, and bulk electronic properties
    to the beginner workflow.
    """

    import os
    from aiida import load_profile
    from teros.core.workgraph import build_core_workgraph
    from teros.core.builders import get_electronic_properties_defaults  # NEW

    def main():
        load_profile(profile='your_profile_name')

        script_dir = os.path.dirname(os.path.abspath(__file__))
        structures_dir = os.path.join(script_dir, 'structures')

        # ... (All beginner parameters: bulk_parameters, metal_parameters, etc.) ...

        # NEW: Electronic properties configuration
        ep_defaults = get_electronic_properties_defaults(
            energy_cutoff=520,
            electronic_convergence=1e-5,
            ncore=4,
            ispin=2,
            kpoints_mesh_density=0.3,
            band_kpoints_distance=0.2,
            dos_kpoints_distance=0.2,
            line_density=0.2,
            nedos=2000,
            band_mode="seekpath-aiida",
        )

        # Build workflow with new features enabled
        wg = build_core_workgraph(
            # ... all beginner parameters ...

            # NEW: Enable all intermediate features
            compute_relaxation_energy=True,
            compute_cleavage=True,
            compute_thermodynamics=True,
            thermodynamics_sampling=10,

            # NEW: Bulk electronic properties
            compute_electronic_properties_bulk=True,
            bands_parameters=ep_defaults,
            band_settings=ep_defaults['band_settings'],
            bands_options=bulk_options,

            name='Ag2O_Intermediate_Workflow',
        )

        wg.submit(wait=False)
        print(f"✓ Workflow submitted! PK: {wg.pk}")

        return wg

    if __name__ == '__main__':
        main()

.. tip::
   See ``examples/complete/complete_ag2o_example.py`` for the complete working script with all parameters.
