===============================
Advanced: Complete Workflow
===============================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

This workflow demonstrates the complete capabilities of PS-TEROS for production-quality calculations on ternary oxide systems. All features are enabled, providing publication-ready data.

**Material**: Ag₃PO₄ (silver phosphate, body-centered cubic)

**What's Included**:

1. All features from :doc:`intermediate workflow <intermediate-with-features>`
2. **Ternary oxide support**: Additional nonmetal reference (phosphorus)
3. **Slab electronic properties**: DOS and band structure for selected slab terminations
4. **Per-slab configuration**: Independent parameters for heterogeneous terminations
5. **High-resolution thermodynamics**: Denser chemical potential sampling

**Prerequisites**: Understand :doc:`beginner <beginner-surface-energy>` and :doc:`intermediate <intermediate-with-features>` workflows

**Runtime**: ~8-12 hours on typical cluster (40 cores)

Why Ternary Oxides?
====================

**Complexity**: Ternary oxides (AxByOz) have more degrees of freedom than binary oxides:

* Two chemical potentials: μ_B and μ_O (μ_A determined by equilibrium)
* More diverse terminations (metal-rich, nonmetal-rich, oxygen-rich)
* Richer surface chemistry and catalytic potential

**Example: Ag₃PO₄**:

.. math::

   \Delta H_f = E_{Ag_3PO_4} - 3 \mu_{Ag} - \mu_P - 4 \mu_O

Surface energy depends on both μ_P and μ_O, creating a 2D phase diagram.

Slab Electronic Properties
===========================

**Why calculate for slabs?**

* Band structure shows surface state emergence
* Gap narrowing/widening at surfaces compared to bulk
* Orbital contributions from surface vs. bulk atoms
* Identify metallic surface states on insulating bulk

**Considerations**:

* More expensive than bulk (larger systems, lower symmetry)
* Requires denser k-point sampling in plane parallel to surface
* Not all terminations need electronic properties (select key ones)

New Code Elements
=================

1. Nonmetal Reference
---------------------

Add phosphorus parameters alongside metal and oxygen:

.. code-block:: python

    # Phosphorus reference parameters
    nonmetal_filename = 'P.cif'  # Black phosphorus or other polymorph

    nonmetal_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,        # Semiconductor smearing
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,          # Full relaxation for bulk P
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    nonmetal_options = bulk_options.copy()

    nonmetal_potential_mapping = {
        'P': 'P',
    }

Include in ``build_core_workgraph()``:

.. code-block:: python

    wg = build_core_workgraph(
        # ... existing parameters ...

        # Nonmetal reference (NEW for ternary)
        nonmetal_name=nonmetal_filename,
        nonmetal_potential_mapping=nonmetal_potential_mapping,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,

        # ... rest of parameters ...
    )

2. Slab Electronic Properties Configuration
--------------------------------------------

Import the slab electronic properties builder:

.. code-block:: python

    from teros.core.builders import (
        get_electronic_properties_defaults,
        get_slab_electronic_properties_defaults,  # NEW
    )

Configure with denser k-points for 2D systems:

.. code-block:: python

    # Slab electronic properties (denser than bulk)
    slab_ep_defaults = get_slab_electronic_properties_defaults(
        energy_cutoff=slab_parameters['ENCUT'],
        electronic_convergence=1e-5,
        ncore=4,
        ispin=2,
        lasph=True,
        lreal="Auto",
        kpoints_mesh_density=0.25,    # Denser than bulk (0.3)
        band_kpoints_distance=0.15,   # Denser path sampling
        dos_kpoints_distance=0.2,
        line_density=0.15,            # More points along paths
        nedos=2000,
        sigma_bands=0.01,
        symprec=1e-4,
        band_mode="seekpath-aiida",
    )

3. Per-Slab Configuration
-------------------------

Define which slabs to calculate and customize per-slab:

.. code-block:: python

    # Enable slab electronic properties
    compute_electronic_properties_slabs = True

    # Configure specific terminations
    slab_electronic_properties = {
        'term_0': {  # Example: Ag-rich termination
            'bands_parameters': slab_ep_defaults,
            'bands_options': {
                'resources': {
                    'num_machines': 1,
                    'num_cores_per_machine': 40,
                },
                'queue_name': 'par40',
            },
            'band_settings': slab_ep_defaults['band_settings'],
        },
        'term_1': {  # Example: PO₄-rich termination
            'bands_parameters': slab_ep_defaults,
            'bands_options': {
                'resources': {
                    'num_machines': 1,
                    'num_cores_per_machine': 40,
                },
                'queue_name': 'par40',
            },
            'band_settings': slab_ep_defaults['band_settings'],
        },
        # Add more terminations as needed
        # Can customize parameters per termination:
        # 'term_2': {
        #     'bands_parameters': custom_params_for_term2,
        #     'bands_options': higher_memory_options,
        #     'band_settings': custom_band_settings,
        # },
    }

Include in ``build_core_workgraph()``:

.. code-block:: python

    wg = build_core_workgraph(
        # ... all previous parameters ...

        # Slab electronic properties (NEW)
        compute_electronic_properties_slabs=compute_electronic_properties_slabs,
        slab_electronic_properties=slab_electronic_properties,
        slab_bands_parameters=slab_ep_defaults,
        slab_band_settings=slab_ep_defaults['band_settings'],
        slab_bands_options=slab_options,

        name='Ag3PO4_Advanced_Complete',
    )

Complete Workflow Structure
============================

The advanced workflow includes all steps:

.. code-block:: text

    ┌──────────────────────────────────────┐
    │ Parallel Bulk/Reference Relaxation:  │
    │  - Bulk (Ag₃PO₄)                    │
    │  - Metal (Ag)                        │
    │  - Nonmetal (P)  ← NEW              │
    │  - Oxygen (O₂)                       │
    └──────────────────────────────────────┘
             ↓
    Formation Enthalpy (ternary formula)
             ↓
    Bulk Electronic Properties (DOS + Bands)
             ↓
    Generate Slab Terminations (pymatgen)
             ↓
    ┌──────────────────────────────────────┐
    │ Parallel Slab Calculations:          │
    │  - Unrelaxed SCF (all terms)         │
    │  - Full Relaxation (all terms)       │
    └──────────────────────────────────────┘
             ↓
    Relaxation Energies (E_relax)
             ↓
    Cleavage Energies (complementary pairs)
             ↓
    ┌──────────────────────────────────────┐
    │ Surface Thermodynamics:              │
    │  - 2D sampling (μ_P, μ_O)           │
    │  - γ(μ_P, μ_O) for each termination │
    └──────────────────────────────────────┘
             ↓
    ┌──────────────────────────────────────┐
    │ Slab Electronic Properties:  ← NEW   │
    │  For selected terminations:          │
    │   - SCF (LWAVE=True, LCHARG=True)   │
    │   - Band structure (non-SCF)         │
    │   - DOS (non-SCF)                    │
    └──────────────────────────────────────┘

New Outputs
===========

Nonmetal Energy
---------------

.. code-block:: python

    from aiida import load_profile
    from aiida.orm import load_node

    load_profile()
    wg = load_node(<YOUR_PK>)

    # Access nonmetal energy
    nonmetal_energy = wg.outputs.nonmetal_energy.value  # eV

    print(f"Bulk: {wg.outputs.bulk_energy.value:.3f} eV")
    print(f"Metal (Ag): {wg.outputs.metal_energy.value:.3f} eV")
    print(f"Nonmetal (P): {nonmetal_energy:.3f} eV")
    print(f"Oxygen (O₂): {wg.outputs.oxygen_energy.value:.3f} eV")

Formation enthalpy uses all four energies for ternary formula.

Slab Electronic Properties
---------------------------

.. code-block:: python

    # Access slab band structures
    slab_bands = wg.outputs.slab_bands.get_dict()

    print("Slab band structures available:")
    for term_label, bands_pk in slab_bands.items():
        bands = load_node(bands_pk)
        print(f"  {term_label}: {bands.pk}")

    # Access slab DOS
    slab_dos = wg.outputs.slab_dos.get_dict()

    # Access slab primitive structures (used for bands)
    slab_primitive_structures = wg.outputs.slab_primitive_structures.get_dict()

    # Access slab seekpath parameters
    slab_seekpath_parameters = wg.outputs.slab_seekpath_parameters.get_dict()

    for term_label, seekpath_pk in slab_seekpath_parameters.items():
        seekpath = load_node(seekpath_pk).get_dict()
        print(f"\n{term_label} high-symmetry path:")
        print(f"  Path: {' → '.join(seekpath['path'])}")
        print(f"  Points: {list(seekpath['point_coords'].keys())}")

2D Surface Energy Phase Diagram
--------------------------------

For ternary oxides, surface energy depends on two chemical potentials:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # Get surface energies for one termination
    surface_energies = wg.outputs.surface_energies.get_dict()

    term_label = 'term_0'  # Choose termination
    energies = load_node(surface_energies[term_label]).get_dict()

    # Extract 2D grid
    mu_P = np.array(energies['mu_nonmetal'])  # or 'mu_P'
    mu_O = np.array(energies['mu_O'])
    gamma = np.array(energies['gamma'])

    # Reshape to 2D grid (assuming NxN sampling)
    N = int(np.sqrt(len(gamma)))
    mu_P_grid = mu_P.reshape((N, N))
    mu_O_grid = mu_O.reshape((N, N))
    gamma_grid = gamma.reshape((N, N))

    # Contour plot
    fig, ax = plt.subplots(figsize=(10, 8))
    contour = ax.contourf(mu_P_grid, mu_O_grid, gamma_grid,
                          levels=20, cmap='viridis')
    ax.contour(mu_P_grid, mu_O_grid, gamma_grid,
               levels=10, colors='black', linewidths=0.5, alpha=0.3)

    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label('Surface Energy γ (eV/Ų)', rotation=270, labelpad=20)

    ax.set_xlabel('μ_P (eV)')
    ax.set_ylabel('μ_O (eV)')
    ax.set_title(f'Surface Energy Phase Diagram - {term_label}')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'phase_diagram_{term_label}.png', dpi=300)

Compare Bulk vs. Slab Band Structures
--------------------------------------

.. code-block:: python

    import matplotlib.pyplot as plt

    # Get bulk bands
    bulk_bands = wg.outputs.bulk_bands
    bulk_array, bulk_labels, bulk_k = bulk_bands.get_bands()

    # Get slab bands (example: term_0)
    slab_bands_dict = wg.outputs.slab_bands.get_dict()
    slab_bands = load_node(slab_bands_dict['term_0'])
    slab_array, slab_labels, slab_k = slab_bands.get_bands()

    # Plot comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Bulk
    for band in bulk_array.T:
        ax1.plot(bulk_k, band, 'k-', linewidth=0.5)
    ax1.set_xlabel('k-path')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title('Bulk Ag₃PO₄')
    ax1.grid(True, alpha=0.3)

    # Slab
    for band in slab_array.T:
        ax2.plot(slab_k, band, 'b-', linewidth=0.5)
    ax2.set_xlabel('k-path')
    ax2.set_ylabel('Energy (eV)')
    ax2.set_title('Slab term_0')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('bulk_vs_slab_bands.png', dpi=300)

**Look for**:

* Surface states appearing in the gap
* Band edge shifts (band bending at surface)
* Metallization of surface (gap closure)

Production Best Practices
==========================

Parameter Convergence Testing
------------------------------

Before production runs, test convergence:

**ENCUT convergence**:

.. code-block:: python

    encuts = [400, 450, 500, 520, 550, 600]
    # Run small test calculations at each ENCUT
    # Plot formation enthalpy vs. ENCUT
    # Choose ENCUT where ΔH_f changes < 0.01 eV

**K-points convergence**:

.. code-block:: python

    kpoints_spacings = [0.5, 0.4, 0.35, 0.3, 0.25]
    # Test on bulk and smallest slab
    # Choose spacing where energies converge to < 0.01 eV

**Slab thickness convergence**:

.. code-block:: python

    slab_thicknesses = [10, 15, 20, 25]  # Angstroms
    # Compare surface energies
    # Choose thickness where γ converges to < 0.01 eV/Ų

Data Management
---------------

**Large outputs**: Advanced workflows generate significant data:

* Band structure files: ~100-500 MB per structure
* Trajectory data (if saving): ~1-5 GB per relaxation
* Total: ~10-50 GB per complete workflow

**Organization**:

.. code-block:: bash

    project/
    ├── calculations/
    │   ├── bulk/
    │   ├── references/
    │   ├── slabs/
    │   └── electronic_properties/
    ├── results/
    │   ├── formation_enthalpy.dat
    │   ├── surface_energies/
    │   ├── band_structures/
    │   └── figures/
    └── analysis/
        └── notebooks/

**Archiving**:

.. code-block:: bash

    # Export AiiDA data
    verdi export create workflow_<PK>.aiida -N <PK>

    # Compress band structure data
    tar -czf band_structures.tar.gz results/band_structures/

Validation Checklist
====================

Before publishing results, verify:

**✓ Energy Convergence**:

* All bulk/reference relaxations converged (check NSW, forces)
* Electronic SCF converged (check NELM, EDIFF)
* Slab relaxations converged (may take NSW=200-300)

**✓ Physical Sanity**:

* Formation enthalpy negative (for stable compounds)
* Surface energies positive
* Relaxation energies negative
* Band gap reasonable (compare with literature)

**✓ Consistency**:

* Same ENCUT across all calculations
* Consistent pseudopotentials (same family)
* K-point density appropriate for system size

**✓ Thermodynamic Limits**:

* Chemical potential ranges respect stability limits
* No unphysical regions in phase diagram

Common Advanced Issues
======================

Issue: Slab band structure fails
---------------------------------

**Symptoms**: SCF succeeds but band calculation fails

**Solution**:

1. Check LWAVE=True and LCHARG=True in SCF
2. Verify SCF converged fully (not stopped at NELM)
3. Try reducing ``band_kpoints_distance`` (denser path)
4. Check symmetry detection (``symprec`` parameter)
5. Ensure primitive cell detection worked (check ``slab_primitive_structures``)

Issue: Surface energy phase diagram has gaps
---------------------------------------------

**Symptoms**: Missing data points in 2D grid

**Solution**:

1. Some (μ_P, μ_O) combinations may be outside stability region
2. Increase ``thermodynamics_sampling`` for denser grid
3. Check that all slabs successfully relaxed
4. Verify formation enthalpy is correct

Issue: Different terminations show same bands
----------------------------------------------

**Symptoms**: term_0 and term_1 have identical band structures

**Solution**:

1. Check that slabs are actually different (inspect structures)
2. Verify correct slab was used for each calculation
3. Differences may be subtle - compare DOS carefully
4. Some terminations may be nearly equivalent by symmetry

Computational Requirements
==========================

**Minimum**:

* 40 cores per node
* 4-8 GB RAM per core
* 100 GB scratch space
* ~10-20 hours walltime

**Recommended**:

* 80-120 cores total (2-3 nodes)
* 8-16 GB RAM per core
* 500 GB scratch space
* ~12-24 hours walltime

**Parallel Efficiency**:

* Bulk/reference: Good scaling to 40-80 cores
* Slabs: May benefit from 80-120 cores (depends on size)
* Electronic properties: Band calculation serial, DOS parallel

Next Steps
==========

You've completed the full PS-TEROS advanced workflow:

✓ Ternary oxide calculations
✓ Complete electronic structure analysis
✓ Production-quality data generation
✓ Publication-ready phase diagrams

**What's Next?**

* :doc:`AIMD Workflow <aimd-molecular-dynamics>` - Add finite-temperature dynamics
* :doc:`How-To: Restart Calculations <../how-to/restart-calculations>` - Handle failures gracefully
* :doc:`API Reference </api/index>` - Deep dive into parameters

**Publishing Results**

When preparing manuscripts:

1. Report all convergence tests (ENCUT, k-points, slab thickness)
2. Include formation enthalpy comparison with experiment/database
3. Show surface energy phase diagrams with stability regions marked
4. Compare band structures/DOS with experimental data
5. Cite PS-TEROS and AiiDA properly

**Example Citation**:

    Calculations were performed using PS-TEROS (version X.Y.Z), built on
    AiiDA [cite] and AiiDA-WorkGraph [cite], using VASP [cite] as the
    DFT engine with PBE functionals [cite] and PAW pseudopotentials [cite].

Complete Example
================

See ``examples/complete/complete_ag3po4_example.py`` for the full working script with all parameters and detailed comments.

.. tip::
   The complete example includes helpful print statements showing workflow progress, expected outputs, and monitoring commands. Use it as a template for your own materials.
