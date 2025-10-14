==========================================
How to Configure Electronic Properties
==========================================

This guide shows you how to configure and compute electronic properties (DOS and band structure) for bulk and slab structures in PS-TEROS.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

PS-TEROS can compute electronic properties alongside structural calculations:

* **Bulk electronic properties** - DOS and band structure of bulk oxide
* **Slab electronic properties** - DOS and band structure for each surface termination

**Key concepts:**

* **DOS (Density of States)** - Distribution of electronic states vs energy
* **Band structure** - E(k) along high-symmetry paths in Brillouin zone
* **Seekpath integration** - Automatic k-point path selection
* **Builder functions** - `get_electronic_properties_defaults()` for easy setup
* **Per-slab configuration** - Different DOS/bands settings for each termination

When to Compute Electronic Properties
======================================

Compute bulk electronic properties when:

* Characterizing oxide semiconductors (band gap, VBM/CBM positions)
* Comparing with experimental optical/transport data
* Understanding bulk electronic structure before surface analysis

Compute slab electronic properties when:

* Studying surface states and band bending
* Analyzing work functions
* Comparing termination-dependent electronic structure
* Understanding surface reactivity via DOS at Fermi level

**Cost considerations:**

* DOS calculation: ~2-3x SCF cost (needs dense k-mesh)
* Band structure: ~1x SCF cost per path segment
* Total: Electronic properties add ~30-50% to workflow time

Quick Start: Bulk Electronic Properties
========================================

Using Builder Function (Recommended)
-------------------------------------

Easiest way - use the builder with sensible defaults:

.. code-block:: python

   from teros.core.workgraph import build_core_workgraph
   from teros.core.builders import get_electronic_properties_defaults

   # Get default electronic properties configuration
   ep_defaults = get_electronic_properties_defaults(
       energy_cutoff=520,  # Match your bulk ENCUT
       kpoints_mesh_density=0.3,  # Dense k-mesh for DOS
       band_kpoints_distance=0.2,  # Band structure k-spacing
       nedos=2000,  # DOS energy resolution
       band_mode="seekpath-aiida",  # Automatic path selection
   )

   # Build workflow with electronic properties
   wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       code_label='vasp@cluster',
       potential_family='PBE.54',

       # Standard bulk/slab parameters
       bulk_parameters=bulk_parameters,
       slab_parameters=slab_parameters,
       bulk_options=bulk_options,
       slab_options=slab_options,

       # Enable electronic properties
       compute_electronic_properties_bulk=True,
       bands_parameters=ep_defaults,
       band_settings=ep_defaults['band_settings'],
       bands_options=bulk_options,  # Use same resources as bulk

       # Other features
       compute_thermodynamics=True,
       thermodynamics_sampling=10,
   )

   wg.submit(wait=True)

**What the builder provides:**

* DOS parameters (NEDOS, ISMEAR, SIGMA)
* Band structure k-point path (automatic via seekpath)
* Sensible defaults based on your ENCUT

Manual Configuration
--------------------

For full control, configure parameters manually:

.. code-block:: python

   from aiida import orm

   # DOS parameters
   bands_parameters = {
       'kpoints_mesh_density': 0.3,  # Dense mesh (smaller = denser)
       'nedos': 2000,  # Energy points
       'ismear_dos': 0,  # Gaussian smearing for DOS
       'sigma_dos': 0.05,  # Narrow smearing (eV)
       'lorbit': 11,  # Project onto atoms
   }

   # Band structure settings
   band_settings = orm.Dict(dict={
       'band_mode': 'seekpath-aiida',  # Auto high-symmetry path
       'band_kpoints_distance': 0.2,  # k-spacing along path (Å⁻¹)
       'ismear_bands': 0,  # Tetrahedron for accurate eigenvalues
       'sigma_bands': 0.01,  # Small smearing
   })

   # Use in workflow
   wg = build_core_workgraph(
       # ... structure parameters ...
       compute_electronic_properties_bulk=True,
       bands_parameters=bands_parameters,
       band_settings=band_settings,
       bands_options=bulk_options,
   )

Quick Start: Slab Electronic Properties
========================================

Per-Slab Configuration
----------------------

Configure electronic properties for specific slabs:

.. code-block:: python

   from teros.core.builders import get_slab_electronic_properties_defaults

   # Get slab-specific defaults
   slab_ep_defaults = get_slab_electronic_properties_defaults(
       energy_cutoff=520,
       kpoints_mesh_density=0.25,  # Slightly sparser for large slabs
       band_kpoints_distance=0.2,
       nedos=2000,
       band_mode="seekpath-aiida",
   )

   # Options for slab electronic properties (more resources)
   slab_ep_options = {
       'queue_name': 'normal',
       'resources': {'num_machines': 2, 'num_mpiprocs_per_machine': 48},
       'max_wallclock_seconds': 7200,  # 2 hours
   }

   # Configure per termination
   slab_electronic_properties = {
       'ag2o_111_term1': {
           'bands_parameters': slab_ep_defaults,
           'bands_options': slab_ep_options,
           'band_settings': slab_ep_defaults['band_settings'],
       },
       'ag2o_111_term2': {
           'bands_parameters': slab_ep_defaults,
           'bands_options': slab_ep_options,
           'band_settings': slab_ep_defaults['band_settings'],
       },
   }

   # Enable in workflow
   wg = build_core_workgraph(
       # ... other parameters ...
       compute_electronic_properties_slab=True,
       slab_electronic_properties=slab_electronic_properties,
   )

**Why per-slab configuration:**

* Different terminations may need different k-meshes
* Can enable only for interesting terminations (saves compute)
* Adjust DOS resolution based on slab complexity

All Slabs with Same Settings
-----------------------------

To apply same configuration to all slabs:

.. code-block:: python

   # Get default settings
   slab_ep_defaults = get_slab_electronic_properties_defaults(
       energy_cutoff=520,
       kpoints_mesh_density=0.25,
       band_kpoints_distance=0.2,
       nedos=2000,
   )

   # Create configuration for all slabs
   slab_names = ['ag2o_111_term1', 'ag2o_111_term2', 'ag2o_100_term1']
   slab_electronic_properties = {
       slab_name: {
           'bands_parameters': slab_ep_defaults,
           'bands_options': slab_ep_options,
           'band_settings': slab_ep_defaults['band_settings'],
       }
       for slab_name in slab_names
   }

   # Or use a helper pattern
   def apply_to_all_slabs(slab_list, ep_config):
       return {slab: ep_config for slab in slab_list}

Parameter Reference
===================

K-Point Mesh Density (DOS)
---------------------------

Controls DOS k-point mesh via density parameter:

.. code-block:: python

   bands_parameters = {
       'kpoints_mesh_density': 0.3,  # Standard
   }

**Typical values:**

* 0.5 - Coarse (quick test, ~10³ k-points)
* 0.3 - Standard (production, ~10⁴ k-points)
* 0.2 - Dense (high accuracy, ~10⁵ k-points)
* 0.1 - Very dense (converged DOS, expensive)

**Relationship to k-mesh:**

Density = 2π / (k-spacing in Å⁻¹)

For 0.3 density on cubic cell (a=4Å): ~(2π/0.3)/4 ≈ 5 k-points per direction

**Convergence check:**

.. code-block:: python

   densities = [0.5, 0.4, 0.3, 0.25, 0.2]
   for density in densities:
       ep_config = get_electronic_properties_defaults(
           energy_cutoff=520,
           kpoints_mesh_density=density,
       )
       # Run and compare band gaps

When DOS shape stops changing, you're converged.

Band Structure K-Point Spacing
-------------------------------

Controls path resolution:

.. code-block:: python

   band_settings = {
       'band_kpoints_distance': 0.2,  # Å⁻¹ spacing
   }

**Typical values:**

* 0.3 - Coarse (fast, may miss features)
* 0.2 - Standard (good resolution)
* 0.1 - Dense (publication-quality)
* 0.05 - Very dense (converged dispersion)

**Practical guidance:**

* Start with 0.2
* Refine to 0.1 if bands look jagged
* Use 0.05 only for final publication figures

DOS Energy Resolution (NEDOS)
------------------------------

Number of energy grid points for DOS:

.. code-block:: python

   bands_parameters = {
       'nedos': 2000,  # Standard
   }

**Typical values:**

* 1000 - Coarse (quick preview)
* 2000 - Standard (most analyses)
* 5000 - High resolution (narrow peaks, semiconductors)
* 10000 - Very high (ultra-sharp features)

**Trade-off:**

Higher NEDOS = better resolution but larger output files

Smearing for DOS vs Bands
--------------------------

Different smearing for DOS and bands:

.. code-block:: python

   bands_parameters = {
       'ismear_dos': 0,     # Gaussian for DOS
       'sigma_dos': 0.05,   # 0.05 eV smearing
   }

   band_settings = {
       'ismear_bands': 0,      # Tetrahedron or Gaussian
       'sigma_bands': 0.01,    # Narrow smearing
   }

**Guidelines:**

* **DOS**: ISMEAR=0 (Gaussian), SIGMA=0.05-0.1 eV
* **Bands**: ISMEAR=0 (Gaussian), SIGMA=0.01 eV or ISMEAR=-5 (tetrahedron)
* **Metals**: ISMEAR=1 or 2 (Methfessel-Paxton)
* **Semiconductors**: ISMEAR=0 (Gaussian)

Band Mode Selection
-------------------

How to generate k-point path:

.. code-block:: python

   band_settings = {
       'band_mode': 'seekpath-aiida',  # Recommended
   }

**Options:**

* ``seekpath-aiida`` - Automatic high-symmetry path (works for any structure)
* ``manual`` - Provide explicit k-point list (advanced users)

**Seekpath advantages:**

* Automatically identifies Brillouin zone symmetry
* Generates standardized paths (Γ-X-M-Γ, etc.)
* Works for complex structures without manual effort
* Consistent with literature conventions

Advanced: Manual K-Point Path
------------------------------

For custom paths (rare):

.. code-block:: python

   from aiida import orm

   # Define explicit path
   kpoint_path = orm.Dict(dict={
       'path': [
           {'from': 'GAMMA', 'to': 'X', 'num_points': 20},
           {'from': 'X', 'to': 'M', 'num_points': 15},
           {'from': 'M', 'to': 'GAMMA', 'num_points': 20},
       ],
       'labels': {'GAMMA': [0,0,0], 'X': [0.5,0,0], 'M': [0.5,0.5,0]},
   })

   band_settings = {
       'band_mode': 'manual',
       'kpoint_path': kpoint_path,
   }

Analyzing Results
=================

Bulk Electronic Properties
--------------------------

Access DOS and band structure data:

.. code-block:: python

   from aiida import orm

   # Load workflow
   wg_node = orm.load_node(<WG_PK>)

   # Get DOS calculation
   dos_calc = wg_node.outputs.dos_bulk
   dos_data = dos_calc.outputs.dos

   # Get band structure calculation
   bands_calc = wg_node.outputs.bands_bulk
   bands_data = bands_calc.outputs.bands

**Extract DOS:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   # Get DOS data
   energies = dos_data.get_x()[1]  # Energy grid (eV)
   total_dos = dos_data.get_y()[0][1]  # Total DOS

   # Get Fermi energy
   fermi = dos_data.get_fermi_energy()

   # Plot
   plt.figure(figsize=(8, 6))
   plt.plot(energies - fermi, total_dos)
   plt.axvline(0, color='k', linestyle='--', label='Fermi level')
   plt.xlabel('E - E_F (eV)')
   plt.ylabel('DOS (states/eV)')
   plt.xlim(-10, 10)
   plt.legend()
   plt.savefig('bulk_dos.png')

**Extract band structure:**

.. code-block:: python

   # Get band structure
   bands = bands_data.get_bands()  # Shape: (n_bands, n_kpoints)
   kpoints = bands_data.get_kpoints()  # k-point coordinates

   # Plot
   plt.figure(figsize=(8, 6))
   for band in bands:
       plt.plot(band - fermi, 'b-', linewidth=0.5)

   plt.axhline(0, color='k', linestyle='--', label='Fermi level')
   plt.ylabel('E - E_F (eV)')
   plt.xlabel('k-path')
   plt.ylim(-10, 10)
   plt.legend()
   plt.savefig('bulk_bands.png')

**Extract band gap:**

.. code-block:: python

   # Get band gap from DOS calculation
   dos_calc_output = dos_calc.outputs.output_parameters

   if 'band_gap' in dos_calc_output:
       band_gap = dos_calc_output['band_gap']
       print(f"Band gap: {band_gap:.3f} eV")
   else:
       print("Metallic system (no band gap)")

Slab Electronic Properties
---------------------------

Access per-slab DOS and bands:

.. code-block:: python

   # List all slab DOS calculations
   for key in wg_node.outputs:
       if 'dos_slab_' in key:
           slab_name = key.replace('dos_slab_', '')
           dos_calc = wg_node.outputs[key]
           print(f"Slab: {slab_name}")

           # Analyze DOS
           dos_data = dos_calc.outputs.dos
           fermi = dos_data.get_fermi_energy()
           print(f"  Fermi energy: {fermi:.3f} eV")

**Compare terminations:**

.. code-block:: python

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(2, 1, figsize=(8, 10))

   # Plot DOS for different terminations
   for term_name in ['ag2o_111_term1', 'ag2o_111_term2']:
       dos_calc = wg_node.outputs[f'dos_slab_{term_name}']
       dos_data = dos_calc.outputs.dos

       energies = dos_data.get_x()[1]
       total_dos = dos_data.get_y()[0][1]
       fermi = dos_data.get_fermi_energy()

       axes[0].plot(energies - fermi, total_dos, label=term_name)

   axes[0].set_xlabel('E - E_F (eV)')
   axes[0].set_ylabel('DOS')
   axes[0].legend()
   axes[0].set_title('DOS Comparison')

   # Similarly for bands
   for term_name in ['ag2o_111_term1', 'ag2o_111_term2']:
       bands_calc = wg_node.outputs[f'bands_slab_{term_name}']
       bands_data = bands_calc.outputs.bands

       bands = bands_data.get_bands()
       fermi = bands_data.get_fermi_energy()

       # Plot first few bands
       for band in bands[:5]:
           axes[1].plot(band - fermi, linewidth=0.5)

   axes[1].set_ylabel('E - E_F (eV)')
   axes[1].set_title('Band Structure Comparison')

   plt.tight_layout()
   plt.savefig('slab_comparison.png')

Common Issues and Solutions
============================

Issue 1: DOS Calculation Runs Forever
--------------------------------------

**Symptoms:**

* DOS calculation doesn't finish within expected time
* High CPU usage but slow progress

**Cause:**

K-mesh too dense for available resources.

**Solution:**

.. code-block:: python

   # Reduce k-mesh density
   bands_parameters = {
       'kpoints_mesh_density': 0.4,  # Was 0.2
       'nedos': 1000,  # Was 2000
   }

   # Or add more computational resources
   bands_options = {
       'resources': {'num_machines': 2, 'num_mpiprocs_per_machine': 48},
       'max_wallclock_seconds': 7200,  # 2 hours
   }

Issue 2: Band Structure Looks Wrong
------------------------------------

**Symptoms:**

* Discontinuous bands
* Missing high-symmetry points
* Unexpected band crossings

**Possible causes:**

1. **k-spacing too coarse:**

.. code-block:: python

   # Increase resolution
   band_settings = {
       'band_kpoints_distance': 0.1,  # Was 0.3
   }

2. **Wrong band mode:**

.. code-block:: python

   # Use seekpath instead of manual
   band_settings = {
       'band_mode': 'seekpath-aiida',  # Was 'manual'
   }

3. **Electronic convergence issues:**

Check if SCF converged:

.. code-block:: bash

   verdi calcjob outputcat <BANDS_PK> OUTCAR | grep "aborting loop"

If not converged, increase NELM or adjust electronic parameters.

Issue 3: Band Gap Incorrect
----------------------------

**Symptoms:**

* Band gap differs from expected/literature values
* Metal predicted for known semiconductor

**Causes:**

1. **DFT functional limitation (PBE underestimates gaps):**

Expected - compare with PBE literature values, not experimental.

2. **k-mesh not converged:**

.. code-block:: python

   # Increase k-mesh density
   bands_parameters = {
       'kpoints_mesh_density': 0.2,  # Was 0.5
   }

3. **Smearing too large:**

.. code-block:: python

   # Reduce smearing
   bands_parameters = {
       'ismear_dos': 0,
       'sigma_dos': 0.01,  # Was 0.1
   }

4. **Spin polarization needed:**

For magnetic systems:

.. code-block:: python

   bulk_parameters = {
       # ... other parameters ...
       'ISPIN': 2,  # Enable spin polarization
       'MAGMOM': '2.0 2.0 0.0',  # Initial magnetic moments
   }

Issue 4: Slab Electronic Properties Fail
-----------------------------------------

**Symptoms:**

* Slab electronic properties calculations fail while bulk succeeds
* Out of memory errors

**Cause:**

Slabs are much larger than bulk, need more resources.

**Solution:**

.. code-block:: python

   # Significantly increase resources for slabs
   slab_ep_options = {
       'queue_name': 'large',
       'resources': {'num_machines': 4, 'num_mpiprocs_per_machine': 48},
       'max_wallclock_seconds': 14400,  # 4 hours
   }

   # Reduce k-mesh for slabs
   slab_ep_defaults = get_slab_electronic_properties_defaults(
       energy_cutoff=520,
       kpoints_mesh_density=0.4,  # Coarser than bulk
       band_kpoints_distance=0.25,
       nedos=1000,  # Fewer points
   )

Best Practices
==============

Convergence Testing
-------------------

Always test k-mesh convergence:

.. code-block:: python

   # Test sequence
   test_densities = [0.5, 0.4, 0.3, 0.25, 0.2]

   for density in test_densities:
       ep_config = get_electronic_properties_defaults(
           energy_cutoff=520,
           kpoints_mesh_density=density,
       )

       # Run single DOS calculation
       # Compare band gap and DOS shape

   # Use converged density in production

**Convergence criteria:**

* Band gap changes <0.01 eV between consecutive densities
* DOS shape visually similar

Computational Efficiency
------------------------

**1. Start coarse, refine if needed:**

.. code-block:: python

   # Phase 1: Quick check
   quick_ep = get_electronic_properties_defaults(
       energy_cutoff=400,  # Lower than production
       kpoints_mesh_density=0.5,
       nedos=1000,
   )

   # Phase 2: Production (if Phase 1 looks good)
   production_ep = get_electronic_properties_defaults(
       energy_cutoff=520,
       kpoints_mesh_density=0.3,
       nedos=2000,
   )

**2. Selective slab calculations:**

Don't compute electronic properties for all slabs - only interesting ones:

.. code-block:: python

   # Only compute for most stable termination
   slab_electronic_properties = {
       'ag2o_111_term1': {  # Most stable
           'bands_parameters': slab_ep_defaults,
           'bands_options': slab_ep_options,
           'band_settings': slab_ep_defaults['band_settings'],
       },
       # Skip other terminations
   }

**3. Reuse bulk electronic properties:**

If bulk doesn't change, reuse previous calculation:

.. code-block:: python

   # Load previous DOS/bands from earlier workflow
   previous_wg = orm.load_node(<PREVIOUS_WG_PK>)
   previous_dos = previous_wg.outputs.dos_bulk
   previous_bands = previous_wg.outputs.bands_bulk

   # Use in new workflow (implementation depends on caching)

Documentation
-------------

Document your electronic properties setup:

.. code-block:: text

   # Electronic Properties Configuration

   ## DOS Parameters
   - K-mesh density: 0.3 (converged to 0.01 eV in band gap)
   - NEDOS: 2000
   - Smearing: Gaussian, σ=0.05 eV

   ## Band Structure
   - K-path: Seekpath automatic
   - K-spacing: 0.2 Å⁻¹
   - Segments: Γ-X-M-Γ-R

   ## Results Summary
   - Bulk band gap: 1.45 eV (PBE, indirect)
   - Slab work function: 4.8 eV (term1)

Related Guides
==============

* :doc:`../workflows/intermediate-with-features` - Using bulk electronic properties in workflows
* :doc:`../workflows/advanced-complete` - Slab electronic properties example
* :doc:`custom-slabs` - Electronic properties of custom surfaces
* :doc:`restart-calculations` - Troubleshooting failed DOS/bands calculations

External Resources
==================

* `VASP DOS calculation <https://www.vasp.at/wiki/index.php/DOSCAR>`_
* `VASP band structure <https://www.vasp.at/wiki/index.php/Band-structure_calculation>`_
* `Seekpath documentation <https://github.com/giovannipizzi/seekpath>`_
* `K-point convergence guide <https://www.vasp.at/wiki/index.php/KPOINTS>`_
