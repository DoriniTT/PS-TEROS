========================================
How to Use Custom Slab Structures
========================================

This guide shows you how to use manually prepared slab structures in PS-TEROS instead of auto-generated ones.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

PS-TEROS can automatically generate slabs from bulk structures, but sometimes you need manual control:

* Non-standard surface terminations
* Pre-relaxed slab geometries
* Slabs with defects or adsorbates
* Surface reconstructions
* Specific atomic arrangements

**Key concepts:**

* **input_slabs** - Dictionary mapping slab names to structure data
* **Structure formats** - AiiDA StructureData or file paths (POSCAR, CIF)
* **Slab naming** - How PS-TEROS identifies and tracks individual slabs
* **Validation** - Ensuring slab quality before expensive calculations

When to Use Custom Slabs
=========================

Use custom slabs when:

**Automatic generation insufficient:**

* Surface reconstructs (e.g., Au(100) 5×20)
* Polar surface terminations need manual compensation
* You need specific stoichiometry per slab

**Continuing previous work:**

* Starting from previously relaxed structures
* Adding adsorbates to converged surfaces
* Testing small modifications to existing slabs

**Special requirements:**

* Comparing with literature structures
* Non-standard vacuum thickness
* Pre-defined defect configurations

**Don't use custom slabs if:**

* Standard surfaces (low-index, symmetric terminations)
* No special requirements → Auto-generation is faster and less error-prone

Preparing Custom Slabs
=======================

Method 1: Load from Files
--------------------------

Most common approach - prepare POSCAR files with ASE, pymatgen, or VESTA:

**Step 1: Create slab directory**

.. code-block:: bash

   mkdir custom_slabs
   cd custom_slabs

**Step 2: Generate slab structures**

Using ASE:

.. code-block:: python

   from ase.build import bulk, surface
   from ase.io import write

   # Create Ag2O bulk
   bulk_ag2o = bulk('Ag2O', 'fcc', a=4.72)

   # Generate (111) slab with 5 layers
   slab = surface(bulk_ag2o, (1,1,1), layers=5, vacuum=15.0)

   # Save as POSCAR
   write('ag2o_111_term1.vasp', slab)

**Repeat for multiple terminations:**

.. code-block:: python

   # Termination 2: Shift and remove top layer
   slab2 = surface(bulk_ag2o, (1,1,1), layers=5, vacuum=15.0)
   # ... custom modifications ...
   write('ag2o_111_term2.vasp', slab2)

   # Different surface
   slab_100 = surface(bulk_ag2o, (1,0,0), layers=7, vacuum=15.0)
   write('ag2o_100_term1.vasp', slab_100)

**Step 3: Verify structures**

Always check before running calculations:

.. code-block:: bash

   # Check with ASE
   python -c "from ase.io import read; print(read('ag2o_111_term1.vasp'))"

   # Or use VESTA for visual inspection
   vesta ag2o_111_term1.vasp

**Validation checklist:**

1. ☐ Vacuum thickness ≥15 Å
2. ☐ No atomic overlaps (minimum distance >1.5 Å for metals, >1.0 Å for oxides)
3. ☐ Slab thickness sufficient (typically 5-10 layers)
4. ☐ Cell vectors reasonable (no extreme angles or ratios)
5. ☐ Stoichiometry correct for your surface termination

Method 2: From Previous Calculations
-------------------------------------

Reuse structures from earlier PS-TEROS or VASP runs:

.. code-block:: python

   from aiida import load_profile, orm
   load_profile()

   # Load previous workflow
   previous_wg = orm.load_node(<PK>)

   # Extract relaxed slab structures
   relaxed_slabs = {}
   for label, structure in previous_wg.outputs.items():
       if 'relaxed_slab_' in label:
           slab_name = label.replace('relaxed_slab_', '')
           relaxed_slabs[slab_name] = structure

   # Save to files for inspection
   for name, structure in relaxed_slabs.items():
       structure.get_ase().write(f'{name}.vasp')

Method 3: Modify Auto-Generated Slabs
--------------------------------------

Generate with PS-TEROS, modify, then reuse:

**Step 1: Generate initial slabs**

.. code-block:: python

   from teros.core.slabs import generate_slabs_scatter
   from aiida import orm

   # Load bulk structure
   bulk = orm.load_node(<BULK_PK>)

   # Generate slabs normally
   slabs = generate_slabs_scatter(
       structures=[bulk],
       max_index=2,
       min_slab_size=15,
       min_vacuum_size=15,
   )

**Step 2: Extract and modify**

.. code-block:: python

   # Save generated slabs
   for label, slab_data in slabs.items():
       structure = slab_data['structure']
       ase_structure = structure.get_ase()

       # Apply modifications
       # Example: Add adsorbate
       from ase import Atoms
       h2o = Atoms('OH2', positions=[[0,0,20], [0.76,0.59,20], [-0.76,0.59,20]])
       ase_structure += h2o

       # Save modified structure
       ase_structure.write(f'modified_{label}.vasp')

Using Custom Slabs in PS-TEROS
===============================

Option 1: From File Directory
------------------------------

Most straightforward - point to directory with POSCAR files:

.. code-block:: python

   from teros.core.workgraph import build_core_workgraph

   # Load bulk structure (still needed for reference calculations)
   structures_dir = '/path/to/structures'
   bulk_filename = 'Ag2O.vasp'

   # Directory with custom slabs
   custom_slabs_dir = '/path/to/custom_slabs'

   wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       slab_dir=custom_slabs_dir,  # Use custom slabs instead of generating
       code_label='vasp@cluster',
       potential_family='PBE.54',

       # Standard parameters
       bulk_parameters=bulk_parameters,
       slab_parameters=slab_parameters,
       bulk_options=bulk_options,
       slab_options=slab_options,

       compute_thermodynamics=True,
       thermodynamics_sampling=10,
   )

   wg.submit(wait=True)

**Slab naming convention:**

PS-TEROS expects files named: ``{material}_{hkl}_term{n}.vasp``

Examples:
- ``ag2o_111_term1.vasp``
- ``ag2o_111_term2.vasp``
- ``ag2o_100_term1.vasp``
- ``ag3po4_001_term1.vasp``

**If your files have different names**, rename them or use Option 2.

Option 2: From AiiDA Nodes
---------------------------

Use structures already in AiiDA database:

.. code-block:: python

   from aiida import orm

   # Load or create AiiDA StructureData nodes
   slab1 = orm.load_node(<PK1>)
   slab2 = orm.load_node(<PK2>)
   slab3 = orm.load_node(<PK3>)

   # Create input_slabs dictionary
   input_slabs = {
       'ag2o_111_term1': slab1,
       'ag2o_111_term2': slab2,
       'ag2o_100_term1': slab3,
   }

   # Pass to workflow
   wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       input_slabs=input_slabs,  # Provide pre-loaded structures
       code_label='vasp@cluster',
       potential_family='PBE.54',

       bulk_parameters=bulk_parameters,
       slab_parameters=slab_parameters,
       bulk_options=bulk_options,
       slab_options=slab_options,

       compute_thermodynamics=True,
       thermodynamics_sampling=10,
   )

   wg.submit(wait=True)

Option 3: Mixed Auto + Custom
------------------------------

Combine auto-generated slabs with custom ones:

.. code-block:: python

   from teros.core.slabs import generate_slabs_scatter
   from aiida import orm

   # Auto-generate some surfaces
   bulk = orm.StructureData(ase=bulk_ase)
   auto_slabs = generate_slabs_scatter(
       structures=[bulk],
       max_index=1,  # Only (100), (110), (111)
       min_slab_size=15,
   )

   # Load custom reconstructed surface
   custom_slab = orm.StructureData(ase=read('ag2o_111_reconstructed.vasp'))

   # Combine
   all_slabs = {**auto_slabs}  # Start with auto-generated
   all_slabs['ag2o_111_reconstructed'] = custom_slab  # Add custom

   # Use in workflow
   wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       input_slabs=all_slabs,  # Mixed set
       # ... rest of parameters ...
   )

Validation and Quality Checks
==============================

Before Running Calculations
----------------------------

**Check 1: Structure integrity**

.. code-block:: python

   from ase.io import read
   from ase.geometry.analysis import Analysis

   slab = read('ag2o_111_term1.vasp')

   # Check atomic distances
   analysis = Analysis(slab)
   bonds = analysis.get_bonds('Ag', 'O', unique=True)

   print(f"Number of Ag-O bonds: {len(bonds[0])}")
   print(f"Shortest Ag-O distance: {min(bonds[1]):.3f} Å")

   # Should be ~2.0-2.5 Å for Ag-O
   assert min(bonds[1]) > 1.5, "Atoms too close - check structure!"

**Check 2: Vacuum and periodicity**

.. code-block:: python

   # Check cell parameters
   cell = slab.get_cell()
   print(f"Cell vectors:\n{cell}")

   # Check vacuum
   positions = slab.get_positions()
   z_extent = positions[:, 2].max() - positions[:, 2].min()
   vacuum = cell[2, 2] - z_extent

   print(f"Slab thickness: {z_extent:.2f} Å")
   print(f"Vacuum thickness: {vacuum:.2f} Å")

   assert vacuum >= 15.0, "Vacuum too small - increase to ≥15 Å"

**Check 3: Stoichiometry**

.. code-block:: python

   from collections import Counter

   # Count atoms
   symbols = slab.get_chemical_symbols()
   composition = Counter(symbols)

   print(f"Composition: {dict(composition)}")

   # For Ag2O, expect Ag:O ratio ~ 2:1
   ag_count = composition['Ag']
   o_count = composition['O']
   ratio = ag_count / o_count

   print(f"Ag:O ratio = {ratio:.2f} (expect ~2.0 for stoichiometric Ag2O)")

After First Test Calculation
-----------------------------

Always run a single test before full workflow:

.. code-block:: python

   # Test one slab first
   test_input = {
       'ag2o_111_term1': orm.load_node(<CUSTOM_SLAB_PK>)
   }

   test_wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       input_slabs=test_input,  # Just one slab

       # ... parameters ...

       compute_thermodynamics=False,  # Skip for test
   )

   test_wg.submit(wait=True)

**Check test results:**

.. code-block:: bash

   # Wait for completion
   sleep 600

   # Check if successful
   verdi process show <WG_PK>

   # If successful, check output structure
   verdi shell
   >>> node = load_node(<WG_PK>)
   >>> relaxed = node.outputs.relaxed_slab_ag2o_111_term1
   >>> relaxed.get_ase().get_positions()

**Red flags:**

* Large forces (>0.5 eV/Å) after relaxation
* Significant structural changes (>0.5 Å atom displacements)
* Unexpected energy (compare with similar systems)
* Unusual magnetic moments or charges

Common Issues and Solutions
============================

Issue 1: Structure Won't Relax
-------------------------------

**Symptoms:**

* NSW limit reached without convergence
* Forces remain large (>0.1 eV/Å)
* Energy oscillates

**Solutions:**

.. code-block:: python

   # More gradual relaxation
   slab_parameters = {
       'IBRION': 1,      # RMM-DIIS (more stable)
       'POTIM': 0.2,     # Smaller steps
       'EDIFFG': -0.03,  # Relaxed criterion initially
       'NSW': 300,       # More steps allowed
   }

   # Or: Pre-relax with cheaper settings
   prerelax_parameters = {
       'ENCUT': 400,     # Lower cutoff
       'PREC': 'Normal',  # Lower precision
       'EDIFFG': -0.05,  # Relaxed criterion
   }

Issue 2: Slab Distorts Unreasonably
------------------------------------

**Symptoms:**

* Surface reconstruction not in literature
* Atoms move >1 Å from initial positions
* Bulk-like layers change significantly

**Possible causes:**

1. Initial structure has strain or defects
2. Slab too thin (bulk layers can't stabilize)
3. Wrong pseudopotentials
4. Vacuum too small (interaction with periodic image)

**Solutions:**

.. code-block:: python

   # Fix bottom layers
   from ase.constraints import FixAtoms

   slab = read('ag2o_111_term1.vasp')

   # Fix bottom 2 layers (assuming 5 total)
   z_positions = slab.get_positions()[:, 2]
   z_sorted = sorted(set(z_positions))
   fix_layers = z_sorted[:2]  # Bottom 2 unique z values

   constraints = [FixAtoms(indices=[i for i, z in enumerate(z_positions)
                                     if z in fix_layers])]
   slab.set_constraint(constraints)

   # Save constrained structure
   write('ag2o_111_term1_constrained.vasp', slab)

Then inform PS-TEROS about constraints via VASP selective dynamics:

.. code-block:: python

   slab_parameters = {
       'IBRION': 2,
       'ISIF': 2,  # Relax ions only, not cell
       # Constraints read from POSCAR automatically
   }

Issue 3: Energy Much Higher Than Expected
------------------------------------------

**Possible causes:**

* Slab termination is unstable (polar, undercoordinated)
* Wrong stoichiometry (too metal-rich or oxygen-rich)
* Vacuum not large enough
* Structure still has overlapping atoms

**Diagnosis:**

.. code-block:: bash

   # Compare with auto-generated slab energy
   verdi shell
   >>> auto_node = load_node(<AUTO_GENERATED_SLAB_PK>)
   >>> custom_node = load_node(<CUSTOM_SLAB_PK>)
   >>>
   >>> auto_energy = auto_node.outputs.output_parameters['energy']
   >>> custom_energy = custom_node.outputs.output_parameters['energy']
   >>>
   >>> print(f"Energy difference: {custom_energy - auto_energy:.3f} eV")

If difference >0.5 eV/atom, investigate structure quality.

Best Practices
==============

Slab Construction Guidelines
-----------------------------

**1. Vacuum thickness**

* Minimum: 15 Å
* Recommended: 20 Å for charged surfaces or adsorbates
* Check: No charge density in vacuum region (CHGCAR)

**2. Slab thickness**

* Binary oxides (Ag₂O): 5-7 layers minimum
* Ternary oxides (Ag₃PO₄): 7-10 layers
* Rule of thumb: Bulk-like behavior in center layers

**3. Surface termination**

* Prefer stoichiometric or nearly stoichiometric
* Avoid highly polar surfaces (unless compensated)
* Check coordination: Surface atoms should have >50% of bulk coordination

**4. Lattice parameters**

* Use experimental or DFT-relaxed bulk values
* Don't mix different functionals/settings between bulk and slab
* For heterostructures, consider strain effects

Workflow Integration
--------------------

**1. Consistent parameters**

Use same DFT settings for bulk and slabs:

.. code-block:: python

   # Shared settings
   shared_params = {
       'PREC': 'Accurate',
       'ENCUT': 520,
       'EDIFF': 1e-6,
       'ISMEAR': 0,
       'SIGMA': 0.05,
   }

   bulk_parameters = {
       **shared_params,
       'ISIF': 3,    # Relax cell
       'IBRION': 2,
   }

   slab_parameters = {
       **shared_params,
       'ISIF': 2,    # Relax ions only
       'IBRION': 2,
   }

**2. Naming consistency**

Keep names descriptive and systematic:

.. code-block:: text

   Good:
   - ag2o_111_term1.vasp
   - ag2o_111_term2.vasp
   - ag2o_111_oh_adsorbed.vasp

   Bad:
   - slab1.vasp
   - test.vasp
   - final_FINAL_v3.vasp

**3. Documentation**

Create a README in your slab directory:

.. code-block:: text

   custom_slabs/README.md:

   # Custom Ag2O Slabs

   ## ag2o_111_term1.vasp
   - Standard (111) termination, Ag-rich
   - 5 layers, 15 Å vacuum
   - Generated with ASE
   - Reference: Smith et al. J. Phys. Chem. C 2020

   ## ag2o_111_reconstructed.vasp
   - Missing-row reconstruction
   - Modified from DFT relaxation of term1
   - 7 layers, 20 Å vacuum

Related Guides
==============

* :doc:`restart-calculations` - What to do if custom slabs fail to converge
* :doc:`electronic-properties` - Electronic structure of custom surfaces
* :doc:`../workflows/intermediate-with-features` - Using custom slabs in full workflows

External Resources
==================

* `ASE surface builder <https://wiki.fysik.dtu.dk/ase/ase/build/surface.html>`_
* `pymatgen Slab class <https://pymatgen.org/pymatgen.core.surface.html>`_
* `VESTA structure editor <http://jp-minerals.org/vesta/en/>`_
