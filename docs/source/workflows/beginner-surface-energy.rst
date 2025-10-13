========================================
Beginner: Surface Energy Calculation
========================================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

This workflow demonstrates the essential features of PS-TEROS for calculating surface energies of binary oxide materials. You'll learn the complete workflow from structure input to surface energy calculation.

**Material**: Ag₂O (silver oxide, cuprite structure)

**What You'll Calculate**:

1. Relaxed bulk structure (Ag₂O)
2. Relaxed reference structures (Ag metal, O₂ molecule)
3. Formation enthalpy of Ag₂O
4. Multiple slab terminations for (111) surface
5. Relaxed slab structures
6. Surface energies as function of oxygen chemical potential

**Prerequisites**:

* AiiDA profile configured and daemon running
* VASP code set up in AiiDA
* Basic understanding of DFT calculations
* ~2-4 hours of cluster time

**Estimated Runtime**: 2-4 hours on typical cluster (40 cores)

Understanding the Workflow
===========================

Before diving into code, let's understand what PS-TEROS does:

.. code-block:: text

    Input Structures
         ↓
    ┌────────────────────────────────┐
    │ Parallel Relaxation:           │
    │  - Bulk (Ag₂O)                │
    │  - Metal reference (Ag)        │
    │  - Oxygen reference (O₂)       │
    └────────────────────────────────┘
         ↓
    Formation Enthalpy (ΔH_f)
         ↓
    Generate Slab Terminations
    (pymatgen SlabGenerator)
         ↓
    ┌────────────────────────────────┐
    │ Parallel Slab Relaxation:      │
    │  - term_0, term_1, term_2, ... │
    └────────────────────────────────┘
         ↓
    Surface Energy Calculation
    γ(μ_O) for each termination

**Key Concepts**:

* **Scatter-gather pattern**: All parallel operations (bulk/metal/O₂ relaxations, slab relaxations) run simultaneously
* **Provenance tracking**: AiiDA tracks all inputs, outputs, and dependencies
* **Automatic slab generation**: pymatgen generates symmetrically distinct terminations
* **Chemical potential dependence**: Surface energy depends on O₂ reservoir conditions

Step-by-Step Workflow
======================

1. Setup and Imports
--------------------

Create a new Python script ``beginner_ag2o.py``:

.. code-block:: python

    #!/usr/bin/env python
    """
    Beginner PS-TEROS workflow: Basic surface energy calculation for Ag₂O.

    This script demonstrates core functionality without advanced features.
    """

    import os
    from aiida import load_profile
    from teros.core.workgraph import build_core_workgraph

    def main():
        # Load your AiiDA profile
        load_profile(profile='your_profile_name')  # Replace with your profile

        # Define paths
        script_dir = os.path.dirname(os.path.abspath(__file__))
        structures_dir = os.path.join(script_dir, 'structures')

**Directory Structure**:

.. code-block:: text

    your_project/
    ├── beginner_ag2o.py
    └── structures/
        ├── ag2o.cif      # Bulk Ag₂O structure
        ├── Ag.cif        # Ag metal reference
        └── O2.cif        # O₂ molecule reference

**Where to get structures**:

* Materials Project (https://materialsproject.org)
* Your own optimized structures
* Literature structures

2. Define Calculation Parameters
---------------------------------

**Bulk Parameters (Ag₂O)**:

.. code-block:: python

    # Structure filenames
    bulk_filename = 'ag2o.cif'
    metal_filename = 'Ag.cif'
    oxygen_filename = 'O2.cif'

    # VASP code and potentials
    code_label = 'vasp@your_computer'  # Your VASP code label
    potential_family = 'PBE'           # Pseudopotential family

    # Bulk relaxation INCAR parameters
    bulk_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,              # Plane-wave cutoff (eV)
        'EDIFF': 1e-6,            # Electronic convergence
        'ISMEAR': 0,               # Gaussian smearing (semiconductors)
        'SIGMA': 0.05,             # Smearing width (eV)
        'IBRION': 2,               # Conjugate gradient relaxation
        'ISIF': 3,                 # Relax cell + atoms
        'NSW': 100,                # Max ionic steps
        'EDIFFG': -0.01,          # Force convergence (eV/Å)
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options
    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'your_queue',
    }

    # Pseudopotential mapping
    bulk_potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
    }

.. note::
   **ISIF=3** allows full cell relaxation (shape + volume + atoms). Critical for bulk but NOT used for slabs.

**Metal Reference Parameters (Ag)**:

.. code-block:: python

    metal_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 1,           # Methfessel-Paxton for metals
        'SIGMA': 0.2,          # Wider smearing for metals
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    metal_options = bulk_options.copy()  # Same resources

    metal_potential_mapping = {
        'Ag': 'Ag',
    }

.. important::
   **ISMEAR=1** (Methfessel-Paxton) is appropriate for metals. Do NOT use ISMEAR=0 for metals.

**Oxygen Reference Parameters (O₂)**:

.. code-block:: python

    oxygen_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,             # Relax atoms only (NOT cell for molecules)
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': False,        # Use reciprocal space for small systems
        'LWAVE': False,
        'LCHARG': False,
    }

    oxygen_options = bulk_options.copy()

    oxygen_potential_mapping = {
        'O': 'O',
    }

.. warning::
   **LREAL=False** for O₂ molecule! Small systems need reciprocal space projection. **ISIF=2** keeps cell fixed while relaxing atoms.

**Slab Parameters**:

.. code-block:: python

    # Slab generation
    miller_indices = [1, 1, 1]    # (111) surface
    min_slab_thickness = 15.0     # Minimum slab thickness (Å)
    min_vacuum_thickness = 15.0   # Minimum vacuum gap (Å)

    # Slab relaxation INCAR
    slab_parameters = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'NELM': 100,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.1,
        'IBRION': 2,
        'ISIF': 2,                # Relax atoms only (NOT cell for slabs)
        'NSW': 100,
        'EDIFFG': -0.1,          # Slightly relaxed for surfaces
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        'LASPH': True,           # Aspherical contributions
    }

    slab_options = bulk_options.copy()

.. note::
   **ISIF=2** for slabs: Only atoms relax, cell is fixed. **EDIFFG=-0.1** is slightly relaxed compared to bulk (-0.01) because surface convergence is slower.

3. Build and Submit Workflow
-----------------------------

.. code-block:: python

    # Build the workflow
    wg = build_core_workgraph(
        # ===== Required Parameters =====
        structures_dir=structures_dir,
        bulk_name=bulk_filename,
        code_label=code_label,
        potential_family=potential_family,

        # ===== Bulk Parameters =====
        bulk_potential_mapping=bulk_potential_mapping,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,

        # ===== Reference Structures =====
        metal_name=metal_filename,
        metal_potential_mapping=metal_potential_mapping,
        metal_parameters=metal_parameters,
        metal_options=metal_options,

        oxygen_name=oxygen_filename,
        oxygen_potential_mapping=oxygen_potential_mapping,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,

        # ===== Slab Generation =====
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,

        # ===== Slab Relaxation =====
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        relax_slabs=True,

        # ===== Feature Flags (Beginner: Basic only) =====
        compute_relaxation_energy=False,      # Not needed for basic workflow
        compute_cleavage=False,                # Not needed for basic workflow
        compute_thermodynamics=True,           # Essential: Surface energies
        thermodynamics_sampling=10,            # Chemical potential grid points

        # ===== Other Settings =====
        kpoints_spacing=0.4,                   # K-points spacing (Å⁻¹)
        clean_workdir=False,                   # Keep files for inspection

        name='Ag2O_Beginner_Workflow',
    )

    print(f"WorkGraph created: {wg.name}")
    print(f"Total tasks: {len(wg.tasks)}")

    # Submit to AiiDA daemon
    wg.submit(wait=False)

    print(f"\n✓ Workflow submitted!")
    print(f"  PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")

    return wg

if __name__ == '__main__':
    main()

**Run the workflow**:

.. code-block:: bash

    # 1. Activate AiiDA environment
    source ~/envs/aiida/bin/activate

    # 2. Ensure daemon is running
    verdi daemon status
    verdi daemon start  # If not running

    # 3. Run the script
    python beginner_ag2o.py

4. Monitor the Workflow
-----------------------

**Check overall status**:

.. code-block:: bash

    verdi process show <PK>

**Example output**:

.. code-block:: text

    Property     Value
    -----------  ------------------------------------
    type         WorkGraphNode
    state        ProcessState.RUNNING
    pk           12345
    uuid         abc123...
    ctime        2025-10-13 15:30:00
    mtime        2025-10-13 15:35:00

    Inputs      PK    Type
    ----------  ----  -------------
    structures  ...   StructureData
    ...

    Outputs            PK    Type
    -----------------  ----  -------------
    bulk_energy        ...   Float
    metal_energy       ...   Float
    oxygen_energy      ...   Float
    formation_enthalpy ...   Float
    slab_structures    ...   Dict
    slab_energies      ...   Dict
    surface_energies   ...   Dict

**Check for errors**:

.. code-block:: bash

    verdi process report <PK>

**Monitor all running processes**:

.. code-block:: bash

    verdi process list
    verdi process list -a -p1  # More detail

.. tip::
   If a calculation fails, check ``verdi process report <CALC_PK>`` for the specific calculation that failed.

5. Understanding the Outputs
=============================

When the workflow completes (state = ``Finished [0]``), you can access outputs:

**Energy Outputs**:

.. code-block:: python

    from aiida import load_profile
    from aiida.orm import load_node

    load_profile()

    wg = load_node(<YOUR_PK>)

    # Access outputs
    bulk_energy = wg.outputs.bulk_energy.value  # eV
    metal_energy = wg.outputs.metal_energy.value  # eV
    oxygen_energy = wg.outputs.oxygen_energy.value  # eV
    formation_enthalpy = wg.outputs.formation_enthalpy.value  # eV/formula unit

    print(f"Formation enthalpy: {formation_enthalpy:.3f} eV/formula unit")

**Slab Structures**:

.. code-block:: python

    slab_structures = wg.outputs.slab_structures.get_dict()

    print("Generated terminations:")
    for term_label, structure_pk in slab_structures.items():
        structure = load_node(structure_pk)
        print(f"  {term_label}: {structure.get_formula()} ({structure.get_ase().get_global_number_of_atoms()} atoms)")

**Surface Energies**:

.. code-block:: python

    surface_energies = wg.outputs.surface_energies.get_dict()

    # Each termination has surface energies at different μ_O values
    for term_label, energies_pk in surface_energies.items():
        energies = load_node(energies_pk).get_dict()
        print(f"\n{term_label}:")
        print(f"  μ_O range: {energies['mu_O'][0]:.3f} to {energies['mu_O'][-1]:.3f} eV")
        print(f"  γ range: {min(energies['gamma']):.3f} to {max(energies['gamma']):.3f} eV/Ų")

**Visualization**:

.. code-block:: python

    import matplotlib.pyplot as plt

    # Plot surface energy vs. chemical potential
    for term_label, energies_pk in surface_energies.items():
        energies = load_node(energies_pk).get_dict()
        plt.plot(energies['mu_O'], energies['gamma'], label=term_label)

    plt.xlabel('μ_O (eV)')
    plt.ylabel('γ (eV/Ų)')
    plt.legend()
    plt.title('Surface Energy vs. Oxygen Chemical Potential')
    plt.grid(True)
    plt.savefig('surface_energies.png', dpi=300)

Expected Results
================

**Typical Outputs for Ag₂O (111)**:

* **Formation enthalpy**: Around -0.3 to -0.5 eV/formula unit (compare with experimental: ~-0.31 eV)
* **Number of terminations**: 2-4 symmetrically distinct terminations
* **Surface energy range**: 0.5-2.0 eV/Ų depending on termination and μ_O
* **Calculation time**: 2-4 hours on 40-core cluster

**Sanity Checks**:

✓ Formation enthalpy should be **negative** (Ag₂O is stable)
✓ Surface energies should be **positive**
✓ Different terminations show different γ(μ_O) slopes
✓ At least one termination is stable across the μ_O range

Troubleshooting
===============

Problem: Bulk relaxation not converging
----------------------------------------

**Symptoms**: ``NSW`` reaches 100 without ``EDIFFG`` convergence

**Solution**:

1. Check starting structure quality (use Materials Project optimized structures)
2. Increase ``NSW`` to 200
3. Try ``IBRION=1`` (RMM-DIIS) instead of IBRION=2
4. Check if ``EDIFFG=-0.01`` is too strict, try ``-0.02``

Problem: O₂ calculation fails
------------------------------

**Symptoms**: O₂ relaxation gives errors

**Solution**:

1. Ensure ``LREAL=False`` (required for molecules)
2. Check O₂ structure has sufficient vacuum (>10 Å from cell boundaries)
3. Verify ``ISIF=2`` (NOT 3) to keep cell fixed

Problem: No slab terminations generated
----------------------------------------

**Symptoms**: ``slab_structures`` output is empty

**Solution**:

1. Check ``miller_indices`` are valid for your bulk structure
2. Increase ``min_slab_thickness`` (try 20 Å)
3. Verify bulk structure has correct symmetry

Problem: Surface energies seem wrong
-------------------------------------

**Symptoms**: All γ values are negative or unreasonably large

**Solution**:

1. Check formation enthalpy is correct (should be negative for stable oxides)
2. Verify slab relaxations converged (check ``verdi process report <PK>``)
3. Ensure consistent ENCUT across all calculations

Next Steps
==========

Congratulations! You've completed the beginner workflow. You now understand:

✓ How to set up PS-TEROS calculations
✓ Parameter inheritance (bulk → metal → oxygen → slab)
✓ Workflow structure and monitoring
✓ Output interpretation and visualization

**Ready for more?**

* :doc:`Intermediate Workflow <intermediate-with-features>` - Add relaxation energies, cleavage energies, and electronic properties
* :doc:`How-To Guides </how-to/index>` - Task-specific solutions
* :doc:`API Reference </api/index>` - Detailed parameter documentation

**Questions or issues?**

* Check the :doc:`theory page </theory>` for background
* See :doc:`installation </installation>` for setup issues
* Open an issue on `GitHub <https://github.com/your-repo/PS-TEROS/issues>`_

Complete Example Script
========================

Download the complete beginner script: :download:`beginner_ag2o.py <../../../examples/complete/beginner_ag2o.py>`

.. literalinclude:: ../../../examples/complete/beginner_ag2o_simplified.py
   :language: python
   :linenos:
   :caption: Complete beginner workflow (simplified)

.. note::
   This simplified version removes advanced features to focus on core concepts. See ``examples/complete/complete_ag2o_example.py`` for the full-featured version.
