================================
AIMD: Molecular Dynamics
================================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

This workflow demonstrates ab initio molecular dynamics (AIMD) on slab structures using PS-TEROS. AIMD allows you to study finite-temperature behavior, surface dynamics, and structural fluctuations.

**Prerequisites**:

* Complete :doc:`beginner workflow <beginner-surface-energy>` to get relaxed slabs
* Understanding of molecular dynamics concepts
* Sufficient computational resources (~10-20 hours per MD stage)
* VASP or CP2K code configured in AiiDA

**What You'll Learn**:

* Sequential AIMD staging (equilibration → production)
* Automatic restart chaining between stages
* Temperature control and ramping protocols
* Using VASP or CP2K calculators for AIMD
* Applying atomic constraints (fixed atoms)
* Trajectory analysis and visualization

**Materials**: Ag₂O or Ag₃PO₄ slabs from previous workflows

**Runtime**: Variable (depends on MD length, typically 10-50 hours)

**Supported Calculators**:

* **VASP**: Full support (recommended for most users)
* **CP2K**: AIMD support only (slab relaxation still uses VASP)

When to Use AIMD
================

**AIMD is useful for**:

* Studying temperature-dependent surface structure
* Investigating surface reconstruction dynamics
* Sampling multiple local minima
* Calculating temperature-dependent properties (free energy, entropy)
* Understanding diffusion and reaction mechanisms

**AIMD is NOT needed for**:

* Basic surface energy calculations (use :doc:`beginner workflow <beginner-surface-energy>`)
* Comparing different terminations at 0 K
* Quick screening calculations

AIMD Basics
===========

**Key Concepts**:

* **Equilibration**: Initial MD run to reach thermal equilibrium (~100-200 steps)
* **Production**: Main MD run for data collection (~500-1000+ steps)
* **Thermostat**: Controls temperature (Nosé-Hoover, Langevin, etc.)
* **Timestep (POTIM)**: Time interval between MD steps (typically 1-3 fs)

**VASP INCAR Tags**:

* ``IBRION=0``: Molecular dynamics mode
* ``MDALGO``: Thermostat algorithm (2 = Nosé-Hoover)
* ``POTIM``: Timestep in femtoseconds
* ``TEBEG/TEEND``: Initial/final temperature in Kelvin
* ``NSW``: Number of MD steps
* ``SMASS``: Nosé mass parameter

Sequential AIMD Workflow
=========================

PS-TEROS AIMD uses **sequential staging** with automatic restart:

.. code-block:: text

    Relaxed Slabs (from beginner/intermediate workflow)
         ↓
    ┌────────────────────────────────────┐
    │ Stage 1: Equilibration             │
    │  Temperature: 300 K                │
    │  Steps: 100                        │
    │  Purpose: Reach thermal equilibrium│
    └────────────────────────────────────┘
         ↓ (automatic restart via remote_folders)
    ┌────────────────────────────────────┐
    │ Stage 2: Production Run            │
    │  Temperature: 300 K                │
    │  Steps: 500                        │
    │  Purpose: Data collection          │
    │  Reads: WAVECAR from Stage 1       │
    └────────────────────────────────────┘
         ↓ (optional: extend)
    ┌────────────────────────────────────┐
    │ Stage 3: Extended Production       │
    │  Temperature: 300 K                │
    │  Steps: 1000                       │
    │  Purpose: More statistics          │
    │  Reads: WAVECAR from Stage 2       │
    └────────────────────────────────────┘

Each stage runs on **all slabs in parallel** (scatter-gather pattern).

Code Implementation
===================

1. Import AIMD Module
---------------------

.. code-block:: python

    from aiida import load_profile, orm
    from teros.core.aimd import aimd_single_stage_scatter

    load_profile(profile='your_profile_name')

2. Define Base AIMD Parameters
-------------------------------

Parameters that stay **constant** across all stages:

.. code-block:: python

    # Base AIMD parameters (constant for all stages)
    aimd_parameters = {
        'IBRION': 0,           # Molecular dynamics
        'MDALGO': 2,           # Nosé-Hoover thermostat
        'POTIM': 2.0,          # 2 fs timestep
        'SMASS': 3.0,          # Nosé mass (adjust for system)
        'PREC': 'Normal',      # Precision (Normal sufficient for MD)
        'ENCUT': 400,          # Slightly lower than relaxation (400-450)
        'EDIFF': 1e-5,         # Slightly relaxed electronic convergence
        'ISMEAR': 0,           # Gaussian smearing
        'SIGMA': 0.1,
        'LWAVE': True,         # Write WAVECAR for restart
        'LCHARG': True,        # Write CHGCAR for restart
        'LREAL': 'Auto',
        'ALGO': 'Normal',
    }

    # Scheduler options
    aimd_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'your_queue',
    }

    # Code and potentials
    code = orm.load_code('vasp@your_computer')
    potential_family = orm.Str('PBE')
    potential_mapping = {'Ag': 'Ag', 'O': 'O'}  # Adjust for your system
    kpoints_spacing = 0.5  # Coarser than relaxation (faster MD)

.. note::
   **POTIM**: Start with 1-2 fs. Too large = energy drift. Too small = slow convergence.
   **ENCUT**: Can be 10-15% lower than relaxation to speed up MD.
   **EDIFF**: Can be relaxed to 1e-5 (faster, still accurate for MD).

3. Get Relaxed Slabs
--------------------

Load slabs from a previous PS-TEROS workflow:

.. code-block:: python

    # Load previous workflow
    previous_wg = orm.load_node(<PREVIOUS_WORKFLOW_PK>)

    # Extract relaxed slab structures
    slab_structures = previous_wg.outputs.slab_structures.get_dict()

    # Load into dictionary
    slabs = {}
    for term_label, structure_pk in slab_structures.items():
        slabs[term_label] = orm.load_node(structure_pk)

    print(f"Loaded {len(slabs)} slabs: {list(slabs.keys())}")

4. Stage 1: Equilibration
--------------------------

.. code-block:: python

    print("\n=== Stage 1: Equilibration ===")

    # Run equilibration
    stage1 = aimd_single_stage_scatter(
        slabs=slabs,                      # Dict of StructureData
        temperature=300.0,                # Target temperature (K)
        steps=100,                        # Number of MD steps
        code=code,
        aimd_parameters=aimd_parameters,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=aimd_options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=False,
    )

    # Submit to daemon
    stage1.submit(wait=False)

    print(f"Stage 1 submitted: PK = {stage1.pk}")
    print(f"Monitor with: verdi process show {stage1.pk}")

Wait for Stage 1 to complete before proceeding.

5. Stage 2: Production Run with Restart
----------------------------------------

.. code-block:: python

    # After Stage 1 completes, load it
    stage1 = orm.load_node(<STAGE1_PK>)

    print("\n=== Stage 2: Production Run ===")

    # Run production MD with restart from Stage 1
    stage2 = aimd_single_stage_scatter(
        slabs=stage1.outputs.structures,           # Output from Stage 1
        temperature=300.0,
        steps=500,                                  # Longer production run
        code=code,
        aimd_parameters=aimd_parameters,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=aimd_options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=False,
        restart_folders=stage1.outputs.remote_folders,  # Restart from Stage 1
    )

    stage2.submit(wait=False)

    print(f"Stage 2 submitted: PK = {stage2.pk}")

**Key point**: ``restart_folders`` tells VASP to read WAVECAR and CHGCAR from Stage 1.

6. Optional: Stage 3 Extended Production
-----------------------------------------

.. code-block:: python

    # After Stage 2 completes
    stage2 = orm.load_node(<STAGE2_PK>)

    print("\n=== Stage 3: Extended Production ===")

    stage3 = aimd_single_stage_scatter(
        slabs=stage2.outputs.structures,
        temperature=300.0,
        steps=1000,                                 # Even longer
        code=code,
        aimd_parameters=aimd_parameters,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=aimd_options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=False,
        restart_folders=stage2.outputs.remote_folders,  # Restart from Stage 2
    )

    stage3.submit(wait=False)

Advanced: Temperature Ramping
==============================

To gradually heat/cool the system:

.. code-block:: python

    # Temperature ramping from 100 K to 300 K
    temperatures = [100, 150, 200, 250, 300]
    stages = []

    current_slabs = initial_slabs  # From relaxation
    current_folders = {}           # No restart for first stage

    for i, temp in enumerate(temperatures):
        print(f"\n=== Stage {i+1}: {temp} K ===")

        stage = aimd_single_stage_scatter(
            slabs=current_slabs,
            temperature=temp,
            steps=50,                    # Short equilibration at each T
            code=code,
            aimd_parameters=aimd_parameters,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=aimd_options,
            kpoints_spacing=kpoints_spacing,
            clean_workdir=False,
            restart_folders=current_folders,
        )

        stage.submit(wait=False)
        stages.append(stage)

        # Wait for this stage to complete before next
        # (In practice, manually check or use workflow management)
        print(f"  PK: {stage.pk}")

        # Update for next iteration (after stage completes)
        # current_slabs = stage.outputs.structures
        # current_folders = stage.outputs.remote_folders

Using CP2K for AIMD
===================

PS-TEROS supports CP2K as an alternative calculator for AIMD simulations. CP2K is particularly well-suited for large systems and offers efficient Born-Oppenheimer molecular dynamics.

.. note::
   **Current CP2K Support**: AIMD calculations only. Slab relaxation, surface energies, and other features still use VASP.

When to Use CP2K
----------------

**Use CP2K for AIMD when:**

* Working with large systems (>200 atoms)
* Need GPW/GAPW methods for efficient basis sets
* Require advanced metadynamics or path sampling
* Want tight-binding DFT (DFTB) for exploratory runs
* Need specific CP2K features (QMMM, excited states)

**Use VASP when:**

* Standard oxide surface calculations
* Need PAW accuracy
* Want consistent parameters with relaxation
* Familiar with VASP INCAR parameters

CP2K AIMD Workflow
------------------

1. Setup CP2K Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from teros.core.workgraph import build_core_workgraph
    from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

    # Get default CP2K AIMD parameters
    aimd_params = get_aimd_defaults_cp2k(
        cutoff=400,           # Plane wave cutoff (Ry)
        rel_cutoff=60,        # Relative cutoff (Ry)
        timestep=1.0,         # MD timestep (fs)
        eps_scf=1e-6,         # SCF convergence
        max_scf=40,           # Max SCF iterations
        thermostat='NOSE',    # Thermostat type
    )

    # Add KIND section for your system
    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {
            "_": "Ag",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11",
            "POTENTIAL": "GTH-PBE-q11",
        },
        {
            "_": "O",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6",
            "POTENTIAL": "GTH-PBE-q6",
        }
    ]

    # Scheduler options
    aimd_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 48,
        },
        'queue_name': 'normal',
    }

.. tip::
   **CP2K Cutoff Values**: 400 Ry cutoff with 60 Ry relative cutoff provides good accuracy for most oxides. Increase for heavier elements or tighter convergence.

2. Build Workflow with CP2K
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Define AIMD sequence
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    # Build workflow using 'aimd_only' preset
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',  # Use CP2K for AIMD

        # Structure inputs
        structures_dir='/path/to/structures',
        bulk_name='ag2o.cif',

        # VASP for bulk/slab generation
        code_label='VASP@computer',
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=bulk_params,
        bulk_options=bulk_options,

        # Slab generation parameters
        miller_indices=[1, 1, 1],
        min_slab_thickness=15.0,
        min_vacuum_thickness=15.0,

        # CP2K for AIMD
        aimd_code_label='CP2K@computer',
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        name='AIMD_CP2K_Workflow',
    )

    wg.submit(wait=False)
    print(f"Submitted: {wg.pk}")

**Workflow stages:**

1. Bulk relaxation with VASP
2. Slab generation from relaxed bulk
3. AIMD on slabs with CP2K (sequential stages)

3. Using Pre-existing Slabs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For AIMD-only workflows with manual slab structures:

.. code-block:: python

    from aiida import orm

    # Load pre-existing slabs
    input_slabs = {
        'slab_001': orm.load_node(<STRUCTURE_PK_1>),
        'slab_010': orm.load_node(<STRUCTURE_PK_2>),
    }

    # Build AIMD-only workflow (no bulk relaxation/slab generation)
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',

        # Input slabs directly (no bulk needed!)
        input_slabs=input_slabs,

        # CP2K AIMD
        aimd_code_label='CP2K@computer',
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        name='AIMD_CP2K_InputSlabs',
    )

    wg.submit(wait=False)

.. note::
   When using ``input_slabs``, bulk structure parameters are not required.

CP2K-Specific Parameters
-------------------------

Key CP2K AIMD parameters (in ``aimd_parameters``):

**Electronic Structure:**

.. code-block:: python

    'FORCE_EVAL': {
        'METHOD': 'QS',
        'DFT': {
            'BASIS_SET_FILE_NAME': 'BASIS_MOLOPT',  # Auto-provided
            'POTENTIAL_FILE_NAME': 'GTH_POTENTIALS',  # Auto-provided
            'MGRID': {
                'CUTOFF': 400,       # Plane wave cutoff (Ry)
                'REL_CUTOFF': 60,    # Relative cutoff (Ry)
                'NGRIDS': 4,         # Multi-grid levels
            },
            'SCF': {
                'EPS_SCF': 1e-6,     # Convergence threshold
                'MAX_SCF': 40,       # Max iterations
                'SCF_GUESS': 'ATOMIC',
                'OT': {              # Orbital transformation
                    'MINIMIZER': 'DIIS',
                    'PRECONDITIONER': 'FULL_SINGLE_INVERSE',
                },
            },
            'XC': {
                'XC_FUNCTIONAL': {'_': 'PBE'},  # Exchange-correlation
                'VDW_POTENTIAL': {               # Dispersion corrections
                    'POTENTIAL_TYPE': 'PAIR_POTENTIAL',
                    'PAIR_POTENTIAL': {
                        'TYPE': 'DFTD3',
                        'PARAMETER_FILE_NAME': 'dftd3.dat',
                        'REFERENCE_FUNCTIONAL': 'PBE',
                    },
                },
            },
        },
    }

**Molecular Dynamics:**

.. code-block:: python

    'MOTION': {
        'MD': {
            'ENSEMBLE': 'NVT',        # Canonical ensemble
            'TIMESTEP': 1.0,          # fs
            'TEMPERATURE': 300.0,      # Set automatically per stage
            'STEPS': 100,              # Set automatically per stage
            'THERMOSTAT': {
                'TYPE': 'NOSE',        # Nosé-Hoover thermostat
                'REGION': 'GLOBAL',    # Apply to all atoms
            },
        },
    }

.. warning::
   Do not manually set ``TEMPERATURE`` or ``STEPS`` in ``aimd_parameters``. These are automatically set by PS-TEROS based on ``aimd_sequence``.

Restart Between Stages
-----------------------

CP2K restart works similarly to VASP:

.. code-block:: python

    # Stage 1 output provides restart for Stage 2
    stage1_outputs = wg.tasks['aimd_stage_00_300K'].outputs

    # Access outputs
    final_structures = stage1_outputs.structures
    remote_folders = stage1_outputs.remote_folders  # Contains .restart files
    parameters = stage1_outputs.parameters           # Output parameters
    trajectories = stage1_outputs.trajectories       # MD trajectory

**Restart files used:**

* ``aiida-1.restart``: Wave functions and coordinates
* ``WAVECAR`` equivalent in CP2K format

Fixed Atoms Constraints
========================

PS-TEROS supports fixing atoms during AIMD (both VASP and CP2K) to simulate realistic surface conditions where bottom layers remain fixed.

When to Use Fixed Atoms
------------------------

**Use fixed atoms for:**

* Surface slab calculations (fix bottom to represent bulk)
* Preventing artificial drift of entire slab
* Mimicking experimental substrate constraints
* Studying adsorbate dynamics on rigid surfaces

**Typical setup:**

* Fix bottom 2-3 layers (~7-10 Å)
* Let top surface layers relax freely
* Maintain bulk-like structure at bottom

Fixed Atoms Implementation
---------------------------

1. Fixed Atoms with Auto-Generated Slabs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',  # or 'vasp'

        # Structure and workflow setup
        structures_dir='/path/to/structures',
        bulk_name='ag2o.cif',
        miller_indices=[1, 1, 1],
        min_slab_thickness=15.0,
        min_vacuum_thickness=15.0,

        # AIMD configuration
        aimd_code_label='CP2K@computer',
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        # Fixed atoms configuration
        fix_atoms=True,           # Enable fixed atoms
        fix_type='bottom',        # Fix bottom layers
        fix_thickness=7.0,        # Fix bottom 7 Å
        fix_elements=None,        # Fix all elements (or ['Ag', 'O'])
        fix_components='XYZ',     # Fix all components (fully rigid)

        name='AIMD_Fixed_Bottom',
    )

2. Fixed Atoms with Input Slabs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Load pre-existing slabs
    input_slabs = {
        'slab_111': orm.load_node(<STRUCTURE_PK>),
    }

    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',

        # Use input slabs
        input_slabs=input_slabs,

        # AIMD configuration
        aimd_code_label='CP2K@computer',
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        # Fixed atoms
        fix_atoms=True,
        fix_type='bottom',
        fix_thickness=7.0,
        fix_elements=None,
        fix_components='XYZ',

        name='AIMD_InputSlabs_Fixed',
    )

Fixed Atoms Parameters
----------------------

**fix_type**: Where to fix atoms

* ``'bottom'``: Fix atoms from bottom up to ``fix_thickness`` Å
* ``'top'``: Fix atoms from top down to ``fix_thickness`` Å  
* ``'center'``: Fix atoms within ``fix_thickness``/2 Å of slab center

**fix_thickness**: Thickness of fixed region (Angstroms)

* Typical: 7-10 Å for oxide surfaces
* Rule of thumb: 2-3 atomic layers
* Check slab structure to determine appropriate value

**fix_elements**: Which elements to fix (optional)

* ``None``: Fix all elements in the region (default)
* ``['Ag']``: Fix only Ag atoms
* ``['Ag', 'O']``: Fix both Ag and O atoms

**fix_components**: Which Cartesian components to constrain

* ``'XYZ'``: Fully rigid (default)
* ``'XY'``: Fix in-plane only, allow perpendicular movement
* ``'Z'``: Fix perpendicular only, allow in-plane movement

Example: Partial Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Fix oxygen atoms in bottom 5 Å, allow XY movement but not Z
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',

        # ... other parameters ...

        fix_atoms=True,
        fix_type='bottom',
        fix_thickness=5.0,
        fix_elements=['O'],     # Only oxygen
        fix_components='Z',     # Only perpendicular direction

        name='AIMD_Partial_Fix',
    )

Verification: Checking Fixed Atoms
-----------------------------------

For **CP2K**, check the input file:

.. code-block:: bash

    verdi calcjob inputcat <CP2K_CALC_PK> aiida.inp | grep -A 10 CONSTRAINT

Output should show:

.. code-block:: text

    &CONSTRAINT
      &FIXED_ATOMS
        COMPONENTS_TO_FIX XYZ
        LIST 1 2 3 4 5 6 7 8 21 22 23 24 26
      &END FIXED_ATOMS
    &END CONSTRAINT

For **VASP**, check POSCAR:

.. code-block:: bash

    verdi calcjob inputcat <VASP_CALC_PK> POSCAR

Selective dynamics section should show::

    Selective dynamics
    Cartesian
    0.0 0.0 0.0  F F F  # Fixed atom
    1.0 0.0 0.0  T T T  # Free atom

**Legend:**

* ``F F F``: Fixed in X, Y, Z
* ``T T T``: Free in X, Y, Z  
* ``T T F``: Free in XY, fixed in Z

Accessing AIMD Outputs
=======================

Structures and Energies
------------------------

.. code-block:: python

    from aiida.orm import load_node

    stage = load_node(<STAGE_PK>)

    # Final structures (last MD step)
    final_structures = stage.outputs.structures.get_dict()

    for term_label, structure_pk in final_structures.items():
        structure = load_node(structure_pk)
        print(f"{term_label}: {structure.get_formula()}")

    # Final energies
    final_energies = stage.outputs.energies.get_dict()

    for term_label, energy_pk in final_energies.items():
        energy = load_node(energy_pk).value
        print(f"{term_label}: {energy:.3f} eV")

    # Remote folders (for restart)
    remote_folders = stage.outputs.remote_folders.get_dict()

Trajectory Analysis
-------------------

AIMD trajectories are stored in AiiDA ``TrajectoryData``:

.. code-block:: python

    # Access VASP calculation node
    # (Need to navigate through WorkGraph structure)
    # Find the VaspWorkChain for specific slab

    # Example: Get trajectory for term_0
    # (Simplified - actual navigation depends on WorkGraph structure)

    vasp_calc = load_node(<VASP_CALC_PK>)

    if 'trajectory' in vasp_calc.outputs:
        trajectory = vasp_calc.outputs.trajectory

        # Get trajectory data
        positions = trajectory.get_positions()  # Shape: (n_steps, n_atoms, 3)
        cells = trajectory.get_cells()          # Shape: (n_steps, 3, 3)
        energies = trajectory.get_array('energies')  # Shape: (n_steps,)

        print(f"Trajectory length: {len(energies)} steps")
        print(f"Energy range: {min(energies):.3f} to {max(energies):.3f} eV")

.. tip::
   Use ``verdi process show <PK>`` to navigate the WorkGraph and find individual VASP calculations.

Visualization
=============

Plot Energy vs. Time
--------------------

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np

    # Assuming you have energies array from trajectory
    timestep_fs = 2.0  # POTIM value
    times = np.arange(len(energies)) * timestep_fs  # femtoseconds

    plt.figure(figsize=(10, 5))
    plt.plot(times, energies, 'k-', linewidth=0.5)
    plt.xlabel('Time (fs)')
    plt.ylabel('Energy (eV)')
    plt.title('AIMD Energy Evolution')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('aimd_energy.png', dpi=300)

    # Calculate average and std dev (after equilibration)
    eq_step = 50  # Skip first 50 steps
    avg_energy = np.mean(energies[eq_step:])
    std_energy = np.std(energies[eq_step:])

    print(f"Average energy: {avg_energy:.3f} ± {std_energy:.3f} eV")

Animate Trajectory (ASE)
-------------------------

.. code-block:: python

    from ase.io import write

    # Convert trajectory to ASE Atoms list
    # (Requires converting AiiDA TrajectoryData)

    # Example: If you have ASE atoms list
    # write('trajectory.xyz', atoms_list)
    # Then visualize with ASE gui or OVITO

Troubleshooting
===============

Problem: Energy drift
---------------------

**Symptoms**: Energy increases/decreases continuously

**Solution**:

1. Reduce ``POTIM`` (try 1.0-1.5 fs)
2. Tighten ``EDIFF`` (try 1e-6)
3. Check forces after relaxation (should be < 0.05 eV/Å)
4. Ensure thermostat is working (check ``MDALGO``, ``SMASS``)

Problem: Temperature not stabilizing
-------------------------------------

**Symptoms**: Temperature oscillates wildly or doesn't reach target

**Solution**:

1. Adjust ``SMASS`` (larger = slower temperature response)
2. Increase equilibration steps (try 200-300)
3. Check if system is too small (< 50 atoms may have large fluctuations)
4. Try different thermostat (``MDALGO=1`` for Andersen)

Problem: WAVECAR restart fails
-------------------------------

**Symptoms**: Stage 2+ fails immediately

**Solution**:

1. Verify ``LWAVE=True`` and ``LCHARG=True`` in aimd_parameters
2. Check ``clean_workdir=False`` (must keep remote data)
3. Ensure remote_folders from previous stage are correctly passed
4. Check file sizes: WAVECAR should be > 0 bytes

Performance Tips
================

**Speedup Strategies**:

* **Coarser k-points**: Use ``kpoints_spacing=0.5-0.6`` (vs 0.3-0.4 for relaxation)
* **Lower ENCUT**: 10-15% reduction from relaxation value
* **Relaxed EDIFF**: Use 1e-5 instead of 1e-6
* **Fewer atoms**: Test on smaller slabs first
* **Parallel efficiency**: AIMD scales well to 80-120 cores

**Time Estimates** (40 cores, 100-atom slab, ENCUT=400):

* 100 MD steps: ~2-4 hours
* 500 MD steps: ~10-20 hours
* 1000 MD steps: ~20-40 hours

Analysis Workflow
=================

Typical AIMD analysis pipeline:

1. **Check equilibration**: Plot energy vs. time, verify stabilization
2. **Discard equilibration**: Use only production run data
3. **Calculate averages**: Mean energy, structure, properties
4. **RDF/MSD**: Radial distribution functions, mean square displacement
5. **Visualization**: Animate trajectory, identify interesting events

.. code-block:: python

    # Example: Check if equilibrated
    def is_equilibrated(energies, window=50):
        """Check if last `window` steps have stable energy."""
        recent = energies[-window:]
        std = np.std(recent)
        mean = np.mean(recent)
        cv = std / abs(mean)  # Coefficient of variation
        return cv < 0.001  # < 0.1% variation

    if is_equilibrated(stage1_energies):
        print("✓ System equilibrated")
    else:
        print("✗ Need more equilibration steps")

Complete Example Script
========================

.. code-block:: python

    #!/usr/bin/env python
    """
    AIMD workflow: Sequential MD stages with automatic restart.

    Runs equilibration → production on relaxed Ag₂O slabs.
    """

    from aiida import load_profile, orm
    from teros.core.aimd import aimd_single_stage_scatter

    load_profile(profile='your_profile')

    # Base parameters
    aimd_parameters = {
        'IBRION': 0, 'MDALGO': 2, 'POTIM': 2.0, 'SMASS': 3.0,
        'PREC': 'Normal', 'ENCUT': 400, 'EDIFF': 1e-5,
        'ISMEAR': 0, 'SIGMA': 0.1,
        'LWAVE': True, 'LCHARG': True,
        'LREAL': 'Auto', 'ALGO': 'Normal',
    }

    options = {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 40},
        'queue_name': 'your_queue',
    }

    code = orm.load_code('vasp@your_computer')

    # Load relaxed slabs
    previous_wg = orm.load_node(<PREVIOUS_PK>)
    slab_pks = previous_wg.outputs.slab_structures.get_dict()
    slabs = {label: orm.load_node(pk) for label, pk in slab_pks.items()}

    # Stage 1: Equilibration
    stage1 = aimd_single_stage_scatter(
        slabs=slabs, temperature=300, steps=100,
        code=code, aimd_parameters=aimd_parameters,
        potential_family=orm.Str('PBE'),
        potential_mapping={'Ag': 'Ag', 'O': 'O'},
        options=options, kpoints_spacing=0.5, clean_workdir=False,
    )
    stage1.submit(wait=False)
    print(f"Stage 1 PK: {stage1.pk}")

    # After Stage 1 completes, run Stage 2:
    # stage1 = orm.load_node(<STAGE1_PK>)
    # stage2 = aimd_single_stage_scatter(
    #     slabs=stage1.outputs.structures, temperature=300, steps=500,
    #     # ... same parameters ...
    #     restart_folders=stage1.outputs.remote_folders,
    # )
    # stage2.submit(wait=False)

Next Steps
==========

You've learned AIMD with PS-TEROS:

✓ Sequential staging with restart
✓ Temperature control and ramping
✓ Trajectory analysis
✓ Troubleshooting common issues

**For More**:

* :doc:`How-To: AIMD Stages <../how-to/aimd-stages>` - Advanced staging patterns
* :doc:`API: AIMD Module </api/aimd>` - Function reference
* VASP Manual - AIMD tags and theory

**Further Analysis**:

* `MDAnalysis <https://www.mdanalysis.org/>`_ - Python trajectory analysis
* `OVITO <https://www.ovito.org/>`_ - Visualization
* `Pymatgen <https://pymatgen.org/>`_ - Structure analysis tools
