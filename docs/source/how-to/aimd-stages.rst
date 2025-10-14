============================
How to Set Up AIMD Stages
============================

This guide shows you how to design and execute multi-stage AIMD simulations with restart chaining in PS-TEROS.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

Ab initio molecular dynamics (AIMD) simulations often require multiple stages:

* **Equilibration** - Bring system to target temperature/pressure
* **Production** - Collect statistics at equilibrium
* **Annealing** - Temperature ramping protocols
* **Long runs** - Break into chunks to fit walltime limits

**Key concepts:**

* **restart_folders** - Remote directories with WAVECAR for continuation
* **Sequential staging** - Chain outputs from one stage as inputs to next
* **Temperature protocols** - Gradual heating/cooling strategies
* **Checkpointing** - Save progress to recover from failures

When to Use Multiple Stages
============================

Use multi-stage AIMD when:

**Walltime constraints:**

* Single AIMD run exceeds queue limits
* Need >1000 steps total
* Testing before committing long runs

**Physical protocols:**

* Heating from 0 K → target temperature
* Cooling/annealing simulations
* Equilibration → production separation
* NPT equilibration → NVT production

**Practical workflows:**

* Check equilibration before full production
* Adjust parameters between stages based on monitoring
* Create checkpoints for long simulations

**Don't use stages if:**

* Single run fits walltime comfortably
* No temperature/pressure ramping needed
* Simple canonical ensemble (NVT) at fixed T

Designing Stage Sequences
==========================

Pattern 1: Equilibration + Production
--------------------------------------

Most common pattern - equilibrate first, then collect statistics:

**Stage 1: Equilibration (short)**

.. code-block:: python

   from teros.core.aimd import aimd_single_stage_scatter
   from aiida import orm

   # Prepare AIMD parameters
   aimd_parameters = {
       'PREC': 'Normal',
       'ENCUT': 400,
       'ISMEAR': 0,
       'SIGMA': 0.1,
       'ALGO': 'VeryFast',
       'NELM': 60,
       'LREAL': 'Auto',
       'ISYM': 0,

       # MD-specific
       'IBRION': 0,      # MD
       'SMASS': 0,       # Nose thermostat
       'POTIM': 2.0,     # 2 fs timestep
       'NSW': 100,       # 100 steps = 0.2 ps (equilibration)
       'TEBEG': 300.0,
       'TEEND': 300.0,
   }

   # Run equilibration
   equilibration = aimd_single_stage_scatter(
       slabs=slabs,
       temperature=300.0,
       steps=100,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

**Stage 2: Production (long)**

.. code-block:: python

   # Continue from equilibration
   production = aimd_single_stage_scatter(
       slabs=equilibration.outputs.structures,  # Final structures
       temperature=300.0,
       steps=500,  # 500 steps = 1.0 ps (production)
       restart_folders=equilibration.outputs.remote_folders,  # Restart
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

**Why this works:**

* Equilibration brings system to thermal equilibrium
* Production continues from equilibrated state with same electronic structure
* WAVECAR from equilibration provides good initial guess → faster convergence

Pattern 2: Temperature Ramping
-------------------------------

Gradually heat system to avoid thermal shock:

**Stage 1: Low temperature (0 → 100 K)**

.. code-block:: python

   aimd_params_stage1 = {
       **base_aimd_params,
       'NSW': 50,        # 50 steps
       'TEBEG': 0.0,     # Start at 0 K
       'TEEND': 100.0,   # Ramp to 100 K
   }

   stage1 = aimd_single_stage_scatter(
       slabs=slabs,
       temperature=50.0,  # Average temperature for logging
       steps=50,
       code=code,
       aimd_parameters=aimd_params_stage1,
       options=aimd_options,
   )

**Stage 2: Medium temperature (100 → 200 K)**

.. code-block:: python

   aimd_params_stage2 = {
       **base_aimd_params,
       'NSW': 50,
       'TEBEG': 100.0,   # Continue from Stage 1
       'TEEND': 200.0,
   }

   stage2 = aimd_single_stage_scatter(
       slabs=stage1.outputs.structures,
       temperature=150.0,
       steps=50,
       restart_folders=stage1.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_params_stage2,
       options=aimd_options,
   )

**Stage 3: Target temperature (200 → 300 K)**

.. code-block:: python

   aimd_params_stage3 = {
       **base_aimd_params,
       'NSW': 50,
       'TEBEG': 200.0,
       'TEEND': 300.0,
   }

   stage3 = aimd_single_stage_scatter(
       slabs=stage2.outputs.structures,
       temperature=250.0,
       steps=50,
       restart_folders=stage2.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_params_stage3,
       options=aimd_options,
   )

**Stage 4: Equilibration at target**

.. code-block:: python

   aimd_params_stage4 = {
       **base_aimd_params,
       'NSW': 200,
       'TEBEG': 300.0,
       'TEEND': 300.0,   # Hold at 300 K
   }

   stage4_equilibration = aimd_single_stage_scatter(
       slabs=stage3.outputs.structures,
       temperature=300.0,
       steps=200,
       restart_folders=stage3.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_params_stage4,
       options=aimd_options,
   )

**Stage 5: Production**

.. code-block:: python

   aimd_params_production = {
       **base_aimd_params,
       'NSW': 500,
       'TEBEG': 300.0,
       'TEEND': 300.0,
   }

   production = aimd_single_stage_scatter(
       slabs=stage4_equilibration.outputs.structures,
       temperature=300.0,
       steps=500,
       restart_folders=stage4_equilibration.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_params_production,
       options=aimd_options,
   )

**Rationale:**

* Gradual heating prevents structure distortion
* Each stage shorter than walltime limit
* Can monitor temperature stabilization between stages
* WAVECAR restart maintains electronic structure consistency

Pattern 3: Iterative Extension
-------------------------------

Start short, check convergence, extend if needed:

**Stage 1: Initial test (100 steps)**

.. code-block:: python

   stage1 = aimd_single_stage_scatter(
       slabs=slabs,
       temperature=300.0,
       steps=100,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

**Check convergence (manual step)**

.. code-block:: bash

   # Wait for completion
   sleep 300

   # Check if equilibrated
   verdi process show <STAGE1_PK>

   # Analyze temperature/energy drift
   verdi shell
   >>> from aiida import orm
   >>> node = orm.load_node(<CALC_PK>)
   >>> outcar = node.outputs.retrieved.get_object_content('OUTCAR')
   >>> # Parse temperature, energy from OUTCAR...

**If not equilibrated: Stage 2 (another 100 steps)**

.. code-block:: python

   stage2 = aimd_single_stage_scatter(
       slabs=stage1.outputs.structures,
       temperature=300.0,
       steps=100,  # Another 100
       restart_folders=stage1.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

**If equilibrated: Stage 2 (production)**

.. code-block:: python

   stage2_production = aimd_single_stage_scatter(
       slabs=stage1.outputs.structures,
       temperature=300.0,
       steps=500,  # Long production run
       restart_folders=stage1.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

Pattern 4: NPT → NVT Transition
--------------------------------

Equilibrate volume, then fix for production:

**Stage 1: NPT equilibration**

.. code-block:: python

   npt_parameters = {
       **base_aimd_params,
       'ISIF': 3,        # Relax cell shape/volume
       'PSTRESS': 0.0,   # Target 0 GPa
       'LANGEVIN_GAMMA': [10.0, 10.0, 10.0],  # Langevin thermostat
       'LANGEVIN_GAMMA_L': 10.0,              # Barostat coupling
       'NSW': 200,
       'PMASS': 500,     # Mass for barostat
   }

   npt_equilibration = aimd_single_stage_scatter(
       slabs=slabs,
       temperature=300.0,
       steps=200,
       code=code,
       aimd_parameters=npt_parameters,
       options=aimd_options,
   )

**Stage 2: NVT production (fixed volume)**

.. code-block:: python

   nvt_parameters = {
       **base_aimd_params,
       'ISIF': 2,        # Fixed cell
       'NSW': 500,
       # No PSTRESS/PMASS (NVT)
   }

   nvt_production = aimd_single_stage_scatter(
       slabs=npt_equilibration.outputs.structures,  # Equilibrated volume
       temperature=300.0,
       steps=500,
       restart_folders=npt_equilibration.outputs.remote_folders,
       code=code,
       aimd_parameters=nvt_parameters,
       options=aimd_options,
   )

Practical Implementation
========================

Complete Workflow Example
--------------------------

Full workflow with all stages in one script:

.. code-block:: python

   from aiida_workgraph import WorkGraph
   from teros.core.aimd import aimd_single_stage_scatter, prepare_aimd_parameters
   from aiida import orm, load_profile

   load_profile()

   # Load structures
   slabs = [orm.load_node(<PK1>), orm.load_node(<PK2>)]

   # Load code and options
   code = orm.load_code('vasp@cluster')
   aimd_options = {
       'queue_name': 'normal',
       'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 48},
       'max_wallclock_seconds': 3600,  # 1 hour per stage
   }

   # Prepare parameters using builder
   aimd_parameters = prepare_aimd_parameters(
       prec='Normal',
       encut=400,
       potim=2.0,
       lreal=True,
       algo='VeryFast',
   )

   # Build WorkGraph with all stages
   wg = WorkGraph('aimd_multi_stage')

   # Stage 1: Equilibration
   equilibration = wg.add_task(
       aimd_single_stage_scatter,
       name='equilibration',
       slabs=slabs,
       temperature=300.0,
       steps=100,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

   # Stage 2: Production 1
   production1 = wg.add_task(
       aimd_single_stage_scatter,
       name='production_1',
       slabs=equilibration.outputs.structures,
       temperature=300.0,
       steps=500,
       restart_folders=equilibration.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

   # Stage 3: Production 2 (extended)
   production2 = wg.add_task(
       aimd_single_stage_scatter,
       name='production_2',
       slabs=production1.outputs.structures,
       temperature=300.0,
       steps=500,
       restart_folders=production1.outputs.remote_folders,
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

   # Submit workflow
   wg.submit(wait=True)
   print(f"WorkGraph submitted: {wg.pk}")

**Result:**

* Total: 1100 steps = 2.2 ps (with 2 fs timestep)
* Automatic restart chaining between stages
* Each stage fits within 1 hour walltime
* Checkpoints every 500 steps

Monitoring Stages
-----------------

**During execution:**

.. code-block:: bash

   # Check workflow status
   verdi process show <WG_PK>

   # List all stages
   verdi process list -a -p 1  # Past 1 day

   # Monitor specific stage
   watch -n 60 'verdi process show <STAGE_PK>'

**After each stage:**

.. code-block:: python

   from aiida import orm

   # Load stage calculation
   stage_node = orm.load_node(<STAGE_PK>)

   # Check if successful
   if stage_node.is_finished_ok:
       print("Stage completed successfully")

       # Get final energy
       energy = stage_node.outputs.output_parameters['energy']
       print(f"Final energy: {energy:.3f} eV")

       # Get final structure
       structure = stage_node.outputs.structure
       structure.get_ase().write('stage_final.vasp')
   else:
       print(f"Stage failed: {stage_node.exit_status}")
       print(stage_node.exit_message)

**Temperature/energy analysis:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   # Extract OSZICAR data
   oszicar = stage_node.outputs.retrieved.get_object_content('OSZICAR')

   lines = [l for l in oszicar.split('\n') if 'T=' in l]
   temperatures = [float(l.split('T=')[1].split()[0]) for l in lines]
   energies = [float(l.split('E0=')[1].split()[0]) for l in lines]

   # Plot
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

   ax1.plot(temperatures)
   ax1.axhline(300, color='r', linestyle='--', label='Target')
   ax1.set_ylabel('Temperature (K)')
   ax1.legend()

   ax2.plot(energies)
   ax2.set_ylabel('Energy (eV)')
   ax2.set_xlabel('MD Step')

   plt.tight_layout()
   plt.savefig('stage_monitoring.png')

Common Issues and Solutions
============================

Issue 1: WAVECAR Restart Fails
-------------------------------

**Symptoms:**

* "Error reading WAVECAR" in VASP output
* Stage restarts from scratch despite restart_folders

**Causes:**

* WAVECAR file corrupted or incomplete
* Parameter changes incompatible with WAVECAR (ENCUT, NBANDS, KPOINTS)
* File system issues during write

**Solutions:**

.. code-block:: python

   # Make sure parameters don't change between stages
   shared_electronic = {
       'ENCUT': 400,
       'PREC': 'Normal',
       'ALGO': 'VeryFast',
       # Keep these constant!
   }

   # Stage 1
   stage1_params = {**shared_electronic, 'NSW': 100, 'TEBEG': 300.0, 'TEEND': 300.0}

   # Stage 2 - same electronic parameters
   stage2_params = {**shared_electronic, 'NSW': 500, 'TEBEG': 300.0, 'TEEND': 300.0}

**Fallback option:**

.. code-block:: python

   # If restart fails, continue without WAVECAR
   stage2_no_restart = aimd_single_stage_scatter(
       slabs=stage1.outputs.structures,  # Keep structures
       # No restart_folders = fresh electronic structure
       temperature=300.0,
       steps=500,
       code=code,
       aimd_parameters=stage2_params,
       options=aimd_options,
   )

Issue 2: Temperature Not Equilibrating
---------------------------------------

**Symptoms:**

* Temperature oscillates wildly (±50 K or more)
* Doesn't reach target temperature by end of stage

**Solutions:**

**A. Increase thermostat coupling:**

.. code-block:: python

   aimd_parameters = {
       **base_params,
       'SMASS': 0,       # Nose thermostat
       'NBLOCK': 50,     # More frequent rescaling
       'KBLOCK': 10,     # Tighter temperature control
   }

**B. Use Langevin thermostat:**

.. code-block:: python

   aimd_parameters = {
       **base_params,
       'MDALGO': 2,      # Langevin thermostat
       'LANGEVIN_GAMMA': [10.0, 10.0, 10.0],  # Friction coefficients
       'ISIF': 2,        # Fixed cell
   }

**C. Longer equilibration:**

.. code-block:: python

   # Double equilibration time
   equilibration = aimd_single_stage_scatter(
       slabs=slabs,
       temperature=300.0,
       steps=200,  # Was 100
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

Issue 3: System Drifts Between Stages
--------------------------------------

**Symptoms:**

* Energy jumps significantly when new stage starts
* Structure distorts unexpectedly between stages
* Momentum not conserved

**Cause:**

VASP doesn't preserve velocities by default when restarting.

**Solution:**

For true trajectory continuation, you need to manually handle velocity restart. PS-TEROS currently focuses on structural continuation with thermal re-equilibration.

**Workaround:**

Add short re-equilibration at start of each stage:

.. code-block:: python

   # Each production stage starts with mini-equilibration
   def stage_with_reequilibration(previous_stage, total_steps):
       # Short equilibration
       equil = aimd_single_stage_scatter(
           slabs=previous_stage.outputs.structures,
           temperature=300.0,
           steps=50,  # 50 steps to stabilize
           restart_folders=previous_stage.outputs.remote_folders,
           code=code,
           aimd_parameters=aimd_parameters,
           options=aimd_options,
       )

       # Main production
       prod = aimd_single_stage_scatter(
           slabs=equil.outputs.structures,
           temperature=300.0,
           steps=total_steps,
           restart_folders=equil.outputs.remote_folders,
           code=code,
           aimd_parameters=aimd_parameters,
           options=aimd_options,
       )

       return prod

Best Practices
==============

Stage Design Guidelines
-----------------------

**1. Stage duration**

* Equilibration: 0.2-0.5 ps (100-250 steps with 2 fs)
* Production: 1-2 ps per stage (500-1000 steps)
* Stay within 50-80% of walltime limit
* Account for startup overhead (~5-10% of walltime)

**2. Temperature protocols**

* Ramp gradually: <100 K per stage
* Hold at target for at least 0.2 ps before production
* Monitor temperature fluctuations (<±20 K is good)

**3. Parameter consistency**

Keep these constant across stages:
- ENCUT
- PREC
- ALGO
- ISMEAR/SIGMA
- LREAL
- KPOINTS

Can change between stages:
- NSW (number of steps)
- TEBEG/TEEND (temperature)
- ISIF (ensemble)
- SMASS/MDALGO (thermostat)

**4. Checkpointing strategy**

.. code-block:: python

   # Conservative: checkpoint every 100 steps
   stages = [100, 100, 100, 100, 100]  # 5 stages of 100 steps

   # Balanced: ramp up after equilibration
   stages = [100, 500, 500]  # 1 equilibration + 2 production

   # Aggressive: minimize restart overhead
   stages = [200, 1000]  # If walltime allows

**5. Error handling**

Always plan for failures:

.. code-block:: python

   # Get stage status
   if not stage1.is_finished_ok:
       print("Stage 1 failed - diagnose before continuing")
       # Don't automatically launch stage 2

   # Check equilibration quality
   if temperature_fluctuation > 30.0:  # K
       print("Not equilibrated - extend equilibration stage")

Documentation Template
----------------------

When setting up multi-stage AIMD, document your protocol:

.. code-block:: text

   # AIMD Protocol: Ag2O (111) Surface at 300 K

   ## Objective
   1.0 ps production trajectory at 300 K (NVT)

   ## Stages
   1. Equilibration: 0-0.2 ps (100 steps)
      - Purpose: Reach thermal equilibrium
      - Thermostat: Nose (SMASS=0)
      - Temperature: 300 K

   2. Production 1: 0.2-0.7 ps (250 steps)
      - Collect first half of trajectory
      - Continue from equilibration WAVECAR

   3. Production 2: 0.7-1.2 ps (250 steps)
      - Complete trajectory
      - Continue from Production 1 WAVECAR

   ## Parameters
   - Timestep: 2 fs
   - Thermostat: Nose
   - Total time: 1.2 ps
   - Checkpoints: After each stage

   ## Files
   - Configuration: aimd_multi_stage.py
   - Trajectories: stage{1,2,3}_<slab>.traj

Related Guides
==============

* :doc:`restart-calculations` - Recovering from stage failures
* :doc:`../workflows/aimd-molecular-dynamics` - Complete AIMD workflow
* :doc:`electronic-properties` - DOS/bands from AIMD trajectories

External Resources
==================

* `VASP AIMD guide <https://www.vasp.at/wiki/index.php/Category:Molecular_Dynamics>`_
* `Temperature control in VASP <https://www.vasp.at/wiki/index.php/SMASS>`_
* `Thermostats comparison <https://www.vasp.at/wiki/index.php/MDALGO>`_
