=======================================
How to Restart Failed Calculations
=======================================

This guide shows you how to diagnose and recover from failed calculations in PS-TEROS workflows.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

Calculations can fail for many reasons: convergence issues, walltime limits, node crashes, or incorrect parameters. PS-TEROS provides several mechanisms to restart and recover from failures without losing work.

**Key concepts:**

* **Failed nodes** - AiiDA processes that didn't complete successfully
* **Restart folders** - Remote directories containing WAVECAR/CHGCAR for continuation
* **Scatter-gather resilience** - Failed individual slabs don't block entire workflow
* **Parameter adjustment** - Modifying settings based on failure diagnosis

When to Use This Guide
======================

Use this guide when:

* Calculations exceed walltime limits
* Electronic convergence fails (EDIFF not reached)
* Ionic convergence fails (EDDDG not reached)
* Node crashes or network issues interrupt calculations
* You need to adjust parameters mid-workflow

Diagnosing Failures
===================

Step 1: Identify Failed Processes
----------------------------------

Use AiiDA's ``verdi`` commands to check workflow status:

.. code-block:: bash

   # Check overall workflow status
   verdi process show <PK>

   # Get detailed report
   verdi process report <PK>

   # List all failed processes
   verdi process list -S FAILED

**Example output:**

.. code-block:: text

   PK  Created    State        Process Label
   123 1h ago     ⏹ Failed     VaspCalculation
   124 1h ago     ⏹ Failed     VaspCalculation
   125 1h ago     ✓ Finished   VaspCalculation

Two calculations failed while one succeeded. The scatter-gather pattern allows the workflow to continue with successful calculations.

Step 2: Examine Error Messages
-------------------------------

Get detailed information about a specific failed calculation:

.. code-block:: bash

   # Show calculation details
   verdi process show 123

   # View calculation logs
   verdi calcjob outputls 123  # List output files
   verdi calcjob outputcat 123 stdout  # Read VASP output
   verdi calcjob outputcat 123 OUTCAR  # Read OUTCAR

**Common error patterns:**

* **"EDIFF not reached"** → Increase ``NELM`` or adjust ``AMIX``/``BMIX``
* **"ZBRENT: fatal error"** → Electronic structure problem, try different ``ALGO``
* **"Walltime exceeded"** → Increase time limit or adjust ``NPAR``/``KPAR``
* **"POSMAP: wrong number of atoms"** → Structure inconsistency, check input

Step 3: Check Remote Calculations
----------------------------------

Sometimes you need to inspect files on the compute cluster:

.. code-block:: bash

   # Show remote folder location
   verdi calcjob remotecat 123 --path  # Get remote path

   # Copy files locally for inspection
   verdi calcjob gotocomputer 123  # SSH to compute node

Common Failure Scenarios
========================

Scenario 1: Electronic Convergence Failure
-------------------------------------------

**Symptoms:**

* OUTCAR shows "EDIFF not reached in NELM steps"
* Calculation stops before ionic relaxation completes

**Diagnosis:**

Check convergence history in OUTCAR:

.. code-block:: bash

   verdi calcjob outputcat <PK> OUTCAR | grep "DAV:" | tail -20

If energy oscillates or changes slowly, electronic convergence is difficult.

**Solution:**

Adjust electronic convergence parameters:

.. code-block:: python

   # Original parameters (too strict or poor convergence)
   bulk_parameters = {
       'NELM': 60,    # Too few steps
       'ALGO': 'Fast',  # May not converge for difficult systems
       'EDIFF': 1e-6,   # Very strict
   }

   # Improved parameters
   bulk_parameters = {
       'NELM': 100,      # More electronic steps
       'ALGO': 'Normal',  # More robust algorithm
       'EDIFF': 1e-5,     # Slightly relaxed (still accurate)
       'AMIX': 0.2,       # Reduce mixing for stability
       'BMIX': 0.0001,    # Reduce magnetic mixing
   }

**Restart the workflow** with new parameters using the same structure files.

Scenario 2: Walltime Limit Exceeded
------------------------------------

**Symptoms:**

* Calculation status shows FAILED with "walltime exceeded"
* OUTCAR ends mid-iteration

**Diagnosis:**

Check how far the calculation progressed:

.. code-block:: bash

   verdi calcjob outputcat <PK> OUTCAR | grep "Iteration"

**Solution A: Increase Walltime**

Modify job options:

.. code-block:: python

   # Original (too short)
   bulk_options = {
       'queue_name': 'normal',
       'max_wallclock_seconds': 3600,  # 1 hour
       'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 24},
   }

   # Extended walltime
   bulk_options = {
       'queue_name': 'normal',
       'max_wallclock_seconds': 14400,  # 4 hours
       'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 24},
   }

**Solution B: Speed Up Convergence**

Optimize VASP parallelization:

.. code-block:: python

   bulk_parameters = {
       # ... other parameters ...
       'NPAR': 4,   # Parallelize over bands (tune based on cores)
       'LREAL': 'Auto',  # Real-space projection (faster for large cells)
   }

**Solution C: Use Restart Folders**

For AIMD or long relaxations, chain calculations using restart folders:

.. code-block:: python

   # Stage 1: Initial run
   stage1 = aimd_single_stage_scatter(
       slabs=slabs,
       steps=100,  # Partial run
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

   # Stage 2: Continue from Stage 1
   stage2 = aimd_single_stage_scatter(
       slabs=stage1.outputs.structures,  # Final structures
       steps=100,  # Another 100 steps
       restart_folders=stage1.outputs.remote_folders,  # Continue
       code=code,
       aimd_parameters=aimd_parameters,
       options=aimd_options,
   )

This pattern works for AIMD and can be adapted for difficult relaxations.

Scenario 3: Ionic Convergence Issues
-------------------------------------

**Symptoms:**

* Electronic steps converge but structure won't relax
* "EDDDG not reached" warnings in OUTCAR
* Forces remain large after many ionic steps

**Diagnosis:**

Check forces and stress:

.. code-block:: bash

   verdi calcjob outputcat <PK> OUTCAR | grep -A 10 "TOTAL-FORCE"

**Solution:**

Adjust ionic convergence and dynamics:

.. code-block:: python

   # Original (may be too strict or poor algorithm)
   slab_parameters = {
       'EDIFFG': -0.01,  # Very strict force criterion
       'IBRION': 2,      # CG algorithm
       'POTIM': 0.5,     # Large step size
   }

   # Improved
   slab_parameters = {
       'EDIFFG': -0.02,   # Relaxed criterion (still good)
       'IBRION': 1,       # RMM-DIIS (more robust)
       'POTIM': 0.3,      # Smaller, more stable steps
       'NSW': 200,        # More ionic steps allowed
   }

Scenario 4: Memory or Disk Issues
----------------------------------

**Symptoms:**

* Job crashes without clear error
* "No space left on device" in stderr
* Calculation hangs indefinitely

**Solution:**

.. code-block:: python

   # Reduce memory/disk usage
   bulk_parameters = {
       'LWAVE': False,   # Don't write WAVECAR (saves disk)
       'LCHARG': False,  # Don't write CHGCAR (saves disk)
       'NPAR': 4,        # Better memory distribution
   }

   # For slabs (if you need to restart later)
   slab_parameters = {
       'LWAVE': True,    # Keep WAVECAR for potential restart
       'LCHARG': False,  # Save space (CHGCAR less critical)
   }

Restarting Workflows
====================

Option 1: Restart Individual Calculations
------------------------------------------

If only a few slabs failed, you can restart just those:

**Step 1: Identify failed slabs**

.. code-block:: python

   from aiida import load_profile, orm
   load_profile()

   # Load your workflow
   wg_node = orm.load_node(<PK>)

   # Find failed slab calculations
   from aiida.cmdline.utils.query import calc_states
   failed_calcs = orm.QueryBuilder().append(
       orm.WorkChainNode,
       filters={'uuid': wg_node.uuid}
   ).append(
       orm.CalcJobNode,
       filters={'attributes.exit_status': {'!==': 0}},
       edge_filters={'label': {'like': 'slab_%'}}
   ).all()

   print(f"Found {len(failed_calcs)} failed calculations")

**Step 2: Create new workflow with only failed slabs**

Extract the failed slab structures and rerun:

.. code-block:: python

   # Get structures from failed calculations
   failed_structures = []
   for calc in failed_calcs:
       structure = calc[0].inputs.structure
       failed_structures.append(structure)

   # Build new workflow with adjusted parameters
   retry_wg = build_core_workgraph(
       structures=failed_structures,  # Only failed slabs
       bulk_name=bulk_filename,
       code_label=code_label,
       potential_family=potential_family,

       # Adjusted parameters based on diagnosis
       slab_parameters=improved_slab_parameters,
       slab_options=extended_slab_options,

       # Same settings as original
       compute_thermodynamics=True,
       thermodynamics_sampling=10,
   )

   retry_wg.submit(wait=True)

Option 2: Restart Entire Workflow
----------------------------------

If many calculations failed due to systematic parameter issues:

.. code-block:: python

   # Rebuild entire workflow with corrected parameters
   wg = build_core_workgraph(
       structures_dir=structures_dir,
       bulk_name=bulk_filename,
       code_label=code_label,
       potential_family=potential_family,

       # Corrected parameters
       bulk_parameters=corrected_bulk_parameters,
       slab_parameters=corrected_slab_parameters,
       bulk_options=extended_bulk_options,
       slab_options=extended_slab_options,

       # Same feature flags
       compute_relaxation_energy=True,
       compute_cleavage=True,
       compute_thermodynamics=True,
       thermodynamics_sampling=10,
   )

   wg.submit(wait=True)

PS-TEROS will recompute all calculations. If you're using AiiDA's caching and haven't changed parameters for successful calculations, those may be retrieved from cache.

Option 3: Manual Continuation (Advanced)
-----------------------------------------

For very long calculations, you can manually continue from checkpoint files:

.. code-block:: python

   from aiida import orm

   # Get remote folder from failed calculation
   failed_calc = orm.load_node(<FAILED_PK>)
   restart_folder = failed_calc.outputs.remote_folder

   # Create new calculation continuing from this folder
   # (This is what AIMD restart_folders parameter does automatically)

See the :doc:`aimd-stages` guide for detailed examples of checkpoint/restart patterns.

Best Practices
==============

Prevention Strategies
---------------------

**1. Start Conservative**

Use relaxed convergence criteria first, then tighten if needed:

.. code-block:: python

   # Phase 1: Quick test
   test_parameters = {
       'ENCUT': 400,      # Lower cutoff
       'EDIFF': 1e-4,     # Relaxed electronic
       'EDIFFG': -0.05,   # Relaxed ionic
   }

   # Phase 2: Production (after confirming convergence)
   production_parameters = {
       'ENCUT': 520,
       'EDIFF': 1e-6,
       'EDIFFG': -0.02,
   }

**2. Test Before Large Workflows**

Run a single slab calculation before launching scatter-gather over all slabs:

.. code-block:: python

   # Test with one slab first
   from teros.core.slabs import generate_slabs_scatter
   test_slab = generate_slabs_scatter(
       structures=[test_structure],  # Just one
       max_index=1,
       min_slab_size=15,
   )

   # Check if this succeeds before full run

**3. Use Appropriate Resources**

Match computational resources to problem size:

.. code-block:: python

   # Small slabs (<50 atoms)
   small_options = {
       'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 24},
       'max_wallclock_seconds': 3600,  # 1 hour
   }

   # Large slabs (>100 atoms)
   large_options = {
       'resources': {'num_machines': 2, 'num_mpiprocs_per_machine': 48},
       'max_wallclock_seconds': 14400,  # 4 hours
   }

**4. Monitor Early**

Check workflow progress early:

.. code-block:: bash

   # Launch workflow
   python my_workflow.py

   # Wait 5-10 minutes for first calculations
   sleep 600

   # Check for early failures
   verdi process list -S FAILED -p 1  # Past 1 day

If failures appear early, diagnose and fix before wasting compute time.

Troubleshooting Checklist
--------------------------

When calculations fail, work through this checklist:

1. ☐ Run ``verdi process report <PK>`` for error details
2. ☐ Check OUTCAR for convergence history
3. ☐ Verify structure makes sense (reasonable lattice parameters, no overlapping atoms)
4. ☐ Compare parameters with successful examples in ``examples/`` directory
5. ☐ Check if similar system worked with different parameters
6. ☐ Verify pseudopotential family is correct (``verdi data vasp-potcar listfamilies``)
7. ☐ Ensure sufficient walltime and resources
8. ☐ Test with single calculation before large scatter-gather

Related Guides
==============

* :doc:`custom-slabs` - Creating slab structures that converge well
* :doc:`aimd-stages` - Sequential staging patterns for long calculations
* :doc:`electronic-properties` - Convergence considerations for DOS/bands

For deeper debugging, see also:

* `AiiDA troubleshooting <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/faq.html>`_
* `VASP error messages <https://www.vasp.at/wiki/index.php/Category:Errors>`_
