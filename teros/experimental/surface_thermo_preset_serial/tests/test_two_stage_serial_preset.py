#!/usr/bin/env python3
"""
Test script for two-stage serial surface thermodynamics workflow.

Stage 1: Generate slabs from relaxed bulk structure
Stage 2: Run surface thermodynamics calculations with pre-provided slabs

This two-stage approach resolves the architectural issue where we need:
1. Slabs generated from RELAXED bulk (not input bulk)
2. Flat-graph architecture (no nested workgraphs)
3. max_number_jobs concurrency control

The solution: Generate slabs in Stage 1, then use them as input_slabs in Stage 2.
"""

import time
from pathlib import Path
from aiida import orm, load_profile

# Load AiiDA profile
load_profile()

from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph,
    generate_slabs_from_relaxed_bulk_workgraph,
)

# Configuration
structures_dir = '/home/thiagotd/git/PS-TEROS/examples/vasp/structures'
code_label = 'VASP-6.4.1@cluster06'
potential_family = 'PBE'

# VASP options for cluster02
options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 24,
    },
    'max_wallclock_seconds': 3600,
}

print("\n" + "="*80)
print("STAGE 1: Generate slabs from relaxed bulk structure")
print("="*80 + "\n")

# Stage 1: Generate slabs from relaxed bulk
stage1_wg = generate_slabs_from_relaxed_bulk_workgraph(
    structures_dir=structures_dir,
    bulk_name='ag2o.cif',
    code_label=code_label,
    potential_family=potential_family,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    bulk_options=options,

    # Slab generation parameters
    miller_indices=[(1, 0, 0)],  # Single Miller index for testing
    min_slab_thickness=8,
    min_vacuum_thickness=10,
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
)

# Submit Stage 1
stage1_wg.name = "stage1_generate_slabs"
stage1_wg.submit(wait=False)
print(f"Stage 1 submitted: PK={stage1_wg.pk}")
print(f"Monitor with: verdi process show {stage1_wg.pk}")

# Wait for Stage 1 to complete
print("\nWaiting for Stage 1 to complete...")
max_wait_time = 1800  # 30 minutes
start_time = time.time()
check_interval = 30  # Check every 30 seconds

while True:
    stage1_wg = orm.load_node(stage1_wg.pk)

    if stage1_wg.is_finished:
        break

    elapsed = time.time() - start_time
    if elapsed > max_wait_time:
        print(f"\nTimeout: Stage 1 did not complete within {max_wait_time/60:.1f} minutes")
        print(f"Current state: {stage1_wg.process_state}")
        exit(1)

    print(f"  Stage 1 still running... (elapsed: {elapsed:.0f}s, state: {stage1_wg.process_state})")
    time.sleep(check_interval)

# Check Stage 1 results
elapsed = time.time() - start_time
print(f"\nStage 1 completed in {elapsed:.1f}s")
print(f"Exit status: {stage1_wg.exit_status}")
print(f"Exit message: {stage1_wg.exit_message}")

if stage1_wg.exit_status != 0:
    print("\nERROR: Stage 1 failed!")
    print(f"Check with: verdi process report {stage1_wg.pk}")
    exit(1)

# Extract generated slabs from Stage 1 outputs
print("\nExtracting generated slabs from Stage 1 outputs...")

# The slabs are in a dynamic namespace, so we need to iterate through them
slabs_namespace = stage1_wg.outputs.slabs
slab_dict = {}

# Get all slab outputs
for key in slabs_namespace.keys():
    slab_structure = slabs_namespace[key]
    slab_dict[key] = slab_structure
    print(f"  Found slab: {key} (PK={slab_structure.pk})")

print(f"\nTotal slabs generated: {len(slab_dict)}")

if len(slab_dict) == 0:
    print("\nERROR: No slabs were generated in Stage 1!")
    exit(1)

# Also get the relaxed bulk structure
relaxed_bulk = stage1_wg.outputs.relaxed_bulk
print(f"Relaxed bulk structure: PK={relaxed_bulk.pk}")

print("\n" + "="*80)
print("STAGE 2: Run surface thermodynamics with pre-provided slabs")
print("="*80 + "\n")

# Stage 2: Run main workflow with pre-provided slabs
stage2_wg = surface_thermodynamics_serial_workgraph(
    structures_dir=structures_dir,
    bulk_name='ag2o.cif',  # Still needed for reference calculations
    code_label=code_label,
    potential_family=potential_family,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    bulk_options=options,

    # Reference materials
    metal_name='ag.cif',
    metal_potential_mapping={'Ag': 'Ag'},
    metal_options=options,
    oxygen_name='o2.cif',
    oxygen_potential_mapping={'O': 'O'},
    oxygen_options=options,

    # Slab parameters
    slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    slab_options=options,

    # Pre-provided slabs from Stage 1!
    input_slabs=slab_dict,

    # Control flags
    relax_slabs=False,  # No relaxation for quick testing
    compute_thermodynamics=True,  # Compute thermodynamics
    thermodynamics_sampling=50,  # Reduced sampling for testing
    compute_relaxation_energy=False,  # Not needed without relaxation
)

# Configure max_number_jobs for concurrency control
stage2_wg.name = "stage2_surface_thermo"
stage2_wg.max_number_jobs = 2  # Limit concurrent VASP calculations

# Submit Stage 2
stage2_wg.submit(wait=False)
print(f"Stage 2 submitted: PK={stage2_wg.pk}")
print(f"Monitor with: verdi process show {stage2_wg.pk}")

# Wait for Stage 2 to complete
print("\nWaiting for Stage 2 to complete...")
start_time = time.time()

while True:
    stage2_wg = orm.load_node(stage2_wg.pk)

    if stage2_wg.is_finished:
        break

    elapsed = time.time() - start_time
    if elapsed > max_wait_time:
        print(f"\nTimeout: Stage 2 did not complete within {max_wait_time/60:.1f} minutes")
        print(f"Current state: {stage2_wg.process_state}")
        exit(1)

    print(f"  Stage 2 still running... (elapsed: {elapsed:.0f}s, state: {stage2_wg.process_state})")
    time.sleep(check_interval)

# Check Stage 2 results
elapsed = time.time() - start_time
print(f"\nStage 2 completed in {elapsed:.1f}s")
print(f"Exit status: {stage2_wg.exit_status}")
print(f"Exit message: {stage2_wg.exit_message}")

if stage2_wg.exit_status != 0:
    print("\nERROR: Stage 2 failed!")
    print(f"Check with: verdi process report {stage2_wg.pk}")
    exit(1)

print("\n" + "="*80)
print("SUCCESS: Both stages completed successfully!")
print("="*80 + "\n")

print(f"Stage 1 (slab generation): PK={stage1_wg.pk}")
print(f"Stage 2 (thermodynamics):  PK={stage2_wg.pk}")

# Display Stage 2 outputs
print("\nStage 2 outputs:")
for key in stage2_wg.outputs.keys():
    output = stage2_wg.outputs[key]
    if isinstance(output, orm.Dict):
        print(f"  {key}: {output.get_dict()}")
    else:
        print(f"  {key}: {output}")

print("\n" + "="*80)
print("Two-stage workflow demonstration complete!")
print("="*80)
