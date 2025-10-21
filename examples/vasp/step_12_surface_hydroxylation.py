#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 12: Surface Hydroxylation and Vacancy Generation

This script demonstrates the surface_hydroxylation module:
- Generate surface variants with different OH coverages
- Generate oxygen-deficient surfaces (vacancies)
- Coverage-based deduplication algorithm
- Batch VASP relaxation with controlled parallelization
- Result organization and analysis

The module enables studying:
- Hydroxylation effects on surface stability
- Oxygen vacancy formation energies
- Coverage-dependent electronic properties
- Surface reactivity under different conditions

Material: Ag2O (111) surface
Modes: hydroxylation, vacancies, or combined
Coverage: 0-100% with representative sampling

Usage:
    source ~/envs/aiida/bin/activate
    python step_12_surface_hydroxylation.py
"""

import sys
import os
from pathlib import Path
from aiida import orm, load_profile
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

def main():
    """Step 12: Test surface hydroxylation workflow."""

    print("\n" + "="*70)
    print("STEP 12: SURFACE HYDROXYLATION AND VACANCY GENERATION")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   ✓ Daemon is running")

    # Setup structure
    print("\n2. Input structure:")
    print("   OPTION 1: Use relaxed slab from previous workflow")
    print("   OPTION 2: Create test structure")

    # For this example, create a simple test structure
    # In production, you would use: structure_pk = <relaxed_slab_from_step_05>

    print("\n   Creating test Ag2O (111) surface with O adlayer...")
    from ase.build import surface
    from ase import Atom

    # Create Ag2O surface (simplified for demonstration)
    # In production, use the relaxed slab from step 5
    slab = surface('Ag', (1, 1, 1), size=(2, 2, 4), vacuum=10.0)

    # Add oxygen atoms on top (simplified O adlayer)
    import numpy as np
    z_max = max(atom.position[2] for atom in slab)
    top_atoms = [atom for atom in slab if abs(atom.position[2] - z_max) < 0.1]

    for ag_atom in top_atoms:
        o_pos = ag_atom.position.copy()
        o_pos[2] += 2.0  # Place O 2 Å above Ag
        slab.append(Atom('O', position=o_pos))

    slab.center(vacuum=10.0, axis=2)
    structure = orm.StructureData(ase=slab)

    print(f"   ✓ Test structure: {len(slab)} atoms ({sum(1 for a in slab if a.symbol == 'O')} O atoms)")
    print(f"   Structure PK: {structure.pk}")

    # Surface modification parameters
    print("\n3. Surface modification parameters:")

    # Choose mode: 'hydrogen', 'vacancies', or 'combine'
    mode = 'hydrogen'  # Hydroxylation

    surface_params = {
        'mode': mode,                    # Hydroxylation mode
        'species': 'O',                  # Target oxygen atoms
        'z_window': 0.5,                 # Surface detection window (Å)
        'which_surface': 'top',          # Modify top surface only
        'oh_dist': 0.98,                 # O-H bond distance (Å)
        'include_empty': False,          # Don't include pristine surface
        'deduplicate_by_coverage': True, # Enable deduplication
        'coverage_bins': 3,              # Sample 3 representative coverages
    }

    print(f"   Mode: {surface_params['mode']}")
    print(f"   Coverage bins: {surface_params['coverage_bins']}")
    print(f"   Expected structures: ~3-5")

    if mode == 'hydrogen':
        print(f"   → Will add OH groups to surface O atoms")
    elif mode == 'vacancies':
        print(f"   → Will remove O atoms (oxygen vacancies)")
    elif mode == 'combine':
        print(f"   → Will generate both hydroxylation and vacancies")

    # VASP configuration
    print("\n4. VASP configuration:")
    print("   Using LIGHTWEIGHT parameters for demonstration")
    print("   For production, use converged ENCUT, k-points, and tight forces")

    code_label = 'VASP-6.4.1@cluster02'  # Update to your VASP code
    print(f"   VASP code: {code_label}")

    vasp_config = {
        # Direct INCAR parameters (no translation layer)
        'parameters': {
            'PREC': 'Normal',       # Lower precision for speed
            'ENCUT': 400,           # Lower cutoff for demonstration
            'EDIFF': 1e-4,          # Looser electronic convergence
            'ISMEAR': 0,            # Gaussian smearing
            'SIGMA': 0.05,
            'ALGO': 'Fast',         # Faster algorithm
            'LREAL': 'Auto',        # Real-space projection (faster)
            'LWAVE': False,         # Don't write WAVECAR
            'LCHARG': False,        # Don't write CHGCAR
            # Relaxation parameters
            'ISIF': 2,              # Relax positions only
            'NSW': 50,              # Fewer ionic steps for demo
            'IBRION': 2,            # Conjugate gradient
            'EDIFFG': -0.05,        # Loose force convergence (eV/Å)
        },
        # K-points
        'kpoints_spacing': 0.5,     # Coarse k-points (Å⁻¹)

        # Pseudopotentials
        'potential_family': 'PBE',
        'potential_mapping': {},    # Use defaults

        # Cleanup
        'clean_workdir': False,     # Keep files for inspection
    }

    print(f"   ENCUT: {vasp_config['parameters']['ENCUT']} eV")
    print(f"   NSW (max steps): {vasp_config['parameters']['NSW']}")
    print(f"   EDIFFG (force): {vasp_config['parameters']['EDIFFG']} eV/Å")
    print(f"   K-points spacing: {vasp_config['kpoints_spacing']} Å⁻¹")

    # Scheduler options
    options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,  # Update to your cluster
        },
        'queue_name': 'normal',            # Update to your queue
        'max_wallclock_seconds': 3600 * 4, # 4 hours
    }

    # Parallelization control
    max_parallel = 2  # Process 2 structures for demonstration

    print(f"\n5. Parallelization:")
    print(f"   Batch size: {max_parallel} structures")
    print(f"   (Use higher values in production for efficiency)")

    # Build workflow
    print("\n6. Building workflow...")

    wg = build_surface_hydroxylation_workgraph(
        structure=structure,            # or structure_pk=<PK>
        surface_params=surface_params,
        code_label=code_label,
        vasp_config=vasp_config,
        options=options,
        max_parallel_jobs=max_parallel,
        name='Step12_SurfaceHydroxylation_Ag2O',
    )

    print("   ✓ Workflow built successfully")

    # Submit
    print("\n7. Submitting to AiiDA daemon...")
    result = wg.submit()
    pk = result.pk

    print(f"\n{'='*70}")
    print("STEP 12 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkflow PK: {pk}")

    print(f"\nMonitor with:")
    print(f"  verdi process show {pk}")
    print(f"  verdi process report {pk}")
    print(f"  verdi process list")

    print(f"\nWatch progress (updates every 30s):")
    print(f"  watch -n 30 verdi process show {pk}")

    print(f"\nExpected workflow steps:")
    print(f"  1. generate_structures: Creates ~{surface_params['coverage_bins']} variants")
    print(f"  2. relax_slabs_with_semaphore: VASP relaxations (batch={max_parallel})")
    print(f"  3. Returns raw outputs: manifest, structures, energies")

    print(f"\nExpected outputs (raw namespace format):")
    print(f"  - manifest: Dict with variant metadata")
    print(f"  - structures: {{idx_variantname: StructureData}}")
    print(f"  - energies: {{idx_variantname: Float}}")

    print(f"\nDescriptive output keys example:")
    print(f"  - 0_oh_000_3_7572: First structure, 3.76 OH/nm²")
    print(f"  - 1_oh_001_7_5145: Second structure, 7.51 OH/nm²")
    print(f"  (Dots replaced with underscores for AiiDA compatibility)")

    print(f"\nAfter completion, organize results:")
    print(f"  python -c \"")
    print(f"from aiida import orm")
    print(f"from teros.core.surface_hydroxylation import organize_hydroxylation_results")
    print(f"node = orm.load_node({pk})")
    print(f"results = organize_hydroxylation_results(node)")
    print(f"print('Statistics:', results['statistics'])")
    print(f"print('Successful:', len(results['successful_relaxations']), 'structures')")
    print(f"for r in results['successful_relaxations']:")
    print(f"    print(f\\\"{{r['name']}}: {{r['energy']:.6f}} eV (coverage={{r['coverage']:.2f}})\\\")\"")
    print(f"  \"")

    print(f"\nExpected runtime:")
    print(f"  - Structure generation: < 1 minute")
    print(f"  - VASP relaxations: ~10-30 minutes (depends on system size)")

    print(f"\nFor production calculations:")
    print(f"  1. Use relaxed slab from step 5: structure_pk=<PK>")
    print(f"  2. Use converged VASP parameters:")
    print(f"     - ENCUT: 520 eV (or system-appropriate)")
    print(f"     - EDIFFG: -0.02 eV/Å (tight forces)")
    print(f"     - kpoints_spacing: 0.3 Å⁻¹")
    print(f"  3. Increase coverage_bins: 10 (better sampling)")
    print(f"  4. Increase max_parallel_jobs: 5-10 (efficiency)")

    print(f"\nScientific applications:")
    print(f"  - Study hydroxylation coverage effects on stability")
    print(f"  - Calculate oxygen vacancy formation energies")
    print(f"  - Identify optimal surface configurations")
    print(f"  - Analyze coverage-dependent electronic properties")

    print(f"\n{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
