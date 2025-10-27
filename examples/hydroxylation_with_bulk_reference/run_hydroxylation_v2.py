#!/home/thiagotd/envs/aiida/bin/python
"""
Surface Hydroxylation v2 - Ag3PO4 with Bulk and Pristine Reference Calculations

This example demonstrates the updated hydroxylation workflow with automatic
bulk and pristine slab relaxations for surface energy calculations (Section S2).

NEW in v2:
- Automatic bulk relaxation from CIF file (ISIF=3 for cell relaxation)
- Automatic pristine slab relaxation (reference γ₀)
- Complete reference data for surface thermodynamics analysis
- Four new outputs: bulk_structure, bulk_energy, pristine_structure, pristine_energy

Usage:
    source ~/envs/aiida/bin/activate
    python run_hydroxylation_v2.py
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

def main():
    """Run hydroxylation workflow with bulk and pristine reference calculations."""

    print("\n" + "="*70)
    print("SURFACE HYDROXYLATION v2 - Ag3PO4 with Reference Calculations")
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

    # Load slab structure (use your actual slab structure PK or file)
    print("\n2. Loading slab structure...")
    # MODIFY THIS: Use your actual slab structure
    structure_file = Path("st2_ag3po4_110_12A.vasp")

    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        structure = orm.StructureData(ase=atoms)
        print(f"   ✓ Loaded from file: {structure_file}")
    else:
        # Or load from database
        structure_pk = 1234  # MODIFY THIS
        structure = orm.load_node(structure_pk)
        print(f"   ✓ Loaded from PK: {structure_pk}")

    print(f"   Composition: {structure.get_composition()}")
    print(f"   Structure PK: {structure.pk}")

    # Surface modification parameters
    print("\n3. Surface modification parameters:")
    surface_params = {
        'mode': 'combine',               # Hydroxylation + vacancies
        'species': 'O',                  # Target oxygen atoms
        'z_window': 0.5,                 # Surface detection window (Å)
        'which_surface': 'top',          # Modify top surface only
        'oh_dist': 0.98,                 # O-H bond distance (Å)
        'include_empty': False,          # Not needed - pristine calculated automatically
        'deduplicate_by_coverage': True, # Enable deduplication
        'coverage_bins': 3,              # Sample 3 coverages for quick test
    }

    print(f"   Mode: {surface_params['mode']}")
    print(f"   Coverage bins: {surface_params['coverage_bins']}")
    print(f"   Expected structures: ~3-5")

    # Slab VASP parameters (ISIF=2, no cell relaxation)
    print("\n4. Slab VASP Configuration (ISIF=2):")
    code_label = 'VASPGAM-6.5.0@lovelace-parexp'

    slab_builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': True,
                'LCHARG': True,
                'NCORE': 3,
                'KPAR': 4,
                'ISIF': 2,              # Ionic relaxation only
                'NSW': 500,
                'IBRION': 2,
                'EDIFFG': -0.1,
            }
        },
        'kpoints_spacing': 1.0,
        'potential_family': 'PBE',
        'potential_mapping': {
            'Ag': 'Ag',
            'P': 'P',
            'O': 'O',
            'H': 'H',
        },
        'options': {
            'resources': {
                'num_machines': 2,
                'num_cores_per_machine': 48,
            },
            'queue_name': 'parexp',
        },
        'clean_workdir': False,
    }

    print(f"   VASP code: {code_label}")
    print(f"   ISIF: {slab_builder_inputs['parameters']['incar']['ISIF']} (ionic only)")
    print(f"   ENCUT: {slab_builder_inputs['parameters']['incar']['ENCUT']} eV")

    # NEW: Bulk VASP parameters (ISIF=3, cell relaxation)
    print("\n5. Bulk VASP Configuration (ISIF=3):")

    bulk_builder_inputs = {
        'parameters': {
            'incar': {
                'PREC': 'Accurate',
                'ENCUT': 500,
                'EDIFF': 1e-6,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Fast',
                'LREAL': False,
                'LWAVE': False,
                'LCHARG': False,
                'ISIF': 3,              # Cell + ionic relaxation
                'NSW': 500,
                'IBRION': 2,
                'EDIFFG': -0.01,        # Tighter convergence
            }
        },
        'kpoints_spacing': 0.3,  # Denser k-points
        'potential_family': 'PBE',
        'potential_mapping': {
            'Ag': 'Ag',
            'P': 'P',
            'O': 'O',
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 48,
            },
            'queue_name': 'parexp',
        },
        'clean_workdir': False,
    }

    print(f"   ISIF: {bulk_builder_inputs['parameters']['incar']['ISIF']} (cell + ionic)")
    print(f"   ENCUT: {bulk_builder_inputs['parameters']['incar']['ENCUT']} eV")
    print(f"   K-points: {bulk_builder_inputs['kpoints_spacing']} Å⁻¹ (denser than slab)")

    # Parallelization and fixing
    max_parallel = 5
    fix_type = 'bottom'
    fix_thickness = 5.0

    print(f"\n6. Workflow Configuration:")
    print(f"   Batch size: {max_parallel} structures")
    print(f"   Fix type: {fix_type}")
    print(f"   Fix thickness: {fix_thickness} Å")

    # Build workflow with bulk reference
    print("\n7. Building workflow with bulk and pristine reference...")

    wg = build_surface_hydroxylation_workgraph(
        structure=structure,
        surface_params=surface_params,
        code_label=code_label,
        builder_inputs=slab_builder_inputs,

        # NEW: Bulk structure and parameters
        bulk_cif_path='ag3po4.cif',              # Option 1: CIF file
        # bulk_structure_pk=1234,                # Option 2: StructureData PK
        bulk_builder_inputs=bulk_builder_inputs,

        max_parallel_jobs=max_parallel,
        fix_type=fix_type,
        fix_thickness=fix_thickness,
        name='Ag3PO4_Hydroxylation_v2_with_references',
    )

    print("   ✓ Workflow built successfully")

    # Submit
    print("\n8. Submitting to AiiDA daemon...")
    result = wg.submit()
    pk = result.pk

    # Write PK to file
    pk_file = Path("workflow_pk.txt")
    with open(pk_file, 'w') as f:
        f.write(f"{pk}\n")
    print(f"   ✓ Workflow PK written to: {pk_file}")

    print(f"\n{'='*70}")
    print("WORKFLOW SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkflow PK: {pk}")

    print(f"\nMonitor with:")
    print(f"  verdi process show {pk}")
    print(f"  verdi process report {pk}")

    print(f"\nExpected workflow steps:")
    print(f"  1. Bulk relaxation (ISIF=3, cell optimization)")
    print(f"  2. Pristine slab relaxation (reference γ₀)")
    print(f"  3. Generate structure variants (~{surface_params['coverage_bins']})")
    print(f"  4. Relax all variants (batch={max_parallel})")

    print(f"\nNEW outputs in v2:")
    print(f"  - bulk_structure, bulk_energy")
    print(f"  - pristine_structure, pristine_energy")
    print(f"  - Plus existing: manifest, structures, energies")

    print(f"\nAfter completion, analyze results:")
    print(f"  python -c \"")
    print(f"from aiida import orm")
    print(f"from teros.core.surface_hydroxylation import organize_hydroxylation_results")
    print(f"node = orm.load_node({pk})")
    print(f"results = organize_hydroxylation_results(node)")
    print(f"ref = results['reference_data']")
    print(f"print('Bulk energy:', ref['bulk_energy'], 'eV')")
    print(f"print('Pristine energy:', ref['pristine_energy'], 'eV')")
    print(f"print('Successful variants:', len(results['successful_relaxations']))")
    print(f"  \"")

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
