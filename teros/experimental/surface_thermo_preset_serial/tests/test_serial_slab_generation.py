#!/home/thiagotd/envs/aiida/bin/python
"""
Test Script: Serial Preset - Slab Generation and Relaxation

This script tests ONLY the slab generation and relaxation features of the
experimental serial preset module, without running the full surface
thermodynamics workflow.

Features tested:
1. Internal slab generation from miller_indices (using input bulk structure)
2. SCF calculations on generated slabs
3. Relaxation calculations on slabs
4. Relaxation energy calculation (unrelaxed - relaxed)
5. Custom builders for full parameter control
6. Concurrency control with max_number_jobs

Material: Ag2O
Miller indices: (100), (110) - generates multiple terminations
Concurrency limit: 2 jobs

Usage:
    source ~/envs/aiida/bin/activate
    python test_serial_slab_generation.py
"""

import sys
import os
from aiida import load_profile, orm
from aiida.engine import submit
from aiida_workgraph import WorkGraph
from aiida.plugins import WorkflowFactory

# Import the node builders from serial preset
from teros.experimental.surface_thermo_preset_serial.slab_operations import (
    build_scf_slabs_nodes,
    build_relax_slabs_nodes,
    build_energy_extraction_nodes,
    build_relaxation_energy_nodes,
)
from teros.experimental.surface_thermo_preset_serial.utils import (
    prepare_vasp_parameters,
    create_default_scf_parameters,
)

VaspWorkChain = WorkflowFactory('vasp.vasp')


def generate_slabs_from_bulk(bulk_filepath, miller_indices, min_slab_thickness=8, min_vacuum_thickness=10):
    """Generate slabs from bulk structure using pymatgen."""
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.surface import SlabGenerator

    # Load bulk structure
    bulk_ase = read(bulk_filepath)
    bulk_structure_input = orm.StructureData(ase=bulk_ase)

    # Convert to pymatgen
    pmg_structure = AseAtomsAdaptor.get_structure(bulk_ase)

    slab_structures = {}

    for miller in miller_indices:
        print(f"   Generating slabs for Miller index {miller}...")

        slabgen = SlabGenerator(
            pmg_structure,
            miller,
            min_slab_size=min_slab_thickness,
            min_vacuum_size=min_vacuum_thickness,
            lll_reduce=False,
            center_slab=True,
            primitive=True,
            in_unit_planes=False,
            max_normal_search=None,
        )

        slabs = slabgen.get_slabs(symmetrize=False)

        for idx, slab in enumerate(slabs):
            slab_id = f"slab_{''.join(map(str, miller))}_term_{idx}"
            slab_ase = AseAtomsAdaptor.get_atoms(slab)
            slab_structures[slab_id] = orm.StructureData(ase=slab_ase)
            print(f"      Created: {slab_id}")

    print(f"   ✓ Generated {len(slab_structures)} slab structures")
    return slab_structures


def build_slab_test_workgraph(
    slab_structures,
    code,
    potential_family,
    potential_mapping,
    kpoints_spacing,
    scf_parameters,
    relax_parameters,
    options,
    clean_workdir=False,
):
    """Build a simple workgraph to test slab SCF and relaxation."""

    wg = WorkGraph(name="test_slab_generation_relaxation")

    print("\n   Building workgraph nodes...")

    # Phase 1: SCF calculations (unrelaxed energies)
    print("   1. Adding SCF nodes for all slabs...")
    scf_nodes = build_scf_slabs_nodes(
        wg=wg,
        slabs=slab_structures,
        code=code,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        kpoints_spacing=kpoints_spacing,
        parameters=scf_parameters,
        options=options,
        clean_workdir=clean_workdir,
    )
    print(f"      ✓ Added {len(scf_nodes)} SCF nodes")

    # Extract unrelaxed energies
    print("   2. Adding energy extraction nodes for SCF...")
    unrelaxed_energy_nodes = build_energy_extraction_nodes(
        wg=wg,
        vasp_nodes=scf_nodes,
        node_type="scf",
    )
    print(f"      ✓ Added {len(unrelaxed_energy_nodes)} energy extraction nodes")

    # Phase 2: Relaxation calculations
    print("   3. Adding relaxation nodes for all slabs...")
    relax_nodes = build_relax_slabs_nodes(
        wg=wg,
        slabs=slab_structures,
        code=code,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        kpoints_spacing=kpoints_spacing,
        parameters=relax_parameters,
        options=options,
        clean_workdir=clean_workdir,
    )
    print(f"      ✓ Added {len(relax_nodes)} relaxation nodes")

    # Extract relaxed energies
    print("   4. Adding energy extraction nodes for relaxation...")
    relaxed_energy_nodes = build_energy_extraction_nodes(
        wg=wg,
        vasp_nodes=relax_nodes,
        node_type="relaxed",
    )
    print(f"      ✓ Added {len(relaxed_energy_nodes)} energy extraction nodes")

    # Phase 3: Calculate relaxation energies
    print("   5. Adding relaxation energy calculation nodes...")
    relaxation_energy_nodes = build_relaxation_energy_nodes(
        wg=wg,
        unrelaxed_energies=unrelaxed_energy_nodes,
        relaxed_energies=relaxed_energy_nodes,
    )
    print(f"      ✓ Added {len(relaxation_energy_nodes)} relaxation energy nodes")

    print("   ✓ WorkGraph construction complete")

    return wg


def main():
    """Test slab generation and relaxation with serial preset."""

    print("\n" + "="*70)
    print("TEST: SLAB GENERATION AND RELAXATION (SERIAL PRESET)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    bulk_filepath = os.path.join(structures_dir, 'ag2o.cif')

    print(f"\n2. Bulk structure:")
    print(f"   {bulk_filepath}")

    # Code configuration
    code_label = 'VASP-6.5.0@bohr-new'
    code = orm.load_code(code_label)
    potential_family = 'PBE'

    print(f"\n3. VASP Configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # VASP options
    options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # Miller indices for slab generation
    miller_indices = [
        (1, 0, 0),
        (1, 1, 0),
    ]

    print(f"\n4. Generating slabs from bulk structure...")
    print(f"   Miller indices: {miller_indices}")
    print(f"   Min slab thickness: 8 Å")
    print(f"   Min vacuum thickness: 10 Å")

    slab_structures = generate_slabs_from_bulk(
        bulk_filepath=bulk_filepath,
        miller_indices=miller_indices,
        min_slab_thickness=8,
        min_vacuum_thickness=10,
    )

    # VASP parameters - LIGHT for testing
    scf_parameters = {
        'PREC': 'Normal',
        'ENCUT': 300,
        'EDIFF': 1e-3,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 0,  # SCF only
        'ALGO': 'Fast',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    relax_parameters = {
        'PREC': 'Normal',
        'ENCUT': 300,
        'EDIFF': 1e-3,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,  # Ions only
        'NSW': 200,
        'EDIFFG': -0.5,
        'ALGO': 'Fast',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    potential_mapping = {'Ag': 'Ag', 'O': 'O'}
    kpoints_spacing = 1.0

    print(f"\n5. VASP Parameters:")
    print(f"   SCF: NSW=0, ENCUT=300, kpoints={kpoints_spacing}")
    print(f"   Relax: NSW=200, ISIF=2, ENCUT=300, kpoints={kpoints_spacing}")

    # Build workgraph
    print(f"\n6. Building test workgraph...")

    try:
        wg = build_slab_test_workgraph(
            slab_structures=slab_structures,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            kpoints_spacing=kpoints_spacing,
            scf_parameters=scf_parameters,
            relax_parameters=relax_parameters,
            options=options,
            clean_workdir=False,
        )

        print("   ✓ WorkGraph built successfully")

    except Exception as e:
        print(f"   ✗ Error building workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None

    # Set concurrent job limit
    print(f"\n7. Configuring concurrency control...")
    wg.max_number_jobs = 2
    print(f"   ✓ max_number_jobs = 2")
    print(f"   This limits concurrent VASP jobs to 2 at a time")

    # Display graph structure
    print(f"\n8. WorkGraph structure:")
    print(f"   Total slabs: {len(slab_structures)}")
    print(f"   SCF calculations: {len(slab_structures)}")
    print(f"   Relaxation calculations: {len(slab_structures)}")
    print(f"   Energy extractions: {len(slab_structures) * 2}")
    print(f"   Relaxation energy calculations: {len(slab_structures)}")
    print(f"   Total VASP jobs: {len(slab_structures) * 2}")

    # Submit
    print(f"\n9. Submitting to AiiDA daemon...")

    try:
        result = wg.submit()
        pk = result.pk if hasattr(result, 'pk') else result

        print(f"\n{'='*70}")
        print("TEST SUBMITTED SUCCESSFULLY")
        print(f"{'='*70}")
        print(f"\nWorkGraph PK: {pk}")
        print(f"\nMonitor with:")
        print(f"  verdi process show {pk}")
        print(f"  verdi process report {pk}")
        print(f"\nWatch concurrent jobs (should max out at 2):")
        print(f"  watch -n 2 'verdi process list -p 1'")

        print(f"\nExpected execution:")
        print(f"  Phase 1: SCF calculations ({len(slab_structures)} slabs)")
        print(f"           Limited to 2 concurrent VASP jobs")
        print(f"  Phase 2: Relaxation calculations ({len(slab_structures)} slabs)")
        print(f"           Limited to 2 concurrent VASP jobs")
        print(f"  Phase 3: Energy extraction and relaxation energy calculation")

        print(f"\nGenerated slabs:")
        for slab_id in sorted(slab_structures.keys()):
            print(f"  - {slab_id}")

        print(f"\nExpected outputs:")
        print(f"  For each slab:")
        print(f"    - unrelaxed_energy (from SCF)")
        print(f"    - relaxed_energy (from relaxation)")
        print(f"    - relaxed_structure")
        print(f"    - relaxation_energy (unrelaxed - relaxed)")

        print(f"\nTest Features:")
        print(f"  ✓ Internal slab generation from miller_indices")
        print(f"  ✓ Flat-graph architecture (all nodes at same level)")
        print(f"  ✓ max_number_jobs controls ALL VASP calculations")
        print(f"  ✓ No nested sub-workgraphs")
        print(f"  ✓ SCF + Relaxation workflow")
        print(f"  ✓ Relaxation energy calculation")

        print(f"\nThis is a FOCUSED test of slab operations only.")
        print(f"No bulk, reference, or thermodynamics calculations included.")
        print(f"{'='*70}\n")

        return pk

    except Exception as e:
        print(f"\n✗ Error submitting workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == '__main__':
    try:
        pk = main()
        if pk is None:
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
