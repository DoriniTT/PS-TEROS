#!/home/thiagotd/envs/aiida/bin/python
"""
Test: Adsorption Energy with Relaxation - Ag(111) + OH

This script tests the new relaxation workflow:
1. Relax complete system (Ag + OH)
2. Separate relaxed structure
3. SCF calculations on all three components

System: Ag(111) + OH radical
"""

from aiida import load_profile, orm
from ase import Atoms, Atom
from ase.build import fcc111
from teros.core.workgraph import build_core_workgraph


def create_ag_oh_structure():
    """Create Ag(111) slab with OH adsorbate for testing."""
    # Create Ag(111) slab (4 layers, 2x2)
    slab = fcc111('Ag', size=(2, 2, 4), vacuum=10.0)

    # Manually add OH atoms on top site (directly above Ag atom)
    # Get position of first Ag atom at top layer
    top_ag_z = max(atom.position[2] for atom in slab if atom.symbol == 'Ag')
    oh_height = 2.0  # Height above top Ag layer

    # Add O atom
    o_pos = [slab[0].position[0], slab[0].position[1], top_ag_z + oh_height]
    slab.append(Atom('O', position=o_pos))

    # Add H atom (OH bond length ~0.97 Å)
    h_pos = [o_pos[0], o_pos[1], o_pos[2] + 0.97]
    slab.append(Atom('H', position=h_pos))

    return slab


def main():
    """Run Ag + OH adsorption energy calculation with relaxation."""

    print("=" * 70)
    print("ADSORPTION ENERGY WITH RELAXATION TEST")
    print("System: Ag(111) + OH")
    print("=" * 70)

    # 1. Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # 2. Create structure
    print("\n2. Creating Ag(111) + OH structure...")
    ase_structure = create_ag_oh_structure()
    structure = orm.StructureData(ase=ase_structure)

    print(f"   Formula: {structure.get_formula()}")
    print(f"   Total atoms: {len(structure.sites)}")
    print("   ✓ Structure created")

    # 3. VASP Configuration
    print("\n3. VASP Configuration:")
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # Common settings
    common_potential_mapping = {'Ag': 'Ag', 'O': 'O', 'H': 'H'}
    common_options = {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600 * 6,
        'queue_name': 'regular',
    }
    kpoints_spacing = 0.3

    # Relaxation builder inputs
    relax_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 1,
            'SIGMA': 0.2,
            'IBRION': 2,
            'NSW': 50,
            'ISIF': 2,
            'LWAVE': False,
            'LCHARG': False,
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
    }

    # SCF builder inputs
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 1,
            'SIGMA': 0.2,
            # NSW=0, IBRION=-1 set automatically
            'LWAVE': False,
            'LCHARG': False,
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
    }

    print(f"   Code: {code_label}")
    print(f"   Relaxation: ENABLED")
    print(f"   Max ionic steps: {relax_builder_inputs['parameters']['incar']['NSW']}")

    # 4. Build WorkGraph
    print("\n4. Building WorkGraph...")

    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,

        adsorption_structures={'oh_ag': structure},
        adsorption_formulas={'oh_ag': 'OH'},
        adsorption_potential_mapping=common_potential_mapping,

        relax_before_adsorption=True,
        adsorption_relax_builder_inputs=relax_builder_inputs,
        adsorption_scf_builder_inputs=scf_builder_inputs,

        name='Ag111_OH_RelaxTest',
    )

    print("   ✓ WorkGraph built")

    # 5. Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 70}")
    print("CALCULATION SUBMITTED")
    print(f"{'=' * 70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor:")
    print(f"  verdi process status {wg.pk}")
    print(f"\nExpected workflow:")
    print(f"  Phase 1: Relax Ag+OH → relaxed structure")
    print(f"  Phase 2: Separate → Ag substrate, OH molecule, complete")
    print(f"  Phase 3: SCF → 3 single-point calculations")
    print(f"\nOutputs:")
    print(f"  wg.outputs.relaxed_complete_structures['oh_ag']")
    print(f"  wg.outputs.adsorption_energies['oh_ag']")
    print(f"{'=' * 70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
