"""
Test script for serial surface thermodynamics preset.

This script tests the flat-graph implementation where all VASP nodes
exist at the same level, allowing max_concurrent_jobs to control execution.
"""

from aiida import orm
from aiida.engine import submit
from ase.io import read

from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph
)


def main():
    """Run serial surface thermodynamics test."""

    # Configuration
    structures_dir = "/path/to/structures"  # UPDATE THIS PATH

    # Load pre-generated slabs (required for this version)
    # In practice, you would generate these first or provide them
    input_slabs = {
        'slab_100': orm.StructureData(ase=read(f"{structures_dir}/slab_100.cif")),
        'slab_110': orm.StructureData(ase=read(f"{structures_dir}/slab_110.cif")),
        'slab_111': orm.StructureData(ase=read(f"{structures_dir}/slab_111.cif")),
        'slab_210': orm.StructureData(ase=read(f"{structures_dir}/slab_210.cif")),
    }

    # Build workgraph
    wg = surface_thermodynamics_serial_workgraph.build(
        # Structure files
        structures_dir=structures_dir,
        bulk_name="bulk.cif",

        # Code and potentials
        code_label="VASP-6.4.1@cluster02",
        potential_family="PBE.54",
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},

        # K-points
        kpoints_spacing=0.4,

        # Bulk parameters
        bulk_parameters={
            'PREC': 'Accurate',
            'EDIFF': 1e-6,
            'ENCUT': 520,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 3,
            'NSW': 100,
            'LWAVE': False,
            'LCHARG': False,
        },
        bulk_options={
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },

        # Reference materials for thermodynamics
        metal_name="Ag.cif",
        oxygen_name="O2.cif",
        metal_potential_mapping={'Ag': 'Ag'},
        oxygen_potential_mapping={'O': 'O'},

        # Slab parameters
        input_slabs=input_slabs,
        slab_parameters={
            'PREC': 'Accurate',
            'EDIFF': 1e-6,
            'ENCUT': 520,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 2,
            'NSW': 100,
            'LWAVE': False,
            'LCHARG': False,
        },
        slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        slab_kpoints_spacing=0.4,

        # Control flags
        relax_slabs=True,
        compute_thermodynamics=True,
        thermodynamics_sampling=100,
        compute_relaxation_energy=True,
    )

    # CRITICAL: Set max_concurrent_jobs to test concurrency control
    wg.max_concurrent_jobs = 2

    # Submit workflow
    result = submit(wg)

    print(f"Submitted WorkGraph: {result.pk}")
    print(f"Monitor with: verdi process show {result.pk}")
    print(f"Check concurrent jobs with: verdi process list")
    print()
    print("Expected behavior:")
    print("  - Maximum 2 VASP jobs running simultaneously")
    print("  - All VASP calculations complete successfully")
    print("  - Surface energies calculated for all slabs")

    return result


if __name__ == "__main__":
    main()
