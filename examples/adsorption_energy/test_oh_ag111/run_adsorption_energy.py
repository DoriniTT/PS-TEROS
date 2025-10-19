#!/usr/bin/env python
"""
Example: Calculate adsorption energy of OH on Ag(111) surface.

This example demonstrates the adsorption energy module workflow:
1. Load substrate+adsorbate structures
2. Separate into substrate, molecule, complete
3. Relax all three systems with VASP
4. Calculate adsorption energy
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from aiida_workgraph import WorkGraph

# Import PSTEROS core functions
from teros.core import (
    get_structure_from_file,
    compute_adsorption_energies_scatter,
)


def main():
    """Run adsorption energy calculation workflow."""
    load_profile()

    print("=" * 70)
    print("Adsorption Energy Calculation: OH on Ag(111)")
    print("=" * 70)
    print()

    # Configuration
    base_path = Path(__file__).parent

    # Load VASP code (adjust for your setup)
    code_label = 'vasp@localhost'  # Modify to match your VASP code
    try:
        code = orm.load_code(code_label)
    except:
        print(f"ERROR: Could not load VASP code '{code_label}'")
        print("Please create a VASP code in AiiDA first:")
        print(f"  verdi code create core.code.installed")
        sys.exit(1)

    # Load structures
    # For this example, we'll create a simple test structure programmatically
    # In real use, you'd load from CIF files using get_structure_from_file()

    from pymatgen.core import Structure, Lattice
    from ase import Atoms

    # Create simple Ag slab with OH
    lattice = Lattice.from_parameters(a=5.8, b=5.8, c=20.0,
                                     alpha=90, beta=90, gamma=90)

    ag_positions = [[0.0, 0.0, 10.0], [2.9, 0.0, 10.0],
                    [0.0, 2.9, 10.0], [2.9, 2.9, 10.0]]
    oh_positions = [[1.45, 1.45, 12.5], [1.45, 1.45, 13.5]]

    species = ['Ag'] * 4 + ['O', 'H']
    positions = ag_positions + oh_positions

    structure = Structure(lattice, species, positions, coords_are_cartesian=True)

    # Prepare input data
    structures = {
        'site1': orm.StructureData(pymatgen=structure),
    }

    adsorbate_formulas = {
        'site1': 'OH',
    }

    # VASP parameters (minimal for testing, adjust for production)
    parameters = {
        'ENCUT': 400,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'NSW': 50,
        'ISIF': 2,  # Relax ions only, keep cell fixed
        'LWAVE': False,
        'LCHARG': False,
    }

    # Scheduler options
    options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 4,
        },
        'max_wallclock_seconds': 3600,
        'queue_name': 'debug',  # Adjust for your cluster
    }

    # Potential mapping
    potential_mapping = {
        'Ag': 'Ag',
        'O': 'O',
        'H': 'H',
    }

    print("Setting up WorkGraph...")

    # Create WorkGraph
    wg = WorkGraph('adsorption_energy_oh_ag111')

    # Add adsorption energy calculation task
    ads_task = wg.add_task(
        compute_adsorption_energies_scatter,
        name='compute_adsorption_energies',
        structures=structures,
        adsorbate_formulas=adsorbate_formulas,
        code=code,
        potential_family='PBE',  # Adjust for your setup
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    print("Submitting WorkGraph to AiiDA daemon...")
    wg.submit(wait=False)

    print()
    print("=" * 70)
    print(f"WorkGraph submitted! PK: {wg.pk}")
    print()
    print("Monitor progress with:")
    print(f"  verdi process show {wg.pk}")
    print(f"  verdi process report {wg.pk}")
    print()
    print("After completion, check results with:")
    print(f"  verdi process show {wg.pk}")
    print("=" * 70)

    return wg.pk


if __name__ == '__main__':
    pk = main()
    sys.exit(0)
