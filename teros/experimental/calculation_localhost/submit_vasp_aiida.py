#!/usr/bin/env python
"""
AiiDA submission script for VASP calculation
This script sets up and submits a VASP calculation using AiiDA
"""

from aiida import orm, load_profile
from aiida.engine import submit
from aiida.plugins import CalculationFactory, DataFactory
import sys

# Load the AiiDA profile
load_profile()

# Get the VASP calculation class
VaspCalculation = CalculationFactory('vasp.vasp')

# Get data types
StructureData = DataFactory('core.structure')
KpointsData = DataFactory('core.array.kpoints')
PotcarData = DataFactory('vasp.potcar')

def create_structure():
    """Create Si structure from POSCAR"""
    # Si FCC structure with lattice parameter 5.43 Å
    structure = StructureData(cell=[
        [0.0, 2.715, 2.715],
        [2.715, 0.0, 2.715],
        [2.715, 2.715, 0.0]
    ])
    structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Si')
    structure.append_atom(position=(1.35825, 1.35825, 1.35825), symbols='Si')  # 0.25 * 5.43
    return structure

def create_kpoints():
    """Create k-points mesh"""
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([1, 1, 1], offset=[0.0, 0.0, 0.0])
    return kpoints

def load_potcar_for_element(element):
    """Load POTCAR from file for specific element"""
    potcar_file = '/home/thiagotd/git/PS-TEROS/teros/experimental/calculation_localhost/POTCAR'

    try:
        # Use get_or_create_from_file to create PotcarData node
        potcar, created = PotcarData.get_or_create_from_file(potcar_file)
        if created:
            print(f"Created new POTCAR node with PK: {potcar.pk}")
        else:
            print(f"Using existing POTCAR node with PK: {potcar.pk}")
        return potcar
    except Exception as e:
        print(f"Error loading POTCAR: {e}")
        print("You may need to upload POTCAR families first")
        return None

def get_potcar_mapping(structure):
    """Create POTCAR mapping for structure"""
    # Get unique elements from structure
    elements = structure.get_kind_names()

    potcars = {}
    for element in elements:
        potcar = load_potcar_for_element(element)
        if potcar is None:
            raise ValueError(f"Could not load POTCAR for {element}")
        potcars[element] = potcar

    return potcars

def setup_calculation():
    """Set up the VASP calculation"""

    # Get the code (replace with your actual code label)
    try:
        code = orm.load_code('vasp-6.5.1-std@localhost')
    except Exception as e:
        print(f"Error loading code: {e}")
        print("Available codes:")
        for code in orm.QueryBuilder().append(orm.Code).all(flat=True):
            print(f"  - {code.label}@{code.computer.label}")
        sys.exit(1)

    # Create structure
    structure = create_structure()
    print(f"Created structure with {len(structure.sites)} atoms")

    # Create k-points
    kpoints = create_kpoints()
    print("Created k-points mesh: 1x1x1")

    # INCAR parameters
    incar_params = orm.Dict(dict={
        'SYSTEM': 'Si Test via AiiDA',
        'ENCUT': 240,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'EDIFF': 1e-4,
        'NSW': 0,
        'IBRION': -1,
        'LCHARG': False,
        'LWAVE': False,
    })

    # Set up calculation builder
    builder = VaspCalculation.get_builder()
    builder.code = code
    builder.structure = structure
    builder.kpoints = kpoints
    builder.parameters = incar_params

    # Handle POTCAR - load from file
    try:
        potcars = get_potcar_mapping(structure)
        # The builder.potential is a namespace that expects PotcarData nodes directly
        for kind_name, potcar in potcars.items():
            builder.potential[kind_name] = potcar
        print(f"Loaded POTCARs for: {list(potcars.keys())}")
    except Exception as e:
        print(f"Warning: Could not load POTCARs from file: {e}")
        print("Attempting to use POTCAR family approach...")
        # Fallback: use potential_family if available
        # This requires POTCARs to be uploaded as a family
        # Uncomment and modify as needed:
        # builder.potential_family = orm.Str('PBE')
        # builder.potential_mapping = orm.Dict(dict={'Si': 'Si'})

    # Metadata for job submission
    builder.metadata.description = 'Si bulk test calculation via AiiDA'
    builder.metadata.options.resources = {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 1
    }
    builder.metadata.options.max_wallclock_seconds = 3600  # 1 hour
    builder.metadata.options.withmpi = True

    return builder

def main():
    """Main execution"""
    print("Setting up VASP calculation...")

    builder = setup_calculation()

    # Dry run option
    if '--dry-run' in sys.argv:
        print("\nDry run mode - not submitting")
        print(f"Code: {builder.code.label}@{builder.code.computer.label}")
        print(f"Structure: {builder.structure.get_formula()}")
        print(f"INCAR parameters: {builder.parameters.get_dict()}")
        return

    # Submit the calculation
    print("\nSubmitting calculation...")
    node = submit(builder)
    print(f"Submitted calculation with PK: {node.pk}")
    print(f"Monitor with: verdi process status {node.pk}")
    print(f"Check list: verdi process list")

if __name__ == '__main__':
    main()
