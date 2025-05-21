"""
Runner script for TEROS DFT workflow using VASP.

This script combines the builder configuration and workflow execution into a single file.
It sets up the VASP calculation parameters and runs the create_teros_workgraph function for surface thermodynamics.
"""
import os
import sys
from ase.io import read
from aiida.orm import Dict, Bool, Str, Int, Float, List, StructureData, KpointsData, load_code
from aiida.plugins import WorkflowFactory
from aiida import load_profile

# Import the create_teros_workgraph function
from teros import create_teros_workgraph

# Load AiiDA profile
load_profile()

#########################################################
# VASP CONFIGURATION - Edit parameters below as needed
#########################################################

# Workflow name
WORKGRAPH_NAME = "teros_vasp"

# Import the VASP workchain
DFTWORKCHAIN = WorkflowFactory('vasp.vasp')

# Define DFT code to use
CODE = 'VASP'

# Computational resources
CODE_LABEL = "VASP-6.4.3@bohr"  # Adjust to your VASP code label configured in AiiDA (e.g., vasp@mycluster)
NCORE = 2
NUM_CORES = 40
NUM_MACHINES = 1
QUEUE_NAME = 'par40'
MAX_WALLTIME = 3 * 24 * 60 * 60  # 3 days
MAX_ITERATIONS = 5

# Base directory for file paths
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# File paths
STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/bulk/Ag6O8P2_optimized.cif")
AG_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/Ag.cif")
P_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/P.cif")
O2_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/O2.cif")

# VASP potential settings
POTENTIAL_FAMILY = "PBE"
POTENTIAL_MAPPING = {
    'Ag': 'Ag',
    'P': 'P',
    'O': 'O',
}

# Electronic parameters
ALGO = "Fast"
PREC = "Accurate"
ENCUT = 500
KSPACING = 0.3
ISMEAR = 0
SIGMA = 0.01
ISPIN = 2

# Bulk relaxation parameters
ISIF = 3
IBRION = 2
NSW = 500
EDIFFG = -0.01
EDIFF = 1.0e-5
KPOINTS = KpointsData()
KPOINTS.set_kpoints_mesh([3, 3, 3])

# Slab relaxation parameters
SLAB_ISIF = 2
SLAB_EDIFFG = -0.05
SLAB_DIPOLE = False
LDIPOL = True
IDIPOL = 3
KPOINTS_SLAB = KpointsData()
KPOINTS_SLAB.set_kpoints_mesh([2, 2, 1])

# Load VASP code
vasp_code = load_code(CODE_LABEL)

# Set resource parameters based on code/computer label
# This example distinguishes between a generic "cluster" setup and a local "bohr" setup.
# Adjust this logic based on your AiiDA computer configuration.
resource_kwargs = {}
if "cluster" in CODE_LABEL.lower():  # Example: if your code label contains "cluster"
    resource_kwargs["num_mpiprocs_per_machine"] = 16 # Typical for clusters using MPI
    queue_name = None # Or specify your queue if not the default
else:  # Default or "bohr" (assumed local or specific setup)
    resource_kwargs["num_cores_per_machine"] = NUM_CORES # For machines where you specify total cores
    queue_name = QUEUE_NAME

# Set the slab parameters
MILLER_INDICES = (1, 0, 0)   # Miller index for surface orientation
MIN_VACUUM_THICKNESS = 15.0  # Minimum vacuum thickness in Angstroms
MIN_SLAB_THICKNESS = 15.0    # Minimum thickness of the slab in Angstroms
LLL_REDUCE = True            # Whether to reduce the cell using LLL algorithm
CENTER_SLAB = True           # Whether to center the slab in the simulation cell
SYMMETRIZE = True            # Whether to generate symmetrically distinct slab terminations
PRIMITIVE = True             # Whether to use the primitive cell for slab generation
MAX_NORMAL_SEARCH = None     # Max normal search for pymatgen slab generation (None for default behavior).
IN_UNIT_PLANES = False       # Whether to restrict slab generation to unit planes.
SAMPLING = 300 # Number of sampling points for the surface phase diagram (e.g., 300 results in a 300x300 grid).
               # Use a smaller number (e.g., 50-100) for faster testing, especially for initial runs.
               # 300 usually provides a good resolution for publication-quality diagrams.

# Option for manually created slab structures
# Provide paths to your pre-existing slab structure files (e.g., CIF, POSCAR).
# These are examples; uncomment and modify paths if you have your own slab files.
MANUAL_SLABS_PATHS = [
    # os.path.join(BASE_DIR, "structures/slabs/MySlab_term1.cif"),
    # os.path.join(BASE_DIR, "structures/slabs/MySlab_term2.cif"),
]

# Choose how to handle slabs:
# Set to True to use manual slabs, False to use automatic generation
USE_MANUAL_SLABS = len(MANUAL_SLABS_PATHS) > 0

def load_manual_slabs():
    """
    Load manually created slabs from files into AiiDA StructureData objects
    
    Returns:
        dict: Dictionary mapping slab identifiers to StructureData objects
    """
    manual_slabs = {}
    
    if not MANUAL_SLABS_PATHS:
        return None
        
    for i, path in enumerate(MANUAL_SLABS_PATHS):
        if os.path.exists(path):
            # Load the structure from file and store as AiiDA StructureData
            structure = StructureData(ase=read(path))
            manual_slabs[f"s_{i}"] = structure
            print(f"Loaded manual slab {i+1}/{len(MANUAL_SLABS_PATHS)} from {path}")
        else:
            print(f"Warning: Slab file not found at {path}")
            
    return manual_slabs if manual_slabs else None

#########################################################
# BUILDER FUNCTIONS
#########################################################

def builder_bulk_relax():
    """Set up a VASP bulk relaxation calculation"""
    builder = DFTWORKCHAIN.get_builder()
    
    # Load structure from file
    structure = StructureData(ase=read(STRUCTURE_PATH))
    builder.structure = structure
    
    # Set VASP code
    builder.code = vasp_code
    
    # Set pseudopotentials
    builder.potential_family = Str(POTENTIAL_FAMILY)
    builder.potential_mapping = Dict(dict=POTENTIAL_MAPPING)
    
    # Define parameters for electronic structure and ionic relaxation
    parameters = {
        'incar': {
            # Electronic parameters
            'NCORE': NCORE,
            'ALGO': ALGO,
            'PREC': PREC,
            'ENCUT': ENCUT,
            'ISPIN': ISPIN,
            'ISMEAR': ISMEAR,
            'SIGMA': SIGMA,
            'NELM': 60,
            'NELMIN': 6,
            'LREAL': 'Auto',
            'EDIFF': EDIFF,
            
            # Ionic relaxation parameters
            'ISIF': ISIF,
            'IBRION': IBRION,
            'NSW': NSW,
            'EDIFFG': EDIFFG,
        },
    }
    builder.parameters = Dict(dict=parameters)
    
    # Set kpoints
    builder.kpoints = KPOINTS
    
    # Configure computational resources
    options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
    }
    if queue_name:
        options['queue_name'] = queue_name
        
    builder.options = Dict(dict=options)
    
    # Set parser settings
    settings = {
        'parser_settings': {
            'add_energies': True,
            'add_forces': True,
            'add_stress': True,
        }
    }
    builder.settings = Dict(dict=settings)
    
    # Set restart settings
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    
    # Set meta info
    builder.metadata.label = "VASP bulk structure relaxation"
    builder.metadata.description = "VaspWorkChain calculation for bulk optimization"
    
    return builder

def builder_slab_relax():
    """Set up a VASP slab relaxation calculation"""
    builder = DFTWORKCHAIN.get_builder()
    
    # Set VASP code
    builder.code = vasp_code
    
    # Set pseudopotentials
    builder.potential_family = Str(POTENTIAL_FAMILY)
    builder.potential_mapping = Dict(dict=POTENTIAL_MAPPING)
    
    # Define parameters for electronic structure and ionic relaxation with slab-specific settings
    parameters = {
        'incar': {
            # Electronic parameters
            'NCORE': NCORE,
            'PREC': PREC,
            'ALGO': ALGO,
            'ENCUT': ENCUT,
            'ISMEAR': ISMEAR,
            'ISPIN': ISPIN,
            'SIGMA': SIGMA,
            'NELM': 60,
            'NELMIN': 6,
            'LREAL': 'Auto',
            'EDIFF': EDIFF,
            
            # Ionic relaxation parameters for slabs
            'ISIF': SLAB_ISIF,
            'IBRION': IBRION,
            'NSW': NSW,
            'EDIFFG': SLAB_EDIFFG,
        },
    }
    
    # Add dipole correction for slabs
    if SLAB_DIPOLE:
        parameters['incar'].update({
            'LDIPOL': LDIPOL,
            'IDIPOL': IDIPOL,
        })
    
    builder.parameters = Dict(dict=parameters)
    
    # Set kpoints
    builder.kpoints = KPOINTS_SLAB
    
    # Configure computational resources
    options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
    }
    if queue_name:
        options['queue_name'] = queue_name
        
    builder.options = Dict(dict=options)
    
    # Set parser settings
    settings = {
        'parser_settings': {
            'add_energies': True,
            'add_forces': True,
            'add_stress': True,
        }
    }
    builder.settings = Dict(dict=settings)
    
    # Set restart settings
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    
    # Set meta info
    builder.metadata.label = "VASP slab structure relaxation"
    builder.metadata.description = "VaspWorkChain calculation for slab with fixed cell"
    
    return builder

def builder_ag_relax():
    """Set up a VASP relaxation calculation for pure silver"""
    builder = DFTWORKCHAIN.get_builder()
    
    # Load structure from file
    structure = StructureData(ase=read(AG_STRUCTURE_PATH))
    builder.structure = structure
    
    # Set VASP code
    builder.code = vasp_code
    
    # Set pseudopotentials
    builder.potential_family = Str(POTENTIAL_FAMILY)
    builder.potential_mapping = Dict(dict=POTENTIAL_MAPPING)
    
    # Define parameters for electronic structure and ionic relaxation
    parameters = {
        'incar': {
            # Electronic parameters
            'NCORE': NCORE,
            'ALGO': ALGO,
            'PREC': PREC,
            'ENCUT': ENCUT,
            'ISPIN': ISPIN,
            'ISMEAR': 1,     # Metal-specific smearing for Ag
            'SIGMA': 0.1,    # Wider smearing for metal
            'NELM': 60,
            'NELMIN': 6,
            'LREAL': 'Auto',
            'EDIFF': EDIFF,
            
            # Ionic relaxation parameters
            'ISIF': ISIF,
            'IBRION': IBRION,
            'NSW': NSW,
            'EDIFFG': EDIFFG,
        },
    }
    builder.parameters = Dict(dict=parameters)
    
    # Set kpoints (dense mesh for metals)
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([7, 7, 7])
    builder.kpoints = kpoints
    
    # Configure computational resources
    options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
    }
    if queue_name:
        options['queue_name'] = queue_name
        
    builder.options = Dict(dict=options)
    
    # Set parser settings
    settings = {
        'parser_settings': {
            'add_energies': True,
            'add_forces': True,
            'add_stress': True,
        }
    }
    builder.settings = Dict(dict=settings)
    
    # Set restart settings
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    
    # Set meta info
    builder.metadata.label = "VASP pure Ag relaxation"
    builder.metadata.description = "VaspWorkChain calculation for pure silver"
    
    return builder

def builder_o2_relax():
    """Set up a VASP relaxation calculation for O2 molecule"""
    builder = DFTWORKCHAIN.get_builder()
    
    # Load structure from file
    structure = StructureData(ase=read(O2_STRUCTURE_PATH))
    builder.structure = structure
    
    # Set VASP code
    builder.code = vasp_code
    
    # Set pseudopotentials
    builder.potential_family = Str(POTENTIAL_FAMILY)
    builder.potential_mapping = Dict(dict=POTENTIAL_MAPPING)
    
    # Define parameters for electronic structure and ionic relaxation
    parameters = {
        'incar': {
            # Electronic parameters
            'NCORE': NCORE,
            'ALGO': ALGO,
            'PREC': PREC,
            'ENCUT': ENCUT,
            'ISPIN': 2,      # Spin polarization (important for O2)
            'ISMEAR': 0,     # Gaussian smearing for molecule
            'SIGMA': 0.01,   # Small sigma for molecule
            'NELM': 60,
            'NELMIN': 6,
            'LREAL': 'Auto',
            'EDIFF': EDIFF,
            
            # Ionic relaxation parameters
            'ISIF': 2,       # Only relax atomic positions for molecule
            'IBRION': IBRION,
            'NSW': NSW,
            'EDIFFG': EDIFFG,
            
            # Important parameters for molecular calculation
            'LWAVE': True,
            'LCHARG': True,
        },
    }
    builder.parameters = Dict(dict=parameters)
    
    # Set kpoints (gamma-point only for molecule)
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([1, 1, 1])
    builder.kpoints = kpoints
    
    # Configure computational resources
    options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
    }
    if queue_name:
        options['queue_name'] = queue_name
        
    builder.options = Dict(dict=options)
    
    # Set parser settings
    settings = {
        'parser_settings': {
            'add_energies': True,
            'add_forces': True,
            'add_stress': False,  # No stress for molecule
        }
    }
    builder.settings = Dict(dict=settings)
    
    # Set restart settings
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    
    # Set meta info
    builder.metadata.label = "VASP O2 relaxation"
    builder.metadata.description = "VaspWorkChain calculation for O2 molecule"
    
    return builder

def builder_p_relax():
    """Set up a VASP relaxation calculation for pure phosphorus"""
    builder = DFTWORKCHAIN.get_builder()
    
    # Load structure from file
    structure = StructureData(ase=read(P_STRUCTURE_PATH))
    builder.structure = structure
    
    # Set VASP code
    builder.code = vasp_code
    
    # Set pseudopotentials
    builder.potential_family = Str(POTENTIAL_FAMILY)
    builder.potential_mapping = Dict(dict=POTENTIAL_MAPPING)
    
    # Define parameters for electronic structure and ionic relaxation
    parameters = {
        'incar': {
            # Electronic parameters
            'NCORE': NCORE,
            'ALGO': ALGO,
            'PREC': PREC,
            'ENCUT': ENCUT,
            'ISPIN': ISPIN,
            'ISMEAR': ISMEAR,
            'SIGMA': SIGMA,
            'NELM': 60,
            'NELMIN': 6,
            'LREAL': 'Auto',
            'EDIFF': EDIFF,
            
            # Ionic relaxation parameters
            'ISIF': ISIF,
            'IBRION': IBRION,
            'NSW': NSW,
            'EDIFFG': EDIFFG,
        },
    }
    builder.parameters = Dict(dict=parameters)
    
    # Set kpoints (appropriate mesh for phosphorus)
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([5, 5, 5])
    builder.kpoints = kpoints
    
    # Configure computational resources
    options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
    }
    if queue_name:
        options['queue_name'] = queue_name
        
    builder.options = Dict(dict=options)
    
    # Set parser settings
    settings = {
        'parser_settings': {
            'add_energies': True,
            'add_forces': True,
            'add_stress': True,
        }
    }
    builder.settings = Dict(dict=settings)
    
    # Set restart settings
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    
    # Set meta info
    builder.metadata.label = "VASP pure P relaxation"
    builder.metadata.description = "VaspWorkChain calculation for pure phosphorus"
    
    return builder

def get_reference_builders(elements):
    """
    Get reference builders for the specified elements
    
    Args:
        elements: List of elements in the oxide
        
    Returns:
        Dict mapping elements to their builders
    """
    ref_builders = {}
    
    for element in elements:
        if element == 'O':
            # Oxygen reference is typically an O2 molecule.
            # The key 'O2' is used internally by TEROS to identify the oxygen reference energy (halved for E_O).
            ref_builders['O2'] = builder_o2_relax()
        elif element == 'Ag':
            ref_builders['Ag'] = builder_ag_relax()
        elif element == 'P' and os.path.exists(P_STRUCTURE_PATH):
            ref_builders['P'] = builder_p_relax()
        else:
            print(f"Warning: No builder found for element {element}")
    
    return ref_builders

def run_workflow():
    """
    Configure and run the TEROS workflow using VASP
    
    Returns:
        WorkGraph: The submitted workgraph
    """
    # Get the builders
    bulk_builder = builder_bulk_relax()
    slab_builder = builder_slab_relax()
    
    bulk_structure = bulk_builder.structure
    
    elements = list(set(bulk_structure.get_ase().get_chemical_symbols()))
    
    # Get reference builders for all elements in the structure
    reference_builders = get_reference_builders(elements)
    
    # Determine which slab approach to use: manual or automatic
    manual_slabs = load_manual_slabs() if USE_MANUAL_SLABS else None
    
    # Build and submit the workflow using the original approach
    wg = create_teros_workgraph(
        dft_workchain=DFTWORKCHAIN,
        builder_bulk=bulk_builder,
        builder_slab=slab_builder,
        reference_builders=reference_builders,
        workgraph_name=WORKGRAPH_NAME,
        code=CODE,
        sampling=Int(SAMPLING),
        manual_slabs=manual_slabs,
        # Slab generation parameters
        miller_indices=List(list(MILLER_INDICES)),
        min_slab_thickness=Float(MIN_SLAB_THICKNESS),
        min_vacuum_thickness=Float(MIN_VACUUM_THICKNESS),
        lll_reduce=Bool(LLL_REDUCE),
        center_slab=Bool(CENTER_SLAB),
        symmetrize=Bool(SYMMETRIZE),
        primitive=Bool(PRIMITIVE),
        in_unit_planes=Bool(IN_UNIT_PLANES),
        max_normal_search=Int(MAX_NORMAL_SEARCH) if MAX_NORMAL_SEARCH is not None else None,
    )
    
    print(f"Submitting TEROS workflow using {CODE}...")
    print(f"Using {'manual slabs' if manual_slabs else 'automatic slab generation with default parameters'}")
    
    wg.submit(wait=False)  # wait=False submits the workflow and returns immediately (asynchronous).
                           # Set to True to wait for the workflow to complete (synchronous, useful for scripting).
    wg.to_html()           # Generates an HTML representation of the workgraph, useful for visualization and debugging.
                           # The file is typically saved as wg.html in the current directory.
    
    print("\nWorkflow submitted.")
    print(f"WorkGraph PK: {wg.pk}")
    print("You can check the status with: verdi process list -a")
    
    return wg

if __name__ == "__main__":
    wg = run_workflow()
