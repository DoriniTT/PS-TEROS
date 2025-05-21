"""
Runner script for TEROS DFT workflow using Quantum ESPRESSO.

This script combines the builder configuration and workflow execution into a single file.
It sets up the Quantum ESPRESSO calculation parameters and runs the create_teros_workgraph 
function for surface thermodynamics.
"""
import os
import sys
from ase.io import read
from aiida.orm import Dict, Bool, Str, Int, Float, StructureData, load_code
from aiida.plugins import WorkflowFactory
from aiida import load_profile

# Import the create_teros_workgraph function
from teros import create_teros_workgraph

# Import QE-specific types
from aiida_quantumespresso.common.types import SpinType, RelaxType

# Load AiiDA profile
load_profile()

#########################################################
# GLOBAL CONFIGURATION - Edit parameters below as needed
#########################################################

# Import the workchain factory
PwRelaxWorkChain = DFTWORKCHAIN = WorkflowFactory('quantumespresso.pw.relax')

# Workflow name
WORKGRAPH_NAME = "teros_ag3po4_qe"

# Define DFT code to use
CODE = 'QUANTUM_ESPRESSO'

# Computational resources
CODE_LABEL = "QE-PW-6.6@bohr"  # Adjust to your Quantum ESPRESSO `pw.x` code label configured in AiiDA (e.g., qe-pw@mycluster)

# File paths
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/bulk/Ag6O8P2_optimized.cif")
AG_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/Ag.cif")
P_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/P.cif")
O2_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/O2.cif")

# Pseudopotential family name
PSEUDO_FAMILY = "SSSP/1.3/PBE/efficiency"  # Ensure this pseudo family is installed in AiiDA.

# Protocol level (accuracy vs speed) - 'fast', 'moderate', 'precise'
# 'fast': suitable for quick tests, uses lower cutoffs and k-point densities.
# 'moderate': a balance between speed and accuracy.
# 'precise': aims for high accuracy, uses tighter convergence criteria and denser settings.
PROTOCOL = "fast"  

# Key calculation parameters
MAX_FORCE = 1.0e-3          # Target maximum force on atoms for ionic relaxation (Ry/bohr).
MAX_FORCE = 1.0e-3          # Target maximum force on atoms for ionic relaxation (Ry/bohr).
MAX_STEPS = 200             # Maximum number of ionic steps in a relaxation.
SYSTEM_2D = True            # If True, applies settings suitable for 2D systems (slabs) like 'assume_isolated'.

# Computational resource settings
NUM_MACHINES = 1
NUM_CORES_PER_MACHINE = 40
QUEUE_NAME = 'par40'
MAX_WALLTIME = 3 * 24 * 60 * 60  # 3 days

# Load QE pw.x code
qe_code = load_code(CODE_LABEL)

# Set resource parameters based on code/computer label
resource_kwargs = {}
if "cluster" in CODE_LABEL.lower():
    resource_kwargs["num_mpiprocs_per_machine"] = 24
    queue_name = None
else:  # Default or "bohr"
    resource_kwargs["num_cores_per_machine"] = NUM_CORES_PER_MACHINE
    queue_name = QUEUE_NAME

#########################################################
# BUILDER FUNCTIONS
#########################################################

def builder_bulk_relax():
    """
    Set up a Quantum ESPRESSO variable-cell relaxation calculation using the PwRelaxWorkChain
    with protocol-based configuration.
    
    Returns:
        builder: Configured QE relaxation workflow builder
    """
    # Load structure from file
    structure = StructureData(ase=read(STRUCTURE_PATH))
    
    # Get builder with protocol
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=qe_code,
        structure=structure,
        protocol=PROTOCOL,
        relaxation_type=RelaxType.CELL,  # Variable-cell relaxation
        spin_type=SpinType.COLLINEAR,    # Enable spin polarization
        pseudo_family=PSEUDO_FAMILY
    )
    
    # Adjust a few key parameters
    parameters = builder.base.pw.parameters.get_dict()
    
    # Modify convergence parameters for forces
    parameters['CONTROL']['forc_conv_thr'] = MAX_FORCE
    parameters['CONTROL']['nstep'] = MAX_STEPS
    
    # Add scf_must_converge under ELECTRONS
    # This allows the calculation to proceed even if SCF doesn't converge in one step,
    # which can be useful for difficult systems or initial relaxations.
    parameters.setdefault('ELECTRONS', {})
    parameters['ELECTRONS']['scf_must_converge'] = False
    
    # Update parameters
    builder.base.pw.parameters = Dict(parameters)
    
    # Configure computational resources
    builder.base.pw.metadata.options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
        'queue_name': queue_name,
    }
    
    # Set meta info
    builder.metadata.label = "QE bulk structure relaxation"
    builder.metadata.description = "PwRelax calculation for bulk optimization of Ag3PO4"
    builder.clean_workdir = Bool(False)
    
    return builder

def builder_slab_relax():
    """
    Set up a Quantum ESPRESSO slab relaxation with protocol-based configuration.
    The structure is not set here as it will be provided by the map_zone in the workflow.
    
    Returns:
        builder: Configured QE workflow builder for slabs
    """
    # Create a dummy structure to initialize the protocol builder.
    # The actual slab structures will be generated by the TEROS workflow
    # and injected into copies of this builder.
    dummy_structure = StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 20]])
    dummy_structure.append_atom(position=(0, 0, 0), symbols='Ag') # Example atoms
    dummy_structure.append_atom(position=(0, 0, 2), symbols='O')
    dummy_structure.append_atom(position=(2, 2, 0), symbols='P')
    
    # Get builder with protocol
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=qe_code,
        structure=dummy_structure,
        protocol=PROTOCOL, 
        relaxation_type=RelaxType.POSITIONS,  # Only relax atomic positions
        spin_type=SpinType.COLLINEAR,         # Enable spin polarization
        pseudo_family=PSEUDO_FAMILY
    )
    
    # Remove the dummy structure as it will be replaced by the workflow
    delattr(builder, 'structure')
    
    # Adjust parameters for slab calculations
    parameters = builder.base.pw.parameters.get_dict()
    
    # Modify convergence parameters
    parameters['CONTROL']['forc_conv_thr'] = MAX_FORCE
    parameters['CONTROL']['nstep'] = MAX_STEPS
    
    # Configure slab-specific settings
    if SYSTEM_2D:
        # For 2D systems (slabs), 'assume_isolated' helps in treating the system as non-periodic in z.
        parameters['SYSTEM']['assume_isolated'] = 'esm'  # Use Effective Screening Medium method.
        parameters['SYSTEM']['esm_bc'] = 'bc3'  # Boundary condition: vacuum-slab-vacuum.
                                               # Other options: 'pbc' (periodic), 'bc1' (vacuum-slab-solid), etc.
    
    # Add scf_must_converge under ELECTRONS
    parameters.setdefault('ELECTRONS', {})
    parameters['ELECTRONS']['scf_must_converge'] = False
    
    # Update parameters
    builder.base.pw.parameters = Dict(parameters)
    
    # Configure computational resources
    builder.base.pw.metadata.options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'withmpi': True,
        'queue_name': queue_name,
    }
    
    # Set meta info
    builder.metadata.label = "QE slab structure relaxation"
    builder.metadata.description = "PwRelax calculation for slab with fixed cell"
    builder.clean_workdir = Bool(False)
    return builder

def builder_ag_relax():
    """
    Set up a Quantum ESPRESSO relaxation for pure silver.
    
    Returns:
        builder: Configured QE relaxation workflow builder for pure Ag
    """
    # Load structure from file
    structure = StructureData(ase=read(AG_STRUCTURE_PATH))
    
    # Get builder with protocol
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=qe_code,
        structure=structure,
        protocol=PROTOCOL,
        relaxation_type=RelaxType.CELL,  # Variable-cell relaxation
        spin_type=SpinType.NONE,         # Non-spin-polarized for metal
        pseudo_family=PSEUDO_FAMILY
    )
    
    # Adjust parameters for metal calculations
    parameters = builder.base.pw.parameters.get_dict()
    
    # Modify convergence parameters for forces
    parameters['CONTROL']['forc_conv_thr'] = MAX_FORCE
    parameters['CONTROL']['nstep'] = MAX_STEPS
    
    # Metal-specific parameters
    parameters['SYSTEM']['occupations'] = 'smearing'
    parameters['SYSTEM']['smearing'] = 'cold'  # Cold smearing for metals
    parameters['SYSTEM']['degauss'] = 0.01     # Degauss parameter in Ry
    
    # Add scf_must_converge under ELECTRONS
    parameters.setdefault('ELECTRONS', {})
    parameters['ELECTRONS']['scf_must_converge'] = False
    
    # Update parameters
    builder.base.pw.parameters = Dict(parameters)
    
    # Configure computational resources
    builder.base.pw.metadata.options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'queue_name': queue_name,
    }
    
    # Set meta info
    builder.metadata.label = "QE pure Ag relaxation"
    builder.metadata.description = "PwRelax calculation for pure silver"
    builder.clean_workdir = Bool(False)
    
    return builder

def builder_p_relax():
    """
    Set up a Quantum ESPRESSO relaxation for pure phosphorus.
    
    Returns:
        builder: Configured QE relaxation workflow builder for pure P
    """
    # Load structure from file
    structure = StructureData(ase=read(P_STRUCTURE_PATH))
    
    # Get builder with protocol
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=qe_code,
        structure=structure,
        protocol=PROTOCOL,
        relaxation_type=RelaxType.CELL,  # Variable-cell relaxation
        spin_type=SpinType.COLLINEAR,    # Enable spin polarization
        pseudo_family=PSEUDO_FAMILY
    )
    
    # Adjust parameters
    parameters = builder.base.pw.parameters.get_dict()
    
    # Modify convergence parameters for forces
    parameters['CONTROL']['forc_conv_thr'] = MAX_FORCE
    parameters['CONTROL']['nstep'] = MAX_STEPS
    
    # Add scf_must_converge under ELECTRONS
    parameters.setdefault('ELECTRONS', {})
    parameters['ELECTRONS']['scf_must_converge'] = False
    
    # Update parameters
    builder.base.pw.parameters = Dict(parameters)
    
    # Configure computational resources
    builder.base.pw.metadata.options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'queue_name': queue_name,
    }
    
    # Set meta info
    builder.metadata.label = "QE pure P relaxation"
    builder.metadata.description = "PwRelax calculation for pure phosphorus"
    builder.clean_workdir = Bool(False)
    
    return builder

def builder_o2_relax():
    """
    Set up a Quantum ESPRESSO relaxation for O2 molecule.
    
    Returns:
        builder: Configured QE relaxation workflow builder for O2
    """
    # Load structure from file
    structure = StructureData(ase=read(O2_STRUCTURE_PATH))
    
    # Create a builder with protocol
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=qe_code,
        structure=structure,
        protocol=PROTOCOL,
        relaxation_type=RelaxType.POSITIONS,  # Only relax positions
        spin_type=SpinType.COLLINEAR,         # Enable spin polarization (needed for O2)
        pseudo_family=PSEUDO_FAMILY
    )
    
    # Adjust parameters for molecular calculation
    parameters = builder.base.pw.parameters.get_dict()
    
    # Modify convergence parameters
    parameters['CONTROL']['forc_conv_thr'] = MAX_FORCE
    parameters['CONTROL']['nstep'] = MAX_STEPS
    
    # Isolated molecule settings
    parameters['SYSTEM']['assume_isolated'] = 'martyna-tuckerman'  # Isolated molecule
    parameters['SYSTEM']['nosym'] = True  # No symmetry for molecule
    parameters['SYSTEM']['occupations'] = 'smearing'  # Smearing for better convergence
    parameters['SYSTEM']['smearing'] = 'gaussian'
    parameters['SYSTEM']['degauss'] = 0.005
    
    # Set initial magnetization for O (atom index 1, assuming O is the first type if sorted alphabetically or by pseudo order)
    # This is often necessary for O2 to converge to a triplet state.
    # Ensure the atom order matches what QE expects or use element-specific magnetization.
    parameters['SYSTEM']['starting_magnetization(1)'] = 0.5  # Example for oxygen
    parameters['SYSTEM']['nspin'] = 2  # Force spin-polarized calculation
    
    # Add scf_must_converge under ELECTRONS
    parameters.setdefault('ELECTRONS', {})
    parameters['ELECTRONS']['scf_must_converge'] = False
    
    # Update parameters
    builder.base.pw.parameters = Dict(parameters)
    
    # Configure computational resources
    builder.base.pw.metadata.options = {
        'resources': {
            'num_machines': NUM_MACHINES,
            **resource_kwargs,
        },
        'max_wallclock_seconds': MAX_WALLTIME,
        'queue_name': queue_name,
    }
    
    # Set meta info
    builder.metadata.label = "QE O2 molecule relaxation"
    builder.metadata.description = "PwRelax calculation for O2 molecule"
    builder.clean_workdir = Bool(False)
    
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
    Configure and run the TEROS workflow using Quantum ESPRESSO
    
    Returns:
        WorkGraph: The submitted workgraph
    """
    # Get the builders
    bulk_builder = builder_bulk_relax()
    slab_builder = builder_slab_relax()
    
    # Extract bulk structure for element analysis
    bulk_structure = bulk_builder.structure
    elements = list(set(bulk_structure.get_ase().get_chemical_symbols()))
    
    # Get reference builders for all elements in the structure
    reference_builders = get_reference_builders(elements)
    
    # Build and submit the workflow using the functional approach
    wg = create_teros_workgraph(
        dft_workchain=DFTWORKCHAIN,
        builder_bulk=bulk_builder,
        builder_slab=slab_builder,
        reference_builders=reference_builders,
        workgraph_name=WORKGRAPH_NAME,
        code=CODE
    )
    
    print(f"Submitting TEROS workflow using {CODE}...")
    wg.submit(wait=False)  # wait=False submits the workflow and returns immediately (asynchronous).
                           # Set to True to wait for the workflow to complete (synchronous).
    
    # Generate HTML visualization
    wg.to_html()           # Generates an HTML representation of the workgraph (e.g., wg.html).
    
    print("\nWorkflow submitted.")
    print(f"WorkGraph PK: {wg.pk}")
    print("You can check the status with: verdi process list -a")
    
    return wg

if __name__ == "__main__":
    wg = run_workflow()
