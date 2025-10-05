"""
Runner script for TEROS DFT workflow using CP2K.

This script combines the builder configuration and workflow execution into a single file.
It sets up the CP2K calculation parameters and runs the create_teros_workgraph 
function for surface thermodynamics.
"""
import os
import sys
from ase.io import read
from aiida.orm import Int, Dict, Bool, SinglefileData, StructureData, load_code
from aiida.plugins import WorkflowFactory
from aiida.common.extendeddicts import AttributeDict
from aiida import load_profile

# Import the create_teros_workgraph function
from teros import create_teros_workgraph

# Load AiiDA profile
load_profile()

#########################################################
# GLOBAL CONFIGURATION - Edit parameters below as needed
#########################################################

DFTWORKCHAIN = Cp2kWorkChain = WorkflowFactory('cp2k.base')

# Workflow name
WORKGRAPH_NAME = "aimd_ag3po4" # Example name, can be customized for your specific project

# Computational resources
CODE_LABEL = "CP2K-GPU-2022@bohr"  # Adjust to your CP2K code label configured in AiiDA (e.g., cp2k@mycluster)
QUEUE_NAME = 'gpu_a100'           # Specify your queue name if running on a cluster, or set to None for default/local execution.

# Base directory for file paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# File paths
STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/bulk/Ag6O8P2_optimized.cif") # Path to the primary bulk material structure file
AG_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/Ag.cif")    # Path to pure Ag structure file (for reference energy)
P_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/P.cif")      # Path to pure P structure file (for reference energy)
O2_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/O2.cif")    # Path to O2 molecule structure file (for oxygen reference energy)
BASIS_PATH = os.path.join(BASE_DIR, "potentials/BASIS_MOLOPT") # Path to CP2K basis set file (e.g., BASIS_MOLOPT). Ensure this file exists.
PSEUDO_PATH = os.path.join(BASE_DIR, "potentials/GTH_POTENTIALS") # Path to CP2K pseudopotential file (e.g., GTH_POTENTIALS). Ensure this file exists.

# Additional retrieve list for CP2K calculations.
# This tells AiiDA which files to retrieve from the calculation directory in addition to standard output files.
# "aiida-1.restart" is useful for restarting CP2K calculations.
ADDITIONAL_RETRIEVE_LIST = ["aiida-1.restart"]

# Optimization parameters
MAX_FORCE = 1.0e-3          # Target maximum force on atoms for ionic relaxation (Hartree/Bohr).
MAX_STEPS = 200             # Maximum number of ionic steps in a geometry optimization.
OPTIMIZER = "BFGS"          # Type of optimizer for geometry relaxation (e.g., BFGS, CG, LBFGS).
OPTIMIZE_CELL = True        # If True, perform cell optimization for bulk calculations. Set to False for fixed-cell relaxations.
KEEP_SYMMETRY = True        # If True, maintain symmetry during cell optimization. CP2K's ability to strictly enforce this can vary.
KEEP_ANGLES = True          # If True, keep cell angles fixed during cell optimization (relevant for certain crystal symmetries).
PRESSURE = 0.0              # External pressure for cell optimization (GPa). CP2K expects pressure in [bar]. This will be handled by AiiDA plugin or needs manual conversion.

# Computational resource settings
NUM_MACHINES = 1            # Number of machines (nodes) to use for the calculation.
MAX_ITERATIONS = 3          # Max restart iterations for the AiiDA workchain (handles job resubmissions if they fail or run out of time).
MAX_WALLTIME = 3 * 24 * 60 * 60  # Maximum wallclock time in seconds (here: 3 days).

# Slab specific parameters
SLAB_PERIODICITY = "XY"     # Periodicity for slab calculations ("NONE", "X", "Y", "Z", "XY", "XZ", "YZ", "XYZ").
                            # "XY" means periodic in X and Y directions, non-periodic in Z (typical for surfaces).
SLAB_DIPOLE_CORRECTION = False # Whether to apply dipole correction for slabs. Useful for asymmetric slabs.

# Load CP2K code
cp2k_code = load_code(CODE_LABEL)

# Set resource parameters based on code/computer label
resource_kwargs = {}
if "cluster" in CODE_LABEL.lower():
    resource_kwargs["num_mpiprocs_per_machine"] = 24
    queue_name = None
else:
    resource_kwargs["num_cores_per_machine"] = 40
    queue_name = QUEUE_NAME

#########################################################
# BUILDER FUNCTIONS
#########################################################

def builder_bulk_relax():
    # Determine run type based on global OPTIMIZE_CELL setting
    run_type = "CELL_OPT" if OPTIMIZE_CELL else "GEO_OPT"
    
    # Basic CP2K input parameters structure. This dictionary is converted to CP2K input format.
    # Sections (keys like "GLOBAL", "MOTION", "FORCE_EVAL") correspond to CP2K input sections.
    # Keywords within sections match CP2K keywords.
    parameters = {
        "GLOBAL": {  # Global parameters for the CP2K run
            "RUN_TYPE": run_type,       # Type of calculation (e.g., GEO_OPT, CELL_OPT, MD, ENERGY_FORCE)
            "PRINT_LEVEL": "LOW"      # Verbosity of the CP2K output (LOW, MEDIUM, HIGH, DEBUG)
        },
    }
    
    # MOTION section: parameters for geometry optimization or cell optimization
    if run_type == "GEO_OPT": # Geometry optimization (fixed cell)
        parameters["MOTION"] = {
            "GEO_OPT": {
                "OPTIMIZER": OPTIMIZER,     # Algorithm for optimization (e.g., BFGS, CG)
                "MAX_FORCE": MAX_FORCE,     # Convergence criterion for max force on any atom
                "MAX_ITER": MAX_STEPS,      # Maximum number of optimization steps
            },
        }
    else: # Cell optimization (variable cell)
        parameters["MOTION"] = {
            "CELL_OPT": {
                "TYPE": "DIRECT_CELL_OPT", # Type of cell optimization strategy
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
                "EXTERNAL_PRESSURE": PRESSURE, # Target external pressure in [bar] (AiiDA plugin might handle unit conversion if input is GPa)
                "KEEP_SYMMETRY": KEEP_SYMMETRY, # Attempt to preserve crystal symmetry
                "KEEP_ANGLES": KEEP_ANGLES,     # Keep cell angles fixed (e.g., for orthorhombic to tetragonal transitions)
                "PRESSURE_TOLERANCE": 100.0,    # Tolerance for pressure convergence during cell optimization
            },
        }
        
    # FORCE_EVAL section: parameters for calculating forces and energies (the core DFT setup)
    parameters["FORCE_EVAL"] = {
        "METHOD": "QS",  # Using the Quickstep method (Gaussian and Plane Waves DFT)
        "STRESS_TENSOR": "ANALYTICAL", # How to compute stress tensor (ANALYTICAL, NUMERICAL)
        "DFT": { # DFT specific parameters
            "BASIS_SET_FILE_NAME": "BASIS_MOLOPT", # Filename for basis sets. Must match the name of the SinglefileData provided to `builder.cp2k.file['basis']`.
            "POTENTIAL_FILE_NAME": "GTH_POTENTIALS", # Filename for pseudopotentials. Must match `builder.cp2k.file['pseudo']`.
            "CHARGE": 0,        # Total charge of the system
            "MULTIPLICITY": 1,  # Spin multiplicity (2S+1). 1 for closed-shell, 2 for doublet, 3 for triplet (like O2).
            "XC": { # Exchange-Correlation functional settings
                "XC_FUNCTIONAL": {"_": "PBE"} # Defines the XC functional (e.g., PBE, LDA, BLYP). The "_" key is a CP2K convention.
            },
            "POISSON": { # Poisson solver settings for electrostatic interactions
                "PERIODIC": "XYZ",      # Periodicity of the Poisson solver (e.g., XYZ for bulk, XY for slabs, NONE for molecules)
                "PSOLVER": "PERIODIC",  # Type of Poisson solver (e.g., PERIODIC, MT for Martyna-Tuckerman for isolated systems)
            },
            "SCF": { # Self-Consistent Field (SCF) cycle parameters
                "SCF_GUESS": "ATOMIC", # Initial guess for wavefunction (e.g., ATOMIC, RESTART, RANDOM)
                "EPS_SCF": 1.0e-6,     # Convergence criterion for SCF energy difference between iterations
                "MAX_SCF": 80,         # Maximum number of SCF cycles
                "OT": {                # Settings for Orbital Transformation (OT) method (often more robust than DIIS for difficult convergence)
                    "PRECONDITIONER": "FULL_SINGLE_INVERSE", # Preconditioner for OT
                    "MINIMIZER": "DIIS" # Minimization algorithm for OT (can also be CG, SD etc.)
                },
            },
            "QS": { # Quickstep (QS) specific settings
                "EPS_DEFAULT": 1.0e-10, # Default accuracy for QS method (related to grid generation)
            }
        },
        "SUBSYS": { # Subsystem definition (atoms, basis sets, potentials, cell)
            "KIND": [ # Define parameters for each atomic kind (element).
                       # Ensure BASIS_SET and POTENTIAL names are correct and present in the provided files.
                       # The string after "q" usually denotes the number of valence electrons for that pseudo/basis.
                {"_": "Ag", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
                {"_": "P", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q5", "POTENTIAL": "GTH-PBE-q5"},
                {"_": "O", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
            ]
        },
    }
    
    # Get the builder for Cp2kWorkChain from aiida-cp2k plugin
    builder = Cp2kWorkChain.get_builder()
    builder.cp2k.code = cp2k_code
    builder.cp2k.structure = StructureData(ase=read(STRUCTURE_PATH))
    builder.cp2k.parameters = Dict(parameters)
    resources = {"num_machines": NUM_MACHINES}
    resources.update(resource_kwargs)
    builder.cp2k.metadata = {
        'options': {
            'resources': resources,
            'max_wallclock_seconds': MAX_WALLTIME,
            'queue_name': queue_name,
        }
    }
    builder.cp2k.file = {
        "basis": SinglefileData(file=BASIS_PATH),
        "pseudo": SinglefileData(file=PSEUDO_PATH),
    }
    settings = AttributeDict()
    settings.additional_retrieve_list = ADDITIONAL_RETRIEVE_LIST
    builder.cp2k.settings = Dict(dict=settings)
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    return builder

def builder_slab_relax():
    parameters = {
        "GLOBAL": {
            "RUN_TYPE": "GEO_OPT",
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "GEO_OPT": {
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "STRESS_TENSOR": "ANALYTICAL",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "MULTIPLICITY": 1,
                "SURFACE_DIPOLE_CORRECTION": SLAB_DIPOLE_CORRECTION,
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"}
                },
                "POISSON": {
                    "PERIODIC": SLAB_PERIODICITY,
                    "PSOLVER": "MT",
                },
                "SCF": {
                    "SCF_GUESS": "ATOMIC",
                    "EPS_SCF": 1.0e-6,
                    "MAX_SCF": 80,
                    "OT": {
                        "PRECONDITIONER": "FULL_SINGLE_INVERSE",
                        "MINIMIZER": "DIIS"
                    },
                },
                "QS": {
                    "EPS_DEFAULT": 1.0e-10,
                }
            },
            "SUBSYS": {
                "KIND": [
                    {"_": "Ag", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
                    {"_": "P", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q5", "POTENTIAL": "GTH-PBE-q5"},
                    {"_": "O", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
                ]
            },
        },
    }
    if SLAB_DIPOLE_CORRECTION and SLAB_PERIODICITY in ("XY", "XZ", "YZ"):
        non_periodic = 'Z' if SLAB_PERIODICITY == "XY" else ('Y' if SLAB_PERIODICITY == "XZ" else 'X')
        parameters["FORCE_EVAL"]["DFT"]["SURF_DIP_DIR"] = non_periodic
    builder = Cp2kWorkChain.get_builder()
    builder.cp2k.code = cp2k_code
    # Structure is NOT set here for the slab_relax builder.
    # The TEROS workflow will generate slab structures using Pymatgen 
    # and inject them into copies of this builder within the workgraph.
    builder.cp2k.parameters = Dict(parameters) # Pass the CP2K parameters dictionary
    resources = {"num_machines": NUM_MACHINES}
    resources.update(resource_kwargs)
    builder.cp2k.metadata = {
        'options': {
            'resources': resources,
            'max_wallclock_seconds': MAX_WALLTIME,
            'queue_name': queue_name,
        }
    }
    builder.cp2k.file = {
        "basis": SinglefileData(file=BASIS_PATH),
        "pseudo": SinglefileData(file=PSEUDO_PATH),
    }
    settings = AttributeDict()
    settings.additional_retrieve_list = ADDITIONAL_RETRIEVE_LIST
    builder.cp2k.settings = Dict(dict=settings)
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    return builder

def builder_ag_relax():
    run_type = "CELL_OPT"
    parameters = {
        "GLOBAL": {
            "RUN_TYPE": run_type,
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "CELL_OPT": {
                "TYPE": "DIRECT_CELL_OPT",
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
                "EXTERNAL_PRESSURE": PRESSURE,
                "KEEP_SYMMETRY": KEEP_SYMMETRY,
                "KEEP_ANGLES": KEEP_ANGLES,
                "PRESSURE_TOLERANCE": 100.0,
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "STRESS_TENSOR": "ANALYTICAL",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "MULTIPLICITY": 1,
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"}
                },
                "POISSON": {
                    "PERIODIC": "XYZ",
                    "PSOLVER": "PERIODIC",
                },
                "SCF": {
                    "SCF_GUESS": "ATOMIC",
                    "EPS_SCF": 1.0e-6,
                    "MAX_SCF": 80,
                    "SMEAR": {
                        "_": "ON",
                        "METHOD": "FERMI_DIRAC",
                        "ELECTRONIC_TEMPERATURE": 300,
                    },
                    "ADDED_MOS": 50,
                    "MIXING": {
                        "METHOD": "BROYDEN_MIXING",
                        "ALPHA": 0.4,
                        "NBROYDEN": 8
                    }
                },
                "QS": {
                    "EPS_DEFAULT": 1.0e-10,
                }
            },
            "SUBSYS": {
                "KIND": [
                    {"_": "Ag", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
                ]
            },
        },
    }
    builder = Cp2kWorkChain.get_builder()
    builder.cp2k.code = cp2k_code
    builder.cp2k.structure = StructureData(ase=read(AG_STRUCTURE_PATH))
    builder.cp2k.parameters = Dict(parameters)
    resources = {"num_machines": NUM_MACHINES}
    resources.update(resource_kwargs)
    builder.cp2k.metadata = {
        'options': {
            'resources': resources,
            'max_wallclock_seconds': MAX_WALLTIME,
            'queue_name': queue_name,
        }
    }
    builder.cp2k.file = {
        "basis": SinglefileData(file=BASIS_PATH),
        "pseudo": SinglefileData(file=PSEUDO_PATH),
    }
    settings = AttributeDict()
    settings.additional_retrieve_list = ADDITIONAL_RETRIEVE_LIST
    builder.cp2k.settings = Dict(dict=settings)
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    return builder

def builder_p_relax():
    run_type = "CELL_OPT"
    parameters = {
        "GLOBAL": {
            "RUN_TYPE": run_type,
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "CELL_OPT": {
                "TYPE": "DIRECT_CELL_OPT",
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
                "EXTERNAL_PRESSURE": PRESSURE,
                "KEEP_SYMMETRY": KEEP_SYMMETRY,
                "KEEP_ANGLES": KEEP_ANGLES,
                "PRESSURE_TOLERANCE": 100.0,
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "STRESS_TENSOR": "ANALYTICAL",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "MULTIPLICITY": 1,
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"}
                },
                "POISSON": {
                    "PERIODIC": "XYZ",
                    "PSOLVER": "PERIODIC",
                },
                "SCF": {
                    "SCF_GUESS": "ATOMIC",
                    "EPS_SCF": 1.0e-6,
                    "MAX_SCF": 80,
                    "OT": {
                        "PRECONDITIONER": "FULL_SINGLE_INVERSE",
                        "MINIMIZER": "DIIS"
                    },
                },
                "QS": {
                    "EPS_DEFAULT": 1.0e-10,
                }
            },
            "SUBSYS": {
                "KIND": [
                    {"_": "P", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q5", "POTENTIAL": "GTH-PBE-q5"},
                ]
            },
        },
    }
    builder = Cp2kWorkChain.get_builder()
    builder.cp2k.code = cp2k_code
    builder.cp2k.structure = StructureData(ase=read(P_STRUCTURE_PATH))
    builder.cp2k.parameters = Dict(parameters)
    resources = {"num_machines": NUM_MACHINES}
    resources.update(resource_kwargs)
    builder.cp2k.metadata = {
        'options': {
            'resources': resources,
            'max_wallclock_seconds': MAX_WALLTIME,
            'queue_name': queue_name,
        }
    }
    builder.cp2k.file = {
        "basis": SinglefileData(file=BASIS_PATH),
        "pseudo": SinglefileData(file=PSEUDO_PATH),
    }
    settings = AttributeDict()
    settings.additional_retrieve_list = ADDITIONAL_RETRIEVE_LIST
    builder.cp2k.settings = Dict(dict=settings)
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    return builder

def builder_o2_relax():
    run_type = "GEO_OPT"
    parameters = {
        "GLOBAL": {
            "RUN_TYPE": run_type,
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "GEO_OPT": {
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "MULTIPLICITY": 3,
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"}
                },
                "POISSON": {
                    "PERIODIC": "NONE",
                    "PSOLVER": "MT",
                },
                "SCF": {
                    "SCF_GUESS": "ATOMIC",
                    "EPS_SCF": 1.0e-6,
                    "MAX_SCF": 100,
                    "OT": {
                        "PRECONDITIONER": "FULL_SINGLE_INVERSE",
                        "MINIMIZER": "DIIS"
                    },
                },
                "QS": {
                    "EPS_DEFAULT": 1.0e-10,
                }
            },
            "SUBSYS": {
                "KIND": [
                    {"_": "O", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
                ],
                "CELL": {
                    "ABC": "15.0 15.0 15.0",
                    "PERIODIC": "NONE"
                }
            },
        },
    }
    builder = Cp2kWorkChain.get_builder()
    builder.cp2k.code = cp2k_code
    builder.cp2k.structure = StructureData(ase=read(O2_STRUCTURE_PATH))
    builder.cp2k.parameters = Dict(parameters)
    resources = {"num_machines": NUM_MACHINES}
    resources.update(resource_kwargs)
    builder.cp2k.metadata = {
        'options': {
            'resources': resources,
            'max_wallclock_seconds': MAX_WALLTIME,
            'queue_name': queue_name,
        }
    }
    builder.cp2k.file = {
        "basis": SinglefileData(file=BASIS_PATH),
        "pseudo": SinglefileData(file=PSEUDO_PATH),
    }
    settings = AttributeDict()
    settings.additional_retrieve_list = ADDITIONAL_RETRIEVE_LIST
    builder.cp2k.settings = Dict(dict=settings)
    builder.max_iterations = Int(MAX_ITERATIONS)
    builder.clean_workdir = Bool(False)
    return builder

#########################################################
# CP2K WORKGRAPH CONFIGURATION
#########################################################

CODE = 'CP2K'
CUSTOM_WORKGRAPH_NAME = f"{WORKGRAPH_NAME}_teros"

def get_reference_builders(elements):
    ref_builders = {}
    for element in elements:
        if element == 'O':
            # Oxygen reference is typically an O2 molecule (triplet state).
            # The key 'O2' is used internally by TEROS to identify the oxygen reference energy (which is then halved for E_O).
            ref_builders['O2'] = builder_o2_relax()
        elif element == 'Ag':
            ref_builders['Ag'] = builder_ag_relax() # Reference for elemental Ag
        elif element == 'P': # Ensure P_STRUCTURE_PATH points to a valid structure if P is in your material
            ref_builders['P'] = builder_p_relax()   # Reference for elemental P
        else:
            print(f"Warning: No builder found for element {element}")
    return ref_builders

def run_workflow():
    bulk_builder = builder_bulk_relax()
    slab_builder = builder_slab_relax()
    bulk_structure = bulk_builder.cp2k.structure
    elements = list(set(bulk_structure.get_ase().get_chemical_symbols()))
    reference_builders = get_reference_builders(elements)
    wg = create_teros_workgraph(
        dft_workchain=DFTWORKCHAIN,
        builder_bulk=bulk_builder,
        builder_slab=slab_builder,
        reference_builders=reference_builders,
        workgraph_name=CUSTOM_WORKGRAPH_NAME,
        code=CODE
    )
    print(f"Submitting TEROS workflow using {CODE}...")
    wg.submit(wait=False)  # wait=False submits the workflow and returns immediately (asynchronous execution).
                           # Set to True to make the script wait for the workflow to complete (synchronous execution).
    wg.to_html()           # Generates an HTML representation of the workgraph's structure (e.g., wg.html).
                           # This is useful for visualizing the workflow components and their connections.
    print("\nWorkflow submitted.")
    print(f"WorkGraph PK: {wg.pk}")
    print("You can check the status with: verdi process list -a")
    return wg

if __name__ == "__main__":
    wg = run_workflow()
