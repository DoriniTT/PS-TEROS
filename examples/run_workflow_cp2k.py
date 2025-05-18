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
WORKGRAPH_NAME = "aimd_ag3po4"

# Computational resources
CODE_LABEL = "CP2K-GPU-2022@bohr"  # GPU code
QUEUE_NAME = 'gpu_a100'

# Base directory for file paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# File paths
STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/bulk/Ag6O8P2_optimized.cif")
AG_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/Ag.cif")
P_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/P.cif")
O2_STRUCTURE_PATH = os.path.join(BASE_DIR, "structures/pure_elements/O2.cif")
BASIS_PATH = os.path.join(BASE_DIR, "potentials/BASIS_MOLOPT")
PSEUDO_PATH = os.path.join(BASE_DIR, "potentials/GTH_POTENTIALS")

# Additional retrieve list
ADDITIONAL_RETRIEVE_LIST = ["aiida-1.restart"]

# Optimization parameters
MAX_FORCE = 1.0e-3
MAX_STEPS = 200
OPTIMIZER = "BFGS"
OPTIMIZE_CELL = True
KEEP_SYMMETRY = True
KEEP_ANGLES = True
PRESSURE = 0.0

# Computational resource settings
NUM_MACHINES = 1
MAX_ITERATIONS = 3 
MAX_WALLTIME = 3 * 24 * 60 * 60  # 3 days

# Slab specific parameters
SLAB_PERIODICITY = "XY"
SLAB_DIPOLE_CORRECTION = False

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
    run_type = "CELL_OPT" if OPTIMIZE_CELL else "GEO_OPT"
    parameters = {
        "GLOBAL": {
            "RUN_TYPE": run_type,
            "PRINT_LEVEL": "LOW"
        },
    }
    if run_type == "GEO_OPT":
        parameters["MOTION"] = {
            "GEO_OPT": {
                "OPTIMIZER": OPTIMIZER,
                "MAX_FORCE": MAX_FORCE,
                "MAX_ITER": MAX_STEPS,
            },
        }
    else:
        parameters["MOTION"] = {
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
        }
    parameters["FORCE_EVAL"] = {
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
                {"_": "Ag", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
                {"_": "P", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q5", "POTENTIAL": "GTH-PBE-q5"},
                {"_": "O", "BASIS_SET": "ORB DZVP-MOLOPT-SR-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
            ]
        },
    }
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
    # Structure is NOT set here - it will be provided by the map_zone
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
            ref_builders['O2'] = builder_o2_relax()
        elif element == 'Ag':
            ref_builders['Ag'] = builder_ag_relax()
        elif element == 'P':
            ref_builders['P'] = builder_p_relax()
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
    wg.submit(wait=False)
    wg.to_html()
    print("\nWorkflow submitted.")
    print(f"WorkGraph PK: {wg.pk}")
    print("You can check the status with: verdi process list -a")
    return wg

if __name__ == "__main__":
    wg = run_workflow()
