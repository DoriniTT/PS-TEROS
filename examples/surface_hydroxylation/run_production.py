#!/home/thiagotd/envs/aiida/bin/python
"""
Production example for surface_hydroxylation module.

Uses realistic perovskite oxide surface with production VASP settings.

IMPORTANT:
  - Replace <PK_of_relaxed_slab> with actual PK from surface_thermodynamics output
  - Adjust compute resources (account, queue) to match your cluster
  - Verify VASP pseudopotentials are available
  - Review VASP settings for your specific system

Expected runtime: Several hours (depends on structure size and settings)
"""

from aiida import orm, load_profile
from teros.core.surface_hydroxylation import SurfaceHydroxylationWorkGraph

# Load AiiDA profile
load_profile('psteros')

# ============================================================================
# INPUT STRUCTURE
# ============================================================================
# Load production slab structure
# TODO: Replace with actual relaxed surface from surface_thermodynamics
# Example: structure = orm.load_node(12345)

# For testing, you can create a structure manually:
# from ase.build import fcc111
# slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
# # Add O atoms...
# structure = orm.StructureData(ase=slab)

print("ERROR: Please set the structure PK!")
print("Replace <PK_of_relaxed_slab> in this file with an actual PK.")
print("\nExample:")
print("  structure = orm.load_node(12345)")
exit(1)

structure = orm.load_node(<PK_of_relaxed_slab>)

print(f"Production slab: {len(structure.sites)} atoms")

# ============================================================================
# SURFACE GENERATION PARAMETERS
# ============================================================================
# Combined mode: Creates both oxygen-deficient and hydroxylated variants
surface_params = orm.Dict(dict={
    'mode': 'combine',  # 'vacancies', 'hydrogen', or 'combine'
    'species': 'O',     # Target species for modification
    'z_window': 0.5,    # Angstroms - surface identification window
    'which_surface': 'top',  # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,    # Angstroms - O-H bond distance
    'include_empty': False,  # Include unmodified structure
    'deduplicate_by_coverage': True,  # Remove symmetry-equivalent structures
    'coverage_bins': 10  # Number of coverage bins (controls # of structures)
})

print("\nSurface modification parameters:")
print(f"  Mode: {surface_params['mode']}")
print(f"  Coverage bins: {surface_params['coverage_bins']}")
print(f"  Expected structures: ~10-20 (depends on surface)")

# ============================================================================
# VASP CONFIGURATION
# ============================================================================
# Production VASP settings - adjust for your system!
builder_config = orm.Dict(dict={
    'metadata': {
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 16  # Adjust for your cluster
            },
            'max_wallclock_seconds': 3600 * 10,  # 10 hours
            'account': 'your_production_account',  # CHANGE THIS
            'queue_name': 'normal'  # CHANGE THIS
        }
    },
    'kpoints_distance': orm.Float(0.3),  # Converged k-points (test convergence!)
    'parameters': orm.Dict(dict={
        'EDIFF': 1e-6,      # Energy convergence
        'EDIFFG': -0.02,    # Force convergence (eV/Å)
        'ENCUT': 520,       # Plane wave cutoff (eV)
        'ISMEAR': 0,        # Gaussian smearing
        'SIGMA': 0.05,      # Smearing width
        'NSW': 200,         # Max ionic steps
        'IBRION': 2,        # Conjugate gradient
        'ISIF': 2,          # Relax atoms only (fixed cell)
        'LREAL': False,     # Accurate projectors
        'PREC': 'Accurate', # Precision mode
        'ALGO': 'Normal',   # Electronic minimization
        'NELM': 100         # Max electronic steps
    }),
    'relax': orm.Dict(dict={
        'positions': True,
        'shape': False,
        'volume': False
    })
})

print("\nVASP configuration:")
print(f"  Cores: {builder_config['metadata']['options']['resources']['num_mpiprocs_per_machine']}")
print(f"  Walltime: {builder_config['metadata']['options']['max_wallclock_seconds'] / 3600} hours")
print(f"  ENCUT: {builder_config['parameters']['ENCUT']} eV")
print(f"  EDIFFG: {builder_config['parameters']['EDIFFG']} eV/Å")

# ============================================================================
# PARALLELIZATION CONTROL
# ============================================================================
# Limit parallel jobs based on cluster capacity
# Set to number of VASP jobs you want running simultaneously
max_parallel = orm.Int(5)

print(f"\nParallelization: Max {max_parallel.value} VASP jobs at once")

# ============================================================================
# WORKFLOW SETUP AND SUBMISSION
# ============================================================================
print("\nCreating production workflow...")

# Create WorkGraph and add the hydroxylation task
from aiida_workgraph import WorkGraph
wg = WorkGraph(name='surface_hydroxylation_production')

wg.add_task(
    SurfaceHydroxylationWorkGraph,
    name='hydroxylation',
    structure=structure,
    surface_params=surface_params,
    builder_config=builder_config,
    max_parallel_jobs=max_parallel,
)

print("Submitting production workflow...")
result = wg.submit()

print("\n" + "="*70)
print(f"WORKFLOW SUBMITTED: PK = {result.pk}")
print("="*70)
print(f"\nThis will take several hours. Monitor with:")
print(f"  verdi process show {result.pk}")
print(f"  verdi process report {result.pk}")
print(f"\nCheck status periodically:")
print(f"  verdi process list -a  # List all processes")
print(f"\nOnce complete, analyze results:")
print(f"  verdi process show {result.pk}")
print(f"  # Access outputs: node = orm.load_node({result.pk})")
print(f"  # Statistics: node.outputs.statistics.get_dict()")
print(f"  # Results: node.outputs.successful_relaxations")
print("="*70 + "\n")
