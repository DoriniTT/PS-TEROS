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
load_profile('presto')

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
surface_params = {
    'mode': 'combine',  # 'vacancies', 'hydrogen', or 'combine'
    'species': 'O',     # Target species for modification
    'z_window': 0.5,    # Angstroms - surface identification window
    'which_surface': 'top',  # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,    # Angstroms - O-H bond distance
    'include_empty': False,  # Include unmodified structure
    'deduplicate_by_coverage': True,  # Remove symmetry-equivalent structures
    'coverage_bins': 10  # Number of coverage bins (controls # of structures)
}

print("\nSurface modification parameters:")
print(f"  Mode: {surface_params['mode']}")
print(f"  Coverage bins: {surface_params['coverage_bins']}")
print(f"  Expected structures: ~10-20 (depends on surface)")

# ============================================================================
# VASP CONFIGURATION
# ============================================================================
# Production VASP settings - adjust for your system!
# Complete configuration for vasp.v2.relax workflow plugin

# IMPORTANT: Adjust these settings for your cluster and system
code_label = 'your_vasp_code@your_computer'  # CHANGE THIS
try:
    code = orm.load_code(code_label)
except:
    print(f"ERROR: Could not load VASP code '{code_label}'")
    print("Available codes:")
    from aiida.orm import QueryBuilder
    qb = QueryBuilder()
    qb.append(orm.Code, project=['label'])
    for (label,) in qb.iterall():
        print(f"  - {label}")
    exit(1)

# VASP configuration (split into separate components for WorkGraph serialization)
vasp_config = {
    # Relaxation settings
    'relax': {
        'perform': True,
        'positions': True,
        'shape': False,
        'volume': False,
        'force_cutoff': 0.02,  # eV/Angstrom - production convergence
        'steps': 200,          # Max ionic steps
        'algo': 'cg',          # Conjugate gradient
    },

    # Base VASP settings (production quality)
    'base': {
        'PREC': 'Accurate',
        'ENCUT': 520,          # Converged cutoff (test for your system!)
        'EDIFF': 1e-6,         # Energy convergence
        'ISMEAR': 0,           # Gaussian smearing
        'SIGMA': 0.05,         # Smearing width
        'ALGO': 'Normal',
        'LREAL': False,        # Accurate projectors (slower but more accurate)
        'NELM': 100,           # Max electronic steps
        'LWAVE': False,        # Don't write WAVECAR (saves space)
        'LCHARG': False,       # Don't write CHGCAR (saves space)
    },

    # K-points (IMPORTANT: test convergence!)
    'kpoints_spacing': 0.3,  # Angstrom^-1 (typical production value)

    # Pseudopotentials
    'potential_family': 'PBE',  # Or 'PBE.54', etc. - depends on your setup
    'potential_mapping': {},    # Optional: element-specific potentials

    # Cleanup
    'clean_workdir': False,  # Keep files for analysis (set True to save space)
}

# Scheduler options (separate from VASP config for serialization)
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 16,  # Adjust for your cluster
    },
    'queue_name': 'normal',              # CHANGE THIS
    'max_wallclock_seconds': 3600 * 10,  # 10 hours
    'account': 'your_account',           # CHANGE THIS if required
}

print("\nVASP configuration:")
print(f"  Code: {code.label}")
print(f"  Cores: {options['resources']['num_mpiprocs_per_machine']}")
print(f"  Walltime: {options['max_wallclock_seconds'] / 3600} hours")
print(f"  ENCUT: {vasp_config['base']['ENCUT']} eV")
print(f"  Force cutoff: {vasp_config['relax']['force_cutoff']} eV/Å")
print(f"  K-points spacing: {vasp_config['kpoints_spacing']} Å⁻¹")

# ============================================================================
# PARALLELIZATION CONTROL
# ============================================================================
# BATCH APPROACH: Process only the first N structures
# Example: If 15 structures generated, max_parallel=5 processes structures 0-4
# After completion, increase to 10 to process structures 5-9, etc.
max_parallel = 5

print(f"\nParallelization: Processing first {max_parallel} structures per run")

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
    code=code,
    vasp_config=vasp_config,
    options=options,
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
