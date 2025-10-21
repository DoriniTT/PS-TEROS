#!/home/thiagotd/envs/aiida/bin/python
"""
Production example for surface_hydroxylation module.

Uses realistic perovskite oxide surface with production VASP settings.

IMPORTANT:
  - Set the structure PK from a relaxed surface (from surface_thermodynamics or manual relaxation)
  - Adjust compute resources (account, queue, cores) to match your cluster
  - Verify VASP pseudopotentials are available for all elements
  - Review and test VASP settings for your specific system (especially ENCUT and k-points)
  - Test with small max_parallel first, then scale up

Expected runtime: Several hours to days (depends on structure size, # of variants, and settings)
"""

from aiida import orm, load_profile
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

# Load AiiDA profile
load_profile('presto')

# ============================================================================
# INPUT STRUCTURE
# ============================================================================
# Load production slab structure from a previous relaxation
#
# Method 1: Load from surface_thermodynamics output
#   structure = orm.load_node(12345)  # Replace with actual PK
#
# Method 2: Query for recent relaxed structures
#   from aiida.orm import QueryBuilder, StructureData
#   qb = QueryBuilder()
#   qb.append(StructureData, filters={'label': {'like': 'relaxed%'}})
#   qb.order_by({StructureData: {'ctime': 'desc'}})
#   structure = qb.first()[0]
#
# For this script, set the PK here:
STRUCTURE_PK = None  # ← SET THIS to your relaxed surface PK

if STRUCTURE_PK is None:
    print("\n" + "="*70)
    print("ERROR: Please set STRUCTURE_PK in this file!")
    print("="*70)
    print("\nTo find available structures:")
    print("  verdi data core.structure list")
    print("\nOr query for relaxed structures:")
    print("  from aiida.orm import QueryBuilder, StructureData")
    print("  qb = QueryBuilder()")
    print("  qb.append(StructureData)")
    print("  for (pk, label) in qb.iterall():")
    print("      print(f'PK: {pk}, Label: {label}')")
    print("\nThen set: STRUCTURE_PK = <your_pk>")
    print("="*70 + "\n")
    exit(1)

try:
    structure = orm.load_node(STRUCTURE_PK)
except Exception as e:
    print(f"\nERROR: Could not load structure PK {STRUCTURE_PK}")
    print(f"Error: {e}\n")
    exit(1)

print(f"\n{'='*70}")
print(f"Production slab loaded: {len(structure.sites)} atoms")
print(f"Structure PK: {STRUCTURE_PK}")
print(f"{'='*70}")

# ============================================================================
# SURFACE GENERATION PARAMETERS
# ============================================================================
# Combined mode: Creates both oxygen-deficient and hydroxylated variants
#
# Mode options:
#   'vacancies' - Only oxygen vacancy structures
#   'hydrogen'  - Only hydroxylated structures (adds H to surface O)
#   'combine'   - Both vacancies and hydroxylation (comprehensive)
#
# Coverage bins: Controls structure generation
#   - Higher values = more structures (finer coverage sampling)
#   - Typical values: 5-10 for production
#   - With deduplicate_by_coverage=True, duplicates are removed
#
surface_params = {
    'mode': 'combine',  # 'vacancies', 'hydrogen', or 'combine'
    'species': 'O',     # Target species for modification
    'z_window': 0.5,    # Angstroms - surface identification window
    'which_surface': 'top',  # 'top', 'bottom', or 'both'
    'oh_dist': 0.98,    # Angstroms - O-H bond distance (typical for hydroxyl)
    'include_empty': False,  # Include unmodified structure
    'deduplicate_by_coverage': True,  # Remove symmetry-equivalent structures
    'coverage_bins': 10  # Number of coverage bins (controls # of structures)
}

print("\nSurface modification parameters:")
print(f"  Mode: {surface_params['mode']}")
print(f"  Target species: {surface_params['species']}")
print(f"  Coverage bins: {surface_params['coverage_bins']}")
print(f"  Deduplication: {surface_params['deduplicate_by_coverage']}")
print(f"  Expected structures: ~10-20 (depends on surface symmetry)")

# ============================================================================
# VASP CONFIGURATION
# ============================================================================
# Production VASP settings - adjust for your system!
#
# IMPORTANT: This module uses vasp.v2.vasp (NOT vasp.v2.relax)
# The workflow handles relaxation by setting INCAR parameters (ISIF, NSW, IBRION, EDIFFG)
#
# CRITICAL: Test convergence for your system!
#   - ENCUT: Run convergence tests (typically 1.3x max ENMAX from POTCAR)
#   - k-points: Test kpoints_spacing convergence (0.25-0.35 typical)
#   - Force cutoff: 0.02-0.05 eV/Å typical for surfaces
#
# IMPORTANT: Adjust these settings for your cluster and system
code_label = 'your_vasp_code@your_computer'  # ← CHANGE THIS

# Validate code exists (builder function will load it)
try:
    orm.load_code(code_label)
except Exception as e:
    print(f"\n{'='*70}")
    print(f"ERROR: Could not load VASP code '{code_label}'")
    print(f"{'='*70}")
    print(f"Error: {e}\n")
    print("Available codes:")
    from aiida.orm import QueryBuilder
    qb = QueryBuilder()
    qb.append(orm.Code, project=['label', 'description'])
    found_codes = False
    for (label, desc) in qb.iterall():
        print(f"  - {label}")
        if desc:
            print(f"    {desc}")
        found_codes = True
    if not found_codes:
        print("  (No codes found - you need to set up VASP code first)")
        print("\n  Setup VASP code with:")
        print("    verdi code create core.code.installed")
    print(f"\n{'='*70}\n")
    exit(1)

# VASP configuration
# Direct INCAR parameters for vasp.v2.vasp workflow
#
# IMPORTANT: Test convergence for your system!
#   - ENCUT: Typically 1.3x max ENMAX from POTCAR
#   - k-points: Test kpoints_spacing (0.25-0.35 typical for surfaces)
#   - EDIFFG: Force convergence (0.02-0.05 eV/Å typical)
#
vasp_config = {
    # INCAR parameters (production quality - adjust for your system!)
    'parameters': {
        # Electronic structure
        'PREC': 'Accurate',     # Precision level
        'ENCUT': 520,           # Plane-wave cutoff (eV) - TEST CONVERGENCE!
        'EDIFF': 1e-6,          # Electronic convergence (eV)
        'ISMEAR': 0,            # Gaussian smearing (good for surfaces)
        'SIGMA': 0.05,          # Smearing width (eV)
        'ALGO': 'Normal',       # Electronic minimization algorithm
        'LREAL': False,         # Real-space projection (False = more accurate)
        'NELM': 100,            # Max electronic steps
        'LWAVE': False,         # Don't write WAVECAR (saves disk space)
        'LCHARG': False,        # Don't write CHGCAR (saves disk space)

        # Relaxation parameters
        'ISIF': 2,              # Relax positions only (2=pos, 3=pos+cell, 7=all)
        'NSW': 200,             # Max ionic steps
        'IBRION': 2,            # Ionic relaxation: 2=CG, 1=RMM-DIIS
        'EDIFFG': -0.02,        # Force convergence (eV/Å, negative = force criterion)
    },

    # K-points sampling (IMPORTANT: test convergence!)
    # Typical values: 0.25-0.35 Å⁻¹ for surfaces
    # Lower values = denser k-point mesh = more accurate (but slower)
    'kpoints_spacing': 0.3,     # Angstrom^-1

    # Pseudopotentials
    # Check available families with: verdi data core.upf listfamilies
    'potential_family': 'PBE',  # e.g., 'PBE', 'PBE.54', 'PBE_54'
    'potential_mapping': {},    # Optional: {'Element': 'PotentialName'}
                                # Example: {'Ti': 'Ti_pv', 'O': 'O'}

    # Cleanup
    'clean_workdir': False,     # Keep remote files for debugging
                                # Set True in production to save disk space
}

# Scheduler options (separate from VASP config for proper serialization)
# IMPORTANT: Adjust for your cluster configuration
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 16,  # ← ADJUST: cores per calculation
                                          # Typical: 16-32 for medium systems
                                          # Balance: more cores = faster but less throughput
    },
    'queue_name': 'normal',               # ← CHANGE: your cluster queue
    'max_wallclock_seconds': 3600 * 10,   # 10 hours (36000 seconds)
                                          # Adjust based on system size and settings
    'account': 'your_account',            # ← CHANGE: if your cluster requires account
    # Optional settings:
    # 'prepend_text': 'module load vasp/6.4.2',  # Commands before VASP
    # 'custom_scheduler_commands': '#SBATCH --constraint=haswell',
}

print("\nVASP configuration:")
print(f"  Code: {code_label}")
print(f"  Workflow: vasp.v2.vasp")
print(f"  Cores per job: {options['resources']['num_mpiprocs_per_machine']}")
print(f"  Walltime: {options['max_wallclock_seconds'] / 3600:.1f} hours")
print(f"  Queue: {options['queue_name']}")
print(f"\n  INCAR parameters:")
print(f"    ENCUT: {vasp_config['parameters']['ENCUT']} eV")
print(f"    EDIFF: {vasp_config['parameters']['EDIFF']} eV")
print(f"    EDIFFG: {vasp_config['parameters']['EDIFFG']} eV/Å")
print(f"    NSW: {vasp_config['parameters']['NSW']}")
print(f"    ISIF: {vasp_config['parameters']['ISIF']} (positions only)")
print(f"    IBRION: {vasp_config['parameters']['IBRION']} (CG)")
print(f"\n  K-points:")
print(f"    Spacing: {vasp_config['kpoints_spacing']} Å⁻¹")
print(f"\n  Pseudopotentials:")
print(f"    Family: {vasp_config['potential_family']}")

# ============================================================================
# PARALLELIZATION CONTROL
# ============================================================================
# SIMPLE BATCH APPROACH: Process only the first N structures
#
# How it works:
#   - generate_structures creates all variants (e.g., 15 structures)
#   - Only first max_parallel structures are submitted for relaxation
#   - Already-completed structures are automatically skipped (AiiDA caching)
#
# Recommended workflow:
#   1. Start with small max_parallel (e.g., 2-5) to test
#   2. After successful completion, increase to process more
#   3. Scale up based on cluster resources and queue limits
#
# Example progression:
#   Run 1: max_parallel=5  → processes structures 0-4
#   Run 2: max_parallel=10 → processes structures 5-9 (0-4 cached)
#   Run 3: max_parallel=15 → processes structures 10-14 (0-9 cached)
#
max_parallel = 5  # ← START SMALL, scale up after testing

print(f"\n{'='*70}")
print(f"Parallelization: Processing first {max_parallel} structures")
print(f"{'='*70}")
print(f"Note: Increase this value in subsequent runs to process more structures")
print(f"      Already-completed structures are skipped automatically")

# ============================================================================
# WORKFLOW SETUP AND SUBMISSION
# ============================================================================
print(f"\n{'='*70}")
print("Creating production workflow...")
print(f"{'='*70}")

# Build WorkGraph using builder function
wg = build_surface_hydroxylation_workgraph(
    structure=structure,
    surface_params=surface_params,
    code_label=code_label,
    vasp_config=vasp_config,
    options=options,
    max_parallel_jobs=max_parallel,
    name='surface_hydroxylation_production',
)

print("\nSubmitting to AiiDA daemon...")
print("(Make sure daemon is running: verdi daemon status)")

result = wg.submit()

print("\n" + "="*70)
print(f"✓ WORKFLOW SUBMITTED SUCCESSFULLY")
print("="*70)
print(f"\nWorkflow PK: {result.pk}")
print(f"\nExpected behavior:")
print(f"  1. generate_structures creates all structure variants")
print(f"  2. First {max_parallel} structures submitted for VASP relaxation")
print(f"  3. Relaxations run in parallel on cluster")
print(f"  4. collect_results organizes outputs")
print(f"  5. Workflow completes with exit status [0]")
print(f"\nExpected runtime: Several hours (depends on system size and queue)")
print(f"\n{'='*70}")
print("MONITORING COMMANDS")
print("="*70)
print(f"\nCheck workflow status:")
print(f"  verdi process show {result.pk}")
print(f"  verdi process report {result.pk}")
print(f"\nList all running processes:")
print(f"  verdi process list")
print(f"\nWatch progress (updates every 30s):")
print(f"  watch -n 30 verdi process show {result.pk}")
print(f"\nCheck daemon status:")
print(f"  verdi daemon status")
print(f"\n{'='*70}")
print("ANALYZING RESULTS (after completion)")
print("="*70)
print(f"\nView outputs:")
print(f"  verdi process show {result.pk}")
print(f"\nAccess data in Python:")
print(f"  from aiida import orm")
print(f"  node = orm.load_node({result.pk})")
print(f"  ")
print(f"  # Summary statistics")
print(f"  print(node.outputs.statistics.get_dict())")
print(f"  ")
print(f"  # Successful relaxations")
print(f"  results = node.outputs.successful_relaxations.get_dict()")
print(f"  for r in results['results']:")
print(f"      print(f\"{{r['name']}}: {{r['energy']}} eV @ {{r['coverage']}} coverage\")")
print(f"      structure = orm.load_node(r['structure_pk'])")
print(f"  ")
print(f"  # Failed relaxations (if any)")
print(f"  if 'failed_relaxations' in node.outputs:")
print(f"      failed = node.outputs.failed_relaxations.get_dict()")
print(f"      print(f\"Failed: {{failed['results']}}\")")
print(f"\n{'='*70}")
print("SCALING UP")
print("="*70)
print(f"\nAfter successful completion with max_parallel={max_parallel}:")
print(f"  1. Edit this file and increase max_parallel (e.g., to {max_parallel*2})")
print(f"  2. Run again - already-completed structures will be skipped")
print(f"  3. Only new structures will be processed")
print(f"  4. Repeat until all desired structures are processed")
print("="*70 + "\n")
