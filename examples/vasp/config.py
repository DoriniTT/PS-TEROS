"""
Shared configuration for PS-TEROS example scripts.

This file contains common settings that can be imported by all step examples.
Modify these values to match your cluster and AiiDA setup.

Usage in example scripts:
    from config import (
        AIIDA_PROFILE, VASP_CODE, POTENTIAL_FAMILY,
        BULK_PARAMS, SLAB_PARAMS, COMMON_OPTIONS
    )
"""

# ==============================================================================
# AIIDA CONFIGURATION
# ==============================================================================

# AiiDA profile name (check with: verdi profile list)
AIIDA_PROFILE = 'presto'

# ==============================================================================
# VASP CODE CONFIGURATION
# ==============================================================================

# VASP code label (check with: verdi code list)
# Common formats: 'VASP-6.4.1@cluster', 'vasp_std@supercomputer'
VASP_CODE = 'VASP-6.4.1@cluster'

# PAW potential family (check with: verdi data vasp-potcar listfamilies)
POTENTIAL_FAMILY = 'PBE'

# ==============================================================================
# SCHEDULER OPTIONS
# ==============================================================================

# Default scheduler options for VASP calculations
# Adjust num_cores_per_machine and queue_name for your cluster
COMMON_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 24,
    },
    'queue_name': 'normal',
    # 'max_wallclock_seconds': 86400,  # Optional: 24 hours
}

# Options with more cores (for larger calculations)
LARGE_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 40,
    },
    'queue_name': 'par40',
}

# ==============================================================================
# DEFAULT VASP PARAMETERS
# ==============================================================================

# Standard bulk relaxation parameters (ISIF=3: full cell + ions)
BULK_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 3,          # Full relaxation (cell + ions)
    'NSW': 100,
    'EDIFFG': -0.01,
    'ALGO': 'Normal',
    'LREAL': 'Auto',
    'LWAVE': False,
    'LCHARG': False,
}

# Standard slab relaxation parameters (ISIF=2: ions only)
SLAB_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 2,          # Ions only (fixed cell)
    'NSW': 100,
    'EDIFFG': -0.01,
    'ALGO': 'Normal',
    'LREAL': 'Auto',
    'LWAVE': False,
    'LCHARG': False,
}

# Metal-specific parameters (different smearing)
METAL_PARAMS = {
    'PREC': 'Accurate',
    'ENCUT': 520,
    'EDIFF': 1e-6,
    'ISMEAR': 1,        # Methfessel-Paxton for metals
    'SIGMA': 0.2,
    'IBRION': 2,
    'ISIF': 3,
    'NSW': 100,
    'EDIFFG': -0.01,
    'ALGO': 'Normal',
    'LREAL': 'Auto',
    'LWAVE': False,
    'LCHARG': False,
}

# Light parameters for testing (faster but less accurate)
LIGHT_PARAMS = {
    'PREC': 'Normal',
    'ENCUT': 400,
    'EDIFF': 1e-4,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IBRION': 2,
    'ISIF': 2,
    'NSW': 50,
    'EDIFFG': -0.05,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'LWAVE': False,
    'LCHARG': False,
}

# ==============================================================================
# K-POINT SETTINGS
# ==============================================================================

# Standard k-point spacings (in 1/Angstrom)
KPOINTS_FINE = 0.03      # High accuracy
KPOINTS_STANDARD = 0.04  # Standard calculations
KPOINTS_COARSE = 0.06    # Faster calculations
KPOINTS_AIMD = 0.5       # AIMD (Gamma-centered typically)

# ==============================================================================
# POTENTIAL MAPPINGS
# ==============================================================================

# Common element mappings (element -> PAW potential name)
# Most elements map to themselves, but some have variants (e.g., 'Ag_pv')
POTENTIAL_MAPPING_AG2O = {'Ag': 'Ag', 'O': 'O'}
POTENTIAL_MAPPING_AG3PO4 = {'Ag': 'Ag', 'P': 'P', 'O': 'O'}
POTENTIAL_MAPPING_AG = {'Ag': 'Ag'}
POTENTIAL_MAPPING_O2 = {'O': 'O'}

# ==============================================================================
# SLAB GENERATION DEFAULTS
# ==============================================================================

SLAB_DEFAULTS = {
    'min_slab_thickness': 18.0,    # Angstroms
    'min_vacuum_thickness': 15.0,  # Angstroms
    'lll_reduce': True,
    'center_slab': True,
    'symmetrize': True,
    'primitive': True,
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def get_slab_params(base_params=None):
    """Get slab parameters (ISIF=2) from bulk parameters."""
    params = (base_params or BULK_PARAMS).copy()
    params['ISIF'] = 2
    return params


def get_structures_dir():
    """Get the path to the structures directory."""
    import os
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'structures')


def load_profile():
    """Load the configured AiiDA profile."""
    from aiida import load_profile as aiida_load_profile
    aiida_load_profile(profile=AIIDA_PROFILE)
