"""
AIMD Builder

Material-agnostic builder for ab initio molecular dynamics calculations
using VASP.
"""

def get_aimd_defaults(
    energy_cutoff: float = 400,
    electronic_convergence: float = 1e-5,
    timestep: float = 1.0,
    smass: float = 0.0,
    mdalgo: int = 2,
    ismear: int = 0,
    sigma: float = 0.1,
    isym: int = 0,
    ncore: int = 4,
    lreal: str = "Auto",
) -> dict:
    """
    Get default parameters for VASP AIMD calculations.

    Returns sensible defaults for NVT molecular dynamics using Nosé-Hoover
    thermostat (MDALGO=2), similar to CP2K NOSE thermostat.

    Args:
        energy_cutoff: ENCUT value (eV). Default: 400
        electronic_convergence: EDIFF value. Default: 1e-5
        timestep: POTIM timestep in fs. Default: 1.0
        smass: Nosé mass parameter. Default: 0.0 (automatic)
        mdalgo: MD algorithm (2=NVT Nosé-Hoover). Default: 2
        ismear: Smearing method. Default: 0 (Gaussian)
        sigma: Smearing width (eV). Default: 0.1
        isym: Symmetry (0=off for MD). Default: 0
        ncore: Cores per band. Default: 4
        lreal: Projection operators. Default: "Auto"

    Returns:
        Dictionary with AIMD INCAR parameters. Does NOT include:
        - TEBEG/TEEND (set automatically per stage)
        - NSW (set automatically per stage)

    Example:
        >>> aimd_params = get_aimd_defaults(timestep=2.0)
        >>> aimd_sequence = [
        ...     {'temperature': 300, 'steps': 1000},
        ...     {'temperature': 500, 'steps': 1000},
        ... ]
        >>> wg = build_core_workgraph(
        ...     run_aimd=True,
        ...     aimd_sequence=aimd_sequence,
        ...     aimd_parameters=aimd_params,
        ...     **other_params
        ... )
    """
    return {
        'PREC': 'Normal',
        'ENCUT': energy_cutoff,
        'EDIFF': electronic_convergence,
        'ISMEAR': ismear,
        'SIGMA': sigma,
        'IBRION': 0,
        'MDALGO': mdalgo,
        'SMASS': smass,
        'POTIM': timestep,
        'ISYM': isym,
        'LREAL': lreal,
        'NCORE': ncore,
        'LWAVE': True,
        'LCHARG': True,
    }
