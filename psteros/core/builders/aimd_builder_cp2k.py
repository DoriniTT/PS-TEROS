"""
AIMD Builder for CP2K

Material-agnostic builder for ab initio molecular dynamics calculations
using CP2K.
"""


def get_aimd_defaults_cp2k(
    cutoff: float = 400,
    rel_cutoff: float = 60,
    timestep: float = 1.0,
    eps_scf: float = 1e-6,
    max_scf: int = 40,
    thermostat: str = "NOSE",
) -> dict:
    """
    Get default CP2K AIMD parameters.

    Returns sensible defaults for NVT molecular dynamics using NOSE thermostat.

    Args:
        cutoff: CUTOFF value in Ry. Default: 400
        rel_cutoff: REL_CUTOFF value in Ry. Default: 60
        timestep: TIMESTEP in fs. Default: 1.0
        eps_scf: SCF convergence threshold. Default: 1e-6
        max_scf: Maximum SCF iterations. Default: 40
        thermostat: Thermostat type (NOSE, CSVR, etc.). Default: NOSE

    Returns:
        Dictionary with CP2K AIMD parameters. Does NOT include:
        - TEMPERATURE (set automatically per stage)
        - STEPS (set automatically per stage)
        - KIND section (must be added separately for each system)
        - FIXED_ATOMS (use fixed_atoms module if needed)

    Example:
        >>> aimd_params = get_aimd_defaults_cp2k(timestep=2.0)
        >>> # Add KIND section for your system
        >>> aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        ...     {"_": "Ag", "BASIS_SET": "...", "POTENTIAL": "..."},
        ...     {"_": "O", "BASIS_SET": "...", "POTENTIAL": "..."}
        ... ]
        >>> aimd_sequence = [
        ...     {'temperature': 300, 'steps': 1000},
        ...     {'temperature': 500, 'steps': 1000},
        ... ]
        >>> wg = build_core_workgraph(
        ...     calculator='cp2k',
        ...     aimd_sequence=aimd_sequence,
        ...     aimd_parameters=aimd_params,
        ...     **other_params
        ... )
    """
    return {
        "GLOBAL": {
            "RUN_TYPE": "MD",
            "PRINT_LEVEL": "LOW"
        },
        "MOTION": {
            "MD": {
                "ENSEMBLE": "NVT",
                "TIMESTEP": timestep,
                "THERMOSTAT": {
                    "TYPE": thermostat,
                    "REGION": "GLOBAL",
                },
            },
        },
        "FORCE_EVAL": {
            "METHOD": "QS",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "MULTIPLICITY": 1,
                "MGRID": {
                    "CUTOFF": cutoff,
                    "REL_CUTOFF": rel_cutoff,
                    "NGRIDS": 4
                },
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"},
                    "VDW_POTENTIAL": {
                        "POTENTIAL_TYPE": "PAIR_POTENTIAL",
                        "PAIR_POTENTIAL": {
                            "TYPE": "DFTD3",
                            "REFERENCE_FUNCTIONAL": "PBE",
                            "PARAMETER_FILE_NAME": "dftd3.dat"
                        }
                    }
                },
                "POISSON": {
                    "PERIODIC": 'XYZ',
                    "PSOLVER": "PERIODIC",
                },
                "SURFACE_DIPOLE_CORRECTION": ".TRUE.",
                "SCF": {
                    "IGNORE_CONVERGENCE_FAILURE": ".TRUE.",
                    "SCF_GUESS": "ATOMIC",
                    "EPS_SCF": eps_scf,
                    "MAX_SCF": max_scf,
                    "OT": {
                        "PRECONDITIONER": "FULL_SINGLE_INVERSE",
                        "MINIMIZER": "DIIS"
                    },
                    "OUTER_SCF": {
                        "MAX_SCF": 10,
                        "EPS_SCF": eps_scf
                    }
                },
                "QS": {
                    "METHOD": "GPW",
                    "EPS_DEFAULT": 1.0e-12,
                    "EXTRAPOLATION": "ASPC",
                    "EXTRAPOLATION_ORDER": 3
                }
            },
            "SUBSYS": {}  # KIND section must be added by user
        },
    }


def get_basis_molopt_content() -> str:
    """
    Returns hardcoded BASIS_MOLOPT content for common elements (H, O, P, Ag).

    Returns:
        String containing BASIS_MOLOPT file content
    """
    return """# URL: https://cp2k-basis.pierrebeaujean.net/api/basis/DZVP-MOLOPT-PBE-GTH/data?elements=H,O,P,Ag
# BUILD: 06/08/2024 @ 14:19
# FETCHED: 08/08/2025 @ 21:53
# ---
# H [10s5p|2s1p]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/BASIS_MOLOPT_UZH#L38
H  DZVP-MOLOPT-PBE-GTH-q1 DZVP-MOLOPT-GGA-GTH-q1
1
2 0 1 5 2 1
  9.586641744358  0.028188454243 -0.011916098826 -0.104299547500
  2.202359864130  0.137112096750 -0.048501541876  0.333755971437
  0.604094259906  0.421114159535 -0.107830568846  0.500163684961
  0.146497785045  0.825868105506 -0.330879704789  0.775957173703
  0.139525796829  0.347865521314  0.936160668051 -0.159547199117
# O [10s10p5d|2s2p1d]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/BASIS_MOLOPT_UZH#L309
O  DZVP-MOLOPT-PBE-GTH-q6 DZVP-MOLOPT-GGA-GTH-q6
1
2 0 2 5 2 2 1
 10.277042424065 -0.149868256025  0.048759688626 -0.092477972588  0.069652198011 -0.067753697289
  3.564052538806 -0.150283608565  0.038882119766 -0.290527105694  0.215775427853  0.182303508691
  1.313635632268  0.539068797290 -0.125699044360 -0.518398147192  0.367414070471  0.334682703086
  0.488903246714  0.802346121883 -0.075751253852 -0.626250099252  0.267796450585  0.910370319697
  0.155331533694  0.143526566604  0.987204219660 -0.496100795824 -0.861325430904 -0.146247176017
# P [8s8p4d|2s2p1d]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/BASIS_MOLOPT_UZH#L551
P  DZVP-MOLOPT-PBE-GTH-q5 DZVP-MOLOPT-GGA-GTH-q5
1
2 0 2 4 2 2 1
  1.529974428103  0.372140622492  0.294912996621  0.142007115982  0.102489790742  0.258356603330
  0.572156237343 -0.286654463183 -0.354714954359 -0.690107171500 -0.444463997488  0.786129849829
  0.221472852050 -0.855355955676  0.020683118630 -0.704204434492  0.354647334824  0.561466877411
  0.075885075902 -0.218418783347  0.887003852405 -0.087648076394  0.816194134770  0.002582691280
# Ag [10s10p10d5f|2s2p2d1f]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/BASIS_MOLOPT_UZH#L1600
Ag  DZVP-MOLOPT-PBE-GTH-q11 DZVP-MOLOPT-GGA-GTH-q11
1
2 0 3 5 2 2 2 1
  2.495281378251  0.074382968060 -0.050379970111 -0.032009065730 -0.074229383961 -0.403506355709  0.047073978724 -0.010308846444
  1.147928349138 -0.202979901206  0.143630028949  0.054577776230  0.246610290741 -0.676447605213 -0.224013654162 -0.069104090809
  0.461374479372 -0.370367808151  0.058891282271  0.476571057807  0.302330644307 -0.550328606307  0.419980371593 -0.770987410161
  0.172597307993  0.776846732626 -0.450530057203  0.001576332235 -0.899775431851 -0.276707865274 -0.393547443116 -0.367828851362
  0.051247258073  0.461089117702  0.877716786797  0.876854736025 -0.180758725397 -0.013131667027  0.785072493958 -0.515168614919
"""


def get_gth_potentials_content() -> str:
    """
    Returns hardcoded GTH_POTENTIALS content for common elements (H, O, P, Ag).

    Returns:
        String containing GTH_POTENTIALS file content
    """
    return """# URL: https://cp2k-basis.pierrebeaujean.net/api/pseudopotentials/GTH-PBE/data?elements=H,O,P,Ag
# BUILD: 06/08/2024 @ 14:19
# FETCHED: 08/08/2025 @ 21:53
# ---
# H [0|1s]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/POTENTIAL_UZH#L985
H  GTH-PBE-q1 GTH-GGA-q1
1 0 0 0
      0.20059317   2    -4.17806832     0.72440924
     0
# O [2|2s4p]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/POTENTIAL_UZH#L1037
O  GTH-PBE-q6 GTH-GGA-q6
2 4 0 0
      0.24446328   2   -16.67548222     2.48908598
     1
      0.22097111   1    18.33446866
# P [10|2s3p]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/POTENTIAL_UZH#L1103
P  GTH-PBE-q5 GTH-GGA-q5
2 3 0 0
      0.43012670   1    -5.86287518
     2
      0.39637888   2    11.00906771    -3.47035684
                                       4.48022987
      0.44829366   1     3.05605781
# Ag [36|1s10d]
# SOURCE: https://github.com/cp2k/cp2k/raw/ac0226eb549c7ef1ea50d0597d545f29d4c8fc87/data/POTENTIAL_UZH#L1831
Ag  GTH-PBE-q11 GTH-GGA-q11
1 0 10 0
      0.57261106   1    -0.09385803
     3
      0.52717235   3     9.59197051    -5.27424307     0.99704951
                                       8.43011316    -2.57436732
                                                      2.02769689
      0.62063433   2     3.90685517    -1.68549606
                                       2.06129973
      0.39996164   2    -2.69173033    -0.43354834
                                       0.39191109
"""


def prepare_aimd_parameters_cp2k(
    base_parameters: dict,
    temperature: float,
    steps: int,
) -> dict:
    """
    Inject temperature and steps into CP2K AIMD parameters.

    Similar to VASP's prepare_aimd_parameters but for CP2K.

    Args:
        base_parameters: Base CP2K parameters from get_aimd_defaults_cp2k()
        temperature: Target temperature in K
        steps: Number of MD steps for this stage

    Returns:
        Complete CP2K parameters dict for this AIMD stage
    """
    import copy
    params = copy.deepcopy(base_parameters)
    params['MOTION']['MD']['TEMPERATURE'] = temperature
    params['MOTION']['MD']['STEPS'] = steps
    return params
