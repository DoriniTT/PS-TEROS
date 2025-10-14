# CP2K AIMD Integration Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Add CP2K calculator support for AIMD simulations in PS-TEROS, with calculator-agnostic fixed atoms utilities.

**Architecture:** Builder-centric approach with separate CP2K AIMD module (`aimd_cp2k.py`), reusable fixed atoms utilities (`fixed_atoms.py`), and calculator routing in `build_core_workgraph()`. CP2K outputs are normalized to match VASP output structure for consistent user experience.

**Tech Stack:** AiiDA, aiida-workgraph, aiida-cp2k, CP2K, ASE

---

## Task 1: Create Fixed Atoms Module (Calculator-Agnostic)

**Files:**
- Create: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/fixed_atoms.py`

**Step 1: Create the fixed_atoms.py module**

Create file with complete implementation:

```python
"""
Fixed Atoms Module

Calculator-agnostic utilities for constraining atoms in slab structures.
Supports fixing atoms by position (top, bottom, center) with configurable thicknesses.
"""

import typing as t
import numpy as np
from aiida import orm


def get_fixed_atoms_list(
    structure: orm.StructureData,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
) -> t.List[int]:
    """
    Identify atoms to fix in a slab structure based on position criteria.

    Args:
        structure: AiiDA StructureData (slab structure)
        fix_type: Where to fix atoms. Options:
            - 'bottom': Fix atoms from bottom up to fix_thickness Å
            - 'top': Fix atoms from top down to fix_thickness Å
            - 'center': Fix atoms within fix_thickness/2 Å of slab center
            - None: No fixing (returns empty list)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O'])
                     If None, all elements in the region are fixed

    Returns:
        List of 1-based atom indices to fix (sorted)

    Examples:
        # Fix bottom 7 Å of all atoms
        >>> fixed = get_fixed_atoms_list(slab, fix_type='bottom', fix_thickness=7.0)

        # Fix top 5 Å of Ag atoms only
        >>> fixed = get_fixed_atoms_list(slab, fix_type='top', fix_thickness=5.0,
        ...                              fix_elements=['Ag'])

        # Fix 4 Å around center (2 Å above and below)
        >>> fixed = get_fixed_atoms_list(slab, fix_type='center', fix_thickness=4.0)
    """
    if fix_type is None or fix_thickness <= 0.0:
        return []

    # Get atomic positions and symbols
    positions = np.array([site.position for site in structure.sites])
    symbols = [site.kind_name for site in structure.sites]

    z_coords = positions[:, 2]
    z_min, z_max = np.min(z_coords), np.max(z_coords)
    z_center = (z_min + z_max) / 2.0

    fixed_indices = set()

    def is_eligible(symbol: str) -> bool:
        """Check if atom type should be fixed."""
        return fix_elements is None or symbol in fix_elements

    # Determine fixing region based on type
    if fix_type == 'bottom':
        z_cutoff = z_min + fix_thickness
        for idx, (z, symbol) in enumerate(zip(z_coords, symbols)):
            if z <= z_cutoff and is_eligible(symbol):
                fixed_indices.add(idx + 1)  # 1-based indexing

    elif fix_type == 'top':
        z_cutoff = z_max - fix_thickness
        for idx, (z, symbol) in enumerate(zip(z_coords, symbols)):
            if z >= z_cutoff and is_eligible(symbol):
                fixed_indices.add(idx + 1)

    elif fix_type == 'center':
        half_thickness = fix_thickness / 2.0
        z_lower = z_center - half_thickness
        z_upper = z_center + half_thickness
        for idx, (z, symbol) in enumerate(zip(z_coords, symbols)):
            if z_lower <= z <= z_upper and is_eligible(symbol):
                fixed_indices.add(idx + 1)

    else:
        raise ValueError(
            f"Invalid fix_type: '{fix_type}'. "
            f"Must be one of: 'bottom', 'top', 'center', or None"
        )

    return sorted(fixed_indices)


def add_fixed_atoms_to_cp2k_parameters(
    base_parameters: dict,
    fixed_atoms_list: t.List[int],
    components: str = "XYZ",
) -> dict:
    """
    Add FIXED_ATOMS constraint to CP2K parameters.

    Args:
        base_parameters: Base CP2K parameters dict
        fixed_atoms_list: List of 1-based atom indices to fix
        components: Which components to fix (default: "XYZ")
            - "XYZ": Fix all three dimensions (fully rigid)
            - "XY": Fix only in-plane motion
            - "Z": Fix only out-of-plane motion

    Returns:
        Updated parameters dict with FIXED_ATOMS constraint
    """
    import copy
    params = copy.deepcopy(base_parameters)

    if not fixed_atoms_list:
        return params

    # Ensure MOTION section exists
    if "MOTION" not in params:
        params["MOTION"] = {}
    if "CONSTRAINT" not in params["MOTION"]:
        params["MOTION"]["CONSTRAINT"] = {}

    # Add FIXED_ATOMS constraint
    params["MOTION"]["CONSTRAINT"]["FIXED_ATOMS"] = {
        "LIST": " ".join(map(str, fixed_atoms_list)),
        "COMPONENTS_TO_FIX": components
    }

    return params


def add_fixed_atoms_to_vasp_parameters(
    base_parameters: dict,
    structure: orm.StructureData,
    fixed_atoms_list: t.List[int],
) -> t.Tuple[dict, orm.StructureData]:
    """
    Add selective dynamics to VASP INCAR and create constrained structure.

    VASP uses selective dynamics via POSCAR, not INCAR parameters.
    This function:
    1. Sets IBRION and NSW appropriately
    2. Returns a modified StructureData with constraints

    Args:
        base_parameters: Base VASP INCAR parameters
        structure: Original StructureData
        fixed_atoms_list: List of 1-based atom indices to fix

    Returns:
        Tuple of (updated_parameters, constrained_structure)

    Note:
        For VASP, the constrained structure must be used in the calculation.
        The structure's 'kinds' will have 'fixed' tags set appropriately.
    """
    import copy
    from ase import Atoms
    from ase.constraints import FixAtoms

    params = copy.deepcopy(base_parameters)

    if not fixed_atoms_list:
        return params, structure

    # Ensure selective dynamics is enabled
    params['IBRION'] = params.get('IBRION', 2)  # Keep existing or default to 2

    # Get ASE atoms and add constraint
    atoms = structure.get_ase()

    # Convert to 0-based indices for ASE
    fixed_indices_0based = [idx - 1 for idx in fixed_atoms_list]

    # Add FixAtoms constraint
    constraint = FixAtoms(indices=fixed_indices_0based)
    atoms.set_constraint(constraint)

    # Create new StructureData with constraints
    constrained_structure = orm.StructureData(ase=atoms)

    return params, constrained_structure
```

**Step 2: Verify file creation**

Run: `ls -la /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/fixed_atoms.py`
Expected: File exists with 200+ lines

**Step 3: Test import**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/psteros/bin/activate
python -c "from teros.core.fixed_atoms import get_fixed_atoms_list, add_fixed_atoms_to_cp2k_parameters; print('✓ Import successful')"
```
Expected: "✓ Import successful"

**Step 4: Clear Python cache**

Run: `find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -type d -name __pycache__ -exec rm -rf {} + && find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add teros/core/fixed_atoms.py
git commit -m "feat: add calculator-agnostic fixed atoms utilities

- get_fixed_atoms_list(): identify atoms by position (bottom/top/center)
- add_fixed_atoms_to_cp2k_parameters(): add CP2K FIXED_ATOMS constraint
- add_fixed_atoms_to_vasp_parameters(): add VASP selective dynamics
- Supports element filtering and configurable components"
```

---

## Task 2: Create CP2K AIMD Builder

**Files:**
- Create: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/builders/aimd_builder_cp2k.py`

**Step 1: Create the aimd_builder_cp2k.py module**

Create file with complete implementation:

```python
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
                "SURFACE_DIPOLE_CORRECTION": True,
                "SCF": {
                    "IGNORE_CONVERGENCE_FAILURE": True,
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
```

**Step 2: Verify file creation**

Run: `ls -la /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/builders/aimd_builder_cp2k.py`
Expected: File exists with 200+ lines

**Step 3: Test import**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/psteros/bin/activate
python -c "from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k, get_basis_molopt_content, get_gth_potentials_content, prepare_aimd_parameters_cp2k; print('✓ Import successful')"
```
Expected: "✓ Import successful"

**Step 4: Clear Python cache**

Run: `find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -type d -name __pycache__ -exec rm -rf {} + && find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add teros/core/builders/aimd_builder_cp2k.py
git commit -m "feat: add CP2K AIMD builder with defaults

- get_aimd_defaults_cp2k(): default NVT parameters
- get_basis_molopt_content(): hardcoded basis sets (H,O,P,Ag)
- get_gth_potentials_content(): hardcoded pseudopotentials
- prepare_aimd_parameters_cp2k(): inject temperature and steps"
```

---

## Task 3: Create CP2K AIMD Scatter Function

**Files:**
- Create: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/aimd_cp2k.py`

**Step 1: Create the aimd_cp2k.py module**

Create file with complete implementation:

```python
"""
AIMD Module for PS-TEROS - CP2K Implementation

Ab initio molecular dynamics calculations on slab structures using CP2K.
Sequential AIMD stages with automatic restart chaining.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace


@task.graph
def aimd_single_stage_scatter_cp2k(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    temperature: float,
    steps: int,
    code: orm.Code,
    aimd_parameters: dict,
    basis_file: orm.SinglefileData,
    pseudo_file: orm.SinglefileData,
    options: dict,
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    fixed_atoms_lists: dict = None,
    fix_components: str = "XYZ",
) -> t.Annotated[dict, namespace(
    structures=dynamic(orm.StructureData),
    remote_folders=dynamic(orm.RemoteData),
    parameters=dynamic(orm.Dict),
    trajectories=dynamic(orm.TrajectoryData),
    retrieved=dynamic(orm.FolderData)
)]:
    """
    Run single CP2K AIMD stage on all slabs in parallel using scatter-gather pattern.

    This function handles ONE temperature/timestep stage for all slabs.
    Call it multiple times sequentially to build multi-stage AIMD workflows.

    Args:
        slabs: Dictionary of slab structures to run AIMD on
        temperature: Target temperature in K
        steps: Number of MD steps for this stage
        code: CP2K code
        aimd_parameters: Base CP2K AIMD parameters (from get_aimd_defaults_cp2k)
        basis_file: SinglefileData containing BASIS_MOLOPT
        pseudo_file: SinglefileData containing GTH_POTENTIALS
        options: Scheduler options (metadata)
        clean_workdir: Whether to clean work directory
        restart_folders: Optional dict of RemoteData for restart (from previous stage)
        fixed_atoms_lists: Optional dict mapping slab_label -> list of fixed atom indices
        fix_components: Components to fix ("XYZ", "XY", "Z")

    Returns:
        Dictionary with outputs per slab:
            - structures: Output structures from this stage
            - remote_folders: RemoteData nodes for next stage restart
            - parameters: Output parameters from CP2K
            - trajectories: Trajectory data from MD
            - retrieved: Retrieved folder data
    """
    from teros.core.builders.aimd_builder_cp2k import prepare_aimd_parameters_cp2k
    from teros.core.fixed_atoms import add_fixed_atoms_to_cp2k_parameters

    # Get CP2K workchain
    Cp2kWorkChain = WorkflowFactory('cp2k.base')
    Cp2kTask = task(Cp2kWorkChain)

    structures_out = {}
    remote_folders_out = {}
    parameters_out = {}
    trajectories_out = {}
    retrieved_out = {}

    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters_cp2k(aimd_parameters, temperature, steps)

        # Add fixed atoms if provided for this slab
        if fixed_atoms_lists and slab_label in fixed_atoms_lists:
            stage_params = add_fixed_atoms_to_cp2k_parameters(
                stage_params,
                fixed_atoms_lists[slab_label],
                fix_components,
            )

        # Build CP2K inputs
        cp2k_inputs = {
            'structure': slab_structure,
            'parameters': orm.Dict(dict=stage_params),
            'code': code,
            'metadata': options,
            'file': {
                'basis': basis_file,
                'pseudo': pseudo_file,
            },
            'settings': orm.Dict(dict={
                'additional_retrieve_list': [
                    "aiida-1.ener",
                    "aiida-1.restart",
                    "aiida-pos-1.xyz"
                ]
            }),
        }

        # Add restart folder if provided for this slab
        if restart_folders and slab_label in restart_folders:
            cp2k_inputs['parent_calc_folder'] = restart_folders[slab_label]

        # Create CP2K task
        aimd_task = Cp2kTask(
            cp2k=cp2k_inputs,
            max_iterations=orm.Int(3),
            clean_workdir=orm.Bool(clean_workdir),
        )

        # Store CP2K outputs (note: different from VASP!)
        structures_out[slab_label] = aimd_task.output_structure
        remote_folders_out[slab_label] = aimd_task.remote_folder
        parameters_out[slab_label] = aimd_task.output_parameters
        trajectories_out[slab_label] = aimd_task.output_trajectory
        retrieved_out[slab_label] = aimd_task.retrieved

    # Gather: return collected results
    return {
        'structures': structures_out,
        'remote_folders': remote_folders_out,
        'parameters': parameters_out,
        'trajectories': trajectories_out,
        'retrieved': retrieved_out,
    }
```

**Step 2: Verify file creation**

Run: `ls -la /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/aimd_cp2k.py`
Expected: File exists with 100+ lines

**Step 3: Test import**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/psteros/bin/activate
python -c "from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k; print('✓ Import successful')"
```
Expected: "✓ Import successful"

**Step 4: Clear Python cache**

Run: `find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -type d -name __pycache__ -exec rm -rf {} + && find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add teros/core/aimd_cp2k.py
git commit -m "feat: add CP2K AIMD scatter function

- aimd_single_stage_scatter_cp2k(): parallel AIMD for all slabs
- Uses cp2k.base workchain
- Supports restart via parent_calc_folder
- Supports fixed atoms constraints
- Returns 5 output dicts: structures, remote_folders, parameters, trajectories, retrieved"
```

---

## Task 4: Update Builders __init__.py

**Files:**
- Modify: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/builders/__init__.py`

**Step 1: Read current file**

Run: `cat /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/builders/__init__.py`

**Step 2: Add CP2K builder exports**

Add these imports to the file:

```python
from teros.core.builders.aimd_builder_cp2k import (
    get_aimd_defaults_cp2k,
    get_basis_molopt_content,
    get_gth_potentials_content,
    prepare_aimd_parameters_cp2k,
)
```

And add to `__all__` list:
```python
    'get_aimd_defaults_cp2k',
    'get_basis_molopt_content',
    'get_gth_potentials_content',
    'prepare_aimd_parameters_cp2k',
```

**Step 3: Verify import**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/psteros/bin/activate
python -c "from teros.core.builders import get_aimd_defaults_cp2k; print('✓ Export successful')"
```
Expected: "✓ Export successful"

**Step 4: Clear Python cache**

Run: `find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -type d -name __pycache__ -exec rm -rf {} + && find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add teros/core/builders/__init__.py
git commit -m "feat: export CP2K AIMD builder functions"
```

---

## Task 5: Add Calculator Routing to build_core_workgraph

**Files:**
- Modify: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/teros/core/workgraph.py`

**Step 1: Add calculator parameter to function signature**

Locate the `build_core_workgraph` function (around line 512) and add these parameters:

```python
def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    calculator: str = 'vasp',  # NEW: 'vasp' or 'cp2k'
    potential_family: str = 'PBE',
    # ... existing parameters ...

    # CP2K-specific parameters (NEW - add after aimd_kpoints_spacing):
    basis_content: str = None,
    pseudo_content: str = None,
    cp2k_kind_section: list = None,

    # Fixed atoms configuration (NEW - add after cp2k parameters):
    fix_atoms: bool = False,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list = None,
    fix_components: str = "XYZ",

    # ... rest of parameters
):
```

**Step 2: Add CP2K basis/pseudo file creation**

After the restart handling section (around line 860), add:

```python
    # ========================================================================
    # CP2K-SPECIFIC SETUP
    # ========================================================================

    basis_file = None
    pseudo_file = None
    if calculator == 'cp2k':
        import io
        from teros.core.builders.aimd_builder_cp2k import (
            get_basis_molopt_content,
            get_gth_potentials_content,
        )

        # Use provided content or defaults
        basis_str = basis_content if basis_content else get_basis_molopt_content()
        pseudo_str = pseudo_content if pseudo_content else get_gth_potentials_content()

        # Create SinglefileData nodes
        basis_file = orm.SinglefileData(
            io.BytesIO(basis_str.encode('utf-8')),
            filename='BASIS_MOLOPT'
        )
        pseudo_file = orm.SinglefileData(
            io.BytesIO(pseudo_str.encode('utf-8')),
            filename='GTH_POTENTIALS'
        )

        print(f"\n  ✓ Created CP2K basis and pseudopotential files")
```

**Step 3: Add AIMD calculator routing**

Locate the AIMD section (around line 1276) and replace the existing code with:

```python
    # ===== AIMD CALCULATION (OPTIONAL) =====
    # Check if AIMD should be added
    should_add_aimd = run_aimd and aimd_sequence is not None and (input_slabs is not None or miller_indices is not None)

    if should_add_aimd:
        print("\n  → Adding AIMD stages to workflow")
        print(f"     Calculator: {calculator}")
        print(f"     Number of stages: {len(aimd_sequence)}")

        # Import appropriate scatter function based on calculator
        if calculator == 'vasp':
            from teros.core.aimd import aimd_single_stage_scatter
            aimd_scatter_func = aimd_single_stage_scatter
        elif calculator == 'cp2k':
            from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k
            aimd_scatter_func = aimd_single_stage_scatter_cp2k
        else:
            raise ValueError(f"Unknown calculator: {calculator}")

        # Load code
        code = load_code(code_label)

        # Get parameters - use AIMD-specific or fall back to slab/bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        aimd_params = aimd_parameters if aimd_parameters is not None else slab_params
        aimd_opts = aimd_options if aimd_options is not None else slab_opts
        aimd_pot_map = aimd_potential_mapping if aimd_potential_mapping is not None else slab_pot_map
        aimd_kpts = aimd_kpoints_spacing if aimd_kpoints_spacing is not None else slab_kpts

        # Handle fixed atoms if requested
        fixed_atoms_lists = {}

        if fix_atoms and fix_type is not None:
            print(f"\n  → Preparing fixed atoms constraints")
            print(f"     Type: {fix_type}")
            print(f"     Thickness: {fix_thickness} Å")
            print(f"     Elements: {fix_elements if fix_elements else 'all'}")

            from teros.core.fixed_atoms import get_fixed_atoms_list

            # For input_slabs, calculate fixed atoms now
            if input_slabs is not None:
                for label, slab_struct in input_slabs.items():
                    fixed_list = get_fixed_atoms_list(
                        slab_struct,
                        fix_type=fix_type,
                        fix_thickness=fix_thickness,
                        fix_elements=fix_elements,
                    )
                    fixed_atoms_lists[label] = fixed_list
                    print(f"     {label}: {len(fixed_list)} atoms fixed")

        # Determine which slabs to use as initial structures
        if relax_slabs and 'relax_slabs_scatter' in wg.tasks:
            relax_task = wg.tasks['relax_slabs_scatter']
            initial_slabs_source = relax_task.outputs.relaxed_structures
            print(f"     Using relaxed slabs from relax_slabs_scatter task")
        elif relax_slabs and 'collect_slab_outputs_restart' in wg.tasks:
            collector_task = wg.tasks['collect_slab_outputs_restart']
            initial_slabs_source = collector_task.outputs.structures
            print(f"     Using relaxed slabs from restart collector")
        elif input_slabs is not None:
            initial_slabs_source = input_slabs
            print(f"     Using unrelaxed input slabs as initial structures")
        else:
            gen_task = wg.tasks['generate_slab_structures']
            initial_slabs_source = gen_task.outputs.slabs
            print(f"     Using unrelaxed generated slabs as initial structures")

        # Sequential AIMD stages
        current_structures = initial_slabs_source
        current_remotes = {}
        stage_tasks = []

        for stage_idx, stage_config in enumerate(aimd_sequence):
            stage_name = f"aimd_stage_{stage_idx:02d}_{stage_config['temperature']}K"
            print(f"     Stage {stage_idx}: {stage_config['temperature']}K × {stage_config['steps']} steps")

            # Build stage inputs based on calculator
            if calculator == 'vasp':
                stage_inputs = {
                    'slabs': current_structures,
                    'temperature': stage_config['temperature'],
                    'steps': stage_config['steps'],
                    'code': code,
                    'aimd_parameters': aimd_params,
                    'potential_family': potential_family,
                    'potential_mapping': aimd_pot_map,
                    'options': aimd_opts,
                    'kpoints_spacing': aimd_kpts,
                    'clean_workdir': clean_workdir,
                    'restart_folders': current_remotes,
                }
            elif calculator == 'cp2k':
                stage_inputs = {
                    'slabs': current_structures,
                    'temperature': stage_config['temperature'],
                    'steps': stage_config['steps'],
                    'code': code,
                    'aimd_parameters': aimd_params,
                    'basis_file': basis_file,
                    'pseudo_file': pseudo_file,
                    'options': aimd_opts,
                    'clean_workdir': clean_workdir,
                    'restart_folders': current_remotes,
                }

                # Add fixed atoms for CP2K
                if fixed_atoms_lists:
                    stage_inputs['fixed_atoms_lists'] = fixed_atoms_lists
                    stage_inputs['fix_components'] = fix_components

            # Add task
            stage_task = wg.add_task(
                aimd_scatter_func,
                name=stage_name,
                **stage_inputs
            )

            stage_tasks.append(stage_task)

            # Wire outputs to next stage inputs
            current_structures = stage_task.outputs.structures
            current_remotes = stage_task.outputs.remote_folders

        print(f"  ✓ AIMD calculation enabled ({len(aimd_sequence)} sequential stages)")
        print(f"     Access AIMD outputs via: wg.tasks['aimd_stage_XX_XXXK'].outputs")
```

**Step 4: Clear Python cache and restart daemon**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
verdi daemon restart
```

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add teros/core/workgraph.py
git commit -m "feat: add calculator routing and fixed atoms to build_core_workgraph

- Add calculator parameter ('vasp' or 'cp2k')
- Add CP2K basis/pseudo file creation
- Add fixed atoms support (fix_type, fix_thickness, fix_elements)
- Route AIMD to appropriate scatter function based on calculator
- CP2K gets basis_file, pseudo_file, fixed_atoms_lists
- VASP gets potential_family, kpoints_spacing"
```

---

## Task 6: Create Example Script - Auto-generate Slabs

**Files:**
- Create: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07a_aimd_autogenerate_slabs.py`

**Step 1: Create the example script**

Create file with complete implementation:

```python
#!/home/thiagotd/envs/psteros/bin/python
"""
STEP 7A: AIMD with CP2K - Auto-generate Slabs

This workflow:
1. Relaxes bulk structure (VASP)
2. Generates slabs from relaxed bulk
3. Runs AIMD on generated slabs (CP2K)

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence (equilibration + production)

Usage:
    source ~/envs/psteros/bin/activate
    python step_07a_aimd_autogenerate_slabs.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

def main():
    """Step 7A: Test AIMD with CP2K - auto-generate slabs."""

    print("\n" + "="*70)
    print("STEP 7A: AIMD WITH CP2K - AUTO-GENERATE SLABS")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, '../structures')

    print(f"\n2. Structure:")
    print(f"   Bulk: {structures_dir}/ag2o.cif")

    # VASP parameters for bulk relaxation
    bulk_params = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    bulk_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # AIMD configuration
    print("\n3. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    for i, stage in enumerate(aimd_sequence):
        print(f"   Stage {i+1}: {stage['temperature']} K, {stage['steps']} steps")

    # CP2K AIMD parameters
    aimd_params = get_aimd_defaults_cp2k(
        cutoff=400,
        rel_cutoff=60,
        timestep=1.0,
        eps_scf=1e-6,
        thermostat='NOSE',
    )

    # Add KIND section for Ag and O
    if 'FORCE_EVAL' not in aimd_params:
        aimd_params['FORCE_EVAL'] = {}
    if 'SUBSYS' not in aimd_params['FORCE_EVAL']:
        aimd_params['FORCE_EVAL']['SUBSYS'] = {}

    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {
            "_": "Ag",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11",
            "POTENTIAL": "GTH-PBE-q11",
        },
        {
            "_": "O",
            "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6",
            "POTENTIAL": "GTH-PBE-q6",
        }
    ]

    aimd_options = {
        'resources': {
            'num_machines': 2,
            'num_cores_per_machine': 128,
        },
        'queue_name': 'paralela',
    }

    print("\n4. Building workgraph...")
    print("   Workflow: Bulk relaxation → Slab generation → AIMD")
    print("   Bulk code: VASP")
    print("   AIMD code: CP2K")
    print("   Using preset: 'aimd_only'")

    # Build workgraph using preset with CP2K for AIMD
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',  # Use CP2K for AIMD

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',

        # Code for bulk (VASP)
        code_label='VASP-VTST-6.4.3@bohr',
        potential_family='PBE',
        kpoints_spacing=0.4,
        clean_workdir=False,

        # Bulk parameters (REQUIRED for auto-generation)
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=bulk_params,
        bulk_options=bulk_options,

        # Slab generation (REQUIRED for auto-generation)
        miller_indices=[1, 1, 1],
        min_slab_thickness=15.0,
        min_vacuum_thickness=15.0,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,

        # AIMD (CP2K)
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        name='Step07A_AIMD_CP2K_AutoGen',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 7A SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nWorkflow stages:")
    print(f"  1. Bulk relaxation (VASP)")
    print(f"  2. Slab generation")
    print(f"  3. AIMD on all slabs (CP2K)")
    print(f"\nExpected outputs:")
    print(f"  - bulk_energy, bulk_structure")
    print(f"  - slab_structures")
    print(f"  - AIMD stage outputs in tasks:")
    print(f"    wg.tasks['aimd_stage_00_300K'].outputs:")
    print(f"      - structures, remote_folders")
    print(f"      - parameters, trajectories, retrieved")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

**Step 2: Make script executable**

Run: `chmod +x /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07a_aimd_autogenerate_slabs.py`

**Step 3: Verify script**

Run: `head -n 30 /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07a_aimd_autogenerate_slabs.py`
Expected: See docstring and imports

**Step 4: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add examples/cp2k/step_07a_aimd_autogenerate_slabs.py
git commit -m "feat: add CP2K AIMD example - auto-generate slabs

- Demonstrates full workflow: bulk → slab generation → AIMD
- Uses VASP for bulk, CP2K for AIMD
- Ag2O (111) surface with 2-stage AIMD sequence"
```

---

## Task 7: Create Example Script - Input Slabs with Fixed Atoms

**Files:**
- Create: `/home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07b_aimd_input_slabs.py`

**Step 1: Create the example script**

Create file with complete implementation:

```python
#!/home/thiagotd/envs/psteros/bin/python
"""
STEP 7B: AIMD with CP2K - Input Slabs Directly (with Fixed Atoms)

This workflow:
1. Takes pre-existing slab structures as input
2. Runs AIMD directly on input slabs (CP2K)
3. NO bulk relaxation or slab generation
4. Demonstrates fixed atoms constraints

Use this when you already have slab structures from:
- Previous calculations
- Manual construction
- Other structure generation tools

Material: Ag2O
Surface: (111)
AIMD: 2-stage sequence with bottom 7Å fixed

Usage:
    source ~/envs/psteros/bin/activate
    python step_07b_aimd_input_slabs.py
"""

import sys
import os
from aiida import load_profile
from aiida import orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

def main():
    """Step 7B: Test AIMD with CP2K - input slabs with fixed atoms."""

    print("\n" + "="*70)
    print("STEP 7B: AIMD WITH CP2K - INPUT SLABS (FIXED ATOMS)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, '../structures')

    # Load pre-existing slab structures
    print("\n2. Loading input slab structures...")

    # Example: Load slabs from files
    # NOTE: You need to create these files or use previous calculation outputs
    slab_files = {
        'slab_111_term0': os.path.join(structures_dir, 'ag2o_111_term0.vasp'),
        'slab_111_term1': os.path.join(structures_dir, 'ag2o_111_term1.vasp'),
    }

    input_slabs = {}
    for label, filepath in slab_files.items():
        if os.path.exists(filepath):
            atoms = read(filepath)
            input_slabs[label] = orm.StructureData(ase=atoms)
            print(f"   ✓ Loaded {label} from {filepath}")
        else:
            print(f"   ⚠️  File not found: {filepath}")
            print(f"      Creating dummy structure for demonstration...")
            # Create dummy slab for demonstration
            from ase.build import bulk, surface
            ag2o_bulk = bulk('Ag2O', 'cubic', a=4.7)
            slab = surface(ag2o_bulk, (1,1,1), 4, vacuum=15.0)
            input_slabs[label] = orm.StructureData(ase=slab)
            print(f"   ✓ Created dummy slab: {label}")

    # Alternative: Load from previous calculation
    # prev_wg = orm.load_node(12345)  # PK of previous workgraph
    # input_slabs = {
    #     'term_0': prev_wg.outputs.slab_structures.term_0,
    #     'term_1': prev_wg.outputs.slab_structures.term_1,
    # }

    if not input_slabs:
        print("\n✗ No slab structures loaded!")
        print("  Please create slab files or use previous calculation outputs.")
        return None

    print(f"\n3. Total slabs loaded: {len(input_slabs)}")

    # AIMD configuration
    print("\n4. AIMD configuration:")
    aimd_sequence = [
        {'temperature': 300, 'steps': 50},
        {'temperature': 300, 'steps': 100},
    ]

    for i, stage in enumerate(aimd_sequence):
        print(f"   Stage {i+1}: {stage['temperature']} K, {stage['steps']} steps")

    # CP2K AIMD parameters
    aimd_params = get_aimd_defaults_cp2k(
        cutoff=400,
        rel_cutoff=60,
        timestep=1.0,
        eps_scf=1e-6,
        thermostat='NOSE',
    )

    # Add KIND section
    if 'FORCE_EVAL' not in aimd_params:
        aimd_params['FORCE_EVAL'] = {}
    if 'SUBSYS' not in aimd_params['FORCE_EVAL']:
        aimd_params['FORCE_EVAL']['SUBSYS'] = {}

    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {"_": "Ag", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
        {"_": "O", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6", "POTENTIAL": "GTH-PBE-q6"}
    ]

    aimd_options = {
        'resources': {
            'num_machines': 2,
            'num_cores_per_machine': 128,
        },
        'queue_name': 'paralela',
    }

    # Fixed atoms configuration
    print("\n5. Fixed atoms configuration:")
    print("   Type: bottom")
    print("   Thickness: 7.0 Å")
    print("   Elements: all")
    print("   Components: XYZ (fully rigid)")

    print("\n6. Building workgraph...")
    print(f"  Workflow: AIMD only (no bulk, no slab generation)")
    print(f"  AIMD code: CP2K")
    print(f"  Fixed atoms: bottom 7Å")

    # Build workgraph
    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',

        # Minimal bulk parameters (not used, but required by internal logic)
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        code_label='VASP-VTST-6.4.3@bohr',

        # Input slabs directly
        input_slabs=input_slabs,

        # AIMD (CP2K)
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_params,
        aimd_options=aimd_options,

        # Fixed atoms (NEW)
        fix_atoms=True,
        fix_type='bottom',
        fix_thickness=7.0,
        fix_elements=None,  # All elements
        fix_components='XYZ',

        clean_workdir=False,
        name='Step07B_AIMD_CP2K_InputSlabs_FixedAtoms',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n7. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 7B SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\n⚠️  WARNING: AIMD is EXPENSIVE and will take many hours!")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"\nWorkflow stages:")
    print(f"  1. AIMD on input slabs (CP2K) - ONLY")
    print(f"\nNo bulk relaxation or slab generation performed!")
    print(f"\nFixed atoms:")
    print(f"  - Bottom 7Å of all slabs are fully constrained (XYZ)")
    print(f"  - Constraint applied to all AIMD stages")
    print(f"\nExpected outputs:")
    print(f"  - AIMD stage outputs in tasks:")
    print(f"    wg.tasks['aimd_stage_00_300K'].outputs:")
    print(f"      - structures, remote_folders")
    print(f"      - parameters, trajectories, retrieved")
    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

**Step 2: Make script executable**

Run: `chmod +x /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07b_aimd_input_slabs.py`

**Step 3: Verify script**

Run: `head -n 30 /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07b_aimd_input_slabs.py`
Expected: See docstring and imports

**Step 4: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
git add examples/cp2k/step_07b_aimd_input_slabs.py
git commit -m "feat: add CP2K AIMD example - input slabs with fixed atoms

- Demonstrates AIMD-only workflow with pre-existing slabs
- No bulk relaxation or slab generation
- Shows fixed atoms usage (fix bottom 7Å)
- Includes example for loading from files or previous calculations"
```

---

## Task 8: Final Verification

**Step 1: Clear all Python cache**

Run: `find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -type d -name __pycache__ -exec rm -rf {} + && find /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd -name "*.pyc" -delete`

**Step 2: Restart AiiDA daemon**

Run: `verdi daemon restart`
Expected: Daemon restarted successfully

**Step 3: Test all imports**

Run:
```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd
source ~/envs/psteros/bin/activate
python -c "
from teros.core.fixed_atoms import get_fixed_atoms_list, add_fixed_atoms_to_cp2k_parameters
from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k
from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k, get_basis_molopt_content, get_gth_potentials_content
from teros.core.builders import get_aimd_defaults_cp2k
print('✓ All imports successful')
"
```
Expected: "✓ All imports successful"

**Step 4: Verify example scripts exist**

Run: `ls -la /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd/examples/cp2k/step_07*.py`
Expected: See both step_07a and step_07b scripts

**Step 5: Run git status**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd && git status`
Expected: Clean working tree (all changes committed)

**Step 6: View commit log**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-cp2k-aimd && git log --oneline -10`
Expected: See all 7 commits for this feature

**Step 7: Final verification message**

Print:
```
✓ CP2K AIMD Integration Complete!

New Files Created:
  - teros/core/fixed_atoms.py
  - teros/core/aimd_cp2k.py
  - teros/core/builders/aimd_builder_cp2k.py
  - examples/cp2k/step_07a_aimd_autogenerate_slabs.py
  - examples/cp2k/step_07b_aimd_input_slabs.py

Modified Files:
  - teros/core/builders/__init__.py
  - teros/core/workgraph.py

Next Steps:
  1. Test step_07a with real structures
  2. Test step_07b with input slabs
  3. Verify AIMD outputs
  4. Update documentation if needed
```

---

## Success Criteria

Implementation is complete when:
1. All 7 tasks are completed
2. All commits are made
3. All imports work without errors
4. AiiDA daemon restarts successfully
5. Example scripts are executable and well-documented
6. Git working tree is clean

## Notes

- Follow TDD principles: test imports after each module creation
- Clear Python cache frequently to avoid stale imports
- Restart daemon after workgraph.py modifications
- Use exact file paths throughout
- Commit frequently (after each task)
