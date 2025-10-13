# AIMD Slabs Sequential Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Add ab initio molecular dynamics (AIMD) capability to PS-TEROS as a parallel, independent analysis on slab structures with sequential temperature stages.

**Architecture:** Create new `teros/core/aimd.py` module with scatter-gather pattern (similar to `slabs.py`, `thermodynamics.py`). Each slab gets sequential AIMD stages where each stage restarts from the previous one using VASP's restart_folder. AIMD runs in parallel to traditional relaxation as an optional analysis branch.

**Tech Stack:** AiiDA, AiiDA-WorkGraph, aiida-vasp, VASP (IBRION=0, MDALGO=2 for NVT), Python 3.10+

---

## Task 1: Create AIMD Core Module with Helper Function

**Files:**
- Create: `teros/core/aimd.py`

**Step 1: Create the file with imports and helper function**

Create `teros/core/aimd.py`:

```python
"""
AIMD Module for PS-TEROS

Ab initio molecular dynamics calculations on slab structures.
Sequential AIMD stages with automatic restart chaining.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace
from teros.core.slabs import extract_total_energy


def prepare_aimd_parameters(
    base_parameters: dict,
    temperature: float,
    steps: int
) -> dict:
    """
    Prepare INCAR parameters for a single AIMD stage.

    Takes base AIMD parameters and injects stage-specific values:
    - TEBEG and TEEND (both set to temperature for isothermal)
    - NSW (number of MD steps for this stage)

    Args:
        base_parameters: Base AIMD INCAR dict (IBRION=0, MDALGO, POTIM, etc.)
        temperature: Target temperature in K
        steps: Number of MD steps

    Returns:
        Complete INCAR dict for this stage
    """
    params = base_parameters.copy()
    params['TEBEG'] = temperature
    params['TEEND'] = temperature  # Isothermal
    params['NSW'] = steps
    return params
```

**Step 2: Verify imports work**

Run:
```bash
source ~/envs/psteros/bin/activate
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-aimd-slabs
python -c "from teros.core.aimd import prepare_aimd_parameters; print('Import successful')"
```

Expected: "Import successful"

**Step 3: Commit**

```bash
git add teros/core/aimd.py
git commit -m "feat(aimd): add core module with prepare_aimd_parameters helper"
```

---

## Task 2: Implement Sequential AIMD for Single Slab

**Files:**
- Modify: `teros/core/aimd.py`

**Step 1: Add get_settings helper reference**

At the top of `teros/core/aimd.py`, add to imports:

```python
def get_settings():
    """
    Parser settings for aiida-vasp (copied from workgraph.py for consistency).
    """
    return {
        'parser_settings': {
            'add_energy': True,
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }
```

**Step 2: Add aimd_sequential_slab function**

Add to `teros/core/aimd.py` after `prepare_aimd_parameters`:

```python
@task.graph(outputs=[dynamic(namespace())])
def aimd_sequential_slab(
    slab_label: str,
    structure: orm.StructureData,
    aimd_sequence: list,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
) -> t.Annotated[dict, dynamic(namespace())]:
    """
    Run sequential AIMD stages for a single slab termination.

    Creates a chain of VASP calculations where each stage:
    1. Uses structure from previous stage (or initial for first)
    2. Uses restart_folder from previous stage for WAVECAR/CHGCAR
    3. Runs MD at specified temperature for specified steps

    Args:
        slab_label: Label for this slab (e.g., 'term_0')
        structure: Initial slab structure
        aimd_sequence: List of dicts with 'temperature' and 'steps' keys
        code: VASP code
        aimd_parameters: Base AIMD INCAR parameters
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory

    Returns:
        Dictionary with flattened outputs for all stages:
            - stage_00_300K_structure, stage_00_300K_trajectory, etc.
            - final_structure, final_remote, final_trajectory
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    prev_structure = structure
    prev_remote = None
    all_outputs = {}

    # Loop through AIMD sequence
    for stage_idx, stage_config in enumerate(aimd_sequence):
        temp = stage_config['temperature']
        steps = stage_config['steps']

        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters(
            aimd_parameters, temp, steps
        )

        # Build VASP inputs
        vasp_inputs = {
            'structure': prev_structure,
            'code': code,
            'parameters': {'incar': stage_params},
            'options': options,
            'potential_family': potential_family,
            'potential_mapping': potential_mapping,
            'kpoints_spacing': kpoints_spacing,
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }

        # Add restart folder if available
        if prev_remote is not None:
            vasp_inputs['restart_folder'] = prev_remote

        # Create VASP task
        aimd_task = VaspTask(**vasp_inputs)

        # Store outputs with descriptive keys
        stage_prefix = f"stage_{stage_idx:02d}_{temp:03.0f}K"
        all_outputs[f"{stage_prefix}_structure"] = aimd_task.structure
        all_outputs[f"{stage_prefix}_trajectory"] = aimd_task.trajectory
        all_outputs[f"{stage_prefix}_energy"] = extract_total_energy(
            energies=aimd_task.misc
        ).result
        all_outputs[f"{stage_prefix}_remote"] = aimd_task.remote_folder
        all_outputs[f"{stage_prefix}_retrieved"] = aimd_task.retrieved

        # Update for next stage
        prev_structure = aimd_task.structure
        prev_remote = aimd_task.remote_folder

    # Add final outputs
    all_outputs['final_structure'] = prev_structure
    all_outputs['final_remote'] = prev_remote
    all_outputs['final_trajectory'] = aimd_task.trajectory

    return all_outputs
```

**Step 3: Verify syntax**

Run:
```bash
python -c "from teros.core.aimd import aimd_sequential_slab; print('Import successful')"
```

Expected: "Import successful"

**Step 4: Commit**

```bash
git add teros/core/aimd.py
git commit -m "feat(aimd): add aimd_sequential_slab for single slab AIMD chain"
```

---

## Task 3: Implement Scatter Function for All Slabs

**Files:**
- Modify: `teros/core/aimd.py`

**Step 1: Add aimd_slabs_scatter function**

Add to `teros/core/aimd.py` after `aimd_sequential_slab`:

```python
def aimd_slabs_scatter(
    slabs: dict,
    aimd_sequence: list,
    code: orm.Code,
    aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
) -> t.Annotated[dict, dynamic(namespace())]:
    """
    Run AIMD on all slabs in parallel using scatter-gather pattern.

    Each slab gets the same AIMD sequence but runs independently in parallel.

    Args:
        slabs: Dictionary of slab structures (from slab generation)
        aimd_sequence: List of {'temperature': K, 'steps': N} dicts
        code: VASP code
        aimd_parameters: Base AIMD INCAR parameters
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory

    Returns:
        Dynamic namespace with nested structure:
            - term_0: {stage_00_300K_structure, ..., final_structure}
            - term_1: {stage_00_300K_structure, ..., final_structure}
            - ...
    """
    aimd_outputs = {}

    for slab_label, slab_structure in slabs.items():
        # Run sequential AIMD for this slab
        slab_aimd_result = aimd_sequential_slab(
            slab_label=slab_label,
            structure=slab_structure,
            aimd_sequence=aimd_sequence,
            code=code,
            aimd_parameters=aimd_parameters,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            kpoints_spacing=kpoints_spacing,
            clean_workdir=clean_workdir,
        )

        # Store outputs for this slab
        aimd_outputs[slab_label] = slab_aimd_result

    return aimd_outputs
```

**Step 2: Verify syntax**

Run:
```bash
python -c "from teros.core.aimd import aimd_slabs_scatter; print('Import successful')"
```

Expected: "Import successful"

**Step 3: Commit**

```bash
git add teros/core/aimd.py
git commit -m "feat(aimd): add aimd_slabs_scatter for parallel AIMD on all slabs"
```

---

## Task 4: Create AIMD Builder for Default Parameters

**Files:**
- Create: `teros/core/builders/aimd_builder.py`
- Modify: `teros/core/builders/__init__.py`

**Step 1: Create aimd_builder.py**

Create `teros/core/builders/aimd_builder.py`:

```python
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
```

**Step 2: Update builders __init__.py**

Modify `teros/core/builders/__init__.py` to add:

```python
from teros.core.builders.aimd_builder import get_aimd_defaults
```

**Step 3: Verify imports**

Run:
```bash
python -c "from teros.core.builders import get_aimd_defaults; print(get_aimd_defaults())"
```

Expected: Dictionary with AIMD parameters printed

**Step 4: Commit**

```bash
git add teros/core/builders/aimd_builder.py teros/core/builders/__init__.py
git commit -m "feat(aimd): add aimd_builder with get_aimd_defaults"
```

---

## Task 5: Integrate AIMD into core_workgraph Function

**Files:**
- Modify: `teros/core/workgraph.py`

**Step 1: Add aimd import at top of file**

In `teros/core/workgraph.py`, add to imports (around line 13-27):

```python
from teros.core.aimd import aimd_slabs_scatter
```

**Step 2: Add aimd_results to outputs decorator**

In `@task.graph(outputs=[...])` decorator for `core_workgraph` (around line 59-68), add:

```python
'aimd_results',
```

**Step 3: Add AIMD parameters to function signature**

In `core_workgraph()` function signature (around line 69-112), add these parameters after `compute_cleavage`:

```python
    run_aimd: bool = False,
    aimd_sequence: list = None,
    aimd_parameters: dict = None,
    aimd_options: dict = None,
    aimd_potential_mapping: dict = None,
    aimd_kpoints_spacing: float = None,
```

**Step 4: Add AIMD section in function body**

After the cleavage calculation section (around line 414), add:

```python
    # ===== AIMD CALCULATION (OPTIONAL) =====
    aimd_outputs = {}
    if run_aimd and slab_namespace is not None:
        # Use AIMD-specific parameters or fall back to slab parameters
        aimd_params = aimd_parameters if aimd_parameters is not None else slab_params
        aimd_opts = aimd_options if aimd_options is not None else slab_opts
        aimd_pot_map = aimd_potential_mapping if aimd_potential_mapping is not None else slab_pot_map
        aimd_kpts = aimd_kpoints_spacing if aimd_kpoints_spacing is not None else slab_kpts

        # Run AIMD on all slabs in parallel
        aimd_results = aimd_slabs_scatter(
            slabs=slab_namespace,
            aimd_sequence=aimd_sequence,
            code=code,
            aimd_parameters=aimd_params,
            potential_family=potential_family,
            potential_mapping=aimd_pot_map,
            options=aimd_opts,
            kpoints_spacing=aimd_kpts,
            clean_workdir=clean_workdir,
        )

        aimd_outputs = aimd_results
```

**Step 5: Add aimd_results to return statement**

In the return statement at the end of `core_workgraph()` (around line 416-438), add:

```python
        'aimd_results': aimd_outputs,
```

**Step 6: Verify syntax**

Run:
```bash
python -c "from teros.core.workgraph import core_workgraph; print('Import successful')"
```

Expected: "Import successful"

**Step 7: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat(aimd): integrate AIMD into core_workgraph as optional branch"
```

---

## Task 6: Add AIMD Parameters to build_core_workgraph

**Files:**
- Modify: `teros/core/workgraph.py`

**Step 1: Add AIMD parameters to function signature**

In `build_core_workgraph()` function signature (around line 441-488), add after `compute_electronic_properties_bulk` parameters:

```python
    run_aimd: bool = False,
    aimd_sequence: list = None,
    aimd_parameters: dict = None,
    aimd_options: dict = None,
    aimd_potential_mapping: dict = None,
    aimd_kpoints_spacing: float = None,
```

**Step 2: Update docstring**

In `build_core_workgraph()` docstring (around line 489-579), add after electronic properties description:

```python
        run_aimd: Run AIMD on generated/input slabs. Default: False
        aimd_sequence: List of AIMD stages [{'temperature': K, 'steps': N}, ...]. Default: None
        aimd_parameters: AIMD INCAR parameters (use get_aimd_defaults()). Default: None
        aimd_options: Scheduler options for AIMD. Default: None (uses slab_options)
        aimd_potential_mapping: Potential mapping for AIMD. Default: None (uses slab_potential_mapping)
        aimd_kpoints_spacing: K-points spacing for AIMD. Default: None (uses slab_kpoints_spacing)
```

**Step 3: Pass AIMD parameters to core_workgraph.build()**

In the `core_workgraph.build()` call (around line 635-677), add after `compute_cleavage`:

```python
        run_aimd=run_aimd,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        aimd_options=aimd_options,
        aimd_potential_mapping=aimd_potential_mapping,
        aimd_kpoints_spacing=aimd_kpoints_spacing,
```

**Step 4: Verify syntax**

Run:
```bash
python -c "from teros.core.workgraph import build_core_workgraph; print('Import successful')"
```

Expected: "Import successful"

**Step 5: Clear Python cache**

Run:
```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

**Step 6: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat(aimd): add AIMD parameters to build_core_workgraph"
```

---

## Task 7: Create Complete AIMD Example for Ag2O

**Files:**
- Create: `examples/complete/complete_ag2o_aimd_example.py`

**Step 1: Create example file**

Copy `examples/complete/complete_ag2o_example.py` to `examples/complete/complete_ag2o_aimd_example.py` and modify:

Add after line 37 (after imports):

```python
from teros.core.builders import get_aimd_defaults
```

Add after line 313 (after electronic properties section):

```python
    # ===== AIMD PARAMETERS =====
    print("\n" + "=" * 80)
    print("AIMD PARAMETERS")
    print("=" * 80)

    # Get default AIMD parameters
    aimd_parameters = get_aimd_defaults(
        energy_cutoff=400,
        timestep=1.0,  # 1 fs
        mdalgo=2,  # NVT Nosé-Hoover
        ismear=0,
        sigma=0.1,
    )

    # Define AIMD sequence - shorter for testing
    aimd_sequence = [
        {'temperature': 300, 'steps': 100},  # Short test run
        {'temperature': 300, 'steps': 100},  # Another short run
    ]

    print(f"  AIMD sequence: {len(aimd_sequence)} stages")
    print(f"  Timestep: {aimd_parameters['POTIM']} fs")
    print(f"  Thermostat: Nosé-Hoover (MDALGO={aimd_parameters['MDALGO']})")
    for i, stage in enumerate(aimd_sequence):
        print(f"    Stage {i}: {stage['temperature']}K × {stage['steps']} steps")
```

Update `build_core_workgraph()` call (around line 320) to add:

```python
        # AIMD parameters (NEW)
        run_aimd=True,
        aimd_sequence=aimd_sequence,
        aimd_parameters=aimd_parameters,
        aimd_options=slab_options,
```

Update workflow name:

```python
        name='Ag2O_Complete_Workflow_with_AIMD',
```

Update expected outputs section (around line 460):

```python
    print("\n8. AIMD Results:")
    print("   - aimd_results: Nested namespace with all AIMD outputs")
    print("   - For each slab (term_0, term_1, ...):")
    print("     - stage_00_300K_structure, _trajectory, _energy, _remote, _retrieved")
    print("     - stage_01_300K_structure, _trajectory, _energy, _remote, _retrieved")
    print("     - final_structure, final_remote, final_trajectory")
```

**Step 2: Verify syntax**

Run:
```bash
python -c "import sys; sys.path.insert(0, '/home/thiagotd/git/PS-TEROS/.worktree/feature-aimd-slabs'); exec(open('examples/complete/complete_ag2o_aimd_example.py').read().replace('wg.submit', '# wg.submit').replace('if __name__', 'if False'))"
```

Expected: Script runs without syntax errors (won't submit)

**Step 3: Commit**

```bash
git add examples/complete/complete_ag2o_aimd_example.py
git commit -m "feat(aimd): add complete Ag2O example with AIMD"
```

---

## Task 8: Create Complete AIMD Example for Ag3PO4

**Files:**
- Create: `examples/complete/complete_ag3po4_aimd_example.py`

**Step 1: Create example file**

Copy `examples/complete/complete_ag3po4_example.py` to `examples/complete/complete_ag3po4_aimd_example.py` and apply same modifications as Task 7:

- Add `get_aimd_defaults` import
- Add AIMD parameters section
- Add AIMD sequence (2 short stages for testing)
- Update `build_core_workgraph()` call with AIMD parameters
- Update workflow name to include "_with_AIMD"
- Update expected outputs section

**Step 2: Verify syntax**

Run:
```bash
python -c "import sys; sys.path.insert(0, '/home/thiagotd/git/PS-TEROS/.worktree/feature-aimd-slabs'); exec(open('examples/complete/complete_ag3po4_aimd_example.py').read().replace('wg.submit', '# wg.submit').replace('if __name__', 'if False'))"
```

Expected: Script runs without syntax errors

**Step 3: Commit**

```bash
git add examples/complete/complete_ag3po4_aimd_example.py
git commit -m "feat(aimd): add complete Ag3PO4 example with AIMD"
```

---

## Task 9: Restart Daemon and Test Binary Case

**Files:**
- Test: `examples/complete/complete_ag2o_aimd_example.py`

**Step 1: Clear cache and restart daemon**

Run:
```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
verdi daemon restart
sleep 5
```

**Step 2: Submit Ag2O AIMD workflow**

Run:
```bash
source ~/envs/psteros/bin/activate && /home/thiagotd/envs/psteros/bin/python /home/thiagotd/git/PS-TEROS/.worktree/feature-aimd-slabs/examples/complete/complete_ag2o_aimd_example.py
```

Expected: Workflow submits successfully, prints PK

**Step 3: Wait and check status**

Run:
```bash
sleep 30
verdi process list
```

Expected: See running AIMD processes

**Step 4: After completion, check outputs**

Run (replace PK):
```bash
verdi process show <PK>
```

Expected: See `aimd_results` in outputs with nested structure

**Step 5: Verify AIMD outputs exist**

Run (replace PK):
```bash
verdi process show <PK> | grep -A 20 aimd_results
```

Expected: See term_0, term_1 with stage outputs

---

## Task 10: Test Ternary Case

**Files:**
- Test: `examples/complete/complete_ag3po4_aimd_example.py`

**Step 1: Submit Ag3PO4 AIMD workflow**

Run:
```bash
source ~/envs/psteros/bin/activate && /home/thiagotd/envs/psteros/bin/python /home/thiagotd/git/PS-TEROS/.worktree/feature-aimd-slabs/examples/complete/complete_ag3po4_aimd_example.py
```

Expected: Workflow submits successfully

**Step 2: Wait and check status**

Run:
```bash
sleep 30
verdi process list
```

Expected: See running AIMD processes for ternary case

**Step 3: After completion, verify outputs**

Run (replace PK):
```bash
verdi process show <PK>
verdi process report <PK>
```

Expected: Workflow completes with [0], aimd_results present

**Step 4: Check AIMD chaining worked**

Verify that stage 1 used restart_folder from stage 0 by checking task connections in HTML visualization or by inspecting the workgraph structure.

---

## Task 11: Create Documentation

**Files:**
- Create: `examples/complete/AIMD_README.md`

**Step 1: Create documentation**

Create `examples/complete/AIMD_README.md`:

```markdown
# AIMD (Ab Initio Molecular Dynamics) in PS-TEROS

## Overview

AIMD capability allows running sequential molecular dynamics simulations on slab structures as a parallel, independent analysis alongside traditional relaxation calculations.

## Features

- Sequential AIMD stages with automatic restart chaining
- Flexible temperature/timestep configuration
- Runs on all generated or input slabs
- Full provenance tracking of all stages
- Compatible with binary and ternary oxides

## Usage

```python
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_aimd_defaults

# Define AIMD parameters
aimd_parameters = get_aimd_defaults(
    energy_cutoff=400,
    timestep=1.0,  # 1 fs
)

# Define AIMD sequence
aimd_sequence = [
    {'temperature': 300, 'steps': 1000},
    {'temperature': 300, 'steps': 1000},
    {'temperature': 500, 'steps': 1000},
]

# Build workflow
wg = build_core_workgraph(
    # ... standard parameters ...
    run_aimd=True,
    aimd_sequence=aimd_sequence,
    aimd_parameters=aimd_parameters,
    aimd_options=slab_options,
)
```

## Output Structure

```
aimd_results/
  term_0/
    stage_00_300K_structure
    stage_00_300K_trajectory
    stage_00_300K_energy
    stage_00_300K_remote
    stage_00_300K_retrieved
    stage_01_300K_structure
    ...
    final_structure
    final_remote
    final_trajectory
  term_1/
    ...
```

## Examples

- Binary oxide: `examples/complete/complete_ag2o_aimd_example.py`
- Ternary oxide: `examples/complete/complete_ag3po4_aimd_example.py`

## VASP Parameters

AIMD uses NVT ensemble (Nosé-Hoover thermostat):
- `IBRION = 0` (MD mode)
- `MDALGO = 2` (NVT)
- `POTIM` (timestep in fs)
- `TEBEG/TEEND` (set automatically per stage)
- `NSW` (set automatically per stage)

## Notes

- Each stage restarts from previous using WAVECAR/CHGCAR
- Fixed atoms support coming in future update
- External restart capability coming in future update
```

**Step 2: Commit**

```bash
git add examples/complete/AIMD_README.md
git commit -m "docs(aimd): add AIMD documentation"
```

---

## Testing Checklist

After implementation:

- [ ] Binary oxide (Ag2O) AIMD completes successfully
- [ ] Ternary oxide (Ag3PO4) AIMD completes successfully
- [ ] AIMD runs in parallel with relaxation (both can run simultaneously)
- [ ] Each AIMD stage uses restart_folder from previous stage
- [ ] Output structure matches design (nested namespace with all stages)
- [ ] Works with auto-generated slabs
- [ ] Works with manually input slabs
- [ ] All outputs accessible via `verdi process show`
- [ ] HTML visualization shows AIMD task chain
- [ ] No errors in `verdi process report`

---

## Future Enhancements

Not in this implementation:
- Fixed atoms support (selective dynamics for slabs)
- External restart capability (restart AIMD from previous workflow PK)
- Stage-level restart (restart from specific stage)
- AIMD-specific analysis tools (MSD, RDF, etc.)
