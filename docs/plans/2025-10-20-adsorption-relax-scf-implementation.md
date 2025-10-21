# Adsorption Energy Relaxation + SCF Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add optional relaxation step before adsorption energy calculations with full VASP builder input support

**Architecture:** Sequential 3-phase workflow: (1) optional relax complete system using vasp.v2.relax, (2) separate relaxed structure into substrate/molecule/complete, (3) SCF calculations using vasp.v2.vasp. Maintains backward compatibility with existing parameter dict API.

**Tech Stack:** AiiDA-WorkGraph, aiida-vasp (v2.relax, v2.vasp), pymatgen, NetworkX

---

## Task 1: Add `_build_vasp_inputs()` Helper Function

**Files:**
- Modify: `teros/core/adsorption_energy.py` (add after line 23, before `parse_formula()`)

**Context:** This helper constructs VASP WorkChain inputs from either new-style builder_inputs or old-style parameters dict, with automatic SCF enforcement.

**Step 1: Write failing test for new-style builder inputs**

Create test file `teros/core/tests/test_adsorption_energy_builder.py`:

```python
"""Tests for VASP builder input construction in adsorption_energy module."""

import pytest
from aiida import orm
from teros.core.adsorption_energy import _build_vasp_inputs


def test_build_vasp_inputs_from_builder_dict():
    """Test building inputs from new-style builder_inputs dict."""
    # Create mock structure and code
    from ase import Atoms
    ase_struct = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure = orm.StructureData(ase=ase_struct)
    code = orm.Code.get_from_string('VASP-VTST-6.4.3@bohr')

    # Builder inputs
    builder_inputs = {
        'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 520}},
        'options': {'resources': {'num_machines': 1}},
        'potential_family': 'PBE',
        'potential_mapping': {'H': 'H'},
        'kpoints_spacing': 0.3,
    }

    # Build inputs
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=builder_inputs,
    )

    # Verify structure and code are set
    assert result['structure'] == structure
    assert result['code'] == code

    # Verify builder_inputs are preserved
    assert result['parameters']['incar']['PREC'] == 'Accurate'
    assert result['parameters']['incar']['ENCUT'] == 520
    assert result['options']['resources']['num_machines'] == 1
    assert result['potential_family'] == 'PBE'
    assert result['potential_mapping']['H'] == 'H'
    assert result['kpoints_spacing'] == 0.3
```

**Step 2: Run test to verify it fails**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy
source ~/envs/aiida/bin/activate
pytest teros/core/tests/test_adsorption_energy_builder.py::test_build_vasp_inputs_from_builder_dict -v
```

Expected: `ImportError: cannot import name '_build_vasp_inputs'` or `AttributeError`

**Step 3: Write failing test for old-style parameters**

Add to `teros/core/tests/test_adsorption_energy_builder.py`:

```python
def test_build_vasp_inputs_from_parameters_dict():
    """Test building inputs from old-style parameters dict (backward compat)."""
    from ase import Atoms
    ase_struct = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure = orm.StructureData(ase=ase_struct)
    code = orm.Code.get_from_string('VASP-VTST-6.4.3@bohr')

    # Old-style parameters
    parameters = {'PREC': 'Accurate', 'ENCUT': 520, 'NSW': 100}
    options = {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40}}
    potential_mapping = {'H': 'H'}

    # Build inputs
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=None,
        parameters=parameters,
        options=options,
        potential_family='PBE',
        potential_mapping=potential_mapping,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    # Verify structure and code are set
    assert result['structure'] == structure
    assert result['code'] == code

    # Verify parameters are wrapped in 'incar' dict
    assert result['parameters']['incar']['PREC'] == 'Accurate'
    assert result['parameters']['incar']['ENCUT'] == 520
    assert result['parameters']['incar']['NSW'] == 100

    # Verify options
    assert result['options']['resources']['num_machines'] == 1
    assert result['options']['resources']['num_mpiprocs_per_machine'] == 40

    # Verify other fields
    assert result['potential_family'] == 'PBE'
    assert result['potential_mapping']['H'] == 'H'
    assert result['kpoints_spacing'] == 0.3
    assert result['clean_workdir'] == True

    # Verify settings dict is created
    assert 'settings' in result
    assert isinstance(result['settings'], orm.Dict)
```

**Step 4: Write failing test for force_scf enforcement**

Add to `teros/core/tests/test_adsorption_energy_builder.py`:

```python
def test_build_vasp_inputs_force_scf():
    """Test that force_scf=True enforces NSW=0 and IBRION=-1."""
    from ase import Atoms
    ase_struct = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure = orm.StructureData(ase=ase_struct)
    code = orm.Code.get_from_string('VASP-VTST-6.4.3@bohr')

    # Builder inputs with relaxation settings
    builder_inputs = {
        'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ENCUT': 520}},
        'options': {'resources': {'num_machines': 1}},
        'potential_family': 'PBE',
        'potential_mapping': {'H': 'H'},
    }

    # Build inputs with force_scf=True
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=builder_inputs,
        force_scf=True,
    )

    # Verify NSW=0 and IBRION=-1 are enforced
    assert result['parameters']['incar']['NSW'] == 0
    assert result['parameters']['incar']['IBRION'] == -1

    # Verify other parameters are preserved
    assert result['parameters']['incar']['ENCUT'] == 520
```

**Step 5: Run tests to verify they fail**

```bash
pytest teros/core/tests/test_adsorption_energy_builder.py -v
```

Expected: All 3 tests FAIL (function not defined yet)

**Step 6: Implement `_build_vasp_inputs()` helper function**

Add to `teros/core/adsorption_energy.py` after line 23:

```python
def _build_vasp_inputs(
    structure: orm.StructureData,
    code: orm.Code,
    builder_inputs: dict | None = None,
    parameters: dict | None = None,
    options: dict | None = None,
    potential_family: str | None = None,
    potential_mapping: dict | None = None,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    force_scf: bool = False,
) -> dict:
    """
    Build VASP WorkChain inputs from either builder_inputs or old-style parameters.

    This helper provides backward compatibility by accepting either:
    1. New-style: builder_inputs dict (full control, recommended)
    2. Old-style: parameters/options/etc dicts (automatic conversion)

    Args:
        structure: Structure to calculate
        code: AiiDA code for VASP
        builder_inputs: New-style builder dict (takes priority if provided)
        parameters: Old-style INCAR dict (fallback)
        options: Old-style scheduler dict (fallback)
        potential_family: Pseudopotential family name
        potential_mapping: Element → potential mapping
        kpoints_spacing: K-points spacing in Angstrom^-1
        clean_workdir: Whether to clean remote working directory
        force_scf: If True, enforce NSW=0 and IBRION=-1 for single-point calculation

    Returns:
        Complete input dictionary for VASP WorkChain

    Raises:
        ValueError: If neither builder_inputs nor parameters is provided
    """
    if builder_inputs is not None:
        # Use new-style builder inputs
        # Make a deep copy to avoid modifying original
        import copy
        inputs = copy.deepcopy(builder_inputs)

        # Override structure and code (always set by workflow)
        inputs['structure'] = structure
        inputs['code'] = code

    elif parameters is not None:
        # Construct from old-style parameters
        inputs = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping) if potential_mapping is not None else {},
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict={
                'parser_settings': {
                    'add_trajectory': True,
                    'add_structure': True,
                    'add_kpoints': True,
                }
            }),
        }

        if kpoints_spacing is not None:
            inputs['kpoints_spacing'] = kpoints_spacing

    else:
        raise ValueError(
            "Must provide either 'builder_inputs' or 'parameters'. "
            "Use builder_inputs for full control (recommended) or parameters for backward compatibility."
        )

    # Force SCF mode if requested
    if force_scf:
        # Ensure parameters dict exists and has 'incar' key
        if 'parameters' not in inputs:
            inputs['parameters'] = {'incar': {}}
        if 'incar' not in inputs['parameters']:
            inputs['parameters']['incar'] = {}

        # Override NSW and IBRION for single-point calculation
        inputs['parameters']['incar']['NSW'] = 0
        inputs['parameters']['incar']['IBRION'] = -1

    return inputs
```

**Step 7: Run tests to verify they pass**

```bash
pytest teros/core/tests/test_adsorption_energy_builder.py -v
```

Expected: All 3 tests PASS

**Step 8: Commit**

```bash
git add teros/core/adsorption_energy.py teros/core/tests/test_adsorption_energy_builder.py
git commit -m "feat: add _build_vasp_inputs helper for flexible VASP input construction"
```

---

## Task 2: Update `compute_adsorption_energies_scatter()` Signature

**Files:**
- Modify: `teros/core/adsorption_energy.py:280-297` (function signature and docstring)

**Context:** Add new parameters for relaxation control and SCF builder inputs while maintaining backward compatibility.

**Step 1: Update function signature**

Modify `teros/core/adsorption_energy.py` line 280:

```python
@task.graph
def compute_adsorption_energies_scatter(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    adsorbate_formulas: t.Annotated[dict[str, str], dict],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],

    # NEW: Relaxation control
    relax_before_adsorption: bool = False,
    relax_builder_inputs: dict | None = None,

    # NEW: SCF control
    scf_builder_inputs: dict | None = None,

    # DEPRECATED (kept for backward compatibility)
    parameters: t.Mapping[str, t.Any] | None = None,
    options: t.Mapping[str, t.Any] | None = None,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    relaxed_complete_structures=dynamic(orm.StructureData),  # NEW
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
)]:
```

**Step 2: Update docstring**

Replace docstring at line 299-325:

```python
    """
    Scatter-gather workflow for calculating adsorption energies.

    This workflow supports two modes:
    1. Direct SCF: separate → SCF (3N jobs)
    2. Relax + SCF: relax → separate → SCF (N+3N jobs)

    Workflow phases:
    - Phase 1 (optional): Relax complete structures using vasp.v2.relax
    - Phase 2: Separate relaxed (or original) structures into substrate/molecule/complete
    - Phase 3: SCF calculations using vasp.v2.vasp for all three components

    Args:
        structures: Dynamic namespace of complete structures (substrate + adsorbate)
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'system1': 'OOH', 'system2': 'OH'}
        code: AiiDA code for VASP
        potential_family: Pseudopotential family name
        potential_mapping: Element to potential mapping

        relax_before_adsorption: If True, relax complete structures before separation
        relax_builder_inputs: Full builder dict for vasp.v2.relax (NSW, IBRION, ISIF)
        scf_builder_inputs: Full builder dict for vasp.v2.vasp (NSW=0 enforced)

        parameters: [DEPRECATED] Old-style INCAR parameters dict (for backward compat)
        options: [DEPRECATED] Old-style scheduler options dict
        kpoints_spacing: [DEPRECATED] K-points spacing in Angstrom^-1
        clean_workdir: Whether to clean remote working directories

    Returns:
        Dictionary with namespaces:
        - relaxed_complete_structures: Relaxed structures (empty if relax=False)
        - separated_structures: Dict of separated systems for each input
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies of complete systems
        - adsorption_energies: Final adsorption energies (E_ads = E_complete - E_substrate - E_molecule)

    Example (new API):
        >>> relax_inputs = {
        ...     'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ENCUT': 520}},
        ...     'options': {'resources': {'num_machines': 1}},
        ...     'potential_family': 'PBE',
        ...     'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
        ... }
        >>> scf_inputs = {
        ...     'parameters': {'incar': {'ENCUT': 520}},  # NSW=0 added automatically
        ...     'options': {'resources': {'num_machines': 1}},
        ...     'potential_family': 'PBE',
        ...     'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
        ... }
        >>> results = compute_adsorption_energies_scatter(
        ...     structures={'oh_lamno3': structure},
        ...     adsorbate_formulas={'oh_lamno3': 'OH'},
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={...},
        ...     relax_before_adsorption=True,
        ...     relax_builder_inputs=relax_inputs,
        ...     scf_builder_inputs=scf_inputs,
        ... )

    Example (old API, backward compatible):
        >>> results = compute_adsorption_energies_scatter(
        ...     structures={'oh_ag': structure},
        ...     adsorbate_formulas={'oh_ag': 'OH'},
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={'Ag': 'Ag', 'O': 'O', 'H': 'H'},
        ...     parameters={'PREC': 'Accurate', 'ENCUT': 520},
        ...     options={'resources': {'num_machines': 1}},
        ...     kpoints_spacing=0.3,
        ... )
    """
```

**Step 3: Verify syntax**

```bash
python3 -m py_compile teros/core/adsorption_energy.py
```

Expected: No syntax errors

**Step 4: Commit**

```bash
git add teros/core/adsorption_energy.py
git commit -m "feat: update compute_adsorption_energies_scatter signature with relax/scf params"
```

---

## Task 3: Implement Phase 1 (Optional Relaxation)

**Files:**
- Modify: `teros/core/adsorption_energy.py:360-361` (after validation, before Phase 2)

**Context:** Add relaxation logic that runs before structure separation when `relax_before_adsorption=True`.

**Step 1: Add Phase 1 implementation**

Insert after line 360 (after input validation, before Phase 2 separation loop):

```python
    # Get workflow plugins
    vasp_relax = WorkflowFactory('vasp.v2.relax')
    vasp_scf = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_relax)
    scf_task_cls = task(vasp_scf)

    # ===== PHASE 1: OPTIONAL RELAXATION =====
    # Relax complete structures before separation if requested
    if relax_before_adsorption:
        relaxed_structures: dict[str, orm.StructureData] = {}

        for key, structure in structures.items():
            # Build relaxation inputs
            relax_inputs = _build_vasp_inputs(
                structure=structure,
                code=code,
                builder_inputs=relax_builder_inputs,
                parameters=parameters,
                options=options,
                potential_family=potential_family,
                potential_mapping=potential_mapping,
                kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
                force_scf=False,  # Allow relaxation (NSW > 0)
            )

            # Run relaxation
            relax_task = relax_task_cls(**relax_inputs)
            relaxed_structures[key] = relax_task.outputs.structure

        # Use relaxed structures for separation
        structures_to_separate = relaxed_structures
    else:
        # No relaxation: use original structures
        relaxed_structures = {}
        structures_to_separate = structures
```

**Step 2: Verify syntax**

```bash
python3 -m py_compile teros/core/adsorption_energy.py
```

Expected: No syntax errors

**Step 3: Commit**

```bash
git add teros/core/adsorption_energy.py
git commit -m "feat: implement Phase 1 optional relaxation for complete structures"
```

---

## Task 4: Update Phase 2 (Separation) to Use Conditional Input

**Files:**
- Modify: `teros/core/adsorption_energy.py:364-377` (separation loop)

**Context:** Update separation phase to use `structures_to_separate` instead of `structures` directly.

**Step 1: Update separation loop**

Replace the separation loop (currently at lines 364-377):

```python
    # ===== PHASE 2: STRUCTURE SEPARATION =====
    # Separate relaxed (if Phase 1 ran) or original structures
    separated_dict: dict[str, dict] = {}
    for key, structure in structures_to_separate.items():
        adsorbate_str = adsorbate_formulas[key]

        separated = separate_adsorbate_structure(
            structure=structure,
            adsorbate_formula=orm.Str(adsorbate_str)
        )

        separated_dict[key] = {
            'substrate': separated.substrate,
            'molecule': separated.molecule,
            'complete': separated.complete,
        }
```

**Step 2: Verify syntax**

```bash
python3 -m py_compile teros/core/adsorption_energy.py
```

Expected: No syntax errors

**Step 3: Commit**

```bash
git add teros/core/adsorption_energy.py
git commit -m "refactor: use structures_to_separate in Phase 2 separation"
```

---

## Task 5: Update Phase 3 (SCF) to Use New Builder Inputs

**Files:**
- Modify: `teros/core/adsorption_energy.py:379-413` (VASP calculations loop)

**Context:** Replace current VASP task creation with new builder-based approach using `vasp.v2.vasp` plugin explicitly.

**Step 1: Update SCF calculations section**

Replace the VASP calculations section (lines 379-413):

```python
    # ===== PHASE 3: SCF CALCULATIONS =====
    # Single-point calculations for substrate, molecule, and complete systems
    substrate_energies: dict[str, orm.Float] = {}
    molecule_energies: dict[str, orm.Float] = {}
    complete_energies: dict[str, orm.Float] = {}

    for key, separated in separated_dict.items():
        # Helper function to create SCF inputs for each structure
        def create_scf_inputs(struct: orm.StructureData) -> dict:
            return _build_vasp_inputs(
                structure=struct,
                code=code,
                builder_inputs=scf_builder_inputs,
                parameters=parameters,
                options=options,
                potential_family=potential_family,
                potential_mapping=potential_mapping,
                kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
                force_scf=True,  # Enforce NSW=0, IBRION=-1
            )

        # Substrate SCF
        substrate_calc = scf_task_cls(**create_scf_inputs(separated['substrate']))
        substrate_energies[key] = extract_total_energy(energies=substrate_calc.misc).result

        # Molecule SCF
        molecule_calc = scf_task_cls(**create_scf_inputs(separated['molecule']))
        molecule_energies[key] = extract_total_energy(energies=molecule_calc.misc).result

        # Complete system SCF
        complete_calc = scf_task_cls(**create_scf_inputs(separated['complete']))
        complete_energies[key] = extract_total_energy(energies=complete_calc.misc).result
```

**Step 2: Remove old VASP workflow setup**

Delete lines 335-336 (old vasp_wc and vasp_task_cls setup):

```python
    # Get VASP workflow
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    vasp_task_cls = task(vasp_wc)
```

These are now defined earlier in Phase 1 implementation.

**Step 3: Verify syntax**

```bash
python3 -m py_compile teros/core/adsorption_energy.py
```

Expected: No syntax errors

**Step 4: Commit**

```bash
git add teros/core/adsorption_energy.py
git commit -m "feat: update Phase 3 SCF to use new builder inputs with force_scf"
```

---

## Task 6: Update Return Statement with Relaxed Structures

**Files:**
- Modify: `teros/core/adsorption_energy.py:396-413` (Phase 4 and return statement)

**Context:** Add `relaxed_complete_structures` to the output namespace.

**Step 1: Update Phase 4 comment**

Change the comment before adsorption energy calculation (around line 396):

```python
    # ===== PHASE 4: ADSORPTION ENERGY CALCULATION =====
    # Calculate E_ads = E_complete - E_substrate - E_molecule
    adsorption_energies: dict[str, orm.Float] = {}
    for key in structures.keys():
        E_ads = calculate_adsorption_energy(
            E_complete=complete_energies[key],
            E_substrate=substrate_energies[key],
            E_molecule=molecule_energies[key],
        ).result

        adsorption_energies[key] = E_ads
```

**Step 2: Update return statement**

Replace the return statement (around line 407):

```python
    # Return all results
    return {
        'relaxed_complete_structures': relaxed_structures,  # NEW: Empty dict if relax=False
        'separated_structures': separated_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': complete_energies,
        'adsorption_energies': adsorption_energies,
    }
```

**Step 3: Verify syntax**

```bash
python3 -m py_compile teros/core/adsorption_energy.py
```

Expected: No syntax errors

**Step 4: Commit**

```bash
git add teros/core/adsorption_energy.py
git commit -m "feat: add relaxed_complete_structures to output namespace"
```

---

## Task 7: Integrate into `build_core_workgraph()`

**Files:**
- Modify: `teros/core/workgraph.py:149-160` (add parameters to function signature)
- Modify: `teros/core/workgraph.py:1556-1603` (update adsorption task creation)

**Context:** Add new parameters to the main workgraph builder and connect the new outputs.

**Step 1: Add parameters to function signature**

Find the `build_core_workgraph()` function signature (around line 149) and add after `adsorption_kpoints_spacing`:

```python
    # Adsorption energy: Relaxation control (NEW)
    relax_before_adsorption: bool = False,
    adsorption_relax_builder_inputs: dict = None,
    adsorption_scf_builder_inputs: dict = None,
```

**Step 2: Update adsorption task creation**

Find the adsorption energy section (around line 1556-1603) and replace with:

```python
    # ===== ADSORPTION ENERGY CALCULATION (OPTIONAL) =====
    # Check if adsorption energy should be computed
    should_add_adsorption = run_adsorption_energy and adsorption_structures is not None and adsorption_formulas is not None

    if should_add_adsorption:
        print("\n  → Adding adsorption energy calculation")
        print(f"     Number of structures: {len(adsorption_structures)}")
        print(f"     Structure keys: {list(adsorption_structures.keys())}")
        print(f"     Adsorbate formulas: {adsorption_formulas}")
        if relax_before_adsorption:
            print(f"     Relaxation: ENABLED (relax → separate → SCF)")
        else:
            print(f"     Relaxation: DISABLED (separate → SCF)")

        from aiida.orm import Code, load_code

        # Load code
        code = load_code(code_label)

        # Get parameters - use adsorption-specific or fall back to slab/bulk parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
        slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

        ads_params = adsorption_parameters if adsorption_parameters is not None else slab_params
        ads_opts = adsorption_options if adsorption_options is not None else slab_opts
        ads_pot_map = adsorption_potential_mapping if adsorption_potential_mapping is not None else slab_pot_map
        ads_kpts = adsorption_kpoints_spacing if adsorption_kpoints_spacing is not None else slab_kpts

        # Backward compatibility: construct builder inputs from old-style if needed
        scf_builder_inputs = adsorption_scf_builder_inputs
        if scf_builder_inputs is None and ads_params is not None:
            # Auto-construct from old-style parameters
            scf_builder_inputs = {
                'parameters': {'incar': dict(ads_params)},
                'options': dict(ads_opts),
                'potential_family': potential_family,
                'potential_mapping': dict(ads_pot_map),
                'clean_workdir': clean_workdir,
                'settings': orm.Dict(dict={
                    'parser_settings': {
                        'add_trajectory': True,
                        'add_structure': True,
                        'add_kpoints': True,
                    }
                }),
            }
            if ads_kpts is not None:
                scf_builder_inputs['kpoints_spacing'] = ads_kpts

        # Add adsorption energy scatter task
        adsorption_task = wg.add_task(
            compute_adsorption_energies_scatter,
            name='compute_adsorption_energies_scatter',
            structures=adsorption_structures,
            adsorbate_formulas=adsorption_formulas,
            code=code,
            potential_family=potential_family,
            potential_mapping=ads_pot_map,

            # NEW: Relaxation parameters
            relax_before_adsorption=relax_before_adsorption,
            relax_builder_inputs=adsorption_relax_builder_inputs,
            scf_builder_inputs=scf_builder_inputs,

            # OLD: Backward compatibility fallback
            parameters=ads_params,
            options=ads_opts,
            kpoints_spacing=ads_kpts,
            clean_workdir=clean_workdir,
        )

        # Connect outputs
        wg.outputs.relaxed_complete_structures = adsorption_task.outputs.relaxed_complete_structures  # NEW
        wg.outputs.separated_structures = adsorption_task.outputs.separated_structures
        wg.outputs.substrate_energies = adsorption_task.outputs.substrate_energies
        wg.outputs.molecule_energies = adsorption_task.outputs.molecule_energies
        wg.outputs.complete_energies = adsorption_task.outputs.complete_energies
        wg.outputs.adsorption_energies = adsorption_task.outputs.adsorption_energies

        print(f"  ✓ Adsorption energy calculation enabled")
        if relax_before_adsorption:
            print(f"     Access relaxed structures via: wg.outputs.relaxed_complete_structures")
        print(f"     Access adsorption energies via: wg.outputs.adsorption_energies")
```

**Step 3: Verify syntax**

```bash
python3 -m py_compile teros/core/workgraph.py
```

Expected: No syntax errors

**Step 4: Test backward compatibility with existing example**

```bash
source ~/envs/aiida/bin/activate
python /home/thiagotd/git/PS-TEROS/.worktree/feature-adsorption-energy/examples/vasp/step_12_adsorption_energy.py
```

Expected: Script runs without errors (backward compatible)

**Step 5: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat: integrate relaxation support into build_core_workgraph with backward compat"
```

---

## Task 8: Update LaMnO3 Example to Use New API

**Files:**
- Modify: `teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py`

**Context:** Convert the LaMnO3 example to use new builder-based API with relaxation enabled.

**Step 1: Read current file**

```bash
cat teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py
```

**Step 2: Update to new builder-based API**

Replace the VASP configuration section (lines 65-120):

```python
    # 3. VASP Configuration
    print("\n3. VASP Configuration:")
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # Adsorption structures and formulas (dicts with matching keys)
    adsorption_structures = {
        'oh_lamno3': complete_structure,
    }

    adsorption_formulas = {
        'oh_lamno3': 'OH',  # The adsorbate to identify and separate
    }

    # Common settings
    common_potential_mapping = {
        'La': 'La',
        'Mn': 'Mn_pv',  # Use _pv for transition metals
        'O': 'O',
        'H': 'H',
    }

    common_options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'max_wallclock_seconds': 3600 * 24,  # 24 hours
        'queue_name': 'regular',
    }

    kpoints_spacing = 0.3  # Angstrom^-1

    # Relaxation inputs (Phase 1: relax complete system)
    relax_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 520,           # Cutoff energy
            'EDIFF': 1e-6,          # Electronic convergence
            'ISMEAR': 0,            # Gaussian smearing
            'SIGMA': 0.05,          # Small smearing width
            'ISPIN': 2,             # Spin-polarized
            'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments (La, Mn, O, H)
            'LORBIT': 11,           # DOSCAR and PROCAR
            'IBRION': 2,            # Conjugate gradient relaxation
            'NSW': 100,             # Max ionic steps
            'ISIF': 2,              # Relax ions only (keep cell fixed)
            'LWAVE': False,         # Don't write WAVECAR
            'LCHARG': False,        # Don't write CHGCAR
            'NCORE': 8,             # Parallelization
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
        'clean_workdir': False,
        'settings': orm.Dict(dict={
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }),
    }

    # SCF inputs (Phase 3: single-point calculations)
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 520,           # Cutoff energy
            'EDIFF': 1e-6,          # Electronic convergence
            'ISMEAR': 0,            # Gaussian smearing
            'SIGMA': 0.05,          # Small smearing width
            'ISPIN': 2,             # Spin-polarized
            'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',  # Initial magnetic moments
            'LORBIT': 11,           # DOSCAR and PROCAR
            # NSW=0 and IBRION=-1 will be set automatically by workflow
            'LWAVE': False,         # Don't write WAVECAR
            'LCHARG': False,        # Don't write CHGCAR
            'NCORE': 8,             # Parallelization
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
        'clean_workdir': False,
        'settings': orm.Dict(dict={
            'parser_settings': {
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }),
    }

    print(f"   K-points spacing: {kpoints_spacing} Å⁻¹")
    print(f"   Spin-polarized: Yes")
    print(f"   Relaxation: ENABLED (relax → separate → SCF)")
    print(f"   Relax max ionic steps: {relax_builder_inputs['parameters']['incar']['NSW']}")
```

**Step 3: Update workflow building section**

Replace the workflow building section (lines 125-155):

```python
    # 4. Build WorkGraph
    print("\n4. Building WorkGraph...")
    print("   Using preset: 'adsorption_energy'")
    print("   Workflow phases:")
    print("     Phase 1: Relax complete system (LaMnO3 + OH)")
    print("     Phase 2: Separate relaxed structure")
    print("     Phase 3: SCF calculations (3 single-point calculations)")
    print("       - Substrate (LaMnO3 from relaxed)")
    print("       - Molecule (OH from relaxed)")
    print("       - Complete (LaMnO3+OH from relaxed)")
    print("")
    print("   Formula: E_ads = E_complete - E_substrate - E_molecule")
    print("   Negative E_ads = favorable (exothermic) adsorption")
    print("   Positive E_ads = unfavorable (endothermic) adsorption")

    # Build workgraph using adsorption_energy preset with relaxation
    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,  # Keep files for analysis

        # Adsorption energy specific parameters
        adsorption_structures=adsorption_structures,
        adsorption_formulas=adsorption_formulas,

        # NEW: Relaxation + SCF builder inputs
        relax_before_adsorption=True,
        adsorption_relax_builder_inputs=relax_builder_inputs,
        adsorption_scf_builder_inputs=scf_builder_inputs,

        name='LaMnO3_OH_RelaxThenAdsorption',
    )

    print("   ✓ WorkGraph built successfully")
```

**Step 4: Update output information section**

Replace the output information section (lines 170-194):

```python
    print(f"\nExpected outputs:")
    print(f"  ")
    print(f"  1. Relaxed structure:")
    print(f"     - relaxed_complete_structures['oh_lamno3']: Relaxed La32Mn28O89H")
    print(f"  ")
    print(f"  2. Separated structures:")
    print(f"     - separated_structures['oh_lamno3']:")
    print(f"       * substrate: La32Mn28O88 (149 atoms)")
    print(f"       * molecule: OH (2 atoms)")
    print(f"       * complete: La32Mn28O89H (150 atoms)")
    print(f"  ")
    print(f"  3. Individual energies:")
    print(f"     - substrate_energies['oh_lamno3']: E(LaMnO3) from SCF")
    print(f"     - molecule_energies['oh_lamno3']: E(OH) from SCF")
    print(f"     - complete_energies['oh_lamno3']: E(LaMnO3+OH) from SCF")
    print(f"  ")
    print(f"  4. Adsorption energy:")
    print(f"     - adsorption_energies['oh_lamno3']: E_ads (eV)")
    print(f"  ")
    print(f"Expected E_ads for OH/LaMnO3:")
    print(f"  Literature range: -2 to -4 eV (DFT-PBE+U)")
    print(f"  (Negative = favorable adsorption)")
    print(f"\nWorkflow details:")
    print(f"  Phase 1: Relaxation of complete system")
    print(f"    - 1 VASP relaxation (NSW=100, IBRION=2)")
    print(f"  Phase 2: Structure separation (automatic)")
    print(f"  Phase 3: SCF calculations")
    print(f"    - 3 VASP SCF (NSW=0, single-point)")
    print(f"  Total VASP jobs: 4 (1 relax + 3 SCF)")
    print(f"\nCell handling:")
    print(f"  All systems use the SAME simulation cell")
    print(f"  This eliminates basis set superposition error (BSSE)")
    print(f"\nEstimated runtime:")
    print(f"  Phase 1 (relax): ~12-18 hours (150 atoms, ionic convergence)")
    print(f"  Phase 3 (SCF): ~2-4 hours each (3 parallel jobs)")
    print(f"  Total: ~12-18 hours (relax is bottleneck, SCF runs after)")
```

**Step 5: Verify syntax**

```bash
python3 -m py_compile teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py
```

Expected: No syntax errors

**Step 6: Commit**

```bash
git add teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py
git commit -m "refactor: update LaMnO3 example to use builder API with relaxation"
```

---

## Task 9: Create Test Example with Relaxation

**Files:**
- Create: `examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py`
- Create: `examples/adsorption_energy/test_relax_oh_ag111/README.md`

**Context:** Create a simple test case for the new relaxation workflow using Ag+OH system.

**Step 1: Create directory**

```bash
mkdir -p examples/adsorption_energy/test_relax_oh_ag111
```

**Step 2: Create test script**

Create `examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py`:

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
Test: Adsorption Energy with Relaxation - Ag(111) + OH

This script tests the new relaxation workflow:
1. Relax complete system (Ag + OH)
2. Separate relaxed structure
3. SCF calculations on all three components

System: Ag(111) + OH radical
"""

from aiida import load_profile, orm
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from teros.core.workgraph import build_core_workgraph


def create_ag_oh_structure():
    """Create Ag(111) slab with OH adsorbate for testing."""
    # Create Ag(111) slab (4 layers, 2x2)
    slab = fcc111('Ag', size=(2, 2, 4), vacuum=10.0)

    # Add OH on top site (directly above Ag atom)
    add_adsorbate(slab, 'OH', height=2.0, position='ontop')

    return slab


def main():
    """Run Ag + OH adsorption energy calculation with relaxation."""

    print("=" * 70)
    print("ADSORPTION ENERGY WITH RELAXATION TEST")
    print("System: Ag(111) + OH")
    print("=" * 70)

    # 1. Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # 2. Create structure
    print("\n2. Creating Ag(111) + OH structure...")
    ase_structure = create_ag_oh_structure()
    structure = orm.StructureData(ase=ase_structure)

    print(f"   Formula: {structure.get_formula()}")
    print(f"   Total atoms: {len(structure.sites)}")
    print("   ✓ Structure created")

    # 3. VASP Configuration
    print("\n3. VASP Configuration:")
    code_label = 'VASP-VTST-6.4.3@bohr'
    potential_family = 'PBE'

    # Common settings
    common_potential_mapping = {'Ag': 'Ag', 'O': 'O', 'H': 'H'}
    common_options = {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600 * 6,
        'queue_name': 'regular',
    }
    kpoints_spacing = 0.3

    # Relaxation builder inputs
    relax_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 1,
            'SIGMA': 0.2,
            'IBRION': 2,
            'NSW': 50,
            'ISIF': 2,
            'LWAVE': False,
            'LCHARG': False,
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
    }

    # SCF builder inputs
    scf_builder_inputs = {
        'parameters': {'incar': {
            'PREC': 'Accurate',
            'ENCUT': 400,
            'EDIFF': 1e-5,
            'ISMEAR': 1,
            'SIGMA': 0.2,
            # NSW=0, IBRION=-1 set automatically
            'LWAVE': False,
            'LCHARG': False,
        }},
        'options': common_options,
        'potential_family': potential_family,
        'potential_mapping': common_potential_mapping,
        'kpoints_spacing': kpoints_spacing,
    }

    print(f"   Code: {code_label}")
    print(f"   Relaxation: ENABLED")
    print(f"   Max ionic steps: {relax_builder_inputs['parameters']['incar']['NSW']}")

    # 4. Build WorkGraph
    print("\n4. Building WorkGraph...")

    wg = build_core_workgraph(
        workflow_preset='adsorption_energy',
        code_label=code_label,
        potential_family=potential_family,
        clean_workdir=False,

        adsorption_structures={'oh_ag': structure},
        adsorption_formulas={'oh_ag': 'OH'},

        relax_before_adsorption=True,
        adsorption_relax_builder_inputs=relax_builder_inputs,
        adsorption_scf_builder_inputs=scf_builder_inputs,

        name='Ag111_OH_RelaxTest',
    )

    print("   ✓ WorkGraph built")

    # 5. Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'=' * 70}")
    print("CALCULATION SUBMITTED")
    print(f"{'=' * 70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor:")
    print(f"  verdi process status {wg.pk}")
    print(f"\nExpected workflow:")
    print(f"  Phase 1: Relax Ag+OH → relaxed structure")
    print(f"  Phase 2: Separate → Ag substrate, OH molecule, complete")
    print(f"  Phase 3: SCF → 3 single-point calculations")
    print(f"\nOutputs:")
    print(f"  wg.outputs.relaxed_complete_structures['oh_ag']")
    print(f"  wg.outputs.adsorption_energies['oh_ag']")
    print(f"{'=' * 70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
```

**Step 3: Create README**

Create `examples/adsorption_energy/test_relax_oh_ag111/README.md`:

```markdown
# Test: Adsorption Energy with Relaxation

This test validates the new relaxation workflow for adsorption energy calculations.

## System

- **Substrate:** Ag(111) surface (4 layers, 2x2 supercell)
- **Adsorbate:** OH radical
- **Adsorption site:** Top site (directly above Ag atom)

## Workflow

**Phase 1: Relaxation**
- Relax complete Ag+OH system
- Plugin: `vasp.v2.relax`
- Settings: NSW=50, IBRION=2, ISIF=2

**Phase 2: Separation**
- Separate relaxed structure using connectivity analysis
- Components: Ag substrate, OH molecule, complete system

**Phase 3: SCF Calculations**
- Single-point calculations for all three components
- Plugin: `vasp.v2.vasp`
- Settings: NSW=0 (enforced), IBRION=-1 (enforced)

## Running

```bash
source ~/envs/aiida/bin/activate
python run_relax_adsorption.py
```

## Expected Results

- Relaxation should converge in < 50 steps
- Relaxed structure should show OH bond length ~0.97 Å
- Adsorption energy: E_ads ≈ -2.5 to -3.5 eV (exothermic)

## Validation

Check that:
1. Relaxation completed successfully
2. Structure separation identified OH correctly
3. All 3 SCF calculations completed
4. Adsorption energy is negative (favorable)

Access outputs:
```python
wg = load_node(PK)
relaxed = wg.outputs.relaxed_complete_structures['oh_ag']
E_ads = wg.outputs.adsorption_energies['oh_ag']
```
```

**Step 4: Make script executable**

```bash
chmod +x examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py
```

**Step 5: Verify syntax**

```bash
python3 -m py_compile examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py
```

Expected: No syntax errors

**Step 6: Commit**

```bash
git add examples/adsorption_energy/test_relax_oh_ag111/
git commit -m "test: add Ag+OH relaxation workflow test example"
```

---

## Task 10: Update Documentation

**Files:**
- Modify: `examples/vasp/README.md` (Step 12 section)
- Modify: `docs/adsorption_energy_module.md`

**Context:** Document the new relaxation + SCF workflow with examples.

**Step 1: Update examples/vasp/README.md**

Find Step 12 section and add after the existing content:

```markdown

### New Feature: Relaxation + SCF Workflow

**Optional relaxation before adsorption energy calculation:**

The adsorption energy module now supports optional relaxation of the complete system before structure separation and SCF calculations.

**Workflow comparison:**

| Mode | Phases | VASP Jobs | When to Use |
|------|--------|-----------|-------------|
| **Direct SCF** | separate → SCF | 3N | Input geometries are pre-relaxed or well-optimized |
| **Relax + SCF** | relax → separate → SCF | N + 3N | Input geometries need optimization before energy calculation |

**Example with relaxation:**

```python
# Relaxation builder inputs
relax_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'NSW': 100,
        'IBRION': 2,
        'ISIF': 2,
        # ... other INCAR tags
    }},
    'options': {'resources': {'num_machines': 1}},
    'potential_family': 'PBE',
    'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    'kpoints_spacing': 0.3,
}

# SCF builder inputs (NSW=0 enforced automatically)
scf_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        # ... other INCAR tags
    }},
    'options': {'resources': {'num_machines': 1}},
    'potential_family': 'PBE',
    'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    'kpoints_spacing': 0.3,
}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',

    adsorption_structures={'oh_lamno3': structure},
    adsorption_formulas={'oh_lamno3': 'OH'},

    # Enable relaxation
    relax_before_adsorption=True,
    adsorption_relax_builder_inputs=relax_inputs,
    adsorption_scf_builder_inputs=scf_inputs,
)

# Access relaxed structure
relaxed = wg.outputs.relaxed_complete_structures['oh_lamno3']
```

**Test examples:**
- Without relaxation: `examples/vasp/step_12_adsorption_energy.py`
- With relaxation: `examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py`
```

**Step 2: Update docs/adsorption_energy_module.md**

Add a new section after "Key Improvement":

```markdown

## Workflow Modes

### Mode 1: Direct SCF (Default)

**When to use:** Input structures are pre-relaxed or well-optimized

**Workflow:**
```
Input structure → Separate → SCF (3 jobs) → E_ads
```

**Example:**
```python
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={'system': structure},
    adsorption_formulas={'system': 'OH'},
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    adsorption_options={'resources': {'num_machines': 1}},
)
```

**VASP jobs:** 3N (substrate SCF, molecule SCF, complete SCF for each structure)

---

### Mode 2: Relax + SCF (New)

**When to use:** Input structures need geometry optimization before energy calculation

**Workflow:**
```
Input structure → Relax → Separate relaxed → SCF (3 jobs) → E_ads
```

**Example:**
```python
relax_inputs = {
    'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, ...}},
    'options': {'resources': {'num_machines': 1}},
    'potential_family': 'PBE',
    'potential_mapping': {...},
    'kpoints_spacing': 0.3,
}

scf_inputs = {
    'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 520, ...}},
    # NSW=0, IBRION=-1 enforced automatically
    'options': {'resources': {'num_machines': 1}},
    'potential_family': 'PBE',
    'potential_mapping': {...},
    'kpoints_spacing': 0.3,
}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={'system': structure},
    adsorption_formulas={'system': 'OH'},

    relax_before_adsorption=True,
    adsorption_relax_builder_inputs=relax_inputs,
    adsorption_scf_builder_inputs=scf_inputs,
)

# Access relaxed structure
relaxed = wg.outputs.relaxed_complete_structures['system']
```

**VASP jobs:** N + 3N (relax complete for each, then substrate/molecule/complete SCF)

**Key advantages:**
- Geometrically correct: Separates RELAXED structure, not input structure
- Relaxes only complete system (efficient)
- SCF is fast for all three components
- Eliminates geometry-related errors in E_ads

---

## API Reference

### New Parameters (build_core_workgraph)

**`relax_before_adsorption`** (bool, default=False)
- If True, relax complete structures before separation and SCF
- Uses `vasp.v2.relax` plugin

**`adsorption_relax_builder_inputs`** (dict, optional)
- Full builder dictionary for relaxation calculations
- Must include: parameters, options, potential_family, potential_mapping
- Example: `{'parameters': {'incar': {'NSW': 100, ...}}, 'options': {...}, ...}`

**`adsorption_scf_builder_inputs`** (dict, optional)
- Full builder dictionary for SCF calculations
- NSW=0 and IBRION=-1 enforced automatically
- Example: `{'parameters': {'incar': {'ENCUT': 520, ...}}, 'options': {...}, ...}`

### Backward Compatibility

Old parameter style still works (no relaxation):

```python
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    adsorption_structures={...},
    adsorption_formulas={...},
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520},  # Still works!
    adsorption_options={...},
    adsorption_potential_mapping={...},
)
```

Migration is optional and gradual.

---

## Testing

### Test Cases

**Test 1: Backward compatibility**
- File: `examples/vasp/step_12_adsorption_energy.py`
- Tests: Old API still works, no relaxation

**Test 2: Relaxation workflow**
- File: `examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py`
- Tests: New API with relaxation enabled

**Test 3: Production example**
- File: `teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py`
- Tests: Perovskite oxide with relaxation (real-world case)

### Validation Checklist

- [ ] Old scripts run unchanged (backward compatibility)
- [ ] Relaxation completes successfully
- [ ] Separation works on relaxed structures
- [ ] SCF calculations have NSW=0, IBRION=-1
- [ ] E_ads values are physically reasonable
- [ ] Relaxed structures show expected geometry changes
```

**Step 3: Commit**

```bash
git add examples/vasp/README.md docs/adsorption_energy_module.md
git commit -m "docs: document relaxation + SCF workflow with examples and API reference"
```

---

## Final Verification

**Step 1: Run all tests**

```bash
# Test backward compatibility
source ~/envs/aiida/bin/activate
verdi daemon restart
sleep 5

python examples/vasp/step_12_adsorption_energy.py
# Wait for completion
sleep 30
verdi process show <PK>  # Should show success [0]
```

**Step 2: Test new relaxation workflow**

```bash
python examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py
# Wait for completion (longer due to relaxation)
sleep 60
verdi process show <PK>  # Should show success [0]
```

**Step 3: Verify outputs**

```python
from aiida import load_profile, orm
load_profile('psteros')

wg = orm.load_node(PK)

# Check relaxed structure exists
assert 'relaxed_complete_structures' in wg.outputs
assert 'oh_ag' in wg.outputs.relaxed_complete_structures

# Check adsorption energy is negative (favorable)
E_ads = wg.outputs.adsorption_energies['oh_ag'].value
assert E_ads < 0, f"Expected negative E_ads, got {E_ads}"

print(f"✓ All tests passed!")
print(f"✓ E_ads = {E_ads:.3f} eV")
```

**Step 4: Final commit**

```bash
git add -A
git commit -m "chore: final verification and cleanup for relaxation + SCF workflow"
```

---

## Implementation Complete

**Summary of changes:**
1. ✅ Added `_build_vasp_inputs()` helper with backward compatibility
2. ✅ Updated `compute_adsorption_energies_scatter()` signature
3. ✅ Implemented Phase 1 (optional relaxation)
4. ✅ Updated Phase 2 (conditional structure input)
5. ✅ Modified Phase 3 (SCF with builder inputs)
6. ✅ Added `relaxed_complete_structures` to outputs
7. ✅ Integrated into `build_core_workgraph()`
8. ✅ Updated LaMnO3 example to new API
9. ✅ Created test example with relaxation
10. ✅ Updated documentation

**Files modified:** 4
**Files created:** 3
**Tests added:** 3 unit tests + 2 integration examples
**Backward compatible:** Yes ✓

**Next steps:**
- Run production test on LaMnO3 system
- Monitor computational cost (relax vs direct SCF)
- Consider adding DFT+U support in future enhancement
