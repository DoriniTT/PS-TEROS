# Adsorption Energy: Relaxation + SCF Design

**Date:** 2025-10-20
**Status:** Design Approved
**Implementation Branch:** feature-adsorption-energy

## Overview

Enhance the adsorption energy module to support relaxation of the complete system before performing SCF calculations on the separated components (substrate, molecule, complete).

**Current behavior:**
- Single-step workflow: separate → SCF (all 3 components)
- Uses original input geometries for all calculations
- Simple parameter dict API

**New behavior:**
- Optional 3-phase workflow: relax → separate → SCF
- Relaxes complete system before separation (geometrically correct)
- Builder-based API for full control over VASP settings
- Backward compatible with existing scripts

## Requirements

1. **Relaxation control:** Add `relax_before_adsorption` flag
2. **Plugin selection:** Use `vasp.v2.relax` for relaxation, `vasp.v2.vasp` for SCF
3. **Builder inputs:** Accept full builder dicts instead of just INCAR parameters
4. **Backward compatibility:** Keep existing parameter dict API working
5. **Output exposure:** Expose relaxed structures for user inspection

## Architecture

### Overall Workflow (Sequential 3-Phase)

```
Phase 1: Optional Relaxation (if relax_before_adsorption=True)
  ├─ Input: Complete structures (substrate + adsorbate)
  ├─ Plugin: vasp.v2.relax
  ├─ Config: relax_builder_inputs (NSW, IBRION, ISIF)
  ├─ Parallelization: Across all structures in scatter
  └─ Output: relaxed_complete_structures

Phase 2: Structure Separation
  ├─ Input: Relaxed structures (Phase 1) OR original structures
  ├─ Function: separate_adsorbate_structure() [existing]
  ├─ Algorithm: Connectivity analysis with subgraph approach
  └─ Output: separated_structures (substrate, molecule, complete)

Phase 3: SCF Calculations
  ├─ Input: Separated structures (substrate, molecule, complete)
  ├─ Plugin: vasp.v2.vasp
  ├─ Config: scf_builder_inputs (NSW=0, IBRION=-1 enforced)
  ├─ Parallelization: 3N jobs (substrate, molecule, complete for each)
  └─ Output: Energies → E_ads calculation
```

**Why sequential?**
- Separation must use RELAXED geometry (geometrically correct)
- Relaxation changes bond lengths/angles that affect adsorbate identification
- Sequential phases ensure data flows correctly: relax → separate → SCF

**Why relax only complete system?**
- Efficient: 1 relaxation instead of 3
- Substrate geometry comes from separated relaxed structure
- Molecule geometry comes from separated relaxed structure
- SCF is fast for all three components

## API Design

### `compute_adsorption_energies_scatter()` Signature

```python
@task.graph
def compute_adsorption_energies_scatter(
    # Existing parameters
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
)]
```

### Builder Input Structure

Both `relax_builder_inputs` and `scf_builder_inputs` follow the same structure:

```python
builder_inputs = {
    'structure': orm.StructureData,        # Set by workflow automatically
    'code': orm.Code,                      # Set by workflow automatically
    'parameters': {'incar': {              # VASP INCAR parameters
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'NSW': 100,                        # For relax; forced to 0 for SCF
        'IBRION': 2,                       # For relax; forced to -1 for SCF
        'ISIF': 2,
        'ISPIN': 2,
        'MAGMOM': '...',
        # ... other INCAR tags
    }},
    'options': {                           # Scheduler options
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 40,
        },
        'max_wallclock_seconds': 3600 * 24,
        'queue_name': 'regular',
    },
    'potential_family': 'PBE',
    'potential_mapping': {
        'La': 'La',
        'Mn': 'Mn_pv',
        'O': 'O',
        'H': 'H',
    },
    'kpoints_spacing': 0.3,                # Angstrom^-1
    'clean_workdir': True,
    'settings': orm.Dict(dict={
        'parser_settings': {
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }),
}
```

**Key differences between relax and SCF:**
- Relax: `NSW > 0`, `IBRION` in [1, 2, 3], `ISIF` for cell/ion control
- SCF: `NSW = 0`, `IBRION = -1` (forced by workflow)

### Backward Compatibility Strategy

**Priority order for constructing inputs:**

```python
# For SCF calculations
if scf_builder_inputs is not None:
    # Use new-style builder inputs
    inputs = scf_builder_inputs
elif parameters is not None:
    # Construct from old-style parameters
    inputs = {
        'parameters': {'incar': dict(parameters)},
        'options': dict(options),
        'potential_family': potential_family,
        'potential_mapping': dict(potential_mapping),
        'kpoints_spacing': kpoints_spacing,
        'clean_workdir': clean_workdir,
    }
else:
    raise ValueError("Must provide either scf_builder_inputs or parameters")
```

**Result:** Existing scripts using `adsorption_parameters` continue to work unchanged.

## Implementation Details

### Helper Function: `_build_vasp_inputs()`

Location: `teros/core/adsorption_energy.py`

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
    Build VASP inputs from either builder_inputs or old-style parameters.

    Args:
        structure: Structure to calculate
        code: VASP code
        builder_inputs: New-style builder dict (priority)
        parameters: Old-style INCAR dict (fallback)
        options: Old-style scheduler dict (fallback)
        potential_family: Pseudopotential family
        potential_mapping: Element → potential mapping
        kpoints_spacing: K-points spacing (Angstrom^-1)
        clean_workdir: Clean remote directory
        force_scf: If True, force NSW=0 and IBRION=-1

    Returns:
        Complete input dict for VASP WorkChain
    """
    if builder_inputs is not None:
        # Use builder_inputs, merge structure/code
        inputs = dict(builder_inputs)
        inputs['structure'] = structure
        inputs['code'] = code
    else:
        # Construct from old-style parameters
        inputs = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
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

    # Force SCF mode if requested
    if force_scf:
        inputs['parameters']['incar']['NSW'] = 0
        inputs['parameters']['incar']['IBRION'] = -1

    return inputs
```

### Workflow Implementation

Inside `compute_adsorption_energies_scatter()`:

```python
# Get workflow plugins
vasp_relax = WorkflowFactory('vasp.v2.relax')
vasp_scf = WorkflowFactory('vasp.v2.vasp')
relax_task_cls = task(vasp_relax)
scf_task_cls = task(vasp_scf)

# Phase 1: Optional Relaxation
if relax_before_adsorption:
    relaxed_structures: dict[str, orm.StructureData] = {}

    for key, structure in structures.items():
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
            force_scf=False,
        )

        relax_task = relax_task_cls(**relax_inputs)
        relaxed_structures[key] = relax_task.outputs.structure

    structures_to_separate = relaxed_structures
else:
    relaxed_structures = {}
    structures_to_separate = structures

# Phase 2: Separation (existing logic)
separated_dict: dict[str, dict] = {}
for key, structure in structures_to_separate.items():
    separated = separate_adsorbate_structure(
        structure=structure,
        adsorbate_formula=orm.Str(adsorbate_formulas[key])
    )
    separated_dict[key] = {
        'substrate': separated.substrate,
        'molecule': separated.molecule,
        'complete': separated.complete,
    }

# Phase 3: SCF calculations
substrate_energies: dict[str, orm.Float] = {}
molecule_energies: dict[str, orm.Float] = {}
complete_energies: dict[str, orm.Float] = {}

for key, separated in separated_dict.items():
    # Helper to create SCF inputs
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
            force_scf=True,  # Ensures NSW=0, IBRION=-1
        )

    # Substrate SCF
    substrate_calc = scf_task_cls(**create_scf_inputs(separated['substrate']))
    substrate_energies[key] = extract_total_energy(energies=substrate_calc.misc).result

    # Molecule SCF
    molecule_calc = scf_task_cls(**create_scf_inputs(separated['molecule']))
    molecule_energies[key] = extract_total_energy(energies=molecule_calc.misc).result

    # Complete SCF
    complete_calc = scf_task_cls(**create_scf_inputs(separated['complete']))
    complete_energies[key] = extract_total_energy(energies=complete_calc.misc).result

# Phase 4: Calculate adsorption energies (existing logic)
adsorption_energies: dict[str, orm.Float] = {}
for key in structures.keys():
    E_ads = calculate_adsorption_energy(
        E_complete=complete_energies[key],
        E_substrate=substrate_energies[key],
        E_molecule=molecule_energies[key],
    ).result
    adsorption_energies[key] = E_ads

# Return all outputs
return {
    'relaxed_complete_structures': relaxed_structures,  # NEW (empty dict if relax=False)
    'separated_structures': separated_dict,
    'substrate_energies': substrate_energies,
    'molecule_energies': molecule_energies,
    'complete_energies': complete_energies,
    'adsorption_energies': adsorption_energies,
}
```

## Integration with `build_core_workgraph()`

### New Parameters

```python
def build_core_workgraph(
    # ... existing 80+ parameters ...

    # Adsorption energy (existing)
    adsorption_structures: dict | None = None,
    adsorption_formulas: dict | None = None,

    # NEW: Relaxation control
    relax_before_adsorption: bool = False,
    adsorption_relax_builder_inputs: dict | None = None,
    adsorption_scf_builder_inputs: dict | None = None,

    # DEPRECATED (kept for backward compatibility)
    adsorption_parameters: dict | None = None,
    adsorption_options: dict | None = None,
    adsorption_potential_mapping: dict | None = None,
    adsorption_kpoints_spacing: float | None = None,
):
```

### Workgraph Integration

```python
if should_add_adsorption:
    # Backward compatibility: construct builder inputs from old-style if needed
    if adsorption_scf_builder_inputs is None and adsorption_parameters is not None:
        scf_builder_inputs = {
            'parameters': {'incar': dict(adsorption_parameters)},
            'options': dict(adsorption_options) if adsorption_options else dict(slab_opts),
            'potential_family': potential_family,
            'potential_mapping': dict(adsorption_potential_mapping) if adsorption_potential_mapping else dict(slab_pot_map),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict={
                'parser_settings': {
                    'add_trajectory': True,
                    'add_structure': True,
                    'add_kpoints': True,
                }
            }),
        }
        if adsorption_kpoints_spacing is not None:
            scf_builder_inputs['kpoints_spacing'] = adsorption_kpoints_spacing
    else:
        scf_builder_inputs = adsorption_scf_builder_inputs

    # Add adsorption task
    adsorption_task = wg.add_task(
        compute_adsorption_energies_scatter,
        name='compute_adsorption_energies_scatter',
        structures=adsorption_structures,
        adsorbate_formulas=adsorption_formulas,
        code=code,
        potential_family=potential_family,

        # NEW: Relaxation parameters
        relax_before_adsorption=relax_before_adsorption,
        relax_builder_inputs=adsorption_relax_builder_inputs,
        scf_builder_inputs=scf_builder_inputs,

        # OLD: Fallback parameters (for backward compat)
        parameters=adsorption_parameters,
        options=adsorption_options,
        potential_mapping=adsorption_potential_mapping,
        kpoints_spacing=adsorption_kpoints_spacing,
        clean_workdir=clean_workdir,
    )

    # Connect outputs
    wg.outputs.relaxed_complete_structures = adsorption_task.outputs.relaxed_complete_structures  # NEW
    wg.outputs.separated_structures = adsorption_task.outputs.separated_structures
    wg.outputs.substrate_energies = adsorption_task.outputs.substrate_energies
    wg.outputs.molecule_energies = adsorption_task.outputs.molecule_energies
    wg.outputs.complete_energies = adsorption_task.outputs.complete_energies
    wg.outputs.adsorption_energies = adsorption_task.outputs.adsorption_energies
```

## User-Facing Examples

### Example 1: New API with Relaxation

```python
# LaMnO3 + OH with relaxation before adsorption energy

relax_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'NSW': 100,
        'IBRION': 2,
        'ISIF': 2,
        'ISPIN': 2,
        'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',
        'LORBIT': 11,
    }},
    'options': {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600 * 24,
        'queue_name': 'regular',
    },
    'potential_family': 'PBE',
    'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    'kpoints_spacing': 0.3,
}

scf_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        # NSW=0, IBRION=-1 will be set automatically
        'ISPIN': 2,
        'MAGMOM': '32*0.6 28*4.0 89*0.6 1*0.0',
    }},
    'options': {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600 * 6,
        'queue_name': 'regular',
    },
    'potential_family': 'PBE',
    'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    'kpoints_spacing': 0.3,
}

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    clean_workdir=False,

    adsorption_structures={'oh_lamno3': structure},
    adsorption_formulas={'oh_lamno3': 'OH'},

    # NEW: Enable relaxation
    relax_before_adsorption=True,
    adsorption_relax_builder_inputs=relax_inputs,
    adsorption_scf_builder_inputs=scf_inputs,

    name='LaMnO3_OH_RelaxThenAdsorption',
)

wg.submit(wait=False)

# Access relaxed structure
print(wg.outputs.relaxed_complete_structures['oh_lamno3'])

# Access adsorption energy
print(wg.outputs.adsorption_energies['oh_lamno3'])
```

### Example 2: Backward Compatible (Old API)

```python
# Existing scripts continue to work (no relaxation, old parameter style)

wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',

    adsorption_structures={'oh_ag': structure},
    adsorption_formulas={'oh_ag': 'OH'},

    # OLD API still works
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6},
    adsorption_options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40}},
    adsorption_potential_mapping={'Ag': 'Ag', 'O': 'O', 'H': 'H'},
    adsorption_kpoints_spacing=0.3,
)

wg.submit(wait=False)
```

### Example 3: Mixed (Relax with Old-Style Params)

```python
# Use relaxation but keep old parameter style
wg = build_core_workgraph(
    workflow_preset='adsorption_energy',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',

    adsorption_structures={'oh_lamno3': structure},
    adsorption_formulas={'oh_lamno3': 'OH'},

    # Enable relaxation (uses same params for both relax and SCF)
    relax_before_adsorption=True,

    # Old-style params (used for both relax and SCF)
    adsorption_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'NSW': 100, 'IBRION': 2},
    adsorption_options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40}},
    adsorption_potential_mapping={'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
    adsorption_kpoints_spacing=0.3,
)
```

## Migration Path

**Phase 1: Implement with backward compatibility (this design)**
- Keep `adsorption_parameters`, `adsorption_options`, etc.
- Add new `relax_before_adsorption`, `adsorption_relax_builder_inputs`, `adsorption_scf_builder_inputs`
- Auto-convert old style to new style internally
- **Result:** All existing scripts work unchanged

**Phase 2: Documentation update**
- Update examples to use new builder-based API
- Show both old and new styles in docs
- Recommend new style for new projects

**Phase 3: Deprecation (future, optional)**
- Add deprecation warnings for old-style params
- Encourage migration to builder inputs
- Keep old style working for 1-2 major versions

## Testing Strategy

### Unit Tests
1. `test_build_vasp_inputs_new_style()` - Test builder input construction
2. `test_build_vasp_inputs_old_style()` - Test backward compatibility
3. `test_build_vasp_inputs_force_scf()` - Test NSW=0, IBRION=-1 enforcement

### Integration Tests
1. **Test 1:** `examples/adsorption_energy/test_oh_ag111/` - Update to test new API
2. **Test 2:** Create `examples/adsorption_energy/test_relax_oh_ag111/` - Test with `relax_before_adsorption=True`
3. **Test 3:** Update `teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py` to use new API

### Validation Criteria
- All existing tests pass (backward compatibility verified)
- New relaxation workflow completes successfully
- Relaxed structures differ from input structures (geometry changed)
- Final E_ads matches expected range for test systems

## File Modifications

### Files to Modify
1. `teros/core/adsorption_energy.py`
   - Add `_build_vasp_inputs()` helper function
   - Modify `compute_adsorption_energies_scatter()` signature and implementation
   - Add Phase 1 (relaxation) logic
   - Update Phase 3 to use `vasp.v2.vasp` plugin

2. `teros/core/workgraph.py`
   - Add new parameters to `build_core_workgraph()`
   - Add backward compatibility logic
   - Update adsorption task connection with new outputs

3. `examples/vasp/step_12_adsorption_energy.py`
   - Add example of new API usage (optional, for demonstration)

4. `teros/experimental/adsorption_energy/lamno3/run_lamno3_oh_adsorption.py`
   - Update to use new builder-based API with relaxation

### Files to Create
1. `examples/adsorption_energy/test_relax_oh_ag111/run_relax_adsorption.py`
   - New test for relaxation workflow

2. `examples/adsorption_energy/test_relax_oh_ag111/README.md`
   - Documentation for relaxation test

## Timeline Estimate

- **Phase 1:** Implement `_build_vasp_inputs()` helper - 30 min
- **Phase 2:** Modify `compute_adsorption_energies_scatter()` - 1.5 hours
- **Phase 3:** Update `build_core_workgraph()` integration - 1 hour
- **Phase 4:** Update examples and tests - 1 hour
- **Phase 5:** Testing and validation - 2 hours

**Total:** ~6 hours implementation + testing

## Success Criteria

1. **Backward compatibility:** All existing scripts (Ag+OH, LaMnO3+OH) run unchanged
2. **Relaxation works:** `relax_before_adsorption=True` completes successfully
3. **Builder inputs work:** New API with full builder dicts works correctly
4. **Geometry changes:** Relaxed structures show expected geometry changes
5. **Energy correctness:** E_ads values are physically reasonable
6. **Documentation:** Design doc, code comments, and examples are clear

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Breaking existing scripts | High | Rigorous backward compatibility testing |
| Complex API | Medium | Clear documentation, examples for both APIs |
| WorkGraph plugin compatibility | High | Test with both vasp.v2.relax and vasp.v2.vasp |
| Performance regression | Low | Relaxation is optional, old behavior is default |
| NSW=0 not enforced | Medium | Unit test for force_scf parameter |

## Future Enhancements

1. **Selective relaxation:** Option to relax substrate/molecule separately
2. **Constraint support:** Add VASP constraints (selective dynamics) to builder inputs
3. **Multi-step relaxation:** Coarse → fine relaxation sequence
4. **DFT+U support:** Explicit handling of Hubbard U parameters in builder
5. **Alternative separation:** Support user-provided separation instead of automatic

## References

- AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/
- AiiDA-WorkGraph scatter-gather: https://aiida-workgraph.readthedocs.io/
- Original adsorption energy design: `docs/adsorption_energy_module.md`
- PSTEROS workflow presets: `teros/core/workgraph.py`
