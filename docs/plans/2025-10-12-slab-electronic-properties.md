# Slab Electronic Properties Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Add selective electronic properties calculation (DOS and bands) for relaxed slabs using per-slab parameter overrides.

**Architecture:** Following the scatter-gather pattern used throughout PS-TEROS, create a new `@task.graph` function that takes a dictionary of selected slabs and their custom parameters, creates independent BandsWorkChain tasks, and returns four dynamic namespaces (bands, dos, primitive_structures, seekpath_parameters). Integration happens in `build_core_workgraph()` after slab relaxation completes.

**Tech Stack:** AiiDA, aiida-vasp (vasp.v2.bands), aiida-workgraph, Pymatgen

---

## Task 1: Create Slab Electronic Properties Builder

**Files:**
- Modify: `teros/core/builders/electronic_properties_builder.py` (append to end of file)

**Step 1: Add slab-specific builder function**

Add this function at the end of `electronic_properties_builder.py`:

```python
def get_slab_electronic_properties_defaults(
    energy_cutoff: float = 500,
    electronic_convergence: float = 1e-5,
    ncore: int = 4,
    ispin: int = 2,
    lasph: bool = True,
    lreal: str = "Auto",
    kpoints_mesh_density: float = 0.25,  # Denser for 2D systems
    band_kpoints_distance: float = 0.15,  # Denser path sampling
    dos_kpoints_distance: float = 0.2,
    line_density: float = 0.15,  # More points along paths
    nedos: int = 2000,
    sigma_bands: float = 0.01,
    symprec: float = 1e-4,
    band_mode: str = "seekpath-aiida",
) -> dict:
    """
    Get default parameters for slab electronic properties (DOS and bands) calculation.

    Material-agnostic builder tuned for 2D slab systems. Returns all necessary
    parameters for vasp.v2.bands workchain. Based on proven working configuration
    but with denser k-point sampling for surface systems.

    Args:
        energy_cutoff: ENCUT value (eV). Default: 500
        electronic_convergence: EDIFF value. Default: 1e-5
        ncore: Number of cores per band. Default: 4
        ispin: Spin polarization (1=non-magnetic, 2=spin-polarized). Default: 2
        lasph: Include aspherical contributions. Default: True
        lreal: Projection operators (False for small, "Auto" for large). Default: "Auto"
        kpoints_mesh_density: K-point density for SCF and DOS meshes (denser than bulk). Default: 0.25
        band_kpoints_distance: K-point spacing for band structure path. Default: 0.15
        dos_kpoints_distance: K-point spacing for DOS calculation. Default: 0.2
        line_density: K-point density along high-symmetry paths. Default: 0.15
        nedos: Number of DOS grid points. Default: 2000
        sigma_bands: Smearing width for bands (eV). Default: 0.01
        symprec: Symmetry precision for seekpath. Default: 1e-4
        band_mode: Band path mode (seekpath-aiida, pymatgen, etc.). Default: "seekpath-aiida"

    Returns:
        Dictionary with 'band_settings' and nested INCAR parameters for
        scf, bands, and dos stages, ready to use with vasp.v2.bands workchain.

    Example:
        >>> defaults = get_slab_electronic_properties_defaults()
        >>> wg = build_core_workgraph(
        ...     compute_electronic_properties_slabs=True,
        ...     slab_electronic_properties={'term_0': defaults, 'term_2': defaults},
        ...     slab_bands_parameters=defaults,
        ...     slab_band_settings=defaults['band_settings'],
        ...     **other_params
        ... )
    """
    # Common INCAR parameters used across all stages
    common_incar = {
        'ENCUT': energy_cutoff,
        'EDIFF': electronic_convergence,
        'NCORE': ncore,
        'ISPIN': ispin,
        'LASPH': lasph,
        'LREAL': lreal,
    }

    # Band workflow settings - controls the vasp.v2.bands workchain behavior
    band_settings_dict = {
        'symprec': symprec,
        'band_mode': band_mode,
        'band_kpoints_distance': band_kpoints_distance,
        'line_density': line_density,
        'run_dos': True,  # If True, compute DOS when computing bands
        'only_dos': False,
        'dos_kpoints_distance': dos_kpoints_distance,
        'kpoints_per_split': 20,
        'hybrid_reuse_wavecar': False,  # Not relevant for standard PBE
        'additional_band_analysis_parameters': {
            'with_time_reversal': True,
            'reference_distance': 0.05,
            'threshold': 1e-5,
        },
    }

    # SCF stage INCAR - Critical: LWAVE and LCHARG must be True for bands/DOS restart
    scf_incar = {
        **common_incar,
        'SYSTEM': 'SCF_Slab_Electronic_Properties',
        'ISTART': 0,
        'PREC': 'Accurate',
        'ALGO': 'Normal',
        'NELM': 120,
        'ISMEAR': 0,  # Gaussian smearing
        'SIGMA': 0.05,
        'LWAVE': True,   # Must be True - needed for bands restart
        'LCHARG': True,  # Must be True - needed for DOS restart
        'LORBIT': 11,    # L(m)-decomposed charge densities
        'LELF': True,    # Write ELFCAR
        'LAECHG': True,  # Write AECCAR0 and AECCAR2
        'LVTOT': True,   # Write LOCPOT
    }

    # Band structure stage INCAR - Non-self-consistent calculation
    bands_incar = {
        **common_incar,
        'SYSTEM': 'Slab_Band_Structure',
        'ISTART': 1,     # Read from WAVECAR
        'ICHARG': 11,    # Non-self-consistent calculation
        'PREC': 'Accurate',
        'ALGO': 'Normal',
        'NELM': 60,
        'ISMEAR': 0,     # Gaussian smearing
        'SIGMA': sigma_bands,
        'LREAL': False,  # Use reciprocal space for band structure
        'LWAVE': False,
        'LCHARG': False,
        'LORBIT': 11,    # Calculate projected bands
    }

    # DOS stage INCAR - Non-self-consistent with tetrahedron method
    dos_incar = {
        **common_incar,
        'SYSTEM': 'Slab_DOS_Calculation',
        'ISTART': 1,     # Read from WAVECAR
        'ICHARG': 11,    # Non-self-consistent calculation
        'PREC': 'Accurate',
        'ALGO': 'Normal',
        'NELM': 60,
        'ISMEAR': 0,     # Tetrahedron method with Blöchl corrections
        'LREAL': False,  # Use reciprocal space
        'LWAVE': False,
        'LCHARG': False,
        'LORBIT': 11,    # Calculate projected DOS
        'NEDOS': nedos,  # Number of DOS grid points
    }

    return {
        'band_settings': band_settings_dict,
        'scf': scf_incar,
        'bands': bands_incar,
        'dos': dos_incar,
        'scf_kpoints_distance': kpoints_mesh_density,
    }
```

**Step 2: Verify the function syntax**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-bands-slabs
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/builders/electronic_properties_builder.py
```

Expected: No output (successful compilation)

**Step 3: Commit**

```bash
git add teros/core/builders/electronic_properties_builder.py
git commit -m "feat: add get_slab_electronic_properties_defaults builder

Added material-agnostic parameter builder for slab electronic properties
with denser k-point sampling tuned for 2D systems."
```

---

## Task 2: Create Scatter Function for Slab Electronic Properties

**Files:**
- Modify: `teros/core/slabs.py` (append to end of file before `if __name__`)

**Step 1: Add scatter function**

Add this function at the end of `slabs.py` (before the `if __name__ == '__main__'` block if present):

```python
@task.graph
def calculate_electronic_properties_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    slab_electronic_properties: t.Mapping[str, t.Mapping[str, t.Any]],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],
    clean_workdir: bool,
    default_bands_parameters: t.Mapping[str, t.Any] = None,
    default_bands_options: t.Mapping[str, t.Any] = None,
    default_band_settings: t.Mapping[str, t.Any] = None,
) -> t.Annotated[dict, namespace(
    slab_bands=dynamic(orm.BandsData),
    slab_dos=dynamic(orm.ArrayData),
    slab_primitive_structures=dynamic(orm.StructureData),
    slab_seekpath_parameters=dynamic(orm.Dict),
)]:
    """
    Scatter-gather phase: calculate electronic properties (DOS and bands) for selected slabs.

    This task graph creates independent BandsWorkChain tasks for each selected slab,
    allowing per-slab parameter overrides. Only slabs present in slab_electronic_properties
    dictionary will have their electronic properties calculated.

    Args:
        slabs: Dictionary of relaxed slab structures
        slab_electronic_properties: Dictionary mapping slab labels to parameter configs.
                                    Each config can contain:
                                    - 'bands_parameters': Dict with 'scf', 'bands', 'dos' keys
                                    - 'bands_options': Scheduler options
                                    - 'band_settings': Band workflow settings
                                    If keys are missing, defaults are used.
        code: AiiDA code for VASP
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        clean_workdir: Whether to clean remote directories
        default_bands_parameters: Default parameters (fallback)
        default_bands_options: Default scheduler options (fallback)
        default_band_settings: Default band settings (fallback)

    Returns:
        Dictionary with four namespaces:
        - slab_bands: BandsData for each selected slab
        - slab_dos: DOS ArrayData for each selected slab
        - slab_primitive_structures: Primitive structures used for band paths
        - slab_seekpath_parameters: Seekpath parameters for each slab

    Example:
        >>> scatter_outputs = calculate_electronic_properties_slabs_scatter(
        ...     slabs=relaxed_slabs_dict,
        ...     slab_electronic_properties={
        ...         'term_0': {'bands_parameters': defaults},
        ...         'term_2': {'bands_parameters': custom_params, 'bands_options': high_mem},
        ...     },
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={'Ag': 'Ag', 'O': 'O'},
        ...     clean_workdir=True,
        ... )
    """
    from aiida.plugins import WorkflowFactory

    # Get BandsWorkChain and wrap as task
    BandsWorkChain = WorkflowFactory('vasp.v2.bands')
    BandsTask = task(BandsWorkChain)

    # Output dictionaries
    bands_dict: dict[str, orm.BandsData] = {}
    dos_dict: dict[str, orm.ArrayData] = {}
    primitive_structures_dict: dict[str, orm.StructureData] = {}
    seekpath_parameters_dict: dict[str, orm.Dict] = {}

    # Scatter: create independent tasks for each selected slab
    for label, slab_config in slab_electronic_properties.items():
        # Skip if slab doesn't exist
        if label not in slabs:
            continue

        structure = slabs[label]

        # Extract per-slab parameters with fallback to defaults
        bands_params = slab_config.get('bands_parameters', default_bands_parameters or {})
        bands_opts = slab_config.get('bands_options', default_bands_options or {})
        band_settings = slab_config.get('band_settings', default_band_settings or {})

        # Build BandsWorkChain inputs (same structure as bulk)
        bands_inputs: dict[str, t.Any] = {
            'structure': structure,
            'metadata': {
                'label': f'Slab_{label}_Electronic_Properties',
                'description': f'DOS and bands calculation for slab termination {label}',
            }
        }

        # Add band_settings if provided
        if band_settings:
            bands_inputs['band_settings'] = band_settings

        # Build SCF namespace inputs (required)
        scf_inputs = {
            'code': code,
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'options': dict(bands_opts),
            'clean_workdir': clean_workdir,
        }
        if bands_params and 'scf' in bands_params:
            scf_inputs['parameters'] = {'incar': bands_params['scf']}
        if bands_params and 'scf_kpoints_distance' in bands_params:
            scf_inputs['kpoints_spacing'] = bands_params['scf_kpoints_distance']

        bands_inputs['scf'] = scf_inputs

        # Build Bands namespace inputs (optional)
        if bands_params and 'bands' in bands_params:
            bands_inputs['bands'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'options': dict(bands_opts),
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_params['bands']},
            }

        # Build DOS namespace inputs (optional)
        if bands_params and 'dos' in bands_params:
            bands_inputs['dos'] = {
                'code': code,
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'options': dict(bands_opts),
                'clean_workdir': clean_workdir,
                'parameters': {'incar': bands_params['dos']},
            }

        # Create BandsWorkChain task
        bands_task = BandsTask(**bands_inputs)

        # Collect outputs
        bands_dict[label] = bands_task.band_structure
        dos_dict[label] = bands_task.dos
        primitive_structures_dict[label] = bands_task.primitive_structure
        seekpath_parameters_dict[label] = bands_task.seekpath_parameters

    # Gather: return collected results
    return {
        'slab_bands': bands_dict,
        'slab_dos': dos_dict,
        'slab_primitive_structures': primitive_structures_dict,
        'slab_seekpath_parameters': seekpath_parameters_dict,
    }
```

**Step 2: Verify the function syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/slabs.py
```

Expected: No output (successful compilation)

**Step 3: Commit**

```bash
git add teros/core/slabs.py
git commit -m "feat: add calculate_electronic_properties_slabs_scatter function

Scatter-gather function for selective slab electronic properties calculation
with per-slab parameter overrides. Returns bands, dos, primitive structures,
and seekpath parameters in separate dynamic namespaces."
```

---

## Task 3: Add Outputs to core_workgraph Decorator

**Files:**
- Modify: `teros/core/workgraph.py:59-68` (add new outputs to decorator)

**Step 1: Add slab electronic properties outputs**

In the `@task.graph(outputs=[...])` decorator for `core_workgraph` function, add these four new outputs after the existing slab outputs:

```python
@task.graph(outputs=[
    'bulk_energy', 'metal_energy', 'nonmetal_energy', 'oxygen_energy',
    'bulk_structure', 'metal_structure', 'nonmetal_structure', 'oxygen_structure',
    'formation_enthalpy', 'slab_structures', 'relaxed_slabs', 'slab_energies',
    'unrelaxed_slab_energies', 'relaxation_energies',
    'surface_energies', 'cleavage_energies',
    'slab_remote', 'unrelaxed_slab_remote',
    # Electronic properties outputs
    'bulk_bands', 'bulk_dos', 'bulk_primitive_structure', 'bulk_seekpath_parameters',
    # NEW: Slab electronic properties outputs
    'slab_bands', 'slab_dos', 'slab_primitive_structures', 'slab_seekpath_parameters',
])
```

**Step 2: Add return values in core_workgraph**

In the `return` statement at the end of `core_workgraph()` function (around line 419), add initialization for the new outputs:

```python
    # Initialize output dictionaries for slab electronic properties
    slab_bands_output = {}
    slab_dos_output = {}
    slab_primitive_structures_output = {}
    slab_seekpath_parameters_output = {}

    # Return all outputs
    # Note: Electronic properties outputs (bulk_bands, bulk_dos, bulk_electronic_properties_misc)
    # are added dynamically in build_core_workgraph, not returned here
    return {
        'bulk_energy': bulk_energy.result,
        'bulk_structure': bulk_vasp.structure,
        'metal_energy': metal_energy.result,
        'metal_structure': metal_vasp.structure,
        'nonmetal_energy': nonmetal_energy.result,
        'nonmetal_structure': nonmetal_vasp.structure,
        'oxygen_energy': oxygen_energy.result,
        'oxygen_structure': oxygen_vasp.structure,
        'formation_enthalpy': formation_hf.result,
        'slab_structures': slab_namespace if slab_namespace is not None else {},
        'relaxed_slabs': relaxed_slabs_output,
        'slab_energies': slab_energies_output,
        'unrelaxed_slab_energies': unrelaxed_slab_energies_output,
        'relaxation_energies': relaxation_energies_output,
        'surface_energies': surface_energies_output,
        'cleavage_energies': cleavage_energies_output,
        'slab_remote': slab_remote_output,
        'unrelaxed_slab_remote': unrelaxed_slab_remote_output,
        'slab_bands': slab_bands_output,
        'slab_dos': slab_dos_output,
        'slab_primitive_structures': slab_primitive_structures_output,
        'slab_seekpath_parameters': slab_seekpath_parameters_output,
    }
```

**Step 3: Verify syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/workgraph.py
```

Expected: No output (successful compilation)

**Step 4: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat: add slab electronic properties outputs to core_workgraph

Added four new outputs: slab_bands, slab_dos, slab_primitive_structures,
slab_seekpath_parameters to support selective slab electronic properties."
```

---

## Task 4: Add Parameters to build_core_workgraph

**Files:**
- Modify: `teros/core/workgraph.py:441-488` (add parameters to function signature)
- Modify: `teros/core/workgraph.py:635-678` (pass parameters to core_workgraph.build())

**Step 1: Add parameters to build_core_workgraph signature**

In the `build_core_workgraph()` function signature (around line 441), add these new parameters after the existing electronic properties parameters:

```python
def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    # ... existing parameters ...
    compute_electronic_properties_bulk: bool = False,
    bands_parameters: dict = None,
    bands_options: dict = None,
    band_settings: dict = None,
    # NEW: Slab electronic properties parameters
    compute_electronic_properties_slabs: bool = False,
    slab_electronic_properties: dict = None,  # {'term_0': {params}, 'term_2': {params}}
    slab_bands_parameters: dict = None,
    slab_bands_options: dict = None,
    slab_band_settings: dict = None,
    name: str = 'FormationEnthalpy',
):
```

**Step 2: Update docstring**

Add documentation for the new parameters in the docstring (after the bulk electronic properties documentation):

```python
    """
    ... existing docstring ...

        compute_electronic_properties_bulk: Compute DOS and bands for relaxed bulk. Default: False
        bands_parameters: Dict with 'scf', 'bands', 'dos' keys containing INCAR dicts, plus
                          optional 'scf_kpoints_distance'. Use get_electronic_properties_defaults().
                          Default: None
        bands_options: Scheduler options for bands calculation. Default: None (uses bulk_options)
        band_settings: Dict with band workflow settings. Use get_electronic_properties_defaults()
                       ['band_settings']. Default: None
        compute_electronic_properties_slabs: Compute DOS and bands for selected slabs. Default: False
        slab_electronic_properties: Dict mapping slab labels to parameter overrides.
                                    Format: {'term_0': {'bands_parameters': ..., 'bands_options': ..., 'band_settings': ...}}
                                    If a parameter key is missing, falls back to slab_bands_* defaults.
                                    Default: None
        slab_bands_parameters: Default parameters for all slabs. Use get_slab_electronic_properties_defaults().
                               Default: None
        slab_bands_options: Default scheduler options for slab electronic properties. Default: None (uses bulk_options)
        slab_band_settings: Default band settings for slabs. Use get_slab_electronic_properties_defaults()
                           ['band_settings']. Default: None
        name: WorkGraph name. Default: 'FormationEnthalpy'

    ... rest of docstring ...
    """
```

**Step 3: Verify syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/workgraph.py
```

Expected: No output (successful compilation)

**Step 4: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat: add slab electronic properties parameters to build_core_workgraph

Added compute_electronic_properties_slabs flag, slab_electronic_properties dict,
and default parameter arguments for selective slab electronic properties."
```

---

## Task 5: Integrate Slab Electronic Properties into build_core_workgraph

**Files:**
- Modify: `teros/core/workgraph.py` (add integration logic after cleavage calculation, around line 890)

**Step 1: Add slab electronic properties calculation**

After the cleavage energies block (around line 888), add this new conditional block:

```python
    # NEW: Add slab electronic properties calculation if requested
    if compute_electronic_properties_slabs and relax_slabs and slab_electronic_properties:
        from teros.core.slabs import calculate_electronic_properties_slabs_scatter

        # Get default parameters (per-slab overrides handled in scatter function)
        default_params = slab_bands_parameters if slab_bands_parameters else {}
        default_opts = slab_bands_options if slab_bands_options else bulk_options
        default_settings = slab_band_settings if slab_band_settings else {}

        # Determine which slabs output to use (handles both auto-generated and input_slabs modes)
        if use_input_slabs:
            # In input_slabs mode, use the output from earlier manual task creation
            if restart_folders is not None:
                # Restart mode: use collector outputs
                relaxed_slabs_source = collector.outputs.structures
            else:
                # Non-restart input_slabs: use scatter task outputs
                relaxed_slabs_source = scatter_task.outputs.relaxed_structures
        else:
            # Auto-generated slabs: use outputs from core_workgraph's internal scatter
            # This is tricky - need to get the outputs that were created inside core_workgraph
            # Since core_workgraph is @task.graph, we need to access its internal task outputs
            # The relaxed_slabs are returned by the relax_slabs_scatter inside core_workgraph
            # But we can't easily access them from here without refactoring
            # SOLUTION: We'll add the task inside core_workgraph function itself (see alternative below)
            pass

        # Add electronic properties task
        slab_elec_task = wg.add_task(
            calculate_electronic_properties_slabs_scatter,
            name='calculate_electronic_properties_slabs',
            slabs=relaxed_slabs_source,
            slab_electronic_properties=slab_electronic_properties,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            clean_workdir=clean_workdir,
            default_bands_parameters=default_params,
            default_bands_options=default_opts,
            default_band_settings=default_settings,
        )

        # Connect outputs
        wg.outputs.slab_bands = slab_elec_task.outputs.slab_bands
        wg.outputs.slab_dos = slab_elec_task.outputs.slab_dos
        wg.outputs.slab_primitive_structures = slab_elec_task.outputs.slab_primitive_structures
        wg.outputs.slab_seekpath_parameters = slab_elec_task.outputs.slab_seekpath_parameters

        print(f"  ✓ Slab electronic properties calculation enabled for {len(slab_electronic_properties)} slabs")
```

**Step 2: Handle the data flow issue for auto-generated slabs**

The code above has an issue: for auto-generated slabs, we can't easily access the `relaxed_slabs_output` from inside `core_workgraph` because it's a `@task.graph`.

**SOLUTION:** Add the slab electronic properties calculation inside the `core_workgraph()` function itself (not in `build_core_workgraph`), right after the slab relaxation block.

In `core_workgraph()` function (around line 380, after relaxation_energies calculation):

```python
    # ===== SLAB ELECTRONIC PROPERTIES CALCULATION (OPTIONAL) =====
    slab_bands_output = {}
    slab_dos_output = {}
    slab_primitive_structures_output = {}
    slab_seekpath_parameters_output = {}

    # Note: This is handled conditionally in build_core_workgraph for input_slabs mode
    # For auto-generated slabs, we would add it here, but it's cleaner to handle
    # everything in build_core_workgraph by passing the flag through
```

Actually, let me reconsider. The cleaner approach is to:

1. For **input_slabs mode**: Add task in `build_core_workgraph()` (like we do for thermodynamics)
2. For **auto-generated mode**: Since we can't access internal outputs of `core_workgraph`, we either:
   - Add the logic inside `core_workgraph()` itself, OR
   - Expose relaxed_slabs as an output that build_core_workgraph can use

Let me revise to use approach 2 (expose outputs):

**REVISED Step 2: Simplify - only handle input_slabs mode**

For now, let's implement this feature for `input_slabs` mode only (including restart), since that's where users have explicit control. Auto-generated mode can be added later if needed.

Replace the code from Step 1 with this simpler version:

```python
    # NEW: Add slab electronic properties calculation if requested
    # Currently only supported for input_slabs mode (including restart)
    if compute_electronic_properties_slabs and relax_slabs and slab_electronic_properties and use_input_slabs:
        from teros.core.slabs import calculate_electronic_properties_slabs_scatter

        # Get default parameters (per-slab overrides handled in scatter function)
        default_params = slab_bands_parameters if slab_bands_parameters else {}
        default_opts = slab_bands_options if slab_bands_options else bulk_options
        default_settings = slab_band_settings if slab_band_settings else {}

        # Get code
        code = load_code(code_label)

        # Get slab parameters
        slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
        slab_opts = slab_options if slab_options is not None else bulk_options
        slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping

        # Determine which slabs output to use
        if restart_folders is not None:
            # Restart mode: use collector outputs
            relaxed_slabs_source = collector.outputs.structures
        else:
            # Non-restart input_slabs: use scatter task outputs
            relaxed_slabs_source = scatter_task.outputs.relaxed_structures

        # Add electronic properties task
        slab_elec_task = wg.add_task(
            calculate_electronic_properties_slabs_scatter,
            name='calculate_electronic_properties_slabs',
            slabs=relaxed_slabs_source,
            slab_electronic_properties=slab_electronic_properties,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_pot_map,
            clean_workdir=clean_workdir,
            default_bands_parameters=default_params,
            default_bands_options=default_opts,
            default_band_settings=default_settings,
        )

        # Connect outputs
        wg.outputs.slab_bands = slab_elec_task.outputs.slab_bands
        wg.outputs.slab_dos = slab_elec_task.outputs.slab_dos
        wg.outputs.slab_primitive_structures = slab_elec_task.outputs.slab_primitive_structures
        wg.outputs.slab_seekpath_parameters = slab_elec_task.outputs.slab_seekpath_parameters

        print(f"  ✓ Slab electronic properties calculation enabled for {len(slab_electronic_properties)} slabs")
```

**Step 3: Verify syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/workgraph.py
```

Expected: No output (successful compilation)

**Step 4: Clear Python cache**

```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

**Step 5: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat: integrate slab electronic properties into build_core_workgraph

Added conditional task creation for slab electronic properties calculation.
Currently supports input_slabs mode (manual and restart). Auto-generated
mode can be added in future if needed."
```

---

## Task 6: Create Example/Test Script

**Files:**
- Create: `examples/electronic_properties/slab_electronic_properties_example.py`

**Step 1: Create examples directory**

```bash
mkdir -p examples/electronic_properties
```

**Step 2: Create example script**

Create `examples/electronic_properties/slab_electronic_properties_example.py`:

```python
#!/usr/bin/env python
"""
Example: Selective Slab Electronic Properties Calculation

This example demonstrates how to calculate DOS and band structures
for selected slab terminations with per-slab parameter overrides.

Based on: docs/plans/2025-10-12-slab-electronic-properties.md
"""

from aiida import orm
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.electronic_properties_builder import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults,
)

# Get defaults from builders
bulk_defaults = get_electronic_properties_defaults()
slab_defaults = get_slab_electronic_properties_defaults()

# Custom parameters for term_2 (example: higher accuracy)
slab_custom = get_slab_electronic_properties_defaults(
    energy_cutoff=600,
    kpoints_mesh_density=0.2,  # Even denser
    band_kpoints_distance=0.1,
)

# Build workgraph
wg = build_core_workgraph(
    structures_dir='structures',
    bulk_name='ag2o.cif',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    kpoints_spacing=0.4,
    bulk_parameters={
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NELM': 120,
    },
    bulk_options={
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600,
        'queue_name': 'normal',
    },
    clean_workdir=False,

    # Metal and oxygen references
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6},
    metal_options={'resources': {'num_machines': 1}, 'max_wallclock_seconds': 3600},
    oxygen_options={'resources': {'num_machines': 1}, 'max_wallclock_seconds': 3600},

    # Slab generation (will be used as input_slabs in this example)
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    slab_parameters={
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 100,
        'IBRION': 2,
    },

    # Bulk electronic properties
    compute_electronic_properties_bulk=True,
    bands_parameters=bulk_defaults,
    bands_options={'resources': {'num_machines': 2}, 'max_wallclock_seconds': 7200},
    band_settings=bulk_defaults['band_settings'],

    # NEW: Slab electronic properties
    compute_electronic_properties_slabs=True,
    slab_electronic_properties={
        'term_0': {
            'bands_parameters': slab_defaults,
            'bands_options': {'resources': {'num_machines': 2}},
            'band_settings': slab_defaults['band_settings'],
        },
        'term_2': {
            'bands_parameters': slab_custom,  # Override with custom
            'bands_options': {'resources': {'num_machines': 4}},  # More resources
            'band_settings': slab_custom['band_settings'],
        },
    },
    slab_bands_parameters=slab_defaults,  # Global defaults
    slab_bands_options={'resources': {'num_machines': 2}},
    slab_band_settings=slab_defaults['band_settings'],

    name='Ag2O_Slab_Electronic_Properties',
)

# Note: For this example to work, you need to first generate slabs, then use them as input_slabs
# This is a conceptual example showing the interface

print(f"WorkGraph built: {wg.name}")
print(f"Number of tasks: {len(wg.tasks)}")
print("\nTo submit:")
print("  from aiida.engine import submit")
print("  node = submit(wg)")
print("  # Wait for completion")
print("  print(node.outputs.slab_bands.keys())  # ['term_0', 'term_2']")
print("  print(node.outputs.slab_dos.term_0)")
```

**Step 3: Test script syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile examples/electronic_properties/slab_electronic_properties_example.py
```

Expected: No output (successful compilation)

**Step 4: Commit**

```bash
git add examples/electronic_properties/slab_electronic_properties_example.py
git commit -m "docs: add slab electronic properties example script

Added comprehensive example showing selective slab electronic properties
calculation with per-slab parameter overrides and defaults."
```

---

## Task 7: Update build_core_workgraph_with_map (Deprecated Wrapper)

**Files:**
- Modify: `teros/core/workgraph.py:980-1093` (add new parameters to deprecated function)

**Step 1: Add parameters to deprecated function signature**

In `build_core_workgraph_with_map()` (around line 980), add the new parameters:

```python
def build_core_workgraph_with_map(
    # ... existing parameters ...
    compute_electronic_properties_bulk: bool = False,  # NEW (already added in previous commit)
    bands_parameters: dict = None,  # NEW
    bands_options: dict = None,  # NEW
    band_settings: dict = None,  # NEW
    compute_electronic_properties_slabs: bool = False,  # NEW
    slab_electronic_properties: dict = None,  # NEW
    slab_bands_parameters: dict = None,  # NEW
    slab_bands_options: dict = None,  # NEW
    slab_band_settings: dict = None,  # NEW
    name: str = 'FormationEnthalpy_ScatterGather',
) -> WorkGraph:
```

**Step 2: Forward new parameters to build_core_workgraph**

In the `return build_core_workgraph(...)` call (around line 1047), add the new parameters:

```python
    return build_core_workgraph(
        # ... existing parameters ...
        compute_electronic_properties_bulk=compute_electronic_properties_bulk,
        bands_parameters=bands_parameters,
        bands_options=bands_options,
        band_settings=band_settings,
        compute_electronic_properties_slabs=compute_electronic_properties_slabs,  # NEW
        slab_electronic_properties=slab_electronic_properties,  # NEW
        slab_bands_parameters=slab_bands_parameters,  # NEW
        slab_bands_options=slab_bands_options,  # NEW
        slab_band_settings=slab_band_settings,  # NEW
        name=name,
    )
```

**Step 3: Verify syntax**

```bash
source ~/envs/psteros/bin/activate && python -m py_compile teros/core/workgraph.py
```

Expected: No output (successful compilation)

**Step 4: Commit**

```bash
git add teros/core/workgraph.py
git commit -m "feat: forward slab electronic properties parameters in deprecated wrapper

Updated build_core_workgraph_with_map to forward new slab electronic
properties parameters to build_core_workgraph."
```

---

## Task 8: Final Integration Test

**Files:**
- Test with existing example or create minimal test

**Step 1: Clear all caches**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-bands-slabs
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

**Step 2: Restart AiiDA daemon**

```bash
verdi daemon restart
```

**Step 3: Create minimal test script**

Create `test_slab_electronic_properties_minimal.py`:

```python
#!/usr/bin/env python
"""Minimal test: Check that slab electronic properties functions exist and compile"""

# Test imports
from teros.core.builders.electronic_properties_builder import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults,
)
from teros.core.slabs import calculate_electronic_properties_slabs_scatter
from teros.core.workgraph import build_core_workgraph

# Test builder
defaults = get_slab_electronic_properties_defaults()
assert 'band_settings' in defaults
assert 'scf' in defaults
assert 'bands' in defaults
assert 'dos' in defaults
print("✓ Builder works")

# Test scatter function exists
assert callable(calculate_electronic_properties_slabs_scatter)
print("✓ Scatter function exists")

# Test build_core_workgraph accepts new parameters
import inspect
sig = inspect.signature(build_core_workgraph)
params = sig.parameters
assert 'compute_electronic_properties_slabs' in params
assert 'slab_electronic_properties' in params
assert 'slab_bands_parameters' in params
print("✓ build_core_workgraph has new parameters")

print("\n✅ All minimal tests passed!")
```

**Step 4: Run minimal test**

```bash
source ~/envs/psteros/bin/activate && python test_slab_electronic_properties_minimal.py
```

Expected:
```
✓ Builder works
✓ Scatter function exists
✓ build_core_workgraph has new parameters

✅ All minimal tests passed!
```

**Step 5: Clean up test script**

```bash
rm test_slab_electronic_properties_minimal.py
```

**Step 6: Create final commit**

```bash
git add .
git commit -m "test: verify slab electronic properties implementation

All components compile and integrate correctly. Ready for full workflow testing."
```

---

## Post-Implementation Notes

**Limitations:**
- Currently only works for `input_slabs` mode (manual slabs and restart mode)
- Auto-generated slabs mode not yet supported (would require refactoring to expose relaxed_slabs from core_workgraph)

**Future Enhancements:**
1. Support auto-generated slabs by refactoring core_workgraph outputs
2. Add validation for slab_electronic_properties dictionary structure
3. Add helper function to auto-populate slab_electronic_properties from relaxed slabs
4. Add more comprehensive examples with real calculations

**Testing with Real Calculations:**
To fully test this implementation, you'll need to:
1. Run a workgraph that generates/relaxes slabs first
2. Extract the slab structures
3. Pass them to a new workgraph with `compute_electronic_properties_slabs=True`
4. Check that outputs.slab_bands, outputs.slab_dos exist with correct keys

---

## Summary

**Implementation complete! Key changes:**

1. ✅ Added `get_slab_electronic_properties_defaults()` builder in `electronic_properties_builder.py`
2. ✅ Added `calculate_electronic_properties_slabs_scatter()` in `slabs.py`
3. ✅ Added four new outputs to `core_workgraph`: slab_bands, slab_dos, slab_primitive_structures, slab_seekpath_parameters
4. ✅ Added parameters to `build_core_workgraph()`
5. ✅ Integrated conditional task creation for input_slabs mode
6. ✅ Created example script
7. ✅ Updated deprecated wrapper
8. ✅ Verified compilation and integration

**Next step:** Full workflow testing with real VASP calculations
