# Bulk Electronic Properties (DOS and Bands) Implementation Plan

> **For Claude:** Use `${CLAUDE_PLUGIN_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Add optional DOS and band structure calculation for relaxed bulk structures using AiiDA-VASP's vasp.v2.bands workchain.

**Architecture:** Create a material-agnostic builder for electronic properties parameters, add conditional task in core_workgraph that uses vasp.v2.bands with relaxed bulk structure as input. The bands workchain handles SCF→bands→DOS chain internally.

**Tech Stack:** AiiDA-VASP vasp.v2.bands workchain, AiiDA WorkGraph, Python

---

## Task 1: Create Electronic Properties Builder

**Files:**
- Create: `teros/core/builders/electronic_properties_builder.py`
- Modify: `teros/core/builders/__init__.py`

**Step 1: Create the builder file with function signature**

Create `teros/core/builders/electronic_properties_builder.py`:

```python
"""
Electronic Properties Builder

Material-agnostic builder for DOS and band structure calculations
using AiiDA-VASP's vasp.v2.bands workchain.
"""

def get_electronic_properties_defaults(
    energy_cutoff: float = 500,
    electronic_convergence: float = 1e-5,
    ncore: int = 4,
    ispin: int = 2,
    lasph: bool = True,
    lreal: str = "Auto",
    kpoints_mesh_density: float = 0.3,
    band_kpoints_distance: float = 0.2,
    dos_kpoints_distance: float = 0.2,
    line_density: float = 0.2,
    nedos: int = 2000,
    sigma_bands: float = 0.01,
    symprec: float = 1e-4,
    band_mode: str = "seekpath-aiida",
) -> dict:
    """
    Get default parameters for electronic properties (DOS and bands) calculation.

    Material-agnostic builder that returns all necessary parameters for
    vasp.v2.bands workchain. Based on proven working configuration.

    Args:
        energy_cutoff: ENCUT value (eV). Default: 500
        electronic_convergence: EDIFF value. Default: 1e-5
        ncore: Number of cores per band. Default: 4
        ispin: Spin polarization (1=non-magnetic, 2=spin-polarized). Default: 2
        lasph: Include aspherical contributions. Default: True
        lreal: Projection operators (False for small, "Auto" for large). Default: "Auto"
        kpoints_mesh_density: K-point density for SCF and DOS meshes. Default: 0.3
        band_kpoints_distance: K-point spacing for band structure path. Default: 0.2
        dos_kpoints_distance: K-point spacing for DOS calculation. Default: 0.2
        line_density: K-point density along high-symmetry paths. Default: 0.2
        nedos: Number of DOS grid points. Default: 2000
        sigma_bands: Smearing width for bands (eV). Default: 0.01
        symprec: Symmetry precision for seekpath. Default: 1e-4
        band_mode: Band path mode (seekpath-aiida, pymatgen, etc.). Default: "seekpath-aiida"

    Returns:
        Dictionary with 'band_settings' and nested INCAR parameters for
        scf, bands, and dos stages, ready to use with vasp.v2.bands workchain.

    Example:
        >>> defaults = get_electronic_properties_defaults()
        >>> wg = build_core_workgraph(
        ...     compute_electronic_properties_bulk=True,
        ...     bands_parameters=defaults,
        ...     band_settings=defaults['band_settings'],
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
        'run_dos': True,  # Always compute DOS when computing bands
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
        'SYSTEM': 'SCF_for_Electronic_Properties',
        'ISTART': 0,
        'ICHARG': 2,
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
        'SYSTEM': 'Band_Structure',
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
        'SYSTEM': 'DOS_Calculation',
        'ISTART': 1,     # Read from WAVECAR
        'ICHARG': 11,    # Non-self-consistent calculation
        'PREC': 'Accurate',
        'ALGO': 'Normal',
        'NELM': 60,
        'ISMEAR': -5,    # Tetrahedron method with Blöchl corrections
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

**Step 2: Update __init__.py to expose the builder**

Modify `teros/core/builders/__init__.py`:

```python
"""
PS-TEROS Builders Module

Contains default parameter builders for various material systems.
"""

from .default_ag2o_builders import get_ag2o_defaults
from .default_ag3po4_builders import get_ag3po4_defaults
from .electronic_properties_builder import get_electronic_properties_defaults

__all__ = [
    'get_ag2o_defaults',
    'get_ag3po4_defaults',
    'get_electronic_properties_defaults',
]
```

**Step 3: Verify the builder can be imported**

Run:
```bash
source ~/envs/psteros/bin/activate && python -c "from teros.core.builders import get_electronic_properties_defaults; print(get_electronic_properties_defaults()['band_settings'])"
```

Expected: Should print the band_settings dictionary without errors.

**Step 4: Commit the builder**

```bash
git add teros/core/builders/electronic_properties_builder.py teros/core/builders/__init__.py
git commit -m "feat: add electronic properties builder for DOS and bands"
```

---

## Task 2: Add Electronic Properties Task to core_workgraph

**Files:**
- Modify: `teros/core/workgraph.py:58-65` (add outputs to @task.graph decorator)
- Modify: `teros/core/workgraph.py:66-108` (add parameters to core_workgraph function)
- Modify: `teros/core/workgraph.py:403-423` (add electronic properties calculation and outputs)

**Step 1: Add outputs to @task.graph decorator**

In `teros/core/workgraph.py`, modify the `@task.graph` decorator (around line 58):

```python
@task.graph(outputs=[
    'bulk_energy', 'metal_energy', 'nonmetal_energy', 'oxygen_energy',
    'bulk_structure', 'metal_structure', 'nonmetal_structure', 'oxygen_structure',
    'formation_enthalpy', 'slab_structures', 'relaxed_slabs', 'slab_energies',
    'unrelaxed_slab_energies', 'relaxation_energies',
    'surface_energies', 'cleavage_energies',
    'slab_remote', 'unrelaxed_slab_remote',
    'bulk_bands', 'bulk_dos', 'bulk_electronic_properties_misc',  # NEW: Electronic properties outputs
])
```

**Step 2: Add parameters to core_workgraph function signature**

In `teros/core/workgraph.py`, add parameters to `core_workgraph` function (around line 100):

```python
def core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
    bulk_parameters: dict,
    bulk_options: dict,
    clean_workdir: bool,
    metal_name: str = None,
    oxygen_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    miller_indices: list = None,
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,
    input_slabs: dict = None,
    use_input_slabs: bool = False,
    compute_cleavage: bool = True,
    compute_electronic_properties_bulk: bool = False,  # NEW: Enable DOS and bands
    bands_parameters: dict = None,  # NEW: INCAR parameters for scf/bands/dos stages
    bands_options: dict = None,  # NEW: Scheduler options for bands calculation
    band_settings: dict = None,  # NEW: Band workflow settings (symprec, band_mode, etc.)
):
```

**Step 3: Add electronic properties calculation after bulk relaxation**

In `teros/core/workgraph.py`, add the electronic properties calculation after the bulk relaxation section and before slab generation (around line 291, just before the SLAB GENERATION comment):

```python
    # ===== ELECTRONIC PROPERTIES CALCULATION (OPTIONAL) =====
    bulk_bands_output = None
    bulk_dos_output = None
    bulk_electronic_misc_output = {}

    if compute_electronic_properties_bulk:
        # Get BandsWorkChain and wrap it as a task
        BandsWorkChain = WorkflowFactory('vasp.v2.bands')
        BandsTask = task(BandsWorkChain)

        # Use bands-specific options or fall back to bulk options
        bands_opts = bands_options if bands_options is not None else bulk_options

        # Build bands workchain task - uses relaxed bulk structure
        # Note: vasp.v2.bands handles SCF → bands → DOS chain internally
        bands_calc = BandsTask(
            structure=bulk_vasp.structure,  # Use relaxed structure from bulk calculation
            code=code,
            potential_family=potential_family,
            potential_mapping=bulk_potential_mapping,
            options=bands_opts,
            clean_workdir=clean_workdir,
        )

        # Extract outputs
        bulk_bands_output = bands_calc.outputs.bands
        bulk_dos_output = bands_calc.outputs.dos
        bulk_electronic_misc_output = bands_calc.outputs.misc
```

**Step 4: Add electronic properties to return dictionary**

In `teros/core/workgraph.py`, modify the return statement (around line 403):

```python
    # Return all outputs
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
        'bulk_bands': bulk_bands_output,  # NEW: Band structure data
        'bulk_dos': bulk_dos_output,  # NEW: DOS data
        'bulk_electronic_properties_misc': bulk_electronic_misc_output,  # NEW: Misc outputs
    }
```

**Step 5: Update docstring for core_workgraph**

In `teros/core/workgraph.py`, update the docstring (around line 109):

Add to the "Calculations controlled by boolean flags" section:
```
        - compute_electronic_properties_bulk: Calculate DOS and band structure for relaxed bulk
```

Add to Args section:
```
        compute_electronic_properties_bulk: Whether to compute DOS and bands for bulk. Default: False
        bands_parameters: Dictionary with 'scf', 'bands', 'dos' keys containing INCAR dicts.
                          Use get_electronic_properties_defaults() builder. Default: None
        bands_options: Scheduler options for bands calculation. Default: None (uses bulk_options)
        band_settings: Dict with band workflow settings (symprec, band_mode, kpoints_distance, etc.).
                       Use get_electronic_properties_defaults()['band_settings']. Default: None
```

Add to Returns section:
```
            - bulk_bands: BandsData with electronic band structure (if compute_electronic_properties_bulk=True)
            - bulk_dos: ArrayData with density of states (if compute_electronic_properties_bulk=True)
            - bulk_electronic_properties_misc: Dict with additional outputs (if compute_electronic_properties_bulk=True)
```

**Step 6: Commit core_workgraph modifications**

```bash
git add teros/core/workgraph.py
git commit -m "feat: add electronic properties calculation to core_workgraph"
```

---

## Task 3: Add Parameter Handling to build_core_workgraph

**Files:**
- Modify: `teros/core/workgraph.py:426-468` (add parameters to build_core_workgraph signature)
- Modify: `teros/core/workgraph.py:862-865` (add electronic properties task when needed)

**Step 1: Add parameters to build_core_workgraph function signature**

In `teros/core/workgraph.py`, modify `build_core_workgraph` function signature (around line 426):

```python
def build_core_workgraph(
    structures_dir: str,
    bulk_name: str,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,
    metal_name: str = None,
    oxygen_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    miller_indices: list = None,
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,
    input_slabs: dict = None,
    compute_cleavage: bool = True,
    restart_from_node: int = None,
    compute_electronic_properties_bulk: bool = False,  # NEW: Enable DOS and bands
    bands_parameters: dict = None,  # NEW: INCAR parameters
    bands_options: dict = None,  # NEW: Scheduler options
    band_settings: dict = None,  # NEW: Band workflow settings
    name: str = 'FormationEnthalpy',
):
```

**Step 2: Pass electronic properties parameters to core_workgraph.build()**

In `teros/core/workgraph.py`, modify the `core_workgraph.build()` call (around line 608):

```python
    # Build the workgraph (pass None for input_slabs to avoid serialization issues)
    # Note: parameters will be wrapped with {'incar': ...} inside the graph
    wg = core_workgraph.build(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping or {},
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters or {},
        bulk_options=bulk_options or {},
        clean_workdir=clean_workdir,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        metal_potential_mapping=metal_potential_mapping,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        oxygen_potential_mapping=oxygen_potential_mapping,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        nonmetal_name=nonmetal_name,
        nonmetal_potential_mapping=nonmetal_potential_mapping,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        miller_indices=miller_indices if not use_input_slabs else [0, 0, 1],  # Dummy value
        min_slab_thickness=min_slab_thickness if not use_input_slabs else 10.0,  # Dummy value
        min_vacuum_thickness=min_vacuum_thickness if not use_input_slabs else 10.0,  # Dummy value
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        relax_slabs=relax_slabs,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        compute_relaxation_energy=compute_relaxation_energy,
        input_slabs=None,  # Always pass None to avoid serialization
        use_input_slabs=use_input_slabs,  # Pass the flag
        compute_cleavage=compute_cleavage,
        compute_electronic_properties_bulk=False,  # NEW: Will be handled manually below
        bands_parameters=None,  # NEW: Will be handled manually below
        bands_options=None,  # NEW: Will be handled manually below
        band_settings=None,  # NEW: Will be handled manually below
    )
```

**Step 3: Add manual electronic properties task when enabled**

In `teros/core/workgraph.py`, add electronic properties handling after the slab handling section (around line 861, just before setting wg.name):

```python
    # NEW: Add electronic properties calculation if requested
    if compute_electronic_properties_bulk:
        from aiida.plugins import WorkflowFactory

        BandsWC = WorkflowFactory('vasp.v2.bands')

        # Get bulk task from workgraph
        bulk_vasp_task = wg.tasks['VaspWorkChain']

        # Get code
        code = load_code(code_label)

        # Use bands-specific options or fall back to bulk options
        bands_opts = bands_options if bands_options is not None else bulk_options

        # Prepare basic inputs
        bands_inputs = {
            'structure': bulk_vasp_task.outputs.structure,
            'code': code,
            'potential_family': potential_family,
            'potential_mapping': bulk_potential_mapping,
            'options': bands_opts,
            'clean_workdir': clean_workdir,
        }

        # Add band_settings if provided
        if band_settings:
            bands_inputs['band_settings'] = orm.Dict(dict=band_settings)

        # Add nested namespace parameters (scf, bands, dos)
        # bands_parameters should be a dict with keys: 'scf', 'bands', 'dos'
        if bands_parameters:
            for namespace in ['scf', 'bands', 'dos']:
                if namespace in bands_parameters:
                    # Use dot notation for nested namespace
                    bands_inputs[f'{namespace}.parameters'] = orm.Dict(
                        dict={'incar': bands_parameters[namespace]}
                    )

            # Add SCF k-points if provided
            if 'scf_kpoints_distance' in bands_parameters:
                from aiida.orm import KpointsData
                kpoints = KpointsData()
                kpoints.set_cell_from_structure(bulk_vasp_task.outputs.structure)
                kpoints.set_kpoints_mesh_from_density(
                    bands_parameters['scf_kpoints_distance']
                )
                bands_inputs['scf.kpoints'] = kpoints

        # Add the bands task
        bands_task = wg.add_task(
            BandsWC,
            name='BandsWorkChain_bulk',
            **bands_inputs
        )

        # Connect outputs
        wg.outputs.bulk_bands = bands_task.outputs.bands
        wg.outputs.bulk_dos = bands_task.outputs.dos
        wg.outputs.bulk_electronic_properties_misc = bands_task.outputs.misc

        print(f"  ✓ Electronic properties calculation enabled (DOS and bands)")
```

**Step 4: Update build_core_workgraph docstring**

In `teros/core/workgraph.py`, update the docstring (around line 470):

Add to main description:
```
    - Electronic properties calculation (controlled by compute_electronic_properties_bulk flag, default: False)
```

Add to Args section:
```
        compute_electronic_properties_bulk: Compute DOS and bands for relaxed bulk. Default: False
        bands_parameters: Dict with 'scf', 'bands', 'dos' keys containing INCAR dicts, plus
                          optional 'scf_kpoints_distance'. Use get_electronic_properties_defaults().
                          Default: None
        bands_options: Scheduler options for bands calculation. Default: None (uses bulk_options)
        band_settings: Dict with band workflow settings. Use get_electronic_properties_defaults()
                       ['band_settings']. Default: None
```

**Step 5: Commit build_core_workgraph modifications**

```bash
git add teros/core/workgraph.py
git commit -m "feat: add electronic properties parameter handling to build_core_workgraph"
```

---

## Task 4: Update build_core_workgraph_with_map Wrapper

**Files:**
- Modify: `teros/core/workgraph.py:868-972` (add parameters to wrapper function)

**Step 1: Add parameters to build_core_workgraph_with_map signature**

In `teros/core/workgraph.py`, modify `build_core_workgraph_with_map` function signature (around line 868):

```python
def build_core_workgraph_with_map(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    oxygen_name: str,
    miller_indices: list = None,
    min_slab_thickness: float = 18.0,
    min_vacuum_thickness: float = 15.0,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    metal_potential_mapping: dict = None,
    oxygen_potential_mapping: dict = None,
    nonmetal_name: str = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    clean_workdir: bool = False,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = False,
    center_slab: bool = True,
    symmetrize: bool = False,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = False,
    compute_thermodynamics: bool = False,
    thermodynamics_sampling: int = 100,
    input_slabs: dict = None,
    compute_cleavage: bool = True,
    compute_electronic_properties_bulk: bool = False,  # NEW
    bands_parameters: dict = None,  # NEW
    bands_options: dict = None,  # NEW
    band_settings: dict = None,  # NEW
    name: str = 'FormationEnthalpy_ScatterGather',
) -> WorkGraph:
```

**Step 2: Forward new parameters to build_core_workgraph**

In `teros/core/workgraph.py`, modify the return statement in `build_core_workgraph_with_map` (around line 931):

```python
    # Simply forward to build_core_workgraph which now uses scatter-gather
    return build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        metal_name=metal_name,
        oxygen_name=oxygen_name,
        nonmetal_name=nonmetal_name,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping,
        metal_potential_mapping=metal_potential_mapping,
        nonmetal_potential_mapping=nonmetal_potential_mapping,
        oxygen_potential_mapping=oxygen_potential_mapping,
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=clean_workdir,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        relax_slabs=relax_slabs,
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        input_slabs=input_slabs,
        compute_cleavage=compute_cleavage,
        compute_electronic_properties_bulk=compute_electronic_properties_bulk,  # NEW
        bands_parameters=bands_parameters,  # NEW
        bands_options=bands_options,  # NEW
        band_settings=band_settings,  # NEW
        name=name,
    )
```

**Step 3: Commit wrapper modifications**

```bash
git add teros/core/workgraph.py
git commit -m "feat: forward electronic properties parameters in build_core_workgraph_with_map"
```

---

## Task 5: Clear Python Cache and Restart Daemon

**Step 1: Clear Python cache**

Run:
```bash
find /home/thiagotd/git/PS-TEROS/.worktree/feature-dos-bands -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null ; find /home/thiagotd/git/PS-TEROS/.worktree/feature-dos-bands -name "*.pyc" -delete 2>/dev/null ; echo "Cache cleared"
```

Expected: "Cache cleared"

**Step 2: Restart AiiDA daemon**

Run:
```bash
verdi daemon restart
```

Expected: Daemon restarts successfully

**Step 3: Verify daemon status**

Run:
```bash
verdi daemon status
```

Expected: Daemon is running

**Step 4: Test import of modified modules**

Run:
```bash
source ~/envs/psteros/bin/activate && python -c "
from teros.core.workgraph import build_core_workgraph
from teros.core.builders import get_electronic_properties_defaults
print('✓ Imports successful')
defaults = get_electronic_properties_defaults()
print('✓ Builder returns:', list(defaults.keys()))
"
```

Expected: Should print imports successful and show keys: band_settings, scf, bands, dos, scf_kpoints_distance

**Step 5: Final commit**

```bash
git add -A
git commit -m "chore: clear cache and verify electronic properties implementation"
```

---

## Summary

This implementation adds optional DOS and band structure calculation for relaxed bulk structures in PS-TEROS:

1. **Builder**: `get_electronic_properties_defaults()` provides material-agnostic defaults based on proven working parameters
2. **Core workflow**: `compute_electronic_properties_bulk=False` flag enables the calculation
3. **Integration**: Uses `vasp.v2.bands` workchain which handles SCF→bands→DOS chain internally
4. **Outputs**: `bulk_bands`, `bulk_dos`, `bulk_electronic_properties_misc` available after workflow completion

**Next Steps:**
- Create an example script in `examples/` to demonstrate usage
- Test with actual structure files
- Validate outputs with AiiDA tools

**Testing Commands:**
```bash
# After creating example script
verdi daemon restart
source ~/envs/psteros/bin/activate && python examples/electronic_properties/electronic_properties_example.py
sleep 30
verdi process show <PK>
verdi process report <PK>
```
