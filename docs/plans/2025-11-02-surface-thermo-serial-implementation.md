# Surface Thermodynamics Serial Preset Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a serial/flat variant of surface_thermodynamics preset where all VASP nodes exist at the same graph level, allowing max_concurrent_jobs to control concurrent execution.

**Architecture:** Remove @task.graph decorators from scatter operations and replace with regular Python functions that directly add nodes to the main WorkGraph. All VASP calculation nodes will be siblings in the graph hierarchy, not nested in sub-workgraphs.

**Tech Stack:** AiiDA, AiiDA-WorkGraph, VASP, Python 3.x

---

## Task 1: Create Module Structure

**Files:**
- Create: `teros/experimental/surface_thermo_preset_serial/__init__.py`
- Create: `teros/experimental/surface_thermo_preset_serial/workgraph.py`
- Create: `teros/experimental/surface_thermo_preset_serial/slab_operations.py`
- Create: `teros/experimental/surface_thermo_preset_serial/thermodynamics_operations.py`
- Create: `teros/experimental/surface_thermo_preset_serial/utils.py`

**Step 1: Create directory**

```bash
mkdir -p teros/experimental/surface_thermo_preset_serial
```

**Step 2: Create __init__.py with exports**

Create: `teros/experimental/surface_thermo_preset_serial/__init__.py`

```python
"""
Serial Surface Thermodynamics Preset

Experimental flat-graph implementation where all VASP nodes exist at the same
graph level, allowing max_concurrent_jobs to control concurrent execution.
"""

from .workgraph import surface_thermodynamics_serial_workgraph

__all__ = ['surface_thermodynamics_serial_workgraph']
```

**Step 3: Create empty module files**

Create empty files for now:
- `teros/experimental/surface_thermo_preset_serial/workgraph.py`
- `teros/experimental/surface_thermo_preset_serial/slab_operations.py`
- `teros/experimental/surface_thermo_preset_serial/thermodynamics_operations.py`
- `teros/experimental/surface_thermo_preset_serial/utils.py`

**Step 4: Commit module structure**

```bash
git add teros/experimental/surface_thermo_preset_serial/
git commit -m "feat: create surface_thermo_preset_serial module structure"
```

---

## Task 2: Implement Utilities Module

**Files:**
- Modify: `teros/experimental/surface_thermo_preset_serial/utils.py`

**Step 1: Implement parameter preparation utilities**

Write: `teros/experimental/surface_thermo_preset_serial/utils.py`

```python
"""Utility functions for the serial surface thermodynamics preset."""

from aiida import orm


def prepare_vasp_parameters(
    base_parameters: dict,
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    options: dict = None,
    clean_workdir: bool = False,
) -> dict:
    """
    Prepare standardized parameters for VASP WorkChain.

    Args:
        base_parameters: INCAR-like parameters dict
        code: VASP code to use
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        options: Computer options (num_machines, etc.)
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary ready for VaspWorkChain inputs
    """
    params = {
        'code': code,
        'parameters': orm.Dict(dict=base_parameters),
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict=potential_mapping),
        'options': orm.Dict(dict=options or {}),
        'kpoints_spacing': orm.Float(kpoints_spacing),
        'clean_workdir': orm.Bool(clean_workdir),
        'settings': orm.Dict(dict={
            'parser_settings': {
                'add_energy': True,
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }),
    }
    return params


def create_default_bulk_parameters() -> dict:
    """Create default VASP parameters for bulk relaxation."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'LWAVE': False,
        'LCHARG': False,
    }


def create_default_slab_parameters() -> dict:
    """Create default VASP parameters for slab calculations."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,  # Relax ions only, not cell
        'NSW': 100,
        'LWAVE': False,
        'LCHARG': False,
    }


def create_default_scf_parameters() -> dict:
    """Create default VASP parameters for SCF (single-point) calculations."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 0,  # No ionic relaxation
        'LWAVE': False,
        'LCHARG': False,
    }
```

**Step 2: Commit utilities**

```bash
git add teros/experimental/surface_thermo_preset_serial/utils.py
git commit -m "feat: add parameter preparation utilities"
```

---

## Task 3: Implement Slab Operations Module

**Files:**
- Modify: `teros/experimental/surface_thermo_preset_serial/slab_operations.py`

**Step 1: Implement node builder for slab relaxation**

Write: `teros/experimental/surface_thermo_preset_serial/slab_operations.py`

```python
"""Node builders for slab-related operations."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

VaspWorkChain = WorkflowFactory('vasp.vasp')


def build_scf_slabs_nodes(
    wg,
    slabs: dict[str, orm.StructureData],
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    parameters: dict,
    options: dict,
    clean_workdir: bool = False,
) -> dict[str, t.Any]:
    """
    Add VASP SCF (single-point) calculation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        slabs: Dictionary of {slab_id: StructureData}
        code: VASP code
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        parameters: VASP INCAR parameters
        options: Computer options
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary of {slab_id: vasp_task_node}
    """
    from .utils import prepare_vasp_parameters

    scf_nodes = {}

    for slab_id, slab_structure in slabs.items():
        # Prepare parameters
        vasp_params = prepare_vasp_parameters(
            base_parameters=parameters,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            kpoints_spacing=kpoints_spacing,
            options=options,
            clean_workdir=clean_workdir,
        )

        # Add VASP SCF node
        scf_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"scf_slab_{slab_id}",
            structure=slab_structure,
            **vasp_params,
        )

    return scf_nodes


def build_relax_slabs_nodes(
    wg,
    slabs: dict[str, orm.StructureData],
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    parameters: dict,
    options: dict,
    clean_workdir: bool = False,
) -> dict[str, t.Any]:
    """
    Add VASP relaxation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        slabs: Dictionary of {slab_id: StructureData}
        code: VASP code
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        parameters: VASP INCAR parameters
        options: Computer options
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary of {slab_id: vasp_task_node}
    """
    from .utils import prepare_vasp_parameters

    relax_nodes = {}

    for slab_id, slab_structure in slabs.items():
        # Prepare parameters
        vasp_params = prepare_vasp_parameters(
            base_parameters=parameters,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            kpoints_spacing=kpoints_spacing,
            options=options,
            clean_workdir=clean_workdir,
        )

        # Add VASP relax node
        relax_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"relax_slab_{slab_id}",
            structure=slab_structure,
            **vasp_params,
        )

    return relax_nodes


def build_energy_extraction_nodes(
    wg,
    vasp_nodes: dict[str, t.Any],
    node_type: str = "relaxed",
) -> dict[str, t.Any]:
    """
    Add energy extraction nodes for VASP calculations.

    Args:
        wg: WorkGraph instance to add nodes to
        vasp_nodes: Dictionary of {slab_id: vasp_task_node}
        node_type: Type of calculation ("relaxed" or "scf")

    Returns:
        Dictionary of {slab_id: energy_extraction_node}
    """
    # Import the existing extract_total_energy calcfunction
    from teros.core.slabs import extract_total_energy

    energy_nodes = {}

    for slab_id, vasp_node in vasp_nodes.items():
        energy_nodes[slab_id] = wg.add_task(
            extract_total_energy,
            name=f"extract_energy_{node_type}_{slab_id}",
            misc=vasp_node.outputs.misc,
        )

    return energy_nodes


def build_relaxation_energy_nodes(
    wg,
    unrelaxed_energies: dict[str, t.Any],
    relaxed_energies: dict[str, t.Any],
) -> dict[str, t.Any]:
    """
    Add relaxation energy calculation nodes.

    Args:
        wg: WorkGraph instance to add nodes to
        unrelaxed_energies: Dictionary of {slab_id: unrelaxed_energy_node}
        relaxed_energies: Dictionary of {slab_id: relaxed_energy_node}

    Returns:
        Dictionary of {slab_id: relaxation_energy_node}
    """
    relaxation_nodes = {}

    for slab_id in unrelaxed_energies.keys():
        relaxation_nodes[slab_id] = wg.add_task(
            calculate_relaxation_energy,
            name=f"relaxation_energy_{slab_id}",
            unrelaxed_energy=unrelaxed_energies[slab_id].outputs.result,
            relaxed_energy=relaxed_energies[slab_id].outputs.result,
        )

    return relaxation_nodes


@task.calcfunction
def calculate_relaxation_energy(
    unrelaxed_energy: orm.Float,
    relaxed_energy: orm.Float,
) -> orm.Float:
    """
    Calculate relaxation energy (unrelaxed - relaxed).

    Args:
        unrelaxed_energy: Energy before relaxation
        relaxed_energy: Energy after relaxation

    Returns:
        Relaxation energy (positive means structure lowered energy)
    """
    return orm.Float(unrelaxed_energy.value - relaxed_energy.value)
```

**Step 2: Commit slab operations**

```bash
git add teros/experimental/surface_thermo_preset_serial/slab_operations.py
git commit -m "feat: add slab operation node builders"
```

---

## Task 4: Implement Thermodynamics Operations Module

**Files:**
- Modify: `teros/experimental/surface_thermo_preset_serial/thermodynamics_operations.py`

**Step 1: Implement node builder for surface energies**

Write: `teros/experimental/surface_thermo_preset_serial/thermodynamics_operations.py`

```python
"""Node builders for thermodynamics calculations."""

import typing as t
from aiida import orm
from aiida_workgraph import task


def build_surface_energy_nodes(
    wg,
    bulk_structure: orm.StructureData,
    bulk_energy: t.Any,
    slab_structures: dict[str, orm.StructureData],
    slab_energies: dict[str, t.Any],
    reference_energies: t.Any,
    formation_enthalpy: t.Any,
    oxide_type: t.Any,
    sampling: int = 100,
) -> dict[str, t.Any]:
    """
    Add surface energy calculation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        bulk_structure: Bulk structure
        bulk_energy: Bulk energy node output
        slab_structures: Dictionary of {slab_id: StructureData}
        slab_energies: Dictionary of {slab_id: energy_node}
        reference_energies: Reference energies Dict node
        formation_enthalpy: Formation enthalpy Dict node
        oxide_type: Oxide type Str node ('binary' or 'ternary')
        sampling: Number of sampling points for chemical potential grid

    Returns:
        Dictionary of {slab_id: surface_energy_node}
    """
    # Import existing calcfunctions
    from teros.core.thermodynamics import (
        calculate_surface_energy_binary,
        calculate_surface_energy_ternary,
    )

    surface_energy_nodes = {}

    for slab_id, slab_structure in slab_structures.items():
        # Add conditional node based on oxide type
        # For simplicity, we'll create both and the workflow will use the right one

        # Binary surface energy
        binary_node = wg.add_task(
            calculate_surface_energy_binary,
            name=f"surface_energy_binary_{slab_id}",
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energies[slab_id].outputs.result,
            reference_energies=reference_energies,
            sampling=orm.Int(sampling),
        )

        # Ternary surface energy
        ternary_node = wg.add_task(
            calculate_surface_energy_ternary,
            name=f"surface_energy_ternary_{slab_id}",
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energies[slab_id].outputs.result,
            reference_energies=reference_energies,
            formation_enthalpy=formation_enthalpy,
            sampling=orm.Int(sampling),
        )

        # Store both for now (we'll handle selection in main workgraph)
        surface_energy_nodes[slab_id] = {
            'binary': binary_node,
            'ternary': ternary_node,
        }

    return surface_energy_nodes


@task.calcfunction
def select_surface_energy_by_oxide_type(
    oxide_type: orm.Str,
    binary_result: orm.Dict,
    ternary_result: orm.Dict,
) -> orm.Dict:
    """
    Select the appropriate surface energy result based on oxide type.

    Args:
        oxide_type: 'binary' or 'ternary'
        binary_result: Result from binary calculation
        ternary_result: Result from ternary calculation

    Returns:
        The appropriate result based on oxide type
    """
    if oxide_type.value == 'binary':
        return binary_result
    elif oxide_type.value == 'ternary':
        return ternary_result
    else:
        raise ValueError(f"Unknown oxide type: {oxide_type.value}")
```

**Step 2: Commit thermodynamics operations**

```bash
git add teros/experimental/surface_thermo_preset_serial/thermodynamics_operations.py
git commit -m "feat: add thermodynamics operation node builders"
```

---

## Task 5: Implement Main Workgraph

**Files:**
- Modify: `teros/experimental/surface_thermo_preset_serial/workgraph.py`

**Step 1: Implement main workgraph function (Part 1: Setup and bulk)**

Write: `teros/experimental/surface_thermo_preset_serial/workgraph.py`

```python
"""Main workgraph for serial surface thermodynamics preset."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph
from ase.io import read

from teros.core.slabs import generate_slab_structures, extract_total_energy
from teros.core.hf import calculate_formation_enthalpy
from teros.core.thermodynamics import identify_oxide_type

from .slab_operations import (
    build_scf_slabs_nodes,
    build_relax_slabs_nodes,
    build_energy_extraction_nodes,
    build_relaxation_energy_nodes,
)
from .thermodynamics_operations import (
    build_surface_energy_nodes,
    select_surface_energy_by_oxide_type,
)
from .utils import (
    prepare_vasp_parameters,
    create_default_bulk_parameters,
    create_default_slab_parameters,
    create_default_scf_parameters,
)

VaspWorkChain = WorkflowFactory('vasp.vasp')


@task.graph(outputs=[
    'bulk_energy', 'bulk_structure',
    'metal_energy', 'oxygen_energy', 'nonmetal_energy',
    'metal_structure', 'oxygen_structure', 'nonmetal_structure',
    'formation_enthalpy', 'reference_energies',
    'slab_structures', 'relaxed_slabs',
    'slab_energies', 'unrelaxed_slab_energies', 'relaxation_energies',
    'surface_energies', 'oxide_type',
])
def surface_thermodynamics_serial_workgraph(
    structures_dir: str = None,
    bulk_name: str = None,
    code_label: str = None,
    potential_family: str = None,
    bulk_potential_mapping: dict = None,
    kpoints_spacing: float = 0.4,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    clean_workdir: bool = False,
    # Reference materials
    metal_name: str = None,
    oxygen_name: str = None,
    nonmetal_name: str = None,
    metal_potential_mapping: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    oxygen_potential_mapping: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    nonmetal_potential_mapping: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    # Slab parameters
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
    # Control flags
    relax_slabs: bool = True,
    compute_thermodynamics: bool = True,
    thermodynamics_sampling: int = 100,
    compute_relaxation_energy: bool = True,
    input_slabs: dict = None,
):
    """
    Serial surface thermodynamics workflow with flat graph structure.

    All VASP calculation nodes are added directly to the main graph,
    allowing max_concurrent_jobs to control concurrent execution.

    Args:
        structures_dir: Directory containing structure files
        bulk_name: Bulk structure filename
        code_label: VASP code label
        potential_family: Potential family name
        bulk_potential_mapping: Element to potential mapping for bulk
        kpoints_spacing: K-points spacing for bulk
        bulk_parameters: VASP parameters for bulk
        bulk_options: Computer options for bulk
        clean_workdir: Whether to clean working directory
        metal_name: Metal reference structure filename
        oxygen_name: Oxygen reference structure filename
        nonmetal_name: Optional nonmetal reference structure filename
        metal_potential_mapping: Element to potential mapping for metal
        metal_parameters: VASP parameters for metal
        metal_options: Computer options for metal
        oxygen_potential_mapping: Element to potential mapping for oxygen
        oxygen_parameters: VASP parameters for oxygen
        oxygen_options: Computer options for oxygen
        nonmetal_potential_mapping: Element to potential mapping for nonmetal
        nonmetal_parameters: VASP parameters for nonmetal
        nonmetal_options: Computer options for nonmetal
        miller_indices: List of Miller indices to generate slabs
        min_slab_thickness: Minimum slab thickness (Angstroms)
        min_vacuum_thickness: Minimum vacuum thickness (Angstroms)
        slab_parameters: VASP parameters for slabs
        slab_options: Computer options for slabs
        slab_potential_mapping: Element to potential mapping for slabs
        slab_kpoints_spacing: K-points spacing for slabs
        lll_reduce: Whether to LLL reduce slab
        center_slab: Whether to center slab in cell
        symmetrize: Whether to symmetrize slab
        primitive: Whether to use primitive slab
        in_unit_planes: Whether Miller indices are in unit planes
        max_normal_search: Maximum normal search range
        relax_slabs: Whether to relax slabs
        compute_thermodynamics: Whether to compute thermodynamics
        thermodynamics_sampling: Number of chemical potential sampling points
        compute_relaxation_energy: Whether to compute relaxation energy
        input_slabs: Optional pre-generated slabs dict

    Returns:
        Dictionary of outputs
    """
    # Create WorkGraph instance
    wg = WorkGraph(name="surface_thermodynamics_serial")

    # Load code
    code = orm.load_code(code_label)

    # Set default parameters
    if bulk_parameters is None:
        bulk_parameters = create_default_bulk_parameters()
    if slab_parameters is None:
        slab_parameters = create_default_slab_parameters()
    if slab_kpoints_spacing is None:
        slab_kpoints_spacing = kpoints_spacing

    # =========================================================================
    # PHASE 1: Bulk calculations
    # =========================================================================

    # Load bulk structure
    bulk_filepath = f"{structures_dir}/{bulk_name}"
    bulk_structure = orm.StructureData(ase=read(bulk_filepath))

    # Add bulk relaxation node
    bulk_vasp_params = prepare_vasp_parameters(
        base_parameters=bulk_parameters,
        code=code,
        potential_family=potential_family,
        potential_mapping=bulk_potential_mapping,
        kpoints_spacing=kpoints_spacing,
        options=bulk_options,
        clean_workdir=clean_workdir,
    )

    bulk_node = wg.add_task(
        VaspWorkChain,
        name="bulk_relax",
        structure=bulk_structure,
        **bulk_vasp_params,
    )

    # Extract bulk energy
    bulk_energy_node = wg.add_task(
        extract_total_energy,
        name="extract_bulk_energy",
        misc=bulk_node.outputs.misc,
    )

    # =========================================================================
    # PHASE 2: Reference calculations (for thermodynamics)
    # =========================================================================

    reference_nodes = {}
    reference_energy_nodes = {}

    if compute_thermodynamics:
        # Metal reference
        if metal_name:
            metal_filepath = f"{structures_dir}/{metal_name}"
            metal_structure = orm.StructureData(ase=read(metal_filepath))

            metal_vasp_params = prepare_vasp_parameters(
                base_parameters=metal_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=metal_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=metal_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['metal'] = wg.add_task(
                VaspWorkChain,
                name="metal_relax",
                structure=metal_structure,
                **metal_vasp_params,
            )

            reference_energy_nodes['metal'] = wg.add_task(
                extract_total_energy,
                name="extract_metal_energy",
                misc=reference_nodes['metal'].outputs.misc,
            )

        # Oxygen reference
        if oxygen_name:
            oxygen_filepath = f"{structures_dir}/{oxygen_name}"
            oxygen_structure = orm.StructureData(ase=read(oxygen_filepath))

            oxygen_vasp_params = prepare_vasp_parameters(
                base_parameters=oxygen_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=oxygen_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=oxygen_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['oxygen'] = wg.add_task(
                VaspWorkChain,
                name="oxygen_relax",
                structure=oxygen_structure,
                **oxygen_vasp_params,
            )

            reference_energy_nodes['oxygen'] = wg.add_task(
                extract_total_energy,
                name="extract_oxygen_energy",
                misc=reference_nodes['oxygen'].outputs.misc,
            )

        # Optional nonmetal reference
        if nonmetal_name:
            nonmetal_filepath = f"{structures_dir}/{nonmetal_name}"
            nonmetal_structure = orm.StructureData(ase=read(nonmetal_filepath))

            nonmetal_vasp_params = prepare_vasp_parameters(
                base_parameters=nonmetal_parameters or bulk_parameters,
                code=code,
                potential_family=potential_family,
                potential_mapping=nonmetal_potential_mapping or bulk_potential_mapping,
                kpoints_spacing=kpoints_spacing,
                options=nonmetal_options or bulk_options,
                clean_workdir=clean_workdir,
            )

            reference_nodes['nonmetal'] = wg.add_task(
                VaspWorkChain,
                name="nonmetal_relax",
                structure=nonmetal_structure,
                **nonmetal_vasp_params,
            )

            reference_energy_nodes['nonmetal'] = wg.add_task(
                extract_total_energy,
                name="extract_nonmetal_energy",
                misc=reference_nodes['nonmetal'].outputs.misc,
            )

    # Continue in next step...
```

**Step 2: Continue main workgraph (Part 2: Slabs and thermodynamics)**

Append to `teros/experimental/surface_thermo_preset_serial/workgraph.py`:

```python
    # =========================================================================
    # Prepare reference energies dict
    # =========================================================================

    if compute_thermodynamics:
        # Build reference energies dict
        reference_energies_node = wg.add_task(
            build_reference_energies_dict,
            name="build_reference_energies",
            metal_energy=reference_energy_nodes.get('metal').outputs.result if 'metal' in reference_energy_nodes else None,
            oxygen_energy=reference_energy_nodes.get('oxygen').outputs.result if 'oxygen' in reference_energy_nodes else None,
            nonmetal_energy=reference_energy_nodes.get('nonmetal').outputs.result if 'nonmetal' in reference_energy_nodes else None,
            metal_structure=reference_nodes.get('metal').outputs.structure if 'metal' in reference_nodes else None,
            oxygen_structure=reference_nodes.get('oxygen').outputs.structure if 'oxygen' in reference_nodes else None,
            nonmetal_structure=reference_nodes.get('nonmetal').outputs.structure if 'nonmetal' in reference_nodes else None,
        )

        # Calculate formation enthalpy
        formation_enthalpy_node = wg.add_task(
            calculate_formation_enthalpy,
            name="formation_enthalpy",
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            reference_energies=reference_energies_node.outputs.result,
        )

        # Identify oxide type
        oxide_type_node = wg.add_task(
            identify_oxide_type,
            name="identify_oxide_type",
            bulk_structure=bulk_node.outputs.structure,
        )

    # =========================================================================
    # PHASE 3: Slab generation
    # =========================================================================

    if input_slabs:
        # Use provided slabs
        slab_structures = input_slabs
    else:
        # Generate slabs
        slab_gen_node = wg.add_task(
            generate_slab_structures,
            name="generate_slabs",
            bulk_structure=bulk_node.outputs.structure,
            miller_indices=orm.List(list=miller_indices or [[1,0,0], [1,1,0], [1,1,1]]),
            min_slab_thickness=orm.Float(min_slab_thickness),
            min_vacuum_thickness=orm.Float(min_vacuum_thickness),
            lll_reduce=orm.Bool(lll_reduce),
            center_slab=orm.Bool(center_slab),
            symmetrize=orm.Bool(symmetrize),
            primitive=orm.Bool(primitive),
            in_unit_planes=orm.Bool(in_unit_planes),
            max_normal_search=orm.Int(max_normal_search) if max_normal_search else None,
        )
        slab_structures = slab_gen_node.outputs.slabs

    # For direct node addition, we need to handle slabs as a known dict
    # In practice, we'll need to wait for slab generation or use input_slabs
    # For this implementation, let's assume input_slabs is provided as a dict

    if not input_slabs:
        raise NotImplementedError(
            "Dynamic slab generation not yet supported in serial mode. "
            "Please provide input_slabs as a dictionary."
        )

    # =========================================================================
    # PHASE 4: Slab calculations
    # =========================================================================

    # Build SCF (unrelaxed) nodes
    scf_params = create_default_scf_parameters()
    scf_nodes = build_scf_slabs_nodes(
        wg=wg,
        slabs=input_slabs,
        code=code,
        potential_family=potential_family,
        potential_mapping=slab_potential_mapping or bulk_potential_mapping,
        kpoints_spacing=slab_kpoints_spacing,
        parameters=scf_params,
        options=slab_options or bulk_options,
        clean_workdir=clean_workdir,
    )

    # Extract unrelaxed energies
    unrelaxed_energy_nodes = build_energy_extraction_nodes(
        wg=wg,
        vasp_nodes=scf_nodes,
        node_type="scf",
    )

    # Build relaxation nodes (if requested)
    if relax_slabs:
        relax_nodes = build_relax_slabs_nodes(
            wg=wg,
            slabs=input_slabs,
            code=code,
            potential_family=potential_family,
            potential_mapping=slab_potential_mapping or bulk_potential_mapping,
            kpoints_spacing=slab_kpoints_spacing,
            parameters=slab_parameters,
            options=slab_options or bulk_options,
            clean_workdir=clean_workdir,
        )

        # Extract relaxed energies
        relaxed_energy_nodes = build_energy_extraction_nodes(
            wg=wg,
            vasp_nodes=relax_nodes,
            node_type="relaxed",
        )

        # Calculate relaxation energies (if requested)
        if compute_relaxation_energy:
            relaxation_energy_nodes = build_relaxation_energy_nodes(
                wg=wg,
                unrelaxed_energies=unrelaxed_energy_nodes,
                relaxed_energies=relaxed_energy_nodes,
            )

    # =========================================================================
    # PHASE 5: Thermodynamics calculations
    # =========================================================================

    if compute_thermodynamics:
        # Use relaxed energies if available, otherwise SCF energies
        energy_nodes_for_thermo = relaxed_energy_nodes if relax_slabs else unrelaxed_energy_nodes

        # Build surface energy nodes
        surface_energy_nodes = build_surface_energy_nodes(
            wg=wg,
            bulk_structure=bulk_node.outputs.structure,
            bulk_energy=bulk_energy_node.outputs.result,
            slab_structures=input_slabs,
            slab_energies=energy_nodes_for_thermo,
            reference_energies=reference_energies_node.outputs.result,
            formation_enthalpy=formation_enthalpy_node.outputs.result,
            oxide_type=oxide_type_node.outputs.result,
            sampling=thermodynamics_sampling,
        )

        # Select appropriate surface energy based on oxide type
        final_surface_energies = {}
        for slab_id, surf_eng_nodes in surface_energy_nodes.items():
            final_surface_energies[slab_id] = wg.add_task(
                select_surface_energy_by_oxide_type,
                name=f"select_surface_energy_{slab_id}",
                oxide_type=oxide_type_node.outputs.result,
                binary_result=surf_eng_nodes['binary'].outputs.result,
                ternary_result=surf_eng_nodes['ternary'].outputs.result,
            )

    # =========================================================================
    # Collect outputs
    # =========================================================================

    outputs = {
        'bulk_energy': bulk_energy_node.outputs.result,
        'bulk_structure': bulk_node.outputs.structure,
        'slab_structures': input_slabs,
    }

    if compute_thermodynamics:
        outputs['metal_energy'] = reference_energy_nodes.get('metal').outputs.result if 'metal' in reference_energy_nodes else None
        outputs['oxygen_energy'] = reference_energy_nodes.get('oxygen').outputs.result if 'oxygen' in reference_energy_nodes else None
        outputs['nonmetal_energy'] = reference_energy_nodes.get('nonmetal').outputs.result if 'nonmetal' in reference_energy_nodes else None
        outputs['metal_structure'] = reference_nodes.get('metal').outputs.structure if 'metal' in reference_nodes else None
        outputs['oxygen_structure'] = reference_nodes.get('oxygen').outputs.structure if 'oxygen' in reference_nodes else None
        outputs['nonmetal_structure'] = reference_nodes.get('nonmetal').outputs.structure if 'nonmetal' in reference_nodes else None
        outputs['formation_enthalpy'] = formation_enthalpy_node.outputs.result
        outputs['reference_energies'] = reference_energies_node.outputs.result
        outputs['oxide_type'] = oxide_type_node.outputs.result

        # Surface energies
        outputs['surface_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in final_surface_energies.items()
        }

    if relax_slabs:
        outputs['relaxed_slabs'] = {
            slab_id: node.outputs.structure
            for slab_id, node in relax_nodes.items()
        }
        outputs['slab_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in relaxed_energy_nodes.items()
        }

    outputs['unrelaxed_slab_energies'] = {
        slab_id: node.outputs.result
        for slab_id, node in unrelaxed_energy_nodes.items()
    }

    if compute_relaxation_energy and relax_slabs:
        outputs['relaxation_energies'] = {
            slab_id: node.outputs.result
            for slab_id, node in relaxation_energy_nodes.items()
        }

    return outputs


@task.calcfunction
def build_reference_energies_dict(
    metal_energy: orm.Float = None,
    oxygen_energy: orm.Float = None,
    nonmetal_energy: orm.Float = None,
    metal_structure: orm.StructureData = None,
    oxygen_structure: orm.StructureData = None,
    nonmetal_structure: orm.StructureData = None,
) -> orm.Dict:
    """
    Build reference energies dictionary for thermodynamics calculations.

    Args:
        metal_energy: Total energy of metal reference
        oxygen_energy: Total energy of oxygen reference
        nonmetal_energy: Total energy of nonmetal reference
        metal_structure: Metal reference structure
        oxygen_structure: Oxygen reference structure
        nonmetal_structure: Nonmetal reference structure

    Returns:
        Dictionary with '*_energy_per_atom' keys
    """
    ref_dict = {}

    if metal_energy and metal_structure:
        metal_ase = metal_structure.get_ase()
        num_metal_atoms = len(metal_ase)
        ref_dict['metal_energy_per_atom'] = metal_energy.value / num_metal_atoms

    if oxygen_energy and oxygen_structure:
        oxygen_ase = oxygen_structure.get_ase()
        num_oxygen_atoms = len(oxygen_ase)
        ref_dict['oxygen_energy_per_atom'] = oxygen_energy.value / num_oxygen_atoms

    if nonmetal_energy and nonmetal_structure:
        nonmetal_ase = nonmetal_structure.get_ase()
        num_nonmetal_atoms = len(nonmetal_ase)
        ref_dict['nonmetal_energy_per_atom'] = nonmetal_energy.value / num_nonmetal_atoms

    return orm.Dict(dict=ref_dict)
```

**Step 3: Commit main workgraph**

```bash
git add teros/experimental/surface_thermo_preset_serial/workgraph.py
git commit -m "feat: implement main serial surface thermodynamics workgraph"
```

---

## Task 6: Create Test Script

**Files:**
- Create: `examples/vasp/test_surface_thermo_serial/test_serial_preset.py`
- Create: `examples/vasp/test_surface_thermo_serial/README.md`

**Step 1: Create test directory**

```bash
mkdir -p examples/vasp/test_surface_thermo_serial
```

**Step 2: Create test script**

Write: `examples/vasp/test_surface_thermo_serial/test_serial_preset.py`

```python
"""
Test script for serial surface thermodynamics preset.

This script tests the flat-graph implementation where all VASP nodes
exist at the same level, allowing max_concurrent_jobs to control execution.
"""

from aiida import orm
from aiida.engine import submit
from ase.io import read

from teros.experimental.surface_thermo_preset_serial import (
    surface_thermodynamics_serial_workgraph
)


def main():
    """Run serial surface thermodynamics test."""

    # Configuration
    structures_dir = "/path/to/structures"  # UPDATE THIS PATH

    # Load pre-generated slabs (required for this version)
    # In practice, you would generate these first or provide them
    input_slabs = {
        'slab_100': orm.StructureData(ase=read(f"{structures_dir}/slab_100.cif")),
        'slab_110': orm.StructureData(ase=read(f"{structures_dir}/slab_110.cif")),
        'slab_111': orm.StructureData(ase=read(f"{structures_dir}/slab_111.cif")),
        'slab_210': orm.StructureData(ase=read(f"{structures_dir}/slab_210.cif")),
    }

    # Build workgraph
    wg = surface_thermodynamics_serial_workgraph.build(
        # Structure files
        structures_dir=structures_dir,
        bulk_name="bulk.cif",

        # Code and potentials
        code_label="VASP-6.4.1@cluster02",
        potential_family="PBE.54",
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},

        # K-points
        kpoints_spacing=0.4,

        # Bulk parameters
        bulk_parameters={
            'PREC': 'Accurate',
            'EDIFF': 1e-6,
            'ENCUT': 520,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 3,
            'NSW': 100,
            'LWAVE': False,
            'LCHARG': False,
        },
        bulk_options={
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },

        # Reference materials for thermodynamics
        metal_name="Ag.cif",
        oxygen_name="O2.cif",
        metal_potential_mapping={'Ag': 'Ag'},
        oxygen_potential_mapping={'O': 'O'},

        # Slab parameters
        input_slabs=input_slabs,
        slab_parameters={
            'PREC': 'Accurate',
            'EDIFF': 1e-6,
            'ENCUT': 520,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'ISIF': 2,
            'NSW': 100,
            'LWAVE': False,
            'LCHARG': False,
        },
        slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        slab_kpoints_spacing=0.4,

        # Control flags
        relax_slabs=True,
        compute_thermodynamics=True,
        thermodynamics_sampling=100,
        compute_relaxation_energy=True,
    )

    # CRITICAL: Set max_concurrent_jobs to test concurrency control
    wg.max_concurrent_jobs = 2

    # Submit workflow
    result = submit(wg)

    print(f"Submitted WorkGraph: {result.pk}")
    print(f"Monitor with: verdi process show {result.pk}")
    print(f"Check concurrent jobs with: verdi process list")
    print()
    print("Expected behavior:")
    print("  - Maximum 2 VASP jobs running simultaneously")
    print("  - All VASP calculations complete successfully")
    print("  - Surface energies calculated for all slabs")

    return result


if __name__ == "__main__":
    main()
```

**Step 3: Create README**

Write: `examples/vasp/test_surface_thermo_serial/README.md`

```markdown
# Serial Surface Thermodynamics Test

This directory contains test scripts for the experimental serial surface thermodynamics preset.

## Purpose

Test that the flat-graph implementation correctly limits concurrent VASP jobs using `max_concurrent_jobs`.

## Prerequisites

1. Pre-generated slab structures (or modify script to use miller_indices)
2. VASP code configured: `VASP-6.4.1@cluster02`
3. Potential family: `PBE.54`
4. Structure files in a directory

## Running the Test

1. Update `structures_dir` path in `test_serial_preset.py`
2. Ensure AiiDA daemon is running: `verdi daemon start`
3. Run the test:
   ```bash
   source ~/envs/aiida/bin/activate
   python test_serial_preset.py
   ```

## Verification

Monitor execution:
```bash
verdi process list
```

Expected behavior:
- Maximum 2 VASP jobs running simultaneously (due to `max_concurrent_jobs=2`)
- All jobs complete successfully
- Main workgraph returns exit code [0]

Check results:
```bash
verdi process show <PK>
verdi process report <PK>
```

## Success Criteria

1. `verdi process list` never shows more than 2 running VASP jobs
2. All calculations complete successfully
3. Surface energies calculated for all slabs
4. Main node exits with [0]
```

**Step 4: Commit test files**

```bash
git add examples/vasp/test_surface_thermo_serial/
git commit -m "feat: add test script for serial surface thermodynamics preset"
```

---

## Task 7: Update Documentation

**Files:**
- Modify: `docs/plans/2025-11-02-surface-thermodynamics-serial-design.md`

**Step 1: Add implementation status to design doc**

Add to the end of `docs/plans/2025-11-02-surface-thermodynamics-serial-design.md`:

```markdown

---

## Implementation Status

**Date:** 2025-11-02
**Status:** Implemented

### Completed Components

- [x] Module structure created
- [x] Utilities module (parameter preparation)
- [x] Slab operations module (node builders)
- [x] Thermodynamics operations module (node builders)
- [x] Main workgraph implementation
- [x] Test script created

### Known Limitations

1. **Dynamic slab generation not supported**: Current implementation requires `input_slabs` to be provided as a dictionary. Dynamic slab generation from `miller_indices` requires additional work to handle the dynamic namespace properly.

2. **Simplified oxide type handling**: Currently creates both binary and ternary surface energy nodes and selects the appropriate one. Could be optimized to create only the needed type.

### Testing

Test script location: `examples/vasp/test_surface_thermo_serial/test_serial_preset.py`

Verification steps:
1. Run test with `max_concurrent_jobs=2`
2. Monitor with `verdi process list` to confirm max 2 concurrent VASP jobs
3. Verify successful completion with exit code [0]
4. Check surface energies calculated correctly

### Next Steps

1. Test with real data
2. Add support for dynamic slab generation
3. Optimize oxide type handling
4. Add support for optional features (cleavage, electronic properties, AIMD)
5. If successful, consider migrating other presets to flat structure
```

**Step 2: Commit documentation update**

```bash
git add docs/plans/2025-11-02-surface-thermodynamics-serial-design.md
git commit -m "docs: update implementation status in design document"
```

---

## Final Checklist

After completing all tasks:

- [ ] Module structure created in `teros/experimental/surface_thermo_preset_serial/`
- [ ] Utilities module implemented
- [ ] Slab operations module implemented
- [ ] Thermodynamics operations module implemented
- [ ] Main workgraph implemented
- [ ] Test script created
- [ ] Documentation updated
- [ ] All changes committed

## Testing Instructions

1. Clear Python cache:
   ```bash
   find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
   ```

2. Restart AiiDA daemon:
   ```bash
   verdi daemon restart
   ```

3. Update test script paths and run:
   ```bash
   source ~/envs/aiida/bin/activate
   python examples/vasp/test_surface_thermo_serial/test_serial_preset.py
   ```

4. Monitor execution:
   ```bash
   watch -n 2 'verdi process list'
   ```

5. Verify max 2 VASP jobs running simultaneously

6. Check final results:
   ```bash
   verdi process show <PK>
   verdi process report <PK>
   ```

## Success Criteria

- Main node returns exit code [0]
- Surface energies calculated for all slabs
- Maximum 2 VASP jobs ran concurrently (verified via `verdi process list`)
- All intermediate calculations (bulk, references, slabs) completed successfully
