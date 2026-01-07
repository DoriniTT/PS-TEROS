"""WorkGraph builder for custom VASP calculations."""

import re
import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task, namespace, dynamic

from .tasks import extract_total_energy
from ..fixed_atoms import get_fixed_atoms_list


def _sanitize_key(label: str, index: int) -> str:
    """
    Sanitize a label to be a valid Python identifier for AiiDA link labels.

    Args:
        label: The original label string
        index: Fallback index if label is empty or invalid

    Returns:
        A valid Python identifier string
    """
    if not label:
        return f'structure_{index}'

    # Replace common problematic characters
    key = label.replace('.', '_').replace(' ', '_').replace('-', '_')

    # Remove any remaining non-alphanumeric characters (except underscore)
    key = re.sub(r'[^a-zA-Z0-9_]', '', key)

    # Ensure it doesn't start with a digit (prepend 's_' if it does)
    if key and key[0].isdigit():
        key = f's_{key}'

    # If empty after sanitization, use fallback
    if not key:
        return f'structure_{index}'

    return key


@task.graph
def CustomCalculationMultipleWorkGraph(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.InstalledCode,
    builder_inputs_list: list,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list = None,
) -> t.Annotated[dict, namespace(**{
    'structures': dynamic(orm.StructureData),
    'energies': dynamic(orm.Float),
    'misc': dynamic(orm.Dict),
})]:
    """
    WorkGraph task for multiple custom VASP calculations.

    Args:
        structures: Dict of StructureData keyed by sanitized label
        code: VASP InstalledCode
        builder_inputs_list: List of builder input dicts (same order as structures)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix

    Returns:
        Dictionary with namespace outputs:
            - structures: Dict {key: StructureData} for relaxed structures
            - energies: Dict {key: Float} for energies
            - misc: Dict {key: Dict} for misc outputs
    """
    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    structures_out = {}
    energies_out = {}
    misc_out = {}

    # Process each structure
    keys = list(structures.keys())
    for i, key in enumerate(keys):
        struct = structures[key]
        inputs = builder_inputs_list[i]

        # Prepare builder inputs
        prepared_inputs = _prepare_builder_inputs(
            inputs, struct, fix_type, fix_thickness, fix_elements
        )

        # Add VaspTask
        vasp_task = VaspTask(
            structure=struct,
            code=code,
            **prepared_inputs
        )

        # Extract energy
        energy_task = extract_total_energy(misc=vasp_task.misc)

        # Collect outputs
        structures_out[key] = vasp_task.structure
        energies_out[key] = energy_task.result
        misc_out[key] = vasp_task.misc

    return {
        'structures': structures_out,
        'energies': energies_out,
        'misc': misc_out,
    }


def build_custom_calculation_workgraph(
    structure,
    code_label,
    builder_inputs,
    name='custom_calc',
    max_concurrent_jobs=None,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
):
    """
    Build a WorkGraph for custom VASP calculations.

    Args:
        structure: StructureData or list of StructureData
        code_label: str, VASP code label (e.g., 'VASP-6.4.1@cluster02')
        builder_inputs: dict or list of dicts with VASP builder parameters:
            - parameters: Dict with nested 'incar' dict
            - options: Dict with resources, queue, etc.
            - kpoints_spacing: float
            - potential_family: str
            - potential_mapping: dict
            - clean_workdir: bool
            - settings: Dict (optional)
        name: str, WorkGraph name
        max_concurrent_jobs: int, optional - Maximum number of concurrent VASP calculations
                           (only applies to multiple structure workflows). Default: None (unlimited)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None). Default: None (no fixing)
        fix_thickness: Thickness in Angstroms for fixing region. Default: 0.0
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O']).
                     If None, all elements in the region are fixed. Default: None

    Returns:
        WorkGraph ready to submit
    """
    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Detect single vs multiple structures
    is_single = isinstance(structure, orm.StructureData)

    if is_single:
        # Single structure workflow - use plain WorkGraph
        wg = WorkGraph(name=name)

        # Prepare builder inputs (convert plain dicts to AiiDA types)
        prepared_inputs = _prepare_builder_inputs(
            builder_inputs, structure, fix_type, fix_thickness, fix_elements
        )

        # Add VaspTask
        vasp_task = wg.add_task(
            VaspTask,
            name='vasp_calc',
            structure=structure,
            code=code,
            **prepared_inputs
        )

        # Extract energy
        energy_task = wg.add_task(
            extract_total_energy,
            name='extract_energy',
            misc=vasp_task.outputs.misc
        )

        # Set WorkGraph outputs (use 'any' socket type identifier)
        wg.add_output('any', 'energy', from_socket=energy_task.outputs.result)
        wg.add_output('any', 'structure', from_socket=vasp_task.outputs.structure)
        wg.add_output('any', 'misc', from_socket=vasp_task.outputs.misc)

        return wg
    else:
        # Multiple structures workflow - use @task.graph pattern
        structures = structure  # Rename for clarity

        # Handle builder_inputs: single dict or list of dicts
        if isinstance(builder_inputs, dict):
            # Same builder inputs for all structures
            builder_inputs_list = [builder_inputs] * len(structures)
        elif isinstance(builder_inputs, list):
            # Different builder inputs for each structure
            if len(builder_inputs) != len(structures):
                raise ValueError(
                    f"builder_inputs list length ({len(builder_inputs)}) must match "
                    f"structures list length ({len(structures)})"
                )
            builder_inputs_list = builder_inputs
        else:
            raise TypeError("builder_inputs must be dict or list of dicts")

        # Prepare structures dict with sanitized keys
        structures_dict = {}
        for i, struct in enumerate(structures):
            label = struct.label if hasattr(struct, 'label') else ''
            key = _sanitize_key(label, i)
            structures_dict[key] = struct

        # Build WorkGraph using @task.graph function
        wg = CustomCalculationMultipleWorkGraph.build(
            structures=structures_dict,
            code=code,
            builder_inputs_list=builder_inputs_list,
            fix_type=fix_type,
            fix_thickness=fix_thickness,
            fix_elements=fix_elements,
        )

        # Set name and max concurrent jobs
        wg.name = name
        if max_concurrent_jobs is not None:
            wg.max_number_jobs = max_concurrent_jobs

        return wg


def _prepare_builder_inputs(
    builder_inputs,
    structure: orm.StructureData = None,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
):
    """
    Prepare builder inputs by converting plain dicts to AiiDA types.

    Also applies selective dynamics if fix_type is specified.

    Args:
        builder_inputs: dict with VASP builder parameters
        structure: StructureData to apply selective dynamics to (optional)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix

    Returns:
        dict with AiiDA-compatible types
    """
    prepared = {}

    # Convert parameters dict to orm.Dict
    if 'parameters' in builder_inputs:
        if isinstance(builder_inputs['parameters'], dict):
            prepared['parameters'] = orm.Dict(dict=builder_inputs['parameters'])
        else:
            prepared['parameters'] = builder_inputs['parameters']

    # Convert options dict to orm.Dict
    if 'options' in builder_inputs:
        if isinstance(builder_inputs['options'], dict):
            prepared['options'] = orm.Dict(dict=builder_inputs['options'])
        else:
            prepared['options'] = builder_inputs['options']

    # Convert potential_mapping dict to orm.Dict
    if 'potential_mapping' in builder_inputs:
        if isinstance(builder_inputs['potential_mapping'], dict):
            prepared['potential_mapping'] = orm.Dict(dict=builder_inputs['potential_mapping'])
        else:
            prepared['potential_mapping'] = builder_inputs['potential_mapping']

    # Convert settings dict to orm.Dict if present
    if 'settings' in builder_inputs:
        if isinstance(builder_inputs['settings'], dict):
            prepared['settings'] = orm.Dict(dict=builder_inputs['settings'])
        else:
            prepared['settings'] = builder_inputs['settings']

    # Handle kpoints_spacing - ensure it's a float
    if 'kpoints_spacing' in builder_inputs:
        kps = builder_inputs['kpoints_spacing']
        # Keep as plain Python float (don't wrap in orm.Float)
        if isinstance(kps, (int, float)):
            prepared['kpoints_spacing'] = float(kps)
        else:
            prepared['kpoints_spacing'] = kps

    # Copy string/bool values directly
    for key in ('potential_family', 'clean_workdir'):
        if key in builder_inputs:
            prepared[key] = builder_inputs[key]

    # Apply selective dynamics if fix_type is specified
    if structure is not None and fix_type is not None and fix_thickness > 0.0:
        # Get list of atoms to fix (1-based indices)
        fixed_atoms_list = get_fixed_atoms_list(
            structure=structure,
            fix_type=fix_type,
            fix_thickness=fix_thickness,
            fix_elements=fix_elements,
        )

        if fixed_atoms_list:
            # Create positions_dof array for selective dynamics
            # True = relax, False = fix
            num_atoms = len(structure.sites)
            positions_dof = []

            for i in range(1, num_atoms + 1):  # 1-based indexing
                if i in fixed_atoms_list:
                    positions_dof.append([False, False, False])  # Fix atom
                else:
                    positions_dof.append([True, True, True])  # Relax atom

            # Add selective dynamics to VASP inputs
            prepared['dynamics'] = orm.Dict(dict={
                'positions_dof': positions_dof
            })

    return prepared


def get_custom_results(workgraph):
    """
    Extract results from completed custom calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_custom_calculation_workgraph

    Returns:
        dict with:
            - energies: dict (keyed by structure label) or single float
            - structures: dict of StructureData (keyed by structure label), or single StructureData
            - misc: dict (keyed by structure label) or single dict

    For multiple structures, results are organized by structure key (derived from
    structure labels, or 'structure_0', 'structure_1', etc. if no label):
        {
            'energies': {'my_struct_1': -123.45, 'my_struct_2': -234.56, ...},
            'structures': {'my_struct_1': <StructureData>, 'my_struct_2': <StructureData>, ...},
            'misc': {'my_struct_1': {...}, 'my_struct_2': {...}, ...}
        }

    Example:
        >>> results = get_custom_results(wg)
        >>> for key, energy in results['energies'].items():
        ...     print(f"{key}: {energy:.3f} eV")
        >>> # Access structures directly
        >>> for key, struct in results['structures'].items():
        ...     print(f"{key}: {struct.get_formula()}")
    """
    results = {}

    # Check if single or multiple structure workflow
    if hasattr(workgraph.outputs, 'energy'):
        # Single structure
        results['energies'] = workgraph.outputs.energy.value
        results['structures'] = workgraph.outputs.structure
        results['misc'] = workgraph.outputs.misc.get_dict()
    elif hasattr(workgraph.outputs, 'energies'):
        # Multiple structures with namespaced outputs
        # Access energies namespace
        energies = {}
        for key in workgraph.outputs.energies:
            energies[key] = workgraph.outputs.energies[key].value
        results['energies'] = energies

        # Access structures namespace
        structures = {}
        for key in workgraph.outputs.structures:
            structures[key] = workgraph.outputs.structures[key]
        results['structures'] = structures

        # Access misc namespace
        misc_dict = {}
        for key in workgraph.outputs.misc:
            misc_dict[key] = workgraph.outputs.misc[key].get_dict()
        results['misc'] = misc_dict
    else:
        # Legacy: Multiple structures - access via task outputs
        i = 0
        energies = {}
        structures = {}
        misc_dict = {}

        # Find all energy extraction tasks
        while f'extract_energy_{i}' in workgraph.tasks:
            energy_task = workgraph.tasks[f'extract_energy_{i}']
            vasp_task = workgraph.tasks[f'vasp_calc_{i}']

            # Get structure to determine key (use same sanitization as workgraph builder)
            struct = vasp_task.outputs['structure'].value
            label = struct.label if hasattr(struct, 'label') else ''
            key = _sanitize_key(label, i)

            # Access outputs from completed tasks
            energies[key] = energy_task.outputs['result'].value
            structures[key] = struct
            misc_dict[key] = vasp_task.outputs['misc'].value.get_dict()
            i += 1

        if not energies:
            raise ValueError("WorkGraph does not have expected outputs (energy or extract_energy_0 tasks)")

        results['energies'] = energies
        results['structures'] = structures
        results['misc'] = misc_dict

    return results
