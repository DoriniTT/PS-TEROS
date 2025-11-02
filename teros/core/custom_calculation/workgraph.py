"""WorkGraph builder for custom VASP calculations."""

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task

from .tasks import extract_total_energy, extract_relaxed_structure


def build_custom_calculation_workgraph(
    structure,
    code_label,
    builder_inputs,
    name='custom_calc'
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

    Returns:
        WorkGraph ready to submit
    """
    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # Detect single vs multiple structures
    is_single = isinstance(structure, orm.StructureData)

    if is_single:
        # Single structure workflow
        # Prepare builder inputs (convert plain dicts to AiiDA types)
        prepared_inputs = _prepare_builder_inputs(builder_inputs)

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

        # Extract structure
        structure_task = wg.add_task(
            extract_relaxed_structure,
            name='extract_structure',
            misc=vasp_task.outputs.misc
        )

        # Set WorkGraph outputs
        wg.add_output('energy', energy_task.outputs.result)
        wg.add_output('structure', structure_task.outputs.result)
        wg.add_output('misc', vasp_task.outputs.misc)
    else:
        # Multiple structures workflow
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

        # Create VaspTask for each structure
        vasp_tasks = []
        energy_tasks = []
        structure_tasks = []

        for i, (struct, inputs) in enumerate(zip(structures, builder_inputs_list)):
            # Prepare builder inputs
            prepared_inputs = _prepare_builder_inputs(inputs)

            # Add VaspTask
            vasp_task = wg.add_task(
                VaspTask,
                name=f'vasp_calc_{i}',
                structure=struct,
                code=code,
                **prepared_inputs
            )
            vasp_tasks.append(vasp_task)

            # Extract energy
            energy_task = wg.add_task(
                extract_total_energy,
                name=f'extract_energy_{i}',
                misc=vasp_task.outputs.misc
            )
            energy_tasks.append(energy_task)

            # Extract structure
            structure_task = wg.add_task(
                extract_relaxed_structure,
                name=f'extract_structure_{i}',
                misc=vasp_task.outputs.misc
            )
            structure_tasks.append(structure_task)

        # For multiple structures, outputs are accessed via task names
        # (WorkGraph doesn't support dynamic list outputs with plain functions)
        # Users can access: wg.tasks['vasp_calc_0'], wg.tasks['extract_energy_0'], etc.

    return wg


def _prepare_builder_inputs(builder_inputs):
    """
    Prepare builder inputs by converting plain dicts to AiiDA types.

    Args:
        builder_inputs: dict with VASP builder parameters

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

    return prepared


def get_custom_results(workgraph):
    """
    Extract results from completed custom calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_custom_calculation_workgraph

    Returns:
        dict with:
            - energies: list of floats or single float
            - structures: list of StructureData or single StructureData
            - misc: list of dicts or single dict
    """
    results = {}

    # Check if single or multiple structure workflow
    if hasattr(workgraph.outputs, 'energy'):
        # Single structure
        results['energies'] = workgraph.outputs.energy.value
        results['structures'] = workgraph.outputs.structure
        results['misc'] = workgraph.outputs.misc.get_dict()
    else:
        # Multiple structures - access via task outputs
        i = 0
        energies = []
        structures = []
        misc_list = []

        # Find all energy extraction tasks
        while f'extract_energy_{i}' in workgraph.tasks:
            energy_task = workgraph.tasks[f'extract_energy_{i}']
            structure_task = workgraph.tasks[f'extract_structure_{i}']
            vasp_task = workgraph.tasks[f'vasp_calc_{i}']

            # Access outputs from completed tasks
            energies.append(energy_task.outputs['result'].value)
            structures.append(structure_task.outputs['result'].value)
            misc_list.append(vasp_task.outputs['misc'].value.get_dict())
            i += 1

        if not energies:
            raise ValueError("WorkGraph does not have expected outputs (energy or extract_energy_0 tasks)")

        results['energies'] = energies
        results['structures'] = structures
        results['misc'] = misc_list

    return results
