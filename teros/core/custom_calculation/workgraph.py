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
        # Multiple structures - implement in next task
        raise NotImplementedError("Multiple structure support not yet implemented")

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

    # Copy scalar values directly
    for key in ('kpoints_spacing', 'potential_family', 'clean_workdir'):
        if key in builder_inputs:
            prepared[key] = builder_inputs[key]

    return prepared


def get_custom_results(workgraph):
    """
    Extract results from completed custom calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_custom_calculation_workgraph

    Returns:
        dict with energies, structures, misc outputs
    """
    raise NotImplementedError("Function not yet implemented - placeholder for Task 1")
