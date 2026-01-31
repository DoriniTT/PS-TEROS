"""Slab generation brick for the lego module.

Generates slab structures at multiple thicknesses from a bulk structure.
Pure structure generation — no VASP calculations.

Reuses generate_thickness_series from teros.core.convergence.
"""

from aiida import orm
from aiida.common.links import LinkType


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate a slab_gen stage configuration.

    Args:
        stage: Stage configuration dict.
        stage_names: Set of stage names defined so far (before this stage).

    Raises:
        ValueError: If validation fails.
    """
    name = stage['name']

    # structure_from is required
    if 'structure_from' not in stage:
        raise ValueError(
            f"Stage '{name}': slab_gen stages require 'structure_from' "
            f"(name of a previous stage providing bulk structure)"
        )
    structure_from = stage['structure_from']
    if structure_from not in stage_names:
        raise ValueError(
            f"Stage '{name}' structure_from='{structure_from}' must reference "
            f"a previous stage name"
        )

    # miller_indices is required
    if 'miller_indices' not in stage:
        raise ValueError(
            f"Stage '{name}': slab_gen stages require 'miller_indices' "
            f"(e.g., [1, 1, 0])"
        )
    miller = stage['miller_indices']
    if not isinstance(miller, (list, tuple)) or len(miller) != 3:
        raise ValueError(
            f"Stage '{name}': miller_indices must be a list of 3 integers, "
            f"got: {miller}"
        )

    # layer_counts is required with at least 2 values
    if 'layer_counts' not in stage:
        raise ValueError(
            f"Stage '{name}': slab_gen stages require 'layer_counts' "
            f"(e.g., [3, 5, 7, 9, 11])"
        )
    layers = stage['layer_counts']
    if not isinstance(layers, (list, tuple)) or len(layers) < 2:
        raise ValueError(
            f"Stage '{name}': layer_counts must have at least 2 values, "
            f"got: {layers}"
        )


def create_stage_tasks(wg, stage, stage_name, context):
    """Create slab generation tasks in the WorkGraph.

    Args:
        wg: WorkGraph to add tasks to.
        stage: Stage configuration dict.
        stage_name: Unique stage identifier.
        context: Dict with shared context.

    Returns:
        Dict with task references for later stages.
    """
    from teros.core.convergence import generate_thickness_series
    from . import resolve_structure_from

    # Resolve bulk structure from referenced stage
    structure_from = stage['structure_from']
    input_structure = resolve_structure_from(structure_from, context)

    # Stage configuration
    miller_indices = stage['miller_indices']
    layer_counts = stage['layer_counts']
    min_vacuum_thickness = stage.get('min_vacuum_thickness', 15.0)
    termination_index = stage.get('termination_index', 0)
    lll_reduce = stage.get('lll_reduce', True)
    center_slab = stage.get('center_slab', True)
    primitive = stage.get('primitive', True)

    miller_list = orm.List(list=[int(m) for m in miller_indices])
    layers_list = orm.List(list=[int(n) for n in layer_counts])

    slab_gen_task = wg.add_task(
        generate_thickness_series,
        name=f'generate_slabs_{stage_name}',
        bulk_structure=input_structure,
        miller_indices=miller_list,
        layer_counts=layers_list,
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        primitive=orm.Bool(primitive),
        termination_index=orm.Int(termination_index),
    )

    return {
        'generate_slabs': slab_gen_task,
        'slabs': slab_gen_task.outputs.slabs,
        'structure': input_structure,
    }


def expose_stage_outputs(wg, stage_name, stage_tasks_result):
    """Expose slab generation stage outputs on the WorkGraph.

    Slab generation is an intermediate step; it does not expose outputs
    directly on the WorkGraph. Downstream stages (e.g., batch) consume
    the slabs namespace via stage_tasks.

    Args:
        wg: WorkGraph instance.
        stage_name: Unique stage identifier.
        stage_tasks_result: Dict returned by create_stage_tasks.
    """
    # No outputs exposed — slabs are consumed by the next stage (batch)
    pass


def get_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a slab generation stage.

    Args:
        wg_node: The WorkGraph ProcessNode.
        wg_pk: WorkGraph PK.
        stage_name: Name of the slab_gen stage.

    Returns:
        Dict with keys: slab_count, miller_indices, layer_counts,
        pk, stage, type.
    """
    result = {
        'slab_count': 0,
        'miller_indices': None,
        'layer_counts': None,
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'slab_gen',
    }

    # Try to find the generate_slabs calcfunction via link traversal
    gen_task_name = f'generate_slabs_{stage_name}'

    if hasattr(wg_node, 'base'):
        called_calc = wg_node.base.links.get_outgoing(
            link_type=LinkType.CALL_CALC
        )
        for link in called_calc.all():
            if gen_task_name in link.link_label:
                child_node = link.node
                # Count CREATE links (each slab is a StructureData)
                created = child_node.base.links.get_outgoing(
                    link_type=LinkType.CREATE
                )
                slab_count = sum(
                    1 for out_link in created.all()
                    if out_link.link_label.startswith('slabs__')
                )
                result['slab_count'] = slab_count

                # Extract miller_indices and layer_counts from inputs
                inputs = child_node.base.links.get_incoming()
                for inp_link in inputs.all():
                    if inp_link.link_label == 'miller_indices':
                        if hasattr(inp_link.node, 'get_list'):
                            result['miller_indices'] = inp_link.node.get_list()
                    if inp_link.link_label == 'layer_counts':
                        if hasattr(inp_link.node, 'get_list'):
                            result['layer_counts'] = inp_link.node.get_list()
                break

    return result


def print_stage_results(index: int, stage_name: str, stage_result: dict) -> None:
    """Print formatted results for a slab generation stage.

    Args:
        index: 1-based stage index.
        stage_name: Name of the stage.
        stage_result: Result dict from get_stage_results.
    """
    print(f"  [{index}] {stage_name} (SLAB GENERATION)")

    miller = stage_result.get('miller_indices')
    if miller:
        print(f"      Miller indices: ({', '.join(str(m) for m in miller)})")

    layer_counts = stage_result.get('layer_counts')
    if layer_counts:
        print(f"      Layer counts: {layer_counts}")

    slab_count = stage_result.get('slab_count', 0)
    print(f"      Slabs generated: {slab_count}")
