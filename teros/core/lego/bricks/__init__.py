"""Brick registry and shared helpers for the lego module.

Each brick module exports exactly 5 functions:
    validate_stage(stage, stage_names) -> None
    create_stage_tasks(wg, stage, stage_name, context) -> dict
    expose_stage_outputs(wg, stage_name, stage_tasks_result) -> None
    get_stage_results(wg_node, wg_pk, stage_name) -> dict
    print_stage_results(index, stage_name, stage_result) -> None
"""

from . import vasp, dos, batch, bader, thickness


BRICK_REGISTRY = {
    'vasp': vasp,
    'dos': dos,
    'batch': batch,
    'bader': bader,
    'thickness': thickness,
}

VALID_BRICK_TYPES = tuple(BRICK_REGISTRY.keys())


def get_brick_module(brick_type: str):
    """Look up a brick module by type string.

    Args:
        brick_type: One of 'vasp', 'dos', 'batch', 'bader'.

    Returns:
        The brick module.

    Raises:
        ValueError: If the brick type is unknown.
    """
    try:
        return BRICK_REGISTRY[brick_type]
    except KeyError:
        raise ValueError(
            f"Unknown brick type '{brick_type}'. "
            f"Must be one of {VALID_BRICK_TYPES}"
        )


def resolve_structure_from(structure_from: str, context: dict):
    """Resolve a structure socket from a previous stage.

    Args:
        structure_from: Name of the stage to get structure from.
        context: The context dict passed to create_stage_tasks.

    Returns:
        Structure socket (StructureData or task output socket).
    """
    stage_tasks = context['stage_tasks']
    stage_types = context['stage_types']

    ref_stage_type = stage_types.get(structure_from, 'vasp')
    if ref_stage_type == 'thickness':
        # Thickness stages expose the recommended slab structure
        return stage_tasks[structure_from]['structure']
    elif ref_stage_type in ('dos', 'batch', 'bader'):
        # DOS/batch/bader stages don't modify structure, use their input structure
        return stage_tasks[structure_from]['structure']
    else:
        return stage_tasks[structure_from]['vasp'].outputs.structure
