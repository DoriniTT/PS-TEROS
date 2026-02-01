"""Brick registry and shared helpers for the lego module.

Each brick module exports a PORTS dict plus 5 functions:
    PORTS: dict with 'inputs' and 'outputs' port declarations (from connections.py)
    validate_stage(stage, stage_names) -> None
    create_stage_tasks(wg, stage, stage_name, context) -> dict
    expose_stage_outputs(wg, stage_name, stage_tasks_result) -> None
    get_stage_results(wg_node, wg_pk, stage_name) -> dict
    print_stage_results(index, stage_name, stage_result) -> None

Port declarations and validation logic live in connections.py (pure Python,
no AiiDA dependency) so they can be imported in tier1 tests.
"""

from . import vasp, dos, batch, bader, convergence
from .connections import (  # noqa: F401
    PORT_TYPES,
    ALL_PORTS,
    validate_connections,
    _validate_port_types,
    _evaluate_conditional,
    get_brick_info,
)


BRICK_REGISTRY = {
    'vasp': vasp,
    'dos': dos,
    'batch': batch,
    'bader': bader,
    'convergence': convergence,
}

VALID_BRICK_TYPES = tuple(BRICK_REGISTRY.keys())


def get_brick_module(brick_type: str):
    """Look up a brick module by type string.

    Args:
        brick_type: One of 'vasp', 'dos', 'batch', 'bader', 'convergence'.

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

    Only VASP stages produce a meaningful structure output. Referencing
    a non-VASP stage (dos, batch, bader, convergence) raises an error
    because those bricks don't produce structures.

    Args:
        structure_from: Name of the stage to get structure from.
        context: The context dict passed to create_stage_tasks.

    Returns:
        Structure socket (StructureData or task output socket).

    Raises:
        ValueError: If the referenced stage doesn't produce a structure.
    """
    stage_tasks = context['stage_tasks']
    stage_types = context['stage_types']

    ref_stage_type = stage_types.get(structure_from, 'vasp')
    if ref_stage_type == 'vasp':
        return stage_tasks[structure_from]['vasp'].outputs.structure
    else:
        # Non-VASP bricks don't produce structures
        raise ValueError(
            f"structure_from='{structure_from}' references a '{ref_stage_type}' "
            f"stage, which doesn't produce a structure output. "
            f"Point to a VASP stage instead."
        )
