"""Main workgraph builder for AIMD module."""
import typing as t
from aiida import orm
from aiida_workgraph import WorkGraph
from .utils import (
    validate_stage_sequence,
    validate_supercell_spec,
    merge_builder_inputs,
)
from .tasks import create_supercell


def build_aimd_workgraph(
    structures: dict[str, t.Union[orm.StructureData, int]],
    aimd_stages: list[dict],
    code_label: str,
    builder_inputs: dict,
    supercell_specs: dict[str, list[int]] = None,
    structure_overrides: dict[str, dict] = None,
    stage_overrides: dict[int, dict] = None,
    matrix_overrides: dict[tuple, dict] = None,
    max_concurrent_jobs: int = None,
    name: str = 'AIMDWorkGraph',
) -> WorkGraph:
    """
    Build AIMD workgraph with sequential stages.

    Args:
        structures: {name: StructureData or PK} - input structures
        aimd_stages: [{'temperature': K, 'steps': N}, ...] - sequential stages
        code_label: VASP code label (e.g., 'VASP6.5.0@cluster02')
        builder_inputs: Default builder config for all (structure, stage) combinations
        supercell_specs: {structure_name: [nx, ny, nz]} - optional supercell per structure
        structure_overrides: {structure_name: builder_inputs} - override per structure
        stage_overrides: {stage_idx: builder_inputs} - override per stage (0-indexed)
        matrix_overrides: {(structure_name, stage_idx): builder_inputs} - specific overrides
        max_concurrent_jobs: Limit parallel VASP calculations (None = unlimited)
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    Override priority: matrix_overrides > stage_overrides > structure_overrides > builder_inputs

    Example:
        wg = build_aimd_workgraph(
            structures={'slab1': structure1, 'slab2': pk2},
            aimd_stages=[
                {'temperature': 300, 'steps': 100},
                {'temperature': 300, 'steps': 500},
            ],
            code_label='VASP6.5.0@cluster02',
            builder_inputs={
                'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400}},
                'kpoints_spacing': 0.5,
                'potential_family': 'PBE',
                'potential_mapping': {},
                'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
                'clean_workdir': False,
            },
            supercell_specs={'slab1': [2, 2, 1]},
            max_concurrent_jobs=4,
        )
    """
    # Validate inputs
    validate_stage_sequence(aimd_stages)

    if supercell_specs:
        for struct_name, spec in supercell_specs.items():
            if struct_name not in structures:
                raise ValueError(f"supercell_specs key '{struct_name}' not in structures")
            validate_supercell_spec(spec)

    # Initialize defaults
    if structure_overrides is None:
        structure_overrides = {}
    if stage_overrides is None:
        stage_overrides = {}
    if matrix_overrides is None:
        matrix_overrides = {}

    # Create workgraph
    wg = WorkGraph(name=name)

    if max_concurrent_jobs:
        wg.max_number_jobs = max_concurrent_jobs

    # 1. Load and prepare structures
    prepared_structures = {}
    supercell_outputs = {}

    for struct_name, struct_input in structures.items():
        # Load structure if PK
        if isinstance(struct_input, int):
            struct = orm.load_node(struct_input)
            if not isinstance(struct, orm.StructureData):
                raise ValueError(
                    f"PK {struct_input} for '{struct_name}' is not a StructureData node"
                )
        else:
            struct = struct_input

        # Create supercell if requested
        if supercell_specs and struct_name in supercell_specs:
            sc_task = wg.add_task(
                create_supercell,
                structure=struct,
                spec=supercell_specs[struct_name],
                name=f'create_supercell_{struct_name}',
            )
            prepared_structures[struct_name] = sc_task.outputs.result
            supercell_outputs[struct_name] = sc_task.outputs.result
        else:
            prepared_structures[struct_name] = struct

    # TODO: Stage loop implementation in next task

    return wg
