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


def _get_builder_for_structure_stage(
    struct_name: str,
    stage_idx: int,
    base_builder: dict,
    structure_overrides: dict,
    stage_overrides: dict,
    matrix_overrides: dict,
) -> dict:
    """
    Merge builder inputs for specific (structure, stage) combination.

    Priority: matrix > stage > structure > base

    Args:
        struct_name: Structure name
        stage_idx: Stage index (0-based)
        base_builder: Default builder inputs
        structure_overrides: Per-structure overrides
        stage_overrides: Per-stage overrides
        matrix_overrides: Per-(structure,stage) overrides

    Returns:
        Merged builder inputs for this specific combination
    """
    result = base_builder.copy()

    # Apply structure override
    if struct_name in structure_overrides:
        result = merge_builder_inputs(result, structure_overrides[struct_name])

    # Apply stage override
    if stage_idx in stage_overrides:
        result = merge_builder_inputs(result, stage_overrides[stage_idx])

    # Apply matrix override
    matrix_key = (struct_name, stage_idx)
    if matrix_key in matrix_overrides:
        result = merge_builder_inputs(result, matrix_overrides[matrix_key])

    return result


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
        structure_overrides: NOT IMPLEMENTED - reserved for future use
        stage_overrides: NOT IMPLEMENTED - reserved for future use
        matrix_overrides: NOT IMPLEMENTED - reserved for future use
        max_concurrent_jobs: Limit parallel VASP calculations (None = unlimited)
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    CURRENT LIMITATIONS:
        The structure_overrides, stage_overrides, and matrix_overrides parameters
        are accepted but NOT FUNCTIONAL in this version. All structures in all stages
        use the same builder_inputs.

        This is because aimd_single_stage_scatter() currently accepts a single
        aimd_parameters dict for all structures, not per-structure parameters.

        To implement full override support, aimd_single_stage_scatter() must be
        modified to accept per-structure builder configurations.

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

    # 2. Run sequential AIMD stages
    # Import here to avoid circular import (aimd_single_stage_scatter is in parent aimd.py)
    from . import aimd_single_stage_scatter

    stage_results = {}
    current_structures = prepared_structures
    current_remote_folders = None

    for stage_idx, stage_config in enumerate(aimd_stages):
        temperature = stage_config['temperature']
        steps = stage_config['steps']

        # TODO: Implement per-structure override system
        # Currently, all structures use the same builder_inputs.
        # To implement overrides, need to either:
        #   1. Modify aimd_single_stage_scatter to accept per-structure parameters, OR
        #   2. Call aimd_single_stage_scatter separately for each structure
        #
        # The helper function _get_builder_for_structure_stage() is ready to use
        # once the underlying scatter function supports per-structure configuration.

        # Load code
        code = orm.load_code(code_label)

        # Prepare AIMD parameters from builder_inputs
        # Extract base AIMD INCAR parameters
        if 'parameters' in builder_inputs and 'incar' in builder_inputs['parameters']:
            aimd_base_params = builder_inputs['parameters']['incar'].copy()
        else:
            aimd_base_params = {}

        # Create stage task
        stage_task = wg.add_task(
            aimd_single_stage_scatter,
            slabs=current_structures,
            temperature=temperature,
            steps=steps,
            code=code,
            aimd_parameters=aimd_base_params,
            potential_family=builder_inputs.get('potential_family', 'PBE'),
            potential_mapping=builder_inputs.get('potential_mapping', {}),
            options=builder_inputs.get('options', {}),
            kpoints_spacing=builder_inputs.get('kpoints_spacing', 0.5),
            clean_workdir=builder_inputs.get('clean_workdir', False),
            restart_folders=current_remote_folders if current_remote_folders else {},
            max_number_jobs=max_concurrent_jobs,
            name=f'stage_{stage_idx}_aimd',
        )

        # Store results for this stage
        stage_results[stage_idx] = {
            'structures': stage_task.outputs.structures,
            'energies': stage_task.outputs.energies,
            'remote_folders': stage_task.outputs.remote_folders,
        }

        # Update for next stage
        current_structures = stage_task.outputs.structures
        current_remote_folders = stage_task.outputs.remote_folders

    # 3. Set workgraph outputs
    wg.add_output('results', stage_results)
    if supercell_outputs:
        wg.add_output('supercells', supercell_outputs)

    return wg
