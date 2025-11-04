"""Main workgraph builder for MLFF module."""
import typing as t
from aiida import orm
from aiida_workgraph import WorkGraph


def build_mlff_workgraph(
    structures: dict[str, t.Union[orm.StructureData, int]],
    training_steps: int,
    production_steps: int,
    temperature: float,
    code_label: str,
    builder_inputs: dict,
    name: str = 'MLFFWorkGraph',
) -> WorkGraph:
    """
    Build MLFF workgraph with training and production phases.

    Simple two-stage workflow:
    - Stage 0: Training (ML_ISTART=0, DFT+ML learning)
    - Stage 1: Production (ML_ISTART=2, pure ML predictions)

    Args:
        structures: {name: StructureData or PK}
        training_steps: NSW for training phase (e.g., 200)
        production_steps: NSW for production phase (e.g., 5000)
        temperature: Temperature in K (both phases)
        code_label: VASP code label (e.g., 'VASP6.5.0@cluster02')
        builder_inputs: Base config with parameters, kpoints, potentials, options
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    Example:
        wg = build_mlff_workgraph(
            structures={'si_bulk': structure},
            training_steps=200,
            production_steps=5000,
            temperature=300,
            code_label='VASP6.5.0@cluster02',
            builder_inputs={
                'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400, ...}},
                'kpoints_spacing': 0.5,
                'potential_family': 'PBE',
                'potential_mapping': {'Si': 'Si'},
                'options': {'resources': {...}},
                'clean_workdir': False,
            },
        )
    """
    # Import AIMD function (reuse existing infrastructure)
    from teros.core.aimd_functions import aimd_single_stage_scatter

    # Prepare structures
    prepared_structures = {}
    for struct_name, struct_input in structures.items():
        if isinstance(struct_input, int):
            struct = orm.load_node(struct_input)
            if not isinstance(struct, orm.StructureData):
                raise ValueError(
                    f"PK {struct_input} for '{struct_name}' is not a StructureData"
                )
        else:
            struct = struct_input
        prepared_structures[struct_name] = struct

    # Create MLFF stages
    mlff_stages = [
        # Stage 0: Training
        {
            'TEBEG': temperature,
            'TEEND': temperature,
            'NSW': training_steps,
            'ML_LMLFF': True,
            'ML_ISTART': 0,  # Start training from scratch
        },
        # Stage 1: Production
        {
            'TEBEG': temperature,
            'TEEND': temperature,
            'NSW': production_steps,
            'ML_LMLFF': True,
            'ML_ISTART': 2,  # Use trained model
        },
    ]

    # Extract base INCAR
    if 'parameters' in builder_inputs and 'incar' in builder_inputs['parameters']:
        base_incar = builder_inputs['parameters']['incar'].copy()
    else:
        base_incar = {}

    # Create workgraph
    wg = WorkGraph(name=name)
    wg.max_number_jobs = 1  # Sequential execution

    # Load code
    code = orm.load_code(code_label)

    # Build sequential stages
    current_structures = prepared_structures
    current_remote_folders = {}

    for stage_idx, stage_config in enumerate(mlff_stages):
        # Merge base INCAR with stage MLFF parameters
        stage_params = {**base_incar, **stage_config}

        # Create stage task
        stage_task = wg.add_task(
            aimd_single_stage_scatter,
            slabs=current_structures,
            stage_config=stage_config,
            code=code,
            base_aimd_parameters=base_incar,
            potential_family=builder_inputs.get('potential_family', 'PBE'),
            potential_mapping=builder_inputs.get('potential_mapping', {}),
            options=builder_inputs.get('options', {}),
            kpoints_spacing=builder_inputs.get('kpoints_spacing', 0.5),
            clean_workdir=builder_inputs.get('clean_workdir', False),
            restart_folders=current_remote_folders,
            max_number_jobs=1,
            name=f'stage_{stage_idx}_mlff',
        )

        # Update for next stage
        current_structures = stage_task.outputs.structures
        current_remote_folders = stage_task.outputs.remote_folders

    return wg
