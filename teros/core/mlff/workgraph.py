"""Main workgraph builder for MLFF module."""
import typing as t
from aiida import orm
from aiida_workgraph import WorkGraph


def build_mlff_workgraph(
    structures: dict[str, t.Union[orm.StructureData, int]],
    training_steps: int,
    temperature: float,
    code_label: str,
    potential_family: str,
    potential_mapping: dict,
    production_steps: int = None,
    kpoints_spacing: float = 0.5,
    encut: float = 400,
    prec: str = 'Normal',
    options: dict = None,
    name: str = 'MLFFWorkGraph',
) -> WorkGraph:
    """
    Build MLFF workgraph for VASP.

    Simplified workflow based on VASP Liquid Si MLFF example:
    - Training: ML_ISTART=0, learns on-the-fly during MD
    - Production (optional): ML_ISTART=2, uses trained model

    Args:
        structures: {name: StructureData or PK}
        training_steps: NSW for training (e.g., 100-500)
        temperature: Temperature in K
        code_label: VASP code (e.g., 'VASP6.5.0@cluster02')
        potential_family: Potential family (e.g., 'PBE')
        potential_mapping: Element mapping (e.g., {'Si': 'Si'})
        production_steps: Optional NSW for production with trained model
        kpoints_spacing: K-points spacing (default: 0.5)
        encut: Energy cutoff in eV (default: 400)
        prec: VASP precision (default: 'Normal')
        options: Scheduler options (default: cluster02 settings)
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    Example (training only):
        wg = build_mlff_workgraph(
            structures={'si': structure},
            training_steps=200,
            temperature=300,
            code_label='VASP6.5.0@cluster02',
            potential_family='PBE',
            potential_mapping={'Si': 'Si'},
        )

    Example (training + production):
        wg = build_mlff_workgraph(
            structures={'si': structure},
            training_steps=200,
            production_steps=5000,
            temperature=300,
            code_label='VASP6.5.0@cluster02',
            potential_family='PBE',
            potential_mapping={'Si': 'Si'},
        )
    """
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

    # Default options for cluster02
    if options is None:
        options = {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        }

    # Base INCAR for MD+MLFF
    base_incar = {
        # DFT accuracy
        'PREC': prec,
        'ENCUT': encut,
        'EDIFF': 1e-5,
        'ALGO': 'Normal',
        'LREAL': 'Auto',

        # Electronic structure
        'ISMEAR': 0,      # Gaussian smearing for MD
        'SIGMA': 0.05,

        # MD settings
        'IBRION': 0,      # Molecular dynamics
        'MDALGO': 2,      # Nosé-Hoover thermostat
        'POTIM': 1.0,     # 1 fs timestep
        'SMASS': 0.0,     # Automatic Nosé mass

        # Output (needed for restart and ML)
        'LWAVE': True,
        'LCHARG': True,
    }

    # Build stages
    stages = []

    # Stage 0: Training (always included)
    stages.append({
        'TEBEG': temperature,
        'TEEND': temperature,
        'NSW': training_steps,
        'ML_LMLFF': True,
        'ML_ISTART': 0,   # Train from scratch
    })

    # Stage 1: Production (optional)
    if production_steps:
        stages.append({
            'TEBEG': temperature,
            'TEEND': temperature,
            'NSW': production_steps,
            'ML_LMLFF': True,
            'ML_ISTART': 2,   # Use trained model
        })

    # Create workgraph
    wg = WorkGraph(name=name)
    wg.max_number_jobs = 1
    code = orm.load_code(code_label)

    # Build stages sequentially
    current_structures = prepared_structures
    current_remote_folders = {}

    for stage_idx, stage_config in enumerate(stages):
        stage_task = wg.add_task(
            aimd_single_stage_scatter,
            slabs=current_structures,
            stage_config=stage_config,
            code=code,
            base_aimd_parameters=base_incar,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            kpoints_spacing=kpoints_spacing,
            clean_workdir=False,  # Keep ML files!
            restart_folders=current_remote_folders,
            max_number_jobs=1,
            name=f'mlff_stage_{stage_idx}',
        )

        current_structures = stage_task.outputs.structures
        current_remote_folders = stage_task.outputs.remote_folders

    return wg
