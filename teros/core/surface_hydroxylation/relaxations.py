"""Child WorkGraph for parallel VASP relaxations with semaphore limiting."""

import typing as t
from aiida import orm
from aiida.engine import calcfunction
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic


@calcfunction
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP energies output.

    Args:
        energies: Dictionary containing energy outputs from VASP (from misc output)

    Returns:
        Total energy as Float
    """
    energy_dict = energies.get_dict()
    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in VASP outputs. Available keys: {available}')


@calcfunction
def get_exit_info(exit_status: orm.Int, exit_message: orm.Str) -> dict:
    """
    Package exit status and message for output.

    Args:
        exit_status: Exit status from VASP
        exit_message: Exit message from VASP

    Returns:
        dict with exit_status and error
    """
    from aiida.orm import Int, Str

    # Create new nodes to avoid provenance cycles
    return {
        'exit_status': Int(exit_status.value),
        'error': Str(exit_message.value if exit_message.value else '')
    }


@task.graph
def relax_slabs_with_semaphore(
    structures: dict[str, orm.StructureData],
    code_pk: int,
    vasp_config: dict,
    options: dict,
    max_parallel: int,
) -> t.Annotated[dict, namespace(**{
    'structures': dynamic(orm.StructureData),
    'energies': dynamic(orm.Float),
    'exit_statuses': dynamic(orm.Int),
    'errors': dynamic(orm.Str),
})]:
    """
    Relax multiple structures in parallel with batch-based limiting.

    This is the core relaxation WorkGraph that receives structures from generate_structures
    (in structure_0, structure_1, ... format) and returns results compatible with
    collect_results (in structure_N, energy_N, exit_status_N, error_N format).

    PARALLELIZATION APPROACH:
    This implementation uses a simple batch approach: it processes only the first
    max_parallel structures from the input. This allows manual control over how
    many calculations run at once.

    Example workflow:
    1. Generate 15 structures, set max_parallel=5 → relaxes structures 0-4
    2. After completion, set max_parallel=10 → structures 0-4 already done, relaxes 5-9
    3. Repeat until all structures are processed

    Args:
        structures: Dict of structures to relax (e.g., {'structure_0': StructureData, ...})
        code_pk: PK of AiiDA Code for VASP (as int)
        vasp_config: VASP configuration as plain dict containing:
            - relax: Dict with relaxation settings (positions, shape, volume, etc.)
            - base: Dict with base VASP settings (force_cutoff, etc.)
            - kpoints_distance: Float (Angstrom, typical: 0.3-0.5)
            - potential_family: Str (e.g., 'PBE.54')
            - potential_mapping: Dict (optional, element->potential mapping)
            - clean_workdir: Bool (default: False)
        options: Scheduler options as plain dict (resources, queue, time, etc.)
        max_parallel: Maximum number of structures to process (limits to first N structures)

    Returns:
        Dictionary with namespace outputs:
            - structures: Dict {idx: StructureData} for successful relaxations
            - energies: Dict {idx: Float} for successful relaxations
            - exit_statuses: Dict {idx: Int} for all relaxations
            - errors: Dict {idx: Str} for all relaxations

        Each dict uses string index matching input (e.g., '0', '1', '2', ...)
        Only processed structures appear in outputs.

    Note:
        Uses vasp.v2.relax workflow plugin (not vasp.v2.vasp).
        All calculations run in parallel (scatter-gather pattern).
    """
    # Get VASP relax workchain
    VaspRelaxWorkChain = WorkflowFactory('vasp.v2.relax')
    VaspTask = task(VaspRelaxWorkChain)

    # Load code from PK
    code = orm.load_node(code_pk)

    # Extract config from plain dict
    relax = vasp_config.get('relax', {})
    base = vasp_config.get('base', {})
    kpoints_distance = vasp_config.get('kpoints_distance', 0.5)
    potential_family = vasp_config.get('potential_family', 'PBE')
    potential_mapping = vasp_config.get('potential_mapping', {})
    clean_workdir = vasp_config.get('clean_workdir', False)

    # Output dictionaries
    structures_out = {}
    energies_out = {}
    exit_statuses_out = {}
    errors_out = {}

    # SIMPLE BATCH APPROACH: Only process first max_parallel structures
    # Filter for structure keys only (ignore 'manifest' and other keys)
    structure_keys = [k for k in structures.keys() if k.startswith('structure_')]
    # Sort keys to ensure consistent ordering (structure_0, structure_1, ...)
    sorted_keys = sorted(structure_keys, key=lambda k: int(k.split('_')[1]))
    selected_keys = sorted_keys[:max_parallel]

    # Process selected structures
    for key in selected_keys:
        structure = structures[key]

        # Extract index from key (e.g., 'structure_0' -> 0)
        idx = int(key.split('_')[1])

        # Build VASP relax workflow inputs
        vasp_inputs = {
            'structure': structure,
            'code': code,
            'relax': relax,
            'base': base,
            'kpoints_distance': kpoints_distance,
            'potential_family': potential_family,
            'potential_mapping': potential_mapping,
            'options': options,
            'clean_workdir': clean_workdir,
        }

        # Create VASP relaxation task
        # All tasks created in this loop will run in parallel (scatter-gather pattern)
        vasp_task = VaspTask(**vasp_inputs)

        # Extract energy from VASP misc output
        energy = extract_total_energy(energies=vasp_task.misc)

        # Get exit information
        exit_info = get_exit_info(
            exit_status=vasp_task.exit_status,
            exit_message=vasp_task.exit_message
        )

        # Map to output dictionaries using string index (AiiDA Dict requires string keys)
        idx_key = str(idx)
        structures_out[idx_key] = vasp_task.structure
        energies_out[idx_key] = energy.result
        exit_statuses_out[idx_key] = exit_info.exit_status
        errors_out[idx_key] = exit_info.error

    return {
        'structures': structures_out,
        'energies': energies_out,
        'exit_statuses': exit_statuses_out,
        'errors': errors_out,
    }
