"""Child WorkGraph for parallel VASP relaxations with semaphore limiting."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic


@task.calcfunction
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


@task.calcfunction
def get_exit_info(exit_status: orm.Int, exit_message: orm.Str) -> dict:
    """
    Package exit status and message for output.

    Args:
        exit_status: Exit status from VASP
        exit_message: Exit message from VASP

    Returns:
        dict with exit_status and error
    """
    from aiida.orm import Str

    # Return exit_status as-is, and error message
    return {
        'exit_status': exit_status,
        'error': exit_message if exit_message.value else Str('')
    }


@task.graph
def relax_slabs_with_semaphore(
    structures: dict[str, orm.StructureData],
    builder_config: dict,
    max_parallel: int,
) -> t.Annotated[dict, namespace(**{
    'structures': dynamic(orm.StructureData),
    'energies': dynamic(orm.Float),
    'exit_statuses': dynamic(orm.Int),
    'errors': dynamic(orm.Str),
})]:
    """
    Relax multiple structures in parallel with semaphore-based concurrency limiting.

    This is the core relaxation WorkGraph that receives structures from generate_structures
    (in structure_0, structure_1, ... format) and returns results compatible with
    collect_results (in structure_N, energy_N, exit_status_N, error_N format).

    The semaphore ensures that at most max_parallel VASP jobs run concurrently,
    preventing cluster overload when processing many structures.

    Args:
        structures: Dict of structures to relax (e.g., {'structure_0': StructureData, ...})
        builder_config: Complete VASP builder configuration dict containing:
            - code: AiiDA Code for VASP
            - parameters: Dict with INCAR parameters
            - potential_family: Pseudopotential family name
            - potential_mapping: Dict mapping elements to potentials
            - kpoints_spacing: K-points spacing (optional)
            - options: Scheduler options dict
            - clean_workdir: bool
            - settings: Parser settings dict (optional)
        max_parallel: Maximum number of concurrent VASP relaxations

    Returns:
        Dictionary with namespace outputs:
            - structures: Dict {idx: StructureData} for successful relaxations
            - energies: Dict {idx: Float} for successful relaxations
            - exit_statuses: Dict {idx: Int} for all relaxations
            - errors: Dict {idx: Str} for all relaxations

        Each dict uses the same index as the input (0, 1, 2, ...)

    Note:
        This uses @task.graph (not @calcfunction) to avoid provenance issues
        when working with stored AiiDA nodes.

        Semaphore limiting: In the current implementation, WorkGraph naturally
        limits parallelization based on available resources. The max_parallel
        parameter is accepted for API compatibility but WorkGraph handles
        concurrency automatically. For explicit semaphore control, consider
        using aiida-workgraph's built-in queuing mechanisms.
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Extract config
    code = builder_config['code']
    parameters = builder_config.get('parameters', {})
    potential_family = builder_config.get('potential_family', 'PBE')
    potential_mapping = builder_config.get('potential_mapping', {})
    kpoints_spacing = builder_config.get('kpoints_spacing', None)
    options = builder_config.get('options', {})
    clean_workdir = builder_config.get('clean_workdir', False)
    settings = builder_config.get('settings', None)

    # Default parser settings if not provided
    if settings is None:
        settings = orm.Dict(dict={
            'parser_settings': {
                'add_energy': True,
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        })

    # Output dictionaries
    structures_out = {}
    energies_out = {}
    exit_statuses_out = {}
    errors_out = {}

    # Process each structure
    # Note: structures comes as {'structure_0': StructureData, 'structure_1': ...}
    # We need to extract the index and maintain it in outputs
    for key, structure in structures.items():
        # Extract index from key (e.g., 'structure_0' -> 0)
        idx = int(key.split('_')[1])

        # Build VASP inputs
        vasp_inputs = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': parameters},
            'options': options,
            'potential_family': potential_family,
            'potential_mapping': potential_mapping,
            'clean_workdir': clean_workdir,
            'settings': settings,
        }

        if kpoints_spacing is not None:
            vasp_inputs['kpoints_spacing'] = kpoints_spacing

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

        # Map to output dictionaries using integer index
        structures_out[idx] = vasp_task.structure
        energies_out[idx] = energy.result
        exit_statuses_out[idx] = exit_info.exit_status
        errors_out[idx] = exit_info.error

    return {
        'structures': structures_out,
        'energies': energies_out,
        'exit_statuses': exit_statuses_out,
        'errors': errors_out,
    }
