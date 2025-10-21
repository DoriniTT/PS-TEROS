"""Child WorkGraph for parallel VASP relaxations with batch limiting."""

import typing as t
from aiida import orm
from aiida.engine import calcfunction
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic


def get_settings():
    """Parser settings for aiida-vasp."""
    return {
        'parser_settings': {
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }


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


@task.graph
def relax_slabs_with_semaphore(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code_pk: int,
    vasp_config: dict,
    options: dict,
    max_parallel: int,
) -> t.Annotated[dict, namespace(**{
    'structures': dynamic(orm.StructureData),
    'energies': dynamic(orm.Float),
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
            - parameters: Dict with INCAR parameters (ENCUT, EDIFF, ISIF, NSW, IBRION, EDIFFG, etc.)
            - kpoints_spacing: Float (Angstrom, typical: 0.3-0.5)
            - potential_family: Str (e.g., 'PBE', 'PBE.54')
            - potential_mapping: Dict (optional, element->potential mapping)
            - clean_workdir: Bool (default: False)
        options: Scheduler options as plain dict (resources, queue, time, etc.)
        max_parallel: Maximum number of structures to process (limits to first N structures)

    Returns:
        Dictionary with namespace outputs:
            - structures: Dict {key: StructureData} for successful relaxations
            - energies: Dict {key: Float} for successful relaxations

        Each dict uses descriptive keys from input (e.g., '0_oh_000_3_7572', '1_oh_001_7_5145', ...)
        Note: Dots are replaced with underscores for AiiDA link label compatibility.
        Only successfully relaxed structures appear in outputs.
        Failed relaxations will not have entries in the output dicts.

    Note:
        Uses vasp.v2.vasp workflow plugin.
        All calculations run in parallel (scatter-gather pattern).
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Load code from PK
    code = orm.load_node(code_pk)

    # Extract VASP configuration
    incar_params = vasp_config.get('parameters', {})
    kpoints_spacing = vasp_config.get('kpoints_spacing', 0.5)
    potential_family = vasp_config.get('potential_family', 'PBE')
    potential_mapping = vasp_config.get('potential_mapping', {})
    clean_workdir = vasp_config.get('clean_workdir', False)

    # Output dictionaries (only successful relaxations will be added)
    structures_out = {}
    energies_out = {}

    # SIMPLE BATCH APPROACH: Only process first max_parallel structures
    # structures is a dict with descriptive keys: '0_oh_000_3.7572', '1_oh_001_7.5145', etc.
    # Sort keys by extracting the numeric prefix to ensure consistent ordering
    sorted_keys = sorted(structures.keys(), key=lambda k: int(k.split('_')[0]))
    # Extract value from max_parallel if it's an AiiDA node
    limit = max_parallel.value if hasattr(max_parallel, 'value') else int(max_parallel)
    selected_keys = sorted_keys[:limit]

    # Process selected structures
    for key in selected_keys:
        structure = structures[key]

        # Extract index from key (e.g., "0" from "0_oh_000_3.7572")
        idx = int(key.split('_')[0])

        # Build VASP inputs following slabs.py pattern
        vasp_inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(incar_params)},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': orm.Dict(dict=potential_mapping) if potential_mapping else orm.Dict(dict={}),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            vasp_inputs['kpoints_spacing'] = kpoints_spacing

        # Create VASP task
        # All tasks created in this loop will run in parallel (scatter-gather pattern)
        vasp_task = VaspTask(**vasp_inputs)

        # Extract energy from VASP misc output
        energy = extract_total_energy(energies=vasp_task.misc)

        # Map to output dictionaries using descriptive key from input
        # Note: Only successful relaxations will populate these outputs
        # Failed calculations will be detected in collect_results by absence of data
        # Use the same descriptive key from input (e.g., "0_oh_000_3.7572")
        structures_out[key] = vasp_task.structure
        energies_out[key] = energy.result

    return {
        'structures': structures_out,
        'energies': energies_out,
    }
