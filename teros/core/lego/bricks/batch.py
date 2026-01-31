"""Batch brick for the lego module.

Handles batch stages in two modes:

1. Single-structure mode (structure_from + calculations):
   One structure × N INCAR variations. Each entry in 'calculations'
   produces a separate VASP calculation on the same structure.

2. Multi-structure mode (structures_from + base_incar):
   N structures × 1 INCAR. Takes a namespace of structures from a
   previous stage (e.g., slab_gen) and runs the same VASP calculation
   on each. Returns relaxed_structures and energies namespaces.
"""

import typing as t

from aiida import orm
from aiida.common.links import LinkType
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace, get_current_graph

from ..tasks import extract_energy
from ...utils import deep_merge_dicts, get_vasp_parser_settings, extract_max_jobs_value


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate a batch stage configuration.

    Supports two modes:
    - Single-structure: requires structure_from + base_incar + calculations
    - Multi-structure: requires structures_from + base_incar

    Args:
        stage: Stage configuration dict.
        stage_names: Set of stage names defined so far (before this stage).

    Raises:
        ValueError: If validation fails.
    """
    name = stage['name']

    has_structure_from = 'structure_from' in stage
    has_structures_from = 'structures_from' in stage

    if not has_structure_from and not has_structures_from:
        raise ValueError(
            f"Stage '{name}': batch stages require 'structure_from' "
            f"or 'structures_from'"
        )

    if has_structure_from and has_structures_from:
        raise ValueError(
            f"Stage '{name}': batch stages cannot have both 'structure_from' "
            f"and 'structures_from'"
        )

    if 'base_incar' not in stage:
        raise ValueError(f"Stage '{name}': batch stages require 'base_incar'")

    if has_structure_from:
        # Single-structure mode: need calculations dict
        if 'calculations' not in stage or not stage['calculations']:
            raise ValueError(
                f"Stage '{name}': single-structure batch stages require "
                f"non-empty 'calculations' dict"
            )
        structure_from = stage['structure_from']
        if structure_from not in stage_names:
            raise ValueError(
                f"Stage '{name}' structure_from='{structure_from}' must "
                f"reference a previous stage name"
            )

    if has_structures_from:
        # Multi-structure mode: calculations not required
        structures_from = stage['structures_from']
        if structures_from not in stage_names:
            raise ValueError(
                f"Stage '{name}' structures_from='{structures_from}' must "
                f"reference a previous stage name"
            )


def create_stage_tasks(wg, stage, stage_name, context):
    """Create batch stage tasks in the WorkGraph.

    Dispatches to single-structure or multi-structure mode.

    Args:
        wg: WorkGraph to add tasks to.
        stage: Stage configuration dict.
        stage_name: Unique stage identifier.
        context: Dict with shared context.

    Returns:
        Dict with task references for later stages.
    """
    if 'structures_from' in stage:
        return _create_multi_structure_tasks(wg, stage, stage_name, context)
    else:
        return _create_single_structure_tasks(wg, stage, stage_name, context)


def _create_single_structure_tasks(wg, stage, stage_name, context):
    """Create single-structure batch tasks (original behavior).

    One structure × N INCAR variations from 'calculations' dict.
    """
    from . import resolve_structure_from
    from ..workgraph import _prepare_builder_inputs

    code = context['code']
    potential_family = context['potential_family']
    potential_mapping = context['potential_mapping']
    options = context['options']
    base_kpoints_spacing = context['base_kpoints_spacing']
    clean_workdir = context['clean_workdir']

    # Resolve structure from referenced stage
    structure_from = stage['structure_from']
    input_structure = resolve_structure_from(structure_from, context)

    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    base_incar = stage['base_incar']
    calculations = stage['calculations']

    # Stage-level defaults
    stage_kpoints_spacing = stage.get('kpoints_spacing', base_kpoints_spacing)
    stage_kpoints_mesh = stage.get('kpoints', None)
    stage_retrieve = stage.get('retrieve', None)

    calc_tasks = {}
    energy_tasks = {}

    for calc_label, calc_config in calculations.items():
        # Deep-merge base_incar with per-calculation incar overrides
        calc_incar_overrides = calc_config.get('incar', {})
        if calc_incar_overrides:
            merged_incar = deep_merge_dicts(base_incar, calc_incar_overrides)
        else:
            merged_incar = dict(base_incar)

        # Per-calc kpoints or fall back to stage-level
        calc_kpoints_mesh = calc_config.get('kpoints', stage_kpoints_mesh)
        calc_kpoints_spacing = calc_config.get('kpoints_spacing', stage_kpoints_spacing)

        # Per-calc retrieve or fall back to stage-level
        calc_retrieve = calc_config.get('retrieve', stage_retrieve)

        # Prepare builder inputs
        builder_inputs = _prepare_builder_inputs(
            incar=merged_incar,
            kpoints_spacing=calc_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            retrieve=calc_retrieve,
            restart_folder=None,
            clean_workdir=clean_workdir,
            kpoints_mesh=calc_kpoints_mesh,
        )

        # Add VASP task
        vasp_task_name = f'vasp_{stage_name}_{calc_label}'
        vasp_task = wg.add_task(
            VaspTask,
            name=vasp_task_name,
            structure=input_structure,
            code=code,
            **builder_inputs
        )

        # Add energy extraction task
        energy_task_name = f'energy_{stage_name}_{calc_label}'
        energy_task = wg.add_task(
            extract_energy,
            name=energy_task_name,
            misc=vasp_task.outputs.misc,
            retrieved=vasp_task.outputs.retrieved,
        )

        calc_tasks[calc_label] = vasp_task
        energy_tasks[calc_label] = energy_task

    return {
        'calc_tasks': calc_tasks,
        'energy_tasks': energy_tasks,
        'structure': input_structure,
        'mode': 'single',
    }


def _create_multi_structure_tasks(wg, stage, stage_name, context):
    """Create multi-structure batch tasks.

    N structures × 1 INCAR. Uses a @task.graph function to iterate
    over the dynamic namespace of structures at runtime.
    """
    stage_tasks = context['stage_tasks']
    stage_types = context['stage_types']
    code = context['code']
    potential_family = context['potential_family']
    potential_mapping = context['potential_mapping']
    options = context['options']
    base_kpoints_spacing = context['base_kpoints_spacing']
    clean_workdir = context['clean_workdir']

    # Resolve structures namespace from referenced stage
    structures_from = stage['structures_from']
    ref_stage_type = stage_types.get(structures_from, 'vasp')

    if ref_stage_type == 'slab_gen':
        input_structures = stage_tasks[structures_from]['slabs']
    else:
        raise ValueError(
            f"Batch stage '{stage_name}' structures_from='{structures_from}' "
            f"must reference a slab_gen stage (got type='{ref_stage_type}')"
        )

    base_incar = stage['base_incar']
    kpoints_spacing = stage.get('kpoints_spacing', base_kpoints_spacing)
    max_concurrent_jobs = stage.get('max_concurrent_jobs', None)

    batch_graph_task = wg.add_task(
        _run_multi_structure_batch,
        name=f'batch_{stage_name}',
        structures=input_structures,
        code_pk=code.pk,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        incar_params=base_incar,
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
        max_concurrent_jobs=max_concurrent_jobs,
    )

    return {
        'batch_graph': batch_graph_task,
        'relaxed_structures': batch_graph_task.outputs.relaxed_structures,
        'energies': batch_graph_task.outputs.energies,
        'structure': input_structures,
        'mode': 'multi',
    }


def expose_stage_outputs(wg, stage_name, stage_tasks_result):
    """Expose batch stage outputs on the WorkGraph.

    Args:
        wg: WorkGraph instance.
        stage_name: Unique stage identifier.
        stage_tasks_result: Dict returned by create_stage_tasks.
    """
    mode = stage_tasks_result.get('mode', 'single')

    if mode == 'single':
        for calc_label, vasp_task in stage_tasks_result['calc_tasks'].items():
            energy_task = stage_tasks_result['energy_tasks'][calc_label]
            setattr(wg.outputs, f'{stage_name}_{calc_label}_energy',
                    energy_task.outputs.result)
            setattr(wg.outputs, f'{stage_name}_{calc_label}_misc',
                    vasp_task.outputs.misc)
            setattr(wg.outputs, f'{stage_name}_{calc_label}_remote',
                    vasp_task.outputs.remote_folder)
            setattr(wg.outputs, f'{stage_name}_{calc_label}_retrieved',
                    vasp_task.outputs.retrieved)
    else:
        # Multi-structure mode: no individual outputs exposed.
        # Downstream stages (e.g., thickness) consume the namespace sockets.
        pass


def get_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a batch stage in a sequential workflow.

    Args:
        wg_node: The WorkGraph ProcessNode.
        wg_pk: WorkGraph PK.
        stage_name: Name of the batch stage.

    Returns:
        Dict with keys: calculations, pk, stage, type.
    """
    from ..results import _extract_energy_from_misc

    result = {
        'calculations': {},
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'batch',
    }

    energy_suffix = '_energy'
    stage_prefix = f'{stage_name}_'

    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        # Find all output attributes matching the pattern
        calc_labels = []
        for attr_name in dir(outputs):
            if attr_name.startswith(stage_prefix) and attr_name.endswith(energy_suffix):
                calc_label = attr_name[len(stage_prefix):-len(energy_suffix)]
                if calc_label:
                    calc_labels.append(calc_label)

        for calc_label in calc_labels:
            calc_result = {
                'energy': None,
                'misc': None,
                'remote': None,
                'files': None,
            }

            # Energy
            energy_attr = f'{stage_name}_{calc_label}_energy'
            if hasattr(outputs, energy_attr):
                energy_node = getattr(outputs, energy_attr)
                if hasattr(energy_node, 'value'):
                    calc_result['energy'] = energy_node.value
                else:
                    calc_result['energy'] = float(energy_node)

            # Misc
            misc_attr = f'{stage_name}_{calc_label}_misc'
            if hasattr(outputs, misc_attr):
                misc_node = getattr(outputs, misc_attr)
                if hasattr(misc_node, 'get_dict'):
                    calc_result['misc'] = misc_node.get_dict()

            # Remote folder
            remote_attr = f'{stage_name}_{calc_label}_remote'
            if hasattr(outputs, remote_attr):
                calc_result['remote'] = getattr(outputs, remote_attr)

            # Retrieved files
            retrieved_attr = f'{stage_name}_{calc_label}_retrieved'
            if hasattr(outputs, retrieved_attr):
                calc_result['files'] = getattr(outputs, retrieved_attr)

            # Extract energy from misc if not found directly
            if calc_result['energy'] is None and calc_result['misc'] is not None:
                calc_result['energy'] = _extract_energy_from_misc(calc_result['misc'])

            result['calculations'][calc_label] = calc_result

    # Fallback: traverse links if outputs not found
    if not result['calculations']:
        _extract_batch_stage_from_workgraph(wg_node, stage_name, result)

    return result


def _extract_batch_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """Extract batch stage results by traversing WorkGraph links.

    Args:
        wg_node: The WorkGraph ProcessNode.
        stage_name: Name of the batch stage.
        result: Result dict to populate (modified in place).
    """
    from ..results import _extract_energy_from_misc

    if not hasattr(wg_node, 'base'):
        return

    vasp_prefix = f'vasp_{stage_name}_'

    # Collect calc labels from CALL_WORK links
    calc_labels = set()
    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        link_label = link.link_label
        if link_label.startswith(vasp_prefix):
            calc_label = link_label[len(vasp_prefix):]
            calc_labels.add(calc_label)

    for calc_label in calc_labels:
        calc_result = {
            'energy': None,
            'misc': None,
            'remote': None,
            'files': None,
        }

        vasp_task_name = f'vasp_{stage_name}_{calc_label}'
        energy_task_name = f'energy_{stage_name}_{calc_label}'

        # Find VASP task outputs
        for link in called.all():
            if link.link_label == vasp_task_name or vasp_task_name in link.link_label:
                child_node = link.node
                if hasattr(child_node, 'outputs'):
                    outputs = child_node.outputs
                    if calc_result['misc'] is None and hasattr(outputs, 'misc'):
                        misc = outputs.misc
                        if hasattr(misc, 'get_dict'):
                            calc_result['misc'] = misc.get_dict()
                    if calc_result['remote'] is None and hasattr(outputs, 'remote_folder'):
                        calc_result['remote'] = outputs.remote_folder
                    if calc_result['files'] is None and hasattr(outputs, 'retrieved'):
                        calc_result['files'] = outputs.retrieved

        # Find energy task outputs
        called_calc = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
        for link in called_calc.all():
            if link.link_label == energy_task_name or energy_task_name in link.link_label:
                created = link.node.base.links.get_outgoing(link_type=LinkType.CREATE)
                for out_link in created.all():
                    if out_link.link_label == 'result':
                        energy_node = out_link.node
                        if hasattr(energy_node, 'value'):
                            calc_result['energy'] = energy_node.value
                        break

        # Extract energy from misc if not found
        if calc_result['energy'] is None and calc_result['misc'] is not None:
            calc_result['energy'] = _extract_energy_from_misc(calc_result['misc'])

        result['calculations'][calc_label] = calc_result


def print_stage_results(index: int, stage_name: str, stage_result: dict) -> None:
    """Print formatted results for a batch stage.

    Args:
        index: 1-based stage index.
        stage_name: Name of the stage.
        stage_result: Result dict from get_stage_results.
    """
    print(f"  [{index}] {stage_name} (BATCH)")

    calculations = stage_result.get('calculations', {})
    if calculations:
        for calc_label, calc_result in calculations.items():
            energy_str = (
                f"{calc_result['energy']:.6f} eV"
                if calc_result['energy'] is not None
                else "N/A"
            )
            print(f"      [{calc_label}] Energy: {energy_str}")

            if calc_result.get('misc') is not None:
                misc = calc_result['misc']
                run_status = misc.get('run_status', 'N/A')
                print(f"        Status: {run_status}")

            if calc_result.get('remote') is not None:
                print(f"        Remote folder: PK {calc_result['remote'].pk}")

            if calc_result.get('files') is not None:
                files = calc_result['files'].list_object_names()
                print(f"        Retrieved: {', '.join(files)}")
    else:
        print("      (No calculation results found)")


# ─── Multi-structure graph task ────────────────────────────────────────────


@task.graph
def _run_multi_structure_batch(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code_pk: int,
    potential_family: str,
    potential_mapping: dict,
    incar_params: dict,
    options: dict,
    kpoints_spacing: float = 0.03,
    clean_workdir: bool = True,
    max_concurrent_jobs: int = None,
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
)]:
    """Run the same VASP calculation on multiple structures in parallel.

    This graph task iterates over a dynamic namespace of structures
    (e.g., from slab generation) and runs a VaspWorkChain on each.
    Concurrency is controlled via max_concurrent_jobs.

    Args:
        structures: Dynamic namespace of input structures keyed by label.
        code_pk: PK of AiiDA code for VASP.
        potential_family: Pseudopotential family.
        potential_mapping: Element to potential mapping.
        incar_params: INCAR parameters applied to all calculations.
        options: Scheduler options.
        kpoints_spacing: K-points spacing in A^-1.
        clean_workdir: Whether to clean remote directories.
        max_concurrent_jobs: Maximum concurrent VASP calculations.

    Returns:
        Dict with relaxed_structures and energies namespaces.
    """
    # Set max_concurrent_jobs on this workgraph to control concurrency
    if max_concurrent_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = extract_max_jobs_value(max_concurrent_jobs)
        wg.max_number_jobs = max_jobs_value

    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    code = orm.load_node(code_pk)
    settings = orm.Dict(dict=get_vasp_parser_settings(add_energy=True))

    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}

    for label, structure in structures.items():
        vasp_inputs = {
            'structure': structure,
            'code': code,
            'parameters': orm.Dict(dict={'incar': dict(incar_params)}),
            'options': orm.Dict(dict=dict(options)),
            'potential_family': potential_family,
            'potential_mapping': orm.Dict(dict=dict(potential_mapping)),
            'clean_workdir': clean_workdir,
            'settings': settings,
        }

        kpts_val = kpoints_spacing
        if hasattr(kpoints_spacing, 'value'):
            kpts_val = kpoints_spacing.value
        vasp_inputs['kpoints_spacing'] = float(kpts_val)

        vasp_calc = VaspTask(**vasp_inputs)
        relaxed[label] = vasp_calc.structure
        energies_ns[label] = extract_energy(
            misc=vasp_calc.misc,
            retrieved=vasp_calc.retrieved,
        ).result

    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
    }
