"""Thickness convergence brick for the lego module.

Handles slab thickness convergence stages: generates slabs at multiple
thicknesses, relaxes them in parallel, computes surface energies, and
determines the converged slab thickness.

Reuses calcfunctions from teros.core.convergence to avoid code duplication.
"""

import typing as t

from aiida import orm
from aiida.common.links import LinkType
from aiida_workgraph import task, dynamic


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate a thickness convergence stage configuration.

    Args:
        stage: Stage configuration dict.
        stage_names: Set of stage names defined so far (before this stage).

    Raises:
        ValueError: If validation fails.
    """
    name = stage['name']

    # structure_from is required and must reference a previous VASP stage
    if 'structure_from' not in stage:
        raise ValueError(
            f"Stage '{name}': thickness stages require 'structure_from' "
            f"(name of a previous VASP stage providing bulk structure and energy)"
        )
    structure_from = stage['structure_from']
    if structure_from not in stage_names:
        raise ValueError(
            f"Stage '{name}' structure_from='{structure_from}' must reference "
            f"a previous stage name"
        )

    # miller_indices is required
    if 'miller_indices' not in stage:
        raise ValueError(
            f"Stage '{name}': thickness stages require 'miller_indices' "
            f"(e.g., [1, 1, 1])"
        )
    miller = stage['miller_indices']
    if not isinstance(miller, (list, tuple)) or len(miller) != 3:
        raise ValueError(
            f"Stage '{name}': miller_indices must be a list of 3 integers, "
            f"got: {miller}"
        )

    # layer_counts is required with at least 2 values
    if 'layer_counts' not in stage:
        raise ValueError(
            f"Stage '{name}': thickness stages require 'layer_counts' "
            f"(e.g., [3, 5, 7, 9, 11])"
        )
    layers = stage['layer_counts']
    if not isinstance(layers, (list, tuple)) or len(layers) < 2:
        raise ValueError(
            f"Stage '{name}': layer_counts must have at least 2 values, "
            f"got: {layers}"
        )

    # slab_incar is required
    if 'slab_incar' not in stage:
        raise ValueError(
            f"Stage '{name}': thickness stages require 'slab_incar' "
            f"(INCAR parameters for slab relaxations)"
        )


def create_stage_tasks(wg, stage, stage_name, context):
    """Create thickness convergence stage tasks in the WorkGraph.

    This creates a sub-workflow within the sequential WorkGraph:
    1. Generate slabs at multiple thicknesses (from bulk structure)
    2. Relax all slabs in parallel
    3. Compute surface energy for each thickness
    4. Analyze convergence and determine recommended thickness
    5. Select the recommended slab structure

    Args:
        wg: WorkGraph to add tasks to.
        stage: Stage configuration dict.
        stage_name: Unique stage identifier.
        context: Dict with shared context.

    Returns:
        Dict with task references for later stages.
    """
    from teros.core.convergence import (
        generate_thickness_series,
        relax_thickness_series,
        compute_surface_energies,
        gather_surface_energies,
    )

    stage_tasks = context['stage_tasks']
    stage_types = context['stage_types']
    code = context['code']
    potential_family = context['potential_family']
    potential_mapping = context['potential_mapping']
    options = context['options']
    base_kpoints_spacing = context['base_kpoints_spacing']
    clean_workdir = context['clean_workdir']

    # --- Resolve bulk structure and energy from referenced VASP stage ---
    structure_from = stage['structure_from']
    ref_stage_type = stage_types.get(structure_from, 'vasp')
    if ref_stage_type != 'vasp':
        raise ValueError(
            f"Thickness stage '{stage_name}' structure_from='{structure_from}' "
            f"must reference a VASP stage (got type='{ref_stage_type}')"
        )

    bulk_structure = stage_tasks[structure_from]['vasp'].outputs.structure
    bulk_energy = stage_tasks[structure_from]['energy'].outputs.result

    # --- Stage configuration ---
    miller_indices = stage['miller_indices']
    layer_counts = stage['layer_counts']
    slab_incar = stage['slab_incar']
    convergence_threshold = stage.get('convergence_threshold', 0.01)
    min_vacuum_thickness = stage.get('min_vacuum_thickness', 15.0)
    termination_index = stage.get('termination_index', 0)
    lll_reduce = stage.get('lll_reduce', True)
    center_slab = stage.get('center_slab', True)
    primitive = stage.get('primitive', True)
    slab_kpoints_spacing = stage.get('slab_kpoints_spacing', base_kpoints_spacing)
    max_concurrent_jobs = stage.get('max_concurrent_jobs', None)

    # --- 1. Generate slabs at multiple thicknesses ---
    miller_list = orm.List(list=[int(m) for m in miller_indices])
    layers_list = orm.List(list=[int(n) for n in layer_counts])

    slab_gen_task = wg.add_task(
        generate_thickness_series,
        name=f'generate_slabs_{stage_name}',
        bulk_structure=bulk_structure,
        miller_indices=miller_list,
        layer_counts=layers_list,
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        primitive=orm.Bool(primitive),
        termination_index=orm.Int(termination_index),
    )

    # --- 2. Relax all slabs in parallel ---
    relax_task = wg.add_task(
        relax_thickness_series,
        name=f'relax_slabs_{stage_name}',
        slabs=slab_gen_task.outputs.slabs,
        code_pk=code.pk,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=slab_incar,
        options=options,
        kpoints_spacing=slab_kpoints_spacing,
        clean_workdir=clean_workdir,
        max_number_jobs=max_concurrent_jobs,
    )

    # --- 3. Compute surface energies ---
    surface_energy_task = wg.add_task(
        compute_surface_energies,
        name=f'surface_energies_{stage_name}',
        slabs=relax_task.outputs.relaxed_structures,
        energies=relax_task.outputs.energies,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
    )

    # --- 4. Gather results and analyze convergence ---
    gather_task = wg.add_task(
        gather_surface_energies,
        name=f'gather_{stage_name}',
        surface_energies=surface_energy_task.outputs.surface_energies,
        miller_indices=miller_list,
        convergence_threshold=orm.Float(convergence_threshold),
    )

    # --- 5. Select recommended slab structure ---
    select_task = wg.add_task(
        select_recommended_slab,
        name=f'select_{stage_name}',
        convergence_results=gather_task.outputs.result,
        relaxed_structures=relax_task.outputs.relaxed_structures,
    )

    return {
        'generate_slabs': slab_gen_task,
        'relax_slabs': relax_task,
        'surface_energies': surface_energy_task,
        'gather': gather_task,
        'select': select_task,
        'structure': select_task.outputs.recommended_structure,
    }


def expose_stage_outputs(wg, stage_name, stage_tasks_result):
    """Expose thickness convergence stage outputs on the WorkGraph.

    Exposes:
        - {stage_name}_convergence_results: Dict with all slab data,
          surface energies, and convergence analysis (recommended thickness)
        - {stage_name}_recommended_structure: StructureData of the slab
          at the recommended (converged) thickness

    Args:
        wg: WorkGraph instance.
        stage_name: Unique stage identifier.
        stage_tasks_result: Dict returned by create_stage_tasks.
    """
    gather_task = stage_tasks_result['gather']
    select_task = stage_tasks_result['select']

    setattr(
        wg.outputs,
        f'{stage_name}_convergence_results',
        gather_task.outputs.result,
    )
    setattr(
        wg.outputs,
        f'{stage_name}_recommended_structure',
        select_task.outputs.recommended_structure,
    )


def get_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a thickness convergence stage.

    Args:
        wg_node: The WorkGraph ProcessNode.
        wg_pk: WorkGraph PK.
        stage_name: Name of the thickness stage.

    Returns:
        Dict with keys: convergence_results, recommended_structure,
        recommended_layers, converged, surface_energies, pk, stage, type.
    """
    result = {
        'convergence_results': None,
        'recommended_structure': None,
        'recommended_layers': None,
        'converged': False,
        'surface_energies': {},
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'thickness',
    }

    # Try direct output access
    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        conv_attr = f'{stage_name}_convergence_results'
        if hasattr(outputs, conv_attr):
            conv_node = getattr(outputs, conv_attr)
            if hasattr(conv_node, 'get_dict'):
                conv_data = conv_node.get_dict()
                result['convergence_results'] = conv_data
                summary = conv_data.get('summary', {})
                result['recommended_layers'] = summary.get('recommended_layers')
                result['converged'] = summary.get('converged', False)
                thicknesses = summary.get('thicknesses', [])
                energies = summary.get('surface_energies_J_m2', [])
                result['surface_energies'] = dict(zip(thicknesses, energies))

        struct_attr = f'{stage_name}_recommended_structure'
        if hasattr(outputs, struct_attr):
            result['recommended_structure'] = getattr(outputs, struct_attr)

    # Fallback: traverse links
    if result['convergence_results'] is None:
        _extract_thickness_stage_from_workgraph(wg_node, stage_name, result)

    return result


def _extract_thickness_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """Extract thickness stage results by traversing WorkGraph links.

    Args:
        wg_node: The WorkGraph ProcessNode.
        stage_name: Name of the thickness stage.
        result: Result dict to populate (modified in place).
    """
    if not hasattr(wg_node, 'base'):
        return

    gather_task_name = f'gather_{stage_name}'
    select_task_name = f'select_{stage_name}'

    # Look for gather calcfunction (convergence results)
    called_calc = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
    for link in called_calc.all():
        link_label = link.link_label
        child_node = link.node

        if gather_task_name in link_label:
            created = child_node.base.links.get_outgoing(
                link_type=LinkType.CREATE
            )
            for out_link in created.all():
                if out_link.link_label == 'result':
                    conv_node = out_link.node
                    if hasattr(conv_node, 'get_dict'):
                        conv_data = conv_node.get_dict()
                        result['convergence_results'] = conv_data
                        summary = conv_data.get('summary', {})
                        result['recommended_layers'] = summary.get(
                            'recommended_layers'
                        )
                        result['converged'] = summary.get('converged', False)
                        thicknesses = summary.get('thicknesses', [])
                        energies = summary.get('surface_energies_J_m2', [])
                        result['surface_energies'] = dict(
                            zip(thicknesses, energies)
                        )
            break

    # Look for select calcfunction (recommended structure)
    if result['recommended_structure'] is None:
        for link in called_calc.all():
            link_label = link.link_label
            child_node = link.node

            if select_task_name in link_label:
                created = child_node.base.links.get_outgoing(
                    link_type=LinkType.CREATE
                )
                for out_link in created.all():
                    if out_link.link_label == 'result':
                        result['recommended_structure'] = out_link.node
                break

    # Also check CALL_WORK links for graph tasks
    if result['convergence_results'] is None or result['recommended_structure'] is None:
        called_work = wg_node.base.links.get_outgoing(
            link_type=LinkType.CALL_WORK
        )
        for link in called_work.all():
            link_label = link.link_label
            child_node = link.node

            if gather_task_name in link_label and result['convergence_results'] is None:
                child_calcs = child_node.base.links.get_outgoing(
                    link_type=LinkType.CALL_CALC
                )
                for calc_link in child_calcs.all():
                    created = calc_link.node.base.links.get_outgoing(
                        link_type=LinkType.CREATE
                    )
                    for out_link in created.all():
                        if out_link.link_label == 'result' and hasattr(
                            out_link.node, 'get_dict'
                        ):
                            conv_data = out_link.node.get_dict()
                            result['convergence_results'] = conv_data
                            summary = conv_data.get('summary', {})
                            result['recommended_layers'] = summary.get(
                                'recommended_layers'
                            )
                            result['converged'] = summary.get(
                                'converged', False
                            )
                            thicknesses = summary.get('thicknesses', [])
                            energies = summary.get(
                                'surface_energies_J_m2', []
                            )
                            result['surface_energies'] = dict(
                                zip(thicknesses, energies)
                            )

            if select_task_name in link_label and result['recommended_structure'] is None:
                child_calcs = child_node.base.links.get_outgoing(
                    link_type=LinkType.CALL_CALC
                )
                for calc_link in child_calcs.all():
                    created = calc_link.node.base.links.get_outgoing(
                        link_type=LinkType.CREATE
                    )
                    for out_link in created.all():
                        if out_link.link_label == 'result':
                            result['recommended_structure'] = out_link.node


def print_stage_results(index: int, stage_name: str, stage_result: dict) -> None:
    """Print formatted results for a thickness convergence stage.

    Args:
        index: 1-based stage index.
        stage_name: Name of the stage.
        stage_result: Result dict from get_stage_results.
    """
    print(f"  [{index}] {stage_name} (THICKNESS CONVERGENCE)")

    conv = stage_result.get('convergence_results')
    if conv is not None:
        summary = conv.get('summary', {})
        miller = conv.get('miller_indices', [])
        threshold = summary.get('convergence_threshold', 0.01)
        print(f"      Miller indices: ({', '.join(str(m) for m in miller)})")
        print(f"      Threshold: {threshold} J/m\u00b2")

        thicknesses = summary.get('thicknesses', [])
        energies = summary.get('surface_energies_J_m2', [])

        if thicknesses and energies:
            print("      Layers  \u03b3 (J/m\u00b2)  \u0394\u03b3 (mJ/m\u00b2)")
            print("      ------  ---------  ----------")
            for i, (n, gamma) in enumerate(zip(thicknesses, energies)):
                if i == 0:
                    delta_str = "       -"
                else:
                    delta = abs(energies[i] - energies[i - 1]) * 1000
                    delta_str = f"{delta:8.1f}"
                print(f"      {n:>6d}  {gamma:9.4f}  {delta_str}")

        converged = summary.get('converged', False)
        recommended = summary.get('recommended_layers')
        if converged:
            print(f"      Converged: YES at {recommended} layers")
        else:
            max_tested = summary.get('max_tested_layers', '?')
            print(f"      Converged: NO (tested up to {max_tested} layers)")

    struct = stage_result.get('recommended_structure')
    if struct is not None and hasattr(struct, 'get_formula'):
        formula = struct.get_formula()
        n_atoms = len(struct.sites)
        print(
            f"      Recommended structure: {formula} "
            f"({n_atoms} atoms, PK: {struct.pk})"
        )


# ─── Calcfunction tasks ─────────────────────────────────────────────────────


@task.graph(outputs=['recommended_structure'])
def select_recommended_slab(
    convergence_results: orm.Dict,
    relaxed_structures: t.Annotated[
        dict[str, orm.StructureData], dynamic(orm.StructureData)
    ],
) -> dict:
    """Select the relaxed slab structure at the recommended thickness.

    This graph task iterates over the dynamic namespace of relaxed structures
    and passes them to a calcfunction that picks the one matching the
    recommended layer count from the convergence analysis.

    Args:
        convergence_results: Dict with convergence analysis (from
            gather_surface_energies), containing summary.recommended_layers.
        relaxed_structures: Dynamic namespace of relaxed slab StructureData
            keyed by layer count (e.g., layers_3, layers_5, ...).

    Returns:
        Dict with 'recommended_structure' key pointing to the selected
        StructureData.
    """
    kwargs = {'convergence_results': convergence_results}
    for key, struct in relaxed_structures.items():
        kwargs[key] = struct

    picked = _pick_recommended_structure(**kwargs)
    return {'recommended_structure': picked.result}


@task.calcfunction
def _pick_recommended_structure(
    convergence_results: orm.Dict, **structures
) -> orm.StructureData:
    """Pick the relaxed structure at the recommended layer count.

    Args:
        convergence_results: Dict with convergence analysis containing
            summary.recommended_layers.
        **structures: Relaxed slab structures keyed by 'layers_N'.

    Returns:
        StructureData at the recommended thickness.
    """
    conv_data = convergence_results.get_dict()
    summary = conv_data.get('summary', {})
    recommended = summary.get('recommended_layers')

    if recommended is None:
        # If convergence was not reached, return the thickest slab
        available_keys = [k for k in structures if k.startswith('layers_')]
        if not available_keys:
            raise ValueError(
                "No slab structures found in inputs. "
                f"Available keys: {list(structures.keys())}"
            )
        layer_counts = [int(k.split('_')[1]) for k in available_keys]
        recommended = max(layer_counts)

    key = f'layers_{recommended}'
    if key not in structures:
        available = [k for k in structures if k.startswith('layers_')]
        raise ValueError(
            f"Recommended key '{key}' not found in structures. "
            f"Available: {available}"
        )

    # Return a copy: returning an input node directly would create a RETURN
    # link instead of CREATE, which breaks provenance link traversal.
    return orm.StructureData(ase=structures[key].get_ase())
