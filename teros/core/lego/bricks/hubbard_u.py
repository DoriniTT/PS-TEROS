"""Hubbard U brick for the lego module.

Handles Hubbard U parameter calculation stages using the linear response method.
Wraps the existing teros.core.u_calculation module into the lego brick interface.

Reference: https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U
"""

from aiida import orm
from aiida.common.links import LinkType
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

from teros.core.u_calculation.tasks import (
    extract_d_electron_occupation,
    calculate_occupation_response,
    calculate_hubbard_u_linear_regression,
    gather_responses,
    compile_u_calculation_summary,
)
from teros.core.u_calculation.utils import (
    prepare_ground_state_incar,
    prepare_response_incar,
    get_species_order_from_structure,
    DEFAULT_POTENTIAL_VALUES,
)
from teros.core.utils import get_vasp_parser_settings


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate a Hubbard U stage configuration.

    Args:
        stage: Stage configuration dict.
        stage_names: Set of stage names defined so far (before this stage).

    Raises:
        ValueError: If validation fails.
    """
    name = stage['name']

    if 'target_species' not in stage:
        raise ValueError(
            f"Stage '{name}': hubbard_u stages require 'target_species' "
            f"(element symbol for U calculation, e.g., 'Ni', 'Fe')"
        )

    # structure_from is required (can be 'input' or a previous stage name)
    if 'structure_from' not in stage:
        raise ValueError(
            f"Stage '{name}': hubbard_u stages require 'structure_from' "
            f"('input' or a previous stage name)"
        )

    structure_from = stage['structure_from']
    if structure_from != 'input' and structure_from not in stage_names:
        raise ValueError(
            f"Stage '{name}' structure_from='{structure_from}' must be 'input' "
            f"or a previous stage name"
        )

    # Validate potential_values if provided
    potential_values = stage.get('potential_values', None)
    if potential_values is not None:
        if not isinstance(potential_values, (list, tuple)):
            raise ValueError(
                f"Stage '{name}' potential_values must be a list of floats"
            )
        if len(potential_values) < 2:
            raise ValueError(
                f"Stage '{name}' potential_values needs at least 2 values "
                f"for linear regression"
            )
        if 0.0 in potential_values:
            raise ValueError(
                f"Stage '{name}' potential_values must not include 0.0 "
                f"(ground state has LDAU=False, response has LDAU=True)"
            )

    # Validate ldaul if provided
    ldaul = stage.get('ldaul', 2)
    if ldaul not in (2, 3):
        raise ValueError(
            f"Stage '{name}' ldaul={ldaul} must be 2 (d-electrons) or "
            f"3 (f-electrons)"
        )


def create_stage_tasks(wg, stage, stage_name, context):
    """Create Hubbard U stage tasks in the WorkGraph.

    This creates the full linear response workflow:
    1. Ground state calculation (no +U)
    2. For each potential V: NSCF response + SCF response
    3. Gather responses and calculate U via linear regression
    4. Compile summary

    Args:
        wg: WorkGraph to add tasks to.
        stage: Stage configuration dict.
        stage_name: Unique stage identifier.
        context: Dict with shared context.

    Returns:
        Dict with task references for later stages.
    """
    from ..workgraph import _prepare_builder_inputs

    code = context['code']
    potential_family = context['potential_family']
    potential_mapping = context['potential_mapping']
    options = context['options']
    base_kpoints_spacing = context['base_kpoints_spacing']
    clean_workdir = context['clean_workdir']
    input_structure = context['input_structure']

    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Stage configuration
    target_species = stage['target_species']
    potential_values = stage.get('potential_values', DEFAULT_POTENTIAL_VALUES)
    ldaul = stage.get('ldaul', 2)
    ldauj = stage.get('ldauj', 0.0)
    stage_kpoints_spacing = stage.get('kpoints_spacing', base_kpoints_spacing)

    # INCAR: the brick uses a single 'incar' key as base parameters
    # for both ground state and response calculations
    base_incar = stage.get('incar', {
        'encut': 520,
        'ediff': 1e-6,
        'ismear': 0,
        'sigma': 0.05,
        'prec': 'Accurate',
        'algo': 'Normal',
        'nelm': 100,
    })

    # Resolve input structure
    structure_from = stage['structure_from']
    if structure_from == 'input':
        stage_structure = input_structure
    else:
        from . import resolve_structure_from
        stage_structure = resolve_structure_from(structure_from, context)

    # Get species order for LDAU arrays
    # Note: we need the actual StructureData to get species order.
    # If stage_structure is a socket (output of previous task), we need
    # the species from the input structure as a proxy (same composition).
    if isinstance(stage_structure, orm.StructureData):
        all_species = get_species_order_from_structure(stage_structure)
    else:
        # Use input structure as proxy for species order
        all_species = get_species_order_from_structure(input_structure)

    lmaxmix = 4 if ldaul == 2 else 6

    # Parser settings that request orbital data
    settings = get_vasp_parser_settings(
        add_energy=True,
        add_trajectory=True,
        add_structure=True,
        add_kpoints=True,
    )

    # =========================================================================
    # STEP 1: Ground State Calculation
    # =========================================================================
    gs_incar = prepare_ground_state_incar(
        base_params=base_incar,
        lmaxmix=lmaxmix,
    )

    gs_builder_inputs = _prepare_builder_inputs(
        incar=gs_incar,
        kpoints_spacing=stage_kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        options=options,
        retrieve=['OUTCAR'],
        restart_folder=None,
        clean_workdir=False,  # MUST keep CHGCAR/WAVECAR for responses
    )
    # Add parser settings
    if 'settings' in gs_builder_inputs:
        existing = gs_builder_inputs['settings'].get_dict()
        existing.update(settings)
        gs_builder_inputs['settings'] = orm.Dict(dict=existing)
    else:
        gs_builder_inputs['settings'] = orm.Dict(dict=settings)

    ground_state = wg.add_task(
        VaspTask,
        name=f'gs_{stage_name}',
        structure=stage_structure,
        code=code,
        **gs_builder_inputs,
    )

    # Extract ground state d-electron occupation
    gs_occupation = wg.add_task(
        extract_d_electron_occupation,
        name=f'gs_occ_{stage_name}',
        retrieved=ground_state.outputs.retrieved,
        target_species=orm.Str(target_species),
        structure=stage_structure,
    )

    # =========================================================================
    # STEP 2-3: Response Calculations for Each Potential Value
    # =========================================================================
    response_tasks = {}

    for i, V in enumerate(potential_values):
        label = f'V_{i}'
        V_str = f'{V:+.2f}'.replace('.', 'p').replace('-', 'm').replace('+', 'p')

        # ----- Non-SCF Response (ICHARG=11) -----
        nscf_incar = prepare_response_incar(
            base_params=base_incar,
            potential_value=V,
            target_species=target_species,
            all_species=all_species,
            ldaul=ldaul,
            ldauj=ldauj,
            is_scf=False,
            lmaxmix=lmaxmix,
        )

        nscf_builder_inputs = _prepare_builder_inputs(
            incar=nscf_incar,
            kpoints_spacing=stage_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            retrieve=['OUTCAR'],
            restart_folder=None,
            clean_workdir=clean_workdir,
        )
        if 'settings' in nscf_builder_inputs:
            existing = nscf_builder_inputs['settings'].get_dict()
            existing.update(settings)
            nscf_builder_inputs['settings'] = orm.Dict(dict=existing)
        else:
            nscf_builder_inputs['settings'] = orm.Dict(dict=settings)

        nscf_task = wg.add_task(
            VaspTask,
            name=f'nscf_{V_str}_{stage_name}',
            structure=stage_structure,
            code=code,
            restart={'folder': ground_state.outputs.remote_folder},
            **nscf_builder_inputs,
        )

        nscf_occ = wg.add_task(
            extract_d_electron_occupation,
            name=f'nscf_occ_{V_str}_{stage_name}',
            retrieved=nscf_task.outputs.retrieved,
            target_species=orm.Str(target_species),
            structure=stage_structure,
        )

        # ----- SCF Response -----
        scf_incar = prepare_response_incar(
            base_params=base_incar,
            potential_value=V,
            target_species=target_species,
            all_species=all_species,
            ldaul=ldaul,
            ldauj=ldauj,
            is_scf=True,
            lmaxmix=lmaxmix,
        )

        scf_builder_inputs = _prepare_builder_inputs(
            incar=scf_incar,
            kpoints_spacing=stage_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
            retrieve=['OUTCAR'],
            restart_folder=None,
            clean_workdir=clean_workdir,
        )
        if 'settings' in scf_builder_inputs:
            existing = scf_builder_inputs['settings'].get_dict()
            existing.update(settings)
            scf_builder_inputs['settings'] = orm.Dict(dict=existing)
        else:
            scf_builder_inputs['settings'] = orm.Dict(dict=settings)

        scf_task = wg.add_task(
            VaspTask,
            name=f'scf_{V_str}_{stage_name}',
            structure=stage_structure,
            code=code,
            restart={'folder': ground_state.outputs.remote_folder},
            **scf_builder_inputs,
        )

        scf_occ = wg.add_task(
            extract_d_electron_occupation,
            name=f'scf_occ_{V_str}_{stage_name}',
            retrieved=scf_task.outputs.retrieved,
            target_species=orm.Str(target_species),
            structure=stage_structure,
        )

        # Calculate response for this potential
        response = wg.add_task(
            calculate_occupation_response,
            name=f'response_{V_str}_{stage_name}',
            ground_state_occupation=gs_occupation.outputs.result,
            nscf_occupation=nscf_occ.outputs.result,
            scf_occupation=scf_occ.outputs.result,
            potential_value=orm.Float(V),
        )

        response_tasks[label] = response

    # =========================================================================
    # STEP 4: Gather Responses and Calculate U
    # =========================================================================
    gather_kwargs = {
        label: resp_task.outputs.result
        for label, resp_task in response_tasks.items()
    }
    gathered = wg.add_task(
        gather_responses,
        name=f'gather_{stage_name}',
        **gather_kwargs,
    )

    calc_u = wg.add_task(
        calculate_hubbard_u_linear_regression,
        name=f'calc_u_{stage_name}',
        responses=gathered.outputs.result,
    )

    # =========================================================================
    # STEP 5: Compile Summary
    # =========================================================================
    summary = wg.add_task(
        compile_u_calculation_summary,
        name=f'summary_{stage_name}',
        hubbard_u_result=calc_u.outputs.result,
        ground_state_occupation=gs_occupation.outputs.result,
        structure=stage_structure,
        target_species=orm.Str(target_species),
        ldaul=orm.Int(ldaul),
    )

    return {
        'ground_state': ground_state,
        'gs_occupation': gs_occupation,
        'calc_u': calc_u,
        'summary': summary,
        'gathered': gathered,
        'structure': stage_structure,
    }


def expose_stage_outputs(wg, stage_name, stage_tasks_result):
    """Expose Hubbard U stage outputs on the WorkGraph.

    Args:
        wg: WorkGraph instance.
        stage_name: Unique stage identifier.
        stage_tasks_result: Dict returned by create_stage_tasks.
    """
    summary = stage_tasks_result['summary']
    calc_u = stage_tasks_result['calc_u']
    ground_state = stage_tasks_result['ground_state']
    gs_occupation = stage_tasks_result['gs_occupation']
    gathered = stage_tasks_result['gathered']

    setattr(wg.outputs, f'{stage_name}_summary',
            summary.outputs.result)
    setattr(wg.outputs, f'{stage_name}_hubbard_u_result',
            calc_u.outputs.result)
    setattr(wg.outputs, f'{stage_name}_ground_state_occupation',
            gs_occupation.outputs.result)
    setattr(wg.outputs, f'{stage_name}_ground_state_misc',
            ground_state.outputs.misc)
    setattr(wg.outputs, f'{stage_name}_ground_state_remote',
            ground_state.outputs.remote_folder)
    setattr(wg.outputs, f'{stage_name}_all_responses',
            gathered.outputs.result)


def get_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a Hubbard U stage in a sequential workflow.

    Args:
        wg_node: The WorkGraph ProcessNode.
        wg_pk: WorkGraph PK.
        stage_name: Name of the Hubbard U stage.

    Returns:
        Dict with keys: summary, hubbard_u_eV, target_species, chi_r2,
        chi_0_r2, response_data, pk, stage, type.
    """
    result = {
        'summary': None,
        'hubbard_u_eV': None,
        'target_species': None,
        'chi_r2': None,
        'chi_0_r2': None,
        'response_data': None,
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'hubbard_u',
    }

    # Try to access via WorkGraph outputs (exposed outputs)
    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        summary_attr = f'{stage_name}_summary'
        if hasattr(outputs, summary_attr):
            summary_node = getattr(outputs, summary_attr)
            if hasattr(summary_node, 'get_dict'):
                summary_dict = summary_node.get_dict()
                result['summary'] = summary_dict
                # Extract key values
                if 'summary' in summary_dict:
                    result['hubbard_u_eV'] = summary_dict['summary'].get(
                        'hubbard_u_eV')
                    result['target_species'] = summary_dict['summary'].get(
                        'target_species')
                if 'linear_fit' in summary_dict:
                    lf = summary_dict['linear_fit']
                    result['chi_r2'] = lf.get('chi_scf', {}).get('r_squared')
                    result['chi_0_r2'] = lf.get('chi_0_nscf', {}).get(
                        'r_squared')
                if 'response_data' in summary_dict:
                    result['response_data'] = summary_dict['response_data']

    # Fallback: traverse links
    if result['summary'] is None:
        _extract_hubbard_u_stage_from_workgraph(wg_node, stage_name, result)

    return result


def _extract_hubbard_u_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """Extract Hubbard U stage results by traversing WorkGraph links.

    Args:
        wg_node: The WorkGraph ProcessNode.
        stage_name: Name of the Hubbard U stage.
        result: Result dict to populate (modified in place).
    """
    if not hasattr(wg_node, 'base'):
        return

    summary_task_name = f'summary_{stage_name}'

    called_calc = wg_node.base.links.get_outgoing(
        link_type=LinkType.CALL_CALC)
    for link in called_calc.all():
        child_node = link.node
        link_label = link.link_label

        if summary_task_name in link_label or link_label == summary_task_name:
            created = child_node.base.links.get_outgoing(
                link_type=LinkType.CREATE)
            for out_link in created.all():
                out_label = out_link.link_label
                out_node = out_link.node

                if out_label == 'result' and hasattr(out_node, 'get_dict'):
                    summary_dict = out_node.get_dict()
                    result['summary'] = summary_dict
                    if 'summary' in summary_dict:
                        result['hubbard_u_eV'] = summary_dict['summary'].get(
                            'hubbard_u_eV')
                        result['target_species'] = summary_dict[
                            'summary'].get('target_species')
                    if 'linear_fit' in summary_dict:
                        lf = summary_dict['linear_fit']
                        result['chi_r2'] = lf.get('chi_scf', {}).get(
                            'r_squared')
                        result['chi_0_r2'] = lf.get('chi_0_nscf', {}).get(
                            'r_squared')
                    if 'response_data' in summary_dict:
                        result['response_data'] = summary_dict['response_data']
            break


def print_stage_results(
    index: int, stage_name: str, stage_result: dict
) -> None:
    """Print formatted results for a Hubbard U stage.

    Args:
        index: 1-based stage index.
        stage_name: Name of the stage.
        stage_result: Result dict from get_stage_results.
    """
    print(f"  [{index}] {stage_name} (HUBBARD U)")

    if stage_result['hubbard_u_eV'] is not None:
        print(f"      U = {stage_result['hubbard_u_eV']:.3f} eV")

    if stage_result['target_species'] is not None:
        print(f"      Target: {stage_result['target_species']}")

    if stage_result['chi_r2'] is not None:
        print(f"      SCF fit R\u00b2: {stage_result['chi_r2']:.6f}")

    if stage_result['chi_0_r2'] is not None:
        print(f"      NSCF fit R\u00b2: {stage_result['chi_0_r2']:.6f}")

    if stage_result['response_data'] is not None:
        rd = stage_result['response_data']
        potentials = rd.get('potential_values_eV', [])
        if potentials:
            print(f"      Potentials: {potentials}")

    if stage_result['summary'] is not None:
        summary = stage_result['summary']
        gs = summary.get('ground_state', {})
        avg_d = gs.get('average_d_per_atom')
        if avg_d is not None:
            species = stage_result['target_species'] or '?'
            print(f"      Avg d-occupation per {species}: {avg_d:.3f}")
