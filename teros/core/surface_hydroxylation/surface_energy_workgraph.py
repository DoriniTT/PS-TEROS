"""WorkGraph integration for surface energy calculations."""
import typing as t
from aiida.orm import Dict, Dict as AiidaDict
from aiida_workgraph import task, namespace


def _calculate_all_surface_energies_impl(
    structures_dict,
    energies_dict,
    bulk_structure,
    bulk_energy,
    temperature,
    pressures,
):
    """
    Implementation of surface energy calculations for all 3 reactions.

    This is the actual implementation function that can be called directly in tests.
    The @task-decorated version wraps this for WorkGraph execution.

    Args:
        structures_dict: Dict of {name: StructureData}
        energies_dict: Dict of {name: Float}
        bulk_structure: StructureData for bulk
        bulk_energy: Float for bulk energy
        temperature: Float for temperature in K
        pressures: Dict of partial pressures

    Returns:
        dict: Results for all 3 reactions with keys:
            - reaction1_results
            - reaction2_results
            - reaction3_results
    """
    from teros.core.surface_hydroxylation import calculate_surface_energies

    results = {}
    for reaction_num in [1, 2, 3]:
        try:
            results[f'reaction{reaction_num}_results'] = AiidaDict(dict=calculate_surface_energies(
                structures_dict=structures_dict,
                energies_dict=energies_dict,
                bulk_structure=bulk_structure,
                bulk_energy=bulk_energy,
                temperature=temperature,
                pressures=pressures,
                which_reaction=reaction_num,
            ))
        except Exception as e:
            # Log error but continue with other reactions
            results[f'reaction{reaction_num}_results'] = AiidaDict(dict={
                'error': str(e),
                'surface_energies': {},
                'formation_energies': {},
                'stoichiometry': {},
                'reference_data': {},
            })

    return results


@task()
def calculate_all_surface_energies(
    structures_dict,
    energies_dict,
    bulk_structure,
    bulk_energy,
    temperature,
    pressures,
) -> t.Annotated[dict, namespace(
    reaction1_results=Dict,
    reaction2_results=Dict,
    reaction3_results=Dict,
)]:
    """
    Calculate surface energies for all 3 reactions (WorkGraph task).

    Args:
        structures_dict: Dict of {name: StructureData}
        energies_dict: Dict of {name: Float}
        bulk_structure: StructureData for bulk
        bulk_energy: Float for bulk energy
        temperature: Float for temperature in K
        pressures: Dict of partial pressures

    Returns:
        dict: Results for all 3 reactions with keys:
            - reaction1_results
            - reaction2_results
            - reaction3_results
    """
    return _calculate_all_surface_energies_impl(
        structures_dict=structures_dict,
        energies_dict=energies_dict,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
        temperature=temperature,
        pressures=pressures,
    )


def create_surface_energy_task(
    aggregate_results_task,
    bulk_relax_task,
    pristine_relax_task,
    temperature=298.0,
    pressures=None,
):
    """
    Create WorkGraph task for surface energy calculations.

    Args:
        aggregate_results_task: WorkGraph task with structures/energies dict
        bulk_relax_task: WorkGraph task with bulk structure/energy
        pristine_relax_task: WorkGraph task with pristine structure/energy
        temperature: Temperature in Kelvin (default: 298.0)
        pressures: Dict of partial pressures in bar (default: standard)

    Returns:
        WorkGraph task with outputs reaction1_results, reaction2_results, reaction3_results
    """
    if pressures is None:
        pressures = {'H2O': 0.023, 'O2': 0.21, 'H2': 1.0}

    # Create WorkGraph task that depends on upstream tasks
    task_instance = calculate_all_surface_energies(
        structures_dict=aggregate_results_task.outputs.structures,
        energies_dict=aggregate_results_task.outputs.energies,
        bulk_structure=bulk_relax_task.outputs.structure,
        bulk_energy=bulk_relax_task.outputs.energy,
        temperature=temperature,
        pressures=pressures,
    )

    return task_instance
