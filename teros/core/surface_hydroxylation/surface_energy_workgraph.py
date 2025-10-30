"""WorkGraph integration for surface energy calculations."""
import typing as t
from aiida import orm
from aiida.orm import Dict, Dict as AiidaDict
from aiida_workgraph import task, namespace, dynamic


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
        structures_dict: Dict of {name: StructureData or Atoms}
        energies_dict: Dict of {name: Float}
        bulk_structure: StructureData or Atoms for bulk
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
    from aiida.orm import StructureData
    from ase import Atoms

    # Convert ASE Atoms to AiiDA StructureData if needed
    # (WorkGraph automatically converts StructureData → Atoms when passing between tasks)
    if isinstance(bulk_structure, Atoms):
        bulk_structure = StructureData(ase=bulk_structure)

    # Convert structures_dict entries from Atoms to StructureData if needed
    converted_structures = {}
    for name, structure in structures_dict.items():
        if isinstance(structure, Atoms):
            converted_structures[name] = StructureData(ase=structure)
        else:
            converted_structures[name] = structure

    results = {}
    for reaction_num in [1, 2, 3]:
        try:
            # Return unwrapped dict - @task.graph namespace annotation handles wrapping
            results[f'reaction{reaction_num}_results'] = calculate_surface_energies(
                structures_dict=converted_structures,
                energies_dict=energies_dict,
                bulk_structure=bulk_structure,
                bulk_energy=bulk_energy,
                temperature=temperature,
                pressures=pressures,
                which_reaction=reaction_num,
            )
        except Exception as e:
            # Log error but continue with other reactions
            results[f'reaction{reaction_num}_results'] = {
                'error': str(e),
                'surface_energies': {},
                'formation_energies': {},
                'stoichiometry': {},
                'reference_data': {},
            }

    return results


@task.calcfunction
def calculate_surface_energy_single_reaction(
    structure_pks: orm.Dict,
    energies_dict: orm.Dict,
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    temperature: orm.Float,
    pressures: orm.Dict,
    reaction_number: orm.Int,
) -> orm.Dict:
    """
    Calculate surface energy for a single reaction (calcfunction for provenance).

    Args:
        structure_pks: Dict node with {name: structure_pk}
        energies_dict: Dict node with {name: energy_value}
        bulk_structure: StructureData for bulk
        bulk_energy: Float for bulk energy
        temperature: Float for temperature in K
        pressures: Dict of partial pressures
        reaction_number: Int (1, 2, or 3)

    Returns:
        Dict node with surface energy results for the specified reaction
    """
    from teros.core.surface_hydroxylation import calculate_surface_energies

    # Load structures from PKs
    pks_dict = structure_pks.get_dict()
    structures_loaded = {}
    for name, pk in pks_dict.items():
        structures_loaded[name] = orm.load_node(pk)

    # Extract values from orm nodes
    energies_plain = energies_dict.get_dict()
    temperature_value = temperature.value
    pressures_plain = pressures.get_dict()
    bulk_energy_value = bulk_energy.value
    reaction_num = reaction_number.value

    try:
        result = calculate_surface_energies(
            structures_dict=structures_loaded,
            energies_dict=energies_plain,
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy_value,
            temperature=temperature_value,
            pressures=pressures_plain,
            which_reaction=reaction_num,
        )
        return orm.Dict(dict=result)
    except Exception as e:
        # Return error dict if calculation fails
        return orm.Dict(dict={
            'error': str(e),
            'surface_energies': {},
            'formation_energies': {},
            'stoichiometry': {},
            'reference_data': {},
        })


@task.graph
def calculate_all_surface_energies(
    structures_dict: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies_dict: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    temperature: orm.Float,
    pressures: orm.Dict,
) -> t.Annotated[dict, namespace(
    reaction1_results=orm.Dict,
    reaction2_results=orm.Dict,
    reaction3_results=orm.Dict,
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
    # Convert dynamic namespace structures to dict of PKs (JSON-serializable)
    structure_pks_dict = {}
    for key, structure in structures_dict.items():
        structure_pks_dict[key] = structure.pk

    # Convert energies to plain dict with values
    energies_as_dict = {}
    for key, energy in energies_dict.items():
        energies_as_dict[key] = energy.value

    # Wrap in orm.Dict for calcfunction
    structure_pks_node = orm.Dict(dict=structure_pks_dict)
    energies_dict_node = orm.Dict(dict=energies_as_dict)

    # Call calcfunction for each reaction and extract result with .result
    # Following the pattern from thermodynamics.py and slabs.py
    reaction1 = calculate_surface_energy_single_reaction(
        structure_pks=structure_pks_node,
        energies_dict=energies_dict_node,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
        temperature=temperature,
        pressures=pressures,
        reaction_number=orm.Int(1),
    ).result

    reaction2 = calculate_surface_energy_single_reaction(
        structure_pks=structure_pks_node,
        energies_dict=energies_dict_node,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
        temperature=temperature,
        pressures=pressures,
        reaction_number=orm.Int(2),
    ).result

    reaction3 = calculate_surface_energy_single_reaction(
        structure_pks=structure_pks_node,
        energies_dict=energies_dict_node,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
        temperature=temperature,
        pressures=pressures,
        reaction_number=orm.Int(3),
    ).result

    # Return plain dict containing orm.Dict nodes - namespace annotation handles wrapping
    return {
        'reaction1_results': reaction1,
        'reaction2_results': reaction2,
        'reaction3_results': reaction3,
    }


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
