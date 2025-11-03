"""Node builders for thermodynamics calculations."""

import typing as t
from aiida import orm
from aiida_workgraph import task


def build_surface_energy_nodes(
    wg,
    bulk_structure: orm.StructureData,
    bulk_energy: t.Any,
    slab_structures: dict[str, orm.StructureData],
    slab_energies: dict[str, t.Any],
    reference_energies: t.Any,
    formation_enthalpy: t.Any,
    oxide_type: t.Any,
    sampling: int = 100,
) -> dict[str, t.Any]:
    """
    Add surface energy calculation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        bulk_structure: Bulk structure
        bulk_energy: Bulk energy node output
        slab_structures: Dictionary of {slab_id: StructureData}
        slab_energies: Dictionary of {slab_id: energy_node}
        reference_energies: Reference energies Dict node
        formation_enthalpy: Formation enthalpy Dict node
        oxide_type: Oxide type Str node ('binary' or 'ternary')
        sampling: Number of sampling points for chemical potential grid

    Returns:
        Dictionary of {slab_id: surface_energy_node}
    """
    # Import existing calcfunctions
    from teros.core.thermodynamics import (
        calculate_surface_energy_binary,
        calculate_surface_energy_ternary,
    )

    surface_energy_nodes = {}

    for slab_id, slab_structure in slab_structures.items():
        # Add conditional node based on oxide type
        # For simplicity, we'll create both and the workflow will use the right one

        # Binary surface energy
        binary_node = wg.add_task(
            calculate_surface_energy_binary,
            name=f"surface_energy_binary_{slab_id}",
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energies[slab_id].outputs.result,
            reference_energies=reference_energies,
            formation_enthalpy=formation_enthalpy,
            sampling=orm.Int(sampling),
        )

        # Ternary surface energy
        ternary_node = wg.add_task(
            calculate_surface_energy_ternary,
            name=f"surface_energy_ternary_{slab_id}",
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energies[slab_id].outputs.result,
            reference_energies=reference_energies,
            formation_enthalpy=formation_enthalpy,
            sampling=orm.Int(sampling),
        )

        # Store both for now (we'll handle selection in main workgraph)
        surface_energy_nodes[slab_id] = {
            'binary': binary_node,
            'ternary': ternary_node,
        }

    return surface_energy_nodes


@task.calcfunction
def select_surface_energy_by_oxide_type(
    oxide_type: orm.Str,
    binary_result: orm.Dict,
    ternary_result: orm.Dict,
) -> orm.Dict:
    """
    Select the appropriate surface energy result based on oxide type.

    Args:
        oxide_type: 'binary' or 'ternary'
        binary_result: Result from binary calculation
        ternary_result: Result from ternary calculation

    Returns:
        The appropriate result based on oxide type
    """
    if oxide_type.value == 'binary':
        return binary_result
    elif oxide_type.value == 'ternary':
        return ternary_result
    else:
        raise ValueError(f"Unknown oxide type: {oxide_type.value}")
