"""Node builders for slab-related operations."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

VaspWorkChain = WorkflowFactory('vasp.vasp')


def build_scf_slabs_nodes(
    wg,
    slabs: dict[str, orm.StructureData],
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    parameters: dict,
    options: dict,
    clean_workdir: bool = False,
) -> dict[str, t.Any]:
    """
    Add VASP SCF (single-point) calculation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        slabs: Dictionary of {slab_id: StructureData}
        code: VASP code
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        parameters: VASP INCAR parameters
        options: Computer options
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary of {slab_id: vasp_task_node}
    """
    from .utils import prepare_vasp_parameters

    scf_nodes = {}

    for slab_id, slab_structure in slabs.items():
        # Prepare parameters
        vasp_params = prepare_vasp_parameters(
            base_parameters=parameters,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            kpoints_spacing=kpoints_spacing,
            options=options,
            clean_workdir=clean_workdir,
        )

        # Add VASP SCF node
        scf_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"scf_slab_{slab_id}",
            structure=slab_structure,
            **vasp_params,
        )

    return scf_nodes


def build_relax_slabs_nodes(
    wg,
    slabs: dict[str, orm.StructureData],
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    parameters: dict,
    options: dict,
    clean_workdir: bool = False,
) -> dict[str, t.Any]:
    """
    Add VASP relaxation nodes for each slab.

    Args:
        wg: WorkGraph instance to add nodes to
        slabs: Dictionary of {slab_id: StructureData}
        code: VASP code
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        parameters: VASP INCAR parameters
        options: Computer options
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary of {slab_id: vasp_task_node}
    """
    from .utils import prepare_vasp_parameters

    relax_nodes = {}

    for slab_id, slab_structure in slabs.items():
        # Prepare parameters
        vasp_params = prepare_vasp_parameters(
            base_parameters=parameters,
            code=code,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            kpoints_spacing=kpoints_spacing,
            options=options,
            clean_workdir=clean_workdir,
        )

        # Add VASP relax node
        relax_nodes[slab_id] = wg.add_task(
            VaspWorkChain,
            name=f"relax_slab_{slab_id}",
            structure=slab_structure,
            **vasp_params,
        )

    return relax_nodes


def build_energy_extraction_nodes(
    wg,
    vasp_nodes: dict[str, t.Any],
    node_type: str = "relaxed",
) -> dict[str, t.Any]:
    """
    Add energy extraction nodes for VASP calculations.

    Args:
        wg: WorkGraph instance to add nodes to
        vasp_nodes: Dictionary of {slab_id: vasp_task_node}
        node_type: Type of calculation ("relaxed" or "scf")

    Returns:
        Dictionary of {slab_id: energy_extraction_node}
    """
    # Import the existing extract_total_energy calcfunction
    from teros.core.slabs import extract_total_energy

    energy_nodes = {}

    for slab_id, vasp_node in vasp_nodes.items():
        energy_nodes[slab_id] = wg.add_task(
            extract_total_energy,
            name=f"extract_energy_{node_type}_{slab_id}",
            misc=vasp_node.outputs.misc,
        )

    return energy_nodes


def build_relaxation_energy_nodes(
    wg,
    unrelaxed_energies: dict[str, t.Any],
    relaxed_energies: dict[str, t.Any],
) -> dict[str, t.Any]:
    """
    Add relaxation energy calculation nodes.

    Args:
        wg: WorkGraph instance to add nodes to
        unrelaxed_energies: Dictionary of {slab_id: unrelaxed_energy_node}
        relaxed_energies: Dictionary of {slab_id: relaxed_energy_node}

    Returns:
        Dictionary of {slab_id: relaxation_energy_node}
    """
    relaxation_nodes = {}

    for slab_id in unrelaxed_energies.keys():
        relaxation_nodes[slab_id] = wg.add_task(
            calculate_relaxation_energy,
            name=f"relaxation_energy_{slab_id}",
            unrelaxed_energy=unrelaxed_energies[slab_id].outputs.result,
            relaxed_energy=relaxed_energies[slab_id].outputs.result,
        )

    return relaxation_nodes


@task.calcfunction
def calculate_relaxation_energy(
    unrelaxed_energy: orm.Float,
    relaxed_energy: orm.Float,
) -> orm.Float:
    """
    Calculate relaxation energy (unrelaxed - relaxed).

    Args:
        unrelaxed_energy: Energy before relaxation
        relaxed_energy: Energy after relaxation

    Returns:
        Relaxation energy (positive means structure lowered energy)
    """
    return orm.Float(unrelaxed_energy.value - relaxed_energy.value)
