"""Helper tasks for custom VASP calculations."""

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def extract_total_energy(misc: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP misc output.

    Args:
        misc: VASP misc output Dict containing energy data

    Returns:
        Total energy as Float (eV)
    """
    misc_dict = misc.get_dict()

    # Navigate to total_energies if present
    energy_dict = misc_dict
    if 'total_energies' in misc_dict:
        energy_dict = misc_dict['total_energies']

    # Try multiple keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(float(energy_dict[key]))

    # If no recognized key found, raise error with available keys
    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in misc output. Available keys: {available}')


@task.calcfunction
def extract_relaxed_structure(structure: orm.StructureData) -> orm.StructureData:
    """
    Pass through relaxed structure from VASP output.

    In aiida-vasp, the relaxed structure is output directly as 'structure'
    (StructureData), not inside the 'misc' Dict.

    Args:
        structure: VASP output StructureData (relaxed structure)

    Returns:
        The same StructureData (for WorkGraph compatibility)
    """
    return structure


@task.calcfunction
def gather_energies(**kwargs) -> orm.Dict:
    """
    Gather energies from multiple calculations into a single Dict.

    Collects Float nodes from dynamic kwargs and returns a Dict with
    the energy values keyed by their input names.

    Args:
        **kwargs: Dynamic inputs where keys are structure identifiers
                  and values are orm.Float energy nodes

    Returns:
        orm.Dict with structure: {key: energy_value, ...}
    """
    energies = {}
    for key, energy_node in kwargs.items():
        if isinstance(energy_node, orm.Float):
            energies[key] = energy_node.value
        elif hasattr(energy_node, 'value'):
            energies[key] = float(energy_node.value)
        else:
            energies[key] = float(energy_node)
    return orm.Dict(dict=energies)


@task.calcfunction
def gather_misc(**kwargs) -> orm.Dict:
    """
    Gather misc outputs from multiple calculations into a single Dict.

    Collects Dict nodes from dynamic kwargs and returns a Dict with
    the misc data keyed by their input names.

    Args:
        **kwargs: Dynamic inputs where keys are structure identifiers
                  and values are orm.Dict misc nodes

    Returns:
        orm.Dict with structure: {key: misc_dict, ...}
    """
    misc_data = {}
    for key, misc_node in kwargs.items():
        if isinstance(misc_node, orm.Dict):
            misc_data[key] = misc_node.get_dict()
        elif isinstance(misc_node, dict):
            misc_data[key] = misc_node
        else:
            misc_data[key] = dict(misc_node)
    return orm.Dict(dict=misc_data)


@task.calcfunction
def gather_structure_pks(**kwargs) -> orm.Dict:
    """
    Gather structure PKs from multiple calculations into a single Dict.

    Collects StructureData nodes from dynamic kwargs and returns a Dict with
    their PKs keyed by the input names. This provides a lightweight reference
    to all structures without duplicating the data.

    Args:
        **kwargs: Dynamic inputs where keys are structure identifiers
                  and values are orm.StructureData nodes

    Returns:
        orm.Dict with structure: {key: pk, ...}
    """
    structure_pks = {}
    for key, struct_node in kwargs.items():
        if isinstance(struct_node, orm.StructureData):
            structure_pks[key] = struct_node.pk
        elif hasattr(struct_node, 'pk'):
            structure_pks[key] = struct_node.pk
        else:
            # Try to get pk from node
            structure_pks[key] = int(struct_node)
    return orm.Dict(dict=structure_pks)

