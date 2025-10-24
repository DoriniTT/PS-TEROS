"""Child WorkGraph for parallel VASP relaxations with batch limiting."""

import typing as t
from aiida import orm
from aiida.engine import calcfunction
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic

from ..fixed_atoms import get_fixed_atoms_list


def get_settings():
    """Parser settings for aiida-vasp."""
    return {
        'parser_settings': {
            'add_trajectory': True,
            'add_structure': True,
            'add_kpoints': True,
        }
    }


def deep_merge_dicts(base: dict, override: dict) -> dict:
    """
    Deep merge override dict into base dict.

    For nested dicts, recursively merges. For other values, override replaces base.

    Args:
        base: Base dictionary
        override: Override dictionary (values take precedence)

    Returns:
        Merged dictionary

    Example:
        >>> base = {'a': 1, 'b': {'c': 2, 'd': 3}}
        >>> override = {'b': {'c': 99}, 'e': 5}
        >>> deep_merge_dicts(base, override)
        {'a': 1, 'b': {'c': 99, 'd': 3}, 'e': 5}
    """
    import copy
    result = copy.deepcopy(base)

    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            # Recursively merge nested dicts
            result[key] = deep_merge_dicts(result[key], value)
        else:
            # Override value
            result[key] = copy.deepcopy(value)

    return result


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
    builder_inputs: dict,
    max_parallel: int,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    structure_specific_builder_inputs: dict = None,
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
        builder_inputs: Complete builder configuration dict for vasp.v2.vasp WorkChain.
            Must contain all builder parameters (see workgraph.py for details).
            The 'code' and 'structure' keys will be set automatically for each structure.
        max_parallel: Maximum number of structures to process (limits to first N structures)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None). Default: None (no fixing)
        fix_thickness: Thickness in Angstroms for fixing region. Default: 0.0
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O']).
                     If None, all elements in the region are fixed
        structure_specific_builder_inputs: Optional dict to override builder_inputs for specific
                     structure indices. Keys are integer indices (0, 1, 2, ...) matching structure
                     order. Values are PARTIAL builder_inputs dicts that will be MERGED into the
                     default builder_inputs. This allows you to override only specific parameters
                     (e.g., change 'Algo' or increase 'ENCUT') while keeping all other parameters
                     from the default builder.
                     Example: {0: {'parameters': {'incar': {'ALGO': 'Normal'}}, 'kpoints_spacing': 0.2}}
                              only changes ALGO and kpoints_spacing for structure 0, keeping all
                              other parameters (potential_mapping, options, etc.) from default builder.
                     Structures not listed (e.g., 1, 3) will use the default builder_inputs.
                     Default: None (use default builder_inputs for all)

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
        idx_str = str(idx)

        # Check if there's a structure-specific builder for this index
        # Note: structure_specific_builder_inputs uses string keys (required by AiiDA)
        import copy
        if structure_specific_builder_inputs is not None and idx_str in structure_specific_builder_inputs:
            # Deep merge structure-specific overrides into default builder
            # This allows partial overrides (e.g., only changing 'Algo' or 'kpoints_spacing')
            # while keeping all other required parameters from the default builder
            vasp_inputs = deep_merge_dicts(builder_inputs, structure_specific_builder_inputs[idx_str])
        else:
            # Use default builder for all other structures
            vasp_inputs = copy.deepcopy(builder_inputs)

        # Override structure and code (these are set by the workflow)
        vasp_inputs['structure'] = structure
        vasp_inputs['code'] = code

        # Convert plain dicts to orm.Dict where needed
        if 'potential_mapping' in vasp_inputs and not isinstance(vasp_inputs['potential_mapping'], orm.Dict):
            vasp_inputs['potential_mapping'] = orm.Dict(dict=vasp_inputs['potential_mapping'])

        if 'settings' in vasp_inputs and not isinstance(vasp_inputs['settings'], orm.Dict):
            vasp_inputs['settings'] = orm.Dict(dict=vasp_inputs['settings'])
        elif 'settings' not in vasp_inputs:
            # Ensure settings are present (add defaults if not provided by user)
            vasp_inputs['settings'] = orm.Dict(dict=get_settings())

        # Apply selective dynamics if requested
        if fix_type is not None and fix_thickness > 0.0:
            # Get list of atoms to fix (1-based indices)
            fixed_atoms_list = get_fixed_atoms_list(
                structure=structure,
                fix_type=fix_type,
                fix_thickness=fix_thickness,
                fix_elements=fix_elements,
            )

            if fixed_atoms_list:
                # Create positions_dof array for selective dynamics
                # True = relax, False = fix
                num_atoms = len(structure.sites)
                positions_dof = []

                for i in range(1, num_atoms + 1):  # 1-based indexing
                    if i in fixed_atoms_list:
                        positions_dof.append([False, False, False])  # Fix atom
                    else:
                        positions_dof.append([True, True, True])  # Relax atom

                # Add selective dynamics to VASP inputs (plain dict, not orm.Dict)
                vasp_inputs['dynamics'] = {
                    'positions_dof': positions_dof
                }

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
