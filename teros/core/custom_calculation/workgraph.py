"""WorkGraph builder for custom VASP calculations."""

import re
import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task, namespace, dynamic

from .tasks import extract_total_energy
from ..fixed_atoms import get_fixed_atoms_list


def _sanitize_key(label: str, index: int) -> str:
    """
    Sanitize a label to be a valid Python identifier for AiiDA link labels.

    Args:
        label: The original label string
        index: Fallback index if label is empty or invalid

    Returns:
        A valid Python identifier string
    """
    if not label:
        return f'structure_{index}'

    # Replace common problematic characters
    key = label.replace('.', '_').replace(' ', '_').replace('-', '_')

    # Remove any remaining non-alphanumeric characters (except underscore)
    key = re.sub(r'[^a-zA-Z0-9_]', '', key)

    # Ensure it doesn't start with a digit (prepend 's_' if it does)
    if key and key[0].isdigit():
        key = f's_{key}'

    # If empty after sanitization, use fallback
    if not key:
        return f'structure_{index}'

    return key


@task.graph
def CustomCalculationMultipleWorkGraph(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.InstalledCode,
    builder_inputs_list: list,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: list = None,
) -> t.Annotated[dict, namespace(**{
    'structures': dynamic(orm.StructureData),
    'energies': dynamic(orm.Float),
    'misc': dynamic(orm.Dict),
})]:
    """
    WorkGraph task for multiple custom VASP calculations.

    Args:
        structures: Dict of StructureData keyed by sanitized label
        code: VASP InstalledCode
        builder_inputs_list: List of builder input dicts (same order as structures)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix

    Returns:
        Dictionary with namespace outputs:
            - structures: Dict {key: StructureData} for relaxed structures
            - energies: Dict {key: Float} for energies
            - misc: Dict {key: Dict} for misc outputs
    """
    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    structures_out = {}
    energies_out = {}
    misc_out = {}

    # Process each structure
    keys = list(structures.keys())
    for i, key in enumerate(keys):
        struct = structures[key]
        inputs = builder_inputs_list[i]

        # Prepare builder inputs
        prepared_inputs = _prepare_builder_inputs(
            inputs, struct, fix_type, fix_thickness, fix_elements
        )

        # Add VaspTask
        vasp_task = VaspTask(
            structure=struct,
            code=code,
            **prepared_inputs
        )

        # Extract energy
        energy_task = extract_total_energy(misc=vasp_task.misc)

        # Collect outputs
        structures_out[key] = vasp_task.structure
        energies_out[key] = energy_task.result
        misc_out[key] = vasp_task.misc

    return {
        'structures': structures_out,
        'energies': energies_out,
        'misc': misc_out,
    }


def build_custom_calculation_workgraph(
    structure,
    code_label,
    builder_inputs,
    name='custom_calc',
    max_concurrent_jobs=None,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
):
    """
    Build a WorkGraph for custom VASP calculations.

    Args:
        structure: StructureData or list of StructureData
        code_label: str, VASP code label (e.g., 'VASP-6.4.1@cluster02')
        builder_inputs: dict or list of dicts with VASP builder parameters:
            - parameters: Dict with nested 'incar' dict
            - options: Dict with resources, queue, etc.
            - kpoints_spacing: float
            - potential_family: str
            - potential_mapping: dict
            - clean_workdir: bool
            - settings: Dict (optional)
        name: str, WorkGraph name
        max_concurrent_jobs: int, optional - Maximum number of concurrent VASP calculations
                           (only applies to multiple structure workflows). Default: None (unlimited)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None). Default: None (no fixing)
        fix_thickness: Thickness in Angstroms for fixing region. Default: 0.0
        fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O']).
                     If None, all elements in the region are fixed. Default: None

    Returns:
        WorkGraph ready to submit
    """
    # Load code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Detect single vs multiple structures
    is_single = isinstance(structure, orm.StructureData)

    if is_single:
        # Single structure workflow - use plain WorkGraph
        wg = WorkGraph(name=name)

        # Prepare builder inputs (convert plain dicts to AiiDA types)
        prepared_inputs = _prepare_builder_inputs(
            builder_inputs, structure, fix_type, fix_thickness, fix_elements
        )

        # Add VaspTask
        vasp_task = wg.add_task(
            VaspTask,
            name='vasp_calc',
            structure=structure,
            code=code,
            **prepared_inputs
        )

        # Extract energy
        energy_task = wg.add_task(
            extract_total_energy,
            name='extract_energy',
            misc=vasp_task.outputs.misc
        )

        # Set WorkGraph outputs (use 'any' socket type identifier)
        wg.add_output('any', 'energy', from_socket=energy_task.outputs.result)
        wg.add_output('any', 'structure', from_socket=vasp_task.outputs.structure)
        wg.add_output('any', 'misc', from_socket=vasp_task.outputs.misc)

        return wg
    else:
        # Multiple structures workflow - use @task.graph pattern
        structures = structure  # Rename for clarity

        # Handle builder_inputs: single dict or list of dicts
        if isinstance(builder_inputs, dict):
            # Same builder inputs for all structures
            builder_inputs_list = [builder_inputs] * len(structures)
        elif isinstance(builder_inputs, list):
            # Different builder inputs for each structure
            if len(builder_inputs) != len(structures):
                raise ValueError(
                    f"builder_inputs list length ({len(builder_inputs)}) must match "
                    f"structures list length ({len(structures)})"
                )
            builder_inputs_list = builder_inputs
        else:
            raise TypeError("builder_inputs must be dict or list of dicts")

        # Prepare structures dict with sanitized keys
        structures_dict = {}
        for i, struct in enumerate(structures):
            label = struct.label if hasattr(struct, 'label') else ''
            key = _sanitize_key(label, i)
            structures_dict[key] = struct

        # Build WorkGraph using @task.graph function
        wg = CustomCalculationMultipleWorkGraph.build(
            structures=structures_dict,
            code=code,
            builder_inputs_list=builder_inputs_list,
            fix_type=fix_type,
            fix_thickness=fix_thickness,
            fix_elements=fix_elements,
        )

        # Set name and max concurrent jobs
        wg.name = name
        if max_concurrent_jobs is not None:
            wg.max_number_jobs = max_concurrent_jobs

        return wg


def _prepare_builder_inputs(
    builder_inputs,
    structure: orm.StructureData = None,
    fix_type: str = None,
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
):
    """
    Prepare builder inputs by converting plain dicts to AiiDA types.

    Also applies selective dynamics if fix_type is specified.

    Args:
        builder_inputs: dict with VASP builder parameters
        structure: StructureData to apply selective dynamics to (optional)
        fix_type: Where to fix atoms ('bottom'/'top'/'center'/None)
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to fix

    Returns:
        dict with AiiDA-compatible types
    """
    prepared = {}

    # Convert parameters dict to orm.Dict
    if 'parameters' in builder_inputs:
        if isinstance(builder_inputs['parameters'], dict):
            prepared['parameters'] = orm.Dict(dict=builder_inputs['parameters'])
        else:
            prepared['parameters'] = builder_inputs['parameters']

    # Convert options dict to orm.Dict
    if 'options' in builder_inputs:
        if isinstance(builder_inputs['options'], dict):
            prepared['options'] = orm.Dict(dict=builder_inputs['options'])
        else:
            prepared['options'] = builder_inputs['options']

    # Convert potential_mapping dict to orm.Dict
    if 'potential_mapping' in builder_inputs:
        if isinstance(builder_inputs['potential_mapping'], dict):
            prepared['potential_mapping'] = orm.Dict(dict=builder_inputs['potential_mapping'])
        else:
            prepared['potential_mapping'] = builder_inputs['potential_mapping']

    # Convert settings dict to orm.Dict if present
    if 'settings' in builder_inputs:
        if isinstance(builder_inputs['settings'], dict):
            prepared['settings'] = orm.Dict(dict=builder_inputs['settings'])
        else:
            prepared['settings'] = builder_inputs['settings']

    # Handle kpoints_spacing - ensure it's a float
    if 'kpoints_spacing' in builder_inputs:
        kps = builder_inputs['kpoints_spacing']
        # Keep as plain Python float (don't wrap in orm.Float)
        if isinstance(kps, (int, float)):
            prepared['kpoints_spacing'] = float(kps)
        else:
            prepared['kpoints_spacing'] = kps

    # Copy string/bool values directly
    for key in ('potential_family', 'clean_workdir'):
        if key in builder_inputs:
            prepared[key] = builder_inputs[key]

    # Apply selective dynamics if fix_type is specified
    if structure is not None and fix_type is not None and fix_thickness > 0.0:
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

            # Add selective dynamics to VASP inputs
            prepared['dynamics'] = orm.Dict(dict={
                'positions_dof': positions_dof
            })

    return prepared


def get_custom_results(workgraph):
    """
    Extract results from completed custom calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_custom_calculation_workgraph

    Returns:
        dict with:
            - energies: dict (keyed by structure label) or single float
            - structures: dict of StructureData (keyed by structure label), or single StructureData
            - misc: dict (keyed by structure label) or single dict

    For multiple structures, results are organized by structure key (derived from
    structure labels, or 'structure_0', 'structure_1', etc. if no label):
        {
            'energies': {'my_struct_1': -123.45, 'my_struct_2': -234.56, ...},
            'structures': {'my_struct_1': <StructureData>, 'my_struct_2': <StructureData>, ...},
            'misc': {'my_struct_1': {...}, 'my_struct_2': {...}, ...}
        }

    Example:
        >>> results = get_custom_results(wg)
        >>> for key, energy in results['energies'].items():
        ...     print(f"{key}: {energy:.3f} eV")
        >>> # Access structures directly
        >>> for key, struct in results['structures'].items():
        ...     print(f"{key}: {struct.get_formula()}")
    """
    results = {}

    # Check if single or multiple structure workflow
    if hasattr(workgraph.outputs, 'energy'):
        # Single structure
        results['energies'] = workgraph.outputs.energy.value
        results['structures'] = workgraph.outputs.structure
        results['misc'] = workgraph.outputs.misc.get_dict()
    elif hasattr(workgraph.outputs, 'energies'):
        # Multiple structures with namespaced outputs
        # Access energies namespace
        energies = {}
        for key in workgraph.outputs.energies:
            energies[key] = workgraph.outputs.energies[key].value
        results['energies'] = energies

        # Access structures namespace
        structures = {}
        for key in workgraph.outputs.structures:
            structures[key] = workgraph.outputs.structures[key]
        results['structures'] = structures

        # Access misc namespace
        misc_dict = {}
        for key in workgraph.outputs.misc:
            misc_dict[key] = workgraph.outputs.misc[key].get_dict()
        results['misc'] = misc_dict
    else:
        # Legacy: Multiple structures - access via task outputs
        i = 0
        energies = {}
        structures = {}
        misc_dict = {}

        # Find all energy extraction tasks
        while f'extract_energy_{i}' in workgraph.tasks:
            energy_task = workgraph.tasks[f'extract_energy_{i}']
            vasp_task = workgraph.tasks[f'vasp_calc_{i}']

            # Get structure to determine key (use same sanitization as workgraph builder)
            struct = vasp_task.outputs['structure'].value
            label = struct.label if hasattr(struct, 'label') else ''
            key = _sanitize_key(label, i)

            # Access outputs from completed tasks
            energies[key] = energy_task.outputs['result'].value
            structures[key] = struct
            misc_dict[key] = vasp_task.outputs['misc'].value.get_dict()
            i += 1

        if not energies:
            raise ValueError("WorkGraph does not have expected outputs (energy or extract_energy_0 tasks)")

        results['energies'] = energies
        results['structures'] = structures
        results['misc'] = misc_dict

    return results


# =============================================================================
# DOS CALCULATION FUNCTIONS
# =============================================================================

def build_dos_calculation_workgraph(
    structure: t.Union[orm.StructureData, t.List[orm.StructureData]],
    code_label: str,
    scf_inputs: dict,
    dos_inputs: dict,
    dos_kpoints_distance: float = 0.2,
    potential_family: str = 'PBE',
    potential_mapping: dict = None,
    options: dict = None,
    name: str = 'dos_calc',
    max_concurrent_jobs: int = None,
) -> WorkGraph:
    """
    Build a WorkGraph for DOS calculations using vasp.v2.bands workflow.

    This workflow:
    1. Performs SCF calculation to get charge density (CHGCAR)
    2. Performs non-SCF DOS calculation using the CHGCAR

    Args:
        structure: StructureData or list of StructureData
        code_label: str, VASP code label (e.g., 'VASP-6.4.1@cluster')
        scf_inputs: dict with SCF calculation parameters:
            - parameters: Dict with nested 'incar' dict for SCF
            - kpoints_spacing: float, k-point spacing for SCF
            - settings: Dict (optional) for additional retrieve list, etc.
        dos_inputs: dict with DOS calculation parameters:
            - parameters: Dict with nested 'incar' dict for DOS
            - settings: Dict (optional) for additional retrieve list, etc.
        dos_kpoints_distance: float, k-point spacing for DOS (default 0.2)
        potential_family: str, POTCAR family (default 'PBE')
        potential_mapping: dict, element to POTCAR mapping
        options: dict, calculation options (resources, queue, etc.)
        name: str, WorkGraph name
        max_concurrent_jobs: int, maximum concurrent calculations

    Returns:
        WorkGraph ready to submit

    Example:
        >>> scf_inputs = {
        ...     'parameters': {'incar': {'ENCUT': 500, 'EDIFF': 1e-5, ...}},
        ...     'kpoints_spacing': 0.3,
        ... }
        >>> dos_inputs = {
        ...     'parameters': {'incar': {'ICHARG': 11, 'ISMEAR': -5, 'NEDOS': 2001, ...}},
        ... }
        >>> wg = build_dos_calculation_workgraph(
        ...     structure=my_structure,
        ...     code_label='VASP-6.4.1@cluster',
        ...     scf_inputs=scf_inputs,
        ...     dos_inputs=dos_inputs,
        ... )
        >>> wg.submit()
    """
    # Load code
    code = orm.load_code(code_label)

    # Get bands workchain (used for DOS-only calculations)
    BandsWorkChain = WorkflowFactory('vasp.v2.bands')
    BandsTask = task(BandsWorkChain)

    # Detect single vs multiple structures
    is_single = isinstance(structure, orm.StructureData)

    if is_single:
        structures = [structure]
    else:
        structures = structure

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # Process each structure
    for i, struct in enumerate(structures):
        label = struct.label if hasattr(struct, 'label') and struct.label else ''
        key = _sanitize_key(label, i)

        # Prepare inputs for the bands workflow (DOS-only mode)
        prepared_inputs = _prepare_dos_builder_inputs(
            structure=struct,
            code=code,
            scf_inputs=scf_inputs,
            dos_inputs=dos_inputs,
            dos_kpoints_distance=dos_kpoints_distance,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            options=options,
        )

        # Add BandsTask (in DOS-only mode)
        task_name = f'dos_calc_{key}' if not is_single else 'dos_calc'
        wg.add_task(
            BandsTask,
            name=task_name,
            **prepared_inputs
        )

    # Set max concurrent jobs if specified
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg


def _prepare_dos_builder_inputs(
    structure: orm.StructureData,
    code: orm.InstalledCode,
    scf_inputs: dict,
    dos_inputs: dict,
    dos_kpoints_distance: float,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
) -> dict:
    """
    Prepare inputs for vasp.v2.bands workflow in DOS-only mode.

    Args:
        structure: Input structure
        code: VASP code
        scf_inputs: SCF calculation parameters
        dos_inputs: DOS calculation parameters
        dos_kpoints_distance: k-point spacing for DOS
        potential_family: POTCAR family
        potential_mapping: Element to POTCAR mapping
        options: Calculation options

    Returns:
        dict with prepared inputs for the bands workflow
    """
    inputs = {}

    # Structure
    inputs['structure'] = structure

    # Band settings for DOS-only mode
    band_settings = {
        'only_dos': True,
        'run_dos': True,
        'dos_kpoints_distance': dos_kpoints_distance,
    }
    inputs['band_settings'] = orm.Dict(dict=band_settings)

    # SCF namespace inputs
    scf_namespace = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict=potential_mapping or {}),
    }

    # SCF parameters
    if 'parameters' in scf_inputs:
        if isinstance(scf_inputs['parameters'], dict):
            scf_namespace['parameters'] = orm.Dict(dict=scf_inputs['parameters'])
        else:
            scf_namespace['parameters'] = scf_inputs['parameters']

    # SCF kpoints spacing
    if 'kpoints_spacing' in scf_inputs:
        scf_namespace['kpoints_spacing'] = float(scf_inputs['kpoints_spacing'])

    # SCF settings (for additional retrieve list, etc.)
    if 'settings' in scf_inputs:
        if isinstance(scf_inputs['settings'], dict):
            scf_namespace['settings'] = orm.Dict(dict=scf_inputs['settings'])
        else:
            scf_namespace['settings'] = scf_inputs['settings']

    # Options for SCF
    if options:
        if isinstance(options, dict):
            scf_namespace['options'] = orm.Dict(dict=options)
        else:
            scf_namespace['options'] = options

    inputs['scf'] = scf_namespace

    # DOS namespace inputs
    dos_namespace = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict=potential_mapping or {}),
    }

    # DOS parameters
    if 'parameters' in dos_inputs:
        if isinstance(dos_inputs['parameters'], dict):
            dos_namespace['parameters'] = orm.Dict(dict=dos_inputs['parameters'])
        else:
            dos_namespace['parameters'] = dos_inputs['parameters']

    # DOS settings (for additional retrieve list, etc.)
    if 'settings' in dos_inputs:
        if isinstance(dos_inputs['settings'], dict):
            dos_namespace['settings'] = orm.Dict(dict=dos_inputs['settings'])
        else:
            dos_namespace['settings'] = dos_inputs['settings']

    # Options for DOS
    if options:
        if isinstance(options, dict):
            dos_namespace['options'] = orm.Dict(dict=options)
        else:
            dos_namespace['options'] = options

    inputs['dos'] = dos_namespace

    return inputs


def get_dos_results(workgraph) -> dict:
    """
    Extract results from completed DOS calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_dos_calculation_workgraph

    Returns:
        dict with:
            - dos: dict of DOS data nodes (keyed by structure label)
            - projectors: dict of projector data nodes for projected DOS (keyed by structure label)
            - primitive_structure: dict of primitive structures used

    Example:
        >>> results = get_dos_results(wg)
        >>> for key, dos_node in results['dos'].items():
        ...     print(f"{key}: DOS data available")
        >>> # Access projected DOS data
        >>> for key, proj_node in results['projectors'].items():
        ...     print(f"{key}: Projected DOS data available")
    """
    results = {
        'dos': {},
        'projectors': {},
        'primitive_structure': {},
    }

    # Find all DOS calculation tasks
    for task_name, task_obj in workgraph.tasks.items():
        if task_name.startswith('dos_calc'):
            # Extract key from task name
            if task_name == 'dos_calc':
                key = 'structure_0'
            else:
                key = task_name.replace('dos_calc_', '')

            # Get outputs - check if task has completed and has outputs
            try:
                if hasattr(task_obj, 'outputs'):
                    outputs = task_obj.outputs

                    # Get DOS output
                    if hasattr(outputs, 'dos'):
                        dos_output = outputs.dos
                        if hasattr(dos_output, 'value') and dos_output.value is not None:
                            results['dos'][key] = dos_output.value

                    # Get projectors output (for projected DOS)
                    if hasattr(outputs, 'projectors'):
                        proj_output = outputs.projectors
                        if hasattr(proj_output, 'value') and proj_output.value is not None:
                            results['projectors'][key] = proj_output.value

                    # Get primitive structure output
                    if hasattr(outputs, 'primitive_structure'):
                        prim_output = outputs.primitive_structure
                        if hasattr(prim_output, 'value') and prim_output.value is not None:
                            results['primitive_structure'][key] = prim_output.value

            except Exception:
                # Task may not have completed or outputs not available
                pass

    return results
