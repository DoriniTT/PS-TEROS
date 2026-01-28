"""Result extraction functions for the explorer module."""

import typing as t
from aiida import orm
from aiida.common.links import LinkType


def get_results(pk: int) -> dict:
    """
    Extract results from a completed explorer calculation.

    Args:
        pk: PK of the calculation (WorkGraph or VaspWorkChain)

    Returns:
        dict with:
            - energy: float (eV)
            - structure: StructureData (relaxed, if NSW > 0)
            - misc: dict (parsed VASP outputs)
            - files: FolderData (retrieved files)
            - pk: int (the original PK)
    """
    node = orm.load_node(pk)

    result = {
        'energy': None,
        'structure': None,
        'misc': None,
        'files': None,
        'pk': pk,
    }

    # Try to get outputs directly
    if hasattr(node, 'outputs'):
        outputs = node.outputs

        # Energy (might be exposed as workgraph output or in misc)
        if hasattr(outputs, 'energy'):
            energy_node = outputs.energy
            if hasattr(energy_node, 'value'):
                result['energy'] = energy_node.value
            else:
                result['energy'] = float(energy_node)

        # Structure
        if hasattr(outputs, 'structure'):
            result['structure'] = outputs.structure

        # Misc
        if hasattr(outputs, 'misc'):
            misc_node = outputs.misc
            if hasattr(misc_node, 'get_dict'):
                result['misc'] = misc_node.get_dict()
            else:
                result['misc'] = dict(misc_node)

        # Retrieved files
        if hasattr(outputs, 'retrieved'):
            result['files'] = outputs.retrieved
        elif hasattr(outputs, 'files'):
            result['files'] = outputs.files

    # For WorkGraph nodes, traverse to find VASP outputs
    if result['energy'] is None or result['misc'] is None:
        _extract_from_workgraph(node, result)

    # Extract energy from misc if not found directly
    if result['energy'] is None and result['misc'] is not None:
        result['energy'] = _extract_energy_from_misc(result['misc'])

    return result


def _extract_from_workgraph(node, result: dict) -> None:
    """
    Extract results by traversing WorkGraph links to find VASP outputs.

    Args:
        node: The WorkGraph node
        result: Result dict to populate (modified in place)
    """
    if not hasattr(node, 'base'):
        return

    # Try both CALL_WORK (for WorkChains) and CALL_CALC (for calcfunctions)
    for link_type in [LinkType.CALL_WORK, LinkType.CALL_CALC]:
        called = node.base.links.get_outgoing(link_type=link_type)
        for link in called.all():
            child_node = link.node

            # Check if this is a VASP calc
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get misc
                if result['misc'] is None and hasattr(outputs, 'misc'):
                    misc_node = outputs.misc
                    if hasattr(misc_node, 'get_dict'):
                        result['misc'] = misc_node.get_dict()

                # Get structure
                if result['structure'] is None and hasattr(outputs, 'structure'):
                    result['structure'] = outputs.structure

                # Get retrieved
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved

            # Recurse into child WorkGraphs
            if hasattr(child_node, 'base'):
                _extract_from_workgraph(child_node, result)


def _extract_energy_from_misc(misc: dict) -> t.Optional[float]:
    """
    Extract energy from misc dict.

    Args:
        misc: VASP misc output dict

    Returns:
        Energy in eV, or None if not found
    """
    # Navigate to total_energies if present
    energy_dict = misc
    if 'total_energies' in misc:
        energy_dict = misc['total_energies']

    # Try multiple keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return float(energy_dict[key])

    return None


def get_energy(pk: int) -> float:
    """
    Quick shortcut to get just the energy from a calculation.

    Args:
        pk: PK of the calculation

    Returns:
        Energy in eV

    Raises:
        ValueError: If energy cannot be extracted
    """
    results = get_results(pk)
    if results['energy'] is None:
        raise ValueError(f"Could not extract energy from calculation PK {pk}")
    return results['energy']


def get_batch_results(pks: t.Dict[str, int]) -> t.Dict[str, dict]:
    """
    Extract results from multiple calculations.

    Args:
        pks: Dict mapping keys to PKs (e.g., {'clean': 123, 'defect': 456})

    Returns:
        Dict mapping keys to result dicts (from get_results)
    """
    return {key: get_results(pk) for key, pk in pks.items()}


def get_batch_energies(pks: t.Dict[str, int]) -> t.Dict[str, float]:
    """
    Quick shortcut to get energies from multiple calculations.

    Args:
        pks: Dict mapping keys to PKs

    Returns:
        Dict mapping keys to energies (eV)
    """
    return {key: get_energy(pk) for key, pk in pks.items()}


def print_results(pk: int) -> None:
    """
    Print a formatted summary of calculation results.

    Args:
        pk: PK of the calculation
    """
    from .utils import get_status

    status = get_status(pk)
    print(f"Calculation PK {pk}")
    print(f"  Status: {status}")

    if status != 'finished':
        print("  (Calculation not finished, results may be incomplete)")
        return

    results = get_results(pk)

    if results['energy'] is not None:
        print(f"  Energy: {results['energy']:.6f} eV")

    if results['structure'] is not None:
        struct = results['structure']
        formula = struct.get_formula()
        print(f"  Structure: {formula} (PK: {struct.pk})")

    if results['files'] is not None:
        files = results['files'].list_object_names()
        print(f"  Retrieved files: {', '.join(files)}")

    if results['misc'] is not None:
        misc = results['misc']
        # Print some useful info from misc
        if 'run_status' in misc:
            print(f"  Run status: {misc['run_status']}")
        if 'maximum_force' in misc:
            print(f"  Max force: {misc['maximum_force']:.4f} eV/A")


def get_dos_results(pk: int) -> dict:
    """
    Extract results from a completed quick_dos calculation (BandsWorkChain).

    Args:
        pk: PK of the BandsWorkChain (from quick_dos)

    Returns:
        dict with:
            - energy: float (SCF energy in eV)
            - structure: StructureData (input structure)
            - scf_misc: dict (SCF VASP outputs)
            - scf_remote: RemoteData (SCF remote folder)
            - scf_retrieved: FolderData (SCF retrieved files)
            - dos_misc: dict (DOS VASP outputs)
            - dos_remote: RemoteData (DOS remote folder)
            - dos: ArrayData (DOS data, if available)
            - projectors: ArrayData (projected DOS, if available)
            - files: FolderData (retrieved files from DOS calculation)
            - pk: int (the original PK)
    """
    node = orm.load_node(pk)

    result = {
        'energy': None,
        'structure': None,
        'scf_misc': None,
        'scf_remote': None,
        'scf_retrieved': None,
        'dos_misc': None,
        'dos_remote': None,
        'dos': None,
        'projectors': None,
        'files': None,
        'pk': pk,
    }

    # Try to get outputs directly from BandsWorkChain
    if hasattr(node, 'outputs'):
        outputs = node.outputs

        # DOS ArrayData (direct BandsWorkChain output)
        if hasattr(outputs, 'dos'):
            result['dos'] = outputs.dos

        # Projectors (direct BandsWorkChain output)
        if hasattr(outputs, 'projectors'):
            result['projectors'] = outputs.projectors

        # SCF outputs (from modified BandsWorkChain)
        if hasattr(outputs, 'scf_misc'):
            misc_node = outputs.scf_misc
            if hasattr(misc_node, 'get_dict'):
                result['scf_misc'] = misc_node.get_dict()
                if result['energy'] is None:
                    result['energy'] = _extract_energy_from_misc(result['scf_misc'])

        if hasattr(outputs, 'scf_remote_folder'):
            result['scf_remote'] = outputs.scf_remote_folder

        if hasattr(outputs, 'scf_retrieved'):
            result['scf_retrieved'] = outputs.scf_retrieved

        # DOS outputs (from modified BandsWorkChain)
        if hasattr(outputs, 'dos_misc'):
            misc_node = outputs.dos_misc
            if hasattr(misc_node, 'get_dict'):
                result['dos_misc'] = misc_node.get_dict()

        if hasattr(outputs, 'dos_remote_folder'):
            result['dos_remote'] = outputs.dos_remote_folder

        if hasattr(outputs, 'dos_retrieved'):
            result['files'] = outputs.dos_retrieved

    # Get input structure
    if hasattr(node, 'inputs') and hasattr(node.inputs, 'structure'):
        result['structure'] = node.inputs.structure

    # Fallback: Traverse links to find SCF and DOS workchain outputs
    if result['scf_misc'] is None or result['dos_misc'] is None:
        _extract_dos_from_bandsworkchain(node, result)

    return result


def _extract_dos_from_bandsworkchain(node, result: dict) -> None:
    """
    Extract DOS results by traversing BandsWorkChain links.

    BandsWorkChain structure:
    - scf_workchain: VaspWorkChain with misc output
    - dos_workchain: VaspWorkChain with misc, dos, projectors, retrieved outputs

    Args:
        node: The BandsWorkChain node
        result: Result dict to populate (modified in place)
    """
    if not hasattr(node, 'base'):
        return

    # Traverse CALL_WORK links to find VaspWorkChain children
    called = node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        if hasattr(child_node, 'outputs'):
            outputs = child_node.outputs

            # Check if this is the SCF calc
            if 'scf' in link_label.lower():
                if hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        misc_dict = misc.get_dict()
                        result['scf_misc'] = misc_dict
                        if result['energy'] is None:
                            result['energy'] = _extract_energy_from_misc(misc_dict)

            # Check if this is the DOS calc
            if 'dos' in link_label.lower():
                if hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['dos_misc'] = misc.get_dict()

                # Get DOS ArrayData from child if not found at parent level
                if result['dos'] is None and hasattr(outputs, 'dos'):
                    result['dos'] = outputs.dos

                # Get projectors from child if not found at parent level
                if result['projectors'] is None and hasattr(outputs, 'projectors'):
                    result['projectors'] = outputs.projectors

                # Get retrieved files
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved

        # Recurse into nested workchains
        if hasattr(child_node, 'base'):
            _extract_dos_from_bandsworkchain(child_node, result)


def print_dos_results(pk: int) -> None:
    """
    Print a formatted summary of DOS calculation results.

    Args:
        pk: PK of the quick_dos BandsWorkChain
    """
    from .utils import get_status

    status = get_status(pk)
    print(f"DOS Calculation PK {pk}")
    print(f"  Status: {status}")

    if status != 'finished':
        print("  (Calculation not finished, results may be incomplete)")
        return

    results = get_dos_results(pk)

    if results['energy'] is not None:
        print(f"  SCF Energy: {results['energy']:.6f} eV")

    if results['structure'] is not None:
        struct = results['structure']
        formula = struct.get_formula()
        print(f"  Structure: {formula} (PK: {struct.pk})")

    if results['scf_misc'] is not None:
        scf_misc = results['scf_misc']
        run_status = scf_misc.get('run_status', {})
        converged = run_status.get('electronic_converged', 'N/A')
        print(f"  SCF converged: {converged}")

    if results['scf_remote'] is not None:
        print(f"  SCF remote_folder: PK {results['scf_remote'].pk}")

    if results['dos'] is not None:
        dos_node = results['dos']
        print(f"  DOS ArrayData: PK {dos_node.pk}")
        # Show available arrays in the DOS node
        if hasattr(dos_node, 'get_arraynames'):
            arrays = dos_node.get_arraynames()
            print(f"    Arrays: {', '.join(arrays)}")

    if results['projectors'] is not None:
        proj_node = results['projectors']
        print(f"  Projectors ArrayData: PK {proj_node.pk}")

    if results['dos_misc'] is not None:
        dos_misc = results['dos_misc']
        # Print some useful DOS info
        run_status = dos_misc.get('run_status', {})
        if run_status:
            print(f"  DOS run status: {run_status}")
        # Print band gap info if available
        band_props = dos_misc.get('band_properties', {})
        band_gap = band_props.get('band_gap', None)
        if band_gap is not None:
            is_direct = band_props.get('is_direct_gap', False)
            gap_type = "direct" if is_direct else "indirect"
            print(f"  Band gap: {band_gap:.4f} eV ({gap_type})")
        fermi = dos_misc.get('fermi_level', None)
        if fermi is not None:
            print(f"  Fermi level: {fermi:.4f} eV")

    if results['dos_remote'] is not None:
        print(f"  DOS remote_folder: PK {results['dos_remote'].pk}")

    if results['files'] is not None:
        files = results['files'].list_object_names()
        print(f"  DOS retrieved files: {', '.join(files)}")


def get_batch_dos_results(batch_result: dict) -> t.Dict[str, dict]:
    """
    Extract results from a quick_dos_batch calculation.

    Args:
        batch_result: Return value from quick_dos_batch containing
                     __workgraph_pk__ and __task_map__

    Returns:
        Dict mapping structure keys to DOS result dicts:
        {
            'key_0': {
                'energy': float (SCF energy in eV),
                'structure': StructureData (input structure),
                'scf_misc': dict (SCF VASP outputs),
                'dos_misc': dict (DOS VASP outputs),
                'dos': ArrayData (DOS data, if available),
                'projectors': ArrayData (projected DOS, if available),
                'files': FolderData (retrieved files from DOS calc),
                'pk': int (WorkGraph PK),
                'key': str (structure key),
            },
            'key_1': {...},
            ...
        }

    Example:
        >>> result = quick_dos_batch(structures={'pristine': s1, 'vacancy': s2}, ...)
        >>> # Wait for completion...
        >>> batch_results = get_batch_dos_results(result)
        >>> for key, dos_result in batch_results.items():
        ...     print(f"{key}: E = {dos_result['energy']:.4f} eV")
    """
    wg_pk = batch_result['__workgraph_pk__']
    task_map = batch_result['__task_map__']

    wg_node = orm.load_node(wg_pk)

    results = {}
    for key, task_info in task_map.items():
        # Initialize result dict for this structure
        result = {
            'energy': None,
            'structure': None,
            'scf_misc': None,
            'dos_misc': None,
            'scf_remote': None,
            'dos_remote': None,
            'files': None,
            'pk': wg_pk,
            'key': key,
        }

        scf_task_name = task_info.get('scf_task')
        dos_task_name = task_info['dos_task']

        # Try to access via WorkGraph outputs (exposed outputs)
        if hasattr(wg_node, 'outputs'):
            outputs = wg_node.outputs

            # SCF outputs
            scf_misc_attr = f'{key}_scf_misc'
            if hasattr(outputs, scf_misc_attr):
                misc_node = getattr(outputs, scf_misc_attr)
                if hasattr(misc_node, 'get_dict'):
                    result['scf_misc'] = misc_node.get_dict()
                    # Extract energy from SCF misc
                    if result['energy'] is None:
                        result['energy'] = _extract_energy_from_misc(result['scf_misc'])

            scf_remote_attr = f'{key}_scf_remote'
            if hasattr(outputs, scf_remote_attr):
                result['scf_remote'] = getattr(outputs, scf_remote_attr)

            # DOS outputs
            dos_misc_attr = f'{key}_dos_misc'
            if hasattr(outputs, dos_misc_attr):
                misc_node = getattr(outputs, dos_misc_attr)
                if hasattr(misc_node, 'get_dict'):
                    result['dos_misc'] = misc_node.get_dict()

            dos_remote_attr = f'{key}_dos_remote'
            if hasattr(outputs, dos_remote_attr):
                result['dos_remote'] = getattr(outputs, dos_remote_attr)

            dos_retrieved_attr = f'{key}_dos_retrieved'
            if hasattr(outputs, dos_retrieved_attr):
                result['files'] = getattr(outputs, dos_retrieved_attr)

        # Fallback: Traverse links to find VaspWorkChain outputs (for stored nodes)
        if result['energy'] is None or result['dos_misc'] is None:
            _extract_batch_dos_from_workgraph(wg_node, scf_task_name, dos_task_name, result)

        # Get input structure from DOS task
        if result['structure'] is None:
            _get_structure_from_task(wg_node, dos_task_name, result)

        results[key] = result

    return results


def _extract_batch_dos_from_workgraph(
    wg_node, scf_task_name: str, dos_task_name: str, result: dict
) -> None:
    """
    Extract DOS results by traversing WorkGraph links to find VaspWorkChain outputs.

    Args:
        wg_node: The WorkGraph node
        scf_task_name: Name of the SCF VaspWorkChain task
        dos_task_name: Name of the DOS VaspWorkChain task
        result: Result dict to populate (modified in place)
    """
    if not hasattr(wg_node, 'base'):
        return

    # Traverse CALL_WORK links to find VaspWorkChains
    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the SCF task
        if scf_task_name and (scf_task_name in link_label or link_label == scf_task_name):
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get SCF misc
                if result['scf_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['scf_misc'] = misc.get_dict()
                        if result['energy'] is None:
                            result['energy'] = _extract_energy_from_misc(result['scf_misc'])

                # Get SCF remote
                if result['scf_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['scf_remote'] = outputs.remote_folder

        # Check if this is the DOS task
        if dos_task_name in link_label or link_label == dos_task_name:
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get DOS misc
                if result['dos_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['dos_misc'] = misc.get_dict()

                # Get DOS remote
                if result['dos_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['dos_remote'] = outputs.remote_folder

                # Get DOS retrieved files
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved


def _get_structure_from_task(wg_node, task_name: str, result: dict) -> None:
    """
    Get the input structure from a task.

    Args:
        wg_node: The WorkGraph node
        task_name: Name of the task
        result: Result dict to populate (modified in place)
    """
    if not hasattr(wg_node, 'base'):
        return

    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        if task_name in link_label or link_label == task_name:
            if hasattr(child_node, 'inputs') and hasattr(child_node.inputs, 'structure'):
                result['structure'] = child_node.inputs.structure
                return


def print_batch_dos_results(batch_result: dict) -> None:
    """
    Print a formatted summary of batch DOS calculation results.

    Args:
        batch_result: Return value from quick_dos_batch
    """
    from .utils import get_status

    wg_pk = batch_result['__workgraph_pk__']
    status = get_status(wg_pk)

    print(f"Batch DOS Calculation - WorkGraph PK {wg_pk}")
    print(f"  Status: {status}")
    print()

    if status != 'finished':
        print("  (Calculation not finished, results may be incomplete)")
        print()

    results = get_batch_dos_results(batch_result)

    for key, dos_result in results.items():
        print(f"  [{key}]")

        if dos_result['energy'] is not None:
            print(f"    SCF Energy: {dos_result['energy']:.6f} eV")

        if dos_result['structure'] is not None:
            struct = dos_result['structure']
            formula = struct.get_formula()
            print(f"    Structure: {formula} (PK: {struct.pk})")

        if dos_result['scf_misc'] is not None:
            scf_misc = dos_result['scf_misc']
            print(f"    SCF misc: Dict (run_status: {scf_misc.get('run_status', 'N/A')})")

        if dos_result['dos_misc'] is not None:
            dos_misc = dos_result['dos_misc']
            print(f"    DOS misc: Dict (run_status: {dos_misc.get('run_status', 'N/A')})")

        if dos_result['scf_remote'] is not None:
            print(f"    SCF remote_folder: PK {dos_result['scf_remote'].pk}")

        if dos_result['dos_remote'] is not None:
            print(f"    DOS remote_folder: PK {dos_result['dos_remote'].pk}")

        if dos_result['files'] is not None:
            files = dos_result['files'].list_object_names()
            print(f"    DOS retrieved files: {', '.join(files)}")

        print()


def get_sequential_results(sequential_result: dict) -> t.Dict[str, dict]:
    """
    Extract results from all stages of a quick_vasp_sequential calculation.

    Args:
        sequential_result: Return value from quick_vasp_sequential containing
                          __workgraph_pk__ and __stage_names__

    Returns:
        Dict mapping stage names to result dicts:
        {
            'stage_name': {
                'energy': float (eV),
                'structure': StructureData (output structure),
                'misc': dict (VASP outputs),
                'remote': RemoteData (for restart),
                'files': FolderData (retrieved files),
                'pk': int (WorkGraph PK),
                'stage': str (stage name),
            },
            ...
        }

    Example:
        >>> result = quick_vasp_sequential(structure=s, stages=stages, ...)
        >>> # Wait for completion...
        >>> results = get_sequential_results(result)
        >>> for stage, data in results.items():
        ...     print(f"{stage}: E = {data['energy']:.4f} eV")
    """
    stage_names = sequential_result['__stage_names__']

    results = {}
    for stage_name in stage_names:
        results[stage_name] = get_stage_results(sequential_result, stage_name)

    return results


def get_stage_results(sequential_result: dict, stage_name: str) -> dict:
    """
    Extract results from a specific stage of a quick_vasp_sequential calculation.

    Args:
        sequential_result: Return value from quick_vasp_sequential
        stage_name: Name of the stage to extract results from

    Returns:
        For VASP stages:
            - energy: float (eV)
            - structure: StructureData (output structure)
            - misc: dict (VASP outputs)
            - remote: RemoteData (for restart)
            - files: FolderData (retrieved files)
            - pk: int (WorkGraph PK)
            - stage: str (stage name)
            - type: 'vasp'

        For DOS stages:
            - energy: float (SCF energy in eV)
            - scf_misc: dict (SCF VASP outputs)
            - scf_remote: RemoteData (SCF remote folder)
            - dos_misc: dict (DOS VASP outputs)
            - dos_remote: RemoteData (DOS remote folder)
            - files: FolderData (DOS retrieved files)
            - pk: int (WorkGraph PK)
            - stage: str (stage name)
            - type: 'dos'

    Example:
        >>> result = quick_vasp_sequential(structure=s, stages=stages, ...)
        >>> # Wait for completion...
        >>> stage_result = get_stage_results(result, 'relax_1x1_fine')
        >>> print(f"Energy: {stage_result['energy']:.4f} eV")
    """
    wg_pk = sequential_result['__workgraph_pk__']
    stage_names = sequential_result.get('__stage_names__', [])
    stage_types = sequential_result.get('__stage_types__', {})

    if stage_name not in stage_names:
        raise ValueError(
            f"Stage '{stage_name}' not found. Available stages: {stage_names}"
        )

    stage_type = stage_types.get(stage_name, 'vasp')
    wg_node = orm.load_node(wg_pk)

    if stage_type == 'dos':
        return _get_dos_stage_results(wg_node, wg_pk, stage_name)
    else:
        return _get_vasp_stage_results(wg_node, wg_pk, stage_name)


def _get_vasp_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a VASP stage."""
    result = {
        'energy': None,
        'structure': None,
        'misc': None,
        'remote': None,
        'files': None,
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'vasp',
    }

    # Try to access via WorkGraph outputs (exposed outputs)
    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        # Energy
        energy_attr = f'{stage_name}_energy'
        if hasattr(outputs, energy_attr):
            energy_node = getattr(outputs, energy_attr)
            if hasattr(energy_node, 'value'):
                result['energy'] = energy_node.value
            else:
                result['energy'] = float(energy_node)

        # Structure
        struct_attr = f'{stage_name}_structure'
        if hasattr(outputs, struct_attr):
            result['structure'] = getattr(outputs, struct_attr)

        # Misc
        misc_attr = f'{stage_name}_misc'
        if hasattr(outputs, misc_attr):
            misc_node = getattr(outputs, misc_attr)
            if hasattr(misc_node, 'get_dict'):
                result['misc'] = misc_node.get_dict()

        # Remote folder
        remote_attr = f'{stage_name}_remote'
        if hasattr(outputs, remote_attr):
            result['remote'] = getattr(outputs, remote_attr)

        # Retrieved files
        retrieved_attr = f'{stage_name}_retrieved'
        if hasattr(outputs, retrieved_attr):
            result['files'] = getattr(outputs, retrieved_attr)

    # Fallback: Traverse links to find VaspWorkChain outputs (for stored nodes)
    if result['energy'] is None or result['misc'] is None:
        _extract_sequential_stage_from_workgraph(wg_node, stage_name, result)

    # Extract energy from misc if not found directly
    if result['energy'] is None and result['misc'] is not None:
        result['energy'] = _extract_energy_from_misc(result['misc'])

    return result


def _get_dos_stage_results(wg_node, wg_pk: int, stage_name: str) -> dict:
    """Extract results from a DOS stage (BandsWorkChain)."""
    result = {
        'energy': None,
        'scf_misc': None,
        'scf_remote': None,
        'scf_retrieved': None,
        'dos_misc': None,
        'dos_remote': None,
        'files': None,  # DOS retrieved files (includes DOSCAR)
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'dos',
    }

    # Try to access via WorkGraph outputs (exposed outputs)
    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        # DOS ArrayData (from BandsWorkChain, if exposed)
        dos_attr = f'{stage_name}_dos'
        if hasattr(outputs, dos_attr):
            result['dos_arraydata'] = getattr(outputs, dos_attr)

        # Projectors ArrayData (from BandsWorkChain, if exposed)
        projectors_attr = f'{stage_name}_projectors'
        if hasattr(outputs, projectors_attr):
            result['projectors'] = getattr(outputs, projectors_attr)

        # SCF outputs (from modified BandsWorkChain)
        scf_misc_attr = f'{stage_name}_scf_misc'
        if hasattr(outputs, scf_misc_attr):
            misc_node = getattr(outputs, scf_misc_attr)
            if hasattr(misc_node, 'get_dict'):
                result['scf_misc'] = misc_node.get_dict()

        scf_remote_attr = f'{stage_name}_scf_remote'
        if hasattr(outputs, scf_remote_attr):
            result['scf_remote'] = getattr(outputs, scf_remote_attr)

        scf_retrieved_attr = f'{stage_name}_scf_retrieved'
        if hasattr(outputs, scf_retrieved_attr):
            result['scf_retrieved'] = getattr(outputs, scf_retrieved_attr)

        # DOS outputs (from modified BandsWorkChain)
        dos_misc_attr = f'{stage_name}_dos_misc'
        if hasattr(outputs, dos_misc_attr):
            misc_node = getattr(outputs, dos_misc_attr)
            if hasattr(misc_node, 'get_dict'):
                result['dos_misc'] = misc_node.get_dict()

        dos_remote_attr = f'{stage_name}_dos_remote'
        if hasattr(outputs, dos_remote_attr):
            result['dos_remote'] = getattr(outputs, dos_remote_attr)

        dos_retrieved_attr = f'{stage_name}_dos_retrieved'
        if hasattr(outputs, dos_retrieved_attr):
            result['files'] = getattr(outputs, dos_retrieved_attr)

    # Fallback: Traverse links to find BandsWorkChain and extract internal outputs
    if result['scf_misc'] is None or result['dos_misc'] is None:
        _extract_from_bands_workchain(wg_node, stage_name, result)

    # Extract energy from scf_misc if not found directly
    if result['energy'] is None and result['scf_misc'] is not None:
        result['energy'] = _extract_energy_from_misc(result['scf_misc'])

    return result


def _extract_from_bands_workchain(wg_node, stage_name: str, result: dict) -> None:
    """
    Extract DOS stage results by traversing to BandsWorkChain and its children.

    BandsWorkChain internally runs SCF and DOS VaspWorkChains. We traverse
    the links to find these and extract their outputs.

    Args:
        wg_node: The WorkGraph node
        stage_name: Name of the DOS stage
        result: Result dict to populate (modified in place)
    """
    if not hasattr(wg_node, 'base'):
        return

    bands_task_name = f'bands_{stage_name}'

    # Find the BandsWorkChain task
    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the BandsWorkChain task
        if bands_task_name in link_label or link_label == bands_task_name:
            # Get DOS and projectors from BandsWorkChain outputs (if available)
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs
                if 'dos_arraydata' not in result and hasattr(outputs, 'dos'):
                    result['dos_arraydata'] = outputs.dos
                if 'projectors' not in result and hasattr(outputs, 'projectors'):
                    result['projectors'] = outputs.projectors

            # Traverse into BandsWorkChain to find SCF and DOS VaspWorkChains
            if hasattr(child_node, 'base'):
                _extract_from_bands_children(child_node, result)


def _extract_from_bands_children(bands_node, result: dict) -> None:
    """
    Extract outputs from BandsWorkChain's child workchains (SCF and DOS).

    Args:
        bands_node: The BandsWorkChain node
        result: Result dict to populate (modified in place)
    """
    if not hasattr(bands_node, 'base'):
        return

    called = bands_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label.lower()

        if hasattr(child_node, 'outputs'):
            outputs = child_node.outputs

            # SCF workchain
            if 'scf' in link_label:
                if result['scf_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['scf_misc'] = misc.get_dict()
                if result['scf_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['scf_remote'] = outputs.remote_folder

            # DOS workchain
            if 'dos' in link_label and 'seekpath' not in link_label:
                if result['dos_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['dos_misc'] = misc.get_dict()
                if result['dos_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['dos_remote'] = outputs.remote_folder
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved


def _extract_dos_sequential_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """
    Extract DOS stage results by traversing WorkGraph links.

    Args:
        wg_node: The WorkGraph node
        stage_name: Name of the DOS stage to extract
        result: Result dict to populate (modified in place)
    """
    if not hasattr(wg_node, 'base'):
        return

    scf_task_name = f'scf_{stage_name}'
    dos_task_name = f'dos_{stage_name}'
    energy_task_name = f'energy_{stage_name}'

    # Traverse CALL_WORK links to find VaspWorkChains
    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the SCF task for this stage
        if scf_task_name in link_label or link_label == scf_task_name:
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get SCF misc
                if result['scf_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['scf_misc'] = misc.get_dict()

                # Get SCF remote folder
                if result['scf_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['scf_remote'] = outputs.remote_folder

        # Check if this is the DOS task for this stage
        if dos_task_name in link_label or link_label == dos_task_name:
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get DOS misc
                if result['dos_misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['dos_misc'] = misc.get_dict()

                # Get DOS remote folder
                if result['dos_remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['dos_remote'] = outputs.remote_folder

                # Get DOS retrieved files
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved

    # Traverse CALL_CALC links to find energy calcfunction
    called_calc = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
    for link in called_calc.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the energy task for this stage
        if energy_task_name in link_label or link_label == energy_task_name:
            # Get the output of the calcfunction
            created = child_node.base.links.get_outgoing(link_type=LinkType.CREATE)
            for out_link in created.all():
                if out_link.link_label == 'result':
                    energy_node = out_link.node
                    if hasattr(energy_node, 'value'):
                        result['energy'] = energy_node.value
                    break


def _extract_sequential_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """
    Extract stage results by traversing WorkGraph links to find VaspWorkChain outputs.

    Args:
        wg_node: The WorkGraph node
        stage_name: Name of the stage to extract
        result: Result dict to populate (modified in place)
    """
    if not hasattr(wg_node, 'base'):
        return

    vasp_task_name = f'vasp_{stage_name}'
    energy_task_name = f'energy_{stage_name}'

    # Traverse CALL_WORK links to find VaspWorkChain
    called = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
    for link in called.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the VASP task for this stage
        if vasp_task_name in link_label or link_label == vasp_task_name:
            if hasattr(child_node, 'outputs'):
                outputs = child_node.outputs

                # Get misc
                if result['misc'] is None and hasattr(outputs, 'misc'):
                    misc = outputs.misc
                    if hasattr(misc, 'get_dict'):
                        result['misc'] = misc.get_dict()

                # Get structure
                if result['structure'] is None and hasattr(outputs, 'structure'):
                    result['structure'] = outputs.structure

                # Get remote folder
                if result['remote'] is None and hasattr(outputs, 'remote_folder'):
                    result['remote'] = outputs.remote_folder

                # Get retrieved files
                if result['files'] is None and hasattr(outputs, 'retrieved'):
                    result['files'] = outputs.retrieved

    # Traverse CALL_CALC links to find energy calcfunction
    called_calc = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
    for link in called_calc.all():
        child_node = link.node
        link_label = link.link_label

        # Check if this is the energy task for this stage
        if energy_task_name in link_label or link_label == energy_task_name:
            # Get the output of the calcfunction
            created = child_node.base.links.get_outgoing(link_type=LinkType.CREATE)
            for out_link in created.all():
                if out_link.link_label == 'result':
                    energy_node = out_link.node
                    if hasattr(energy_node, 'value'):
                        result['energy'] = energy_node.value
                    break


def print_sequential_results(sequential_result: dict) -> None:
    """
    Print a formatted summary of sequential VASP calculation results.

    Args:
        sequential_result: Return value from quick_vasp_sequential
    """
    from .utils import get_status

    wg_pk = sequential_result['__workgraph_pk__']
    stage_names = sequential_result.get('__stage_names__', [])
    stage_types = sequential_result.get('__stage_types__', {})
    status = get_status(wg_pk)

    print(f"Sequential VASP Calculation - WorkGraph PK {wg_pk}")
    print(f"  Status: {status}")
    print(f"  Stages: {len(stage_names)}")
    print()

    if status != 'finished':
        print("  (Calculation not finished, results may be incomplete)")
        print()

    results = get_sequential_results(sequential_result)

    for i, (stage_name, stage_result) in enumerate(results.items(), 1):
        stage_type = stage_types.get(stage_name, 'vasp')

        if stage_type == 'dos':
            # DOS stage output (SCF + DOS VaspWorkChains)
            print(f"  [{i}] {stage_name} (DOS)")

            if stage_result['energy'] is not None:
                print(f"      SCF Energy: {stage_result['energy']:.6f} eV")

            if stage_result['scf_misc'] is not None:
                scf_misc = stage_result['scf_misc']
                run_status = scf_misc.get('run_status', {})
                converged = run_status.get('electronic_converged', 'N/A')
                print(f"      SCF converged: {converged}")

            if stage_result['dos_misc'] is not None:
                dos_misc = stage_result['dos_misc']
                band_props = dos_misc.get('band_properties', {})
                band_gap = band_props.get('band_gap', None)
                if band_gap is not None:
                    is_direct = band_props.get('is_direct_gap', False)
                    gap_type = "direct" if is_direct else "indirect"
                    print(f"      Band gap: {band_gap:.4f} eV ({gap_type})")
                fermi = dos_misc.get('fermi_level', None)
                if fermi is not None:
                    print(f"      Fermi level: {fermi:.4f} eV")

            if stage_result['scf_remote'] is not None:
                print(f"      SCF Remote folder: PK {stage_result['scf_remote'].pk}")

            if stage_result['dos_remote'] is not None:
                print(f"      DOS Remote folder: PK {stage_result['dos_remote'].pk}")

            if stage_result['files'] is not None:
                files = stage_result['files'].list_object_names()
                print(f"      DOS Retrieved: {', '.join(files)}")

        else:
            # VASP stage output
            print(f"  [{i}] {stage_name}")

            if stage_result['energy'] is not None:
                print(f"      Energy: {stage_result['energy']:.6f} eV")

            if stage_result['structure'] is not None:
                struct = stage_result['structure']
                formula = struct.get_formula()
                n_atoms = len(struct.sites)
                print(f"      Structure: {formula} ({n_atoms} atoms, PK: {struct.pk})")

            if stage_result['misc'] is not None:
                misc = stage_result['misc']
                run_status = misc.get('run_status', 'N/A')
                max_force = misc.get('maximum_force', None)
                force_str = f", max_force: {max_force:.4f} eV/Å" if max_force else ""
                print(f"      Status: {run_status}{force_str}")

            if stage_result['remote'] is not None:
                print(f"      Remote folder: PK {stage_result['remote'].pk}")

            if stage_result['files'] is not None:
                files = stage_result['files'].list_object_names()
                print(f"      Retrieved: {', '.join(files)}")

        print()
