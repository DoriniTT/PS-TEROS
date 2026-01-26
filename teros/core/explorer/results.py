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
            - dos_misc: dict (DOS VASP outputs)
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
        'dos_misc': None,
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

    # Get input structure
    if hasattr(node, 'inputs') and hasattr(node.inputs, 'structure'):
        result['structure'] = node.inputs.structure

    # Traverse links to find SCF and DOS workchain outputs
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
        if 'run_status' in dos_misc:
            print(f"  DOS run status: {dos_misc['run_status']}")

    if results['files'] is not None:
        files = results['files'].list_object_names()
        print(f"  Retrieved files: {', '.join(files)}")
