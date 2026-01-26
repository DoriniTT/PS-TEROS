"""Utility functions for the explorer module."""

import os
import typing as t
from pathlib import Path

from aiida import orm
from aiida.engine import ProcessState


# Map AiiDA process states to simple status strings
STATUS_MAP = {
    ProcessState.CREATED: 'waiting',
    ProcessState.WAITING: 'waiting',
    ProcessState.RUNNING: 'running',
    ProcessState.FINISHED: 'finished',
    ProcessState.EXCEPTED: 'excepted',
    ProcessState.KILLED: 'killed',
}


def get_status(pk: int) -> str:
    """
    Get the status of a calculation by PK.

    Args:
        pk: The process PK

    Returns:
        Status string: 'waiting', 'running', 'finished', 'failed', 'excepted', or 'killed'
    """
    node = orm.load_node(pk)

    # Check if it's a process node
    if not hasattr(node, 'process_state'):
        raise ValueError(f"Node {pk} is not a process node")

    state = node.process_state

    # Map process state to simple string
    status = STATUS_MAP.get(state, 'unknown')

    # For finished processes, check exit status
    if status == 'finished':
        exit_status = node.exit_status
        if exit_status is not None and exit_status != 0:
            status = 'failed'

    return status


def get_restart_info(pk: int) -> dict:
    """
    Extract restart information from a previous calculation.

    Args:
        pk: PK of the previous calculation (WorkGraph or VaspWorkChain)

    Returns:
        dict with:
            - structure: StructureData (relaxed if available, else input)
            - remote_folder: RemoteData (for WAVECAR, CHGCAR restart)
    """
    node = orm.load_node(pk)

    result = {
        'structure': None,
        'remote_folder': None,
    }

    # Try to get outputs directly from the node
    if hasattr(node, 'outputs'):
        outputs = node.outputs

        # Get relaxed structure if available (NSW > 0 relaxation)
        if hasattr(outputs, 'structure'):
            result['structure'] = outputs.structure
        elif hasattr(outputs, 'misc'):
            # For WorkGraph, the structure might be nested
            pass

        # Get remote folder for restart
        if hasattr(outputs, 'remote_folder'):
            result['remote_folder'] = outputs.remote_folder

    # If no structure from outputs, try to get input structure
    if result['structure'] is None:
        # Traverse inputs to find structure
        if hasattr(node, 'inputs'):
            inputs = node.inputs
            if hasattr(inputs, 'structure'):
                result['structure'] = inputs.structure

    # For WorkGraph nodes, traverse CALL_CALC links to find the VASP calc
    if result['remote_folder'] is None:
        from aiida.common.links import LinkType
        if hasattr(node, 'base'):
            called = node.base.links.get_outgoing(link_type=LinkType.CALL_WORK)
            for link in called.all():
                # Find VASP workchain
                if 'vasp' in link.link_label.lower() or 'vasp_calc' in link.link_label.lower():
                    vasp_node = link.node
                    if hasattr(vasp_node, 'outputs'):
                        if hasattr(vasp_node.outputs, 'remote_folder'):
                            result['remote_folder'] = vasp_node.outputs.remote_folder
                        if hasattr(vasp_node.outputs, 'structure') and result['structure'] is None:
                            result['structure'] = vasp_node.outputs.structure
                    break

    return result


def export_files(
    pk: int,
    output_dir: str = './',
    files: t.List[str] = None,
) -> t.List[str]:
    """
    Export retrieved files from a calculation to a local directory.

    Args:
        pk: PK of the calculation
        output_dir: Directory to export files to (default: current directory)
        files: List of file names to export (default: all files)

    Returns:
        List of exported file paths
    """
    node = orm.load_node(pk)

    # Find the retrieved FolderData
    retrieved = None

    if hasattr(node, 'outputs'):
        if hasattr(node.outputs, 'retrieved'):
            retrieved = node.outputs.retrieved
        elif hasattr(node.outputs, 'files'):
            retrieved = node.outputs.files

    # For WorkGraph, traverse to find the VASP calc's retrieved folder
    if retrieved is None:
        from aiida.common.links import LinkType
        if hasattr(node, 'base'):
            # Try CALL_WORK first (for WorkGraphs calling WorkChains)
            for link_type in [LinkType.CALL_WORK, LinkType.CALL_CALC]:
                called = node.base.links.get_outgoing(link_type=link_type)
                for link in called.all():
                    child_node = link.node
                    if hasattr(child_node, 'outputs'):
                        if hasattr(child_node.outputs, 'retrieved'):
                            retrieved = child_node.outputs.retrieved
                            break
                if retrieved is not None:
                    break

    if retrieved is None:
        raise ValueError(f"Could not find retrieved files for PK {pk}")

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Get list of files to export
    available_files = retrieved.list_object_names()

    if files is None:
        files_to_export = available_files
    else:
        files_to_export = [f for f in files if f in available_files]
        missing = [f for f in files if f not in available_files]
        if missing:
            print(f"Warning: Files not found in retrieved: {missing}")
            print(f"Available files: {available_files}")

    # Export files
    exported = []
    for filename in files_to_export:
        content = retrieved.get_object_content(filename)
        output_file = output_path / filename
        if isinstance(content, bytes):
            output_file.write_bytes(content)
        else:
            output_file.write_text(content)
        exported.append(str(output_file))

    return exported


def list_calculations(
    name_pattern: str = None,
    limit: int = 20,
    all_users: bool = False,
) -> t.List[dict]:
    """
    List explorer calculations, optionally filtered by name pattern.

    Args:
        name_pattern: Optional glob pattern to filter by WorkGraph label (e.g., 'sno2*')
        limit: Maximum number of results to return
        all_users: Include calculations from all users

    Returns:
        List of dicts with: pk, label, state, ctime
    """
    from aiida.orm import QueryBuilder
    from aiida_workgraph.orm import WorkGraphNode

    qb = QueryBuilder()

    # Query for WorkGraph nodes
    filters = {}

    if name_pattern:
        # Convert glob pattern to SQL LIKE pattern
        like_pattern = name_pattern.replace('*', '%').replace('?', '_')
        filters['label'] = {'like': like_pattern}

    qb.append(
        WorkGraphNode,
        filters=filters if filters else None,
        project=['id', 'label', 'process_state', 'ctime'],
        tag='wg',
    )

    # Order by creation time, most recent first
    qb.order_by({'wg': {'ctime': 'desc'}})
    qb.limit(limit)

    results = []
    for row in qb.all():
        pk, label, state, ctime = row
        status = STATUS_MAP.get(state, 'unknown') if state else 'unknown'
        results.append({
            'pk': pk,
            'label': label or '',
            'state': status,
            'ctime': ctime.strftime('%Y-%m-%d %H:%M') if ctime else '',
        })

    return results


def prepare_restart_settings(
    restart_from: int,
    copy_wavecar: bool = True,
    copy_chgcar: bool = False,
) -> t.Tuple[orm.StructureData, dict]:
    """
    Prepare structure and settings for restarting from a previous calculation.

    Args:
        restart_from: PK of the previous calculation
        copy_wavecar: Whether to copy WAVECAR for restart (sets ISTART=1)
        copy_chgcar: Whether to copy CHGCAR for restart (sets ICHARG=1)

    Returns:
        Tuple of (structure, restart_settings) where restart_settings contains
        'folder' and 'incar_additions' keys
    """
    restart_info = get_restart_info(restart_from)

    if restart_info['structure'] is None:
        raise ValueError(f"Could not find structure from calculation PK {restart_from}")

    structure = restart_info['structure']

    restart_settings = {
        'folder': restart_info['remote_folder'],
        'incar_additions': {},
    }

    # Set INCAR parameters for restart
    if copy_wavecar and restart_info['remote_folder'] is not None:
        restart_settings['incar_additions']['ISTART'] = 1

    if copy_chgcar and restart_info['remote_folder'] is not None:
        restart_settings['incar_additions']['ICHARG'] = 1

    return structure, restart_settings
