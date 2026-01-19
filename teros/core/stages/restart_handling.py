"""
Stage 3: Restart Handling.

This module provides functions for extracting restart data from a previous
workgraph run, including RemoteData restart folders and slab structures.

The main functions are:
- ``handle_restart_from_node``: Extract restart folders and slab structures
  from a previous workgraph node.
- ``prepare_input_slabs``: Prepare input slabs, potentially overriding with
  restart data.

These functions implement the restart handling stage (Stage 3) of the
``build_core_workgraph()`` function.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from aiida import orm

__all__ = [
    "RestartData",
    "handle_restart_from_node",
    "prepare_input_slabs",
]

logger = logging.getLogger(__name__)


@dataclass
class RestartData:
    """Data extracted from a previous workgraph for restart.

    This dataclass holds the restart folders (RemoteData) and slab structures
    (StructureData) extracted from a previous workgraph run.

    Attributes:
        folders: Dictionary mapping slab labels (e.g., 'term_0') to RemoteData
            nodes containing the restart files. None if no restart data.
        slabs: Dictionary mapping slab labels to StructureData nodes. These
            may be relaxed structures (from relaxed_slabs output) or generated
            structures (from slab_structures output). None if no restart data.

    Example:
        >>> restart_data = handle_restart_from_node(19774)
        >>> restart_data.folders.keys()
        dict_keys(['term_0', 'term_1'])
        >>> restart_data.has_restart_data
        True
    """

    folders: Optional[Dict[str, Any]] = None  # Dict[str, orm.RemoteData]
    slabs: Optional[Dict[str, Any]] = None  # Dict[str, orm.StructureData]

    @property
    def has_restart_data(self) -> bool:
        """Check if any restart data is available.

        Returns:
            True if either folders or slabs are available.
        """
        return self.folders is not None or self.slabs is not None


def handle_restart_from_node(
    restart_from_node: Optional[int],
) -> Tuple[RestartData, Optional[Dict[str, Any]]]:
    """Extract restart data from a previous workgraph node.

    This function loads a previous workgraph node and extracts:
    1. RemoteData restart folders from slab_remote outputs
    2. Slab structures (relaxed or generated) from the outputs

    The extracted data can be used to restart slab relaxation calculations
    from where they left off, avoiding redundant computation.

    Args:
        restart_from_node: PK of the previous workgraph to restart from.
            If None, returns empty RestartData and None for input_slabs.

    Returns:
        Tuple of (RestartData, input_slabs_override):
            - RestartData containing extracted folders and slabs
            - input_slabs dictionary to override user-provided slabs,
              or None if no restart slabs were extracted

    Note:
        - Prefers relaxed_slabs output (for refinement workflows)
        - Falls back to slab_structures output (for restarting failed runs)
        - Returns empty RestartData if extraction fails (with logged error)

    Example:
        >>> restart_data, input_slabs = handle_restart_from_node(19774)
        >>> if restart_data.has_restart_data:
        ...     print(f"Loaded {len(restart_data.folders)} restart folders")
        Loaded 2 restart folders
    """
    if restart_from_node is None:
        return RestartData(), None

    from teros.core.slabs import extract_restart_folders_from_node

    restart_folders = None
    restart_slabs = None
    input_slabs_override = None

    logger.info("=" * 70)
    logger.info("RESTART MODE: Loading data from node %s", restart_from_node)
    logger.info("=" * 70)

    try:
        # Load the previous workgraph node
        prev_node = orm.load_node(restart_from_node)

        # Extract restart folders (RemoteData)
        restart_folders = extract_restart_folders_from_node(restart_from_node)
        logger.info("  ✓ Extracted restart folders: %s", list(restart_folders.keys()))

        # Extract slab structures from previous run
        # Prefer relaxed_slabs if available (for refinement),
        # otherwise use slab_structures (generated)
        if hasattr(prev_node.outputs, "relaxed_slabs"):
            restart_slabs = {}
            for label in prev_node.outputs.relaxed_slabs.keys():
                restart_slabs[label] = prev_node.outputs.relaxed_slabs[label]
            logger.info(
                "  ✓ Extracted RELAXED slab structures: %s",
                list(restart_slabs.keys()),
            )
        elif hasattr(prev_node.outputs, "slab_structures"):
            restart_slabs = {}
            for label in prev_node.outputs.slab_structures.keys():
                restart_slabs[label] = prev_node.outputs.slab_structures[label]
            logger.info("  ✓ Extracted slab structures: %s", list(restart_slabs.keys()))
        else:
            logger.warning(
                "Previous node has no slab_structures, will generate new slabs"
            )

        # Override input_slabs with slabs from previous run
        if restart_slabs:
            input_slabs_override = restart_slabs
            logger.info("  → Using slabs from previous run")

        logger.info("=" * 70)

    except ValueError as e:
        logger.error("Error extracting restart data: %s", e)
        logger.info("  → Proceeding without restart")
        restart_folders = None
        restart_slabs = None

    return RestartData(
        folders=restart_folders, slabs=restart_slabs
    ), input_slabs_override


def prepare_input_slabs(
    input_slabs: Optional[Dict[str, Any]],
    restart_data: RestartData,
) -> Tuple[Optional[Dict[str, Any]], bool]:
    """Prepare input slabs, potentially overriding with restart data.

    This function determines the final input_slabs dictionary to use,
    considering both user-provided slabs and restart data.

    Args:
        input_slabs: User-provided input slab structures. May be None.
        restart_data: RestartData from handle_restart_from_node(). The
            restart_data.slabs will override input_slabs if available.

    Returns:
        Tuple of (final_input_slabs, use_input_slabs):
            - final_input_slabs: The slab dictionary to use (may be from
              restart_data.slabs or the original input_slabs)
            - use_input_slabs: Boolean flag indicating whether to use
              input slabs instead of generating new ones

    Note:
        The use_input_slabs flag is True if input_slabs is provided
        (either originally or via restart) and contains at least one slab.
        This flag is used downstream to decide between slab generation
        and using pre-defined slab structures.

    Example:
        >>> # With restart data
        >>> restart_data = RestartData(slabs={'term_0': slab_node})
        >>> final_slabs, use_slabs = prepare_input_slabs(None, restart_data)
        >>> use_slabs
        True

        >>> # Without restart data
        >>> final_slabs, use_slabs = prepare_input_slabs(None, RestartData())
        >>> use_slabs
        False
    """
    # Override input_slabs with restart data if available
    if restart_data.slabs is not None:
        final_slabs = restart_data.slabs
    else:
        final_slabs = input_slabs

    # Special handling for input_slabs: stored nodes can't be passed through @task.graph
    use_input_slabs = final_slabs is not None and len(final_slabs) > 0

    return final_slabs, use_input_slabs
