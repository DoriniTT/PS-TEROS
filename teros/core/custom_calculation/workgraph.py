"""WorkGraph builder for custom VASP calculations."""

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task


def build_custom_calculation_workgraph(
    structure,
    code_label,
    builder_inputs,
    name='custom_calc'
):
    """
    Build a WorkGraph for custom VASP calculations.

    Args:
        structure: StructureData or list of StructureData
        code_label: str, VASP code label (e.g., 'VASP-6.4.1@cluster02')
        builder_inputs: dict or list of dicts with VASP builder parameters
        name: str, WorkGraph name

    Returns:
        WorkGraph ready to submit
    """
    raise NotImplementedError("Function not yet implemented - placeholder for Task 1")


def get_custom_results(workgraph):
    """
    Extract results from completed custom calculation WorkGraph.

    Args:
        workgraph: Completed WorkGraph from build_custom_calculation_workgraph

    Returns:
        dict with energies, structures, misc outputs
    """
    raise NotImplementedError("Function not yet implemented - placeholder for Task 1")
