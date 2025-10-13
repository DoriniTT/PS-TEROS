"""
AIMD Module for PS-TEROS

Ab initio molecular dynamics calculations on slab structures.
Sequential AIMD stages with automatic restart chaining.
"""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, dynamic, namespace
from teros.core.slabs import extract_total_energy


def prepare_aimd_parameters(
    base_parameters: dict,
    temperature: float,
    steps: int
) -> dict:
    """
    Prepare INCAR parameters for a single AIMD stage.

    Takes base AIMD parameters and injects stage-specific values:
    - TEBEG and TEEND (both set to temperature for isothermal)
    - NSW (number of MD steps for this stage)

    Args:
        base_parameters: Base AIMD INCAR dict (IBRION=0, MDALGO, POTIM, etc.)
        temperature: Target temperature in K
        steps: Number of MD steps

    Returns:
        Complete INCAR dict for this stage
    """
    params = base_parameters.copy()
    params['TEBEG'] = temperature
    params['TEEND'] = temperature  # Isothermal
    params['NSW'] = steps
    return params
