"""Public WorkGraph builder for surface-energy calculations."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from psteros.backends import add_qe_task, add_vasp_task
from psteros.backends.qe import as_aiida_structure
from psteros.config import QeCalculationConfig, SurfaceWorkflowConfig


def build_surface_workgraph(
    structures: Mapping[str, Any],
    config: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build or submit a serial AiiDA WorkGraph for labelled structures.

    The same typed surface recipe supports QE (the primary implementation) and
    VASP.  Each label becomes one backend task.  The graph's hard concurrency
    limit is taken from :class:`~psteros.config.ExecutionPolicy`, currently
    fixed at one for the Bohr campaign.

    Parameters
    ----------
    structures
        Mapping from stable calculation labels to pymatgen ``Structure``,
        AiiDA ``StructureData``, or an AiiDA node PK.
    config
        A typed backend and execution recipe.
    submit
        Submit the completed graph when true; otherwise return it unsubmitted.
    """

    if not structures:
        raise ValueError("at least one labelled structure is required")
    invalid = [label for label in structures if not label or not label.replace("_", "").replace("-", "").isalnum()]
    if invalid:
        raise ValueError(f"invalid calculation labels: {invalid}")

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=config.name)
    workgraph.max_number_jobs = config.execution.max_concurrent_jobs
    for label, structure in structures.items():
        override = config.role_overrides.get(label)
        if config.backend == "qe":
            task = add_qe_task(
                workgraph,
                label=label,
                structure=structure,
                config=config.calculation,
                execution=config.execution,
                override=override,
            )
            workgraph.outputs.__setattr__(f"{label}_parameters", task.outputs.output_parameters)
            workgraph.outputs.__setattr__(f"{label}_structure", task.outputs.output_structure)
            workgraph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
        else:
            task = add_vasp_task(
                workgraph,
                label=label,
                structure=structure,
                config=config.calculation,
                execution=config.execution,
                override=override,
            )
            workgraph.outputs.__setattr__(f"{label}_misc", task.outputs.misc)
            workgraph.outputs.__setattr__(f"{label}_structure", task.outputs.structure)
            workgraph.outputs.__setattr__(f"{label}_retrieved", task.outputs.retrieved)
    if submit:
        workgraph.submit()
    return workgraph


def build_qe_relax_static_workgraph(
    structures: Mapping[str, Any],
    relaxation: SurfaceWorkflowConfig,
    static: SurfaceWorkflowConfig,
    *,
    submit: bool = False,
) -> Any:
    """Build a serial QE relaxation-to-static WorkGraph for each structure.

    Each final SCF consumes the relaxed structure emitted by its preceding
    PwBaseWorkChain. Both recipes share one serial execution policy, so a graph
    containing multiple labels cannot release concurrent GPU jobs.
    """

    if not structures:
        raise ValueError("at least one labelled structure is required")
    if relaxation.backend != "qe" or static.backend != "qe":
        raise TypeError("build_qe_relax_static_workgraph requires QE recipes")
    if not isinstance(relaxation.calculation, QeCalculationConfig):
        raise TypeError("relaxation requires QeCalculationConfig")
    if not isinstance(static.calculation, QeCalculationConfig):
        raise TypeError("static requires QeCalculationConfig")
    if relaxation.execution != static.execution:
        raise ValueError("relaxation and static recipes must share one execution policy")

    from aiida_workgraph import WorkGraph

    workgraph = WorkGraph(name=f"{relaxation.name}_relax_static")
    workgraph.max_number_jobs = relaxation.execution.max_concurrent_jobs
    for label, source_structure in structures.items():
        if not label or not label.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"invalid calculation label: {label!r}")
        initial_structure = as_aiida_structure(source_structure)
        relax_task = add_qe_task(
            workgraph,
            label=f"{label}_relax",
            structure=initial_structure,
            config=relaxation.calculation,
            execution=relaxation.execution,
            override=relaxation.role_overrides.get(label),
        )
        static_task = add_qe_task(
            workgraph,
            label=f"{label}_static",
            structure=relax_task.outputs.output_structure,
            pseudo_structure=initial_structure,
            config=static.calculation,
            execution=static.execution,
            override=static.role_overrides.get(label),
        )
        workgraph.outputs.__setattr__(
            f"{label}_relaxed_structure", relax_task.outputs.output_structure
        )
        workgraph.outputs.__setattr__(
            f"{label}_relax_parameters", relax_task.outputs.output_parameters
        )
        workgraph.outputs.__setattr__(
            f"{label}_relax_retrieved", relax_task.outputs.retrieved
        )
        workgraph.outputs.__setattr__(
            f"{label}_static_parameters", static_task.outputs.output_parameters
        )
        workgraph.outputs.__setattr__(
            f"{label}_static_retrieved", static_task.outputs.retrieved
        )
    if submit:
        workgraph.submit()
    return workgraph
