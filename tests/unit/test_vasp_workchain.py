"""The work chain behind the VASP tasks accepts the upper-case INCAR tags of the recipes.

aiida-vasp 5 refuses upper-case keys on a stored ``Dict`` (a WorkGraph stores the
inputs of its tasks), so psteros runs ``PsterosVaspWorkChain``: ``VaspWorkChain``
with a case-tolerant ``parameters`` validator. Nothing else about it changes.
"""

from __future__ import annotations

import pytest

import psteros

pytest.importorskip("aiida_vasp")

from aiida import orm  # noqa: E402
from aiida.common.exceptions import InputValidationError  # noqa: E402
from aiida.plugins import WorkflowFactory  # noqa: E402

from psteros.backends.vasp_workchain import PsterosVaspWorkChain, case_tolerant_parameters_validator  # noqa: E402

pytestmark = [pytest.mark.tier2, pytest.mark.requires_aiida]

VaspWorkChain = WorkflowFactory("vasp.v2.vasp")
UPPER = {"incar": {"ENCUT": 400, "EDIFF": 1e-6, "IBRION": -1, "NSW": 0, "ISMEAR": 0, "LREAL": False}}


def test_aiida_vasp_refuses_the_upper_case_incar_of_a_stored_dict() -> None:
    stored = orm.Dict(UPPER).store()
    with pytest.raises(InputValidationError, match="lower case keys"):
        VaspWorkChain.spec().inputs["parameters"].validate(stored)


def test_psteros_work_chain_accepts_it_without_changing_the_node() -> None:
    stored = orm.Dict(UPPER).store()
    assert PsterosVaspWorkChain.spec().inputs["parameters"].validate(stored) is None
    assert stored.get_dict() == UPPER  # the recipe's Dict is not rewritten
    assert PsterosVaspWorkChain.spec().inputs["parameters"].validate(orm.Dict(UPPER)) is None


def test_lower_case_input_and_selective_dynamics_still_pass() -> None:
    lower = {"incar": {"encut": 400, "nsw": 0}, "dynamics": {"positions_dof": [[True, True, False]]}}
    assert case_tolerant_parameters_validator(orm.Dict(lower).store()) is None
    assert case_tolerant_parameters_validator(None) is None


def test_the_content_is_still_validated() -> None:
    with pytest.raises(InputValidationError, match="incar"):
        case_tolerant_parameters_validator(orm.Dict({"ENCUT": 400}).store())  # no incar namespace
    with pytest.raises(InputValidationError, match="massager"):
        case_tolerant_parameters_validator(orm.Dict({"incar": {"NOT_A_VASP_TAG": 1}}).store())


def test_the_interface_is_that_of_vasp_work_chain() -> None:
    assert issubclass(PsterosVaspWorkChain, VaspWorkChain)
    parent, child = VaspWorkChain.spec(), PsterosVaspWorkChain.spec()
    assert set(child.inputs) == set(parent.inputs)
    assert set(child.outputs) == set(parent.outputs)
    assert {code.status for code in child.exit_codes.values()} == {code.status for code in parent.exit_codes.values()}
    for name, port in parent.inputs.items():
        assert type(child.inputs[name]) is type(port) and child.inputs[name].required == port.required
        assert getattr(child.inputs[name], "valid_type", None) == getattr(port, "valid_type", None)


def test_vasp_tasks_run_the_psteros_work_chain(tmp_path) -> None:
    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(
        computer=computer, filepath_executable="/bin/true", label="vasp", default_calc_job_plugin="vasp.vasp"
    ).store()
    config = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="case",
        execution=psteros.ExecutionPolicy(computer=computer.label),
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer.label}", UPPER["incar"]),
    )
    graph = psteros.build_surface_workgraph({"o2": psteros.triplet_o2_cell()}, config)
    task = graph.tasks["o2_vasp"]
    assert task.inputs.parameters.value.get_dict() == UPPER  # same INCAR, same case as before
    executor = task.get_executor()
    assert (executor.module_path, executor.callable_name) == (PsterosVaspWorkChain.__module__, "PsterosVaspWorkChain")
