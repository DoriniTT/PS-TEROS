"""Tier-2 tests: the QE relaxation stage accepts a vc-relax that ended with exit 501.

aiida-quantumespresso reports 501 when a vc-relax converged but its final SCF
exceeded the thresholds.  aiida-workgraph treats any non-zero exit as failure
and would skip the static SCF, so psteros relaxes with ``PwRelaxStageWorkChain``.
"""

from __future__ import annotations

import io

import pytest

pytest.importorskip("aiida_quantumespresso")

import psteros

pytestmark = [pytest.mark.tier2, pytest.mark.requires_aiida]


def _upf(element: str, z_valence: float) -> bytes:
    return (
        '<UPF version="2.0.1">\n'
        f'<PP_HEADER element="{element}" pseudo_type="US" functional="PBE" z_valence="{z_valence}" />\n'
        "</UPF>\n"
    ).encode()


@pytest.fixture
def qe_inputs(tmp_path):
    """A local code, a pseudo family and fake (header-only) pseudopotentials."""
    from aiida import orm
    from aiida_pseudo.data.pseudo import UpfData
    from aiida_pseudo.groups.family import PseudoPotentialFamily

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    code = orm.InstalledCode(
        computer=computer, filepath_executable="/bin/true", label="pw",
        default_calc_job_plugin="quantumespresso.pw",
    ).store()
    (tmp_path / "pseudos").mkdir()
    for element, z in (("Sn", 14.0), ("O", 6.0)):
        (tmp_path / "pseudos" / f"{element}.upf").write_bytes(_upf(element, z))
    family = PseudoPotentialFamily.create_from_folder(
        tmp_path / "pseudos", f"fake-{tmp_path.name}", pseudo_type=UpfData
    )
    structure = orm.StructureData(pymatgen=psteros.rutile_sno2_bulk())
    return {
        "computer": computer,
        "code": code,
        "family": family,
        "workchain_inputs": {
            "pw": {
                "code": code,
                "structure": structure,
                "parameters": orm.Dict({
                    "CONTROL": {"calculation": "vc-relax"}, "SYSTEM": {"ecutwfc": 30.0}, "ELECTRONS": {},
                }),
                "pseudos": {
                    element: UpfData(io.BytesIO(_upf(element, z)), filename=f"{element}.upf")
                    for element, z in (("Sn", 14.0), ("O", 6.0))
                },
                "metadata": {"options": {
                    "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 1},
                    "max_wallclock_seconds": 60,
                }},
            },
            "kpoints_distance": orm.Float(0.5),
        },
    }


def _inspect_final_scf_warning(process_class, qe_inputs):
    """Run the restart logic of ``process_class`` on a finished exit-501 PwCalculation."""
    from aiida import orm
    from aiida.engine import ProcessState
    from aiida.engine.utils import instantiate_process
    from aiida.manage import get_manager

    process = instantiate_process(get_manager().get_runner(), process_class, **qe_inputs["workchain_inputs"])
    process.setup()
    calculation = orm.CalcJobNode(
        computer=qe_inputs["computer"], process_type="aiida.calculations:quantumespresso.pw"
    )
    calculation.set_process_state(ProcessState.FINISHED)
    calculation.set_exit_status(501)
    calculation.store()
    process.ctx.children = [calculation]
    process.ctx.iteration = 1
    return process, process.inspect_process()


def test_stock_pw_base_workchain_reports_exit_501(qe_inputs) -> None:
    from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

    process, exit_code = _inspect_final_scf_warning(PwBaseWorkChain, qe_inputs)
    assert exit_code is not None and exit_code.status == 501


def test_relax_stage_accepts_the_relaxed_structure(qe_inputs) -> None:
    from psteros.backends.qe_workchains import PwRelaxStageWorkChain

    process, exit_code = _inspect_final_scf_warning(PwRelaxStageWorkChain, qe_inputs)
    # A zero exit code lets the outline continue; the restart loop then stops and
    # ``results`` attaches the outputs, so the work chain finishes with exit 0.
    assert exit_code.status == 0
    assert process.ctx.is_finished
    assert not process.should_run_process()
    assert process.results() is None


def test_relax_static_graph_uses_the_relax_stage_only_for_relaxations(qe_inputs) -> None:
    execution = psteros.ExecutionPolicy(computer="local", queue="q")

    def recipe(calculation: str) -> psteros.SurfaceWorkflowConfig:
        return psteros.SurfaceWorkflowConfig(
            backend="qe",
            calculation=psteros.QeCalculationConfig(
                qe_inputs["code"].full_label, qe_inputs["family"].label,
                {"CONTROL": {"calculation": calculation}, "SYSTEM": {"ecutwfc": 30.0}, "ELECTRONS": {}},
            ),
            execution=execution,
        )

    graph = psteros.build_qe_relax_static_workgraph(
        {"bulk": psteros.rutile_sno2_bulk()}, recipe("vc-relax"), recipe("scf")
    )
    relax = graph.tasks["bulk_relax_qe"].get_executor()
    static = graph.tasks["bulk_static_qe"].get_executor()
    assert (relax.module_path, relax.callable_name) == (
        "psteros.backends.qe_workchains", "PwRelaxStageWorkChain",
    )
    assert static.callable_name == "PwBaseWorkChain"

    single = psteros.build_surface_workgraph({"bulk": psteros.rutile_sno2_bulk()}, recipe("scf"))
    assert single.tasks["bulk_qe"].get_executor().callable_name == "PwBaseWorkChain"
